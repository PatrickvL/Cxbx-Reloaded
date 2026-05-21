// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
// ******************************************************************
// *
// *  This file is part of the Cxbx project.
// *
// *  Cxbx and Cxbe are free software; you can redistribute them
// *  and/or modify them under the terms of the GNU General Public
// *  License as published by the Free Software Foundation; either
// *  version 2 of the license, or (at your option) any later version.
// *
// *  This program is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// *  GNU General Public License for more details.
// *
// *  You should have recieved a copy of the GNU General Public License
// *  along with this program; see the file COPYING.
// *  If not, write to the Free Software Foundation, Inc.,
// *  59 Temple Place - Suite 330, Bostom, MA 02111-1307, USA.
// *
// *  (c) 2002-2003 Aaron Robinson <caustik@caustik.com>
// * 
// *  This file is heavily based on code from XQEMU
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a.c
// *  Copyright (c) 2012 espes
// *  Copyright (c) 2015 Jannik Vogel
// *  Copyright (c) 2018 Matt Borgerson
// *
// *  Contributions for Cxbx-Reloaded
// *  Copyright (c) 2017-2018 Luke Usher <luke.usher@outlook.com>
// *  Copyright (c) 2018 Patrick van Logchem <pvanlogchem@gmail.com>
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::NV2A

#include <core/kernel/exports/xboxkrnl.h> // For PKINTERRUPT, etc.

#ifdef _MSC_VER                         // Check if MS Visual C compiler
#endif

#include <string> // For std::string
#include <distorm.h> // For uint32_t
#include <process.h> // For __beginthreadex(), etc.

#include "core\kernel\init\CxbxKrnl.h" // For XBOX_MEMORY_SIZE, DWORD, etc
#include "core\kernel\support\Emu.h"
#include "core\kernel\support\NativeHandle.h"
#include "core\kernel\exports\EmuKrnl.h"
#include <backends/imgui_impl_win32.h>
#include "core/common/video/RenderBase.hpp"
#include "core\hle\Intercept.hpp"
#include "common/win32/Threads.h"
#include "Logging.h"
#include "Timer.h"

#include "vga.h"
#include "nv2a.h" // For NV2AState
#include "nv2a_int.h" // from https://github.com/espes/xqemu/tree/xbox/hw/xbox
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11_Profiler.h"
#include <cassert>

// glib types
typedef char gchar;
typedef int gint;
typedef unsigned int guint;
typedef unsigned int guint32;
typedef const void *gconstpointer;
typedef gint   gboolean;
typedef void* gpointer;

typedef guint32 GQuark;

typedef struct _GError GError;

struct _GError
{
	GQuark       domain;
	gint         code;
	gchar       *message;
};

static void update_irq(NV2AState *d)
{
	/* PGRAPH - Auto-ack CONTEXT_SWITCH only. The PULLER thread uses
	 * CONTEXT_SWITCH as an internal synchronization mechanism and waits on
	 * interrupt_cond for the ack. We auto-ack it here because the Xbox
	 * miniport ISR may not handle it properly. Other PGRAPH interrupts
	 * (ERROR, NOTIFY) must route through the real ISR so that
	 * D3DDevice_InsertCallback works correctly. */
	if (d->pgraph.pending_interrupts & NV_PGRAPH_INTR_CONTEXT_SWITCH) {
		d->pgraph.pending_interrupts &= ~NV_PGRAPH_INTR_CONTEXT_SWITCH;
		qemu_cond_broadcast(&d->pgraph.interrupt_cond);
	}

	/* Compute live PMC interrupt status from sub-units.
	 * PMC_INTR_0 on real NV2A is a read-only register that reflects live
	 * sub-unit status. We don't need to cache pmc.pending_interrupts for
	 * read purposes (the READ handler computes it live), but we still
	 * need to know if anything is pending for Assert(true/false). */
	bool any_pending = false;
	if (d->pfifo.pending_interrupts & d->pfifo.enabled_interrupts)
		any_pending = true;
	if (d->pgraph.pending_interrupts & d->pgraph.enabled_interrupts)
		any_pending = true;
	if (d->pcrtc.pending_interrupts & d->pcrtc.enabled_interrupts)
		any_pending = true;
	if (d->pvideo.pending_interrupts & d->pvideo.enabled_interrupts)
		any_pending = true;
	if (d->ptimer.pending_interrupts & d->ptimer.enabled_interrupts)
		any_pending = true;

	if (any_pending && d->pmc.enabled_interrupts) {
		HalSystemInterrupts[3].Assert(true);
		// Wake the DPC thread so it can fire the ISR for non-VBlank
		// interrupts (e.g. PGRAPH INTR_ERROR from D3DDevice_InsertCallback).
		extern void KeSignalVBlankPending();
		KeSignalVBlankPending();
	}
	else {
		HalSystemInterrupts[3].Assert(false);
		// PGRAPH INTR_ERROR (InsertCallback) stalls the GPU pipeline on real
		// hardware until the CPU acknowledges it. If the ISR temporarily
		// disabled NV_PMC_INTR_EN_0 (standard ISR prologue), we still need
		// the DPC thread awake to re-fire the ISR once PMC is re-enabled.
		// Without this, the puller blocks forever waiting for an ack that
		// never comes because the DPC thread is asleep.
		if (d->pgraph.pending_interrupts & NV_PGRAPH_INTR_ERROR) {
			extern void KeSignalVBlankPending();
			KeSignalVBlankPending();
		}
	}
}


#include "EmuNV2A_DEBUG.cpp"


#define DEBUG_READ32(DEV)              EmuLog(LOG_LEVEL::DEBUG, "Rd32 NV2A " #DEV "(0x%08X) = 0x%08X [Handled %s]", addr, result, DebugNV_##DEV(addr))
#define DEBUG_READ32_UNHANDLED(DEV)  { EmuLog(LOG_LEVEL::DEBUG, "Rd32 NV2A " #DEV "(0x%08X) = 0x%08X [Unhandled %s]", addr, result, DebugNV_##DEV(addr)); return result; }

#define DEBUG_WRITE32(DEV)             EmuLog(LOG_LEVEL::DEBUG, "Wr32 NV2A " #DEV "(0x%08X, 0x%08X) [Handled %s]", addr, value, DebugNV_##DEV(addr))
#define DEBUG_WRITE32_UNHANDLED(DEV) { EmuLog(LOG_LEVEL::DEBUG, "Wr32 NV2A " #DEV "(0x%08X, 0x%08X) [Unhandled %s]", addr, value, DebugNV_##DEV(addr)); return; }

#define DEVICE_READ32_NAME(DEV) EmuNV2A_##DEV##_Read32
#define DEVICE_READ32(DEV) uint32_t DEVICE_READ32_NAME(DEV)(NV2AState *d, xbox::addr_xt addr)
#define DEVICE_READ32_SWITCH() uint32_t result = 0; switch (addr) 
#define DEVICE_READ32_REG(dev) result = d->dev.regs[RI(addr)]
#define DEVICE_READ32_END(DEV) DEBUG_READ32(DEV); return result

#define DEVICE_WRITE32_NAME(DEV) EmuNV2A_##DEV##_Write32
#define DEVICE_WRITE32(DEV) void DEVICE_WRITE32_NAME(DEV)(NV2AState *d, xbox::addr_xt addr, uint32_t value)
#define DEVICE_WRITE32_REG(dev) d->dev.regs[RI(addr)] = value
#define DEVICE_WRITE32_END(DEV) DEBUG_WRITE32(DEV)

static inline uint32_t ldl_le_p(const void *p)
{
	return *(uint32_t*)p;
}

static inline void stq_le_p(uint64_t *p, uint64_t v)
{
	*p = v;
}

static inline void stl_le_p(uint32_t *p, uint32_t v)
{
	*p = v;
}

static DMAObject nv_dma_load(NV2AState *d, xbox::addr_xt dma_obj_address)
{
	assert(dma_obj_address < d->pramin.ramin_size);

	uint32_t *dma_obj = (uint32_t*)(d->pramin.ramin_ptr + dma_obj_address);
	uint32_t flags = ldl_le_p(dma_obj);
	uint32_t limit = ldl_le_p(dma_obj + 1);
	uint32_t frame = ldl_le_p(dma_obj + 2);

	DMAObject object;
	object.dma_class = GET_MASK(flags, NV_DMA_CLASS);
	object.dma_target = GET_MASK(flags, NV_DMA_TARGET);
	object.address = (frame & NV_DMA_ADDRESS) | GET_MASK(flags, NV_DMA_ADJUST);
	object.limit = limit;

	return object;
}

static void *nv_dma_map(NV2AState *d, xbox::addr_xt dma_obj_address, xbox::addr_xt *len)
{
//	assert(dma_obj_address < d->pramin.ramin_size);

	DMAObject dma = nv_dma_load(d, dma_obj_address);

	/* TODO: Handle targets and classes properly */
	NV2A_DPRINTF("dma_map %x, %x, %x %x\n",
		dma.dma_class, dma.dma_target, dma.address, dma.limit);

	dma.address &= 0x07FFFFFF;

	// assert(dma.address + dma.limit < memory_region_size(d->vram));
	*len = dma.limit;
	return d->vram_ptr + dma.address;
//	return (void*)(PHYSICAL_MAP_BASE  + dma.address);
}

uint32_t NV2ADevice::ResolveDmaBaseAddress(NV2AState *d, xbox::addr_xt dma_obj_address)
{
	if (dma_obj_address == 0 || dma_obj_address >= d->pramin.ramin_size)
		return 0;

	uint32_t *dma_obj = (uint32_t*)(d->pramin.ramin_ptr + dma_obj_address);
	uint32_t flags = ldl_le_p(dma_obj);
	uint32_t frame = ldl_le_p(dma_obj + 2);
	uint32_t address = (frame & NV_DMA_ADDRESS) | GET_MASK(flags, NV_DMA_ADJUST);
	return address & 0x07FFFFFF; // Mask to 128 MB physical address space
}

#include "EmuNV2A_PBUS.cpp"
#include "EmuNV2A_PCRTC.cpp"
#include "EmuNV2A_PFB.cpp"
#include "nv2a_method_table.h"
#include "EmuNV2A_PGRAPH.cpp"
#include "EmuNV2A_PFIFO.cpp"
#include "EmuNV2A_PMC.cpp"
#include "EmuNV2A_PRAMDAC.cpp"
#include "EmuNV2A_PRMCIO.cpp"
#include "EmuNV2A_PRMVIO.cpp"
#include "EmuNV2A_PTIMER.cpp"
#include "EmuNV2A_PVIDEO.cpp"
#include "EmuNV2A_USER.cpp"

#include "EmuNV2A_PRMA.cpp"
#include "EmuNV2A_PCOUNTER.cpp"
#include "EmuNV2A_PVPE.cpp"
#include "EmuNV2A_PTV.cpp"
#include "EmuNV2A_PRMFB.cpp"
#include "EmuNV2A_PSTRAPS.cpp"
#include "EmuNV2A_PRMDIO.cpp"
#include "EmuNV2A_PRAMIN.cpp"

const NV2ABlockInfo regions[] = { // blocktable

// Note : Avoid designated initializers to facilitate C++ builds
#define ENTRY(OFFSET, SIZE, NAME) \
	{ \
        #NAME, OFFSET, SIZE, \
        { DEVICE_READ32_NAME(NAME), DEVICE_WRITE32_NAME(NAME) }, \
    }, \

#define ENTRY_MIRROR(OFFSET, SIZE, NAME, MIRROR) \
	{ \
        #NAME, OFFSET, SIZE, \
        { DEVICE_READ32_NAME(MIRROR), DEVICE_WRITE32_NAME(MIRROR) }, \
    }, \
	
#define ENTRY_END(OFFSET, SIZE, NAME) \
	{ \
        #NAME, OFFSET, SIZE, \
        { nullptr, nullptr }, \
    }, \

	/* card master control */
	ENTRY(0x000000, 0x001000, PMC)
	/* bus control */
	ENTRY(0x001000, 0x001000, PBUS)
	/* MMIO and DMA FIFO submission to PGRAPH and VPE */
	ENTRY(0x002000, 0x002000, PFIFO)
	/* access to BAR0/BAR1 from real mode */
	ENTRY(0x007000, 0x001000, PRMA)
	/* video overlay */
	ENTRY(0x008000, 0x001000, PVIDEO)
	/* time measurement and time-based alarms */
	ENTRY(0x009000, 0x001000, PTIMER)
	/* performance monitoring counters */
	ENTRY(0x00a000, 0x001000, PCOUNTER)
	/* MPEG2 decoding engine */
	ENTRY(0x00b000, 0x001000, PVPE)
	/* TV encoder */
	ENTRY(0x00d000, 0x001000, PTV)
	/* aliases VGA memory window */
	ENTRY(0x0a0000, 0x020000, PRMFB)
	/* aliases VGA sequencer and graphics controller registers */
	ENTRY(0x0c0000, 0x008000, PRMVIO) // Size was 0x001000
	/* memory interface */
	ENTRY(0x100000, 0x001000, PFB)
	/* straps readout / override */
	ENTRY(0x101000, 0x001000, PSTRAPS)
	/* accelerated 2d/3d drawing engine */
	ENTRY(0x400000, 0x002000, PGRAPH)
	/* more CRTC controls */
	ENTRY(0x600000, 0x001000, PCRTC)
	/* aliases VGA CRTC and attribute controller registers */
	ENTRY(0x601000, 0x001000, PRMCIO)
	/* RAMDAC, cursor, and PLL control */
	ENTRY(0x680000, 0x001000, PRAMDAC)
	/* aliases VGA palette registers */
	ENTRY(0x681000, 0x001000, PRMDIO)
	/* RAMIN access */
	ENTRY(0x700000, 0x100000, PRAMIN)
	/* PFIFO MMIO and DMA submission area */
	ENTRY(0x800000, 0x400000, USER) // Size was 0x800000
	/* UREMAP User area mirror - TODO : Confirm */
	ENTRY_MIRROR(0xC00000, 0x400000, UREMAP, USER) // NOTE : Mirror of USER
	/* Terminating entry */
	ENTRY_END(0xFFFFFF, 0x000000, END)
#undef ENTRY
#undef ENTRY_MIRROR
#undef ENTRY_END
};

const NV2ABlockInfo* EmuNV2A_Block(xbox::addr_xt addr)
{
	// Find the block in the block table
	const NV2ABlockInfo* block = &regions[0];
	int i = 0;

	while (block->size > 0) {
		if (addr >= block->offset && addr < block->offset + block->size) {
			return block;
		}

		block = &regions[++i];
	}

	return nullptr;
}

void NV2ADevice::UpdateHostDisplay(NV2AState *d)
{
	g_renderbase->UpdateFPSCounter();
}

void nv2a_vblank_interrupt(void *opaque)
{
	NV2AState *d = static_cast<NV2AState *>(opaque);

	if (!d->exiting) [[likely]] {
		// Increment VBlank counter (used by NV_PCRTC_RASTER FIELD bit for interlace detection)
		d->pcrtc.vblank_count++;

		// Signal that a VBlank occurred. Don't touch pcrtc.pending_interrupts here!
		// The main thread will set/clear it atomically around the ISR call to prevent
		// the timer from re-asserting it while the ISR is processing.
		d->vblank_pending.test_and_set();

		// Wake the main thread so it can dispatch the ISR:
		extern void KeSignalVBlankPending();
		KeSignalVBlankPending();

		// TODO: We should swap here for the purposes of supporting overlays + direct framebuffer access
		// But it causes crashes on AMD hardware for reasons currently unknown...
		//NV2ADevice::UpdateHostDisplay(d);
	}
}

// See NV2ABlockInfo regions[] PRAMIN
#define NV_PRAMIN_ADDR   0x00700000
#define NV_PRAMIN_SIZE              0x100000

// Flat MMIO backing storage — 16 MiB reserved, only engine block pages committed.
uint8_t* g_pNV2AMMIO = nullptr;

static void CxbxAllocateFlatMMIO(NV2AState *d)
{
	// Reserve 16 MiB virtual address space (no physical memory committed yet)
	g_pNV2AMMIO = (uint8_t*)VirtualAlloc(nullptr, NV2A_MMIO_TOTAL_SIZE,
		MEM_RESERVE, PAGE_NOACCESS);
	if (!g_pNV2AMMIO) {
		CxbxrAbort("VirtualAlloc failed to reserve NV2A flat MMIO buffer (16 MiB). Error 0x%08X", GetLastError());
	}

	// Commit only the pages where actual register blocks reside
	struct { uint32_t offset; uint32_t size; } blocks[] = {
		{ NV2A_MMIO_OFF_PMC,     NV_PMC_REGS_BYTES },      // 4 KB
		{ NV2A_MMIO_OFF_PFIFO,   NV_PFIFO_REGS_BYTES },    // 8 KB
		{ NV2A_MMIO_OFF_PVIDEO,  NV_PVIDEO_REGS_BYTES },   // 4 KB
		{ NV2A_MMIO_OFF_PTIMER,  NV_PTIMER_REGS_BYTES },   // 4 KB
		{ NV2A_MMIO_OFF_PFB,     NV_PFB_REGS_BYTES },      // 4 KB
		{ NV2A_MMIO_OFF_PGRAPH,  NV_PGRAPH_REGS_BYTES },   // 8 KB
		{ NV2A_MMIO_OFF_PCRTC,   NV_PCRTC_REGS_BYTES },    // 4 KB
		{ NV2A_MMIO_OFF_PRAMDAC, NV_PRAMDAC_REGS_BYTES },  // 4 KB
	};

	for (auto& blk : blocks) {
		LPVOID ret = VirtualAlloc(g_pNV2AMMIO + blk.offset, blk.size,
			MEM_COMMIT, PAGE_READWRITE);
		if (!ret) {
			CxbxrAbort("VirtualAlloc failed to commit NV2A MMIO block at offset 0x%06X (size %u). Error 0x%08X",
				blk.offset, blk.size, GetLastError());
		}
	}

	// Point each struct's regs pointer into the flat buffer
	d->pmc.regs     = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PMC);
	d->pfifo.regs   = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PFIFO);
	d->pvideo.regs  = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PVIDEO);
	d->ptimer.regs  = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PTIMER);
	d->pfb.regs     = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PFB);
	d->pgraph.regs  = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PGRAPH);
	d->pcrtc.regs   = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PCRTC);
	d->pramdac.regs = (uint32_t*)(g_pNV2AMMIO + NV2A_MMIO_OFF_PRAMDAC);

	printf("[0x%.4X] INIT: NV2A flat MMIO buffer reserved at %p (16 MiB, %zu KB committed)\n",
		GetCurrentThreadId(), g_pNV2AMMIO,
		(NV_PMC_REGS_BYTES + NV_PFIFO_REGS_BYTES + NV_PVIDEO_REGS_BYTES +
		 NV_PTIMER_REGS_BYTES + NV_PFB_REGS_BYTES + NV_PGRAPH_REGS_BYTES +
		 NV_PCRTC_REGS_BYTES + NV_PRAMDAC_REGS_BYTES) / 1024);
}

void CxbxReserveNV2AMemory(NV2AState *d)
{
	// The NV2A memory was reserved already by the loader!

	printf("[0x%.4X] INIT: Reserved %d MiB of Xbox NV2A memory at 0x%.8X to 0x%.8X\n",
		GetCurrentThreadId(), NV2A_SIZE / ONE_MB, NV2A_ADDR, NV2A_ADDR + NV2A_SIZE - 1);

	// Allocate PRAMIN Region (the loader only reserved this region, it still needs to be committed!)
	// We are looping here because memory-reservation happens in 64 KiB increments
	d->pramin.ramin_size = NV_PRAMIN_SIZE;
	d->pramin.ramin_ptr = (uint8_t*)(NV2A_ADDR + NV_PRAMIN_ADDR);
	for (int i = 0; i < 16; i++) {
		LPVOID ret = VirtualAlloc((LPVOID)(NV2A_ADDR + NV_PRAMIN_ADDR + i * 64 * ONE_KB), 64 * ONE_KB, MEM_COMMIT, PAGE_READWRITE);
		if (ret != (LPVOID)(NV2A_ADDR + NV_PRAMIN_ADDR + i * 64 * ONE_KB)) {
			CxbxrAbort("VirtualAlloc failed to commit the memory for the nv2a pramin. The error was 0x%08X", GetLastError());
		}
	}

	printf("[0x%.4X] INIT: Allocated %d MiB of Xbox NV2A PRAMIN memory at 0x%.8x to 0x%.8x\n",
		GetCurrentThreadId(), d->pramin.ramin_size / ONE_MB, (uintptr_t)d->pramin.ramin_ptr, (uintptr_t)d->pramin.ramin_ptr + d->pramin.ramin_size - 1);

	// Allocate flat MMIO backing buffer and point struct regs pointers into it
	CxbxAllocateFlatMMIO(d);
}

/* NV2ADevice */

NV2ADevice::NV2ADevice()
{
	m_nv2a_state = new NV2AState();
}

NV2ADevice::~NV2ADevice()
{
	Reset(); // TODO : Review this

	g_renderbase->DeviceRelease();

	delete m_nv2a_state;
}

// PCI Device functions

void NV2ADevice::Init()
{
	PCIBarRegister r;

	// Register Memory bar :
	r.Raw.type = PCI_BAR_TYPE_MEMORY;
	r.Memory.address = NV2A_ADDR >> 4;
	RegisterBAR(0, NV2A_SIZE, r.value);

	// Register physical memory on bar 1
	r.Memory.address = 0xF0000000 >> 4;
	RegisterBAR(1, g_SystemMaxMemory, r.value);
	
	m_DeviceId = 0x02A5;
	m_VendorId = PCI_VENDOR_ID_NVIDIA;
	m_RevisionAndClassCode = 0x030000A1; // VGA-compatible display controller, rev A1

	NV2AState *d = m_nv2a_state; // glue

	CxbxReserveNV2AMemory(d);

	d->pcrtc.start = 0;

	// Enable PMC hardware interrupts - the miniport normally writes this during init,
	// but since our miniport init may race with VBlank delivery, set it here to ensure
	// the ISR can see pending interrupts from the start.
	d->pmc.enabled_interrupts = NV_PMC_INTR_EN_0_HARDWARE;
	d->pcrtc.enabled_interrupts = NV_PCRTC_INTR_0_VBLANK;
	// Enable all PGRAPH interrupt sources - LoadEngines writes 0xFFFFFFFF to
	// NV_PGRAPH_INTR_EN during D3D init, but the MMIO write may arrive after
	// the first pushbuffer commands. Pre-enable so NV097_NO_OPERATION's
	// interrupt handshake works from the first Swap.
	d->pgraph.enabled_interrupts = 0xFFFFFFFF;

	d->vram_ptr = (uint8_t*)PHYSICAL_MAP_BASE;
	d->vram_size = g_SystemMaxMemory;

	d->pramdac.core_clock_coeff = 0x00011C01; /* 233MHz...? */
	d->pramdac.core_clock_freq = 233333324;
	d->pramdac.memory_clock_coeff = 0;
	d->pramdac.video_clock_coeff = 0x0003C20D; /* 25182Khz...? */

	// Setup the conditions/mutexes
	pgraph_init(d);

	d->vblank_last = get_now();
	d->vblank_cb = nv2a_vblank_interrupt;

    qemu_mutex_init(&d->pfifo.pfifo_lock);
    qemu_cond_init(&d->pfifo.puller_cond);
    qemu_cond_init(&d->pfifo.pusher_cond);
    qemu_cond_init(&d->pfifo.flush_complete_cond);
    d->pfifo.flush_requested = false;

    d->pfifo.regs[RI(NV_PFIFO_CACHE1_STATUS)] |= NV_PFIFO_CACHE1_STATUS_LOW_MARK;

    // Initialize RAMHT/RAMFC registers to match what the read handlers return.
    // On real Xbox hardware these are programmed by the BIOS during boot.
    // The D3D runtime reads these to find the RAMHT/RAMFC locations in PRAMIN
    // and writes entries there; ramht_lookup must use the same values internally.
    d->pfifo.regs[RI(NV_PFIFO_RAMHT)] = 0x03000100;
    d->pfifo.regs[RI(NV_PFIFO_RAMFC)] = 0x00890110;

    // Populate RAMHT with handle→instance mappings.
    // On real Xbox the BIOS creates these entries; we emulate that here.
    // RAMHT is at PRAMIN offset 0x10000 (4KB, 512 entries × 8 bytes).
    // Hash for handles < 2048 with channel_id=0 is just the handle value.
    // Entry format: word0=handle, word1=context (valid|engine|instance).
    {
        struct { uint32_t handle; uint32_t pramin_offset; } ramht_entries[] = {
            { 0x02, 0x60 },  // DMA_TO_MEMORY: error notifier (small, base=0x03FD6020)
            { 0x03, 0x10 },  // DMA_FROM_MEMORY: textures/vertices (all VRAM)
            { 0x04, 0x30 },  // Bidirectional: state context (all VRAM)
            { 0x07, 0x70 },  // Bidirectional: extended notifier (small, base=0x03FD6040)
            { 0x08, 0x90 },  // Bidirectional: semaphore/fence (base=0x03FD6000, limit=0x20)
            { 0x09, 0x20 },  // DMA_TO_MEMORY: color render target (all VRAM)
            { 0x0A, 0x40 },  // Bidirectional: zeta/depth buffer (all VRAM)
            { 0x0B, 0xA0 },  // DMA_FROM_MEMORY: second read context (all VRAM)
            { 0x0C, 0x80 },  // Bidirectional: report/occlusion query (256MB)
            { 0x11, 0x10 },  // DMA_FROM_MEMORY: vertex read (shares with handle 3)
            { 0x19, 0x80 },  // Bidirectional: catch-all init DMA (256MB)
        };

        const uint32_t ramht_base = 0x10000; // PRAMIN offset of RAMHT
        for (auto& e : ramht_entries) {
            uint32_t hash = e.handle; // For handles < 2048 with channel_id=0
            uint32_t context = NV_RAMHT_STATUS | NV_RAMHT_ENGINE_GRAPHICS
                             | (e.pramin_offset >> 4);
            uint8_t *entry_ptr = d->pramin.ramin_ptr + ramht_base + hash * 8;
            *(uint32_t*)(entry_ptr + 0) = e.handle;
            *(uint32_t*)(entry_ptr + 4) = context;
        }
    }

    // FIFO threads are started later by StartFifoThreads(), after the host
    // D3D11 device has been created (the puller thread calls D3D11 APIs).
}

void NV2ADevice::StartFifoThreads()
{
	NV2AState *d = m_nv2a_state;
	d->pfifo.puller_thread = std::thread(pfifo_puller_thread, d);
	d->pfifo.pusher_thread = std::thread(pfifo_pusher_thread, d);
}

void NV2ADevice::Reset()
{
	NV2AState *d = m_nv2a_state; // glue
	if (!d) return;

	d->exiting = true;

	qemu_cond_broadcast(&d->pfifo.puller_cond);
	qemu_cond_broadcast(&d->pfifo.pusher_cond);
	qemu_cond_broadcast(&d->pfifo.flush_complete_cond);
	d->pfifo.puller_thread.join();
	d->pfifo.pusher_thread.join();
	qemu_mutex_destroy(&d->pfifo.pfifo_lock); // Cxbxr addition

	pgraph_destroy(&d->pgraph);
}

uint32_t NV2ADevice::IORead(int barIndex, uint32_t port, unsigned size)
{
	return 0;
}

void NV2ADevice::IOWrite(int barIndex, uint32_t port, uint32_t value, unsigned size)
{
}

uint32_t NV2ADevice::BlockRead(const NV2ABlockInfo* block, uint32_t addr, unsigned size)
{
	switch (size) {
	case sizeof(uint8_t) :
		return block->ops.read(m_nv2a_state, addr - block->offset) & 0xFF;
	case sizeof(uint16_t) :
		assert((addr & 1) == 0); // TODO : What if this fails?	

		return block->ops.read(m_nv2a_state, addr - block->offset) & 0xFFFF;
	case sizeof(uint32_t) :
		assert((addr & 3) == 0); // TODO : What if this fails?	

		return block->ops.read(m_nv2a_state, addr - block->offset);
	default:
		assert(false);

		return 0;
	}
}

uint32_t NV2ADevice::MMIORead(int barIndex, uint32_t addr, unsigned size)
{ 
	switch (barIndex) {
	case 0: {
		// Access NV2A regardless of HLE/LLE mode
		const NV2ABlockInfo* block = EmuNV2A_Block(addr);
		if (block != nullptr) {
			return BlockRead(block, addr, size);
		}
		break;
	}
	case 1: {
		// TODO : access physical memory
		break;
	}
	}

	EmuLog(LOG_LEVEL::WARNING, "NV2ADevice::MMIORead: Unhandled barIndex %d, addr %08X, size %d", barIndex, addr, size);
	return 0;
}

void NV2ADevice::BlockWrite(const NV2ABlockInfo* block, uint32_t addr, uint32_t value, unsigned size)
{
	switch (size) {
	case sizeof(uint8_t) : {
#if 0
		xbox::addr_xt aligned_addr;
		uint32_t aligned_value;
		int shift;
		uint32_t mask;

		aligned_addr = addr & ~3;
		aligned_value = block->ops.read(m_nv2a_state, aligned_addr - block->offset);
		shift = (addr & 3) * 8;
		mask = 0xFF << shift;
		block->ops.write(m_nv2a_state, aligned_addr - block->offset, (aligned_value & ~mask) | (value << shift));
#else
		block->ops.write(m_nv2a_state, addr - block->offset, value);
#endif
		return;
	}
	case sizeof(uint16_t) : {
		assert((addr & 1) == 0); // TODO : What if this fails?

		xbox::addr_xt aligned_addr;
		uint32_t aligned_value;
		int shift;
		uint32_t mask;

		aligned_addr = addr & ~3;
		aligned_value = block->ops.read(m_nv2a_state, aligned_addr - block->offset);
		shift = (addr & 2) * 16;
		mask = 0xFFFF << shift;
		block->ops.write(m_nv2a_state, aligned_addr - block->offset, (aligned_value & ~mask) | (value << shift));
		return;
	}
	case sizeof(uint32_t) :
		assert((addr & 3) == 0); // TODO : What if this fails?	

		block->ops.write(m_nv2a_state, addr - block->offset, value);
		return;
	}
}

void NV2ADevice::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	switch (barIndex) {
	case 0: {
		// Access NV2A regardless of HLE/LLE mode
		const NV2ABlockInfo* block = EmuNV2A_Block(addr);

		if (block != nullptr) {
			BlockWrite(block, addr, value, size);
			return;
		}

		break;
	}
	case 1: {
		// TODO : access physical memory
		break;
	}
	}

	EmuLog(LOG_LEVEL::WARNING, "NV2ADevice::MMIOWrite: Unhandled barIndex %d, addr %08X, value %08X, size %d", barIndex, addr, value, size);
}

int NV2ADevice::GetFrameHeight(NV2AState* d)
{
	// Derive frame_height from hardware registers
	int height = ((int)d->prmcio.cr[NV_CIO_CR_VDE_INDEX])
		| (((int)d->prmcio.cr[NV_CIO_CR_OVL_INDEX] & 0x02) >> 1 << 8)
		| (((int)d->prmcio.cr[NV_CIO_CR_OVL_INDEX] & 0x40) >> 6 << 9)
		| (((int)d->prmcio.cr[NV_CIO_CRE_LSR_INDEX] & 0x02) >> 1 << 10);

	return height++;
}

int NV2ADevice::GetFrameWidth(NV2AState* d)
{
	// Derive bytes per pixel from VGA pixel index register
	int bpp;
	switch (d->prmcio.cr[NV_CIO_CRE_PIXEL_INDEX] & 0x03) {
	case 1:  bpp = 1; break;
	case 2:  bpp = 2; break;
	default: bpp = 4; break;
	}

	// Test case : Arctic Thunder, sets a 16 bit framebuffer (R5G6B5) not via
	// AvSetDisplayMode(), but via VGA control register writes, which implies
	// that it's format argument cannot be used to determine the framebuffer
	// width. Instead, read the framebuffer width from the VGA control registers :
	int width = ((int)d->prmcio.cr[NV_CIO_CR_OFFSET_INDEX])
		| (0x700 & ((int)d->prmcio.cr[NV_CIO_CRE_RPC0_INDEX] << 3))
		| (0x800 & ((int)d->prmcio.cr[NV_CIO_CRE_LSR_INDEX] << 6));
	width *= 8;
	width /= bpp;

	return width;
}

uint64_t NV2ADevice::vblank_tick(uint64_t now)
{
	// PCRTC always fires VBlank at the NTSC rate (~59.94Hz / ~16.67ms).
	// Some PAL games (e.g. Dead or Alive Ultimate) disable PCRTC VBlank and
	// instead use PTIMER to generate VBlank interrupts at 50Hz, suggesting
	// that PCRTC can only trigger VBlanks at the NTSC frequency.
	NV2AState *d = m_nv2a_state;
	// ~59.94Hz in QPC ticks: freq * 16667 / 1000000
	const int64_t vblank_period = HostQPCFrequency * 16667 / 1000000;

	uint64_t next = d->vblank_last + vblank_period;

	if (now >= next) {
		// Use the absolute QPC from HostLastQPC (set by get_now() moments
		// before) for PCRTC_RASTER scanline position and jitter profiling.
		int64_t qpcNow = HostLastQPC.load(std::memory_order_relaxed);

		// Measure VBlank jitter: how late (or early) did we fire vs ideal?
		if (g_bCxbxProfilerEnabled) {
			int64_t lastQPC = d->vblank_last_qpc.load(std::memory_order_acquire);
			if (lastQPC > 0) {
				LONGLONG actualTicks = qpcNow - lastQPC;
				LONGLONG jitterTicks = actualTicks > vblank_period
					? actualTicks - vblank_period : vblank_period - actualTicks;
				InterlockedAdd64(&g_ProfileAccum[PROF_VBLANK_JITTER], jitterTicks);
			}
		}

		d->vblank_last_qpc.store(qpcNow, std::memory_order_release);

		d->vblank_cb(d);
		d->vblank_last = now;
		return now + vblank_period;
	}

	return next;
}

uint64_t NV2ADevice::ptimer_tick(uint64_t now)
{
	// Test case: Dead or Alive Ultimate uses this when in PAL50 mode only
	if (m_nv2a_state->ptimer_active) {
		const uint64_t ptimer_period = m_nv2a_state->ptimer_period;
		uint64_t next = m_nv2a_state->ptimer_last + ptimer_period;

		if (now >= next) {
			if (!m_nv2a_state->exiting) [[likely]] {
				m_nv2a_state->ptimer.pending_interrupts |= NV_PTIMER_INTR_0_ALARM;
				update_irq(m_nv2a_state);

				// Wake main thread to dispatch the interrupt (same as VBlank path).
				// Don't call Trigger() here — it would fire the ISR concurrently
				// with DPCs on the main thread, causing a deadlock.
				extern void KeSignalVBlankPending();
				KeSignalVBlankPending();
			}
			m_nv2a_state->ptimer_last = now;
			return now + ptimer_period;
		}

		return next;
	}

	return -1;
}
