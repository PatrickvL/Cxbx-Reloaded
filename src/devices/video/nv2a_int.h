// Source : https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a_int.h
/*
 * QEMU Geforce NV2A internal definitions
 *
 * Copyright (c) 2012 espes
 * Copyright (c) 2015 Jannik Vogel
 * Copyright (c) 2018 Matt Borgerson
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 or
 * (at your option) version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <https://www.gnu.org/licenses/>.
 */

#ifndef HW_NV2A_INT_H
#define HW_NV2A_INT_H

#undef USE_SHADER_CACHE

#include <queue>
#include <thread>
#include <cstdint>
#include <atomic>

#include "xbox_types.h" // For xbox::addr_xt

#include "qemu-thread.h" // For qemu_mutex, etc

#ifdef USE_SHADER_CACHE
#include "glib_compat.h" // For GHashTable, g_hash_table_new, g_hash_table_lookup, g_hash_table_insert
#endif

#include "swizzle.h"

#include "nv2a_debug.h" // For HWADDR_PRIx, NV2A_DPRINTF, NV2A_DPRINTF_IF, etc.
#include "nv2a_regs.h" // For NV2A_MAX_TEXTURES, etc
#include "core\\hle\\D3D8\\Rendering\\NV2A_PGRAPH_Helpers.h" // For NV2ASurfaceState


typedef xbox::addr_xt hwaddr; // Compatibility; Cxbx uses xbox::addr_xt, xqemu and OpenXbox use hwaddr 
typedef uint32_t value_t; // Compatibility; Cxbx values are uint32_t (xqemu and OpenXbox use uint64_t)

// Register index: convert MMIO byte offset to uint32_t array index
#define RI(byte_offset) ((byte_offset) >> 2)

#define NV_PMC_SIZE                 (0x001000 / 4)
#define _NV_PFIFO_SIZE              (0x002000 / 4) // Underscore prefix to prevent clash with NV_PFIFO_SIZE
#define NV_PVIDEO_SIZE              (0x001000 / 4)
#define NV_PTIMER_SIZE              (0x001000 / 4)
#define NV_PFB_SIZE                 (0x001000 / 4)
#define NV_PGRAPH_SIZE              (0x002000 / 4)
#define NV_PCRTC_SIZE               (0x001000 / 4)
#define NV_PRAMDAC_SIZE             (0x001000 / 4)

// Byte sizes for each register block (used for GPU upload, memcpy, etc.)
#define NV_PMC_REGS_BYTES           (NV_PMC_SIZE * sizeof(uint32_t))       // 4 KB
#define NV_PFIFO_REGS_BYTES         (_NV_PFIFO_SIZE * sizeof(uint32_t))    // 8 KB
#define NV_PVIDEO_REGS_BYTES        (NV_PVIDEO_SIZE * sizeof(uint32_t))    // 4 KB
#define NV_PTIMER_REGS_BYTES        (NV_PTIMER_SIZE * sizeof(uint32_t))    // 4 KB
#define NV_PFB_REGS_BYTES           (NV_PFB_SIZE * sizeof(uint32_t))       // 4 KB
#define NV_PGRAPH_REGS_BYTES        (NV_PGRAPH_SIZE * sizeof(uint32_t))    // 8 KB
#define NV_PCRTC_REGS_BYTES         (NV_PCRTC_SIZE * sizeof(uint32_t))     // 4 KB
#define NV_PRAMDAC_REGS_BYTES       (NV_PRAMDAC_SIZE * sizeof(uint32_t))   // 4 KB

// Flat MMIO backing storage: 16 MiB reserved, only engine block pages committed.
// Block offsets mirror real NV2A MMIO layout (offset from NV2A base 0xFD000000).
#define NV2A_MMIO_TOTAL_SIZE        0x01000000u  // 16 MiB
#define NV2A_MMIO_OFF_PMC           0x000000u
#define NV2A_MMIO_OFF_PFIFO         0x002000u
#define NV2A_MMIO_OFF_PVIDEO        0x008000u
#define NV2A_MMIO_OFF_PTIMER        0x009000u
#define NV2A_MMIO_OFF_PFB           0x100000u
#define NV2A_MMIO_OFF_PGRAPH        0x400000u
#define NV2A_MMIO_OFF_PCRTC         0x600000u
#define NV2A_MMIO_OFF_PRAMDAC       0x680000u

// Global flat MMIO buffer pointer (allocated in CxbxReserveNV2AMemory)
extern uint8_t* g_pNV2AMMIO;

#define VSH_TOKEN_SIZE 4 // Compatibility; TODO : Move this to nv2a_vsh.h
#define MAX(a,b) ((a)>(b) ? (a) : (b)) // Compatibility
#define MIN(a,b) ((a)<(b) ? (a) : (b)) // Compatibility

#define g_free(x) free(x) // Compatibility
#define g_malloc(x) malloc(x) // Compatibility
#define g_malloc0(x) calloc(1, x) // Compatibility
#define g_realloc(x, y) realloc(x, y) // Compatibility

#undef USE_TEXTURE_CACHE

#if __cplusplus >= 201402L
#  define NV2A_CONSTEXPR constexpr
#else
#  define NV2A_CONSTEXPR static
#endif

// GCC implementation of FFS
static int ffs(int valu)
{
	int bit;

	if (valu == 0)
		return 0;

	for (bit = 1; !(valu & 1); bit++)
		valu >>= 1;

	return bit;
}

#define GET_MASK(v, mask) (((v) & (mask)) >> (ffs(mask)-1))

#define SET_MASK(v, mask, val)                            \
    do {                                                  \
        const unsigned int __val = (val);                 \
        const unsigned int __mask = (mask);               \
        (v) &= ~(__mask);                                 \
        (v) |= ((__val) << (ffs(__mask) - 1)) & (__mask); \
    } while (0)

// Power-of-two CASE statements
#define CASE_1(v, step) case (v)
#define CASE_2(v, step) CASE_1(v, step) : CASE_1(v + (step) * 1, step)
#define CASE_4(v, step) CASE_2(v, step) : CASE_2(v + (step) * 2, step)
#define CASE_8(v, step) CASE_4(v, step) : CASE_4(v + (step) * 4, step)
#define CASE_16(v, step) CASE_8(v, step) : CASE_8(v + (step) * 8, step)
#define CASE_32(v, step) CASE_16(v, step) : CASE_16(v + (step) * 16, step)
#define CASE_64(v, step) CASE_32(v, step) : CASE_32(v + (step) * 32, step)

// Non-power-of-two CASE statements
#define CASE_3(v, step) CASE_2(v, step) : CASE_1(v + (step) * 2, step)
#define CASE_6(v, step) CASE_4(v, step) : CASE_2(v + (step) * 4, step)

#define NV2A_DEVICE(obj) \
    OBJECT_CHECK(NV2AState, (obj), "nv2a")

//void reg_log_read(int block, hwaddr addr, uint64_t val);
//void reg_log_write(int block, hwaddr addr, uint64_t val);

enum FIFOEngine {
	ENGINE_SOFTWARE = 0,
	ENGINE_GRAPHICS = 1,
	ENGINE_DVD = 2,
};

typedef struct DMAObject {
	unsigned int dma_class;
	unsigned int dma_target;
	xbox::addr_xt address;
	xbox::addr_xt limit;
} DMAObject;

typedef struct VertexAttribute {
	bool dma_select;
	xbox::addr_xt offset;

	/* inline arrays are packed in order?
	* Need to pass the offset to converted attributes */
	unsigned int inline_array_offset;

	float inline_value[4];

	unsigned int format;
	unsigned int size; /* size of the data type */
	unsigned int count; /* number of components */
	uint32_t stride;

	float *inline_buffer;
	float *inline_buffer_pool; // Persistent allocation reused across draws (avoids malloc/free per draw)
} VertexAttribute;

typedef struct Surface {
	bool draw_dirty;
	bool write_enabled_cache;
} Surface;

typedef struct TextureShape {
	bool cubemap;
	unsigned int dimensionality;
	unsigned int color_format;
	unsigned int levels;
	unsigned int width, height, depth;

	unsigned int min_mipmap_level, max_mipmap_level;
	unsigned int pitch;
} TextureShape;

typedef struct TextureKey {
	TextureShape state;
	uint64_t data_hash;
	uint8_t* texture_data;
	uint8_t* palette_data;
} TextureKey;

// Hardware tessellation patch state (NV097_SET_BEGIN_PATCH / SET_END_PATCH)
#define NV2A_PATCH_MAX_CURVES     16   // max attribute curves per patch
#define NV2A_PATCH_MAX_COEFFS   1024   // max float4 coefficient writes per patch

typedef struct PatchCurve {
	int curveType;          // curve type (0=END, 1=STRIP, 2=LEFT_GUARD, 3=RIGHT_GUARD, etc.)
	int coeffStart;         // index into coefficients array
	int coeffCount;         // number of float4 entries
} PatchCurve;

typedef struct PatchState {
	// Configuration from SET_BEGIN_PATCH0-3
	uint32_t patch0;
	uint32_t patch1;
	uint32_t patch2;
	uint32_t patch3;

	// Swatch config from SET_BEGIN_END_SWATCH
	uint32_t swatch;

	// Curve accumulation
	bool active;                    // inside a BEGIN_PATCH..END_PATCH
	int currentCurveAttr;           // current curve attribute (-1 if none)
	int curveCount;                 // number of completed curves
	PatchCurve curves[NV2A_PATCH_MAX_CURVES];

	// Coefficient buffer (shared across all curves)
	float coefficients[NV2A_PATCH_MAX_COEFFS * 4]; // float4 entries
	int totalCoeffs;                // total float4 entries written
} PatchState;

// Dirty group indices for PGRAPHState::dirty[] array.
// Indexed from NV097MethodEntry::dirty_group (0 = no dirty flag).
enum NV2ADirtyGroup {
	NV2A_DIRTY_NONE = 0,         // no dirty flag (must be 0 — table guard skips this)
	NV2A_DIRTY_PGRAPH = 0,       // alias: any pg->regs[] write (set explicitly, not via table guard)
	NV2A_DIRTY_PROGRAM,          // program_data[] was written
	NV2A_DIRTY_SURFACE,          // surface configuration changed
	NV2A_DIRTY_TEXTURE,          // texture state changed
	NV2A_DIRTY_BLEND,            // blend / color mask state changed
	NV2A_DIRTY_RASTERIZER,       // rasterizer state changed (cull, polygon, etc.)
	NV2A_DIRTY_DEPTH_STENCIL,    // depth/stencil state changed
	NV2A_DIRTY_SHADER,           // shader/combiner program changed
	// --- Values above here (1–7) fit in the 3-bit dirty_group table field ---
	NV2A_DIRTY_LIGHTING,         // ltctxa/ltctxb/ltc1/light write (not table-assignable, set in code)
	NV2A_DIRTY_COUNT
};

typedef struct KelvinState {
	xbox::addr_xt object_instance;
} KelvinState;

// NV2A Transform Engine ("Cheops") internal SRAM state.
// These banks are NOT MMIO-mapped; accessed indirectly via pushbuffer
// methods (auto-incrementing CHEOPS_OFFSET) or RDI debug interface.
typedef struct CheopsState {
	// XFPR: Transform Program RAM (136 × 128-bit instructions)
	uint32_t xfpr[NV2A_MAX_TRANSFORM_PROGRAM_LENGTH][VSH_TOKEN_SIZE];

	// XFCTX: Transform Context RAM (192 × float4) — vertex shader constants,
	// matrices, viewport params, eye position, etc.
	uint32_t xfctx[NV2A_VERTEXSHADER_CONSTANTS][4];
	uint32_t xfctx_dirty[6]; // Bitmap: 192 bits across 6 words

	// LTCTXA: Lighting Context A (26 × float4) — fog, ambient, material color,
	// per-light attenuation/spot params
	uint32_t ltctxa[NV2A_LTCTXA_COUNT][4];
	uint32_t ltctxa_dirty[1]; // Bitmap: 26 bits in 1 word

	// LTCTXB: Lighting Context B (52 × float4) — per-light diffuse/specular/ambient colors
	uint32_t ltctxb[NV2A_LTCTXB_COUNT][4];
	uint32_t ltctxb_dirty[2]; // Bitmap: 52 bits across 2 words

	// LTC1: Lighting Constants 1 (20 × float4) — light range, material power params
	uint32_t ltc1[NV2A_LTC1_COUNT][4];
	uint32_t ltc1_dirty[1]; // Bitmap: 20 bits in 1 word

	// SET_TRANSFORM_DATA (0x1E80): input v0 register for LAUNCH_TRANSFORM_PROGRAM
	uint32_t vertex_state_shader_v0[4];
} CheopsState;

typedef struct ContextSurfaces2DState {
	xbox::addr_xt object_instance;
	xbox::addr_xt dma_image_source;
	xbox::addr_xt dma_image_dest;
	unsigned int color_format;
	unsigned int source_pitch, dest_pitch;
	xbox::addr_xt source_offset, dest_offset;
} ContextSurfaces2DState;

typedef struct ImageBlitState {
	xbox::addr_xt object_instance;
	xbox::addr_xt context_surfaces;
	unsigned int operation;
	unsigned int in_x, in_y;
	unsigned int out_x, out_y;
	unsigned int width, height;
} ImageBlitState;

typedef struct PGRAPHState {
	QemuMutex pgraph_lock;

	uint32_t pending_interrupts;
	uint32_t enabled_interrupts;
	QemuCond interrupt_cond;

	/* subchannels state we're not sure the location of... */
	ContextSurfaces2DState context_surfaces_2d;
	ImageBlitState image_blit;
	KelvinState kelvin;

	QemuCond fifo_access_cond;
	QemuCond flip_3d;

	Surface surface_color, surface_zeta;
	NV2ASurfaceState surface_state;

	xbox::addr_xt dma_semaphore;

	xbox::addr_xt dma_a, dma_b; // PRAMIN offsets for DMA context A/B (texture/palette address resolution)

	xbox::addr_xt dma_report;
	unsigned int zpass_pixel_count_enable;
	bool zpass_pixel_count_active; // true when a D3D11 occlusion query is between Begin/End
	unsigned int zpass_pixel_count_result;

	unsigned int primitive_mode;

	uint32_t clear_surface_flags; // NV097_CLEAR_SURFACE parameter (Z/STENCIL/COLOR mask)

	// NV2A Transform Engine ("Cheops") — internal SRAM banks
	CheopsState xf;

	// Dirty generation counters indexed by NV2ADirtyGroup.  Bumped by
	// nv097_dispatch_method and PGRAPH switch handlers; each consumer
	// independently tracks its own "last seen" value per group.
	uint32_t dirty[NV2A_DIRTY_COUNT];

	// Light geometry — SRAM bank unknown (xemu: "should figure out where
	// these are in lighting context").  Packed contiguously for data-driven
	// dispatch via NV097_TARGET_LIGHT.  Per-light: 12 floats (48 bytes).
	struct LightGeometry {
		float infinite_half_vector[3];
		float infinite_direction[3];
		float local_position[3];
		float local_attenuation[3];
	} light[NV2A_MAX_LIGHTS];
	static_assert(sizeof(LightGeometry) == 12 * sizeof(float), "LightGeometry must be 48 bytes");

	float point_params[8]; // NV097_SET_POINT_PARAMS attenuation coefficients
	float line_width;      // NV097_SET_LINE_WIDTH (float, pixels)

	VertexAttribute vertex_attributes[NV2A_VERTEXSHADER_ATTRIBUTES];
	uint32_t vertex_attributes_generation; // bumped when FORMAT or OFFSET changes

	unsigned int inline_array_length;
	uint32_t inline_array[NV2A_MAX_BATCH_LENGTH];
	unsigned int inline_elements_length;
	uint16_t inline_elements[NV2A_MAX_BATCH_LENGTH]; // Cxbx-Reloaded TODO : Restore uint32_t once D3D11_draw_inline_elements can using that

	unsigned int inline_buffer_length;

	unsigned int draw_arrays_length;
	unsigned int draw_arrays_max_count;
	bool draw_arrays_prevent_connect;  // Don't merge adjacent entries across bracket boundaries

	int32_t draw_arrays_start[1250];
	int32_t draw_arrays_count[1250];

	// Hardware tessellation state
	PatchState patch;

	bool texture_matrix_enable[NV2A_MAX_TEXTURES]; // NV097_SET_TEXTURE_MATRIX_ENABLE per stage

	uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PGRAPH
} PGRAPHState;

typedef struct OverlayState {
	bool video_buffer_use;
	int pitch;
	bool is_transparent;
#ifdef DEBUG
	hwaddr base;
	hwaddr limit;
#endif
	hwaddr offset;
	uint32_t in_height;
	uint32_t in_width;
	int out_x;
	int out_y;
	int out_width;
	int out_height;

	bool covers_framebuffer;
	int old_in_width;
	int old_in_height;
	int old_pitch;
} OverlayState;

typedef struct NV2AState {
	void(* vblank_cb)(void *);
	uint64_t vblank_last;
	std::atomic<int64_t> vblank_last_qpc{0}; // QPC timestamp of last VBlank (for PCRTC_RASTER sync)
	std::atomic_flag vblank_pending = ATOMIC_FLAG_INIT; // Set by timer, consumed by main thread
    // PCIDevice dev;
    // qemu_irq irq;
    bool exiting;
	bool enable_overlay = false;
	bool ptimer_active = false;
	uint64_t ptimer_last;
	uint64_t ptimer_period;

    // VGACommonState vga;
    // GraphicHwOps hw_ops;
    // QEMUTimer *vblank_timer;

    // MemoryRegion *vram;
    // MemoryRegion vram_pci;
    uint8_t *vram_ptr;
    size_t vram_size;
    // MemoryRegion ramin;
	struct {
		uint8_t *ramin_ptr;
		size_t ramin_size;
	} pramin;

    // MemoryRegion mmio;
    // MemoryRegion block_mmio[NV_NUM_BLOCKS];

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PMC
    } pmc;

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PFIFO
		QemuMutex pfifo_lock;
		std::thread puller_thread;
		QemuCond puller_cond;
		std::thread pusher_thread;
		QemuCond pusher_cond;
		// Flush synchronization: HLE thread signals flush_requested, then
		// waits on flush_complete_cond.  The puller signals back when CACHE1
		// is drained (LOW_MARK set) and the pusher has no pending DMA data.
		bool flush_requested;
		QemuCond flush_complete_cond;
    } pfifo;

    struct {
		uint32_t pending_interrupts;
		uint32_t enabled_interrupts;
		//QemuCond interrupt_cond; // pvideo.interrupt_cond not used (yet)
		OverlayState overlays[2]; // NV2A supports 2 video overlays
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PVIDEO
    } pvideo;

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
        uint32_t numerator;
        uint32_t denominator;
        uint32_t alarm_time;
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PTIMER
    } ptimer;

    struct {
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PFB
    } pfb;

    struct PGRAPHState pgraph;

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
        hwaddr start;
        uint32_t vblank_count; // Incremented each VBlank; bit 0 determines interlace field (even/odd)
        uint32_t last_present_vblank; // VBlank count at last present (prevents double-present)
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PCRTC
    } pcrtc;

    struct {
        uint32_t core_clock_coeff;
        uint64_t core_clock_freq;
        uint32_t memory_clock_coeff;
        uint32_t video_clock_coeff;
		uint32_t* regs; // Backed by g_pNV2AMMIO + NV2A_MMIO_OFF_PRAMDAC
    } pramdac;

	// PRMDIO: VGA DAC palette (gamma LUT)
	struct {
		uint32_t write_mode_address; // byte index into palette[] (0..767)
		uint32_t read_mode_address;  // byte index for read-back
		uint8_t  palette[256 * 3];   // 256 entries x RGB (written sequentially)
		bool     dirty;              // set true when palette[] is modified
	} puserdac;

	// PRMCIO (Actually the VGA controller)
	struct {
		uint8_t cr_index;
		uint8_t cr[256]; /* CRT registers */
		uint8_t ar_index;
		uint8_t ar[0x15]; /* Attribute Controller registers (VGA_ATT_C) */
		bool    ar_flip_flop;  /* false=index, true=data */
	} prmcio; // Not in xqemu/openxbox?

	// PRMVIO: VGA Sequencer and Graphics Controller
	struct {
		uint8_t seq_index;
		uint8_t seq[256];   /* Sequencer registers (VGA_SEQ_C used) */
		uint8_t gfx_index;
		uint8_t gfx[256];   /* Graphics Controller registers (VGA_GFX_C used) */
		uint8_t misc_output; /* Misc Output Register */
	} prmvio;
} NV2AState;

typedef value_t(*read_func)(NV2AState *d, hwaddr addr); //, unsigned int size);
typedef void(*write_func)(NV2AState *d, hwaddr addr, value_t val); //, unsigned int size);

typedef struct {
	read_func read;
	write_func write;
} MemoryRegionOps;

#if 0
// Valid after PCI init :
#define NV20_REG_BASE_KERNEL 0xFD000000

typedef volatile DWORD *PPUSH;

typedef struct {
	DWORD Ignored[0x10];
	PPUSH Put; // On Xbox1, this field is only written to by the CPU (the GPU uses this as a trigger to start executing from the given address)
	PPUSH Get; // On Xbox1, this field is only read from by the CPU (the GPU reflects in here where it is/stopped executing)
	PPUSH Reference; // TODO : xbox::addr_xt / void* / DWORD ? 
	DWORD Ignored2[0x7ED];
} Nv2AControlDma;
#endif

#endif
