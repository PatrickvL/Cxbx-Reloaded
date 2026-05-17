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
// *  This file is heavily based on code from XQEMU
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a_pfifo.c
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

typedef struct RAMHTEntry {
	uint32_t handle;
	xbox::addr_xt instance;
	enum FIFOEngine engine;
	unsigned int channel_id : 5;
	bool valid;
} RAMHTEntry;

#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11_Profiler.h"

static RAMHTEntry ramht_lookup(NV2AState *d, uint32_t handle); // forward declaration
static void pfifo_run_puller(NV2AState *d); // forward declaration
static void pfifo_run_pusher(NV2AState *d); // forward declaration

/* PFIFO - MMIO and DMA FIFO submission to PGRAPH and VPE */
DEVICE_READ32(PFIFO)
{
    // Fast path for DMA_GET reads.  When the DMA pusher cannot process
    // (access flags not set in HLE mode), advance GET to PUT so polls
    // like BlockUntilIdle() return immediately.  When the pusher CAN
    // process, drain pending commands inline so the native polling loop
    // sees GET advance and exits naturally.
    if (addr == NV_PFIFO_CACHE1_DMA_GET) {
        uint32_t get_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
        uint32_t put_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];
        if (get_v != put_v) {
            uint32_t push0    = d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH0)];
            uint32_t dma_push = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUSH)];
            bool pusher_can_run = GET_MASK(push0, NV_PFIFO_CACHE1_PUSH0_ACCESS)
                               && GET_MASK(dma_push, NV_PFIFO_CACHE1_DMA_PUSH_ACCESS)
                               && !GET_MASK(dma_push, NV_PFIFO_CACHE1_DMA_PUSH_STATUS);
            if (!pusher_can_run) {
                d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)] = put_v;
                get_v = put_v;
            } else {
                // Drain pending commands inline — enables native polling loops
                // (BlockUntilIdle) to work without patches.
                pfifo_flush_to_pgraph(d);
                get_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
            }
        }
        uint32_t result = get_v;
        DEVICE_READ32_END(PFIFO);
    }

    qemu_mutex_lock(&d->pfifo.pfifo_lock);

	DEVICE_READ32_SWITCH() {
	case NV_PFIFO_RAMHT:
		result = 0x03000100; // = NV_PFIFO_RAMHT_SIZE_4K | NV_PFIFO_RAMHT_BASE_ADDRESS(NumberOfPaddingBytes >> 12) | NV_PFIFO_RAMHT_SEARCH_128
		break;
	case NV_PFIFO_RAMFC:
		result = 0x00890110; // = ? | NV_PFIFO_RAMFC_SIZE_2K | ?
		break;
	case NV_PFIFO_INTR_0:
		result = d->pfifo.pending_interrupts;
		break;
	case NV_PFIFO_INTR_EN_0:
		result = d->pfifo.enabled_interrupts;
		break;
	case NV_PFIFO_RUNOUT_STATUS:
		result = NV_PFIFO_RUNOUT_STATUS_LOW_MARK; /* low mark empty */
		break;
	default:
		DEVICE_READ32_REG(pfifo); // Was : DEBUG_READ32_UNHANDLED(PFIFO);
		break;
	}

    qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	DEVICE_READ32_END(PFIFO);
}

DEVICE_WRITE32(PFIFO)
{
    qemu_mutex_lock(&d->pfifo.pfifo_lock);

	switch(addr) {
		case NV_PFIFO_INTR_0:
			d->pfifo.pending_interrupts &= ~value;
			update_irq(d);
			break;
		case NV_PFIFO_INTR_EN_0:
			d->pfifo.enabled_interrupts = value;
			update_irq(d);
			break;
		default:
			DEVICE_WRITE32_REG(pfifo); // Was : DEBUG_WRITE32_UNHANDLED(PFIFO);
			break;
	}

    qemu_cond_broadcast(&d->pfifo.pusher_cond);
    qemu_cond_broadcast(&d->pfifo.puller_cond);

    qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	DEVICE_WRITE32_END(PFIFO);
}

static void pfifo_run_puller(NV2AState *d)
{
    uint32_t *pull0 = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_PULL0)];
    uint32_t *pull1 = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_PULL1)];
    uint32_t *engine_reg = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_ENGINE)];

    uint32_t *status = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_STATUS)];
    uint32_t *get_reg = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_GET)];
    uint32_t *put_reg = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUT)];

    // Acquire pgraph_lock once for the entire CACHE1 drain rather than per
    // method.  Eliminates N-1 redundant lock/unlock pairs per puller wake.
    // pgraph_handle_method may internally release/reacquire for waits (e.g.
    // NV097_NO_OPERATION interrupt handshake), which is safe because
    // QemuMutex wraps CRITICAL_SECTION (see thread-win32.h) — reentrant
    // by definition on Windows.
    qemu_mutex_lock(&d->pgraph.pgraph_lock);

    while (true) {
        if (!GET_MASK(*pull0, NV_PFIFO_CACHE1_PULL0_ACCESS)) break;

        /* empty cache1 */
        if (*status & NV_PFIFO_CACHE1_STATUS_LOW_MARK) break;

        uint32_t get = *get_reg;
        uint32_t put = *put_reg;

        assert(get < 128*4 && (get % 4) == 0);
        uint32_t method_entry = d->pfifo.regs[RI(NV_PFIFO_CACHE1_METHOD + get*2)];
        uint32_t parameter = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DATA + get*2)];

        uint32_t new_get = (get+4) & 0x1fc;
        *get_reg = new_get;

        if (new_get == put) {
            // set low mark
            *status |= NV_PFIFO_CACHE1_STATUS_LOW_MARK;
        }
        if (*status & NV_PFIFO_CACHE1_STATUS_HIGH_MARK) {
            // unset high mark
            *status &= ~NV_PFIFO_CACHE1_STATUS_HIGH_MARK;
            // signal pusher
            qemu_cond_signal(&d->pfifo.pusher_cond);            
        }

        uint32_t method = method_entry & 0x1FFC;
        uint32_t subchannel = GET_MASK(method_entry, NV_PFIFO_CACHE1_METHOD_SUBCHANNEL);

        // Process pushbuffer methods into PGRAPH register state.
        // Skip object binding (method 0) — Xbox uses a single channel.
        if (method >= 0x180 && method < 0x200) {
            // DMA context binding methods: parameter is a handle that must
            // be resolved via RAMHT to get the PRAMIN instance address.
            RAMHTEntry entry = ramht_lookup(d, parameter);
            if (entry.valid) {
                pgraph_handle_method(d, subchannel, method, entry.instance);
            }
        } else if (method >= 0x100) {
            pgraph_handle_method(d, subchannel, method, parameter);
        }

    }

    qemu_mutex_unlock(&d->pgraph.pgraph_lock);
}

// Defined in HostSync.cpp — marks the current thread as the PFIFO puller
// so CxbxUpdateNativeD3DResources skips pfifo_flush (prevents deadlock).
extern void CxbxSetPullerContext(bool active);

// Forward declaration: auto-present on VBlank for games without explicit FLIP_STALL
#include "nv2a_pgraph_backend.h"
extern bool g_pgraph_explicit_flip_stall_seen;

int pfifo_puller_thread(NV2AState *d)
{
    g_AffinityPolicy->SetAffinityOther();
    CxbxSetThreadName("Cxbx NV2A FIFO puller");
    CxbxSetPullerContext(true);

    qemu_mutex_lock(&d->pfifo.pfifo_lock);
    while (!d->exiting) {
        pfifo_run_puller(d);

        // Present logic — at most one present per wake-up cycle.
        // The two paths are mutually exclusive (else-if) to prevent
        // double-presenting when both auto-present and overlay are relevant.
        // Release pfifo_lock during flip_stall to avoid deadlock: the HLE
        // draw path acquires D3D11 lock → pfifo_lock (via pfifo_flush),
        // while the puller holds pfifo_lock → D3D11 lock (via flip_stall).
        if (g_pgraph_backend.flip_stall) {
            if (!g_pgraph_explicit_flip_stall_seen
                && d->pgraph.surface_color.draw_dirty) {
                // Auto-present fallback for raw pushbuffer games that never
                // issue NV097_FLIP_STALL.
                d->pgraph.surface_color.draw_dirty = false;
                qemu_mutex_unlock(&d->pfifo.pfifo_lock);
                g_pgraph_backend.flip_stall(d);
                qemu_mutex_lock(&d->pfifo.pfifo_lock);
            } else if (d->enable_overlay
                       && !d->pgraph.surface_color.draw_dirty) {
                // PVIDEO overlay present: during FMV playback, no PGRAPH
                // draws or FLIP_STALL commands occur. The VBlank handler
                // wakes us when the overlay is active so video frames are
                // composited and displayed at frame rate.
                qemu_mutex_unlock(&d->pfifo.pfifo_lock);
                g_pgraph_backend.flip_stall(d);
                qemu_mutex_lock(&d->pfifo.pfifo_lock);
            }
        }

        // If the HLE thread is waiting for a PFIFO flush, signal it now
        // that CACHE1 has been drained.
        if (d->pfifo.flush_requested) {
            qemu_cond_signal(&d->pfifo.flush_complete_cond);
        }

        qemu_cond_wait(&d->pfifo.puller_cond, &d->pfifo.pfifo_lock);
    }
    qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	return 0;
}

// ---------------------------------------------------------------------------
// pfifo_submit_pushbuffer  --  submit a block of NV2A push buffer commands
// through a software command parser that dispatches methods directly to
// PGRAPH.  This replaces the old EmuExecutePushBufferRaw() which did the
// same thing but from outside the NV2A module.
//
// The push buffer contains NV2A FIFO commands (method headers + data words)
// as written by the Xbox D3D runtime during BeginPushBuffer recording.
// We parse them using the standard NV4 DMA pusher format and feed each
// method to pgraph_handle_method, just as the PFIFO puller would.
// ---------------------------------------------------------------------------

void pfifo_submit_pushbuffer(NV2AState *d, void *pPushData, uint32_t uSizeInBytes)
{
    if (!pPushData || uSizeInBytes < 4) return;

    // Mark this thread as a puller context so that draw callbacks
    // (HLE_draw_state_update → CxbxUpdateNativeD3DResources) skip
    // pfifo_flush_to_pgraph (registers are already current) and
    // pgraph_lock acquisition (we hold it for the whole buffer below).
    CxbxSetPullerContext(true);

    // Ensure PGRAPH context control has channel ID set (matches EmuExecutePushBufferRaw)
    d->pgraph.regs[RI(NV_PGRAPH_CTX_CONTROL)] |= NV_PGRAPH_CTX_CONTROL_CHID;

    uint32_t *dma_limit = (uint32_t*)((uintptr_t)pPushData + uSizeInBytes);
    uint32_t *dma_put   = dma_limit;
    uint32_t *dma_get   = (uint32_t*)pPushData;

    // DMA pusher state
    struct {
        uint32_t mthd;    // Current method
        uint32_t subc;    // Current subchannel
        uint32_t mcnt;    // Remaining method count
        bool     ni;      // Non-increasing flag
    } state = {};

    bool subr_active = false;
    uint32_t *subr_return = nullptr;

    // Acquire pgraph_lock once for the entire pushbuffer rather than once per
    // dispatched method.  pgraph_handle_method may release and re-acquire it
    // internally (CRITICAL_SECTION is reentrant), but the net effect is that
    // we hold it across all methods, eliminating N-1 redundant lock/unlock pairs.
    qemu_mutex_lock(&d->pgraph.pgraph_lock);

    while (dma_get != dma_put) {
        if (dma_get >= dma_limit) {
            EmuLog(LOG_LEVEL::WARNING, "pfifo_submit_pushbuffer: overran buffer");
            goto done;
        }

        uint32_t word = *dma_get++;

        // Data word of an active method command
        if (state.mcnt) {
            // Dispatch subchannel 0 methods (3D) to PGRAPH,
            // matching the old EmuExecutePushBufferRaw behaviour.
            // pgraph_lock is held for the whole buffer (acquired above).
            if (state.subc == 0) {
                uint32_t method = state.mthd << 2;
                if (method >= 0x180 && method < 0x200) {
                    // DMA context binding: resolve handle via RAMHT
                    RAMHTEntry entry = ramht_lookup(d, word);
                    if (entry.valid) {
                        pgraph_handle_method(d, state.subc, method, entry.instance);
                    }
                } else {
                    pgraph_handle_method(d, state.subc, method, word);
                }

                // Squash repeated BEGIN/DRAW_ARRAYS/END (xemu approach):
                // peek ahead for END, BEGIN(same_mode), DRAW_ARRAYS pattern.
                if (method == NV097_DRAW_ARRAYS && state.mcnt == 1) {
                    PGRAPHState *pg = &d->pgraph;
                    ptrdiff_t remaining_words = dma_put - dma_get;
                    if (remaining_words >= 6 &&
                        pg->inline_elements_length == 0 &&
                        pg->draw_arrays_length > 0 &&
                        pg->draw_arrays_length < (ARRAY_SIZE(pg->draw_arrays_start) - 1)) {
                        uint32_t w0 = dma_get[0];  // expected: END header
                        uint32_t w1 = dma_get[1];  // expected: END param (0)
                        uint32_t w2 = dma_get[2];  // expected: BEGIN header
                        uint32_t w3 = dma_get[3];  // expected: BEGIN param (primitive_mode)
                        uint32_t w4 = dma_get[4];  // expected: DRAW_ARRAYS header
                        if ((w0 & 0x1FFC) == NV097_SET_BEGIN_END &&
                            w1 == NV097_SET_BEGIN_END_OP_END &&
                            (w2 & 0x1FFC) == NV097_SET_BEGIN_END &&
                            w3 == pg->primitive_mode &&
                            (w4 & 0x1FFC) == NV097_DRAW_ARRAYS) {
                            dma_get += 4;  // skip END hdr+param, BEGIN hdr+param
                            pg->draw_arrays_prevent_connect = true;
                        }
                    }
                }
            }

            if (!state.ni)
                state.mthd++;

            state.mcnt--;
            continue;
        }

        // First word of a new command — decode type/instruction/flags
        // Using the NV4 DMA pusher command format.
        uint32_t type = word & 3;
        if (type == 1) {
            // JUMP_LONG
            dma_get = (uint32_t*)(CONTIGUOUS_MEMORY_BASE | (word & 0xFFFFFFFC));
            continue;
        }
        if (type == 2) {
            // CALL
            subr_return = dma_get;
            subr_active = true;
            dma_get = (uint32_t*)(CONTIGUOUS_MEMORY_BASE | (word & 0xFFFFFFFC));
            continue;
        }

        // type == 0: check instruction field (bits 29-31)
        uint32_t instruction = (word >> 29) & 7;
        if (instruction == 1) {
            // JUMP (short form)
            dma_get = (uint32_t*)(CONTIGUOUS_MEMORY_BASE | (word & 0x1FFFFFFC));
            continue;
        }

        // Check flags (bits 16-17)
        uint32_t flags = (word >> 16) & 3;
        if (flags == 2) {
            // RETURN
            if (subr_active) {
                dma_get = subr_return;
                subr_active = false;
            }
            continue;
        }

        // Increasing or non-increasing methods
        state.ni   = (instruction == 2);
        state.mthd = (word >> 2) & 0x7FF;          // method / 4
        state.subc = (word >> 13) & 7;
        state.mcnt = (word >> 18) & 0x7FF;
    }

done:
    qemu_mutex_unlock(&d->pgraph.pgraph_lock);
    CxbxSetPullerContext(false);
}

// ---------------------------------------------------------------------------
// pfifo_flush_to_pgraph  --  drain all pending DMA pushbuffer commands into
// PGRAPH so register state is current before the next draw.  Called from the
// HLE thread (CxbxUpdateNativeD3DResources) before each draw.
//
// The pushbuffer is processed *inline* on the calling thread instead of
// waking the background pusher thread and sleeping until it finishes.
// Inline processing eliminates the two OS thread context switches
// (cond_signal to wake pusher + cond_wait to sleep until done) that the
// old round-trip protocol required — each switch costs ~10-50 µs on Windows,
// so the saving is up to ~100 µs per draw call.
//
// pfifo_lock is already held for the duration so the pusher thread cannot
// be concurrently modifying PFIFO/DMA state.  After the inline run the
// pusher thread will simply find an empty buffer on its next wake.
// ---------------------------------------------------------------------------
void pfifo_flush_to_pgraph(NV2AState *d)
{
    // Fast path: lockless check — if DMA buffer is already drained, skip
    // the mutex entirely.  This is the common case for back-to-back draws
    // within the same state batch (all methods already processed by a
    // prior flush or by the puller thread).  The reads are benign races:
    // if GET/PUT change between read and lock, the locked re-read will
    // catch it.  CACHE1 (puller queue) is also checked: if empty, there
    // are no pending methods anywhere.
    {
        uint32_t get_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
        uint32_t put_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];
        uint32_t cache_status = d->pfifo.regs[RI(NV_PFIFO_CACHE1_STATUS)];
        if (get_v == put_v && (cache_status & NV_PFIFO_CACHE1_STATUS_LOW_MARK)) {
            return; // Nothing pending — no lock needed
        }
    }

    qemu_mutex_lock(&d->pfifo.pfifo_lock);

    uint32_t get_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
    uint32_t put_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];

    if (get_v != put_v) {
        // Check whether the DMA pusher can actually process commands.
        // If not (access flags not set), the commands were submitted by
        // HLE-patched D3D calls that already set PGRAPH state — skip them.
        uint32_t push0    = d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH0)];
        uint32_t dma_push = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUSH)];
        bool pusher_can_run = GET_MASK(push0, NV_PFIFO_CACHE1_PUSH0_ACCESS)
                           && GET_MASK(dma_push, NV_PFIFO_CACHE1_DMA_PUSH_ACCESS)
                           && !GET_MASK(dma_push, NV_PFIFO_CACHE1_DMA_PUSH_STATUS);
        if (pusher_can_run) {
            // Process the DMA push buffer inline on the calling thread.
            // pfifo_lock is already held, exactly as in pfifo_pusher_thread.
            // Mark this thread as puller context so that any draw callback
            // triggered by pgraph_handle_method (e.g. pgraph_draw_arrays)
            // does not attempt a re-entrant pfifo_flush_to_pgraph.
            CxbxSetPullerContext(true);
            pfifo_run_pusher(d);
            CxbxSetPullerContext(false);
        } else {
            // Advance GET past the unprocessable commands.
            d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)] = put_v;
        }
    }

    d->pfifo.flush_requested = false;
    qemu_mutex_unlock(&d->pfifo.pfifo_lock);
}

static void pfifo_run_pusher(NV2AState *d)
{
    uint32_t *push0 = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH0)];
    uint32_t *push1 = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH1)];
    uint32_t *dma_subroutine = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_SUBROUTINE)];
    uint32_t *dma_state = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_STATE)];
    uint32_t *dma_push = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUSH)];
    uint32_t *dma_get = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
    uint32_t *dma_put = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];
    uint32_t *dma_dcount = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_DCOUNT)];

    if (!GET_MASK(*push0, NV_PFIFO_CACHE1_PUSH0_ACCESS)) return;
    if (!GET_MASK(*dma_push, NV_PFIFO_CACHE1_DMA_PUSH_ACCESS)) return;

    /* suspended */
    if (GET_MASK(*dma_push, NV_PFIFO_CACHE1_DMA_PUSH_STATUS)) return;

    // TODO: should we become busy here??
    // NV_PFIFO_CACHE1_DMA_PUSH_STATE _BUSY

    unsigned int channel_id = GET_MASK(*push1,
                                       NV_PFIFO_CACHE1_PUSH1_CHID);


	/* Channel running DMA */
	uint32_t channel_modes = d->pfifo.regs[RI(NV_PFIFO_MODE)];
	assert(channel_modes & (1 << channel_id));

    assert(GET_MASK(*push1, NV_PFIFO_CACHE1_PUSH1_MODE)
            == NV_PFIFO_CACHE1_PUSH1_MODE_DMA);

	/* We're running so there should be no pending errors... */
    assert(GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR)
            == NV_PFIFO_CACHE1_DMA_STATE_ERROR_NONE);

    hwaddr dma_instance =
        GET_MASK(d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_INSTANCE)],
                 NV_PFIFO_CACHE1_DMA_INSTANCE_ADDRESS_MASK) << 4; // TODO : Use NV_PFIFO_CACHE1_DMA_INSTANCE_ADDRESS_MOVE?

    hwaddr dma_len;
    uint8_t *dma = (uint8_t*)nv_dma_map(d, dma_instance, &dma_len);

    // Acquire pgraph_lock once for the entire pushbuffer rather than once per
    // method.  A typical frame has hundreds of methods; eliminating the
    // per-method lock/unlock pair saves ~100 ns × N ≈ tens of µs per frame.
    // pgraph_handle_method may release and re-acquire pgraph_lock internally
    // (e.g. NV097_NO_OPERATION notification, context switch) — CRITICAL_SECTION
    // is reentrant, so this is safe: each internal unlock/relock is balanced and
    // the function returns with the lock held.
    qemu_mutex_lock(&d->pgraph.pgraph_lock);

	/* based on the convenient pseudocode in envytools */
    while (true) {
        uint32_t dma_get_v = *dma_get;
        uint32_t dma_put_v = *dma_put;
        if (dma_get_v == dma_put_v) break;
        if (dma_get_v >= dma_len) {
            assert(false);
            SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR,
                     NV_PFIFO_CACHE1_DMA_STATE_ERROR_PROTECTION);
            break;
        }

        uint32_t word = ldl_le_p((uint32_t*)(dma + dma_get_v));
        dma_get_v += 4;

        uint32_t method_type =
            GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_TYPE);
        uint32_t method_subchannel =
            GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_SUBCHANNEL);
        uint32_t method =
            GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD) << 2;
        uint32_t method_count =
            GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_COUNT);

        uint32_t subroutine_state =
            GET_MASK(*dma_subroutine, NV_PFIFO_CACHE1_DMA_SUBROUTINE_STATE);

        if (method_count) {
            /* data word of methods command */
            d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_DATA_SHADOW)] = word;

            // Handle PFIFO-level class methods that don't go to PGRAPH:
            // REF_CNT (method 0x0050, per envytools hw/fifo/puller.html):
            // updates NV_PFIFO_CACHE1_REF, readable via NV_USER_REF (xemu
            // user.c offset 0x48).
            if (method == 0x0050) {
                d->pfifo.regs[RI(NV_PFIFO_CACHE1_REF)] = word;
            }
            // Bypass CACHE1: dispatch directly to PGRAPH instead of staging in
            // the CACHE1 ring buffer for the puller thread to pick up later.
            // This eliminates a full OS thread wake cycle per command batch.
            // pgraph_lock is held for the whole buffer (acquired above).
            else if (method >= 0x180 && method < 0x200) {
                // DMA context binding methods: parameter is a handle that must
                // be resolved via RAMHT to get the PRAMIN instance address.
                RAMHTEntry entry = ramht_lookup(d, word);
                if (entry.valid) {
                    pgraph_handle_method(d, method_subchannel, method, entry.instance);
                }
            } else if (method >= 0x100) {
                pgraph_handle_method(d, method_subchannel, method, word);
            }

            // Squash repeated BEGIN/DRAW_ARRAYS/END sequences (xemu approach):
            // After the last data word of a DRAW_ARRAYS command, peek ahead in
            // the DMA pushbuffer for the exact pattern END, BEGIN(same_mode),
            // DRAW_ARRAYS.  If found, skip the END/BEGIN pair so the next
            // DRAW_ARRAYS accumulates into the same batch.  Any intervening
            // method (texture switch, render target change, etc.) breaks the
            // pattern and prevents incorrect cross-state merging.
            if (method == NV097_DRAW_ARRAYS && method_count == 1) {
                PGRAPHState *pg = &d->pgraph;
                uint32_t dma_put_v = *dma_put;
                uint32_t remaining = (dma_put_v > dma_get_v) ? (dma_put_v - dma_get_v) : 0;
                if (remaining >= 24 &&  // 6 words: END hdr+param, BEGIN hdr+param, DA hdr+param
                    pg->inline_elements_length == 0 &&
                    pg->draw_arrays_length > 0 &&
                    pg->draw_arrays_length < (ARRAY_SIZE(pg->draw_arrays_start) - 1)) {
                    uint32_t *peek = (uint32_t*)(dma + dma_get_v);
                    uint32_t w0 = ldl_le_p(&peek[0]);  // expected: END header
                    uint32_t w1 = ldl_le_p(&peek[1]);  // expected: END param (0)
                    uint32_t w2 = ldl_le_p(&peek[2]);  // expected: BEGIN header
                    uint32_t w3 = ldl_le_p(&peek[3]);  // expected: BEGIN param (primitive_mode)
                    uint32_t w4 = ldl_le_p(&peek[4]);  // expected: DRAW_ARRAYS header
                    if ((w0 & 0x1FFC) == NV097_SET_BEGIN_END &&
                        w1 == NV097_SET_BEGIN_END_OP_END &&
                        (w2 & 0x1FFC) == NV097_SET_BEGIN_END &&
                        w3 == pg->primitive_mode &&
                        (w4 & 0x1FFC) == NV097_DRAW_ARRAYS) {
                        // Skip END header+param and BEGIN header+param (4 words)
                        dma_get_v += 16;
                        pg->draw_arrays_prevent_connect = true;
                    }
                }
            }

            if (method_type == NV_PFIFO_CACHE1_DMA_STATE_METHOD_TYPE_INC) {
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD,
                         (method + 4) >> 2);
            }
            SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_COUNT,
                     method_count - 1);
            (*dma_dcount)++;
		} else {
			/* no command active - this is the first word of a new one */
            d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_RSVD_SHADOW)] = word;

			/* match all forms */
			if ((word & 0xe0000003) == 0x20000000) {
				/* old jump */
                d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET_JMP_SHADOW)] =
                    dma_get_v;
                dma_get_v = word & 0x1fffffff;
				NV2A_DPRINTF("pb OLD_JMP 0x%08X\n", dma_get_v);
			} else if ((word & 3) == 1) {
				/* jump */
                d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET_JMP_SHADOW)] =
                    dma_get_v;
                dma_get_v = word & 0xfffffffc;
				NV2A_DPRINTF("pb JMP 0x%08X\n", dma_get_v);
			} else if ((word & 3) == 2) {
				/* call */
                if (subroutine_state) {
                    SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR,
                             NV_PFIFO_CACHE1_DMA_STATE_ERROR_CALL);
                    break;
                } else {
                    *dma_subroutine = dma_get_v;
                    SET_MASK(*dma_subroutine,
                             NV_PFIFO_CACHE1_DMA_SUBROUTINE_STATE, 1);
                    dma_get_v = word & 0xfffffffc;
                    NV2A_DPRINTF("pb CALL 0x%08X\n", dma_get_v);
                }
            } else if (word == 0x00020000) {
                /* return */
                if (!subroutine_state) {
                    SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR,
                             NV_PFIFO_CACHE1_DMA_STATE_ERROR_RETURN);
                    // break;
                } else {
                    dma_get_v = *dma_subroutine & 0xfffffffc;
                    SET_MASK(*dma_subroutine,
                             NV_PFIFO_CACHE1_DMA_SUBROUTINE_STATE, 0);
                    NV2A_DPRINTF("pb RET 0x%08X\n", dma_get_v);
                }
            } else if ((word & 0xe0030003) == 0) {
                /* increasing methods */
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD,
                         (word & 0x1fff) >> 2 );
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_SUBCHANNEL,
                         (word >> 13) & 7);
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_COUNT,
                         (word >> 18) & 0x7ff);
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_TYPE,
                         NV_PFIFO_CACHE1_DMA_STATE_METHOD_TYPE_INC);
                *dma_dcount = 0;
            } else if ((word & 0xe0030003) == 0x40000000) {
                /* non-increasing methods */
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD,
                         (word & 0x1fff) >> 2 );
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_SUBCHANNEL,
                         (word >> 13) & 7);
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_COUNT,
                         (word >> 18) & 0x7ff);
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_METHOD_TYPE,
                         NV_PFIFO_CACHE1_DMA_STATE_METHOD_TYPE_NON_INC);
                *dma_dcount = 0;
            } else {
                NV2A_DPRINTF("pb reserved cmd 0x%08X - 0x%08X\n",
                             dma_get_v, word);
                SET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR,
                         NV_PFIFO_CACHE1_DMA_STATE_ERROR_RESERVED_CMD);
                // break;
                assert(false);
            }
        }

        *dma_get = dma_get_v;

        if (GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR)) {
            break;
        }
    }

    // Release the batched pgraph_lock acquired before the loop.
    qemu_mutex_unlock(&d->pgraph.pgraph_lock);

    // NV2A_DPRINTF("DMA pusher done: max 0x%08X, 0x%08X - 0x%08X\n",
    //      dma_len, control->dma_get, control->dma_put);

    uint32_t error = GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR);
    if (error) {
        NV2A_DPRINTF("pb error: %d\n", error);
        EmuLog(LOG_LEVEL::WARNING, "PFIFO DMA pusher error: %d", error);

        SET_MASK(*dma_push, NV_PFIFO_CACHE1_DMA_PUSH_STATUS, 1); /* suspended */

        // d->pfifo.pending_interrupts |= NV_PFIFO_INTR_0_DMA_PUSHER;
        // update_irq(d);
    }
}

int pfifo_pusher_thread(NV2AState *d)
{
    g_AffinityPolicy->SetAffinityOther();
    CxbxSetThreadName("Cxbx NV2A FIFO pusher");
    // Pusher now dispatches methods directly to PGRAPH (bypassing CACHE1), so
    // draw callbacks (CxbxUpdateNativeD3DResources) must not attempt a
    // re-entrant pfifo_flush_to_pgraph which would deadlock on pfifo_lock.
    CxbxSetPullerContext(true);

    qemu_mutex_lock(&d->pfifo.pfifo_lock);
    while (true) {
        {
            CXBX_PROFILE_SCOPE(PROF_PFIFO_PUSHER);
            pfifo_run_pusher(d);
        }

        // flush_requested is no longer set by pfifo_flush_to_pgraph (flush now
        // processes the pushbuffer inline on the calling thread).  The check
        // and signal below are kept as a safety net in case any future code
        // path restores the old protocol, but they are normally dead code.
        if (d->pfifo.flush_requested) {
            d->pfifo.flush_requested = false;
            qemu_cond_signal(&d->pfifo.flush_complete_cond);
        }

        qemu_cond_wait(&d->pfifo.pusher_cond, &d->pfifo.pfifo_lock);

        if (d->exiting) {
            break;
        }
    }
    qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	return 0;
}

unsigned int ramht_size(NV2AState *d)
{
	return 
		1 << (GET_MASK(d->pfifo.regs[RI(NV_PFIFO_RAMHT)], NV_PFIFO_RAMHT_SIZE_MASK) + 12);
}

static uint32_t ramht_hash(NV2AState *d, uint32_t handle)
{
	/* XXX: Think this is different to what nouveau calculates... */
	unsigned int bits = ffs(ramht_size(d)) - 2;

	uint32_t hash = 0;
	while (handle) {
		hash ^= (handle & ((1 << bits) - 1));
		handle >>= bits;
	}

    unsigned int channel_id = GET_MASK(d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH1)],
                                       NV_PFIFO_CACHE1_PUSH1_CHID);
    hash ^= channel_id << (bits - 4);

	return hash;
}

static RAMHTEntry ramht_lookup(NV2AState *d, uint32_t handle)
{
	uint32_t hash = ramht_hash(d, handle);
	assert(hash * 8 < ramht_size(d));

	xbox::addr_xt ramht_address =
		GET_MASK(d->pfifo.regs[RI(NV_PFIFO_RAMHT)],
			NV_PFIFO_RAMHT_BASE_ADDRESS_MASK) << 12;

	uint8_t *entry_ptr = d->pramin.ramin_ptr + ramht_address + hash * 8;

	uint32_t entry_handle = ldl_le_p((uint32_t*)entry_ptr);
	uint32_t entry_context = ldl_le_p((uint32_t*)(entry_ptr + 4));

	RAMHTEntry entry;
	entry.handle = entry_handle;
	entry.instance = (entry_context & NV_RAMHT_INSTANCE) << 4;
	entry.engine = (FIFOEngine)((entry_context & NV_RAMHT_ENGINE) >> 16);
	entry.channel_id = (entry_context & NV_RAMHT_CHID) >> 24;
	entry.valid = entry_context & NV_RAMHT_STATUS;

	return entry;
}
