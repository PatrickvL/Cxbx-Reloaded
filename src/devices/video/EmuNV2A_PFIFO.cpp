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

static RAMHTEntry ramht_lookup(NV2AState *d, uint32_t handle); // forward declaration

/* PFIFO - MMIO and DMA FIFO submission to PGRAPH and VPE */
DEVICE_READ32(PFIFO)
{
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

    // TODO
    // CacheEntry working_cache[NV2A_CACHE1_SIZE];
    // int working_cache_size = 0;
    // pull everything into our own queue

    // TODO think more about locking

    while (true) {
        if (!GET_MASK(*pull0, NV_PFIFO_CACHE1_PULL0_ACCESS)) return;

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

        // Process pushbuffer methods into PGRAPH register state so the
        // RC/VS interpreters can read from it.  Skip object binding
        // (method 0) and object-reference methods (0x180..0x1FF) which
        // require RAMHT lookups.  Skip context-switch / FIFO-wait —
        // Xbox uses a single GPU channel, so it's safe to write directly.
        if (method >= 0x100 && !(method >= 0x180 && method < 0x200)) {
            qemu_mutex_lock(&d->pgraph.pgraph_lock);
            pgraph_handle_method(d, subchannel, method, parameter);
            qemu_mutex_unlock(&d->pgraph.pgraph_lock);
        }

    }
}

// Defined in HostSync.cpp — marks the current thread as the PFIFO puller
// so CxbxUpdateNativeD3DResources skips pfifo_flush (prevents deadlock).
extern void CxbxSetPullerContext(bool active);

int pfifo_puller_thread(NV2AState *d)
{
    g_AffinityPolicy->SetAffinityOther();
    CxbxSetThreadName("Cxbx NV2A FIFO puller");
    CxbxSetPullerContext(true);

    qemu_mutex_lock(&d->pfifo.pfifo_lock);
    while (true) {
        pfifo_run_puller(d);

        // If the HLE thread is waiting for a PFIFO flush, signal it now
        // that CACHE1 has been drained.  The waiter will re-check whether
        // the DMA pusher also needs another cycle.
        if (d->pfifo.flush_requested) {
            qemu_cond_signal(&d->pfifo.flush_complete_cond);
        }

        qemu_cond_wait(&d->pfifo.puller_cond, &d->pfifo.pfifo_lock);

        if (d->exiting) {
            break;
        }
    }
    qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	return NULL;
}

// ---------------------------------------------------------------------------
// pfifo_flush_to_pgraph  --  block until all pending pushbuffer commands have
// been pushed into CACHE1 by the DMA pusher AND pulled/dispatched to PGRAPH
// by the puller.  Called from the HLE thread before each draw.
// ---------------------------------------------------------------------------
void pfifo_flush_to_pgraph(NV2AState *d)
{
    qemu_mutex_lock(&d->pfifo.pfifo_lock);

    while (true) {
        uint32_t status  = d->pfifo.regs[RI(NV_PFIFO_CACHE1_STATUS)];
        bool cache1_empty = (status & NV_PFIFO_CACHE1_STATUS_LOW_MARK) != 0;
        bool dma_idle     = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)]
                         == d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];

        if (cache1_empty && dma_idle)
            break;

        // Tell the puller to signal us after its next drain cycle.
        d->pfifo.flush_requested = true;

        // Wake both threads so the pusher can feed CACHE1 and the puller
        // can drain it.  Harmless if either thread has nothing to do.
        qemu_cond_signal(&d->pfifo.pusher_cond);
        qemu_cond_signal(&d->pfifo.puller_cond);

        // Release pfifo_lock and sleep until the puller finishes a drain.
        qemu_cond_wait(&d->pfifo.flush_complete_cond, &d->pfifo.pfifo_lock);
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

    uint32_t *status = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_STATUS)];
    uint32_t *get_reg = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_GET)];
    uint32_t *put_reg = &d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUT)];

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
            /* full */
            if (*status & NV_PFIFO_CACHE1_STATUS_HIGH_MARK) return;


            /* data word of methods command */
            d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_DATA_SHADOW)] = word;

            uint32_t put = *put_reg;
            uint32_t get = *get_reg;

            assert((method & 3) == 0);
            uint32_t method_entry = 0;
            SET_MASK(method_entry, NV_PFIFO_CACHE1_METHOD_ADDRESS, method >> 2);
            SET_MASK(method_entry, NV_PFIFO_CACHE1_METHOD_TYPE, method_type);
            SET_MASK(method_entry, NV_PFIFO_CACHE1_METHOD_SUBCHANNEL, method_subchannel);

            // NV2A_DPRINTF("push %d 0x%08X 0x%08X - subch %d\n", put/4, method_entry, word, method_subchannel);

            assert(put < 128*4 && (put%4) == 0);
            d->pfifo.regs[RI(NV_PFIFO_CACHE1_METHOD + put*2)] = method_entry;
            d->pfifo.regs[RI(NV_PFIFO_CACHE1_DATA + put*2)] = word;

            uint32_t new_put = (put+4) & 0x1fc;
            *put_reg = new_put;
            if (new_put == get) {
                // set high mark
                *status |= NV_PFIFO_CACHE1_STATUS_HIGH_MARK;
            }
            if (*status & NV_PFIFO_CACHE1_STATUS_LOW_MARK) {
                // unset low mark
                *status &= ~NV_PFIFO_CACHE1_STATUS_LOW_MARK;
                // signal puller
                qemu_cond_signal(&d->pfifo.puller_cond);
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

    // NV2A_DPRINTF("DMA pusher done: max 0x%08X, 0x%08X - 0x%08X\n",
    //      dma_len, control->dma_get, control->dma_put);

    uint32_t error = GET_MASK(*dma_state, NV_PFIFO_CACHE1_DMA_STATE_ERROR);
    if (error) {
        NV2A_DPRINTF("pb error: %d\n", error);
        assert(false);

        SET_MASK(*dma_push, NV_PFIFO_CACHE1_DMA_PUSH_STATUS, 1); /* suspended */

        // d->pfifo.pending_interrupts |= NV_PFIFO_INTR_0_DMA_PUSHER;
        // update_irq(d);
    }
}

int pfifo_pusher_thread(NV2AState *d)
{
    g_AffinityPolicy->SetAffinityOther();
    CxbxSetThreadName("Cxbx NV2A FIFO pusher");

    qemu_mutex_lock(&d->pfifo.pfifo_lock);
    while (true) {
        pfifo_run_pusher(d);
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
