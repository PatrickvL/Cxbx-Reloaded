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
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a_user.c
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

/* USER - PFIFO MMIO and DMA submission area */
DEVICE_READ32(USER)
{
	unsigned int channel_id = addr >> 16;
	assert(channel_id < NV2A_NUM_CHANNELS);

	// Fast path for DMA_GET reads.  When the DMA pusher cannot process
	// (access flags not set in HLE mode), advance GET to PUT so polls
	// like BlockUntilIdle() return immediately.  When the pusher CAN
	// process, drain pending commands inline so the native polling loop
	// sees GET advance and exits naturally.
	if ((addr & 0xFFFF) == NV_USER_DMA_GET) {
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
				// Drain pending commands inline — enables native
				// BlockUntilIdle polling to work without a patch.
				pfifo_flush_to_pgraph(d);
				get_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
			}
		}
		uint32_t result = get_v;
		DEVICE_READ32_END(USER);
	}

	// Fast path for NV_USER_REF reads (reference counter).
	// NV_USER_REF (offset 0x48) is defined in xemu's nv2a_regs.h and
	// handled in xemu's user.c — it maps to NV_PFIFO_CACHE1_REF.
	// The value is updated by REF_CNT (method 0x0050, documented by
	// envytools: hw/fifo/puller.html#syncing-with-host-reference-counter).
	// Since we process commands inline on DMA_PUT writes, REF should
	// already be current here (GET == PUT).  This flush is a safety net
	// that rarely triggers in practice — xemu omits it entirely.
	if ((addr & 0xFFFF) == NV_USER_REF) {
		uint32_t get_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
		uint32_t put_v = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];
		if (get_v != put_v) {
			pfifo_flush_to_pgraph(d);
		}
		uint32_t result = d->pfifo.regs[RI(NV_PFIFO_CACHE1_REF)];
		DEVICE_READ32_END(USER);
	}

	qemu_mutex_lock(&d->pfifo.pfifo_lock);

	uint32_t channel_modes = d->pfifo.regs[RI(NV_PFIFO_MODE)];

	uint32_t result = 0;
	if (channel_modes & (1 << channel_id)) {
		/* DMA Mode */

		unsigned int cur_channel_id =
			GET_MASK(d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH1)],
				NV_PFIFO_CACHE1_PUSH1_CHID);

		if (channel_id == cur_channel_id) {
			switch(addr & 0xFFFF) { // Was DEVICE_READ32_SWITCH()
				case NV_USER_DMA_PUT:
					result = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)];
					break;
				case NV_USER_DMA_GET:
					result = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)];
					break;
				case NV_USER_REF:
					result = d->pfifo.regs[RI(NV_PFIFO_CACHE1_REF)];
					break;
				default:
					DEBUG_READ32_UNHANDLED(USER);
					break;
			}
		} else {
			/* ramfc */
			assert(false);
		}
	} else {
		/* PIO Mode */
		assert(false);
	}

	qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	DEVICE_READ32_END(USER);
}

DEVICE_WRITE32(USER)
{
	unsigned int channel_id = addr >> 16;
	assert(channel_id < NV2A_NUM_CHANNELS);

	qemu_mutex_lock(&d->pfifo.pfifo_lock);

	uint32_t channel_modes = d->pfifo.regs[RI(NV_PFIFO_MODE)];
	if (channel_modes & (1 << channel_id)) {
		/* DMA Mode */
		unsigned int cur_channel_id =
			GET_MASK(d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH1)],
				NV_PFIFO_CACHE1_PUSH1_CHID);

		if (channel_id == cur_channel_id) {
			switch (addr & 0xFFFF) {
			case NV_USER_DMA_PUT: {
				d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUT)] = value;
				// Process commands inline immediately.  Native D3D runtime
				// may wait for completion signals (FLIP_STALL, semaphore)
				// BEFORE polling DMA_GET.  If we don't process here, those
				// signals never fire and the game thread deadlocks.
				{
					uint32_t push0    = d->pfifo.regs[RI(NV_PFIFO_CACHE1_PUSH0)];
					uint32_t dma_push = d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_PUSH)];
					bool pusher_can_run = GET_MASK(push0, NV_PFIFO_CACHE1_PUSH0_ACCESS)
					                   && GET_MASK(dma_push, NV_PFIFO_CACHE1_DMA_PUSH_ACCESS)
					                   && !GET_MASK(dma_push, NV_PFIFO_CACHE1_DMA_PUSH_STATUS);
					if (pusher_can_run) {
						CxbxSetPullerContext(true);
						pfifo_run_pusher(d);
						CxbxSetPullerContext(false);
					}
				}
				break;
			}
			case NV_USER_DMA_GET:
				d->pfifo.regs[RI(NV_PFIFO_CACHE1_DMA_GET)] = value;
				break;
			case NV_USER_REF:
				d->pfifo.regs[RI(NV_PFIFO_CACHE1_REF)] = value;
				break;
			default:
				assert(false);
				break;
			}

            // Kick puller thread (for auto-present fallback on raw-pushbuffer
            // games without explicit FLIP_STALL).  Do NOT signal pusher_cond:
            // command processing is driven exclusively by inline flushes
            // (pfifo_flush_to_pgraph called from DMA_GET reads and before draws).
            // Signaling the pusher would cause it to race for pfifo_lock,
            // introducing intermittent stalls in the game thread.
            qemu_cond_broadcast(&d->pfifo.puller_cond);
		} else {
			/* ramfc */
			assert(false);
		}
	} else {
		/* PIO Mode */
		assert(false);
	}

    qemu_mutex_unlock(&d->pfifo.pfifo_lock);

	DEVICE_WRITE32_END(USER);
}
