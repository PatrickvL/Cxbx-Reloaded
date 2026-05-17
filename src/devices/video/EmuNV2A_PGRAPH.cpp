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
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a_pgraph.c
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

// FIXME
#define qemu_mutex_lock_iothread()
#define qemu_mutex_unlock_iothread()

// Xbox uses 4 KiB pages
#define TARGET_PAGE_BITS 12
#define TARGET_PAGE_SIZE (1 << TARGET_PAGE_BITS)
#define TARGET_PAGE_MASK ~(TARGET_PAGE_SIZE - 1)
#define TARGET_PAGE_ALIGN(addr) (((addr) + TARGET_PAGE_SIZE - 1) & TARGET_PAGE_MASK)

// ---- NV097 method trace infrastructure ----
// Scans push buffer commands written by Xbox D3D API trampolines.
// pgraph_trace_begin/end bracket a trampoline call; pgraph_trace_log_pushbuffer
// decodes the NV097 commands that were written between old/new pPut positions.
static FILE *pgraph_trace_file = nullptr;
static int pgraph_trace_count = 0;
static const int PGRAPH_TRACE_MAX = 50000; // Max trace entries

extern const char *NV2AMethodToString(DWORD dwMethod); // implemented in PushBuffer.cpp

// Cache the Xbox CDevice pPut pointer once, used by all trace calls
static uint32_t **pgraph_trace_pDeviceObj = nullptr;
static bool pgraph_trace_device_looked_up = false;

static void pgraph_trace_ensure_device()
{
	if (pgraph_trace_device_looked_up) return;
	pgraph_trace_device_looked_up = true;
	void *pDeviceGlobal = GetXboxSymbolPointer("D3D_g_pDevice");
	if (pDeviceGlobal && !IsBadReadPtr(pDeviceGlobal, 4)) {
		uint32_t devAddr = *(uint32_t*)pDeviceGlobal;
		if (devAddr && !IsBadReadPtr((void*)devAddr, 4)) {
			pgraph_trace_pDeviceObj = (uint32_t**)devAddr;
		}
	}
}

static void pgraph_trace_ensure_file()
{
	if (!pgraph_trace_file) {
		pgraph_trace_file = fopen("pgraph_trace.txt", "w");
	}
}

uint32_t pgraph_trace_read_pput()
{
	pgraph_trace_ensure_device();
	if (!pgraph_trace_pDeviceObj) return 0;
	// CDevice+0x00 = pPut (virtual address in contiguous Xbox memory)
	return (uint32_t)(uintptr_t)pgraph_trace_pDeviceObj[0];
}

void pgraph_trace_log_pushbuffer(const char *tag, uint32_t pPut_before, uint32_t pPut_after)
{
	if (pgraph_trace_count >= PGRAPH_TRACE_MAX) return;
	if (pPut_before == pPut_after) return; // No commands written
	pgraph_trace_ensure_file();
	if (!pgraph_trace_file) return;

	fprintf(pgraph_trace_file, ">>> %s\n", tag);

	// pPut values are virtual addresses in Xbox contiguous memory (e.g. 0x83Fxxxxx)
	// They can be read as host pointers directly (contiguous memory is identity-mapped)
	const uint32_t *pStart = (const uint32_t*)(uintptr_t)pPut_before;
	const uint32_t *pEnd   = (const uint32_t*)(uintptr_t)pPut_after;

	// Handle ring buffer wrap-around: for now just scan forward
	// (most API calls write small amounts, won't wrap)
	if (pEnd <= pStart) {
		fprintf(pgraph_trace_file, "  (wrap or empty: before=0x%08X after=0x%08X)\n",
			pPut_before, pPut_after);
		fprintf(pgraph_trace_file, "<<< %s\n", tag);
		fflush(pgraph_trace_file);
		return;
	}

	const uint32_t *p = pStart;
	while (p < pEnd && pgraph_trace_count < PGRAPH_TRACE_MAX) {
		uint32_t cmd = *p;

		// Check for PUSH_TYPE: 0 = method submission, non-zero = jump/call
		DWORD pushType = PUSH_TYPE(cmd);
		if (pushType != 0) {
			// Jump/call command — skip it for tracing
			fprintf(pgraph_trace_file, "  [JUMP/CALL 0x%08X]\n", cmd);
			p++;
			pgraph_trace_count++;
			continue;
		}

		DWORD dwMethod, dwSubCh, dwCount;
		D3DPUSH_DECODE(cmd, dwMethod, dwSubCh, dwCount);

		// Non-incrementing flag: bit 30 means all parameters go to same method
		bool bNonInc = (cmd & 0x40000000) != 0;

		for (DWORD i = 0; i < dwCount && (p + 1 + i) < pEnd; i++) {
			uint32_t param = *(p + 1 + i);
			DWORD meth = bNonInc ? dwMethod : (dwMethod + i * 4);
			const char *name = NV2AMethodToString(meth);
			if (name) {
				fprintf(pgraph_trace_file, "  0x%04X %-40s = 0x%08X\n", meth, name, param);
			} else {
				fprintf(pgraph_trace_file, "  0x%04X (subchan=%d)%*s = 0x%08X\n",
					meth, dwSubCh, 32, "", param);
			}
			pgraph_trace_count++;
		}

		p += 1 + dwCount; // Skip command DWORD + parameter DWORDs
	}

	fprintf(pgraph_trace_file, "<<< %s\n", tag);
	fflush(pgraph_trace_file);
}

// Legacy begin/end no-ops (push-buffer-scan approach supersedes PGRAPH-level trace)
void pgraph_trace_begin(const char *) {}
void pgraph_trace_end() {}
void pgraph_trace_close()
{
	if (pgraph_trace_file) {
		fclose(pgraph_trace_file);
		pgraph_trace_file = nullptr;
	}
}
// ---- End trace infrastructure ----

#include "nv2a_pgraph_backend.h"

PgraphBackend g_pgraph_backend = {};

// Set true the first time the title issues an explicit NV097_FLIP_STALL.
// Once observed, the puller's auto-present fallback (intended for raw push
// buffer games that never flip) is disabled to avoid spurious mid-frame
// presents between explicit flips, which causes flicker at half frame rate.
bool g_pgraph_explicit_flip_stall_seen = false;

void pgraph_handle_method(NV2AState *d, unsigned int subchannel, unsigned int method, uint32_t parameter);
static void pgraph_log_method(unsigned int subchannel, unsigned int graphics_class, unsigned int method, uint32_t parameter);
static void pgraph_allocate_inline_buffer_vertices(PGRAPHState *pg, unsigned int attr);
static void pgraph_finish_inline_buffer_vertex(PGRAPHState *pg);
static bool pgraph_get_color_write_enabled(PGRAPHState *pg);
static bool pgraph_get_zeta_write_enabled(PGRAPHState *pg);
static void pgraph_set_surface_dirty(PGRAPHState *pg, bool color, bool zeta);
static void pgraph_apply_anti_aliasing_factor(PGRAPHState *pg, unsigned int *width, unsigned int *height);
static void pgraph_get_surface_dimensions(PGRAPHState *pg, unsigned int *width, unsigned int *height);
static float convert_f16_to_float(uint16_t f16);
static float convert_f24_to_float(uint32_t f24);
static uint8_t* convert_texture_data(const unsigned int color_format, const uint8_t *data, const uint8_t *palette_data, const unsigned int width, const unsigned int height, const unsigned int depth, const unsigned int row_pitch, const unsigned int slice_pitch);
static unsigned int kelvin_map_stencil_op(uint32_t parameter);
static unsigned int kelvin_map_polygon_mode(uint32_t parameter);
static unsigned int kelvin_map_texgen(uint32_t parameter, unsigned int channel);

/* PGRAPH - accelerated 2d/3d drawing engine */

static uint32_t pgraph_rdi_read(PGRAPHState *pg,
                                unsigned int select, unsigned int address)
{
    uint32_t r = 0;
    switch(select) {
    case RDI_INDEX_VTX_CONSTANTS0:
    case RDI_INDEX_VTX_CONSTANTS1:
        assert((address / 4) < NV2A_VERTEXSHADER_CONSTANTS);
        r = pg->xf.xfctx[address / 4][3 - address % 4];
        break;
    default:
        fprintf(stderr, "nv2a: unknown rdi read select 0x%x address 0x%x\n",
                select, address);
        assert(false);
        break;
    }
    return r;
}

static void pgraph_rdi_write(PGRAPHState *pg,
                             unsigned int select, unsigned int address,
                             uint32_t val)
{
    switch(select) {
    case RDI_INDEX_VTX_CONSTANTS0:
    case RDI_INDEX_VTX_CONSTANTS1:
        assert(false); /* Untested */
        assert((address / 4) < NV2A_VERTEXSHADER_CONSTANTS);
        {
            unsigned idx = address / 4;
			unsigned slot = 3 - address % 4;
            if (pg->xf.xfctx[idx][slot] != val) {
                pg->xf.xfctx[idx][slot] = val;
                pg->xf.xfctx_dirty[idx / 32] |= (1u << (idx % 32));
                pg->xf.xfctx_generation++;
            }
        }
        break;
    default:
        NV2A_DPRINTF("unknown rdi write select 0x%x, address 0x%x, val 0x%08x\n",
                     select, address, val);
        break;
    }
}

DEVICE_READ32(PGRAPH)
{
	qemu_mutex_lock(&d->pgraph.pgraph_lock);

    PGRAPHState *pg = &d->pgraph;
	DEVICE_READ32_SWITCH() {
	case NV_PGRAPH_INTR:
		result = pg->pending_interrupts;
		break;
	case NV_PGRAPH_INTR_EN:
		result = pg->enabled_interrupts;
		break;
    case NV_PGRAPH_RDI_DATA: {
        unsigned int select = GET_MASK(pg->regs[RI(NV_PGRAPH_RDI_INDEX)],
                                       NV_PGRAPH_RDI_INDEX_SELECT);
        int address = GET_MASK(pg->regs[RI(NV_PGRAPH_RDI_INDEX)],
                                        NV_PGRAPH_RDI_INDEX_ADDRESS);

        result = pgraph_rdi_read(pg, select, address);

        /* FIXME: Overflow into select? */
        assert(address < GET_MASK(NV_PGRAPH_RDI_INDEX_ADDRESS,
                                  NV_PGRAPH_RDI_INDEX_ADDRESS));
        SET_MASK(pg->regs[RI(NV_PGRAPH_RDI_INDEX)],
                 NV_PGRAPH_RDI_INDEX_ADDRESS, address + 1);
        break;
    }
	default:
		DEVICE_READ32_REG(pgraph); // Was : DEBUG_READ32_UNHANDLED(PGRAPH);
	}

	qemu_mutex_unlock(&pg->pgraph_lock);

//    reg_log_read(NV_PGRAPH, addr, r);

	DEVICE_READ32_END(PGRAPH);
}

DEVICE_WRITE32(PGRAPH)
{
    PGRAPHState *pg = &d->pgraph;
//    reg_log_write(NV_PGRAPH, addr, val);

	qemu_mutex_lock(&pg->pgraph_lock);

	switch (addr) {
	case NV_PGRAPH_INTR:
		pg->pending_interrupts &= ~value;
		qemu_cond_broadcast(&pg->interrupt_cond);
		break;
	case NV_PGRAPH_INTR_EN:
		pg->enabled_interrupts = value;
		break;
	case NV_PGRAPH_INCREMENT:
		if (value & NV_PGRAPH_INCREMENT_READ_3D) {
			SET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
				NV_PGRAPH_SURFACE_READ_3D,
				(GET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
					NV_PGRAPH_SURFACE_READ_3D) + 1)
				% GET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
					NV_PGRAPH_SURFACE_MODULO_3D));
			qemu_cond_broadcast(&pg->flip_3d);

			// For MMIO-only games (no pushbuffer FLIP_STALL), mark surface dirty
			// and wake the puller thread so its auto-present fires.  We can't call
			// g_pgraph_backend.flip_stall directly here because this runs on the DPC/system_events
			// thread, not the puller thread that owns the D3D11 context.
			if (!g_pgraph_explicit_flip_stall_seen) {
				d->pgraph.surface_color.draw_dirty = true;
				qemu_cond_broadcast(&d->pfifo.puller_cond);
			}
		}
		break;
    case NV_PGRAPH_RDI_DATA: {
        unsigned int select = GET_MASK(pg->regs[RI(NV_PGRAPH_RDI_INDEX)],
                                       NV_PGRAPH_RDI_INDEX_SELECT);
        int address = GET_MASK(pg->regs[RI(NV_PGRAPH_RDI_INDEX)],
                                        NV_PGRAPH_RDI_INDEX_ADDRESS);

        pgraph_rdi_write(pg, select, address, value);

        /* FIXME: Overflow into select? */
        assert(address < GET_MASK(NV_PGRAPH_RDI_INDEX_ADDRESS,
                                  NV_PGRAPH_RDI_INDEX_ADDRESS));
        SET_MASK(pg->regs[RI(NV_PGRAPH_RDI_INDEX)],
                 NV_PGRAPH_RDI_INDEX_ADDRESS, address + 1);
        break;
    }
	case NV_PGRAPH_CHANNEL_CTX_TRIGGER: {
		xbox::addr_xt context_address =
			GET_MASK(pg->regs[RI(NV_PGRAPH_CHANNEL_CTX_POINTER)],
				NV_PGRAPH_CHANNEL_CTX_POINTER_INST) << 4;

		if (value & NV_PGRAPH_CHANNEL_CTX_TRIGGER_READ_IN) {
			unsigned pgraph_channel_id =
				GET_MASK(pg->regs[RI(NV_PGRAPH_CTX_USER)], NV_PGRAPH_CTX_USER_CHID);

			NV2A_DPRINTF("PGRAPH: read channel %d context from %" HWADDR_PRIx "\n",
				pgraph_channel_id, context_address);

			uint8_t *context_ptr = d->pramin.ramin_ptr + context_address;
			uint32_t context_user = ldl_le_p((uint32_t*)context_ptr);

			NV2A_DPRINTF("    - CTX_USER = 0x%08X\n", context_user);

			pg->regs[RI(NV_PGRAPH_CTX_USER)] = context_user;
			// pgraph_set_context_user(d, context_user);
		}
		if (value & NV_PGRAPH_CHANNEL_CTX_TRIGGER_WRITE_OUT) {
			/* do stuff ... */
		}

		break;
    }
	case NV_PGRAPH_CTX_CONTROL:
		// CHID bit is managed by hardware channel switch (PFIFO puller), not by
		// CPU MMIO writes.  Preserve it so kernel register init doesn't clear it.
		pg->regs[RI(NV_PGRAPH_CTX_CONTROL)] = value
			| (pg->regs[RI(NV_PGRAPH_CTX_CONTROL)] & NV_PGRAPH_CTX_CONTROL_CHID);
		break;
	default: 
		DEVICE_WRITE32_REG(pgraph); // Was : DEBUG_WRITE32_UNHANDLED(PGRAPH);
		break;
	}

    // events
    switch (addr) {
    case NV_PGRAPH_FIFO:
        qemu_cond_broadcast(&pg->fifo_access_cond);
        break;
    }

	qemu_mutex_unlock(&pg->pgraph_lock);

	DEVICE_WRITE32_END(PGRAPH);
}


void pgraph_handle_method(NV2AState *d,
							unsigned int subchannel,
							unsigned int method,
							uint32_t parameter)
{
	unsigned int slot;

    PGRAPHState *pg = &d->pgraph;

    bool channel_valid =
        d->pgraph.regs[RI(NV_PGRAPH_CTX_CONTROL)] & NV_PGRAPH_CTX_CONTROL_CHID;
    assert(channel_valid);

    unsigned channel_id = GET_MASK(pg->regs[RI(NV_PGRAPH_CTX_USER)], NV_PGRAPH_CTX_USER_CHID);

	ContextSurfaces2DState *context_surfaces_2d = &pg->context_surfaces_2d;
	ImageBlitState *image_blit = &pg->image_blit;
	KelvinState *kelvin = &pg->kelvin;

    assert(subchannel < 8);

	if (method == NV_SET_OBJECT) {
        assert(parameter < d->pramin.ramin_size);
        uint8_t *obj_ptr = d->pramin.ramin_ptr + parameter;

        uint32_t ctx_1 = ldl_le_p((uint32_t*)obj_ptr);
        uint32_t ctx_2 = ldl_le_p((uint32_t*)(obj_ptr+4));
        uint32_t ctx_3 = ldl_le_p((uint32_t*)(obj_ptr+8));
        uint32_t ctx_4 = ldl_le_p((uint32_t*)(obj_ptr+12));
        uint32_t ctx_5 = parameter;

        pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + subchannel * 4)] = ctx_1;
        pg->regs[RI(NV_PGRAPH_CTX_CACHE2 + subchannel * 4)] = ctx_2;
        pg->regs[RI(NV_PGRAPH_CTX_CACHE3 + subchannel * 4)] = ctx_3;
        pg->regs[RI(NV_PGRAPH_CTX_CACHE4 + subchannel * 4)] = ctx_4;
        pg->regs[RI(NV_PGRAPH_CTX_CACHE5 + subchannel * 4)] = ctx_5;
    }

    // is this right?
    pg->regs[RI(NV_PGRAPH_CTX_SWITCH1)] = pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + subchannel * 4)];
    pg->regs[RI(NV_PGRAPH_CTX_SWITCH2)] = pg->regs[RI(NV_PGRAPH_CTX_CACHE2 + subchannel * 4)];
    pg->regs[RI(NV_PGRAPH_CTX_SWITCH3)] = pg->regs[RI(NV_PGRAPH_CTX_CACHE3 + subchannel * 4)];
    pg->regs[RI(NV_PGRAPH_CTX_SWITCH4)] = pg->regs[RI(NV_PGRAPH_CTX_CACHE4 + subchannel * 4)];
    pg->regs[RI(NV_PGRAPH_CTX_SWITCH5)] = pg->regs[RI(NV_PGRAPH_CTX_CACHE5 + subchannel * 4)];

    uint32_t graphics_class = GET_MASK(pg->regs[RI(NV_PGRAPH_CTX_SWITCH1)],
                                       NV_PGRAPH_CTX_SWITCH1_GRCLASS);

	// Logging is slow.. disable for now..
	//pgraph_log_method(subchannel, graphics_class, method, parameter);

    if (subchannel != 0) {
        // catches context switching issues on xbox d3d
        assert(graphics_class != 0x97);
    }

    /* ugly switch for now */
    switch (graphics_class) {

    case NV_CONTEXT_PATTERN: {
		switch (method) {
		case NV044_SET_MONOCHROME_COLOR0:
			pg->regs[RI(NV_PGRAPH_PATT_COLOR0)] = parameter;
			break;
		}
		
		break;
	}

    case NV_CONTEXT_SURFACES_2D: {
		switch (method) {
		case NV062_SET_OBJECT:
			context_surfaces_2d->object_instance = parameter;
			break;
		case NV062_SET_CONTEXT_DMA_IMAGE_SOURCE:
			context_surfaces_2d->dma_image_source = parameter;
			break;
		case NV062_SET_CONTEXT_DMA_IMAGE_DESTIN:
			context_surfaces_2d->dma_image_dest = parameter;
			break;
		case NV062_SET_COLOR_FORMAT:
			context_surfaces_2d->color_format = parameter;
			break;
		case NV062_SET_PITCH:
			context_surfaces_2d->source_pitch = parameter & 0xFFFF;
			context_surfaces_2d->dest_pitch = parameter >> 16;
			break;
		case NV062_SET_OFFSET_SOURCE:
			context_surfaces_2d->source_offset = parameter & 0x07FFFFFF;
			break;
		case NV062_SET_OFFSET_DESTIN:
			context_surfaces_2d->dest_offset = parameter & 0x07FFFFFF;
			break;
		default:
			EmuLog(LOG_LEVEL::WARNING, "Unknown NV_CONTEXT_SURFACES_2D Method: 0x%08X", method);
		}
	
		break; 
	}
	
	case NV_IMAGE_BLIT: {
		switch (method) {
		case NV09F_SET_OBJECT:
			image_blit->object_instance = parameter;
			break;
		case NV09F_SET_CONTEXT_SURFACES:
			image_blit->context_surfaces = parameter;
			break;
		case NV09F_SET_OPERATION:
			image_blit->operation = parameter;
			break;
		case NV09F_CONTROL_POINT_IN:
			image_blit->in_x = parameter & 0xFFFF;
			image_blit->in_y = parameter >> 16;
			break;
		case NV09F_CONTROL_POINT_OUT:
			image_blit->out_x = parameter & 0xFFFF;
			image_blit->out_y = parameter >> 16;
			break;
		case NV09F_SIZE:
			image_blit->width = parameter & 0xFFFF;
			image_blit->height = parameter >> 16;

			/* I guess this kicks it off? */
			if (image_blit->operation == NV09F_SET_OPERATION_SRCCOPY) {

				NV2A_DPRINTF_IF(true, "NV09F_SET_OPERATION_SRCCOPY");

				ContextSurfaces2DState *context_surfaces = context_surfaces_2d;
				assert(context_surfaces->object_instance
					== image_blit->context_surfaces);

				unsigned int bytes_per_pixel;
				switch (context_surfaces->color_format) {
				case NV062_SET_COLOR_FORMAT_LE_Y8:
					bytes_per_pixel = 1;
					break;
				case NV062_SET_COLOR_FORMAT_LE_R5G6B5:
					bytes_per_pixel = 2;
					break;
				case NV062_SET_COLOR_FORMAT_LE_A8R8G8B8:
					bytes_per_pixel = 4;
					break;
				default:
					printf("Unknown blit surface format: 0x%x\n", context_surfaces->color_format);
					assert(false);
					break;
				}

				xbox::addr_xt source_dma_len, dest_dma_len;
				uint8_t *source, *dest;

				source = (uint8_t*)nv_dma_map(d, context_surfaces->dma_image_source,
												&source_dma_len);
				assert(context_surfaces->source_offset < source_dma_len);
				source += context_surfaces->source_offset;

				dest = (uint8_t*)nv_dma_map(d, context_surfaces->dma_image_dest,
												&dest_dma_len);
				assert(context_surfaces->dest_offset < dest_dma_len);
				dest += context_surfaces->dest_offset;

				NV2A_DPRINTF("  - 0x%tx -> 0x%tx\n", source - d->vram_ptr,
														dest - d->vram_ptr);

				unsigned int y;
				for (y = 0; y<image_blit->height; y++) {
					uint8_t *source_row = source
						+ (image_blit->in_y + y) * context_surfaces->source_pitch
						+ image_blit->in_x * bytes_per_pixel;

					uint8_t *dest_row = dest
						+ (image_blit->out_y + y) * context_surfaces->dest_pitch
						+ image_blit->out_x * bytes_per_pixel;

					memmove(dest_row, source_row,
						image_blit->width * bytes_per_pixel);
				}

			} else {
				assert(false);
			}

			break;
		default:
			EmuLog(LOG_LEVEL::WARNING, "Unknown NV_IMAGE_BLIT Method: 0x%08X", method);
		}
		break;
	}

	case NV_KELVIN_PRIMITIVE: {
		// Data-driven dispatch: handles reg copies and masked writes for all
		// table-registered methods. Returns old register value before overwrite.
		uint32_t old_reg = nv097_dispatch_method(pg, method, parameter);

		switch (method) {
		case NV097_SET_OBJECT:
			kelvin->object_instance = parameter;
			break;

		case NV097_NO_OPERATION:
			/* The bios uses nop as a software method call -
			 * it seems to expect a notify interrupt if the parameter isn't 0.
			 * According to a nouveau guy it should still be a nop regardless
			 * of the parameter. It's possible a debug register enables this,
			 * but nothing obvious sticks out. Weird.
			 */
			if (parameter != 0) {
				assert(!(pg->pending_interrupts & NV_PGRAPH_INTR_ERROR));

				SET_MASK(pg->regs[RI(NV_PGRAPH_TRAPPED_ADDR)],
					NV_PGRAPH_TRAPPED_ADDR_CHID, channel_id);
				SET_MASK(pg->regs[RI(NV_PGRAPH_TRAPPED_ADDR)],
					NV_PGRAPH_TRAPPED_ADDR_SUBCH, subchannel);
				SET_MASK(pg->regs[RI(NV_PGRAPH_TRAPPED_ADDR)],
					NV_PGRAPH_TRAPPED_ADDR_MTHD, method);
				pg->regs[RI(NV_PGRAPH_TRAPPED_DATA_LOW)] = parameter;
				pg->regs[RI(NV_PGRAPH_NSOURCE)] = NV_PGRAPH_NSOURCE_NOTIFICATION; /* TODO: check this */
				pg->pending_interrupts |= NV_PGRAPH_INTR_ERROR;

				qemu_mutex_unlock(&pg->pgraph_lock);
				qemu_mutex_lock_iothread();
				update_irq(d);
				qemu_mutex_lock(&pg->pgraph_lock);
				qemu_mutex_unlock_iothread();

				while (pg->pending_interrupts & NV_PGRAPH_INTR_ERROR) {
					// Use timed wait as a safety net: if the DPC signal was lost
					// (due to any unforeseen race), we re-signal after 50ms rather
					// than hanging indefinitely. Normal path returns in <1ms.
					if (qemu_cond_timedwait(&pg->interrupt_cond, &pg->pgraph_lock, 50)) {
						// Timed out — re-signal DPC thread in case the original was lost
						qemu_mutex_unlock(&pg->pgraph_lock);
						update_irq(d);
						qemu_mutex_lock(&pg->pgraph_lock);
					}
				}
			}
			break;

		// NV097_WAIT_FOR_IDLE: On real HW this drains the 3D pipeline before PFIFO
		// continues. In our architecture, D3D11 draw calls execute synchronously on
		// the puller thread, so the pipeline is already idle — intentional no-op.

		case NV097_FLIP_INCREMENT_WRITE: {
			NV2A_DPRINTF("flip increment write %d -> ",
				GET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
					NV_PGRAPH_SURFACE_WRITE_3D));
			SET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
				NV_PGRAPH_SURFACE_WRITE_3D,
				(GET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
					NV_PGRAPH_SURFACE_WRITE_3D) + 1)
				% GET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
					NV_PGRAPH_SURFACE_MODULO_3D));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SURFACE]++;
			NV2A_DPRINTF("%d\n",
				GET_MASK(pg->regs[RI(NV_PGRAPH_SURFACE)],
					NV_PGRAPH_SURFACE_WRITE_3D));

			break;
		}
		case NV097_FLIP_STALL:
			// Title is using explicit flips — disable puller auto-present fallback.
			g_pgraph_explicit_flip_stall_seen = true;

			// Trigger host present via the flip_stall plugin callback
			if (g_pgraph_backend.flip_stall != nullptr) {
				// Clear draw_dirty so the auto-present in the puller loop
				// doesn't fire again after this explicit FLIP_STALL present.
				d->pgraph.surface_color.draw_dirty = false;
				g_pgraph_backend.flip_stall(d);
			}

			// VBlank-gated frame pacing: wait until the next VBlank deadline.
			// SleepPrecise uses adaptive yielding (SwitchToThread with EMA tracking)
			// which donates CPU time to other threads while waiting, then does a
			// final spin for sub-yield precision.
			// Local anchor prevents drift without conflicting with the VBlank
			// interrupt's writes to vblank_last_qpc.
			{
				static int64_t s_flipStallAnchor = 0;
				unsigned int totalLines = pcrtc_get_total_lines(d);
				unsigned int refreshRate = pcrtc_get_refresh_rate(d, totalLines);
				int64_t vblankPeriodTicks = HostQPCFrequency / refreshRate;

				// Seed anchor from the real VBlank timestamp on first call,
				// or reseed if it's fallen too far behind (e.g. after a stall).
				// SleepPrecise handles the "already behind" case internally
				// (returns current QPC immediately), so no pre-check needed.
				int64_t lastVBlank = d->vblank_last_qpc.load(std::memory_order_acquire);
				if (s_flipStallAnchor == 0) {
					LARGE_INTEGER now;
					QueryPerformanceCounter(&now);
					if (now.QuadPart - lastVBlank > vblankPeriodTicks * 2) {
						s_flipStallAnchor = lastVBlank;
					}
				}

				if (s_flipStallAnchor > 0) {
					int64_t nextVBlankQPC = s_flipStallAnchor + vblankPeriodTicks;
					qemu_mutex_unlock(&d->pgraph.pgraph_lock);
					s_flipStallAnchor = SleepPrecise(nextVBlankQPC);
					qemu_mutex_lock(&d->pgraph.pgraph_lock);
				}
			}

			NV2A_DPRINTF("flip stall done\n");
			break;

		case NV097_SET_CONTEXT_DMA_SEMAPHORE:
			pg->dma_semaphore = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_REPORT:
			pg->dma_report = parameter;
			break;

		case NV097_SET_CONTROL0: {
			bool stencil_write_enable =
				parameter & NV097_SET_CONTROL0_STENCIL_WRITE_ENABLE;
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_STENCIL_WRITE_ENABLE,
				stencil_write_enable);

			uint32_t z_format = GET_MASK(parameter, NV097_SET_CONTROL0_Z_FORMAT);
			SET_MASK(pg->regs[RI(NV_PGRAPH_SETUPRASTER)],
				NV_PGRAPH_SETUPRASTER_Z_FORMAT, z_format);

			bool z_perspective =
				parameter & NV097_SET_CONTROL0_Z_PERSPECTIVE_ENABLE;
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_Z_PERSPECTIVE_ENABLE,
				z_perspective);

			int color_space_convert =
				GET_MASK(parameter, NV097_SET_CONTROL0_COLOR_SPACE_CONVERT);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_CSCONVERT,
				color_space_convert);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_DEPTH_STENCIL]++;
			break;
		}

		case NV097_SET_FOG_MODE: {
			/* FIXME: There is also NV_PGRAPH_CSV0_D_FOG_MODE */
			unsigned int mode;
			switch (parameter) {
			case NV097_SET_FOG_MODE_V_LINEAR:
				mode = NV_PGRAPH_CONTROL_3_FOG_MODE_LINEAR; break;
			case NV097_SET_FOG_MODE_V_EXP:
				mode = NV_PGRAPH_CONTROL_3_FOG_MODE_EXP; break;
			case NV097_SET_FOG_MODE_V_EXP2:
				mode = NV_PGRAPH_CONTROL_3_FOG_MODE_EXP2; break;
			case NV097_SET_FOG_MODE_V_EXP_ABS:
				mode = NV_PGRAPH_CONTROL_3_FOG_MODE_EXP_ABS; break;
			case NV097_SET_FOG_MODE_V_EXP2_ABS:
				mode = NV_PGRAPH_CONTROL_3_FOG_MODE_EXP2_ABS; break;
			case NV097_SET_FOG_MODE_V_LINEAR_ABS:
				mode = NV_PGRAPH_CONTROL_3_FOG_MODE_LINEAR_ABS; break;
			default:
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_3)], NV_PGRAPH_CONTROL_3_FOG_MODE,
				mode);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_BLEND]++;
			break;
		}
		case NV097_SET_FOG_GEN_MODE: {
			unsigned int mode;
			switch (parameter) {
			case NV097_SET_FOG_GEN_MODE_V_SPEC_ALPHA:
				mode = NV_PGRAPH_CSV0_D_FOGGENMODE_SPEC_ALPHA; break;
			case NV097_SET_FOG_GEN_MODE_V_RADIAL:
				mode = NV_PGRAPH_CSV0_D_FOGGENMODE_RADIAL; break;
			case NV097_SET_FOG_GEN_MODE_V_PLANAR:
				mode = NV_PGRAPH_CSV0_D_FOGGENMODE_PLANAR; break;
			case NV097_SET_FOG_GEN_MODE_V_ABS_PLANAR:
				mode = NV_PGRAPH_CSV0_D_FOGGENMODE_ABS_PLANAR; break;
			case NV097_SET_FOG_GEN_MODE_V_FOG_X:
				mode = NV_PGRAPH_CSV0_D_FOGGENMODE_FOG_X; break;
			default:
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_FOGGENMODE, mode);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_BLEND]++;
			break;
		}
		case NV097_SET_FOG_COLOR: {
			/* parameter channels are ABGR, PGRAPH channels are ARGB */
			uint8_t alpha = GET_MASK(parameter, NV097_SET_FOG_COLOR_ALPHA);
			uint8_t blue = GET_MASK(parameter, NV097_SET_FOG_COLOR_BLUE);
			uint8_t green = GET_MASK(parameter, NV097_SET_FOG_COLOR_GREEN);
			uint8_t red = GET_MASK(parameter, NV097_SET_FOG_COLOR_RED);
			SET_MASK(pg->regs[RI(NV_PGRAPH_FOGCOLOR)], NV_PGRAPH_FOGCOLOR_ALPHA, alpha);
			SET_MASK(pg->regs[RI(NV_PGRAPH_FOGCOLOR)], NV_PGRAPH_FOGCOLOR_RED, red);
			SET_MASK(pg->regs[RI(NV_PGRAPH_FOGCOLOR)], NV_PGRAPH_FOGCOLOR_GREEN, green);
			SET_MASK(pg->regs[RI(NV_PGRAPH_FOGCOLOR)], NV_PGRAPH_FOGCOLOR_BLUE, blue);
			pg->dirty[NV2A_DIRTY_BLEND]++; // table already bumped dirty[NV2A_DIRTY_PGRAPH]
			break;
		}
		case NV097_SET_BLEND_FUNC_SFACTOR: {
			unsigned int factor;
			switch (parameter) {
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ZERO:
				factor = NV_PGRAPH_BLEND_SFACTOR_ZERO; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_SRC_COLOR:
				factor = NV_PGRAPH_BLEND_SFACTOR_SRC_COLOR; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE_MINUS_SRC_COLOR:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_SRC_COLOR; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_SRC_ALPHA:
				factor = NV_PGRAPH_BLEND_SFACTOR_SRC_ALPHA; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE_MINUS_SRC_ALPHA:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_SRC_ALPHA; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_DST_ALPHA:
				factor = NV_PGRAPH_BLEND_SFACTOR_DST_ALPHA; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE_MINUS_DST_ALPHA:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_DST_ALPHA; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_DST_COLOR:
				factor = NV_PGRAPH_BLEND_SFACTOR_DST_COLOR; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE_MINUS_DST_COLOR:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_DST_COLOR; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_SRC_ALPHA_SATURATE:
				factor = NV_PGRAPH_BLEND_SFACTOR_SRC_ALPHA_SATURATE; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_CONSTANT_COLOR:
				factor = NV_PGRAPH_BLEND_SFACTOR_CONSTANT_COLOR; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE_MINUS_CONSTANT_COLOR:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_CONSTANT_COLOR; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_CONSTANT_ALPHA:
				factor = NV_PGRAPH_BLEND_SFACTOR_CONSTANT_ALPHA; break;
			case NV097_SET_BLEND_FUNC_SFACTOR_V_ONE_MINUS_CONSTANT_ALPHA:
				factor = NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_CONSTANT_ALPHA; break;
			default:
				fprintf(stderr, "Unknown blend source factor: 0x%x\n", parameter);
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_BLEND)], NV_PGRAPH_BLEND_SFACTOR, factor);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_BLEND]++;
			break;
		}

		case NV097_SET_BLEND_FUNC_DFACTOR: {
			unsigned int factor;
			switch (parameter) {
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ZERO:
				factor = NV_PGRAPH_BLEND_DFACTOR_ZERO; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_SRC_COLOR:
				factor = NV_PGRAPH_BLEND_DFACTOR_SRC_COLOR; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE_MINUS_SRC_COLOR:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE_MINUS_SRC_COLOR; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_SRC_ALPHA:
				factor = NV_PGRAPH_BLEND_DFACTOR_SRC_ALPHA; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE_MINUS_SRC_ALPHA:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE_MINUS_SRC_ALPHA; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_DST_ALPHA:
				factor = NV_PGRAPH_BLEND_DFACTOR_DST_ALPHA; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE_MINUS_DST_ALPHA:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE_MINUS_DST_ALPHA; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_DST_COLOR:
				factor = NV_PGRAPH_BLEND_DFACTOR_DST_COLOR; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE_MINUS_DST_COLOR:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE_MINUS_DST_COLOR; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_SRC_ALPHA_SATURATE:
				factor = NV_PGRAPH_BLEND_DFACTOR_SRC_ALPHA_SATURATE; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_CONSTANT_COLOR:
				factor = NV_PGRAPH_BLEND_DFACTOR_CONSTANT_COLOR; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE_MINUS_CONSTANT_COLOR:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE_MINUS_CONSTANT_COLOR; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_CONSTANT_ALPHA:
				factor = NV_PGRAPH_BLEND_DFACTOR_CONSTANT_ALPHA; break;
			case NV097_SET_BLEND_FUNC_DFACTOR_V_ONE_MINUS_CONSTANT_ALPHA:
				factor = NV_PGRAPH_BLEND_DFACTOR_ONE_MINUS_CONSTANT_ALPHA; break;
			default:
				fprintf(stderr, "Unknown blend destination factor: 0x%x\n", parameter);
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_BLEND)], NV_PGRAPH_BLEND_DFACTOR, factor);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_BLEND]++;
			break;
		}

		case NV097_SET_BLEND_EQUATION: {
			unsigned int equation;
			switch (parameter) {
			case NV097_SET_BLEND_EQUATION_V_FUNC_SUBTRACT:
				equation = 0; break;
			case NV097_SET_BLEND_EQUATION_V_FUNC_REVERSE_SUBTRACT:
				equation = 1; break;
			case NV097_SET_BLEND_EQUATION_V_FUNC_ADD:
				equation = 2; break;
			case NV097_SET_BLEND_EQUATION_V_MIN:
				equation = 3; break;
			case NV097_SET_BLEND_EQUATION_V_MAX:
				equation = 4; break;
			case NV097_SET_BLEND_EQUATION_V_FUNC_REVERSE_SUBTRACT_SIGNED:
				equation = 5; break;
			case NV097_SET_BLEND_EQUATION_V_FUNC_ADD_SIGNED:
				equation = 6; break;
			default:
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_BLEND)], NV_PGRAPH_BLEND_EQN, equation);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_BLEND]++;
			break;
		}

		case NV097_SET_COLOR_MASK: {
			// No table entry for this method (NV097 param bits != PGRAPH reg bits),
			// so NV_PGRAPH_CONTROL_0 is unchanged — read current state directly.
			pg->surface_color.write_enabled_cache |= pgraph_get_color_write_enabled(pg);

			bool alpha = parameter & NV097_SET_COLOR_MASK_ALPHA_WRITE_ENABLE;
			bool red = parameter & NV097_SET_COLOR_MASK_RED_WRITE_ENABLE;
			bool green = parameter & NV097_SET_COLOR_MASK_GREEN_WRITE_ENABLE;
			bool blue = parameter & NV097_SET_COLOR_MASK_BLUE_WRITE_ENABLE;
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE, alpha);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE, red);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE, green);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_0)],
				NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE, blue);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_BLEND]++;
			break;
		}
		case NV097_SET_DEPTH_MASK:
			// old_reg is the pre-write NV_PGRAPH_CONTROL_0 value (table already wrote the new ZWRITEENABLE bit)
			pg->surface_zeta.write_enabled_cache |= (old_reg & (NV_PGRAPH_CONTROL_0_ZWRITEENABLE | NV_PGRAPH_CONTROL_0_STENCIL_WRITE_ENABLE)) != 0;
			break;
		case NV097_SET_STENCIL_OP_FAIL:
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_2)],
				NV_PGRAPH_CONTROL_2_STENCIL_OP_FAIL,
				kelvin_map_stencil_op(parameter));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_DEPTH_STENCIL]++;
			break;
		case NV097_SET_STENCIL_OP_ZFAIL:
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_2)],
				NV_PGRAPH_CONTROL_2_STENCIL_OP_ZFAIL,
				kelvin_map_stencil_op(parameter));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_DEPTH_STENCIL]++;
			break;
		case NV097_SET_STENCIL_OP_ZPASS:
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_2)],
				NV_PGRAPH_CONTROL_2_STENCIL_OP_ZPASS,
				kelvin_map_stencil_op(parameter));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_DEPTH_STENCIL]++;
			break;

		case NV097_SET_FRONT_POLYGON_MODE:
			SET_MASK(pg->regs[RI(NV_PGRAPH_SETUPRASTER)],
				NV_PGRAPH_SETUPRASTER_FRONTFACEMODE,
				kelvin_map_polygon_mode(parameter));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_RASTERIZER]++;
			break;
		case NV097_SET_BACK_POLYGON_MODE:
			SET_MASK(pg->regs[RI(NV_PGRAPH_SETUPRASTER)],
				NV_PGRAPH_SETUPRASTER_BACKFACEMODE,
				kelvin_map_polygon_mode(parameter));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_RASTERIZER]++;
			break;
		case NV097_SET_CULL_FACE: {
			unsigned int face;
			switch (parameter) {
			case NV097_SET_CULL_FACE_V_FRONT:
				face = NV_PGRAPH_SETUPRASTER_CULLCTRL_FRONT; break;
			case NV097_SET_CULL_FACE_V_BACK:
				face = NV_PGRAPH_SETUPRASTER_CULLCTRL_BACK; break;
			case NV097_SET_CULL_FACE_V_FRONT_AND_BACK:
				face = NV_PGRAPH_SETUPRASTER_CULLCTRL_FRONT_AND_BACK; break;
			default:
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_SETUPRASTER)],
				NV_PGRAPH_SETUPRASTER_CULLCTRL,
				face);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_RASTERIZER]++;
			break;
		}
		case NV097_SET_FRONT_FACE: {
			bool ccw;
			switch (parameter) {
			case NV097_SET_FRONT_FACE_V_CW:
				ccw = false; break;
			case NV097_SET_FRONT_FACE_V_CCW:
				ccw = true; break;
			default:
				fprintf(stderr, "Unknown front face: 0x%x\n", parameter);
				assert(false);
				break;
			}
			SET_MASK(pg->regs[RI(NV_PGRAPH_SETUPRASTER)],
				NV_PGRAPH_SETUPRASTER_FRONTFACE,
				ccw ? 1 : 0);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_RASTERIZER]++;
			break;
		}
		CASE_4(NV097_SET_TEXGEN_S, 16) : {
			slot = (method - NV097_SET_TEXGEN_S) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_S
				: NV_PGRAPH_CSV1_A_T0_S;
			SET_MASK(pg->regs[RI(reg)], mask, kelvin_map_texgen(parameter, 0));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_TEXTURE]++;
			break;
		}
		CASE_4(NV097_SET_TEXGEN_T, 16) : {
			slot = (method - NV097_SET_TEXGEN_T) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_T
				: NV_PGRAPH_CSV1_A_T0_T;
			SET_MASK(pg->regs[RI(reg)], mask, kelvin_map_texgen(parameter, 1));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_TEXTURE]++;
			break;
		}
		CASE_4(NV097_SET_TEXGEN_R, 16) : {
			slot = (method - NV097_SET_TEXGEN_R) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_R
				: NV_PGRAPH_CSV1_A_T0_R;
			SET_MASK(pg->regs[RI(reg)], mask, kelvin_map_texgen(parameter, 2));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_TEXTURE]++;
			break;
		}
		CASE_4(NV097_SET_TEXGEN_Q, 16) : {
			slot = (method - NV097_SET_TEXGEN_Q) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_Q
				: NV_PGRAPH_CSV1_A_T0_Q;
			SET_MASK(pg->regs[RI(reg)], mask, kelvin_map_texgen(parameter, 3));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_TEXTURE]++;
			break;
		}

		// NV097_SET_PROJECTION_MATRIX, MODEL_VIEW_MATRIX, INVERSE_MODEL_VIEW_MATRIX,
		// COMPOSITE_MATRIX, TEXTURE_MATRIX, TEXGEN_PLANE_S/T/R/Q, FOG_PLANE,
		// VIEWPORT_OFFSET, EYE_POSITION, VIEWPORT_SCALE → xfctx via dispatch table.
		// SCENE_AMBIENT_COLOR, BACK_SCENE_AMBIENT_COLOR, MATERIAL_EMISSION/ALPHA,
		// BACK_MATERIAL_EMISSION/ALPHA, EYE_DIRECTION → ltctxa via dispatch table.
		// SPECULAR_PARAMS, BACK_SPECULAR_PARAMS → ltc1 via dispatch table.
		// TRANSFORM_PROGRAM, TRANSFORM_CONSTANT → indirect SRAM via dispatch table.

		// FOG_PARAMS dual-writes: ltctxa[FOG_K] via dispatch table,
		// regs[FOGPARAM0..2] via this switch (read by HostRender + HostSync for shader fog)
		CASE_3(NV097_SET_FOG_PARAMS, 4) :
			slot = (method - NV097_SET_FOG_PARAMS) / 4;
			pg->regs[RI(NV_PGRAPH_FOGPARAM0) + slot] = parameter;
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SHADER]++;
			break;

		CASE_3(NV097_SET_VERTEX3F, 4) : {
			slot = (method - NV097_SET_VERTEX3F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_POSITION];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_POSITION);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			vertex_attribute->inline_value[3] = 1.0f;
			if (slot == 2) {
				pgraph_finish_inline_buffer_vertex(pg);
			}
			break;
		}

		// NV097_SET_BACK_LIGHT_{AMBIENT,DIFFUSE,SPECULAR}_COLOR → ltctxb via dispatch table.
		// NV097_SET_LIGHT_{AMBIENT,DIFFUSE,SPECULAR}_COLOR → ltctxb via dispatch table.
		// NV097_SET_LIGHT_LOCAL_RANGE → ltc1 via dispatch table.
		// NV097_SET_LIGHT_SPOT_FALLOFF, SPOT_DIRECTION → ltctxa via dispatch table.
		// NV097_SET_LIGHT_INFINITE_HALF_VECTOR, INFINITE_DIRECTION,
		// LOCAL_POSITION, LOCAL_ATTENUATION → light geometry via dispatch table.

		CASE_4(NV097_SET_VERTEX4F, 4): {
			slot = (method - NV097_SET_VERTEX4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_POSITION];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_POSITION);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			if (slot == 3) {
				pgraph_finish_inline_buffer_vertex(pg);
			}
			break;
		}

		CASE_3(NV097_SET_NORMAL3F, 4): {
			slot = (method - NV097_SET_NORMAL3F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_NORMAL];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_NORMAL);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			vertex_attribute->inline_value[3] = 1.0f;
			break;
		}

		CASE_4(NV097_SET_DIFFUSE_COLOR4F, 4): {
			slot = (method - NV097_SET_DIFFUSE_COLOR4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_DIFFUSE];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_DIFFUSE);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			break;
		}

		CASE_3(NV097_SET_DIFFUSE_COLOR3F, 4): {
			slot = (method - NV097_SET_DIFFUSE_COLOR3F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_DIFFUSE];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_DIFFUSE);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			vertex_attribute->inline_value[3] = 1.0f;
			break;
		}

		case NV097_SET_DIFFUSE_COLOR4UB: {
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_DIFFUSE];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_DIFFUSE);
			vertex_attribute->inline_value[0] = (parameter & 0xFF) / 255.0f;
			vertex_attribute->inline_value[1] = ((parameter >> 8) & 0xFF) / 255.0f;
			vertex_attribute->inline_value[2] = ((parameter >> 16) & 0xFF) / 255.0f;
			vertex_attribute->inline_value[3] = ((parameter >> 24) & 0xFF) / 255.0f;
			break;
		}

		CASE_4(NV097_SET_SPECULAR_COLOR4F, 4): {
			slot = (method - NV097_SET_SPECULAR_COLOR4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_SPECULAR];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_SPECULAR);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			break;
		}

		CASE_3(NV097_SET_SPECULAR_COLOR3F, 4): {
			slot = (method - NV097_SET_SPECULAR_COLOR3F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_SPECULAR];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_SPECULAR);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			vertex_attribute->inline_value[3] = 1.0f;
			break;
		}

		case NV097_SET_SPECULAR_COLOR4UB: {
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_SPECULAR];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_SPECULAR);
			vertex_attribute->inline_value[0] = (parameter & 0xFF) / 255.0f;
			vertex_attribute->inline_value[1] = ((parameter >> 8) & 0xFF) / 255.0f;
			vertex_attribute->inline_value[2] = ((parameter >> 16) & 0xFF) / 255.0f;
			vertex_attribute->inline_value[3] = ((parameter >> 24) & 0xFF) / 255.0f;
			break;
		}

		CASE_4(NV097_SET_TEXCOORD0_4F, 4): {
			slot = (method - NV097_SET_TEXCOORD0_4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_TEXTURE0];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_TEXTURE0);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			break;
		}

		CASE_4(NV097_SET_TEXCOORD1_4F, 4): {
			slot = (method - NV097_SET_TEXCOORD1_4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_TEXTURE1];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_TEXTURE1);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			break;
		}

		CASE_4(NV097_SET_TEXCOORD2_4F, 4): {
			slot = (method - NV097_SET_TEXCOORD2_4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_TEXTURE2];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_TEXTURE2);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			break;
		}

		CASE_4(NV097_SET_TEXCOORD3_4F, 4): {
			slot = (method - NV097_SET_TEXCOORD3_4F) / 4;
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_TEXTURE3];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_TEXTURE3);
			vertex_attribute->inline_value[slot] = *(float*)&parameter;
			break;
		}

		case NV097_SET_FOG1F: {
			VertexAttribute *vertex_attribute =
				&pg->vertex_attributes[NV2A_VERTEX_ATTR_FOG];
			pgraph_allocate_inline_buffer_vertices(pg, NV2A_VERTEX_ATTR_FOG);
			vertex_attribute->inline_value[0] = *(float*)&parameter;
			vertex_attribute->inline_value[1] = 0.0f;
			vertex_attribute->inline_value[2] = 0.0f;
			vertex_attribute->inline_value[3] = 1.0f;
			break;
		}

		case NV097_SET_EDGE_FLAG:
			pg->regs[RI(NV_PGRAPH_SETUPRASTER)] =
				(pg->regs[RI(NV_PGRAPH_SETUPRASTER)] & ~(1 << 30))
				| ((parameter ? 1 : 0) << 30);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_RASTERIZER]++;
			break;

		case NV097_SET_LINE_WIDTH:
			pg->line_width = *(float*)&parameter;
			break;

		case NV097_INVALIDATE_VERTEX_CACHE_FILE:
			// This is a hint to flush the post-T&L vertex cache.
			// No action required — our host renderer doesn't cache transformed vertices.
			break;

		CASE_16(NV097_SET_VERTEX_DATA_ARRAY_FORMAT, 4): {

			slot = (method - NV097_SET_VERTEX_DATA_ARRAY_FORMAT) / 4;
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];

			vertex_attribute->format =
				GET_MASK(parameter, NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE);
			vertex_attribute->count =
				GET_MASK(parameter, NV097_SET_VERTEX_DATA_ARRAY_FORMAT_SIZE);
			vertex_attribute->stride =
				GET_MASK(parameter, NV097_SET_VERTEX_DATA_ARRAY_FORMAT_STRIDE);

			pg->vertex_attributes_generation++;

			NV2A_DPRINTF("vertex data array format=%d, count=%d, stride=%d\n",
				vertex_attribute->format,
				vertex_attribute->count,
				vertex_attribute->stride);

			switch (vertex_attribute->format) {
			case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_UB_D3D:
				vertex_attribute->size = 1;
				assert(vertex_attribute->count == 4);
				break;
			case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_UB_OGL:
				vertex_attribute->size = 1;
				break;
			case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_S1:
				vertex_attribute->size = 2;
				break;
			case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_F:
				vertex_attribute->size = 4;
				break;
			case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_S32K:
				vertex_attribute->size = 2;
				break;
			case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_CMP:
				/* 3 signed, normalized components packed in 32-bits. (11,11,10) */
				vertex_attribute->size = 4;
				break;
			default:
				fprintf(stderr, "Unknown vertex type: 0x%x\n", vertex_attribute->format);
				assert(false);
				break;
			}

			break;
		}

		CASE_16(NV097_SET_VERTEX_DATA_ARRAY_OFFSET, 4): {

			slot = (method - NV097_SET_VERTEX_DATA_ARRAY_OFFSET) / 4;

			pg->vertex_attributes[slot].dma_select =
				parameter & 0x80000000;
			pg->vertex_attributes[slot].offset =
				parameter & 0x7fffffff;

			pg->vertex_attributes_generation++;

			break;
		}

		case NV097_CLEAR_REPORT_VALUE:

			/* FIXME: Does this have a value in parameter? Also does this (also?) modify
			 *        the report memory block?
			 */
			pg->zpass_pixel_count_result = 0;

			break;

		case NV097_GET_REPORT: {
			/* FIXME: This was first intended to be watchpoint-based. However,
			 *        qemu / kvm only supports virtual-address watchpoints.
			 *        This'll do for now, but accuracy and performance with other
			 *        approaches could be better
			 */
			// Collect any pending occlusion query result before reading
			if (g_pgraph_backend.zpass_collect != nullptr)
				g_pgraph_backend.zpass_collect(d);

			uint8_t type = GET_MASK(parameter, NV097_GET_REPORT_TYPE);
			assert(type == NV097_GET_REPORT_TYPE_ZPASS_PIXEL_CNT);
			hwaddr offset = GET_MASK(parameter, NV097_GET_REPORT_OFFSET);

			uint64_t timestamp = 0x0011223344556677; /* FIXME: Update timestamp?! */
			uint32_t done = 0;

			hwaddr report_dma_len;
			uint8_t *report_data = (uint8_t*)nv_dma_map(d, pg->dma_report,
														&report_dma_len);
			assert(offset < report_dma_len);
			report_data += offset;

			stq_le_p((uint64_t*)&report_data[0], timestamp);
			stl_le_p((uint32_t*)&report_data[8], pg->zpass_pixel_count_result);
			stl_le_p((uint32_t*)&report_data[12], done);

			break;
		}

		case NV097_SET_BEGIN_END: {
			uint32_t control_0 = pg->regs[RI(NV_PGRAPH_CONTROL_0)];
			uint32_t control_1 = pg->regs[RI(NV_PGRAPH_CONTROL_1)];

			bool depth_test = control_0
				& NV_PGRAPH_CONTROL_0_ZENABLE;
			bool stencil_test = control_1
				& NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE;

			if (parameter == NV097_SET_BEGIN_END_OP_END) {

				if (pg->draw_arrays_length) {
					NV2A_DPRINTF_IF(false, "Draw Arrays");
					assert(pg->inline_buffer_length == 0);
					assert(pg->inline_array_length == 0);
					assert(pg->inline_elements_length == 0);
				} else if (pg->inline_buffer_length) {
					NV2A_DPRINTF_IF(false, "Inline Buffer");
					assert(pg->draw_arrays_length == 0);
					assert(pg->inline_array_length == 0);
					assert(pg->inline_elements_length == 0);
				} else if (pg->inline_array_length) {
					NV2A_DPRINTF_IF(false, "Inline Array");
					assert(pg->draw_arrays_length == 0);
					assert(pg->inline_buffer_length == 0);
					assert(pg->inline_elements_length == 0);
				} else if (pg->inline_elements_length) {
					NV2A_DPRINTF_IF(false, "Inline Elements");
					assert(pg->draw_arrays_length == 0);
					assert(pg->inline_buffer_length == 0);
					assert(pg->inline_array_length == 0);
				} else {
					NV2A_DPRINTF_IF(true, "EMPTY NV097_SET_BEGIN_END");
					assert(false);
				}

				if (g_pgraph_backend.draw != nullptr) {
					g_pgraph_backend.draw(d);
				}

				// End occlusion query and accumulate zpass pixel count
				if (g_pgraph_backend.zpass_end != nullptr) {
					g_pgraph_backend.zpass_end(d);
				}
			} else {

				assert(parameter <= NV097_SET_BEGIN_END_OP_POLYGON);

				pg->primitive_mode = parameter;

				if (g_pgraph_backend.draw_state_update != nullptr) {
					g_pgraph_backend.draw_state_update(d);
				}

				// Begin occlusion query for zpass pixel counting
				if (pg->zpass_pixel_count_enable && g_pgraph_backend.zpass_begin != nullptr) {
					g_pgraph_backend.zpass_begin(d);
				}

				pg->inline_elements_length = 0;
				pg->inline_array_length = 0;
				pg->inline_buffer_length = 0;
				pg->draw_arrays_length = 0;
				pg->draw_arrays_max_count = 0;
			}

			pgraph_set_surface_dirty(pg, true, depth_test || stencil_test);
			break;
		}
		// NV097_SET_TEXTURE_FORMAT: fully handled by nv097_method_table
		// (NV097_REG_DIRECT_RANGE does a full 32-bit copy of the NV097
		// parameter into NV_PGRAPH_TEXFMT0+slot*4).  The Xbox D3D runtime
		// writes pTexture->Format directly as the method argument, so the
		// register value IS the Xbox Format DWORD — bit-for-bit identical.
		// No side effects needed.
		CASE_4(NV097_SET_TEXTURE_FORMAT, 64):
			break;
		CASE_4(NV097_SET_TEXTURE_PALETTE, 64): {
			slot = (method - NV097_SET_TEXTURE_PALETTE) / 64;

			bool dma_select =
				GET_MASK(parameter, NV097_SET_TEXTURE_PALETTE_CONTEXT_DMA) == 1;
			unsigned int length =
				GET_MASK(parameter, NV097_SET_TEXTURE_PALETTE_LENGTH);
			unsigned int offset =
				GET_MASK(parameter, NV097_SET_TEXTURE_PALETTE_OFFSET);

			uint32_t *reg = &pg->regs[RI(NV_PGRAPH_TEXPALETTE0 + slot * 4)];
			SET_MASK(*reg, NV_PGRAPH_TEXPALETTE0_CONTEXT_DMA, dma_select);
			SET_MASK(*reg, NV_PGRAPH_TEXPALETTE0_LENGTH, length);
			SET_MASK(*reg, NV_PGRAPH_TEXPALETTE0_OFFSET, offset);
			pg->dirty[NV2A_DIRTY_TEXTURE]++; // table already bumped dirty[NV2A_DIRTY_PGRAPH]
			// Also wrote: pg->texture_dirty[slot] = true; (field deleted)
			break;
		}

		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x0, 64):
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x4, 64):
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x8, 64):
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0xc, 64):
			slot = (method - NV097_SET_TEXTURE_SET_BUMP_ENV_MAT) / 4;
			assert((slot / 16) > 0); // Stage 0 has no bump env
			// Also wrote: pg->bump_env_matrix[slot/16 - 1][slot%4] = *(float*)&parameter; (field deleted)
			break;

		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_SCALE, 64):
			slot = (method - NV097_SET_TEXTURE_SET_BUMP_ENV_SCALE) / 64;
			assert(slot > 0);
			break;
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_OFFSET, 64):
			slot = (method - NV097_SET_TEXTURE_SET_BUMP_ENV_OFFSET) / 64;
			assert(slot > 0);
			break;

		case NV097_ARRAY_ELEMENT16:
			//LOG_TEST_CASE("NV2A_VB_ELEMENT_U16");	
			// Test-case : Turok (in main menu)	
			// Test-case : Hunter Redeemer	
			// Test-case : Otogi (see https://github.com/Cxbx-Reloaded/Cxbx-Reloaded/pull/1113#issuecomment-385593814)
			assert(pg->inline_elements_length < NV2A_MAX_BATCH_LENGTH);
			pg->inline_elements[
				pg->inline_elements_length++] = parameter & 0xFFFF;
			pg->inline_elements[
				pg->inline_elements_length++] = parameter >> 16;
			break;
		case NV097_ARRAY_ELEMENT32:
			//LOG_TEST_CASE("NV2A_VB_ELEMENT_U32");	
			// Test-case : Turok (in main menu)
			assert(pg->inline_elements_length < NV2A_MAX_BATCH_LENGTH);
			pg->inline_elements[
				pg->inline_elements_length++] = parameter;
			break;
		case NV097_DRAW_ARRAYS: {

			unsigned int start = GET_MASK(parameter, NV097_DRAW_ARRAYS_START_INDEX);
			unsigned int count = GET_MASK(parameter, NV097_DRAW_ARRAYS_COUNT)+1;

			pg->draw_arrays_max_count = MAX(pg->draw_arrays_max_count, start + count);

			assert(pg->draw_arrays_length < ARRAY_SIZE(pg->draw_arrays_start));

			/* Attempt to connect contiguous primitives */
			if (!pg->draw_arrays_prevent_connect && pg->draw_arrays_length > 0) {
				unsigned int last_start =
					pg->draw_arrays_start[pg->draw_arrays_length - 1];
				int32_t* last_count =
					&pg->draw_arrays_count[pg->draw_arrays_length - 1];
				if (start == (last_start + *last_count)) {
					*last_count += count;
					break;
				}
			}

			pg->draw_arrays_start[pg->draw_arrays_length] = start;
			pg->draw_arrays_count[pg->draw_arrays_length] = count;
			pg->draw_arrays_length++;
			pg->draw_arrays_prevent_connect = false;
			break;
		}
		case NV097_INLINE_ARRAY:
			assert(pg->inline_array_length < NV2A_MAX_BATCH_LENGTH);
			pg->inline_array[
				pg->inline_array_length++] = parameter;
			break;
		CASE_32(NV097_SET_VERTEX_DATA2F_M, 4): {
			slot = (method - NV097_SET_VERTEX_DATA2F_M) / 4;
			unsigned int part = slot % 2;
			slot /= 2;
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
			pgraph_allocate_inline_buffer_vertices(pg, slot);
			vertex_attribute->inline_value[part] = *(float*)&parameter;
			/* FIXME: Should these really be set to 0.0 and 1.0 ? Conditions? */
			vertex_attribute->inline_value[2] = 0.0f;
			vertex_attribute->inline_value[3] = 1.0f;
			if ((slot == 0) && (part == 1)) {
				pgraph_finish_inline_buffer_vertex(pg);
			}
			break;
		}
		CASE_64(NV097_SET_VERTEX_DATA4F_M, 4): {
			slot = (method - NV097_SET_VERTEX_DATA4F_M) / 4;
			unsigned int part = slot % 4;
			slot /= 4;
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
			pgraph_allocate_inline_buffer_vertices(pg, slot);
			vertex_attribute->inline_value[part] = *(float*)&parameter;
			if ((slot == 0) && (part == 3)) {
				pgraph_finish_inline_buffer_vertex(pg);
			}
			break;
		}
		CASE_16(NV097_SET_VERTEX_DATA2S, 4): {
			slot = (method - NV097_SET_VERTEX_DATA2S) / 4;
			assert(false); /* FIXME: Untested! */
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
			pgraph_allocate_inline_buffer_vertices(pg, slot);
			vertex_attribute->inline_value[0] = (float)(int16_t)(parameter & 0xFFFF);
			vertex_attribute->inline_value[1] = (float)(int16_t)(parameter >> 16);
			vertex_attribute->inline_value[2] = 0.0f;
			vertex_attribute->inline_value[3] = 1.0f;
			if (slot == 0) {
				pgraph_finish_inline_buffer_vertex(pg);
			}
			break;
		}
		CASE_16(NV097_SET_VERTEX_DATA4UB, 4) : {
			slot = (method - NV097_SET_VERTEX_DATA4UB) / 4;
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
			pgraph_allocate_inline_buffer_vertices(pg, slot);
			vertex_attribute->inline_value[0] = (parameter & 0xFF) / 255.0f;
			vertex_attribute->inline_value[1] = ((parameter >> 8) & 0xFF) / 255.0f;
			vertex_attribute->inline_value[2] = ((parameter >> 16) & 0xFF) / 255.0f;
			vertex_attribute->inline_value[3] = ((parameter >> 24) & 0xFF) / 255.0f;
			if (slot == 0) {
				pgraph_finish_inline_buffer_vertex(pg);
				assert(false); /* FIXME: Untested */
			}
			break;
		}
		CASE_32(NV097_SET_VERTEX_DATA4S_M, 4) : {
			slot = (method - NV097_SET_VERTEX_DATA4S_M) / 4;
			unsigned int part = slot % 2;
			slot /= 2;
			assert(false); /* FIXME: Untested! */
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
			pgraph_allocate_inline_buffer_vertices(pg, slot);
			/* FIXME: Is mapping to [-1,+1] correct? */
			vertex_attribute->inline_value[part * 2 + 0] = ((int16_t)(parameter & 0xFFFF)
														 * 2.0f + 1) / 65535.0f;
			vertex_attribute->inline_value[part * 2 + 1] = ((int16_t)(parameter >> 16)
														 * 2.0f + 1) / 65535.0f;
			if ((slot == 0) && (part == 1)) {
				pgraph_finish_inline_buffer_vertex(pg);
			}
			break;
		}
		case NV097_BACK_END_WRITE_SEMAPHORE_RELEASE: {
			//qemu_mutex_unlock(&pg->pgraph_lock);
			//qemu_mutex_lock_iothread();

			uint32_t semaphore_offset = pg->regs[RI(NV_PGRAPH_SEMAPHOREOFFSET)];

			xbox::addr_xt semaphore_dma_len;
			uint8_t *semaphore_data = (uint8_t*)nv_dma_map(d, pg->dma_semaphore,
				&semaphore_dma_len);
			if (semaphore_offset >= semaphore_dma_len) {
				EmuLog(LOG_LEVEL::WARNING, "Semaphore offset 0x%X >= dma_len 0x%X, skipping release",
					semaphore_offset, semaphore_dma_len);
				break;
			}
			semaphore_data += semaphore_offset;

			stl_le_p((uint32_t*)semaphore_data, parameter);

			//qemu_mutex_lock(&pg->pgraph_lock);
			//qemu_mutex_unlock_iothread();

			break;
		}
		case NV097_CLEAR_SURFACE: {
			pg->clear_surface_flags = parameter;
			if (g_pgraph_backend.draw_clear != nullptr) {
				g_pgraph_backend.draw_clear(d);
			}
			break;
		}

		case NV097_SET_SHADOW_ZSLOPE_THRESHOLD:
			assert(parameter == 0x7F800000); /* FIXME: Unimplemented */
			break;

		case NV097_SET_TRANSFORM_EXECUTION_MODE:
			// Test-case : Whiplash
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_MODE,
				GET_MASK(parameter,
					NV097_SET_TRANSFORM_EXECUTION_MODE_MODE));
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_RANGE_MODE,
				GET_MASK(parameter,
					NV097_SET_TRANSFORM_EXECUTION_MODE_RANGE_MODE));
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SHADER]++;
			break;
		case NV097_SET_TRANSFORM_PROGRAM_LOAD:
			assert(parameter < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
			break;
		case NV097_SET_TRANSFORM_PROGRAM_START:
			assert(parameter < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
			break;
		case NV097_SET_TRANSFORM_CONSTANT_LOAD:
			assert(parameter < NV2A_VERTEXSHADER_CONSTANTS);
			NV2A_DPRINTF("load to %d\n", parameter);
			break;

		CASE_4(NV097_SET_TRANSFORM_DATA, 4): {
			// Stores input v0 components for the next LAUNCH_TRANSFORM_PROGRAM
			slot = (method - NV097_SET_TRANSFORM_DATA) / 4;
			pg->xf.vertex_state_shader_v0[slot] = parameter;
			break;
		}

		case NV097_LAUNCH_TRANSFORM_PROGRAM: {
			// Execute a vertex state shader (XSS) at the given program address.
			// Input v0 was set by prior SET_TRANSFORM_DATA writes; output updates
			// vsh_constants (transform context RAM) in-place.
			unsigned int program_start = parameter;
			assert(program_start < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
			if (g_pgraph_backend.launch_transform_program != nullptr) {
				g_pgraph_backend.launch_transform_program(d, program_start);
			}
			break;
		}

		case NV097_SET_PROVOKING_VERTEX: 
			assert(parameter <= 1);
			// Handled by method table: NV_PGRAPH_CONTROL_3_PROVOKING_VERTEX
			break;

		case NV097_SET_COLOR_MATERIAL: {
			// Compound write: 4 material source fields packed in parameter
			// bits 0-1: emission, 2-3: ambient, 4-5: diffuse, 6-7: specular
			// Values: 0=FROM_MATERIAL, 1=FROM_COLOR1, 2=FROM_COLOR2
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_EMISSION, (parameter >> 0) & 3);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_AMBIENT,  (parameter >> 2) & 3);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_DIFFUSE,  (parameter >> 4) & 3);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_SPECULAR, (parameter >> 6) & 3);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SHADER]++;
			break;
		}

		case NV097_SET_TWO_SIDED_LIGHT_EN:
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_TWO_SIDE_LIGHTING, parameter ? 1 : 0);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SHADER]++;
			break;

		case NV097_SET_POINT_PARAMS_ENABLE:
			// Writes to BOTH CSV0_D and CONTROL_3
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_POINTPARAMSENABLE, parameter ? 1 : 0);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_3)], NV_PGRAPH_CONTROL_3_POINTPARAMSENABLE, parameter ? 1 : 0);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SHADER]++;
			break;

		case NV097_SET_LIGHT_CONTROL: {
			// Extracts 3 fields into CSV0_C
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_SEPARATE_SPECULAR,
				(parameter & NV097_SET_LIGHT_CONTROL_SEPARATE_SPECULAR) ? 1 : 0);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_LOCALEYE,
				(parameter & NV097_SET_LIGHT_CONTROL_LOCALEYE) ? 1 : 0);
			SET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)], NV_PGRAPH_CSV0_C_ALPHA_FROM_MATERIAL_SPECULAR,
				(parameter & NV097_SET_LIGHT_CONTROL_ALPHA_FROM_MATERIAL_SPECULAR) ? 1 : 0);
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_SHADER]++;
			break;
		}

		CASE_8(NV097_SET_POINT_PARAMS, 4): {
			// 8 float point attenuation parameters
			slot = (method - NV097_SET_POINT_PARAMS) / 4;
			pg->point_params[slot] = *(float*)&parameter;
			break;
		}

		case NV097_SET_POINT_SIZE:
			if (parameter > NV097_SET_POINT_SIZE_V_MAX) {
				break;
			}
			pg->regs[RI(NV_PGRAPH_POINTSIZE)] = parameter;
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
			pg->dirty[NV2A_DIRTY_RASTERIZER]++;
			break;

		// TODO: Implement these methods (not table-compatible due to value remapping or multi-reg writes).
		// See xemu pgraph.c for reference implementations.
		//
		// case NV097_SET_SHADE_MODE:
		//     Value remapping: V_FLAT(0x1D00) -> SHADEMODE_FLAT(0),
		//     V_SMOOTH(0x1D01) -> SHADEMODE_SMOOTH(1)
		//     Target: NV_PGRAPH_CONTROL_3_SHADEMODE
		//     break;
		//
		// case NV097_SET_ZMIN_MAX_CONTROL:
		//     Extracts ZCLAMP_EN field, maps CULL->0, CLAMP->1
		//     Target: NV_PGRAPH_ZCOMPRESSOCCLUDE_ZCLAMP_EN
		//     break;

		// TODO: These cases wrote to PGRAPHState fields that have since been deleted.
		// The register writes are handled by method table entries; the struct field
		// writes may need to be restored once replacement PGRAPH register mappings
		// are identified. See xemu pgraph.c for reference implementations.
		//
		// case NV097_SET_CONTEXT_DMA_NOTIFIES:
		//     pg->dma_notifies = parameter;
		//     break;
		case NV097_SET_CONTEXT_DMA_A:
		    pg->dma_a = parameter;
		    break;
		case NV097_SET_CONTEXT_DMA_B:
		    pg->dma_b = parameter;
		    break;
		// case NV097_SET_CONTEXT_DMA_STATE:
		//     pg->dma_state = parameter;
		//     break;
		// case NV097_SET_CONTEXT_DMA_ZETA:
		//     pg->dma_zeta = parameter;
		//     break;
		// case NV097_SET_CONTEXT_DMA_VERTEX_A:
		// case NV097_SET_CONTEXT_DMA_VERTEX_B:
		//     DMA context methods are not needed under HLE: vertex attribute
		//     offsets already store physical byte offsets into contiguous memory,
		//     and both DMA contexts point to the same region.
		//     break;
		//
		CASE_4(NV097_SET_TEXTURE_MATRIX_ENABLE, 4):
		    slot = (method - NV097_SET_TEXTURE_MATRIX_ENABLE) / 4;
		    pg->texture_matrix_enable[slot] = parameter != 0;
		    break;

		case NV097_SET_ZPASS_PIXEL_COUNT_ENABLE:
		    pg->zpass_pixel_count_enable = parameter;
		    if (parameter) {
		        if (g_pgraph_backend.zpass_begin != nullptr)
		            g_pgraph_backend.zpass_begin(d);
		    } else {
		        if (g_pgraph_backend.zpass_end != nullptr)
		            g_pgraph_backend.zpass_end(d);
		    }
		    break;
		//
		// CASE_4(NV097_SET_TEXTURE_OFFSET, 64):
		//     Handled by table: NV_PGRAPH_TEXOFFSET0
		//     Also wrote: pg->texture_dirty[slot] = true;
		//     break;
		// CASE_4(NV097_SET_TEXTURE_IMAGE_RECT, 64):
		//     Handled by table: NV_PGRAPH_TEXIMAGERECT0
		//     Also wrote: pg->texture_dirty[slot] = true;
		//     break;
		//
		// case NV097_SET_TRANSFORM_PROGRAM_CXT_WRITE_EN:
		//     // Test-case: Whiplash
		//     pg->enable_vertex_program_write = parameter;
		//     break;

		// ===== Hardware Tessellation (Patch) Methods =====
		case NV097_SET_BEGIN_PATCH0:
			pg->patch.patch0 = parameter;
			pg->patch.active = true;
			pg->patch.curveCount = 0;
			pg->patch.totalCoeffs = 0;
			pg->patch.currentCurveAttr = -1;
			break;
		case NV097_SET_BEGIN_PATCH1:
			pg->patch.patch1 = parameter;
			break;
		case NV097_SET_BEGIN_PATCH2:
			pg->patch.patch2 = parameter;
			break;
		case NV097_SET_BEGIN_PATCH3:
			pg->patch.patch3 = parameter;
			break;

		case NV097_SET_BEGIN_END_SWATCH:
			if (parameter != 0) {
				pg->patch.swatch = parameter; // Store begin format
			} else {
				// End swatch - draw the accumulated curves and reset for next swatch
				if (g_pgraph_backend.draw_patch != nullptr && pg->patch.active && pg->patch.curveCount > 0) {
					g_pgraph_backend.draw_patch(d);
				}
				// Reset curves for next swatch (keep patch0-3 and active)
				pg->patch.curveCount = 0;
				pg->patch.totalCoeffs = 0;
				pg->patch.currentCurveAttr = -1;
			}
			break;

		case NV097_SET_BEGIN_END_CURVE:
			if (parameter == 0) {
				// End current curve (END_CURVE_DATA)
				if (pg->patch.currentCurveAttr >= 0 && pg->patch.curveCount < NV2A_PATCH_MAX_CURVES) {
					PatchCurve &curve = pg->patch.curves[pg->patch.curveCount];
					curve.curveType = pg->patch.currentCurveAttr;
					curve.coeffCount = pg->patch.totalCoeffs - curve.coeffStart;
					pg->patch.curveCount++;
				}
				pg->patch.currentCurveAttr = -1;
			} else {
				// Begin curve of given type (1=STRIP, 2=LEFT_GUARD, 3=RIGHT_GUARD, etc.)
				pg->patch.currentCurveAttr = (int)parameter;
				if (pg->patch.curveCount < NV2A_PATCH_MAX_CURVES) {
					pg->patch.curves[pg->patch.curveCount].coeffStart = pg->patch.totalCoeffs;
				}
			}
			break;

		CASE_4(NV097_SET_CURVE_COEFFICIENTS, 4): {
			int idx = (method - NV097_SET_CURVE_COEFFICIENTS) / 4;
			if (pg->patch.totalCoeffs < NV2A_PATCH_MAX_COEFFS) {
				int base = pg->patch.totalCoeffs * 4 + idx;
				uint32_t u = parameter;
				float f;
				memcpy(&f, &u, sizeof(f));
				pg->patch.coefficients[base] = f;
				// Advance total count after writing the 4th component
				if (idx == 3) {
					pg->patch.totalCoeffs++;
				}
			}
			break;
		}

		case NV097_SET_END_PATCH:
			// Finalize any open curve
			if (pg->patch.currentCurveAttr >= 0 && pg->patch.curveCount < NV2A_PATCH_MAX_CURVES) {
				PatchCurve &curve = pg->patch.curves[pg->patch.curveCount];
				curve.curveType = pg->patch.currentCurveAttr;
				curve.coeffCount = pg->patch.totalCoeffs - curve.coeffStart;
				pg->patch.curveCount++;
				pg->patch.currentCurveAttr = -1;
			}
			// Dispatch tessellation for any remaining curves (if swatch didn't already draw them)
			if (g_pgraph_backend.draw_patch != nullptr && pg->patch.active && pg->patch.curveCount > 0) {
				g_pgraph_backend.draw_patch(d);
			}
			pg->patch.active = false;
			break;

		default:
			NV2A_DPRINTF_IF(true, "    unhandled  (0x%02x 0x%08x)",
					graphics_class, method);
			break;
		}
		break;
	}

	default:
		NV2A_DPRINTF_IF(true, "Unknown Graphics Class/Method 0x%08X/0x%08X",
						graphics_class, method);
		break;
	}

}

static void pgraph_switch_context(NV2AState *d, unsigned int channel_id)
{
    bool channel_valid =
        d->pgraph.regs[RI(NV_PGRAPH_CTX_CONTROL)] & NV_PGRAPH_CTX_CONTROL_CHID;
    unsigned pgraph_channel_id = GET_MASK(d->pgraph.regs[RI(NV_PGRAPH_CTX_USER)], NV_PGRAPH_CTX_USER_CHID);
	// Cxbx Note : This isn't present in xqemu / OpenXbox : d->pgraph.pgraph_lock.lock();
    bool valid = channel_valid && pgraph_channel_id == channel_id;
	if (!valid) {
        SET_MASK(d->pgraph.regs[RI(NV_PGRAPH_TRAPPED_ADDR)],
                 NV_PGRAPH_TRAPPED_ADDR_CHID, channel_id);

        NV2A_DPRINTF("pgraph switching to ch %d\n", channel_id);

        /* TODO: hardware context switching */
        assert(!(d->pgraph.regs[RI(NV_PGRAPH_DEBUG_3)]
                & NV_PGRAPH_DEBUG_3_HW_CONTEXT_SWITCH));

		qemu_mutex_unlock(&d->pgraph.pgraph_lock);
		qemu_mutex_lock_iothread();
		d->pgraph.pending_interrupts |= NV_PGRAPH_INTR_CONTEXT_SWITCH; // TODO : Should this be done before unlocking pgraph_lock?
		update_irq(d);

		qemu_mutex_lock(&d->pgraph.pgraph_lock);
		qemu_mutex_unlock_iothread();

        // wait for the interrupt to be serviced
		while (d->pgraph.pending_interrupts & NV_PGRAPH_INTR_CONTEXT_SWITCH) {
			qemu_cond_wait(&d->pgraph.interrupt_cond, &d->pgraph.pgraph_lock);
		}
	}
}

static void pgraph_wait_fifo_access(NV2AState *d) {
    while (!(d->pgraph.regs[RI(NV_PGRAPH_FIFO)] & NV_PGRAPH_FIFO_ACCESS)) {
		qemu_cond_wait(&d->pgraph.fifo_access_cond, &d->pgraph.pgraph_lock);
	}
}

static void pgraph_log_method(unsigned int subchannel,
								unsigned int graphics_class,
								unsigned int method, uint32_t parameter) {
	static unsigned int last = 0;
	static unsigned int count = 0;

	extern const char *NV2AMethodToString(DWORD dwMethod); // implemented in PushBuffer.cpp

	if (last == 0x1800 && method != last) {
		const char* method_name = NV2AMethodToString(last); // = 'NV2A_VB_ELEMENT_U16'
		NV2A_DPRINTF_IF(true, "d->pgraph method (%d) 0x%08X %s * %d",
						subchannel, last, method_name, count);
	}
	if (method != 0x1800) {
		// const char* method_name = NV2AMethodToString(method);
		// unsigned int nmethod = 0;
		// switch (graphics_class) {
		// case NV_KELVIN_PRIMITIVE:
		// 	nmethod = method | (0x5c << 16);
		// 	break;
		// case NV_CONTEXT_SURFACES_2D:
		// 	nmethod = method | (0x6d << 16);
		// 	break;
        // case NV_CONTEXT_PATTERN:
        // 	nmethod = method | (0x68 << 16);
        // 	break;
		// default:
		// 	break;
		// }
		// if (method_name) {
		// 	NV2A_DPRINTF("d->pgraph method (%d): %s (0x%x)\n",
		// 		subchannel, method_name, parameter);
		// } else {
			NV2A_DPRINTF("pgraph method (%d): 0x%08X -> 0x%04x (0x%x)\n",
				subchannel, graphics_class, method, parameter);
		// }

	}
	if (method == last) { count++; }
	else { count = 0; }
	last = method;
}

static void pgraph_allocate_inline_buffer_vertices(PGRAPHState *pg,
                                                   unsigned int attr)
{
    unsigned int i;
    VertexAttribute *vertex_attribute = &pg->vertex_attributes[attr];

    if (vertex_attribute->inline_buffer || pg->inline_buffer_length == 0) {
        return;
    }

    /* Reuse persistent pool allocation to avoid malloc/free per draw */
    if (!vertex_attribute->inline_buffer_pool) {
        vertex_attribute->inline_buffer_pool = (float*)g_malloc(NV2A_MAX_BATCH_LENGTH
                                                      * sizeof(float) * 4);
    }
    vertex_attribute->inline_buffer = vertex_attribute->inline_buffer_pool;

    /* Now upload the previous vertex_attribute value */
    for (i = 0; i < pg->inline_buffer_length; i++) {
        memcpy(&vertex_attribute->inline_buffer[i * 4],
               vertex_attribute->inline_value,
               sizeof(float) * 4);
    }
}

static void pgraph_finish_inline_buffer_vertex(PGRAPHState *pg)
{
	unsigned int i;

    assert(pg->inline_buffer_length < NV2A_MAX_BATCH_LENGTH);

    for (i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
        VertexAttribute *vertex_attribute = &pg->vertex_attributes[i];
        if (vertex_attribute->inline_buffer) {
            memcpy(&vertex_attribute->inline_buffer[
                      pg->inline_buffer_length * 4],
                   vertex_attribute->inline_value,
                   sizeof(float) * 4);
        }
    }

    pg->inline_buffer_length++;
}

void pgraph_init(NV2AState *d)
{
    int i;

    PGRAPHState *pg = &d->pgraph;

	nv097_init_method_table();

	qemu_mutex_init(&pg->pgraph_lock);
	qemu_cond_init(&pg->interrupt_cond);
	qemu_cond_init(&pg->fifo_access_cond);
	qemu_cond_init(&pg->flip_3d);

	// Initialize vertex attribute defaults (inline_value / sticky registers).
	// On NV2A hardware reset, diffuse (slot 3) and specular (slot 4) default
	// to white (1,1,1,1).  All others default to (0,0,0,1).  This matches
	// D3D semantics: unset vertex colors are opaque white, positions/texcoords
	// are zero with w=1.
	for (i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
		pg->vertex_attributes[i].inline_value[0] = 0.0f;
		pg->vertex_attributes[i].inline_value[1] = 0.0f;
		pg->vertex_attributes[i].inline_value[2] = 0.0f;
		pg->vertex_attributes[i].inline_value[3] = 1.0f;
	}
	// Diffuse and specular default to opaque white
	pg->vertex_attributes[3].inline_value[0] = 1.0f;
	pg->vertex_attributes[3].inline_value[1] = 1.0f;
	pg->vertex_attributes[3].inline_value[2] = 1.0f;
	pg->vertex_attributes[4].inline_value[0] = 1.0f;
	pg->vertex_attributes[4].inline_value[1] = 1.0f;
	pg->vertex_attributes[4].inline_value[2] = 1.0f;

	// In HLE mode, the puller processes pushbuffer methods into PGRAPH
	// register state (so RC/VS interpreters can read from it) but does
	// not render.  Initialize minimum context so pgraph_handle_method()
	// can execute without tripping assertions:
	// - CTX_CONTROL CHID bit: marks channel as valid
	// - CTX_SWITCH1 GRCLASS: NV_KELVIN_PRIMITIVE (0x97) — Xbox 3D class
	// - CTX_CACHE1: standard Xbox subchannel → object class mapping
		pg->regs[RI(NV_PGRAPH_CTX_CONTROL)] = NV_PGRAPH_CTX_CONTROL_CHID;
		pg->regs[RI(NV_PGRAPH_CTX_SWITCH1)] = NV_KELVIN_PRIMITIVE;
		// Standard Xbox D3D subchannel assignments:
		pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + 0 * 4)] = NV_KELVIN_PRIMITIVE;        // SC 0: 3D (Kelvin)
		pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + 1 * 4)] = NV_CONTEXT_PATTERN;         // SC 1: Pattern
		pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + 2 * 4)] = NV_CONTEXT_SURFACES_2D;     // SC 2: Surfaces2D
		pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + 3 * 4)] = NV_IMAGE_BLIT;              // SC 3: ImageBlit
		pg->regs[RI(NV_PGRAPH_CTX_CACHE1 + 4 * 4)] = NV_MEMORY_TO_MEMORY_FORMAT; // SC 4: MemToMem
}

void pgraph_destroy(PGRAPHState *pg)
{
	qemu_mutex_destroy(&pg->pgraph_lock);
	qemu_cond_destroy(&pg->interrupt_cond);
	qemu_cond_destroy(&pg->fifo_access_cond);
	qemu_cond_destroy(&pg->flip_3d);
}

static bool pgraph_get_color_write_enabled(PGRAPHState *pg)
{
	return pg->regs[RI(NV_PGRAPH_CONTROL_0)] & (
		NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE
		| NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE
		| NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE
		| NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE);
}

static bool pgraph_get_zeta_write_enabled(PGRAPHState *pg)
{
	return pg->regs[RI(NV_PGRAPH_CONTROL_0)] & (
		NV_PGRAPH_CONTROL_0_ZWRITEENABLE
		| NV_PGRAPH_CONTROL_0_STENCIL_WRITE_ENABLE);
}

static void pgraph_set_surface_dirty(PGRAPHState *pg, bool color, bool zeta)
{
    NV2A_DPRINTF("pgraph_set_surface_dirty(%d, %d) -- %d %d\n",
                 color, zeta,
                 pgraph_get_color_write_enabled(pg), pgraph_get_zeta_write_enabled(pg));
    /* FIXME: Does this apply to CLEARs too? */
    color = color && pgraph_get_color_write_enabled(pg);
    zeta = zeta && pgraph_get_zeta_write_enabled(pg);
    pg->surface_color.draw_dirty |= color;
    pg->surface_zeta.draw_dirty |= zeta;
}

static void pgraph_apply_anti_aliasing_factor(PGRAPHState *pg,
                                              unsigned int *width,
                                              unsigned int *height)
{
    auto surf = NV2AGetSurfaceState(pg);
    switch (surf.antiAliasing) {
    case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_1:
        break;
    case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_CORNER_2:
        if (width) { *width *= 2; }
        break;
    case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_SQUARE_OFFSET_4:
        if (width) { *width *= 2; }
        if (height) { *height *= 2; }
        break;
    default:
        assert(false);
        break;
    }
}

static void pgraph_get_surface_dimensions(PGRAPHState *pg,
                                          unsigned int *width,
                                          unsigned int *height)
{
    auto surf = NV2AGetSurfaceState(pg);
    bool swizzle = (surf.surfaceType == NV097_SET_SURFACE_FORMAT_TYPE_SWIZZLE);
    if (swizzle) {
        *width = 1 << surf.logWidth;
        *height = 1 << surf.logHeight;
    } else {
        *width = surf.clipWidth;
        *height = surf.clipHeight;
    }
}

/* 16 bit to [0.0, F16_MAX = 511.9375] */
static float convert_f16_to_float(uint16_t f16) {
    if (f16 == 0x0000) { return 0.0f; }
    uint32_t i = (f16 << 11) + 0x3C000000;
    return *(float*)&i;
}

/* 24 bit to [0.0, F24_MAX] */
static float convert_f24_to_float(uint32_t f24) {
    assert(!(f24 >> 24));
    f24 &= 0xFFFFFF;
    if (f24 == 0x000000) { return 0.0f; }
    uint32_t i = f24 << 7;
    return *(float*)&i;
}

extern void __R6G5B5ToARGBRow_C(const uint8_t* src_r6g5b5, uint8_t* dst_argb, int width);
extern void ____YUY2ToARGBRow_C(const uint8_t* src_yuy2, uint8_t* rgb_buf, int width);
extern void ____UYVYToARGBRow_C(const uint8_t* src_uyvy, uint8_t* rgb_buf, int width);

/* 'converted_format' indicates the format that results when convert_texture_data() returns non-NULL converted_data. */
static const int converted_format = NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8R8G8B8;

static uint8_t* convert_texture_data(const unsigned int color_format,
                                     const uint8_t *data,
                                     const uint8_t *palette_data,
                                     const unsigned int width,
                                     const unsigned int height,
                                     const unsigned int depth,
                                     const unsigned int row_pitch,
                                     const unsigned int slice_pitch)
{
	// Note : Unswizzle is already done when entering here
	switch (color_format) {
	case NV097_SET_TEXTURE_FORMAT_COLOR_SZ_I8_A8R8G8B8: {
		// Test-case : WWE RAW2
		assert(depth == 1); /* FIXME */
		uint8_t* converted_data = (uint8_t*)g_malloc(width * height * 4);
		unsigned int x, y;
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				uint8_t index = data[y * row_pitch + x];
				uint32_t color = *(uint32_t*)(palette_data + index * 4);
				*(uint32_t*)(converted_data + y * width * 4 + x * 4) = color;
			}
		}
		return converted_data;
	}
	case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_X7SY9: {
		assert(false); /* FIXME */
		return NULL;
	}
	case NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_CR8YB8CB8YA8: {
		// Test-case : WWE RAW2
		assert(depth == 1); /* FIXME */
		uint8_t* converted_data = (uint8_t*)g_malloc(width * height * 4);
		unsigned int y;
		for (y = 0; y < height; y++) {
			const uint8_t* line = &data[y * width * 2];
			uint8_t* pixel = &converted_data[(y * width) * 4];
			____YUY2ToARGBRow_C(line, pixel, width);
			// Note : LC_IMAGE_CR8YB8CB8YA8 suggests UYVY format,
			// but for an unknown reason, the actual encoding is YUY2
		}
		return converted_data;
	}	
	case NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_YB8CR8YA8CB8: {
		assert(depth == 1); /* FIXME */
		uint8_t* converted_data = (uint8_t*)g_malloc(width * height * 4);
		unsigned int y;
		for (y = 0; y < height; y++) {
			const uint8_t* line = &data[y * width * 2];
			uint8_t* pixel = &converted_data[(y * width) * 4];
			____UYVYToARGBRow_C(line, pixel, width); // TODO : Validate LC_IMAGE_YB8CR8YA8CB8 indeed requires ____UYVYToARGBRow_C()
		}
		return converted_data;
	}
	case NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_A4V6YB6A4U6YA6: {
		assert(false); /* FIXME */
		return NULL;
	}
	case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8CR8CB8Y8: {
		assert(false); /* FIXME */
		return NULL;
	}
	case NV097_SET_TEXTURE_FORMAT_COLOR_SZ_R6G5B5:
	case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_R6G5B5: {
		assert(depth == 1); /* FIXME */
		uint8_t *converted_data = (uint8_t*)g_malloc(width * height * 4);
		unsigned int y;
		for (y = 0; y < height; y++) {
			uint16_t rgb655 = *(uint16_t*)(data + y * row_pitch);
			int8_t *pixel = (int8_t*)&converted_data[(y * width) * 4];
			__R6G5B5ToARGBRow_C((const uint8_t*)rgb655, (uint8_t*)pixel, width);
		}
		return converted_data;
	}
	default:
        return NULL;
    }
}

static unsigned int kelvin_map_stencil_op(uint32_t parameter)
{
	unsigned int op;
	switch (parameter) {
	case NV097_SET_STENCIL_OP_V_KEEP:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_KEEP; break;
	case NV097_SET_STENCIL_OP_V_ZERO:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_ZERO; break;
	case NV097_SET_STENCIL_OP_V_REPLACE:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_REPLACE; break;
	case NV097_SET_STENCIL_OP_V_INCRSAT:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_INCRSAT; break;
	case NV097_SET_STENCIL_OP_V_DECRSAT:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_DECRSAT; break;
	case NV097_SET_STENCIL_OP_V_INVERT:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_INVERT; break;
	case NV097_SET_STENCIL_OP_V_INCR:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_INCR; break;
	case NV097_SET_STENCIL_OP_V_DECR:
		op = NV_PGRAPH_CONTROL_2_STENCIL_OP_V_DECR; break;
	default:
		assert(false);
		break;
	}
	return op;
}

static unsigned int kelvin_map_polygon_mode(uint32_t parameter)
{
	unsigned int mode;
	switch (parameter) {
	case NV097_SET_FRONT_POLYGON_MODE_V_POINT:
		mode = NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_POINT; break;
	case NV097_SET_FRONT_POLYGON_MODE_V_LINE:
		mode = NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_LINE; break;
	case NV097_SET_FRONT_POLYGON_MODE_V_FILL:
		mode = NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_FILL; break;
	default:
		assert(false);
		break;
	}
	return mode;
}

static unsigned int kelvin_map_texgen(uint32_t parameter, unsigned int channel)
{
	assert(channel < 4);
	unsigned int texgen;
	switch (parameter) {
	case NV097_SET_TEXGEN_S_DISABLE:
		texgen = NV_PGRAPH_CSV1_A_T0_S_DISABLE; break;
	case NV097_SET_TEXGEN_S_EYE_LINEAR:
		texgen = NV_PGRAPH_CSV1_A_T0_S_EYE_LINEAR; break;
	case NV097_SET_TEXGEN_S_OBJECT_LINEAR:
		texgen = NV_PGRAPH_CSV1_A_T0_S_OBJECT_LINEAR; break;
	case NV097_SET_TEXGEN_S_SPHERE_MAP:
		assert(channel < 2);
		texgen = NV_PGRAPH_CSV1_A_T0_S_SPHERE_MAP; break;
	case NV097_SET_TEXGEN_S_REFLECTION_MAP:
		assert(channel < 3);
		texgen = NV_PGRAPH_CSV1_A_T0_S_REFLECTION_MAP; break;
	case NV097_SET_TEXGEN_S_NORMAL_MAP:
		assert(channel < 3);
		texgen = NV_PGRAPH_CSV1_A_T0_S_NORMAL_MAP; break;
	default:
		assert(false);
		break;
	}
	return texgen;
}
