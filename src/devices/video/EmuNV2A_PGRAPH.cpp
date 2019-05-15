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

void (*pgraph_draw_arrays)(NV2AState *d);
void (*pgraph_draw_inline_buffer)(NV2AState *d);
void (*pgraph_draw_inline_array)(NV2AState *d);
void (*pgraph_draw_inline_elements)(NV2AState *d);
void (*pgraph_draw_state_update)(NV2AState *d);
void (*pgraph_draw_clear)(NV2AState *d);
void (*pgraph_update_surface)(NV2AState *d, bool upload, bool color_write, bool zeta_write);
void (*pgraph_report_get)(NV2AState *d);
void (*pgraph_report_clear)(NV2AState *d);

extern void OpenGL_init(PGRAPHState *pg);
extern void OpenGL_destroy(PGRAPHState *pg);

//static void pgraph_set_context_user(NV2AState *d, uint32_t value);
void pgraph_handle_method(NV2AState *d, unsigned int subchannel, unsigned int method, uint32_t parameter);
static void pgraph_log_method(unsigned int subchannel, unsigned int graphics_class, unsigned int method, uint32_t parameter);
static void pgraph_allocate_inline_buffer_vertices(PGRAPHState *pg, unsigned int attr);
static void pgraph_finish_inline_buffer_vertex(PGRAPHState *pg);
static bool pgraph_get_framebuffer_dirty(PGRAPHState *pg);
static bool pgraph_get_color_write_enabled(PGRAPHState *pg);
static bool pgraph_get_zeta_write_enabled(PGRAPHState *pg);
static void pgraph_set_surface_dirty(PGRAPHState *pg, bool color, bool zeta);
static void pgraph_apply_anti_aliasing_factor(PGRAPHState *pg, unsigned int *width, unsigned int *height);
static void pgraph_get_surface_dimensions(PGRAPHState *pg, unsigned int *width, unsigned int *height);

static float convert_f16_to_float(uint16_t f16);
static float convert_f24_to_float(uint32_t f24);
static uint8_t* convert_texture_data(const unsigned int color_format, const uint8_t *data, const uint8_t *palette_data, const unsigned int width, const unsigned int height, const unsigned int depth, const unsigned int row_pitch, const unsigned int slice_pitch);
static guint texture_key_hash(gconstpointer key);
static gboolean texture_key_equal(gconstpointer a, gconstpointer b);
static void texture_key_destroy(gpointer data);
static guint shader_hash(gconstpointer key);
static gboolean shader_equal(gconstpointer a, gconstpointer b);
static unsigned int kelvin_map_stencil_op(uint32_t parameter);
static unsigned int kelvin_map_polygon_mode(uint32_t parameter);
static unsigned int kelvin_map_texgen(uint32_t parameter, unsigned int channel);
static uint64_t fnv_hash(const uint8_t *data, size_t len);
static uint64_t fast_hash(const uint8_t *data, size_t len, unsigned int samples);

/* PGRAPH - accelerated 2d/3d drawing engine */
DEVICE_READ32(PGRAPH)
{
	qemu_mutex_lock(&d->pgraph.pgraph_lock);

	DEVICE_READ32_SWITCH() {
	case NV_PGRAPH_INTR:
		result = d->pgraph.pending_interrupts;
		break;
	case NV_PGRAPH_INTR_EN:
		result = d->pgraph.enabled_interrupts;
		break;
	default:
		DEVICE_READ32_REG(pgraph); // Was : DEBUG_READ32_UNHANDLED(PGRAPH);
	}

	qemu_mutex_unlock(&d->pgraph.pgraph_lock);

//    reg_log_read(NV_PGRAPH, addr, r);

	DEVICE_READ32_END(PGRAPH);
}

DEVICE_WRITE32(PGRAPH)
{
//    reg_log_write(NV_PGRAPH, addr, val);

	qemu_mutex_lock(&d->pgraph.pgraph_lock);

	switch (addr) {
	case NV_PGRAPH_INTR:
		d->pgraph.pending_interrupts &= ~value;
		qemu_cond_broadcast(&d->pgraph.interrupt_cond);
		break;
	case NV_PGRAPH_INTR_EN:
		d->pgraph.enabled_interrupts = value;
		break;
	case NV_PGRAPH_INCREMENT:
		if (value & NV_PGRAPH_INCREMENT_READ_3D) {
			SET_MASK(d->pgraph.regs[NV_PGRAPH_SURFACE],
				NV_PGRAPH_SURFACE_READ_3D,
				(GET_MASK(d->pgraph.regs[NV_PGRAPH_SURFACE],
					NV_PGRAPH_SURFACE_READ_3D) + 1)
				% GET_MASK(d->pgraph.regs[NV_PGRAPH_SURFACE],
					NV_PGRAPH_SURFACE_MODULO_3D));
			qemu_cond_broadcast(&d->pgraph.flip_3d);
		}
		break;
	case NV_PGRAPH_CHANNEL_CTX_TRIGGER: {
		xbaddr context_address =
			GET_MASK(d->pgraph.regs[NV_PGRAPH_CHANNEL_CTX_POINTER], NV_PGRAPH_CHANNEL_CTX_POINTER_INST) << 4;

		if (value & NV_PGRAPH_CHANNEL_CTX_TRIGGER_READ_IN) {
			unsigned pgraph_channel_id =
				GET_MASK(d->pgraph.regs[NV_PGRAPH_CTX_USER], NV_PGRAPH_CTX_USER_CHID);

			NV2A_DPRINTF("PGRAPH: read channel %d context from %" HWADDR_PRIx "\n",
				pgraph_channel_id, context_address);

			uint8_t *context_ptr = d->pramin.ramin_ptr + context_address;
			uint32_t context_user = ldl_le_p((uint32_t*)context_ptr);

			NV2A_DPRINTF("    - CTX_USER = 0x%08X\n", context_user);

			d->pgraph.regs[NV_PGRAPH_CTX_USER] = context_user;
			// pgraph_set_context_user(d, context_user);
		}
		if (value & NV_PGRAPH_CHANNEL_CTX_TRIGGER_WRITE_OUT) {
			/* do stuff ... */
		}

		break;
    }
	default: 
		DEVICE_WRITE32_REG(pgraph); // Was : DEBUG_WRITE32_UNHANDLED(PGRAPH);
		break;
	}

    // events
    switch (addr) {
    case NV_PGRAPH_FIFO:
        qemu_cond_broadcast(&d->pgraph.fifo_access_cond);
        break;
    }

	qemu_mutex_unlock(&d->pgraph.pgraph_lock);

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
        d->pgraph.regs[NV_PGRAPH_CTX_CONTROL] & NV_PGRAPH_CTX_CONTROL_CHID;
    assert(channel_valid);

    unsigned channel_id = GET_MASK(pg->regs[NV_PGRAPH_CTX_USER], NV_PGRAPH_CTX_USER_CHID);

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

        pg->regs[NV_PGRAPH_CTX_CACHE1 + subchannel * 4] = ctx_1;
        pg->regs[NV_PGRAPH_CTX_CACHE2 + subchannel * 4] = ctx_2;
        pg->regs[NV_PGRAPH_CTX_CACHE3 + subchannel * 4] = ctx_3;
        pg->regs[NV_PGRAPH_CTX_CACHE4 + subchannel * 4] = ctx_4;
        pg->regs[NV_PGRAPH_CTX_CACHE5 + subchannel * 4] = ctx_5;
    }

    // is this right?
    pg->regs[NV_PGRAPH_CTX_SWITCH1] = pg->regs[NV_PGRAPH_CTX_CACHE1 + subchannel * 4];
    pg->regs[NV_PGRAPH_CTX_SWITCH2] = pg->regs[NV_PGRAPH_CTX_CACHE2 + subchannel * 4];
    pg->regs[NV_PGRAPH_CTX_SWITCH3] = pg->regs[NV_PGRAPH_CTX_CACHE3 + subchannel * 4];
    pg->regs[NV_PGRAPH_CTX_SWITCH4] = pg->regs[NV_PGRAPH_CTX_CACHE4 + subchannel * 4];
    pg->regs[NV_PGRAPH_CTX_SWITCH5] = pg->regs[NV_PGRAPH_CTX_CACHE5 + subchannel * 4];

    uint32_t graphics_class = GET_MASK(pg->regs[NV_PGRAPH_CTX_SWITCH1],
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
			pg->regs[NV_PGRAPH_PATT_COLOR0] = parameter;
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

				NV2A_DPRINTF("NV09F_SET_OPERATION_SRCCOPY");

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

				xbaddr source_dma_len, dest_dma_len;
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

				SET_MASK(pg->regs[NV_PGRAPH_TRAPPED_ADDR],
					NV_PGRAPH_TRAPPED_ADDR_CHID, channel_id);
				SET_MASK(pg->regs[NV_PGRAPH_TRAPPED_ADDR],
					NV_PGRAPH_TRAPPED_ADDR_SUBCH, subchannel);
				SET_MASK(pg->regs[NV_PGRAPH_TRAPPED_ADDR],
					NV_PGRAPH_TRAPPED_ADDR_MTHD, method);
				pg->regs[NV_PGRAPH_TRAPPED_DATA_LOW] = parameter;
				pg->regs[NV_PGRAPH_NSOURCE] = NV_PGRAPH_NSOURCE_NOTIFICATION; /* TODO: check this */
				pg->pending_interrupts |= NV_PGRAPH_INTR_ERROR;

				qemu_mutex_unlock(&pg->pgraph_lock);
				qemu_mutex_lock_iothread();
				update_irq(d);
				qemu_mutex_lock(&pg->pgraph_lock);
				qemu_mutex_unlock_iothread();

				while (pg->pending_interrupts & NV_PGRAPH_INTR_ERROR) {
					qemu_cond_wait(&pg->interrupt_cond, &pg->pgraph_lock);
				}
			}
			break;

		case NV097_WAIT_FOR_IDLE:
			pgraph_update_surface(d, false, true, true);
			break;


		case NV097_SET_FLIP_READ:
			SET_MASK(pg->regs[NV_PGRAPH_SURFACE], NV_PGRAPH_SURFACE_READ_3D,
				parameter);
			break;
		case NV097_SET_FLIP_WRITE:
			SET_MASK(pg->regs[NV_PGRAPH_SURFACE], NV_PGRAPH_SURFACE_WRITE_3D,
				parameter);
			break;
		case NV097_SET_FLIP_MODULO:
			SET_MASK(pg->regs[NV_PGRAPH_SURFACE], NV_PGRAPH_SURFACE_MODULO_3D,
				parameter);
			break;
		case NV097_FLIP_INCREMENT_WRITE: {
			NV2A_DPRINTF("flip increment write %d -> ",
				GET_MASK(pg->regs[NV_PGRAPH_SURFACE],
					NV_PGRAPH_SURFACE_WRITE_3D));
			SET_MASK(pg->regs[NV_PGRAPH_SURFACE],
				NV_PGRAPH_SURFACE_WRITE_3D,
				(GET_MASK(pg->regs[NV_PGRAPH_SURFACE],
					NV_PGRAPH_SURFACE_WRITE_3D) + 1)
				% GET_MASK(pg->regs[NV_PGRAPH_SURFACE],
					NV_PGRAPH_SURFACE_MODULO_3D));
			NV2A_DPRINTF("%d\n",
				GET_MASK(pg->regs[NV_PGRAPH_SURFACE],
					NV_PGRAPH_SURFACE_WRITE_3D));

#ifdef __APPLE__
			if (glFrameTerminatorGREMEDY) {
				glFrameTerminatorGREMEDY();
			}
#endif // __APPLE__
			break;
		}
		case NV097_FLIP_STALL:
			pgraph_update_surface(d, false, true, true);


			// TODO: Fix this (why does it hang?)
			/* while (true) */ {
				uint32_t surface = pg->regs[NV_PGRAPH_SURFACE];
				NV2A_DPRINTF("flip stall read: %d, write: %d, modulo: %d\n",
					GET_MASK(surface, NV_PGRAPH_SURFACE_READ_3D),
					GET_MASK(surface, NV_PGRAPH_SURFACE_WRITE_3D),
					GET_MASK(surface, NV_PGRAPH_SURFACE_MODULO_3D));

				if (GET_MASK(surface, NV_PGRAPH_SURFACE_READ_3D)
					!= GET_MASK(surface, NV_PGRAPH_SURFACE_WRITE_3D)) {
					break;
				}

				//qemu_cond_wait(&pg->flip_3d, &pg->lock);
			}

			// TODO: Remove this when the AMD crash is solved in vblank_thread
			NV2ADevice::UpdateHostDisplay(d);
			NV2A_DPRINTF("flip stall done\n");
			break;

		case NV097_SET_CONTEXT_DMA_NOTIFIES:
			pg->dma_notifies = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_A:
			pg->dma_a = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_B:
			pg->dma_b = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_STATE:
			pg->dma_state = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_COLOR:
			/* try to get any straggling draws in before the surface's changed :/ */
			pgraph_update_surface(d, false, true, true);

			pg->dma_color = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_ZETA:
			pg->dma_zeta = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_VERTEX_A:
			pg->dma_vertex_a = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_VERTEX_B:
			pg->dma_vertex_b = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_SEMAPHORE:
			pg->dma_semaphore = parameter;
			break;
		case NV097_SET_CONTEXT_DMA_REPORT:
			pg->dma_report = parameter;
			break;

		case NV097_SET_SURFACE_CLIP_HORIZONTAL:
			pgraph_update_surface(d, false, true, true);

			pg->surface_shape.clip_x =
				GET_MASK(parameter, NV097_SET_SURFACE_CLIP_HORIZONTAL_X);
			pg->surface_shape.clip_width =
				GET_MASK(parameter, NV097_SET_SURFACE_CLIP_HORIZONTAL_WIDTH);
			break;
		case NV097_SET_SURFACE_CLIP_VERTICAL:
			pgraph_update_surface(d, false, true, true);
			
			pg->surface_shape.clip_y =
				GET_MASK(parameter, NV097_SET_SURFACE_CLIP_VERTICAL_Y);
			pg->surface_shape.clip_height =
				GET_MASK(parameter, NV097_SET_SURFACE_CLIP_VERTICAL_HEIGHT);
			break;
		case NV097_SET_SURFACE_FORMAT:
			pgraph_update_surface(d, false, true, true);

			pg->surface_shape.color_format =
				GET_MASK(parameter, NV097_SET_SURFACE_FORMAT_COLOR);
			pg->surface_shape.zeta_format =
				GET_MASK(parameter, NV097_SET_SURFACE_FORMAT_ZETA);
			pg->surface_type =
				GET_MASK(parameter, NV097_SET_SURFACE_FORMAT_TYPE);
			pg->surface_shape.anti_aliasing =
				GET_MASK(parameter, NV097_SET_SURFACE_FORMAT_ANTI_ALIASING);
			pg->surface_shape.log_width =
				GET_MASK(parameter, NV097_SET_SURFACE_FORMAT_WIDTH);
			pg->surface_shape.log_height =
				GET_MASK(parameter, NV097_SET_SURFACE_FORMAT_HEIGHT);
			break;
		case NV097_SET_SURFACE_PITCH:
			pgraph_update_surface(d, false, true, true);

			pg->surface_color.pitch =
				GET_MASK(parameter, NV097_SET_SURFACE_PITCH_COLOR);
			pg->surface_zeta.pitch =
				GET_MASK(parameter, NV097_SET_SURFACE_PITCH_ZETA);
			break;
		case NV097_SET_SURFACE_COLOR_OFFSET:
			pgraph_update_surface(d, false, true, true);

			pg->surface_color.offset = parameter;
			break;
		case NV097_SET_SURFACE_ZETA_OFFSET:
			pgraph_update_surface(d, false, true, true);

			pg->surface_zeta.offset = parameter;
			break;

		CASE_8(NV097_SET_COMBINER_ALPHA_ICW, 4) :
			slot = (method - NV097_SET_COMBINER_ALPHA_ICW) / 4;
			pg->regs[NV_PGRAPH_COMBINEALPHAI0 + slot * 4] = parameter;
			break;

		case NV097_SET_COMBINER_SPECULAR_FOG_CW0:
			pg->regs[NV_PGRAPH_COMBINESPECFOG0] = parameter;
			break;

		case NV097_SET_COMBINER_SPECULAR_FOG_CW1:
			pg->regs[NV_PGRAPH_COMBINESPECFOG1] = parameter;
			break;

		CASE_4(NV097_SET_TEXTURE_ADDRESS, 64):
			slot = (method - NV097_SET_TEXTURE_ADDRESS) / 64;
			pg->regs[NV_PGRAPH_TEXADDRESS0 + slot * 4] = parameter;
			break;
		case NV097_SET_CONTROL0: {
			pgraph_update_surface(d, false, true, true);

			bool stencil_write_enable =
				parameter & NV097_SET_CONTROL0_STENCIL_WRITE_ENABLE;
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_STENCIL_WRITE_ENABLE,
				stencil_write_enable);

			uint32_t z_format = GET_MASK(parameter, NV097_SET_CONTROL0_Z_FORMAT);
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_Z_FORMAT, z_format);

			bool z_perspective =
				parameter & NV097_SET_CONTROL0_Z_PERSPECTIVE_ENABLE;
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_Z_PERSPECTIVE_ENABLE,
				z_perspective);

			int color_space_convert =
				GET_MASK(parameter, NV097_SET_CONTROL0_COLOR_SPACE_CONVERT);
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_CSCONVERT,
				color_space_convert);
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
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_3], NV_PGRAPH_CONTROL_3_FOG_MODE,
				mode);
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
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D], NV_PGRAPH_CSV0_D_FOGGENMODE, mode);
			break;
		}
		case NV097_SET_FOG_ENABLE:
			/*
			FIXME: There is also:
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D], NV_PGRAPH_CSV0_D_FOGENABLE,
			parameter);
			*/
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_3], NV_PGRAPH_CONTROL_3_FOGENABLE,
				parameter);
			break;
		case NV097_SET_FOG_COLOR: {
			/* parameter channels are ABGR, PGRAPH channels are ARGB */
			uint8_t alpha = GET_MASK(parameter, NV097_SET_FOG_COLOR_ALPHA);
			uint8_t blue = GET_MASK(parameter, NV097_SET_FOG_COLOR_BLUE);
			uint8_t green = GET_MASK(parameter, NV097_SET_FOG_COLOR_GREEN);
			uint8_t red = GET_MASK(parameter, NV097_SET_FOG_COLOR_RED);
			SET_MASK(pg->regs[NV_PGRAPH_FOGCOLOR], NV_PGRAPH_FOGCOLOR_ALPHA, alpha);
			SET_MASK(pg->regs[NV_PGRAPH_FOGCOLOR], NV_PGRAPH_FOGCOLOR_RED, red);
			SET_MASK(pg->regs[NV_PGRAPH_FOGCOLOR], NV_PGRAPH_FOGCOLOR_GREEN, green);
			SET_MASK(pg->regs[NV_PGRAPH_FOGCOLOR], NV_PGRAPH_FOGCOLOR_BLUE, blue);
			break;
		}
		case NV097_SET_WINDOW_CLIP_TYPE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_WINDOWCLIPTYPE, parameter);
			break;
		CASE_8(NV097_SET_WINDOW_CLIP_HORIZONTAL, 4):
			slot = (method - NV097_SET_WINDOW_CLIP_HORIZONTAL) / 4;
			pg->regs[NV_PGRAPH_WINDOWCLIPX0 + slot * 4] = parameter;
			break;
		CASE_8(NV097_SET_WINDOW_CLIP_VERTICAL, 4):
			slot = (method - NV097_SET_WINDOW_CLIP_VERTICAL) / 4;
			pg->regs[NV_PGRAPH_WINDOWCLIPY0 + slot * 4] = parameter;
			break;
		case NV097_SET_ALPHA_TEST_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_ALPHATESTENABLE, parameter);
			break;
		case NV097_SET_BLEND_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_BLEND], NV_PGRAPH_BLEND_EN, parameter);
			break;
		case NV097_SET_CULL_FACE_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_CULLENABLE,
				parameter);
			break;
		case NV097_SET_DEPTH_TEST_ENABLE:
			// Test-case : Whiplash
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0], NV_PGRAPH_CONTROL_0_ZENABLE,
				parameter);
			break;
		case NV097_SET_DITHER_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_DITHERENABLE, parameter);
			break;
		case NV097_SET_LIGHTING_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_C], NV_PGRAPH_CSV0_C_LIGHTING,
				parameter);
			break;
		case NV097_SET_SKIN_MODE:
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D], NV_PGRAPH_CSV0_D_SKIN,
				parameter);
			break;
		case NV097_SET_STENCIL_TEST_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_1],
				NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE, parameter);
			break;
		case NV097_SET_POLY_OFFSET_POINT_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_POFFSETPOINTENABLE, parameter);
			break;
		case NV097_SET_POLY_OFFSET_LINE_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_POFFSETLINEENABLE, parameter);
			break;
		case NV097_SET_POLY_OFFSET_FILL_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_POFFSETFILLENABLE, parameter);
			break;
		case NV097_SET_ALPHA_FUNC:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_ALPHAFUNC, parameter & 0xF);
			break;
		case NV097_SET_ALPHA_REF:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_ALPHAREF, parameter);
			break;
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
			SET_MASK(pg->regs[NV_PGRAPH_BLEND], NV_PGRAPH_BLEND_SFACTOR, factor);

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
			SET_MASK(pg->regs[NV_PGRAPH_BLEND], NV_PGRAPH_BLEND_DFACTOR, factor);

			break;
		}

		case NV097_SET_BLEND_COLOR:
			pg->regs[NV_PGRAPH_BLENDCOLOR] = parameter;
			break;

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
			SET_MASK(pg->regs[NV_PGRAPH_BLEND], NV_PGRAPH_BLEND_EQN, equation);

			break;
		}

		case NV097_SET_DEPTH_FUNC:
			// Test-case : Whiplash
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0], NV_PGRAPH_CONTROL_0_ZFUNC,
				parameter & 0xF);
			break;

		case NV097_SET_COLOR_MASK: {
			pg->surface_color.write_enabled_cache |= pgraph_get_color_write_enabled(pg);

			bool alpha = parameter & NV097_SET_COLOR_MASK_ALPHA_WRITE_ENABLE;
			bool red = parameter & NV097_SET_COLOR_MASK_RED_WRITE_ENABLE;
			bool green = parameter & NV097_SET_COLOR_MASK_GREEN_WRITE_ENABLE;
			bool blue = parameter & NV097_SET_COLOR_MASK_BLUE_WRITE_ENABLE;
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE, alpha);
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE, red);
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE, green);
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE, blue);
			break;
		}
		case NV097_SET_DEPTH_MASK:
			pg->surface_zeta.write_enabled_cache |= pgraph_get_zeta_write_enabled(pg);

			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_0],
				NV_PGRAPH_CONTROL_0_ZWRITEENABLE, parameter);
			break;
		case NV097_SET_STENCIL_MASK:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_1],
				NV_PGRAPH_CONTROL_1_STENCIL_MASK_WRITE, parameter);
			break;
		case NV097_SET_STENCIL_FUNC:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_1],
				NV_PGRAPH_CONTROL_1_STENCIL_FUNC, parameter & 0xF);
			break;
		case NV097_SET_STENCIL_FUNC_REF:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_1],
				NV_PGRAPH_CONTROL_1_STENCIL_REF, parameter);
			break;
		case NV097_SET_STENCIL_FUNC_MASK:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_1],
				NV_PGRAPH_CONTROL_1_STENCIL_MASK_READ, parameter);
			break;
		case NV097_SET_STENCIL_OP_FAIL:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_2],
				NV_PGRAPH_CONTROL_2_STENCIL_OP_FAIL,
				kelvin_map_stencil_op(parameter));
			break;
		case NV097_SET_STENCIL_OP_ZFAIL:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_2],
				NV_PGRAPH_CONTROL_2_STENCIL_OP_ZFAIL,
				kelvin_map_stencil_op(parameter));
			break;
		case NV097_SET_STENCIL_OP_ZPASS:
			SET_MASK(pg->regs[NV_PGRAPH_CONTROL_2],
				NV_PGRAPH_CONTROL_2_STENCIL_OP_ZPASS,
				kelvin_map_stencil_op(parameter));
			break;

		case NV097_SET_POLYGON_OFFSET_SCALE_FACTOR:
			pg->regs[NV_PGRAPH_ZOFFSETFACTOR] = parameter;
			break;
		case NV097_SET_POLYGON_OFFSET_BIAS:
			pg->regs[NV_PGRAPH_ZOFFSETBIAS] = parameter;
			break;
		case NV097_SET_FRONT_POLYGON_MODE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_FRONTFACEMODE,
				kelvin_map_polygon_mode(parameter));
			break;
		case NV097_SET_BACK_POLYGON_MODE:
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_BACKFACEMODE,
				kelvin_map_polygon_mode(parameter));
			break;
		case NV097_SET_CLIP_MIN:
			pg->regs[NV_PGRAPH_ZCLIPMIN] = parameter;
			break;
		case NV097_SET_CLIP_MAX:
			pg->regs[NV_PGRAPH_ZCLIPMAX] = parameter;
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
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_CULLCTRL,
				face);
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
			SET_MASK(pg->regs[NV_PGRAPH_SETUPRASTER],
				NV_PGRAPH_SETUPRASTER_FRONTFACE,
				ccw ? 1 : 0);
			break;
		}
		case NV097_SET_NORMALIZATION_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_C],
				NV_PGRAPH_CSV0_C_NORMALIZATION_ENABLE,
				parameter);
			break;

		case NV097_SET_LIGHT_ENABLE_MASK:
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D],
				NV_PGRAPH_CSV0_D_LIGHTS,
				parameter);
			break;

		CASE_4(NV097_SET_TEXGEN_S, 16) : {
			slot = (method - NV097_SET_TEXGEN_S) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_S
				: NV_PGRAPH_CSV1_A_T0_S;
			SET_MASK(pg->regs[reg], mask, kelvin_map_texgen(parameter, 0));
			break;
		}
		CASE_4(NV097_SET_TEXGEN_T, 16) : {
			slot = (method - NV097_SET_TEXGEN_T) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_T
				: NV_PGRAPH_CSV1_A_T0_T;
			SET_MASK(pg->regs[reg], mask, kelvin_map_texgen(parameter, 1));
			break;
		}
		CASE_4(NV097_SET_TEXGEN_R, 16) : {
			slot = (method - NV097_SET_TEXGEN_R) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_R
				: NV_PGRAPH_CSV1_A_T0_R;
			SET_MASK(pg->regs[reg], mask, kelvin_map_texgen(parameter, 2));
			break;
		}
		CASE_4(NV097_SET_TEXGEN_Q, 16) : {
			slot = (method - NV097_SET_TEXGEN_Q) / 16;
			unsigned int reg = (slot < 2) ? NV_PGRAPH_CSV1_A
				: NV_PGRAPH_CSV1_B;
			unsigned int mask = (slot % 2) ? NV_PGRAPH_CSV1_A_T1_Q
				: NV_PGRAPH_CSV1_A_T0_Q;
			SET_MASK(pg->regs[reg], mask, kelvin_map_texgen(parameter, 3));
			break;
		}
		CASE_4(NV097_SET_TEXTURE_MATRIX_ENABLE, 4) :
			slot = (method - NV097_SET_TEXTURE_MATRIX_ENABLE) / 4;
			pg->texture_matrix_enable[slot] = parameter;
			break;

		CASE_16(NV097_SET_PROJECTION_MATRIX, 4) : {
			slot = (method - NV097_SET_PROJECTION_MATRIX) / 4;
			// pg->projection_matrix[slot] = *(float*)&parameter;
			unsigned int row = NV_IGRAPH_XF_XFCTX_PMAT0 + slot / 4;
			pg->vsh_constants[row][slot % 4] = parameter;
			pg->vsh_constants_dirty[row] = true;
			break;
		}

		CASE_64(NV097_SET_MODEL_VIEW_MATRIX, 4) : {
			slot = (method - NV097_SET_MODEL_VIEW_MATRIX) / 4;
			unsigned int matnum = slot / 16;
			unsigned int entry = slot % 16;
			unsigned int row = NV_IGRAPH_XF_XFCTX_MMAT0 + matnum * 8 + entry / 4;
			pg->vsh_constants[row][entry % 4] = parameter;
			pg->vsh_constants_dirty[row] = true;
			break;
		}

		CASE_64(NV097_SET_INVERSE_MODEL_VIEW_MATRIX, 4) : {
			slot = (method - NV097_SET_INVERSE_MODEL_VIEW_MATRIX) / 4;
			unsigned int matnum = slot / 16;
			unsigned int entry = slot % 16;
			unsigned int row = NV_IGRAPH_XF_XFCTX_IMMAT0 + matnum * 8 + entry / 4;
			pg->vsh_constants[row][entry % 4] = parameter;
			pg->vsh_constants_dirty[row] = true;
			break;
		}

		CASE_16(NV097_SET_COMPOSITE_MATRIX, 4) : {
			slot = (method - NV097_SET_COMPOSITE_MATRIX) / 4;
			unsigned int row = NV_IGRAPH_XF_XFCTX_CMAT0 + slot / 4;
			pg->vsh_constants[row][slot % 4] = parameter;
			pg->vsh_constants_dirty[row] = true;
			break;
		}

		CASE_64(NV097_SET_TEXTURE_MATRIX, 4) : {
			slot = (method - NV097_SET_TEXTURE_MATRIX) / 4;
			unsigned int tex = slot / 16;
			unsigned int entry = slot % 16;
			unsigned int row = NV_IGRAPH_XF_XFCTX_T0MAT + tex * 8 + entry / 4;
			pg->vsh_constants[row][entry % 4] = parameter;
			pg->vsh_constants_dirty[row] = true;
			break;
		}

		CASE_3(NV097_SET_FOG_PARAMS, 4) :
			slot = (method - NV097_SET_FOG_PARAMS) / 4;
			pg->regs[NV_PGRAPH_FOGPARAM0 + slot * 4] = parameter;
			/* Cxbx note: slot = 2 is right after slot = 1 */
			pg->ltctxa[NV_IGRAPH_XF_LTCTXA_FOG_K][slot] = parameter;
			pg->ltctxa_dirty[NV_IGRAPH_XF_LTCTXA_FOG_K] = true;
			break;

		/* Handles NV097_SET_TEXGEN_PLANE_S,T,R,Q */
		CASE_64(NV097_SET_TEXGEN_PLANE_S, 4) : {
			slot = (method - NV097_SET_TEXGEN_PLANE_S) / 4;
			unsigned int tex = slot / 16;
			unsigned int entry = slot % 16;
			unsigned int row = NV_IGRAPH_XF_XFCTX_TG0MAT + tex * 8 + entry / 4;
			pg->vsh_constants[row][entry % 4] = parameter;
			pg->vsh_constants_dirty[row] = true;
			break;
		}

		case NV097_SET_TEXGEN_VIEW_MODEL:
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D], NV_PGRAPH_CSV0_D_TEXGEN_REF,
				parameter);
			break;

		CASE_4(NV097_SET_FOG_PLANE, 4):
			slot = (method - NV097_SET_FOG_PLANE) / 4;
			pg->vsh_constants[NV_IGRAPH_XF_XFCTX_FOG][slot] = parameter;
			pg->vsh_constants_dirty[NV_IGRAPH_XF_XFCTX_FOG] = true;
			break;

		CASE_3(NV097_SET_SCENE_AMBIENT_COLOR, 4):
			slot = (method - NV097_SET_SCENE_AMBIENT_COLOR) / 4;
			// ??
			pg->ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][slot] = parameter;
			pg->ltctxa_dirty[NV_IGRAPH_XF_LTCTXA_FR_AMB] = true;
			break;

		CASE_4(NV097_SET_VIEWPORT_OFFSET, 4):
			slot = (method - NV097_SET_VIEWPORT_OFFSET) / 4;
			pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPOFF][slot] = parameter;
			pg->vsh_constants_dirty[NV_IGRAPH_XF_XFCTX_VPOFF] = true;
			break;

		CASE_4(NV097_SET_EYE_POSITION, 4):
			slot = (method - NV097_SET_EYE_POSITION) / 4;
			pg->vsh_constants[NV_IGRAPH_XF_XFCTX_EYEP][slot] = parameter;
			pg->vsh_constants_dirty[NV_IGRAPH_XF_XFCTX_EYEP] = true;
			break;

		CASE_8(NV097_SET_COMBINER_FACTOR0, 4):
			slot = (method - NV097_SET_COMBINER_FACTOR0) / 4;
			pg->regs[NV_PGRAPH_COMBINEFACTOR0 + slot * 4] = parameter;
			break;

		CASE_8(NV097_SET_COMBINER_FACTOR1, 4):
			slot = (method - NV097_SET_COMBINER_FACTOR1) / 4;
			pg->regs[NV_PGRAPH_COMBINEFACTOR1 + slot * 4] = parameter;
			break;

		CASE_8(NV097_SET_COMBINER_ALPHA_OCW, 4):
			slot = (method - NV097_SET_COMBINER_ALPHA_OCW) / 4;
			pg->regs[NV_PGRAPH_COMBINEALPHAO0 + slot * 4] = parameter;
			break;

		CASE_8(NV097_SET_COMBINER_COLOR_ICW, 4):
			slot = (method - NV097_SET_COMBINER_COLOR_ICW) / 4;
			pg->regs[NV_PGRAPH_COMBINECOLORI0 + slot * 4] = parameter;
			break;

		CASE_4(NV097_SET_VIEWPORT_SCALE, 4):
			slot = (method - NV097_SET_VIEWPORT_SCALE) / 4;
			pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL][slot] = parameter;
			pg->vsh_constants_dirty[NV_IGRAPH_XF_XFCTX_VPSCL] = true;
			break;

		CASE_32(NV097_SET_TRANSFORM_PROGRAM, 4) : {

			slot = (method - NV097_SET_TRANSFORM_PROGRAM) / 4;

			int program_load = GET_MASK(pg->regs[NV_PGRAPH_CHEOPS_OFFSET],
				NV_PGRAPH_CHEOPS_OFFSET_PROG_LD_PTR);

			assert(program_load < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
			pg->program_data[program_load][slot % 4] = parameter;

			if (slot % 4 == 3) {
				SET_MASK(pg->regs[NV_PGRAPH_CHEOPS_OFFSET],
					NV_PGRAPH_CHEOPS_OFFSET_PROG_LD_PTR, program_load + 1);
			}

			break;
		}

		CASE_32(NV097_SET_TRANSFORM_CONSTANT, 4): {

			slot = (method - NV097_SET_TRANSFORM_CONSTANT) / 4;

			int const_load = GET_MASK(pg->regs[NV_PGRAPH_CHEOPS_OFFSET],
									  NV_PGRAPH_CHEOPS_OFFSET_CONST_LD_PTR);

			assert(const_load < NV2A_VERTEXSHADER_CONSTANTS);
			// VertexShaderConstant *vsh_constant = &pg->vsh_constants[const_load];
			pg->vsh_constants_dirty[const_load] |=
				(parameter != pg->vsh_constants[const_load][slot%4]);
			pg->vsh_constants[const_load][slot%4] = parameter;

			if (slot % 4 == 3) {
				SET_MASK(pg->regs[NV_PGRAPH_CHEOPS_OFFSET],
						 NV_PGRAPH_CHEOPS_OFFSET_CONST_LD_PTR, const_load+1);
			}
			break;
		}

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

		/* Handles NV097_SET_BACK_LIGHT_* */
		CASE_128(NV097_SET_BACK_LIGHT_AMBIENT_COLOR, 4): {
			slot = (method - NV097_SET_BACK_LIGHT_AMBIENT_COLOR) / 4;
			unsigned int part = NV097_SET_BACK_LIGHT_AMBIENT_COLOR / 4 + slot % 16;
			slot /= 16; /* [Light index] */
			assert(slot < 8);
			switch(part * 4) {
			CASE_3(NV097_SET_BACK_LIGHT_AMBIENT_COLOR, 4):
				part -= NV097_SET_BACK_LIGHT_AMBIENT_COLOR / 4;
				pg->ltctxb[NV_IGRAPH_XF_LTCTXB_L0_BAMB + slot*6][part] = parameter;
				pg->ltctxb_dirty[NV_IGRAPH_XF_LTCTXB_L0_BAMB + slot*6] = true;
				break;
			CASE_3(NV097_SET_BACK_LIGHT_DIFFUSE_COLOR, 4):
				part -= NV097_SET_BACK_LIGHT_DIFFUSE_COLOR / 4;
				pg->ltctxb[NV_IGRAPH_XF_LTCTXB_L0_BDIF + slot*6][part] = parameter;
				pg->ltctxb_dirty[NV_IGRAPH_XF_LTCTXB_L0_BDIF + slot*6] = true;
				break;
			CASE_3(NV097_SET_BACK_LIGHT_SPECULAR_COLOR, 4):
				part -= NV097_SET_BACK_LIGHT_SPECULAR_COLOR / 4;
				pg->ltctxb[NV_IGRAPH_XF_LTCTXB_L0_BSPC + slot*6][part] = parameter;
				pg->ltctxb_dirty[NV_IGRAPH_XF_LTCTXB_L0_BSPC + slot*6] = true;
				break;
			default:
				assert(false);
				break;
			}
			break;
		}
		/* Handles all the light source props except for NV097_SET_BACK_LIGHT_* */
		CASE_256(NV097_SET_LIGHT_AMBIENT_COLOR, 4): {
			slot = (method - NV097_SET_LIGHT_AMBIENT_COLOR) / 4;
			unsigned int part = NV097_SET_LIGHT_AMBIENT_COLOR / 4 + slot % 32;
			slot /= 32; /* [Light index] */
			assert(slot < 8);
			switch(part * 4) {
			CASE_3(NV097_SET_LIGHT_AMBIENT_COLOR, 4):
				part -= NV097_SET_LIGHT_AMBIENT_COLOR / 4;
				pg->ltctxb[NV_IGRAPH_XF_LTCTXB_L0_AMB + slot*6][part] = parameter;
				pg->ltctxb_dirty[NV_IGRAPH_XF_LTCTXB_L0_AMB + slot*6] = true;
				break;
			CASE_3(NV097_SET_LIGHT_DIFFUSE_COLOR, 4):
				part -= NV097_SET_LIGHT_DIFFUSE_COLOR / 4;
				pg->ltctxb[NV_IGRAPH_XF_LTCTXB_L0_DIF + slot*6][part] = parameter;
				pg->ltctxb_dirty[NV_IGRAPH_XF_LTCTXB_L0_DIF + slot*6] = true;
				break;
			CASE_3(NV097_SET_LIGHT_SPECULAR_COLOR, 4):
				part -= NV097_SET_LIGHT_SPECULAR_COLOR / 4;
				pg->ltctxb[NV_IGRAPH_XF_LTCTXB_L0_SPC + slot*6][part] = parameter;
				pg->ltctxb_dirty[NV_IGRAPH_XF_LTCTXB_L0_SPC + slot*6] = true;
				break;
			case NV097_SET_LIGHT_LOCAL_RANGE:
				pg->ltc1[NV_IGRAPH_XF_LTC1_r0 + slot][0] = parameter;
				pg->ltc1_dirty[NV_IGRAPH_XF_LTC1_r0 + slot] = true;
				break;
			CASE_3(NV097_SET_LIGHT_INFINITE_HALF_VECTOR, 4):
				part -= NV097_SET_LIGHT_INFINITE_HALF_VECTOR / 4;
				pg->light_infinite_half_vector[slot][part] = *(float*)&parameter;
				break;
			CASE_3(NV097_SET_LIGHT_INFINITE_DIRECTION, 4):
				part -= NV097_SET_LIGHT_INFINITE_DIRECTION / 4;
				pg->light_infinite_direction[slot][part] = *(float*)&parameter;
				break;
			CASE_3(NV097_SET_LIGHT_SPOT_FALLOFF, 4):
				part -= NV097_SET_LIGHT_SPOT_FALLOFF / 4;
				pg->ltctxa[NV_IGRAPH_XF_LTCTXA_L0_K + slot*2][part] = parameter;
				pg->ltctxa_dirty[NV_IGRAPH_XF_LTCTXA_L0_K + slot*2] = true;
				break;
			CASE_4(NV097_SET_LIGHT_SPOT_DIRECTION, 4):
				part -= NV097_SET_LIGHT_SPOT_DIRECTION / 4;
				pg->ltctxa[NV_IGRAPH_XF_LTCTXA_L0_SPT + slot*2][part] = parameter;
				pg->ltctxa_dirty[NV_IGRAPH_XF_LTCTXA_L0_SPT + slot*2] = true;
				break;
			CASE_3(NV097_SET_LIGHT_LOCAL_POSITION, 4):
				part -= NV097_SET_LIGHT_LOCAL_POSITION / 4;
				pg->light_local_position[slot][part] = *(float*)&parameter;
				break;
			CASE_3(NV097_SET_LIGHT_LOCAL_ATTENUATION, 4):
				part -= NV097_SET_LIGHT_LOCAL_ATTENUATION / 4;
				pg->light_local_attenuation[slot][part] = *(float*)&parameter;
				break;
			default:
				assert(false);
				break;
			}
			break;
		}

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

		CASE_16(NV097_SET_VERTEX_DATA_ARRAY_FORMAT, 4): {

			slot = (method - NV097_SET_VERTEX_DATA_ARRAY_FORMAT) / 4;
			VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];

			vertex_attribute->format =
				GET_MASK(parameter, NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE);
			vertex_attribute->count =
				GET_MASK(parameter, NV097_SET_VERTEX_DATA_ARRAY_FORMAT_SIZE);
			vertex_attribute->stride =
				GET_MASK(parameter, NV097_SET_VERTEX_DATA_ARRAY_FORMAT_STRIDE);

			NV2A_DPRINTF("vertex data array format=%d, count=%d, stride=%d\n",
				vertex_attribute->format,
				vertex_attribute->count,
				vertex_attribute->stride);

			{ // TODO : Make this part OpenGL_specific / opengl_enabled :
				vertex_attribute->gl_count = vertex_attribute->count;

				switch (vertex_attribute->format) {
				case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_UB_D3D:
					vertex_attribute->gl_type = GL_UNSIGNED_BYTE;
					vertex_attribute->gl_normalize = GL_TRUE;
					vertex_attribute->size = 1;
					assert(vertex_attribute->count == 4);
					// https://www.opengl.org/registry/specs/ARB/vertex_array_bgra.txt
					vertex_attribute->gl_count = GL_BGRA;
					vertex_attribute->needs_conversion = false;
					break;
				case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_UB_OGL:
					vertex_attribute->gl_type = GL_UNSIGNED_BYTE;
					vertex_attribute->gl_normalize = GL_TRUE;
					vertex_attribute->size = 1;
					vertex_attribute->needs_conversion = false;
					break;
				case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_S1:
					vertex_attribute->gl_type = GL_SHORT;
					vertex_attribute->gl_normalize = GL_TRUE;
					vertex_attribute->size = 2;
					vertex_attribute->needs_conversion = false;
					break;
				case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_F:
					vertex_attribute->gl_type = GL_FLOAT;
					vertex_attribute->gl_normalize = GL_FALSE;
					vertex_attribute->size = 4;
					vertex_attribute->needs_conversion = false;
					break;
				case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_S32K:
					vertex_attribute->gl_type = GL_SHORT;
					vertex_attribute->gl_normalize = GL_FALSE;
					vertex_attribute->size = 2;
					vertex_attribute->needs_conversion = false;
					break;
				case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_CMP:
					/* 3 signed, normalized components packed in 32-bits. (11,11,10) */
					vertex_attribute->size = 4;
					vertex_attribute->gl_type = GL_FLOAT;
					vertex_attribute->gl_normalize = GL_FALSE;
					vertex_attribute->needs_conversion = true;
					vertex_attribute->converted_size = sizeof(float);
					vertex_attribute->converted_count = 3 * vertex_attribute->count;
					break;
				default:
					fprintf(stderr, "Unknown vertex type: 0x%x\n", vertex_attribute->format);
					assert(false);
					break;
				}

				if (vertex_attribute->needs_conversion) {
					vertex_attribute->converted_elements = 0;
				} else {
					if (vertex_attribute->converted_buffer) {
						g_free(vertex_attribute->converted_buffer);
						vertex_attribute->converted_buffer = NULL;
					}
				}
			}

			break;
		}

		CASE_16(NV097_SET_VERTEX_DATA_ARRAY_OFFSET, 4): {

			slot = (method - NV097_SET_VERTEX_DATA_ARRAY_OFFSET) / 4;

			pg->vertex_attributes[slot].dma_select =
				parameter & 0x80000000;
			pg->vertex_attributes[slot].offset =
				parameter & 0x7fffffff;

			pg->vertex_attributes[slot].converted_elements = 0;

			break;
		}

		case NV097_SET_LOGIC_OP_ENABLE:
			SET_MASK(pg->regs[NV_PGRAPH_BLEND],
					 NV_PGRAPH_BLEND_LOGICOP_ENABLE, parameter);
			break;

		case NV097_SET_LOGIC_OP:
			SET_MASK(pg->regs[NV_PGRAPH_BLEND],
					 NV_PGRAPH_BLEND_LOGICOP, parameter & 0xF);
			break;

		case NV097_CLEAR_REPORT_VALUE:

			pgraph_report_clear(d);
			pg->zpass_pixel_count_result = 0;

			break;

		case NV097_SET_ZPASS_PIXEL_COUNT_ENABLE:
			pg->zpass_pixel_count_enable = parameter;
			break;

		case NV097_GET_REPORT: {
			/* FIXME: This was first intended to be watchpoint-based. However,
			 *        qemu / kvm only supports virtual-address watchpoints.
			 *        This'll do for now, but accuracy and performance with other
			 *        approaches could be better
			 */
			uint8_t type = GET_MASK(parameter, NV097_GET_REPORT_TYPE);
			assert(type == NV097_GET_REPORT_TYPE_ZPASS_PIXEL_CNT);
			hwaddr offset = GET_MASK(parameter, NV097_GET_REPORT_OFFSET);

			uint64_t timestamp = 0x0011223344556677; /* FIXME: Update timestamp?! */
			uint32_t done = 0;

			pgraph_report_get(d);

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

		CASE_3(NV097_SET_EYE_DIRECTION, 4):
			slot = (method - NV097_SET_EYE_DIRECTION) / 4;
			pg->ltctxa[NV_IGRAPH_XF_LTCTXA_EYED][slot] = parameter;
			pg->ltctxa_dirty[NV_IGRAPH_XF_LTCTXA_EYED] = true;
			break;

		case NV097_SET_BEGIN_END: {
			uint32_t control_0 = pg->regs[NV_PGRAPH_CONTROL_0];
			uint32_t control_1 = pg->regs[NV_PGRAPH_CONTROL_1];

			bool depth_test = control_0
				& NV_PGRAPH_CONTROL_0_ZENABLE;
			bool stencil_test = control_1
				& NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE;

			if (parameter == NV097_SET_BEGIN_END_OP_END) {

				if (pg->draw_arrays_length) {

					// NV2A_DPRINTF("Draw Arrays");

					assert(pg->inline_buffer_length == 0);
					assert(pg->inline_array_length == 0);
					assert(pg->inline_elements_length == 0);

					pgraph_draw_arrays(d);
				} else if (pg->inline_buffer_length) {

					// NV2A_DPRINTF("Inline Buffer");

					assert(pg->draw_arrays_length == 0);
					assert(pg->inline_array_length == 0);
					assert(pg->inline_elements_length == 0);

					pgraph_draw_inline_buffer(d);
				} else if (pg->inline_array_length) {

					// NV2A_DPRINTF("Inline Array");

					assert(pg->draw_arrays_length == 0);
					assert(pg->inline_buffer_length == 0);
					assert(pg->inline_elements_length == 0);

					pgraph_draw_inline_array(d);
				} else if (pg->inline_elements_length) {

					// NV2A_DPRINTF("Inline Elements");

					assert(pg->draw_arrays_length == 0);
					assert(pg->inline_buffer_length == 0);
					assert(pg->inline_array_length == 0);

					pgraph_draw_inline_elements(d);
				} else {
					NV2A_DPRINTF("EMPTY NV097_SET_BEGIN_END");
					assert(false);
				}
			} else {

				assert(parameter <= NV097_SET_BEGIN_END_OP_POLYGON);

				pg->primitive_mode = parameter;

				pgraph_draw_state_update(d);

				pg->inline_elements_length = 0;
				pg->inline_array_length = 0;
				pg->inline_buffer_length = 0;
				pg->draw_arrays_length = 0;
				pg->draw_arrays_max_count = 0;
			}

			pgraph_set_surface_dirty(pg, true, depth_test || stencil_test);
			break;
		}
		CASE_4(NV097_SET_TEXTURE_OFFSET, 64):
			slot = (method - NV097_SET_TEXTURE_OFFSET) / 64;
			pg->regs[NV_PGRAPH_TEXOFFSET0 + slot * 4] = parameter;
			pg->texture_dirty[slot] = true;
			break;
		CASE_4(NV097_SET_TEXTURE_FORMAT, 64): {
			slot = (method - NV097_SET_TEXTURE_FORMAT) / 64;

			bool dma_select =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_CONTEXT_DMA) == 2;
			bool cubemap =
				parameter & NV097_SET_TEXTURE_FORMAT_CUBEMAP_ENABLE;
			bool border_source =
				parameter & NV097_SET_TEXTURE_FORMAT_BORDER_SOURCE;
			unsigned int dimensionality =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_DIMENSIONALITY);
			unsigned int color_format =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_COLOR);
			unsigned int levels =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_MIPMAP_LEVELS);
			unsigned int log_width =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_BASE_SIZE_U);
			unsigned int log_height =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_BASE_SIZE_V);
			unsigned int log_depth =
				GET_MASK(parameter, NV097_SET_TEXTURE_FORMAT_BASE_SIZE_P);

			uint32_t *reg = &pg->regs[NV_PGRAPH_TEXFMT0 + slot * 4];
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_CONTEXT_DMA, dma_select);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_CUBEMAPENABLE, cubemap);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_BORDER_SOURCE, border_source);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_DIMENSIONALITY, dimensionality);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_COLOR, color_format);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_MIPMAP_LEVELS, levels);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_BASE_SIZE_U, log_width);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_BASE_SIZE_V, log_height);
			SET_MASK(*reg, NV_PGRAPH_TEXFMT0_BASE_SIZE_P, log_depth);

			pg->texture_dirty[slot] = true;
			break;
		}
		CASE_4(NV097_SET_TEXTURE_CONTROL0, 64):
			slot = (method - NV097_SET_TEXTURE_CONTROL0) / 64;
			pg->regs[NV_PGRAPH_TEXCTL0_0 + slot*4] = parameter;
			break;
		CASE_4(NV097_SET_TEXTURE_CONTROL1, 64):
			slot = (method - NV097_SET_TEXTURE_CONTROL1) / 64;
			pg->regs[NV_PGRAPH_TEXCTL1_0 + slot*4] = parameter;
			break;
		CASE_4(NV097_SET_TEXTURE_FILTER, 64):
			slot = (method - NV097_SET_TEXTURE_FILTER) / 64;
			pg->regs[NV_PGRAPH_TEXFILTER0 + slot * 4] = parameter;
			break;
		CASE_4(NV097_SET_TEXTURE_IMAGE_RECT, 64):
			slot = (method - NV097_SET_TEXTURE_IMAGE_RECT) / 64;
			pg->regs[NV_PGRAPH_TEXIMAGERECT0 + slot * 4] = parameter;
			pg->texture_dirty[slot] = true;
			break;
		CASE_4(NV097_SET_TEXTURE_PALETTE, 64): {
			slot = (method - NV097_SET_TEXTURE_PALETTE) / 64;

			bool dma_select =
				GET_MASK(parameter, NV097_SET_TEXTURE_PALETTE_CONTEXT_DMA) == 1;
			unsigned int length =
				GET_MASK(parameter, NV097_SET_TEXTURE_PALETTE_LENGTH);
			unsigned int offset =
				GET_MASK(parameter, NV097_SET_TEXTURE_PALETTE_OFFSET);

			uint32_t *reg = &pg->regs[NV_PGRAPH_TEXPALETTE0 + slot * 4];
			SET_MASK(*reg, NV_PGRAPH_TEXPALETTE0_CONTEXT_DMA, dma_select);
			SET_MASK(*reg, NV_PGRAPH_TEXPALETTE0_LENGTH, length);
			SET_MASK(*reg, NV_PGRAPH_TEXPALETTE0_OFFSET, offset);

			pg->texture_dirty[slot] = true;
			break;
		}

		CASE_4(NV097_SET_TEXTURE_BORDER_COLOR, 64):
			slot = (method - NV097_SET_TEXTURE_BORDER_COLOR) / 64;
			pg->regs[NV_PGRAPH_BORDERCOLOR0 + slot * 4] = parameter;
			break;
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x0, 64):
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x4, 64):
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x8, 64):
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0xc, 64):
			slot = (method - NV097_SET_TEXTURE_SET_BUMP_ENV_MAT) / 4;
			assert((slot / 16) > 0);
			slot -= 16;
			pg->bump_env_matrix[slot / 16][slot % 4] = *(float*)&parameter;
			break;

		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_SCALE, 64):
			slot = (method - NV097_SET_TEXTURE_SET_BUMP_ENV_SCALE) / 64;
			assert(slot > 0);
			slot--;
			pg->regs[NV_PGRAPH_BUMPSCALE1 + slot * 4] = parameter;
			break;
		CASE_4(NV097_SET_TEXTURE_SET_BUMP_ENV_OFFSET, 64):
			slot = (method - NV097_SET_TEXTURE_SET_BUMP_ENV_OFFSET) / 64;
			assert(slot > 0);
			slot--;
			pg->regs[NV_PGRAPH_BUMPOFFSET1 + slot * 4] = parameter;
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

			{ // TODO : Make this OpenGL_specific / opengl_enabled :
				assert(pg->draw_arrays_length < ARRAY_SIZE(pg->gl_draw_arrays_start));

				/* Attempt to connect primitives */
				if (pg->draw_arrays_length > 0) {
					unsigned int last_start =
						pg->gl_draw_arrays_start[pg->draw_arrays_length - 1];
					GLsizei* last_count =
						&pg->gl_draw_arrays_count[pg->draw_arrays_length - 1];
					if (start == (last_start + *last_count)) {
						*last_count += count;
						break;
					}
				}

				pg->gl_draw_arrays_start[pg->draw_arrays_length] = start;
				pg->gl_draw_arrays_count[pg->draw_arrays_length] = count;
				pg->draw_arrays_length++;
			}
			break;
		}
		case NV097_INLINE_ARRAY:
			assert(pg->inline_array_length < NV2A_MAX_BATCH_LENGTH);
			pg->inline_array[
				pg->inline_array_length++] = parameter;
			break;
		CASE_3(NV097_SET_EYE_VECTOR, 4):
			slot = (method - NV097_SET_EYE_VECTOR) / 4;
			pg->regs[NV_PGRAPH_EYEVEC0 + slot * 4] = parameter;
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
			/* FIXME: Is mapping to [-1,+1] correct? */
			vertex_attribute->inline_value[0] = ((int16_t)(parameter & 0xFFFF) * 2.0f + 1)
											 / 65535.0f;
			vertex_attribute->inline_value[1] = ((int16_t)(parameter >> 16) * 2.0f + 1)
											 / 65535.0f;
			/* FIXME: Should these really be set to 0.0 and 1.0 ? Conditions? */
			vertex_attribute->inline_value[2] = 0.0f;
			vertex_attribute->inline_value[3] = 1.0f;
			if (slot == 0) {
				pgraph_finish_inline_buffer_vertex(pg);
				assert(false); /* FIXME: Untested */
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
				assert(false); /* FIXME: Untested */
			}
			break;
		}
		case NV097_SET_SEMAPHORE_OFFSET:
			pg->regs[NV_PGRAPH_SEMAPHOREOFFSET] = parameter;
			break;
		case NV097_BACK_END_WRITE_SEMAPHORE_RELEASE: {
			pgraph_update_surface(d, false, true, true);

			//qemu_mutex_unlock(&pg->pgraph_lock);
			//qemu_mutex_lock_iothread();

			uint32_t semaphore_offset = pg->regs[NV_PGRAPH_SEMAPHOREOFFSET];

			xbaddr semaphore_dma_len;
			uint8_t *semaphore_data = (uint8_t*)nv_dma_map(d, pg->dma_semaphore,
				&semaphore_dma_len);
			assert(semaphore_offset < semaphore_dma_len);
			semaphore_data += semaphore_offset;

			stl_le_p((uint32_t*)semaphore_data, parameter);

			//qemu_mutex_lock(&pg->pgraph_lock);
			//qemu_mutex_unlock_iothread();

			break;
		}
		case NV097_SET_ZSTENCIL_CLEAR_VALUE:
			pg->regs[NV_PGRAPH_ZSTENCILCLEARVALUE] = parameter;
			break;

		case NV097_SET_COLOR_CLEAR_VALUE:
			pg->regs[NV_PGRAPH_COLORCLEARVALUE] = parameter;
			break;

		case NV097_CLEAR_SURFACE: {
			pg->clear_surface = parameter;
			pgraph_draw_clear(d);
			break;
		}

		case NV097_SET_CLEAR_RECT_HORIZONTAL:
			pg->regs[NV_PGRAPH_CLEARRECTX] = parameter;
			break;
		case NV097_SET_CLEAR_RECT_VERTICAL:
			pg->regs[NV_PGRAPH_CLEARRECTY] = parameter;
			break;

		CASE_2(NV097_SET_SPECULAR_FOG_FACTOR, 4) :
			slot = (method - NV097_SET_SPECULAR_FOG_FACTOR) / 4;
			pg->regs[NV_PGRAPH_SPECFOGFACTOR0 + slot * 4] = parameter;
			break;

		case NV097_SET_SHADER_CLIP_PLANE_MODE:
			pg->regs[NV_PGRAPH_SHADERCLIPMODE] = parameter;
			break;

		CASE_8(NV097_SET_COMBINER_COLOR_OCW, 4) :
			slot = (method - NV097_SET_COMBINER_COLOR_OCW) / 4;
			pg->regs[NV_PGRAPH_COMBINECOLORO0 + slot * 4] = parameter;
			break;

		case NV097_SET_COMBINER_CONTROL:
			pg->regs[NV_PGRAPH_COMBINECTL] = parameter;
			break;

		case NV097_SET_SHADOW_ZSLOPE_THRESHOLD:
			pg->regs[NV_PGRAPH_SHADOWZSLOPETHRESHOLD] = parameter;
			assert(parameter == 0x7F800000); /* FIXME: Unimplemented */
			break;

		case NV097_SET_SHADER_STAGE_PROGRAM:
			pg->regs[NV_PGRAPH_SHADERPROG] = parameter;
			break;

		case NV097_SET_SHADER_OTHER_STAGE_INPUT:
			pg->regs[NV_PGRAPH_SHADERCTL] = parameter;
			break;

		case NV097_SET_TRANSFORM_EXECUTION_MODE:
			// Test-case : Whiplash
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D], NV_PGRAPH_CSV0_D_MODE,
				GET_MASK(parameter,
					NV097_SET_TRANSFORM_EXECUTION_MODE_MODE));
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_D], NV_PGRAPH_CSV0_D_RANGE_MODE,
				GET_MASK(parameter,
					NV097_SET_TRANSFORM_EXECUTION_MODE_RANGE_MODE));
			break;
		case NV097_SET_TRANSFORM_PROGRAM_CXT_WRITE_EN:
			// Test-case : Whiplash
			pg->enable_vertex_program_write = parameter;
			break;
		case NV097_SET_TRANSFORM_PROGRAM_LOAD:
			assert(parameter < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
			SET_MASK(pg->regs[NV_PGRAPH_CHEOPS_OFFSET],
				NV_PGRAPH_CHEOPS_OFFSET_PROG_LD_PTR, parameter);
			break;
		case NV097_SET_TRANSFORM_PROGRAM_START:
			assert(parameter < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
			SET_MASK(pg->regs[NV_PGRAPH_CSV0_C],
				NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START, parameter);
			break;
		case NV097_SET_TRANSFORM_CONSTANT_LOAD:
			assert(parameter < NV2A_VERTEXSHADER_CONSTANTS);
			SET_MASK(pg->regs[NV_PGRAPH_CHEOPS_OFFSET],
				NV_PGRAPH_CHEOPS_OFFSET_CONST_LD_PTR, parameter);
			NV2A_DPRINTF("load to %d\n", parameter);
			break;

		case NV097_SET_FLAT_SHADE_OP: 
			assert(parameter <= 1);
			// TODO : value & 1 = first/last? vertex selection for glShaderMode(GL_FLAT)
			break;
		default:
			NV2A_DPRINTF("    unhandled  (0x%02x 0x%08x)",
					graphics_class, method);
			break;
		}
		break;
	}

	default:
		NV2A_DPRINTF("Unknown Graphics Class/Method 0x%08X/0x%08X",
						graphics_class, method);
		break;
	}

}

static void pgraph_switch_context(NV2AState *d, unsigned int channel_id)
{
    bool channel_valid =
        d->pgraph.regs[NV_PGRAPH_CTX_CONTROL] & NV_PGRAPH_CTX_CONTROL_CHID;
    unsigned pgraph_channel_id = GET_MASK(d->pgraph.regs[NV_PGRAPH_CTX_USER], NV_PGRAPH_CTX_USER_CHID);
	// Cxbx Note : This isn't present in xqemu / OpenXbox : d->pgraph.pgraph_lock.lock();
    bool valid = channel_valid && pgraph_channel_id == channel_id;
	if (!valid) {
        SET_MASK(d->pgraph.regs[NV_PGRAPH_TRAPPED_ADDR],
                 NV_PGRAPH_TRAPPED_ADDR_CHID, channel_id);

        NV2A_DPRINTF("pgraph switching to ch %d\n", channel_id);

        /* TODO: hardware context switching */
        assert(!(d->pgraph.regs[NV_PGRAPH_DEBUG_3]
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
    while (!(d->pgraph.regs[NV_PGRAPH_FIFO] & NV_PGRAPH_FIFO_ACCESS)) {
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
		NV2A_DPRINTF("d->pgraph method (%d) 0x%08X %s * %d",
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

    /* Now upload the previous vertex_attribute value */
    vertex_attribute->inline_buffer = (float*)g_malloc(NV2A_MAX_BATCH_LENGTH
                                                  * sizeof(float) * 4);
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
	PGRAPHState *pg = &d->pgraph;

	qemu_mutex_init(&pg->pgraph_lock);
	qemu_cond_init(&pg->interrupt_cond);
	qemu_cond_init(&pg->fifo_access_cond);
	qemu_cond_init(&pg->flip_3d);

	if (pg->opengl_enabled) {
		OpenGL_init(pg);
	}
}

void pgraph_destroy(PGRAPHState *pg)
{
	qemu_mutex_destroy(&pg->pgraph_lock);
	qemu_cond_destroy(&pg->interrupt_cond);
	qemu_cond_destroy(&pg->fifo_access_cond);
	qemu_cond_destroy(&pg->flip_3d);

	if (pg->opengl_enabled) {
		OpenGL_destroy(pg);
	}
}

static bool pgraph_get_framebuffer_dirty(PGRAPHState *pg)
{
    bool shape_changed = memcmp(&pg->surface_shape, &pg->last_surface_shape,
                                sizeof(SurfaceShape)) != 0;
    if (!shape_changed || (!pg->surface_shape.color_format
            && !pg->surface_shape.zeta_format)) {
        return false;
    }
    return true;
}

static bool pgraph_get_color_write_enabled(PGRAPHState *pg)
{
	return pg->regs[NV_PGRAPH_CONTROL_0] & (
		NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE
		| NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE
		| NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE
		| NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE);
}

static bool pgraph_get_zeta_write_enabled(PGRAPHState *pg)
{
	return pg->regs[NV_PGRAPH_CONTROL_0] & (
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
    switch (pg->surface_shape.anti_aliasing) {
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
    bool swizzle = (pg->surface_type == NV097_SET_SURFACE_FORMAT_TYPE_SWIZZLE);
    if (swizzle) {
        *width = 1 << pg->surface_shape.log_width;
        *height = 1 << pg->surface_shape.log_height;
    } else {
        *width = pg->surface_shape.clip_width;
        *height = pg->surface_shape.clip_height;
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

// NOTE: Might want to change guint to guint64 for return.
/* functions for texture LRU cache */
static guint texture_key_hash(gconstpointer key)
{
    const TextureKey *k = (const TextureKey *)key;
    uint64_t state_hash = fnv_hash(
        (const uint8_t*)&k->state, sizeof(TextureShape));
    return guint(state_hash ^ k->data_hash);
}
static gboolean texture_key_equal(gconstpointer a, gconstpointer b)
{
    const TextureKey *ak = (const TextureKey *)a, *bk = (const TextureKey *)b;
    return memcmp(&ak->state, &bk->state, sizeof(TextureShape)) == 0
            && ak->data_hash == bk->data_hash;
}
static void texture_key_destroy(gpointer data)
{
    g_free(data);
}

// NOTE: Might want to change guint to guint64 for return.
/* hash and equality for shader cache hash table */
static guint shader_hash(gconstpointer key)
{
    return (guint)fnv_hash((const uint8_t *)key, sizeof(ShaderState));
}
static gboolean shader_equal(gconstpointer a, gconstpointer b)
{
    const ShaderState *as = (const ShaderState *)a, *bs = (const ShaderState *)b;
    return memcmp(as, bs, sizeof(ShaderState)) == 0;
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

static uint64_t fnv_hash(const uint8_t *data, size_t len)
{
    /* 64 bit Fowler/Noll/Vo FNV-1a hash code */
    uint64_t hval = 0xcbf29ce484222325ULL;
    const uint8_t *dp = data;
    const uint8_t *de = data + len;
    while (dp < de) {
        hval ^= (uint64_t) *dp++;
        hval += (hval << 1) + (hval << 4) + (hval << 5) +
            (hval << 7) + (hval << 8) + (hval << 40);
    }

    return hval;
}

static uint64_t fast_hash(const uint8_t *data, size_t len, unsigned int samples)
{
#ifdef __SSE4_2__
    uint64_t h[4] = {len, 0, 0, 0};
    assert(samples > 0);

    if (len < 8 || len % 8) {
        return fnv_hash(data, len);
    }

    assert(len >= 8 && len % 8 == 0);
    const uint64_t *dp = (const uint64_t*)data;
    const uint64_t *de = dp + (len / 8);
    size_t step = len / 8 / samples;
    if (step == 0) step = 1;

    while (dp < de - step * 3) {
        h[0] = __builtin_ia32_crc32di(h[0], dp[step * 0]);
        h[1] = __builtin_ia32_crc32di(h[1], dp[step * 1]);
        h[2] = __builtin_ia32_crc32di(h[2], dp[step * 2]);
        h[3] = __builtin_ia32_crc32di(h[3], dp[step * 3]);
        dp += step * 4;
    }
    if (dp < de - step * 0)
        h[0] = __builtin_ia32_crc32di(h[0], dp[step * 0]);
    if (dp < de - step * 1)
        h[1] = __builtin_ia32_crc32di(h[1], dp[step * 1]);
    if (dp < de - step * 2)
        h[2] = __builtin_ia32_crc32di(h[2], dp[step * 2]);

    return h[0] + (h[1] << 10) + (h[2] << 21) + (h[3] << 32);
#else
    return fnv_hash(data, len);
#endif
}
