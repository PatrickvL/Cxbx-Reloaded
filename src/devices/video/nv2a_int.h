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

#include "xbox_types.h" // For xbox::addr_xt

#include "qemu-thread.h" // For qemu_mutex, etc

#ifdef USE_SHADER_CACHE
#include "glib_compat.h" // For GHashTable, g_hash_table_new, g_hash_table_lookup, g_hash_table_insert
#endif

#include "swizzle.h"

#include "nv2a_debug.h" // For HWADDR_PRIx, NV2A_DPRINTF, NV2A_GL_DPRINTF, etc.
#include "nv2a_regs.h" // For NV2A_MAX_TEXTURES, etc


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
#define CASE_128(v, step) CASE_64(v, step) : CASE_64(v + (step) * 64, step)
#define CASE_256(v, step) CASE_128(v, step) : CASE_128(v + (step) * 128, step)

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

	bool needs_conversion;
	uint8_t *converted_buffer;
	unsigned int converted_elements;
	unsigned int converted_size;
	unsigned int converted_count;

	float *inline_buffer;

	int32_t gl_count;
	int32_t gl_type;
	int32_t gl_normalize;
} VertexAttribute;

typedef struct Surface {
	bool draw_dirty;
	bool buffer_dirty;
	bool write_enabled_cache;
	unsigned int pitch;

	xbox::addr_xt offset;
} Surface;

typedef struct SurfaceShape {
	unsigned int z_format;
	unsigned int color_format;
	unsigned int zeta_format;
	unsigned int log_width, log_height;
	unsigned int clip_x, clip_y;
	unsigned int clip_width, clip_height;
	unsigned int anti_aliasing;
} SurfaceShape;

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

typedef struct KelvinState {
	xbox::addr_xt object_instance;
} KelvinState;

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
	unsigned int surface_type;
	SurfaceShape surface_shape;
	SurfaceShape last_surface_shape;

	xbox::addr_xt dma_semaphore;

	xbox::addr_xt dma_report;
	unsigned int zpass_pixel_count_result;

	unsigned int primitive_mode;

	uint32_t clear_surface_flags; // NV097_CLEAR_SURFACE parameter (Z/STENCIL/COLOR mask)

	uint32_t program_data[NV2A_MAX_TRANSFORM_PROGRAM_LENGTH][VSH_TOKEN_SIZE]; // XFPR RAM mirror: NV2A Transform Program RAM (on-chip XF SRAM, 136 × 92-bit instructions in 128-bit containers)

	uint32_t vsh_constants[NV2A_VERTEXSHADER_CONSTANTS][4]; // XFCTX RAM mirror: NV2A Transform Context RAM (on-chip XF SRAM, 192 × float4)
	bool vsh_constants_dirty[NV2A_VERTEXSHADER_CONSTANTS];

	/* lighting constant arrays */
	uint32_t ltctxa[NV2A_LTCTXA_COUNT][4];
	bool ltctxa_dirty[NV2A_LTCTXA_COUNT];
	uint32_t ltctxb[NV2A_LTCTXB_COUNT][4];
	bool ltctxb_dirty[NV2A_LTCTXB_COUNT];
	uint32_t ltc1[NV2A_LTC1_COUNT][4];
	bool ltc1_dirty[NV2A_LTC1_COUNT];

	// should figure out where these are in lighting context
	float light_infinite_half_vector[NV2A_MAX_LIGHTS][3];
	float light_infinite_direction[NV2A_MAX_LIGHTS][3];
	float light_local_position[NV2A_MAX_LIGHTS][3];
	float light_local_attenuation[NV2A_MAX_LIGHTS][3];

	float point_params[8]; // NV097_SET_POINT_PARAMS attenuation coefficients
	float line_width;      // NV097_SET_LINE_WIDTH (float, pixels)

	VertexAttribute vertex_attributes[NV2A_VERTEXSHADER_ATTRIBUTES];

	unsigned int inline_array_length;
	uint32_t inline_array[NV2A_MAX_BATCH_LENGTH];
	unsigned int inline_elements_length;
	uint16_t inline_elements[NV2A_MAX_BATCH_LENGTH]; // Cxbx-Reloaded TODO : Restore uint32_t once D3D11_draw_inline_elements can using that

	unsigned int inline_buffer_length;

	unsigned int draw_arrays_length;
	unsigned int draw_arrays_max_count;

	/* FIXME: Unknown size, possibly endless, 1000 will do for now */
	int32_t gl_draw_arrays_start[1000];
	int32_t gl_draw_arrays_count[1000];

	uint32_t regs[NV_PGRAPH_SIZE]; // TODO : union
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
		uint32_t regs[NV_PMC_SIZE]; // Not in xqemu/openxbox? TODO : union
    } pmc;

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
		uint32_t regs[_NV_PFIFO_SIZE]; // TODO : union
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
		uint32_t regs[NV_PVIDEO_SIZE]; // TODO : union
    } pvideo;

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
        uint32_t numerator;
        uint32_t denominator;
        uint32_t alarm_time;
		uint32_t regs[NV_PTIMER_SIZE]; // Not in xqemu/openxbox? TODO : union
    } ptimer;

    struct {
		uint32_t regs[NV_PFB_SIZE]; // TODO : union
    } pfb;

    struct PGRAPHState pgraph;

    struct {
        uint32_t pending_interrupts;
        uint32_t enabled_interrupts;
        hwaddr start;
		uint32_t regs[NV_PCRTC_SIZE]; // Not in xqemu/openxbox? TODO : union
    } pcrtc;

    struct {
        uint32_t core_clock_coeff;
        uint64_t core_clock_freq;
        uint32_t memory_clock_coeff;
        uint32_t video_clock_coeff;
		uint32_t regs[NV_PRAMDAC_SIZE]; // Not in xqemu/openxbox? TODO : union
    } pramdac;

	// PRMCIO (Actually the VGA controller)
	struct {
		uint8_t cr_index;
		uint8_t cr[256]; /* CRT registers */
	} prmcio; // Not in xqemu/openxbox?
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
