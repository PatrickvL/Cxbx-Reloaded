// ******************************************************************
// *  NV2A Method -> Register Dispatch Table
// *
// *  Data-driven NV097 method -> register/array write mapping.
// *  Covers DIRECT_COPY, DIRECT_COPY_RANGED, and MASKED_WRITE patterns.
// *  Methods requiring custom logic fall through (all-zero entry).
// *
// *  Table indexed by method >> 2 (dword index into method space).
// *  Each entry: { target:5, dirty_group:3, shift:7, identity:1,
// *                reg_index:11, width:5 } — 4 bytes.
// *    all zero         -> not data-driven, use switch fallback
// *    width == 0       -> full 32-bit copy: base[reg_index] = parameter
// *    width > 0        -> masked write using mask = ((1<<width)-1) << shift
// *      identity == 0  -> value-shift: SET_MASK(base[ri], mask, parameter)
// *      identity == 1  -> identity:    base[ri] = (old & ~mask) | (param & mask)
// *    target: NV097Target — which array to address:
// *      0 = pg->regs[] (PGRAPH MMIO shadow registers)
// *      1 = Cheops XF SRAM: xfctx (transform context)
// *      2-4 = Cheops XF SRAM: ltctxa, ltctxb, ltc1 (lighting)
// *      5 = pg->light[] (non-SRAM light geometry)
// *      6-7 = Indirect SRAM writes via auto-incrementing CHEOPS_OFFSET
// *    dirty_group: NV2ADirtyGroup to set on write (0 = none)
// *    reg_index: pre-RI'd dword index (byte_offset >> 2) into target array
// *
// *  For a future PFIFO-in-shader: this table fits in an 8KB structured
// *  buffer (2048 entries x 4 bytes) and the dispatch is branchless for
// *  all data-driven methods.  GPU branchless form:
// *    uint ri = entry.reg_index;
// *    uint mask = (width > 0) ? (((1u << width) - 1) << shift) : ~0u;
// *    uint val = (width == 0) ? param
// *             : identity ? (param & mask) : ((param << shift) & mask);
// *    base[ri] = (base[ri] & ~mask) | val;
// ******************************************************************
#pragma once

#include "nv2a_int.h"
#include "nv2a_regs.h"

#include <cstring> // memset

// Target array identifiers for NV097MethodEntry.target.
// Determines which uint32_t array the dispatch writes to.
enum NV097Target : uint8_t {
	NV097_TARGET_PGRAPH = 0,   // pg->regs[] — PGRAPH MMIO shadow registers
	NV097_TARGET_XFCTX  = 1,   // (uint32_t*)pg->xf.xfctx — transform context (VS constants)
	NV097_TARGET_LTCTXA = 2,   // (uint32_t*)pg->xf.ltctxa — lighting context A
	NV097_TARGET_LTCTXB = 3,   // (uint32_t*)pg->xf.ltctxb — lighting context B
	NV097_TARGET_LTC1   = 4,   // (uint32_t*)pg->xf.ltc1 — lighting constants 1
	NV097_TARGET_LIGHT  = 5,   // (uint32_t*)pg->light — light geometry (non-SRAM)
	// --- Indirect targets: SRAM writes via auto-incrementing CHEOPS_OFFSET pointer ---
	// reg_index encodes the float4 component (0..3); the SRAM row is read at
	// runtime from CHEOPS_OFFSET (CONST_LD_PTR / PROG_LD_PTR).
	NV097_TARGET_XFCTX_INDIRECT = 6, // xfctx via CONST_LD_PTR, per-constant dirty
	NV097_TARGET_XFPR_INDIRECT  = 7, // xfpr  via PROG_LD_PTR
	NV097_TARGET_COUNT
};

struct NV097MethodEntry {
	// --- Destination ---
	uint8_t target : 5;       // NV097Target: which array to address (0 = pg->regs[])
	uint8_t dirty_group : 3;  // NV2ADirtyGroup: dirty flag to set on write (0 = none)
	// --- Mask write parameters ---
	uint8_t shift : 7;        // lowest bit position of the mask field
	uint8_t identity : 1;     // 1 = identity copy (param bits at mask position), 0 = value-shift
	// --- Register index and width ---
	uint16_t reg_index : 11;  // pre-RI'd dword index into target array (byte_offset >> 2)
	uint16_t width : 5;       // 0 = full 32-bit copy, >0 = masked write (width bits, max 31)
};

static_assert(sizeof(NV097MethodEntry) == 4, "NV097MethodEntry must be 4 bytes");

// Method space is 0x0000..0x1FFF = 2048 dword slots
#define NV097_METHOD_TABLE_SIZE 2048

// Helper: count the width (number of consecutive 1-bits) of a contiguous mask
static inline uint8_t nv097_mask_width(uint32_t mask)
{
	if (mask == 0) return 0;
	mask >>= (ffs(mask) - 1);
	uint8_t w = 0;
	while (mask & 1) { w++; mask >>= 1; }
	return w;
}

// Helper to register a direct-copy method (full 32-bit write)
#define NV097_REG_DIRECT(tbl, nv097_method, pg_reg, dirty) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.reg_index = (uint16_t)((pg_reg) >> 2); \
		e.shift = 0; \
		e.identity = 0; \
		e.width = 0; \
		e.dirty_group = (dirty); \
	} while (0)

// Helper to register a direct-copy ranged method (variable stride, N slots)
#define NV097_REG_DIRECT_RANGE(tbl, nv097_base, stride, count, pg_base, dirty) \
	do { \
		for (int _i = 0; _i < (count); _i++) { \
			auto &e = (tbl)[((nv097_base) + _i * (stride)) >> 2]; \
			e.reg_index = (uint16_t)(((pg_base) + _i * 4) >> 2); \
			e.shift = 0; \
			e.identity = 0; \
			e.width = 0; \
			e.dirty_group = (dirty); \
		} \
	} while (0)

// Helper to register a value-shift masked write: SET_MASK(reg, mask, parameter)
// Parameter is a small value at bit 0; SET_MASK shifts it into the mask position.
#define NV097_REG_MASKED(tbl, nv097_method, pg_reg, pg_mask, dirty) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.reg_index = (uint16_t)((pg_reg) >> 2); \
		e.shift = (uint8_t)(ffs(pg_mask) - 1); \
		e.identity = 0; \
		e.width = nv097_mask_width(pg_mask); \
		e.dirty_group = (dirty); \
	} while (0)

// Helper to register an identity-copy masked write:
// reg = (reg & ~mask) | (parameter & mask)
// Parameter bits are already at the mask position (no shift needed).
#define NV097_REG_MASKED_IDENTITY(tbl, nv097_method, pg_reg, pg_mask, dirty) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.reg_index = (uint16_t)((pg_reg) >> 2); \
		e.shift = (uint8_t)(ffs(pg_mask) - 1); \
		e.identity = 1; \
		e.width = nv097_mask_width(pg_mask); \
		e.dirty_group = (dirty); \
	} while (0)

// Helper to register a direct-copy write to an XF SRAM bank.
// flat_index is row * 4 + component (index into the uint32_t view of the 2D array).
#define NV097_REG_SRAM_DIRECT(tbl, nv097_method, tgt, flat_idx, dirty) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.target = (tgt); \
		e.reg_index = (uint16_t)(flat_idx); \
		e.shift = 0; \
		e.identity = 0; \
		e.width = 0; \
		e.dirty_group = (dirty); \
	} while (0)

// Build the dispatch table at startup (called once)
static NV097MethodEntry nv097_method_table[NV097_METHOD_TABLE_SIZE];

static void nv097_init_method_table()
{
	// Zero-init: all-zero entry = not data-driven (fallback to switch)
	memset(nv097_method_table, 0, sizeof(nv097_method_table));

	auto &t = nv097_method_table; // shorthand

	// ===== DIRECT_COPY: pg->regs[RI(reg)] = parameter =====
	NV097_REG_DIRECT(t, NV097_SET_SURFACE_CLIP_HORIZONTAL,    NV_PGRAPH_SURFACECLIPX,          NV2A_DIRTY_SURFACE);
	NV097_REG_DIRECT(t, NV097_SET_SURFACE_CLIP_VERTICAL,      NV_PGRAPH_SURFACECLIPY,          NV2A_DIRTY_SURFACE);
	NV097_REG_DIRECT(t, NV097_SET_SURFACE_FORMAT,             NV_PGRAPH_SURFACEFORMAT,         NV2A_DIRTY_SURFACE);
	NV097_REG_DIRECT(t, NV097_SET_SURFACE_PITCH,              NV_PGRAPH_DMA_PITCH,             NV2A_DIRTY_SURFACE);
	NV097_REG_DIRECT(t, NV097_SET_SURFACE_COLOR_OFFSET,       NV_PGRAPH_BOFFSET3,              NV2A_DIRTY_SURFACE);
	NV097_REG_DIRECT(t, NV097_SET_SURFACE_ZETA_OFFSET,        NV_PGRAPH_BOFFSET4,              NV2A_DIRTY_SURFACE);
	NV097_REG_DIRECT(t, NV097_SET_FOG_COLOR,                  NV_PGRAPH_FOGCOLOR,              NV2A_DIRTY_BLEND);
	NV097_REG_DIRECT(t, NV097_SET_SHADOW_ZSLOPE_THRESHOLD,    NV_PGRAPH_SHADOWZSLOPETHRESHOLD, NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_DIRECT(t, NV097_SET_COMBINER_SPECULAR_FOG_CW0,  NV_PGRAPH_COMBINESPECFOG0,       NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT(t, NV097_SET_COMBINER_SPECULAR_FOG_CW1,  NV_PGRAPH_COMBINESPECFOG1,       NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT(t, NV097_SET_BLEND_COLOR,                NV_PGRAPH_BLENDCOLOR,            NV2A_DIRTY_BLEND);
	NV097_REG_DIRECT(t, NV097_SET_POLYGON_OFFSET_SCALE_FACTOR,NV_PGRAPH_ZOFFSETFACTOR,         NV2A_DIRTY_RASTERIZER);
	NV097_REG_DIRECT(t, NV097_SET_POLYGON_OFFSET_BIAS,        NV_PGRAPH_ZOFFSETBIAS,           NV2A_DIRTY_RASTERIZER);
	NV097_REG_DIRECT(t, NV097_SET_CLIP_MIN,                   NV_PGRAPH_ZCLIPMIN,              NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_DIRECT(t, NV097_SET_CLIP_MAX,                   NV_PGRAPH_ZCLIPMAX,              NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_DIRECT(t, NV097_SET_SEMAPHORE_OFFSET,           NV_PGRAPH_SEMAPHOREOFFSET,       NV2A_DIRTY_NONE);
	NV097_REG_DIRECT(t, NV097_SET_ZSTENCIL_CLEAR_VALUE,       NV_PGRAPH_ZSTENCILCLEARVALUE,    NV2A_DIRTY_NONE);
	NV097_REG_DIRECT(t, NV097_SET_COLOR_CLEAR_VALUE,          NV_PGRAPH_COLORCLEARVALUE,       NV2A_DIRTY_NONE);
	NV097_REG_DIRECT(t, NV097_SET_CLEAR_RECT_HORIZONTAL,      NV_PGRAPH_CLEARRECTX,            NV2A_DIRTY_NONE);
	NV097_REG_DIRECT(t, NV097_SET_CLEAR_RECT_VERTICAL,        NV_PGRAPH_CLEARRECTY,            NV2A_DIRTY_NONE);
	NV097_REG_DIRECT(t, NV097_SET_SHADER_CLIP_PLANE_MODE,     NV_PGRAPH_SHADERCLIPMODE,        NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT(t, NV097_SET_COMBINER_CONTROL,           NV_PGRAPH_COMBINECTL,            NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT(t, NV097_SET_SHADER_STAGE_PROGRAM,       NV_PGRAPH_SHADERPROG,            NV2A_DIRTY_SHADER);

	// ===== DIRECT_COPY_RANGED (stride 4) =====
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_ALPHA_ICW,     4, 8, NV_PGRAPH_COMBINEALPHAI0, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_WINDOW_CLIP_HORIZONTAL, 4, 8, NV_PGRAPH_WINDOWCLIPX0,   NV2A_DIRTY_RASTERIZER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_WINDOW_CLIP_VERTICAL,   4, 8, NV_PGRAPH_WINDOWCLIPY0,   NV2A_DIRTY_RASTERIZER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_FACTOR0,       4, 8, NV_PGRAPH_COMBINEFACTOR0, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_FACTOR1,       4, 8, NV_PGRAPH_COMBINEFACTOR1, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_ALPHA_OCW,     4, 8, NV_PGRAPH_COMBINEALPHAO0, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_COLOR_ICW,     4, 8, NV_PGRAPH_COMBINECOLORI0, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_COLOR_OCW,     4, 8, NV_PGRAPH_COMBINECOLORO0, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_SPECULAR_FOG_FACTOR,    4, 2, NV_PGRAPH_SPECFOGFACTOR0, NV2A_DIRTY_SHADER);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_EYE_VECTOR,             4, 3, NV_PGRAPH_EYEVEC0,        NV2A_DIRTY_NONE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COLOR_KEY_COLOR,        4, 4, NV_PGRAPH_COLORKEYCOLOR0, NV2A_DIRTY_TEXTURE);

	// ===== DIRECT_COPY_RANGED (stride 64, texture methods) =====
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_OFFSET,        64, 4, NV_PGRAPH_TEXOFFSET0,    NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_FORMAT,        64, 4, NV_PGRAPH_TEXFMT0,       NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_ADDRESS,       64, 4, NV_PGRAPH_TEXADDRESS0,   NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_CONTROL0,      64, 4, NV_PGRAPH_TEXCTL0_0,     NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_CONTROL1,      64, 4, NV_PGRAPH_TEXCTL1_0,     NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_FILTER,        64, 4, NV_PGRAPH_TEXFILTER0,    NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_IMAGE_RECT,    64, 4, NV_PGRAPH_TEXIMAGERECT0, NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_PALETTE,       64, 4, NV_PGRAPH_TEXPALETTE0,   NV2A_DIRTY_TEXTURE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_BORDER_COLOR,  64, 4, NV_PGRAPH_BORDERCOLOR0,  NV2A_DIRTY_TEXTURE);

	// ===== DIRECT_COPY_RANGED (stride 64, bump env — stages 1-3 only) =====
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_SET_BUMP_ENV_SCALE  + 64, 64, 3, NV_PGRAPH_BUMPSCALE1,  NV2A_DIRTY_NONE);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_SET_BUMP_ENV_OFFSET + 64, 64, 3, NV_PGRAPH_BUMPOFFSET1, NV2A_DIRTY_NONE);

	// ===== DIRECT_COPY: bump env matrix — stages 1-3, 4 entries per stage =====
	// NV_PGRAPH_BUMPMATxy layout: component[stage-1]*4 offset per stage.
	// Stage 0 has no bump env; its methods hit the switch (assert).
	for (int s = 1; s <= 3; s++) {
		NV097_REG_DIRECT(t, NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x0 + s * 64, NV_PGRAPH_BUMPMAT00 + (s - 1) * 4, NV2A_DIRTY_NONE);
		NV097_REG_DIRECT(t, NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x4 + s * 64, NV_PGRAPH_BUMPMAT01 + (s - 1) * 4, NV2A_DIRTY_NONE);
		NV097_REG_DIRECT(t, NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0x8 + s * 64, NV_PGRAPH_BUMPMAT10 + (s - 1) * 4, NV2A_DIRTY_NONE);
		NV097_REG_DIRECT(t, NV097_SET_TEXTURE_SET_BUMP_ENV_MAT + 0xC + s * 64, NV_PGRAPH_BUMPMAT11 + (s - 1) * 4, NV2A_DIRTY_NONE);
	}

	// ===== MASKED_WRITE (value-shift): SET_MASK(regs[RI(reg)], mask, parameter) =====
	NV097_REG_MASKED(t, NV097_SET_FLIP_READ,                    NV_PGRAPH_SURFACE,     NV_PGRAPH_SURFACE_READ_3D,                    NV2A_DIRTY_SURFACE);
	NV097_REG_MASKED(t, NV097_SET_FLIP_WRITE,                   NV_PGRAPH_SURFACE,     NV_PGRAPH_SURFACE_WRITE_3D,                   NV2A_DIRTY_SURFACE);
	NV097_REG_MASKED(t, NV097_SET_FLIP_MODULO,                  NV_PGRAPH_SURFACE,     NV_PGRAPH_SURFACE_MODULO_3D,                  NV2A_DIRTY_SURFACE);
	NV097_REG_MASKED(t, NV097_SET_FOG_ENABLE,                   NV_PGRAPH_CONTROL_3,   NV_PGRAPH_CONTROL_3_FOGENABLE,                NV2A_DIRTY_BLEND);
	NV097_REG_MASKED(t, NV097_SET_WINDOW_CLIP_TYPE,             NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_WINDOWCLIPTYPE,         NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_ALPHA_TEST_ENABLE,            NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ALPHATESTENABLE,          NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_BLEND_ENABLE,                 NV_PGRAPH_BLEND,       NV_PGRAPH_BLEND_EN,                           NV2A_DIRTY_BLEND);
	NV097_REG_MASKED(t, NV097_SET_CULL_FACE_ENABLE,             NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_CULLENABLE,             NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_DEPTH_TEST_ENABLE,            NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ZENABLE,                  NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_DITHER_ENABLE,                NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_DITHERENABLE,             NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_LIGHTING_ENABLE,              NV_PGRAPH_CSV0_C,      NV_PGRAPH_CSV0_C_LIGHTING,                    NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_SKIN_MODE,                    NV_PGRAPH_CSV0_D,      NV_PGRAPH_CSV0_D_SKIN,                        NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_TEST_ENABLE,          NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE,      NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_POLY_OFFSET_POINT_ENABLE,     NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POFFSETPOINTENABLE,     NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_POLY_OFFSET_LINE_ENABLE,      NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POFFSETLINEENABLE,      NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_POLY_OFFSET_FILL_ENABLE,      NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POFFSETFILLENABLE,      NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_DEPTH_MASK,                   NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ZWRITEENABLE,             NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_LINE_SMOOTH_ENABLE,           NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_LINESMOOTHENABLE,       NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_POLY_SMOOTH_ENABLE,           NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POLYSMOOTHENABLE,       NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_ALPHA_FUNC,                   NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ALPHAFUNC,                NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_ALPHA_REF,                    NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ALPHAREF,                 NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_DEPTH_FUNC,                   NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ZFUNC,                    NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_MASK,                 NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_MASK_WRITE,       NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_FUNC,                 NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_FUNC,             NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_FUNC_REF,             NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_REF,              NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_FUNC_MASK,            NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_MASK_READ,        NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_NORMALIZATION_ENABLE,         NV_PGRAPH_CSV0_C,      NV_PGRAPH_CSV0_C_NORMALIZATION_ENABLE,        NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_LIGHT_ENABLE_MASK,            NV_PGRAPH_CSV0_D,      NV_PGRAPH_CSV0_D_LIGHTS,                      NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_TEXGEN_VIEW_MODEL,            NV_PGRAPH_CSV0_D,      NV_PGRAPH_CSV0_D_TEXGEN_REF,                  NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_LOGIC_OP_ENABLE,              NV_PGRAPH_BLEND,       NV_PGRAPH_BLEND_LOGICOP_ENABLE,               NV2A_DIRTY_BLEND);
	NV097_REG_MASKED(t, NV097_SET_LOGIC_OP,                     NV_PGRAPH_BLEND,       NV_PGRAPH_BLEND_LOGICOP,                      NV2A_DIRTY_BLEND);
	NV097_REG_MASKED(t, NV097_SET_SPECULAR_ENABLE,              NV_PGRAPH_CSV0_C,      NV_PGRAPH_CSV0_C_SPECULAR_ENABLE,             NV2A_DIRTY_SHADER);
	NV097_REG_MASKED(t, NV097_SET_POINT_SMOOTH_ENABLE,          NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POINTSMOOTHENABLE,      NV2A_DIRTY_RASTERIZER);
	NV097_REG_MASKED(t, NV097_SET_PROVOKING_VERTEX,             NV_PGRAPH_CONTROL_3,   NV_PGRAPH_CONTROL_3_PROVOKING_VERTEX,         NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_SHADOW_DEPTH_FUNC,            NV_PGRAPH_SHADOWCTL,   NV_PGRAPH_SHADOWCTL_SHADOW_ZFUNC,             NV2A_DIRTY_DEPTH_STENCIL);
	NV097_REG_MASKED(t, NV097_SET_ANTI_ALIASING_CONTROL,        NV_PGRAPH_ANTIALIASING, NV_PGRAPH_ANTIALIASING_ENABLE,               NV2A_DIRTY_RASTERIZER);

	// ===== MASKED_WRITE (transform pipeline — table writes the register, switch still runs for side-effects) =====
	NV097_REG_MASKED(t, NV097_SET_TRANSFORM_PROGRAM_LOAD,       NV_PGRAPH_CHEOPS_OFFSET, NV_PGRAPH_CHEOPS_OFFSET_PROG_LD_PTR,        NV2A_DIRTY_NONE);
	NV097_REG_MASKED(t, NV097_SET_TRANSFORM_PROGRAM_START,      NV_PGRAPH_CSV0_C,        NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START,      NV2A_DIRTY_PROGRAM);
	NV097_REG_MASKED(t, NV097_SET_TRANSFORM_CONSTANT_LOAD,      NV_PGRAPH_CHEOPS_OFFSET, NV_PGRAPH_CHEOPS_OFFSET_CONST_LD_PTR,       NV2A_DIRTY_NONE);

	// ===== MASKED_WRITE (value-shift, mask at bit 0 — equivalent to identity) =====
	NV097_REG_MASKED(t, NV097_SET_DOT_RGBMAPPING,               NV_PGRAPH_SHADERCTL,   NV_PGRAPH_SHADERCTL_DOT_RGBMAPPING,           NV2A_DIRTY_SHADER);

	// ===== MASKED_WRITE (identity-copy: param bits already at mask position) =====
	NV097_REG_MASKED_IDENTITY(t, NV097_SET_SHADER_OTHER_STAGE_INPUT, NV_PGRAPH_SHADERCTL, NV_PGRAPH_SHADERCTL_OTHER_STAGE_INPUT,     NV2A_DIRTY_SHADER);

	// ===== DIRECT SRAM WRITES: xfctx (matrices, viewport, eye, fog) =====

	// Projection matrix: 4×4 = 16 methods → xfctx rows PMAT0..PMAT0+3
	for (int i = 0; i < 16; i++) {
		int row = NV_IGRAPH_XF_XFCTX_PMAT0 + i / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_PROJECTION_MATRIX + i * 4,
			NV097_TARGET_XFCTX, row * 4 + (i % 4), NV2A_DIRTY_NONE);
	}

	// Model-view matrices: 4 matrices × 4×4 = 64 methods
	// SRAM layout: MMAT0 + matnum*8, 4 rows per matrix (stride 8 rows between matrices)
	for (int i = 0; i < 64; i++) {
		int matnum = i / 16;
		int entry = i % 16;
		int row = NV_IGRAPH_XF_XFCTX_MMAT0 + matnum * 8 + entry / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_MODEL_VIEW_MATRIX + i * 4,
			NV097_TARGET_XFCTX, row * 4 + (entry % 4), NV2A_DIRTY_NONE);
	}

	// Inverse model-view matrices: same layout as model-view at IMMAT0
	for (int i = 0; i < 64; i++) {
		int matnum = i / 16;
		int entry = i % 16;
		int row = NV_IGRAPH_XF_XFCTX_IMMAT0 + matnum * 8 + entry / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_INVERSE_MODEL_VIEW_MATRIX + i * 4,
			NV097_TARGET_XFCTX, row * 4 + (entry % 4), NV2A_DIRTY_NONE);
	}

	// Composite matrix: 4×4 = 16 methods → xfctx rows CMAT0..CMAT0+3
	for (int i = 0; i < 16; i++) {
		int row = NV_IGRAPH_XF_XFCTX_CMAT0 + i / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_COMPOSITE_MATRIX + i * 4,
			NV097_TARGET_XFCTX, row * 4 + (i % 4), NV2A_DIRTY_NONE);
	}

	// Texture matrices: 4 textures × 4×4 = 64 methods
	// SRAM layout: T0MAT + tex*8, 4 rows per texture (stride 8 between textures)
	for (int i = 0; i < 64; i++) {
		int tex = i / 16;
		int entry = i % 16;
		int row = NV_IGRAPH_XF_XFCTX_T0MAT + tex * 8 + entry / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_TEXTURE_MATRIX + i * 4,
			NV097_TARGET_XFCTX, row * 4 + (entry % 4), NV2A_DIRTY_NONE);
	}

	// Texgen planes S,T,R,Q: 4 textures × 4×4 = 64 methods
	// SRAM layout: TG0MAT + tex*8, same stride as texture matrices
	for (int i = 0; i < 64; i++) {
		int tex = i / 16;
		int entry = i % 16;
		int row = NV_IGRAPH_XF_XFCTX_TG0MAT + tex * 8 + entry / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_TEXGEN_PLANE_S + i * 4,
			NV097_TARGET_XFCTX, row * 4 + (entry % 4), NV2A_DIRTY_NONE);
	}

	// Fog plane: 4 components → xfctx row FOG
	for (int i = 0; i < 4; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_FOG_PLANE + i * 4,
			NV097_TARGET_XFCTX, NV_IGRAPH_XF_XFCTX_FOG * 4 + i, NV2A_DIRTY_NONE);

	// Viewport offset: 4 components → xfctx row VPOFF
	for (int i = 0; i < 4; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_VIEWPORT_OFFSET + i * 4,
			NV097_TARGET_XFCTX, NV_IGRAPH_XF_XFCTX_VPOFF * 4 + i, NV2A_DIRTY_NONE);

	// Eye position: 4 components → xfctx row EYEP
	for (int i = 0; i < 4; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_EYE_POSITION + i * 4,
			NV097_TARGET_XFCTX, NV_IGRAPH_XF_XFCTX_EYEP * 4 + i, NV2A_DIRTY_NONE);

	// Viewport scale: 4 components → xfctx row VPSCL
	for (int i = 0; i < 4; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_VIEWPORT_SCALE + i * 4,
			NV097_TARGET_XFCTX, NV_IGRAPH_XF_XFCTX_VPSCL * 4 + i, NV2A_DIRTY_NONE);

	// ===== DIRECT SRAM WRITES: ltctxa (fog, ambient, material, per-light spot) =====

	// Fog params: 3 components → ltctxa row FOG_K (switch also writes regs[FOGPARAM0..2])
	for (int i = 0; i < 3; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_FOG_PARAMS + i * 4,
			NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_FOG_K * 4 + i, NV2A_DIRTY_NONE);

	// Scene ambient color: 3 components → ltctxa row FR_AMB
	for (int i = 0; i < 3; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_SCENE_AMBIENT_COLOR + i * 4,
			NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_FR_AMB * 4 + i, NV2A_DIRTY_NONE);

	// Back scene ambient color: 3 components → ltctxa row BR_AMB
	for (int i = 0; i < 3; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_SCENE_AMBIENT_COLOR + i * 4,
			NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_BR_AMB * 4 + i, NV2A_DIRTY_NONE);

	// Material emission: 3 components → ltctxa row CM_COL
	for (int i = 0; i < 3; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_MATERIAL_EMISSION + i * 4,
			NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_CM_COL * 4 + i, NV2A_DIRTY_NONE);

	// Material alpha: 1 component → ltctxa row CM_COL, component 3
	NV097_REG_SRAM_DIRECT(t, NV097_SET_MATERIAL_ALPHA,
		NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_CM_COL * 4 + 3, NV2A_DIRTY_NONE);

	// Back material emission: 3 components → ltctxa row BCM_COL
	for (int i = 0; i < 3; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_MATERIAL_EMISSIONR + i * 4,
			NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_BCM_COL * 4 + i, NV2A_DIRTY_NONE);

	// Back material alpha: 1 component → ltctxa row BCM_COL, component 3
	NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_MATERIAL_ALPHA,
		NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_BCM_COL * 4 + 3, NV2A_DIRTY_NONE);

	// Eye direction: 3 components → ltctxa row EYED
	for (int i = 0; i < 3; i++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_EYE_DIRECTION + i * 4,
			NV097_TARGET_LTCTXA, NV_IGRAPH_XF_LTCTXA_EYED * 4 + i, NV2A_DIRTY_NONE);

	// Per-light spot falloff: 3 components × 8 lights → ltctxa rows L0_K + light*2
	for (int s = 0; s < 8; s++)
		for (int i = 0; i < 3; i++)
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_SPOT_FALLOFF + i * 4 + s * 0x80,
				NV097_TARGET_LTCTXA, (NV_IGRAPH_XF_LTCTXA_L0_K + s * 2) * 4 + i, NV2A_DIRTY_NONE);

	// Per-light spot direction: 4 components × 8 lights → ltctxa rows L0_SPT + light*2
	for (int s = 0; s < 8; s++)
		for (int i = 0; i < 4; i++)
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_SPOT_DIRECTION + i * 4 + s * 0x80,
				NV097_TARGET_LTCTXA, (NV_IGRAPH_XF_LTCTXA_L0_SPT + s * 2) * 4 + i, NV2A_DIRTY_NONE);

	// ===== DIRECT SRAM WRITES: ltctxb (per-light colors, front + back) =====

	// Front light colors: 8 lights × 3 properties (ambient, diffuse, specular) × 3 components
	for (int s = 0; s < 8; s++) {
		for (int i = 0; i < 3; i++) {
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_AMBIENT_COLOR  + i * 4 + s * 0x80,
				NV097_TARGET_LTCTXB, (NV_IGRAPH_XF_LTCTXB_L0_AMB + s * 6) * 4 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_DIFFUSE_COLOR  + i * 4 + s * 0x80,
				NV097_TARGET_LTCTXB, (NV_IGRAPH_XF_LTCTXB_L0_DIF + s * 6) * 4 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_SPECULAR_COLOR + i * 4 + s * 0x80,
				NV097_TARGET_LTCTXB, (NV_IGRAPH_XF_LTCTXB_L0_SPC + s * 6) * 4 + i, NV2A_DIRTY_NONE);
		}
	}

	// Back light colors: 8 lights × 3 properties × 3 components
	for (int s = 0; s < 8; s++) {
		for (int i = 0; i < 3; i++) {
			NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_LIGHT_AMBIENT_COLOR  + i * 4 + s * 0x40,
				NV097_TARGET_LTCTXB, (NV_IGRAPH_XF_LTCTXB_L0_BAMB + s * 6) * 4 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_LIGHT_DIFFUSE_COLOR  + i * 4 + s * 0x40,
				NV097_TARGET_LTCTXB, (NV_IGRAPH_XF_LTCTXB_L0_BDIF + s * 6) * 4 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_LIGHT_SPECULAR_COLOR + i * 4 + s * 0x40,
				NV097_TARGET_LTCTXB, (NV_IGRAPH_XF_LTCTXB_L0_BSPC + s * 6) * 4 + i, NV2A_DIRTY_NONE);
		}
	}

	// ===== DIRECT SRAM WRITES: ltc1 (specular params, per-light range) =====

	// Specular params: 6 components → ltc1 rows l0..l0+1 (2 rows × 3+3 layout)
	for (int i = 0; i < 6; i++) {
		int row = NV_IGRAPH_XF_LTC1_l0 + i / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_SPECULAR_PARAMS + i * 4,
			NV097_TARGET_LTC1, row * 4 + (i % 4), NV2A_DIRTY_NONE);
	}

	// Back specular params: 6 components → ltc1 rows Bl0..Bl0+1
	for (int i = 0; i < 6; i++) {
		int row = NV_IGRAPH_XF_LTC1_Bl0 + i / 4;
		NV097_REG_SRAM_DIRECT(t, NV097_SET_BACK_SPECULAR_PARAMS + i * 4,
			NV097_TARGET_LTC1, row * 4 + (i % 4), NV2A_DIRTY_NONE);
	}

	// Per-light local range: 1 component × 8 lights → ltc1 rows r0 + light, component 0
	for (int s = 0; s < 8; s++)
		NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_LOCAL_RANGE + s * 0x80,
			NV097_TARGET_LTC1, (NV_IGRAPH_XF_LTC1_r0 + s) * 4 + 0, NV2A_DIRTY_NONE);

	// ===== DIRECT WRITES: light geometry (non-SRAM, stored in PGRAPHState) =====
	// Per-light: 4 properties × 3 components = 12 floats.
	// Flat index = light * 12 + property_offset + component.
	for (int s = 0; s < 8; s++) {
		for (int i = 0; i < 3; i++) {
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_INFINITE_HALF_VECTOR + i * 4 + s * 0x80,
				NV097_TARGET_LIGHT, s * 12 + 0 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_INFINITE_DIRECTION + i * 4 + s * 0x80,
				NV097_TARGET_LIGHT, s * 12 + 3 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_LOCAL_POSITION + i * 4 + s * 0x80,
				NV097_TARGET_LIGHT, s * 12 + 6 + i, NV2A_DIRTY_NONE);
			NV097_REG_SRAM_DIRECT(t, NV097_SET_LIGHT_LOCAL_ATTENUATION + i * 4 + s * 0x80,
				NV097_TARGET_LIGHT, s * 12 + 9 + i, NV2A_DIRTY_NONE);
		}
	}

	// ===== INDIRECT SRAM WRITE: transform constants (32 methods → xfctx via CONST_LD_PTR) =====
	// reg_index = component within float4 (0..3); actual SRAM row resolved at
	// dispatch time from NV_PGRAPH_CHEOPS_OFFSET.CONST_LD_PTR, auto-incrementing.
	for (int i = 0; i < 32; i++) {
		auto &e = t[(NV097_SET_TRANSFORM_CONSTANT + i * 4) >> 2];
		e.reg_index = i % 4;
		e.width = 0;
		e.shift = 0;
		e.identity = 0;
		e.target = NV097_TARGET_XFCTX_INDIRECT;
		e.dirty_group = NV2A_DIRTY_NONE;
	}

	// ===== INDIRECT SRAM WRITE: transform program (32 methods → xfpr via PROG_LD_PTR) =====
	for (int i = 0; i < 32; i++) {
		auto &e = t[(NV097_SET_TRANSFORM_PROGRAM + i * 4) >> 2];
		e.reg_index = i % 4;
		e.width = 0;
		e.shift = 0;
		e.identity = 0;
		e.target = NV097_TARGET_XFPR_INDIRECT;
		e.dirty_group = NV2A_DIRTY_PROGRAM;
	}
}

// Attempt data-driven dispatch of an NV097 method.
// Returns the old regs[] value of the written register (before overwrite),
// or 0 if the method has no table entry.
static inline uint32_t nv097_dispatch_method(PGRAPHState *pg, unsigned int method, uint32_t parameter)
{
	// ---- Table-driven register dispatch ----
	unsigned int idx = method >> 2;
	if (idx >= NV097_METHOD_TABLE_SIZE) return 0;

	const NV097MethodEntry &entry = nv097_method_table[idx];
	if (entry.target == 0 && entry.reg_index == 0) return 0; // all-zero = no table entry

	uint8_t target = entry.target;

	// ---- Indirect SRAM writes (auto-incrementing CHEOPS_OFFSET pointer) ----
	// For XFCTX_INDIRECT and XFPR_INDIRECT, reg_index is the float4 component
	// (0..3) and the SRAM row is read from CHEOPS_OFFSET at dispatch time.
	// The load pointer auto-increments after every 4th component.
	if (target >= NV097_TARGET_XFCTX_INDIRECT) {
		uint32_t cheops = pg->regs[RI(NV_PGRAPH_CHEOPS_OFFSET)];
		uint32_t *sram;
		int load_ptr;

		if (target == NV097_TARGET_XFCTX_INDIRECT) {
			sram = (uint32_t*)pg->xf.xfctx;
			load_ptr = GET_MASK(cheops, NV_PGRAPH_CHEOPS_OFFSET_CONST_LD_PTR);
			assert(load_ptr < NV2A_VERTEXSHADER_CONSTANTS);
		} else { // NV097_TARGET_XFPR_INDIRECT
			sram = (uint32_t*)pg->xf.xfpr;
			load_ptr = GET_MASK(cheops, NV_PGRAPH_CHEOPS_OFFSET_PROG_LD_PTR);
			assert(load_ptr < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
		}

		int component = entry.reg_index; // 0..3
		int flat_index = load_ptr * 4 + component;

		uint32_t &slot = sram[flat_index];
		uint32_t old_val = slot;
		slot = parameter;

		if (slot != old_val) {
			// Per-constant dirty tracking for xfctx; global dirty flag for xfpr
			if (target == NV097_TARGET_XFCTX_INDIRECT) {
				pg->xf.xfctx_dirty[load_ptr / 32] |= (1u << (load_ptr % 32));
				pg->xf.xfctx_generation++;
			}
			if (entry.dirty_group)
				pg->dirty[entry.dirty_group]++;
		}

		// Auto-increment load pointer after 4th component
		if (component == 3) {
			if (target == NV097_TARGET_XFCTX_INDIRECT)
				SET_MASK(pg->regs[RI(NV_PGRAPH_CHEOPS_OFFSET)],
						 NV_PGRAPH_CHEOPS_OFFSET_CONST_LD_PTR, load_ptr + 1);
			else
				SET_MASK(pg->regs[RI(NV_PGRAPH_CHEOPS_OFFSET)],
						 NV_PGRAPH_CHEOPS_OFFSET_PROG_LD_PTR, load_ptr + 1);
		}

		return old_val;
	}

	// ---- Direct dispatch: resolve target array base pointer ----
	uint32_t *base;
	auto dirty_group = entry.dirty_group;
	uint32_t *row_dirty = nullptr;
	switch (target) {
	default: // NV097_TARGET_PGRAPH (0) — most common path
		base = pg->regs;
		break;
	case NV097_TARGET_XFCTX:
		base = (uint32_t*)pg->xf.xfctx;
		row_dirty = pg->xf.xfctx_dirty;
		break;
	case NV097_TARGET_LTCTXA:
		base = (uint32_t*)pg->xf.ltctxa;
		dirty_group = NV2A_DIRTY_LIGHTING;
		row_dirty = pg->xf.ltctxa_dirty;
		break;
	case NV097_TARGET_LTCTXB:
		base = (uint32_t*)pg->xf.ltctxb;
		dirty_group = NV2A_DIRTY_LIGHTING;
		row_dirty = pg->xf.ltctxb_dirty;
		break;
	case NV097_TARGET_LTC1:
		base = (uint32_t*)pg->xf.ltc1;
		dirty_group = NV2A_DIRTY_LIGHTING;
		row_dirty = pg->xf.ltc1_dirty;
		break;
	case NV097_TARGET_LIGHT:
		base = (uint32_t*)pg->light;
		dirty_group = NV2A_DIRTY_LIGHTING;
		break;
	}

	uint32_t &reg = base[entry.reg_index];
	uint32_t old_val = reg;

	if (entry.width == 0) {
		// Full 32-bit copy
		reg = parameter;
	} else {
		uint32_t mask = ((1u << entry.width) - 1) << entry.shift;
		if (entry.identity) {
			// Identity copy: parameter bits are already at the mask position
			reg = (reg & ~mask) | (parameter & mask);
		} else {
			// Value-shift: parameter is a small value, shift it into position
			SET_MASK(reg, mask, parameter);
		}
	}

	if (reg != old_val) {
		if (target == NV097_TARGET_PGRAPH)
			pg->dirty[NV2A_DIRTY_PGRAPH]++;
		if (row_dirty) {
			unsigned row_idx = entry.reg_index / 4;
			row_dirty[row_idx / 32] |= (1u << (row_idx % 32));
			if (target == NV097_TARGET_XFCTX)
				pg->xf.xfctx_generation++;
		}
		if (dirty_group)
			pg->dirty[dirty_group]++;
	}

	return old_val;
}
