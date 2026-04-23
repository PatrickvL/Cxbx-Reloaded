// ******************************************************************
// *  NV2A Method -> PGRAPH Register Dispatch Table
// *
// *  Data-driven NV097 method -> NV_PGRAPH register mapping.
// *  Covers DIRECT_COPY, DIRECT_COPY_RANGED, and MASKED_WRITE patterns.
// *  Methods requiring custom logic fall through (pgraph_reg == 0).
// *
// *  Table indexed by method >> 2 (dword index into method space).
// *  Each entry: { pgraph_reg, shift, width } — 4 bytes.
// *    pgraph_reg == 0  -> not data-driven, use switch fallback
// *    width == 0       -> full 32-bit copy: regs[RI(reg)] = parameter
// *    width > 0        -> masked write using mask = ((1<<width)-1) << shift
// *      shift bit 7 clear -> value-shift: SET_MASK(reg, mask, parameter)
// *      shift bit 7 set   -> identity:    reg = (reg & ~mask) | (param & mask)
// *
// *  For a future PFIFO-in-shader: this table fits in an 8KB structured
// *  buffer (2048 entries x 4 bytes) and the dispatch is branchless for
// *  all data-driven methods.  GPU branchless form:
// *    uint s = shift & 0x7F;
// *    uint mask = ((1u << width) - 1) << s;
// *    uint val = (shift & 0x80) ? (param & mask) : ((param << s) & mask);
// *    reg = (reg & ~mask) | val;
// ******************************************************************
#pragma once

#include "nv2a_int.h"
#include "nv2a_regs.h"

#include <cstring> // memset

struct NV097MethodEntry {
	uint16_t pgraph_reg;  // NV_PGRAPH_* byte offset, 0 = custom handler
	uint8_t shift;        // bits 0-6: lowest bit of mask field
	                      // bit 7: identity-copy flag (NV097_MASK_IDENTITY)
	uint8_t width;        // 0 = full 32-bit copy, >0 = masked write (width bits)
};

// Shift bit 7: parameter bits are already at mask position (identity copy)
#define NV097_MASK_IDENTITY 0x80

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
#define NV097_REG_DIRECT(tbl, nv097_method, pg_reg) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.pgraph_reg = (uint16_t)(pg_reg); \
		e.shift = 0; \
		e.width = 0; \
	} while (0)

// Helper to register a direct-copy ranged method (variable stride, N slots)
#define NV097_REG_DIRECT_RANGE(tbl, nv097_base, stride, count, pg_base) \
	do { \
		for (int _i = 0; _i < (count); _i++) { \
			auto &e = (tbl)[((nv097_base) + _i * (stride)) >> 2]; \
			e.pgraph_reg = (uint16_t)((pg_base) + _i * 4); \
			e.shift = 0; \
			e.width = 0; \
		} \
	} while (0)

// Helper to register a value-shift masked write: SET_MASK(reg, mask, parameter)
// Parameter is a small value at bit 0; SET_MASK shifts it into the mask position.
#define NV097_REG_MASKED(tbl, nv097_method, pg_reg, pg_mask) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.pgraph_reg = (uint16_t)(pg_reg); \
		e.shift = (uint8_t)(ffs(pg_mask) - 1); \
		e.width = nv097_mask_width(pg_mask); \
	} while (0)

// Helper to register an identity-copy masked write:
// reg = (reg & ~mask) | (parameter & mask)
// Parameter bits are already at the mask position (no shift needed).
#define NV097_REG_MASKED_IDENTITY(tbl, nv097_method, pg_reg, pg_mask) \
	do { \
		auto &e = (tbl)[(nv097_method) >> 2]; \
		e.pgraph_reg = (uint16_t)(pg_reg); \
		e.shift = (uint8_t)((ffs(pg_mask) - 1) | NV097_MASK_IDENTITY); \
		e.width = nv097_mask_width(pg_mask); \
	} while (0)

// Build the dispatch table at startup (called once)
static NV097MethodEntry nv097_method_table[NV097_METHOD_TABLE_SIZE];

static void nv097_init_method_table()
{
	// Zero-init: all pgraph_reg == 0 means custom (fallback to switch)
	memset(nv097_method_table, 0, sizeof(nv097_method_table));

	auto &t = nv097_method_table; // shorthand

	// ===== DIRECT_COPY: pg->regs[RI(reg)] = parameter =====
	NV097_REG_DIRECT(t, NV097_SET_COMBINER_SPECULAR_FOG_CW0, NV_PGRAPH_COMBINESPECFOG0);
	NV097_REG_DIRECT(t, NV097_SET_COMBINER_SPECULAR_FOG_CW1, NV_PGRAPH_COMBINESPECFOG1);
	NV097_REG_DIRECT(t, NV097_SET_BLEND_COLOR,                NV_PGRAPH_BLENDCOLOR);
	NV097_REG_DIRECT(t, NV097_SET_POLYGON_OFFSET_SCALE_FACTOR,NV_PGRAPH_ZOFFSETFACTOR);
	NV097_REG_DIRECT(t, NV097_SET_POLYGON_OFFSET_BIAS,        NV_PGRAPH_ZOFFSETBIAS);
	NV097_REG_DIRECT(t, NV097_SET_CLIP_MIN,                   NV_PGRAPH_ZCLIPMIN);
	NV097_REG_DIRECT(t, NV097_SET_CLIP_MAX,                   NV_PGRAPH_ZCLIPMAX);
	NV097_REG_DIRECT(t, NV097_SET_SEMAPHORE_OFFSET,           NV_PGRAPH_SEMAPHOREOFFSET);
	NV097_REG_DIRECT(t, NV097_SET_ZSTENCIL_CLEAR_VALUE,       NV_PGRAPH_ZSTENCILCLEARVALUE);
	NV097_REG_DIRECT(t, NV097_SET_COLOR_CLEAR_VALUE,          NV_PGRAPH_COLORCLEARVALUE);
	NV097_REG_DIRECT(t, NV097_SET_CLEAR_RECT_HORIZONTAL,      NV_PGRAPH_CLEARRECTX);
	NV097_REG_DIRECT(t, NV097_SET_CLEAR_RECT_VERTICAL,        NV_PGRAPH_CLEARRECTY);
	NV097_REG_DIRECT(t, NV097_SET_SHADER_CLIP_PLANE_MODE,     NV_PGRAPH_SHADERCLIPMODE);
	NV097_REG_DIRECT(t, NV097_SET_COMBINER_CONTROL,           NV_PGRAPH_COMBINECTL);
	NV097_REG_DIRECT(t, NV097_SET_SHADER_STAGE_PROGRAM,       NV_PGRAPH_SHADERPROG);

	// ===== DIRECT_COPY_RANGED (stride 4) =====
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_ALPHA_ICW,  4, 8, NV_PGRAPH_COMBINEALPHAI0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_WINDOW_CLIP_HORIZONTAL, 4, 8, NV_PGRAPH_WINDOWCLIPX0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_WINDOW_CLIP_VERTICAL,   4, 8, NV_PGRAPH_WINDOWCLIPY0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_FACTOR0,    4, 8, NV_PGRAPH_COMBINEFACTOR0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_FACTOR1,    4, 8, NV_PGRAPH_COMBINEFACTOR1);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_ALPHA_OCW,  4, 8, NV_PGRAPH_COMBINEALPHAO0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_COLOR_ICW,  4, 8, NV_PGRAPH_COMBINECOLORI0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_COMBINER_COLOR_OCW,  4, 8, NV_PGRAPH_COMBINECOLORO0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_SPECULAR_FOG_FACTOR, 4, 2, NV_PGRAPH_SPECFOGFACTOR0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_EYE_VECTOR,          4, 3, NV_PGRAPH_EYEVEC0);

	// ===== DIRECT_COPY_RANGED (stride 64, texture methods) =====
	// NOTE: TEXTURE_IMAGE_RECT omitted — switch also sets pg->texture_dirty[slot]
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_ADDRESS,      64, 4, NV_PGRAPH_TEXADDRESS0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_CONTROL0,     64, 4, NV_PGRAPH_TEXCTL0_0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_CONTROL1,     64, 4, NV_PGRAPH_TEXCTL1_0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_FILTER,       64, 4, NV_PGRAPH_TEXFILTER0);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_BORDER_COLOR, 64, 4, NV_PGRAPH_BORDERCOLOR0);

	// ===== DIRECT_COPY_RANGED (stride 64, bump env — stages 1-3 only) =====
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_SET_BUMP_ENV_SCALE + 64,  64, 3, NV_PGRAPH_BUMPSCALE1);
	NV097_REG_DIRECT_RANGE(t, NV097_SET_TEXTURE_SET_BUMP_ENV_OFFSET + 64, 64, 3, NV_PGRAPH_BUMPOFFSET1);

	// ===== MASKED_WRITE (value-shift): SET_MASK(regs[RI(reg)], mask, parameter) =====
	NV097_REG_MASKED(t, NV097_SET_FLIP_READ,             NV_PGRAPH_SURFACE,     NV_PGRAPH_SURFACE_READ_3D);
	NV097_REG_MASKED(t, NV097_SET_FLIP_WRITE,            NV_PGRAPH_SURFACE,     NV_PGRAPH_SURFACE_WRITE_3D);
	NV097_REG_MASKED(t, NV097_SET_FLIP_MODULO,           NV_PGRAPH_SURFACE,     NV_PGRAPH_SURFACE_MODULO_3D);
	NV097_REG_MASKED(t, NV097_SET_FOG_ENABLE,            NV_PGRAPH_CONTROL_3,   NV_PGRAPH_CONTROL_3_FOGENABLE);
	NV097_REG_MASKED(t, NV097_SET_WINDOW_CLIP_TYPE,      NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_WINDOWCLIPTYPE);
	NV097_REG_MASKED(t, NV097_SET_ALPHA_TEST_ENABLE,     NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ALPHATESTENABLE);
	NV097_REG_MASKED(t, NV097_SET_BLEND_ENABLE,          NV_PGRAPH_BLEND,       NV_PGRAPH_BLEND_EN);
	NV097_REG_MASKED(t, NV097_SET_CULL_FACE_ENABLE,      NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_CULLENABLE);
	NV097_REG_MASKED(t, NV097_SET_DEPTH_TEST_ENABLE,     NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ZENABLE);
	NV097_REG_MASKED(t, NV097_SET_DITHER_ENABLE,         NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_DITHERENABLE);
	NV097_REG_MASKED(t, NV097_SET_LIGHTING_ENABLE,       NV_PGRAPH_CSV0_C,      NV_PGRAPH_CSV0_C_LIGHTING);
	NV097_REG_MASKED(t, NV097_SET_SKIN_MODE,             NV_PGRAPH_CSV0_D,      NV_PGRAPH_CSV0_D_SKIN);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_TEST_ENABLE,   NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE);
	NV097_REG_MASKED(t, NV097_SET_POLY_OFFSET_POINT_ENABLE, NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POFFSETPOINTENABLE);
	NV097_REG_MASKED(t, NV097_SET_POLY_OFFSET_LINE_ENABLE,  NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POFFSETLINEENABLE);
	NV097_REG_MASKED(t, NV097_SET_POLY_OFFSET_FILL_ENABLE,  NV_PGRAPH_SETUPRASTER, NV_PGRAPH_SETUPRASTER_POFFSETFILLENABLE);
	NV097_REG_MASKED(t, NV097_SET_ALPHA_FUNC,            NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ALPHAFUNC);
	NV097_REG_MASKED(t, NV097_SET_ALPHA_REF,             NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ALPHAREF);
	NV097_REG_MASKED(t, NV097_SET_DEPTH_FUNC,            NV_PGRAPH_CONTROL_0,   NV_PGRAPH_CONTROL_0_ZFUNC);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_MASK,          NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_MASK_WRITE);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_FUNC,          NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_FUNC);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_FUNC_REF,      NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_REF);
	NV097_REG_MASKED(t, NV097_SET_STENCIL_FUNC_MASK,     NV_PGRAPH_CONTROL_1,   NV_PGRAPH_CONTROL_1_STENCIL_MASK_READ);
	NV097_REG_MASKED(t, NV097_SET_NORMALIZATION_ENABLE,   NV_PGRAPH_CSV0_C,      NV_PGRAPH_CSV0_C_NORMALIZATION_ENABLE);
	NV097_REG_MASKED(t, NV097_SET_LIGHT_ENABLE_MASK,      NV_PGRAPH_CSV0_D,      NV_PGRAPH_CSV0_D_LIGHTS);
	NV097_REG_MASKED(t, NV097_SET_TEXGEN_VIEW_MODEL,      NV_PGRAPH_CSV0_D,      NV_PGRAPH_CSV0_D_TEXGEN_REF);
	NV097_REG_MASKED(t, NV097_SET_LOGIC_OP_ENABLE,        NV_PGRAPH_BLEND,       NV_PGRAPH_BLEND_LOGICOP_ENABLE);
	NV097_REG_MASKED(t, NV097_SET_LOGIC_OP,               NV_PGRAPH_BLEND,       NV_PGRAPH_BLEND_LOGICOP);

	// ===== MASKED_WRITE (value-shift, mask at bit 0 — equivalent to identity) =====
	NV097_REG_MASKED(t, NV097_SET_DOT_RGBMAPPING,         NV_PGRAPH_SHADERCTL, NV_PGRAPH_SHADERCTL_DOT_RGBMAPPING);

	// ===== MASKED_WRITE (identity-copy: param bits already at mask position) =====
	NV097_REG_MASKED_IDENTITY(t, NV097_SET_SHADER_OTHER_STAGE_INPUT, NV_PGRAPH_SHADERCTL, NV_PGRAPH_SHADERCTL_OTHER_STAGE_INPUT);
}

// Attempt data-driven dispatch of an NV097 method.
// Returns true if handled, false if caller should use the switch fallback.
static inline bool nv097_dispatch_method(PGRAPHState *pg, unsigned int method, uint32_t parameter)
{
	unsigned int idx = method >> 2;
	if (idx >= NV097_METHOD_TABLE_SIZE) return false;

	const NV097MethodEntry &entry = nv097_method_table[idx];
	if (entry.pgraph_reg == 0) return false;

	if (entry.width == 0) {
		// Full 32-bit copy
		pg->regs[RI(entry.pgraph_reg)] = parameter;
	} else {
		uint8_t s = entry.shift & 0x7F;
		uint32_t mask = ((1u << entry.width) - 1) << s;
		if (entry.shift & NV097_MASK_IDENTITY) {
			// Identity copy: parameter bits are already at the mask position
			pg->regs[RI(entry.pgraph_reg)] = (pg->regs[RI(entry.pgraph_reg)] & ~mask) | (parameter & mask);
		} else {
			// Value-shift: parameter is a small value, shift it into position
			SET_MASK(pg->regs[RI(entry.pgraph_reg)], mask, parameter);
		}
	}
	return true;
}
