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
// *  All rights reserved
// *
// ******************************************************************

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h"
#include "devices\video\nv2a.h"        // PGRAPHState, nv2a_regs.h, GET_MASK, RI
#include "core\hle\D3D8\Rendering\NV2A_PGRAPH_Helpers.h"
#include "common/util/hasher.h"
#include <algorithm>                    // std::min, std::sort
#include <unordered_map>
#include <vector>

// Tracked viewport dimensions — updated every time the viewport is set.
// Used by GS constant buffer update to avoid RSGetViewports() per draw.
float g_CurrentViewportWidth = 0.0f;
float g_CurrentViewportHeight = 0.0f;

// ******************************************************************
// * State object cache — avoids redundant Create*State calls.
// * Xbox games use a small number of unique state combinations;
// * the cache typically holds < 50 entries for each type.
// ******************************************************************
namespace {
	struct DescHash {
		template <typename T>
		size_t operator()(const T& desc) const {
			return static_cast<size_t>(ComputeHash(&desc, sizeof(T)));
		}
	};
	struct DescEqual {
		template <typename T>
		bool operator()(const T& a, const T& b) const {
			return std::memcmp(&a, &b, sizeof(T)) == 0;
		}
	};
	std::unordered_map<D3D11_RASTERIZER_DESC, Microsoft::WRL::ComPtr<ID3D11RasterizerState>, DescHash, DescEqual> s_RasterizerCache;
	std::unordered_map<D3D11_DEPTH_STENCIL_DESC, Microsoft::WRL::ComPtr<ID3D11DepthStencilState>, DescHash, DescEqual> s_DepthStencilCache;
	std::unordered_map<D3D11_BLEND_DESC, Microsoft::WRL::ComPtr<ID3D11BlendState>, DescHash, DescEqual> s_BlendCache;
}

// ******************************************************************
// * Unified D3D11 render state mapping
// * Called from ApplySimpleRenderState and ApplyComplexRenderState
// * after Xbox→PC value conversion. Updates D3D11 state descriptors
// * and sets dirty flags for deferred state object recreation.
// ******************************************************************

// Helper: remap color-referencing blend factors to their alpha equivalents
// for use in BlendAlpha slots (D3D11 separates color/alpha factor interpretation)
static D3D11_BLEND RemapBlendForAlpha(D3D11_BLEND blend)
{
	switch (blend) {
	case D3D11_BLEND_SRC_COLOR:      return D3D11_BLEND_SRC_ALPHA;
	case D3D11_BLEND_INV_SRC_COLOR:  return D3D11_BLEND_INV_SRC_ALPHA;
	case D3D11_BLEND_DEST_COLOR:     return D3D11_BLEND_DEST_ALPHA;
	case D3D11_BLEND_INV_DEST_COLOR: return D3D11_BLEND_INV_DEST_ALPHA;
	default:                         return blend;
	}
}

// ******************************************************************
// * Map NV2A PGRAPH blend factor (4-bit) → D3D11_BLEND
// * PGRAPH values 0-10 map to D3D11 values 1-11 (factor + 1).
// * PGRAPH values 12-15 (constant color/alpha) map to D3D11
// * BLEND_FACTOR / INV_BLEND_FACTOR.
// ******************************************************************
static D3D11_BLEND MapPGRAPHBlendFactor(unsigned int factor)
{
	switch (factor) {
	case NV_PGRAPH_BLEND_SFACTOR_ZERO:                    return D3D11_BLEND_ZERO;
	case NV_PGRAPH_BLEND_SFACTOR_ONE:                     return D3D11_BLEND_ONE;
	case NV_PGRAPH_BLEND_SFACTOR_SRC_COLOR:               return D3D11_BLEND_SRC_COLOR;
	case NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_SRC_COLOR:     return D3D11_BLEND_INV_SRC_COLOR;
	case NV_PGRAPH_BLEND_SFACTOR_SRC_ALPHA:               return D3D11_BLEND_SRC_ALPHA;
	case NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_SRC_ALPHA:     return D3D11_BLEND_INV_SRC_ALPHA;
	case NV_PGRAPH_BLEND_SFACTOR_DST_ALPHA:               return D3D11_BLEND_DEST_ALPHA;
	case NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_DST_ALPHA:     return D3D11_BLEND_INV_DEST_ALPHA;
	case NV_PGRAPH_BLEND_SFACTOR_DST_COLOR:               return D3D11_BLEND_DEST_COLOR;
	case NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_DST_COLOR:     return D3D11_BLEND_INV_DEST_COLOR;
	case NV_PGRAPH_BLEND_SFACTOR_SRC_ALPHA_SATURATE:      return D3D11_BLEND_SRC_ALPHA_SAT;
	case NV_PGRAPH_BLEND_SFACTOR_CONSTANT_COLOR:          return D3D11_BLEND_BLEND_FACTOR;
	case NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_CONSTANT_COLOR:return D3D11_BLEND_INV_BLEND_FACTOR;
	case NV_PGRAPH_BLEND_SFACTOR_CONSTANT_ALPHA:          return D3D11_BLEND_BLEND_FACTOR;
	case NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_CONSTANT_ALPHA:return D3D11_BLEND_INV_BLEND_FACTOR;
	default:                                              return D3D11_BLEND_ONE;
	}
}

// ******************************************************************
// * Map NV2A PGRAPH blend equation (3-bit) → D3D11_BLEND_OP
// * PGRAPH encoding (set by pgraph_handle_method NV097_SET_BLEND_EQUATION):
// *   0=SUBTRACT, 1=REV_SUBTRACT, 2=ADD, 3=MIN, 4=MAX,
// *   5=REV_SUBTRACT_SIGNED, 6=ADD_SIGNED
// ******************************************************************
static D3D11_BLEND_OP MapPGRAPHBlendOp(unsigned int eqn)
{
	switch (eqn) {
	case 0:  return D3D11_BLEND_OP_SUBTRACT;
	case 1:  return D3D11_BLEND_OP_REV_SUBTRACT;
	case 2:  return D3D11_BLEND_OP_ADD;
	case 3:  return D3D11_BLEND_OP_MIN;
	case 4:  return D3D11_BLEND_OP_MAX;
	case 5:  return D3D11_BLEND_OP_REV_SUBTRACT; // signed — approximate
	case 6:  return D3D11_BLEND_OP_ADD;           // signed — approximate
	default: return D3D11_BLEND_OP_ADD;
	}
}

// ******************************************************************
// * Read NV2A PGRAPH registers and populate D3D11 state descriptors.
// * This is the PGRAPH-driven replacement for XboxRenderStates.Apply()
// * for blend, depth-stencil, and rasterizer pipeline state.
// *
// * Only sets dirty flags when the relevant PGRAPH registers actually
// * changed since the last call, avoiding expensive state object
// * recreation (CreateBlendState etc.) on every draw.
// ******************************************************************

// Cached PGRAPH register values for change-detection
static uint32_t s_CachedBlendReg = ~0u;
static uint32_t s_CachedBlendColorReg = ~0u;
static uint32_t s_CachedControl0Reg = ~0u;
static uint32_t s_CachedControl1Reg = ~0u;
static uint32_t s_CachedControl2Reg = ~0u;

// Geometry shader binding cache
static ID3D11GeometryShader* s_LastBoundGS = (ID3D11GeometryShader*)~0ull;
static uint32_t s_CachedControl3Reg = ~0u;
static uint32_t s_CachedSetupRasterReg = ~0u;
static uint32_t s_CachedZOffsetBiasReg = ~0u;
static uint32_t s_CachedZOffsetFactorReg = ~0u;
void CxbxD3D11UpdatePipelineStateFromPGRAPH(PGRAPHState *pg)
{
	if (!pg) return;

	// Read CONTROL_0 once — shared between blend (write mask) and depth/stencil blocks
	uint32_t ctrl0 = pg->regs[RI(NV_PGRAPH_CONTROL_0)];

	// ---- Blend state from NV_PGRAPH_BLEND (0x1804) ----
	{
		uint32_t blend = pg->regs[RI(NV_PGRAPH_BLEND)];
		uint32_t bc = pg->regs[RI(NV_PGRAPH_BLENDCOLOR)];

		// Only rebuild blend state if relevant registers changed
		if (blend != s_CachedBlendReg || bc != s_CachedBlendColorReg ||
			(ctrl0 & 0x3C000000u) != (s_CachedControl0Reg & 0x3C000000u)) { // color write mask bits 26-29
			s_CachedBlendReg = blend;
			s_CachedBlendColorReg = bc;

			g_D3D11BlendDesc.RenderTarget[0].BlendEnable = (blend & NV_PGRAPH_BLEND_EN) ? TRUE : FALSE;

			unsigned int sfactor = GET_MASK(blend, NV_PGRAPH_BLEND_SFACTOR);
			unsigned int dfactor = GET_MASK(blend, NV_PGRAPH_BLEND_DFACTOR);
			unsigned int eqn    = GET_MASK(blend, NV_PGRAPH_BLEND_EQN);

			D3D11_BLEND srcBlend  = MapPGRAPHBlendFactor(sfactor);
			D3D11_BLEND destBlend = MapPGRAPHBlendFactor(dfactor);
			D3D11_BLEND_OP blendOp = MapPGRAPHBlendOp(eqn);

			g_D3D11BlendDesc.RenderTarget[0].SrcBlend      = srcBlend;
			g_D3D11BlendDesc.RenderTarget[0].SrcBlendAlpha  = RemapBlendForAlpha(srcBlend);
			g_D3D11BlendDesc.RenderTarget[0].DestBlend     = destBlend;
			g_D3D11BlendDesc.RenderTarget[0].DestBlendAlpha = RemapBlendForAlpha(destBlend);
			g_D3D11BlendDesc.RenderTarget[0].BlendOp       = blendOp;
			g_D3D11BlendDesc.RenderTarget[0].BlendOpAlpha  = blendOp;

			// NV2A BLENDCOLOR is ARGB packed
			float blendR = ((bc >> 16) & 0xFF) / 255.0f;
			float blendG = ((bc >> 8) & 0xFF) / 255.0f;
			float blendB = (bc & 0xFF) / 255.0f;
			float blendA = ((bc >> 24) & 0xFF) / 255.0f;

			// NV2A CONSTANT_ALPHA uses only the alpha for ALL channels.
			// D3D11 BLEND_FACTOR applies per-channel, so when either factor
			// references CONSTANT_ALPHA we must splat alpha into RGB.
			bool usesConstAlpha = (sfactor == NV_PGRAPH_BLEND_SFACTOR_CONSTANT_ALPHA ||
			                       sfactor == NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_CONSTANT_ALPHA ||
			                       dfactor == NV_PGRAPH_BLEND_SFACTOR_CONSTANT_ALPHA ||
			                       dfactor == NV_PGRAPH_BLEND_SFACTOR_ONE_MINUS_CONSTANT_ALPHA);
			if (usesConstAlpha) {
				g_D3D11BlendFactor[0] = blendA;
				g_D3D11BlendFactor[1] = blendA;
				g_D3D11BlendFactor[2] = blendA;
				g_D3D11BlendFactor[3] = blendA;
			} else {
				g_D3D11BlendFactor[0] = blendR;
				g_D3D11BlendFactor[1] = blendG;
				g_D3D11BlendFactor[2] = blendB;
				g_D3D11BlendFactor[3] = blendA;
			}

			// Color write mask
			UINT8 writeMask = 0;
			if (ctrl0 & NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE)   writeMask |= D3D11_COLOR_WRITE_ENABLE_RED;
			if (ctrl0 & NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE) writeMask |= D3D11_COLOR_WRITE_ENABLE_GREEN;
			if (ctrl0 & NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE)  writeMask |= D3D11_COLOR_WRITE_ENABLE_BLUE;
			if (ctrl0 & NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE) writeMask |= D3D11_COLOR_WRITE_ENABLE_ALPHA;
			g_D3D11BlendDesc.RenderTarget[0].RenderTargetWriteMask = writeMask;

			g_bD3D11BlendStateDirty = true;
		}
	}

	// ---- Depth/stencil state from NV_PGRAPH_CONTROL_0/1/2 ----
	{
		uint32_t ctrl1 = pg->regs[RI(NV_PGRAPH_CONTROL_1)];
		uint32_t ctrl2 = pg->regs[RI(NV_PGRAPH_CONTROL_2)];

		// Only rebuild if depth/stencil registers changed
		// ctrl0 depth bits are below bit 26 (write mask is bits 26-29, handled above)
		if (ctrl0 != s_CachedControl0Reg || ctrl1 != s_CachedControl1Reg || ctrl2 != s_CachedControl2Reg) {
			s_CachedControl0Reg = ctrl0;
			s_CachedControl1Reg = ctrl1;
			s_CachedControl2Reg = ctrl2;

			g_D3D11DepthStencilDesc.DepthEnable = (ctrl0 & NV_PGRAPH_CONTROL_0_ZENABLE) ? TRUE : FALSE;
			g_D3D11DepthStencilDesc.DepthWriteMask = (ctrl0 & NV_PGRAPH_CONTROL_0_ZWRITEENABLE)
				? D3D11_DEPTH_WRITE_MASK_ALL : D3D11_DEPTH_WRITE_MASK_ZERO;

			unsigned int zfunc = GET_MASK(ctrl0, NV_PGRAPH_CONTROL_0_ZFUNC);
			g_D3D11DepthStencilDesc.DepthFunc = (D3D11_COMPARISON_FUNC)(zfunc + 1);

			g_D3D11DepthStencilDesc.StencilEnable = (ctrl1 & NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE) ? TRUE : FALSE;

			unsigned int sfunc = GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_FUNC);
			D3D11_COMPARISON_FUNC stencilFunc = (D3D11_COMPARISON_FUNC)(sfunc + 1);
			g_D3D11DepthStencilDesc.FrontFace.StencilFunc = stencilFunc;
			g_D3D11DepthStencilDesc.BackFace.StencilFunc  = stencilFunc;

			g_D3D11StencilRef = GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_REF);
			g_D3D11DepthStencilDesc.StencilReadMask  = (UINT8)GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_MASK_READ);
			g_D3D11DepthStencilDesc.StencilWriteMask = (UINT8)GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_MASK_WRITE);

			D3D11_STENCIL_OP failOp  = (D3D11_STENCIL_OP)GET_MASK(ctrl2, NV_PGRAPH_CONTROL_2_STENCIL_OP_FAIL);
			D3D11_STENCIL_OP zfailOp = (D3D11_STENCIL_OP)GET_MASK(ctrl2, NV_PGRAPH_CONTROL_2_STENCIL_OP_ZFAIL);
			D3D11_STENCIL_OP zpassOp = (D3D11_STENCIL_OP)GET_MASK(ctrl2, NV_PGRAPH_CONTROL_2_STENCIL_OP_ZPASS);

			g_D3D11DepthStencilDesc.FrontFace.StencilFailOp      = failOp;
			g_D3D11DepthStencilDesc.FrontFace.StencilDepthFailOp = zfailOp;
			g_D3D11DepthStencilDesc.FrontFace.StencilPassOp      = zpassOp;
			g_D3D11DepthStencilDesc.BackFace.StencilFailOp       = failOp;
			g_D3D11DepthStencilDesc.BackFace.StencilDepthFailOp  = zfailOp;
			g_D3D11DepthStencilDesc.BackFace.StencilPassOp       = zpassOp;

			g_bD3D11DepthStencilStateDirty = true;
		}
	}

	// DOT_ZW depth override: On NV2A, DOT_ZW replaces the fragment Z value which
	// then flows through the normal depth pipeline. If the game disabled depth
	// (DepthEnable=FALSE) but DOT_ZW is active, we must force depth writing so
	// SV_Depth from the pixel shader actually reaches the depth buffer.
	{
		uint32_t shaderProg = pg->regs[RI(NV_PGRAPH_SHADERPROG)];
		bool hasDotZW = false;
		for (int i = 0; i < 4; i++) {
			if (((shaderProg >> (i * 5)) & 0x1F) == 0x0A) { hasDotZW = true; break; }
		}
		if (hasDotZW) {
			// Force depth enabled with full write — DOT_ZW replaces
			// the fragment Z unconditionally on NV2A hardware.
			g_D3D11DepthStencilDesc.DepthEnable = TRUE;
			g_D3D11DepthStencilDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
			g_bD3D11DepthStencilStateDirty = true;
		}
	}

	// ---- Rasterizer state from NV_PGRAPH_SETUPRASTER (0x1990) ----
	{
		uint32_t setup = pg->regs[RI(NV_PGRAPH_SETUPRASTER)];
		uint32_t zBiasReg = pg->regs[RI(NV_PGRAPH_ZOFFSETBIAS)];
		uint32_t zFactorReg = pg->regs[RI(NV_PGRAPH_ZOFFSETFACTOR)];
		// ---- Point sprite enable from NV_PGRAPH_SETUPRASTER ----
		// D3DRS_POINTSPRITEENABLE → NV097_SET_POINT_SMOOTH_ENABLE →
		// NV_PGRAPH_SETUPRASTER_POINTSMOOTHENABLE. This is distinct from
		// NV_PGRAPH_CONTROL_3_POINTPARAMSENABLE which tracks D3DRS_POINTSCALEENABLE.
		bool bPointSpriteEnabled = (setup & NV_PGRAPH_SETUPRASTER_POINTSMOOTHENABLE) != 0;

		if (s_CachedSetupRasterReg != setup ||
			s_CachedZOffsetBiasReg != zBiasReg ||
			s_CachedZOffsetFactorReg != zFactorReg ||
			g_bPointSpriteEnabled != bPointSpriteEnabled) {
			s_CachedSetupRasterReg = setup;
			s_CachedZOffsetBiasReg = zBiasReg;
			s_CachedZOffsetFactorReg = zFactorReg;
			g_bPointSpriteEnabled = bPointSpriteEnabled;

			// Fill mode: PGRAPH FRONTFACEMODE 0=FILL, 1=POINT, 2=LINE
			unsigned int fillMode = GET_MASK(setup, NV_PGRAPH_SETUPRASTER_FRONTFACEMODE);
			switch (fillMode) {
			case NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_FILL:  g_D3D11RasterizerDesc.FillMode = D3D11_FILL_SOLID; break;
			case NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_LINE:  g_D3D11RasterizerDesc.FillMode = D3D11_FILL_WIREFRAME; break;
			case NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_POINT: g_D3D11RasterizerDesc.FillMode = D3D11_FILL_WIREFRAME; break; // no point fill in D3D11
			default:                                        g_D3D11RasterizerDesc.FillMode = D3D11_FILL_SOLID; break;
			}

			// Front face winding: PGRAPH bit 23 — 0=CW, 1=CCW
			g_D3D11RasterizerDesc.FrontCounterClockwise = (setup & NV_PGRAPH_SETUPRASTER_FRONTFACE) ? TRUE : FALSE;

			// Cull mode
			if (!(setup & NV_PGRAPH_SETUPRASTER_CULLENABLE)) {
				g_D3D11RasterizerDesc.CullMode = D3D11_CULL_NONE;
			} else if (bPointSpriteEnabled) {
				// On real NV2A hardware, point sprites are always front-facing — they are
				// generated by the rasterizer as screen-aligned quads and are never culled.
				// Our GS-expanded quads have a fixed winding order, so we must disable
				// culling when point sprites are active to match hardware behavior.
				g_D3D11RasterizerDesc.CullMode = D3D11_CULL_NONE;
			} else {
				unsigned int cullCtrl = GET_MASK(setup, NV_PGRAPH_SETUPRASTER_CULLCTRL);
				switch (cullCtrl) {
				case NV_PGRAPH_SETUPRASTER_CULLCTRL_FRONT:          g_D3D11RasterizerDesc.CullMode = D3D11_CULL_FRONT; break;
				case NV_PGRAPH_SETUPRASTER_CULLCTRL_BACK:           g_D3D11RasterizerDesc.CullMode = D3D11_CULL_BACK; break;
				case NV_PGRAPH_SETUPRASTER_CULLCTRL_FRONT_AND_BACK: g_D3D11RasterizerDesc.CullMode = D3D11_CULL_NONE; break; // D3D11 can't cull both
				default:                                            g_D3D11RasterizerDesc.CullMode = D3D11_CULL_NONE; break;
				}
			}

			// Line antialiasing
			g_D3D11RasterizerDesc.AntialiasedLineEnable = (setup & NV_PGRAPH_SETUPRASTER_LINESMOOTHENABLE) ? TRUE : FALSE;

			// Depth bias
			float zBias; std::memcpy(&zBias, &zBiasReg, sizeof(float));
			float zFactor; std::memcpy(&zFactor, &zFactorReg, sizeof(float));
			g_D3D11RasterizerDesc.DepthBias = static_cast<INT>(zBias * (float)(1 << 24));
			g_D3D11RasterizerDesc.SlopeScaledDepthBias = zFactor;
			g_D3D11RasterizerDesc.DepthBiasClamp = 0.0f;

			g_bD3D11RasterizerStateDirty = true;
		}
	}
}

// ******************************************************************
// * Map NV2A texture address mode to D3D11 texture address mode.
// * NV2A values: 1=WRAP 2=MIRROR 3=CLAMP_TO_EDGE 4=BORDER 5=CLAMP_OGL
// * D3D11 values: 1=WRAP 2=MIRROR 3=CLAMP 4=BORDER 5=MIRROR_ONCE
// ******************************************************************
static D3D11_TEXTURE_ADDRESS_MODE MapPGRAPHTexAddress(unsigned int addr)
{
	switch (addr) {
	case 1:  return D3D11_TEXTURE_ADDRESS_WRAP;
	case 2:  return D3D11_TEXTURE_ADDRESS_MIRROR;
	case 3:  return D3D11_TEXTURE_ADDRESS_CLAMP;
	case 4:  return D3D11_TEXTURE_ADDRESS_BORDER;
	case 5:  return D3D11_TEXTURE_ADDRESS_CLAMP; // CLAMP_OGL ≈ CLAMP_TO_EDGE
	default: return D3D11_TEXTURE_ADDRESS_WRAP;
	}
}

// ******************************************************************
// * Map NV2A min/mag filter to D3D11 filter components.
// * NV2A min filter: 1=BOX_LOD0(nearest) 2=TENT_LOD0(linear)
// *   3=BOX_NEARESTLOD 4=TENT_NEARESTLOD 5=BOX_TENT_LOD 6=TENT_TENT_LOD
// *   7=CONVOLUTION_2D_LOD0
// * NV2A mag filter: 1=BOX(nearest) 2=TENT(linear) 4=CONVOLUTION_2D
// ******************************************************************
static D3D11_FILTER BuildD3D11Filter(unsigned int minFilter, unsigned int magFilter)
{
	// Decode NV2A filter to (min, mag, mip) triplet
	bool minLinear = false, magLinear = false, mipLinear = false;
	bool anisotropic = false;

	switch (minFilter) {
	case 1: // BOX_LOD0 = nearest, no mip
		minLinear = false; mipLinear = false; break;
	case 2: // TENT_LOD0 = linear, no mip
		minLinear = true; mipLinear = false; break;
	case 3: // BOX_NEARESTLOD = nearest, nearest mip
		minLinear = false; mipLinear = false; break;
	case 4: // TENT_NEARESTLOD = linear, nearest mip
		minLinear = true; mipLinear = false; break;
	case 5: // BOX_TENT_LOD = nearest, linear mip
		minLinear = false; mipLinear = true; break;
	case 6: // TENT_TENT_LOD = linear, linear mip (trilinear)
		minLinear = true; mipLinear = true; break;
	case 7: // CONVOLUTION_2D_LOD0 = anisotropic approx
		anisotropic = true; minLinear = true; mipLinear = true; break;
	default:
		minLinear = false; mipLinear = false; break;
	}

	switch (magFilter) {
	case 1: magLinear = false; break; // BOX = nearest
	case 2: magLinear = true; break;  // TENT = linear
	case 4: anisotropic = true; magLinear = true; break; // CONVOLUTION_2D
	default: magLinear = false; break;
	}

	if (anisotropic) return D3D11_FILTER_ANISOTROPIC;

	return D3D11_ENCODE_BASIC_FILTER(
		minLinear ? D3D11_FILTER_TYPE_LINEAR : D3D11_FILTER_TYPE_POINT,
		magLinear ? D3D11_FILTER_TYPE_LINEAR : D3D11_FILTER_TYPE_POINT,
		mipLinear ? D3D11_FILTER_TYPE_LINEAR : D3D11_FILTER_TYPE_POINT,
		D3D11_FILTER_REDUCTION_TYPE_STANDARD);
}

// ******************************************************************
// * Read PGRAPH texture registers and create D3D11 sampler states.
// * This replaces XboxTextureStates.Apply() for sampler configuration.
// * Registers per stage: TEXADDRESS, TEXFILTER, TEXCTL0, BORDERCOLOR
// ******************************************************************
void CxbxD3D11UpdateSamplersFromPGRAPH(PGRAPHState *pg)
{
	if (!pg) return;

	static ID3D11SamplerState* s_CachedSamplers[4] = {};
	static uint32_t s_CachedTexAddress[4] = {};
	static uint32_t s_CachedTexFilter[4] = {};
	static uint32_t s_CachedTexCtl0[4] = {};
	static uint32_t s_CachedBorderColor[4] = {};

	for (int stage = 0; stage < 4; stage++) {
		uint32_t texAddr   = pg->regs[RI(NV_PGRAPH_TEXADDRESS0 + stage * 4)];
		uint32_t texFilter = pg->regs[RI(NV_PGRAPH_TEXFILTER0 + stage * 4)];
		uint32_t texCtl0   = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];
		uint32_t borderCol = pg->regs[RI(NV_PGRAPH_BORDERCOLOR0 + stage * 4)];

		// Skip if nothing changed
		if (texAddr == s_CachedTexAddress[stage] &&
			texFilter == s_CachedTexFilter[stage] &&
			texCtl0 == s_CachedTexCtl0[stage] &&
			borderCol == s_CachedBorderColor[stage] &&
			s_CachedSamplers[stage] != nullptr) {
			continue;
		}

		s_CachedTexAddress[stage] = texAddr;
		s_CachedTexFilter[stage] = texFilter;
		s_CachedTexCtl0[stage] = texCtl0;
		s_CachedBorderColor[stage] = borderCol;

		// Don't Release the old sampler — it's owned by s_SamplerCache below
		// and may be reused when the same descriptor is requested again.
		s_CachedSamplers[stage] = nullptr;

		// Decode address modes
		unsigned int addrU = GET_MASK(texAddr, NV_PGRAPH_TEXADDRESS0_ADDRU);
		unsigned int addrV = GET_MASK(texAddr, NV_PGRAPH_TEXADDRESS0_ADDRV);
		unsigned int addrP = GET_MASK(texAddr, NV_PGRAPH_TEXADDRESS0_ADDRP);

		// Decode filter modes
		unsigned int minFilter = GET_MASK(texFilter, NV_PGRAPH_TEXFILTER0_MIN);
		unsigned int magFilter = GET_MASK(texFilter, NV_PGRAPH_TEXFILTER0_MAG);

		// LOD bias: 13-bit signed fixed-point (5.8 format: 1 sign + 4 integer, 8 fraction)
		int lodBiasRaw = texFilter & 0x1FFF;
		if (lodBiasRaw & 0x1000) lodBiasRaw |= ~0x1FFF; // sign-extend
		float lodBias = lodBiasRaw / 256.0f;

		// Max anisotropy from TEXCTL0 bits 4-5 (0=1x, 1=2x, 2=4x, 3=8x? or 1=2x)
		unsigned int maxAniso = GET_MASK(texCtl0, NV_PGRAPH_TEXCTL0_0_MAX_ANISOTROPY);
		UINT maxAnisotropy = 1 << maxAniso; // 0→1, 1→2, 2→4, 3→8
		if (maxAnisotropy < 1) maxAnisotropy = 1;

		// LOD clamp from TEXCTL0
		unsigned int minLodRaw = GET_MASK(texCtl0, NV_PGRAPH_TEXCTL0_0_MIN_LOD_CLAMP);
		unsigned int maxLodRaw = GET_MASK(texCtl0, NV_PGRAPH_TEXCTL0_0_MAX_LOD_CLAMP);
		float minLod = minLodRaw / 256.0f;
		float maxLod = maxLodRaw / 256.0f;
		if (maxLod == 0.0f) maxLod = D3D11_FLOAT32_MAX;

		// Border color: ARGB → float4
		float borderColor[4];
		borderColor[0] = ((borderCol >> 16) & 0xFF) / 255.0f; // R
		borderColor[1] = ((borderCol >> 8)  & 0xFF) / 255.0f; // G
		borderColor[2] = (borderCol & 0xFF) / 255.0f;         // B
		borderColor[3] = ((borderCol >> 24) & 0xFF) / 255.0f; // A

		D3D11_SAMPLER_DESC desc = {};
		desc.Filter         = BuildD3D11Filter(minFilter, magFilter);
		desc.AddressU       = MapPGRAPHTexAddress(addrU);
		desc.AddressV       = MapPGRAPHTexAddress(addrV);
		desc.AddressW       = MapPGRAPHTexAddress(addrP);
		desc.MipLODBias     = lodBias;
		desc.MaxAnisotropy  = maxAnisotropy;
		desc.ComparisonFunc = D3D11_COMPARISON_NEVER;
		desc.MinLOD         = minLod;
		desc.MaxLOD         = maxLod;
		std::memcpy(desc.BorderColor, borderColor, sizeof(borderColor));

		// Cache sampler state objects by descriptor to avoid redundant Create calls
		static std::unordered_map<D3D11_SAMPLER_DESC, ID3D11SamplerState*, DescHash, DescEqual> s_SamplerCache;
		auto cacheIt = s_SamplerCache.find(desc);
		if (cacheIt != s_SamplerCache.end()) {
			s_CachedSamplers[stage] = cacheIt->second;
		} else {
			HRESULT hr = g_pD3DDevice->CreateSamplerState(&desc, &s_CachedSamplers[stage]);
			if (SUCCEEDED(hr)) {
				s_SamplerCache[desc] = s_CachedSamplers[stage];
			}
		}
		if (s_CachedSamplers[stage]) {
			g_pD3DDeviceContext->PSSetSamplers(stage, 1, &s_CachedSamplers[stage]);
			g_pD3DDeviceContext->PSSetSamplers(4 + stage, 1, &s_CachedSamplers[stage]);
			g_pD3DDeviceContext->PSSetSamplers(8 + stage, 1, &s_CachedSamplers[stage]);
		}
	}
}

// ******************************************************************
// * Read viewport/scissor from PGRAPH and apply to D3D11.
// * Replaces g_Xbox_Viewport HLE reads with PGRAPH register data.
// *
// * NV2A viewport transform:
// *   vsh_constants[VPSCL] = { Width/2, -Height/2, zScale, 0 }
// *   vsh_constants[VPOFF] = { X+Width/2, Y+Height/2, zOffset, 0 }
// *
// * Scissor: NV2ASurfaceState clip_x/y/width/height (from NV097_SET_SURFACE_CLIP_HORIZONTAL/VERTICAL)
// *   AA factor (NV2ASurfaceState.antiAliasing) is applied first, then the host upscale factor.
// *   This matches xemu's pgraph_apply_anti_aliasing_factor + pgraph_apply_scaling_factor ordering.
// * Depth clip: NV_PGRAPH_ZCLIPMIN / NV_PGRAPH_ZCLIPMAX
// ******************************************************************
void CxbxD3D11UpdateViewportFromPGRAPH(PGRAPHState *pg)
{
	if (!pg) return;

	// Change-detection: viewport/scissor depends on RT dimensions and surface clip state.
	// VPSCL is NOT used here — the D3D11 viewport is always set to full RT size,
	// and the scissor is derived from surface clip registers.
	static uint32_t s_lastPgraphGen = UINT32_MAX;
	static uint32_t s_lastSurfaceGen = UINT32_MAX;

	if (pg->dirty[NV2A_DIRTY_PGRAPH] == s_lastPgraphGen
		&& pg->dirty[NV2A_DIRTY_SURFACE] == s_lastSurfaceGen)
		return;

	s_lastPgraphGen = pg->dirty[NV2A_DIRTY_PGRAPH];
	s_lastSurfaceGen = pg->dirty[NV2A_DIRTY_SURFACE];

	DWORD HostRenderTarget_Width, HostRenderTarget_Height;
	if (!GetHostRenderTargetDimensions(&HostRenderTarget_Width, &HostRenderTarget_Height)) {
		// RT not bound yet — reset dirty tracking so we retry next call.
		s_lastPgraphGen = UINT32_MAX;
		s_lastSurfaceGen = UINT32_MAX;
		return;
	}

	// Scissor from NV2A surface clip registers — matches xemu pgraph_gl/vk_draw_begin.
	// Set for ALL vertex shader modes (PROGRAM and FIXED) to prevent stale scissor state.
	{
		D3D11_VIEWPORT hostViewport;
		hostViewport.TopLeftX = 0;
		hostViewport.TopLeftY = 0;
		hostViewport.Width    = static_cast<float>(HostRenderTarget_Width);
		hostViewport.Height   = static_cast<float>(HostRenderTarget_Height);
		hostViewport.MinDepth = 0.0f;
		hostViewport.MaxDepth = 1.0f;
		CxbxSetViewport(&hostViewport);

		// Apply AA factor first (from NV2ASurfaceState.antiAliasing), then upscale.
		// Matches xemu: pgraph_apply_anti_aliasing_factor then pgraph_apply_scaling_factor.
		auto surf = NV2AGetSurfaceState(pg);
		unsigned int clipX = surf.clipX;
		unsigned int clipY = surf.clipY;
		unsigned int clipW = surf.clipWidth;
		unsigned int clipH = surf.clipHeight;

		// AA factor: matches pgraph_apply_anti_aliasing_factor in EmuNV2A_PGRAPH.cpp
		unsigned int aaFactorX = 1, aaFactorY = 1;
		switch (surf.antiAliasing) {
		case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_CORNER_2:
			aaFactorX = 2;
			aaFactorY = 1;
			break;
		case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_SQUARE_OFFSET_4:
			aaFactorX = 2;
			aaFactorY = 2;
			break;
		default: // NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_1: 1x1
			break;
		}
		clipX *= aaFactorX;
		clipY *= aaFactorY;
		clipW *= aaFactorX;
		clipH *= aaFactorY;

		// Apply host upscale factor
		clipX = static_cast<unsigned int>(clipX * g_RenderUpscaleFactor);
		clipY = static_cast<unsigned int>(clipY * g_RenderUpscaleFactor);
		clipW = static_cast<unsigned int>(clipW * g_RenderUpscaleFactor);
		clipH = static_cast<unsigned int>(clipH * g_RenderUpscaleFactor);

		RECT scissorRect;
		scissorRect.left   = static_cast<LONG>(clipX);
		scissorRect.top    = static_cast<LONG>(clipY);
		scissorRect.right  = std::min(static_cast<LONG>(clipX + clipW), static_cast<LONG>(HostRenderTarget_Width));
		scissorRect.bottom = std::min(static_cast<LONG>(clipY + clipH), static_cast<LONG>(HostRenderTarget_Height));
		CxbxSetScissorRect(&scissorRect);

		g_D3D11RasterizerDesc.ScissorEnable = TRUE;
		g_bD3D11RasterizerStateDirty = true;
	}
}

// ******************************************************************
// * Apply dirty states
// ******************************************************************
void CxbxD3D11ApplyDirtyStates()
{
	LOG_INIT;

	if (g_bD3D11RasterizerStateDirty) {
		auto it = s_RasterizerCache.find(g_D3D11RasterizerDesc);
		if (it != s_RasterizerCache.end()) {
			g_pD3DRasterizerState = it->second;
			g_pD3DDeviceContext->RSSetState(g_pD3DRasterizerState.Get());
		} else {
			HRESULT hr = g_pD3DDevice->CreateRasterizerState(&g_D3D11RasterizerDesc, g_pD3DRasterizerState.ReleaseAndGetAddressOf());
			DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateRasterizerState");
			if (SUCCEEDED(hr)) {
				s_RasterizerCache[g_D3D11RasterizerDesc] = g_pD3DRasterizerState;
				g_pD3DDeviceContext->RSSetState(g_pD3DRasterizerState.Get());
			}
		}
		g_bD3D11RasterizerStateDirty = false;
	}

	if (g_bD3D11DepthStencilStateDirty) {
		auto it = s_DepthStencilCache.find(g_D3D11DepthStencilDesc);
		if (it != s_DepthStencilCache.end()) {
			g_pD3DDepthStencilState = it->second;
			g_pD3DDeviceContext->OMSetDepthStencilState(g_pD3DDepthStencilState.Get(), g_D3D11StencilRef);
		} else {
			HRESULT hr = g_pD3DDevice->CreateDepthStencilState(&g_D3D11DepthStencilDesc, g_pD3DDepthStencilState.ReleaseAndGetAddressOf());
			DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateDepthStencilState");
			if (SUCCEEDED(hr)) {
				s_DepthStencilCache[g_D3D11DepthStencilDesc] = g_pD3DDepthStencilState;
				g_pD3DDeviceContext->OMSetDepthStencilState(g_pD3DDepthStencilState.Get(), g_D3D11StencilRef);
			}
		}
		g_bD3D11DepthStencilStateDirty = false;
	}

	if (g_bD3D11BlendStateDirty) {
		auto it = s_BlendCache.find(g_D3D11BlendDesc);
		if (it != s_BlendCache.end()) {
			g_pD3DBlendState = it->second;
			g_pD3DDeviceContext->OMSetBlendState(g_pD3DBlendState.Get(), g_D3D11BlendFactor, g_D3D11SampleMask);
		} else {
			HRESULT hr = g_pD3DDevice->CreateBlendState(&g_D3D11BlendDesc, g_pD3DBlendState.ReleaseAndGetAddressOf());
			DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateBlendState");
			if (SUCCEEDED(hr)) {
				s_BlendCache[g_D3D11BlendDesc] = g_pD3DBlendState;
				g_pD3DDeviceContext->OMSetBlendState(g_pD3DBlendState.Get(), g_D3D11BlendFactor, g_D3D11SampleMask);
			}
		}
		g_bD3D11BlendStateDirty = false;
	}

	CxbxD3D11FlushVertexShaderConstants();

	// Update GS constant buffer (shared by point sprite and thick line GS)
	// xy = inverse viewport dimensions, z = line width, w = unused
	// Only update when viewport or line width actually changed.
	{
		static float s_LastGSVpWidth = 0.0f;
		static float s_LastGSVpHeight = 0.0f;
		static float s_LastGSLineWidth = 0.0f;

		float vpW = g_CurrentViewportWidth;
		float vpH = g_CurrentViewportHeight;
		if (vpW > 0 && vpH > 0 && g_pD3D11GSConstantBuffer &&
			(vpW != s_LastGSVpWidth || vpH != s_LastGSVpHeight || g_fLineWidth != s_LastGSLineWidth)) {
			s_LastGSVpWidth = vpW;
			s_LastGSVpHeight = vpH;
			s_LastGSLineWidth = g_fLineWidth;
			float gsConstants[4] = { 1.0f / vpW, 1.0f / vpH, g_fLineWidth, 0.0f };
			g_pD3DDeviceContext->UpdateSubresource(g_pD3D11GSConstantBuffer, 0, nullptr, gsConstants, 0, 0);
			g_pD3DDeviceContext->GSSetConstantBuffers(0, 1, &g_pD3D11GSConstantBuffer);
		}
	}

	// Bind or unbind the point sprite geometry shader
	// (Thick line GS is bound at draw time since it depends on primitive type)
	{
		ID3D11GeometryShader* desiredGS = (g_bPointSpriteEnabled && g_pD3D11PointSpriteGS)
			? g_pD3D11PointSpriteGS : nullptr;
		if (desiredGS != s_LastBoundGS) {
			g_pD3DDeviceContext->GSSetShader(desiredGS, nullptr, 0);
			s_LastBoundGS = desiredGS;
		}
	}
}

void CxbxInvalidateGSCache()
{
	s_LastBoundGS = (ID3D11GeometryShader*)~0ull;
}

// ******************************************************************
// * Render target update from PGRAPH surface state
// ******************************************************************

// Track the last PGRAPH surface state we bound, so we only rebind on change.
// Covers offsets, pitches, formats, clip rect, AA, surface type — not just offsets.
static NV2ASurfaceState g_LastBoundSurfaceState = {};

// PGRAPH backbuffer tracking — first color offset bound becomes the backbuffer
static xbox::addr_xt g_PgraphBackBufferOffset = 0;
ID3D11Texture2D* g_pHostPgraphBackBuffer = nullptr;
UINT g_PgraphBackBufferWidth = 0;
UINT g_PgraphBackBufferHeight = 0;

// Implemented after g_PgraphRTCache is defined
void CxbxResetPgraphSurfaceTracking();

// Map NV097 surface color format to DXGI format for host render target creation
static DXGI_FORMAT NV097ColorFormatToDXGI(unsigned int colorFormat)
{
	switch (colorFormat) {
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1R5G5B5_Z1R5G5B5:
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1R5G5B5_O1R5G5B5:
		return DXGI_FORMAT_B5G5R5A1_UNORM;
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_R5G6B5:
		return DXGI_FORMAT_B5G6R5_UNORM;
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_X8R8G8B8_Z8R8G8B8:
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_X8R8G8B8_O8R8G8B8:
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_Z1A7R8G8B8:
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_X1A7R8G8B8_O1A7R8G8B8:
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_A8R8G8B8:
		return DXGI_FORMAT_B8G8R8A8_UNORM;
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_B8:
		return DXGI_FORMAT_R8_UNORM;
	case NV097_SET_SURFACE_FORMAT_COLOR_LE_G8B8:
		return DXGI_FORMAT_R8G8_UNORM;
	default:
		return DXGI_FORMAT_B8G8R8A8_UNORM;
	}
}

// Map NV097 surface zeta format to DXGI format for host depth stencil creation
static DXGI_FORMAT NV097ZetaFormatToDXGI(unsigned int zetaFormat)
{
	switch (zetaFormat) {
	case NV097_SET_SURFACE_FORMAT_ZETA_Z16:
		return DXGI_FORMAT_D16_UNORM;
	case NV097_SET_SURFACE_FORMAT_ZETA_Z24S8:
	default:
		return DXGI_FORMAT_D24_UNORM_S8_UINT;
	}
}

// Unified resource cache: keyed by physical VRAM offset for O(1) lookup.
// Stores render targets, depth stencils, and textures.
// Multiple entries may exist per offset (different format/dimensions).
// Replaces the old RT cache that required O(n) linear scan for offset lookup.
struct ResourceCacheEntry {
	Microsoft::WRL::ComPtr<ID3D11Texture2D> pTexture;
	DXGI_FORMAT format;
	UINT width;
	UINT height;
	uint64_t lastAccessFrame;
	bool isDepthStencil;
};
// Each offset maps to a small vector of entries (typically 1–2).
static std::unordered_map<xbox::addr_xt, std::vector<ResourceCacheEntry>> g_ResourceCache;
static uint64_t g_ResourceCacheFrameCount = 0;
static size_t g_ResourceCacheTotalEntries = 0;

// Elastic eviction watermarks — sized for combined RT + texture population.
static constexpr size_t RESOURCE_CACHE_HIGH_WATERMARK = 256;
static constexpr size_t RESOURCE_CACHE_LOW_WATERMARK = 128;

void CxbxResetPgraphSurfaceTracking()
{
	memset(&g_LastBoundSurfaceState, 0, sizeof(g_LastBoundSurfaceState));
	g_PgraphBackBufferOffset = 0;
	g_pHostPgraphBackBuffer = nullptr;
	g_PgraphBackBufferWidth = 0;
	g_PgraphBackBufferHeight = 0;
	g_ResourceCache.clear();
	g_ResourceCacheTotalEntries = 0;
}

ID3D11Texture2D* CxbxLookupPgraphRTByOffset(xbox::addr_xt offset)
{
	auto it = g_ResourceCache.find(offset);
	if (it != g_ResourceCache.end() && !it->second.empty()) {
		// Return the most recently accessed entry at this offset
		auto& entries = it->second;
		auto* best = &entries[0];
		for (size_t i = 1; i < entries.size(); i++) {
			if (entries[i].lastAccessFrame > best->lastAccessFrame)
				best = &entries[i];
		}
		best->lastAccessFrame = g_ResourceCacheFrameCount;
		return best->pTexture.Get();
	}
	return nullptr;
}

void CxbxPgraphRTCacheEvict()
{
	g_ResourceCacheFrameCount++;

	if (g_ResourceCacheTotalEntries <= RESOURCE_CACHE_HIGH_WATERMARK)
		return;

	// Collect all entries as eviction candidates (exclude pinned backbuffer)
	struct EvictCandidate {
		xbox::addr_xt offset;
		size_t index;
		uint64_t lastAccessFrame;
	};
	std::vector<EvictCandidate> candidates;
	candidates.reserve(g_ResourceCacheTotalEntries);
	for (auto& [offset, entries] : g_ResourceCache) {
		for (size_t i = 0; i < entries.size(); i++) {
			if (entries[i].pTexture.Get() != g_pHostPgraphBackBuffer)
				candidates.push_back({ offset, i, entries[i].lastAccessFrame });
		}
	}

	// Sort by lastAccessFrame ascending (oldest first)
	std::sort(candidates.begin(), candidates.end(),
		[](const EvictCandidate& a, const EvictCandidate& b) {
			return a.lastAccessFrame < b.lastAccessFrame;
		});

	// Evict oldest entries until at or below low watermark
	size_t toEvict = g_ResourceCacheTotalEntries - RESOURCE_CACHE_LOW_WATERMARK;
	size_t evicted = 0;
	for (size_t i = 0; i < candidates.size() && evicted < toEvict; i++) {
		auto mapIt = g_ResourceCache.find(candidates[i].offset);
		if (mapIt == g_ResourceCache.end())
			continue;
		auto& entries = mapIt->second;
		// Find and remove the entry (by matching texture pointer, since indices shift)
		for (auto vecIt = entries.begin(); vecIt != entries.end(); ++vecIt) {
			if (vecIt->lastAccessFrame == candidates[i].lastAccessFrame
				&& vecIt->pTexture.Get() != g_pHostPgraphBackBuffer) {
				entries.erase(vecIt);
				g_ResourceCacheTotalEntries--;
				evicted++;
				break;
			}
		}
		// Remove the offset key if no entries remain
		if (entries.empty())
			g_ResourceCache.erase(mapIt);
	}
}

void CxbxInvalidatePgraphRTBinding()
{
	// Force CxbxD3D11UpdateRenderTargetFromPGRAPH to rebind on the next draw.
	// Must be called after binding a PGRAPH RT as a texture (SRV), because
	// D3D11 automatically unbinds the RTV when the same resource is bound as SRV.
	memset(&g_LastBoundSurfaceState, 0, sizeof(g_LastBoundSurfaceState));
}

// Create a D3D11 render target or depth stencil directly from PGRAPH surface state
static ID3D11Texture2D* CreateHostSurfaceFromPGRAPH(
	xbox::addr_xt offset, DXGI_FORMAT format, UINT width, UINT height, bool isDepthStencil)
{
	if (width == 0 || height == 0)
		return nullptr;

	UINT hostWidth = width * g_RenderUpscaleFactor;
	UINT hostHeight = height * g_RenderUpscaleFactor;

	// Check unified resource cache — reuse if same offset has matching params
	auto it = g_ResourceCache.find(offset);
	if (it != g_ResourceCache.end()) {
		for (auto& entry : it->second) {
			if (entry.format == format
				&& entry.width == hostWidth
				&& entry.height == hostHeight
				&& entry.isDepthStencil == isDepthStencil) {
				entry.lastAccessFrame = g_ResourceCacheFrameCount;
				return entry.pTexture.Get();
			}
		}
	}

	D3D11_TEXTURE2D_DESC desc = {};
	desc.Width = hostWidth;
	desc.Height = hostHeight;
	desc.MipLevels = 1;
	desc.ArraySize = 1;
	desc.Format = format;
	desc.SampleDesc.Count = 1;
	desc.SampleDesc.Quality = 0;
	desc.Usage = D3D11_USAGE_DEFAULT;
	desc.CPUAccessFlags = 0;
	desc.MiscFlags = 0;

	if (isDepthStencil) {
		desc.Format = GetTypelessDepthFormat(format);
		desc.BindFlags = D3D11_BIND_DEPTH_STENCIL | D3D11_BIND_SHADER_RESOURCE;
	} else {
		desc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
	}

	Microsoft::WRL::ComPtr<ID3D11Texture2D> pTexture;
	HRESULT hr = g_pD3DDevice->CreateTexture2D(&desc, nullptr, pTexture.GetAddressOf());
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CreateHostSurfaceFromPGRAPH failed (0x%08X) %ux%u fmt=%u",
			hr, hostWidth, hostHeight, format);
		return nullptr;
	}

	auto* pResult = pTexture.Get();
	g_ResourceCache[offset].push_back({ std::move(pTexture), format, hostWidth, hostHeight,
		g_ResourceCacheFrameCount, isDepthStencil });
	g_ResourceCacheTotalEntries++;

	// Clear newly created depth stencil surfaces to 1.0 (far plane).
	// On real NV2A hardware, newly allocated depth memory contains
	// whatever was there before.  Many games (e.g. MotionBlur) rely on the
	// offscreen RT's depth buffer not being cleared to 0, since they only
	// issue color clears before drawing with depth test LEQUAL.  A D3D11
	// texture starts as all-zeros, causing LEQUAL to reject all fragments.
	if (isDepthStencil && pResult) {
		D3D11_DEPTH_STENCIL_VIEW_DESC dsvDesc = {};
		dsvDesc.Format = GetDepthDSVFormat(format);
		dsvDesc.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2D;
		dsvDesc.Texture2D.MipSlice = 0;
		ID3D11DepthStencilView* pInitDSV = nullptr;
		if (SUCCEEDED(g_pD3DDevice->CreateDepthStencilView(pResult, &dsvDesc, &pInitDSV))) {
			g_pD3DDeviceContext->ClearDepthStencilView(pInitDSV, D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);
			pInitDSV->Release();
		}
	}

	return pResult;
}

void CxbxD3D11UpdateRenderTargetFromPGRAPH(PGRAPHState *pg)
{
	auto surf = NV2AGetSurfaceState(pg);

	// Skip if nothing changed (full state comparison: offsets, pitches, formats, clip, AA)
	if (memcmp(&surf, &g_LastBoundSurfaceState, sizeof(NV2ASurfaceState)) == 0)
		return;

	xbox::addr_xt colorOffset = surf.colorOffset;
	xbox::addr_xt zetaOffset  = surf.zetaOffset;
	xbox::addr_xt prevColorOffset = g_LastBoundSurfaceState.colorOffset;
	xbox::addr_xt prevZetaOffset  = g_LastBoundSurfaceState.zetaOffset;

	UINT rtWidth = surf.clipWidth;
	UINT rtHeight = surf.clipHeight;

	// For swizzled surfaces, actual dimensions are power-of-2 from logWidth/logHeight
	if (surf.surfaceType == 0x2 /*NV097_SET_SURFACE_FORMAT_TYPE_SWIZZLE*/) {
		rtWidth = 1u << surf.logWidth;
		rtHeight = 1u << surf.logHeight;
	}

	// Apply anti-aliasing factor to surface dimensions.
	// On NV2A, AA modes physically expand the surface: 2x doubles width,
	// 4x doubles both. The host RT must match this physical size so that
	// when the game reads the supersampled buffer back as a texture
	// (TEXIMAGERECT reflects the AA-scaled dimensions), our cached RT
	// has the correct dimensions.
	// Matches xemu's pgraph_apply_anti_aliasing_factor.
	switch (surf.antiAliasing) {
	case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_CORNER_2:
		rtWidth *= 2;
		break;
	case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_SQUARE_OFFSET_4:
		rtWidth *= 2;
		rtHeight *= 2;
		break;
	default: // CENTER_1: no scaling
		break;
	}

	// For pitch-linear surfaces, include the clip offset in dimensions.
	// The surface must be large enough to contain the clip region at its offset.
	// Matches xemu's populate_surface_binding_entry.
	if (surf.surfaceType != 0x2 /*NV097_SET_SURFACE_FORMAT_TYPE_SWIZZLE*/) {
		rtWidth += surf.clipX;
		rtHeight += surf.clipY;
	}

	// Color render target (rebind if offset changed, or if format/pitch/clip changed)
	bool colorChanged = (colorOffset != prevColorOffset) ||
		(surf.colorFormat != g_LastBoundSurfaceState.colorFormat) ||
		(surf.colorPitch != g_LastBoundSurfaceState.colorPitch) ||
		(surf.clipWidth != g_LastBoundSurfaceState.clipWidth) ||
		(surf.clipHeight != g_LastBoundSurfaceState.clipHeight);
	if (colorChanged && colorOffset != 0) {
		ID3D11Texture2D *pHostRT = nullptr;
		UINT mipSlice = 0;
		UINT faceIndex = 0;

		// Create host RT directly from PGRAPH surface state
		DXGI_FORMAT colorFmt = NV097ColorFormatToDXGI(surf.colorFormat);
		pHostRT = CreateHostSurfaceFromPGRAPH(colorOffset, colorFmt, rtWidth, rtHeight, false);

		if (pHostRT) {
			CxbxSetRenderTarget(pHostRT, mipSlice, faceIndex);

			// Mark RT pages as GPU-dirty and register for readback.
			// This sets PAGE_NOACCESS so CPU reads trigger a fault → readback.
			uint32_t colorBpp = (colorFmt == DXGI_FORMAT_B8G8R8A8_UNORM) ? 4 :
				(colorFmt == DXGI_FORMAT_R8_UNORM) ? 1 :
				(colorFmt == DXGI_FORMAT_R8G8_UNORM) ? 2 : 2;
			uint32_t colorPitch = surf.colorPitch;
			uint32_t rtSize = colorPitch * rtHeight;
			CxbxPageTrackerMarkGPUDirty(colorOffset, rtSize);
			CxbxPageTrackerRegisterRT(colorOffset, colorPitch,
				rtWidth, rtHeight, colorBpp, surf.surfaceType, pHostRT);
		}

		// Track the backbuffer by matching RT dimensions against presentation parameters.
		// This avoids misidentifying an offscreen surface (shadow map, reflection, etc.)
		// that happens to be rendered before the backbuffer.
		// Also accept half-height surfaces for field rendering (D3DPRESENTFLAG_FIELD),
		// where the Xbox renders 640x240 per field into a 640x480 display.
		// Compare against logical clip dimensions (pre-AA), since presentation
		// parameters reflect the game's logical resolution, not the supersampled size.
		if (g_PgraphBackBufferOffset == 0 && pHostRT) {
			if (surf.clipWidth == g_EmuCDPD.HostPresentationParameters.BackBufferWidth &&
				(surf.clipHeight == g_EmuCDPD.HostPresentationParameters.BackBufferHeight ||
				 surf.clipHeight * 2 == g_EmuCDPD.HostPresentationParameters.BackBufferHeight)) {
				g_PgraphBackBufferOffset = colorOffset;
			}
		}
		if (colorOffset == g_PgraphBackBufferOffset) {
			g_pHostPgraphBackBuffer = pHostRT;
			g_PgraphBackBufferWidth = rtWidth;
			g_PgraphBackBufferHeight = rtHeight;
		}

	}

	// Depth/stencil target (rebind if offset or format changed)
	bool zetaChanged = (zetaOffset != prevZetaOffset) ||
		(surf.zetaFormat != g_LastBoundSurfaceState.zetaFormat) ||
		(surf.zetaPitch != g_LastBoundSurfaceState.zetaPitch);
	if (zetaChanged) {
		if (zetaOffset != 0) {
			ID3D11Texture2D *pHostDS = nullptr;

			// Create host DS directly from PGRAPH state
			DXGI_FORMAT zetaFmt = NV097ZetaFormatToDXGI(surf.zetaFormat);
			pHostDS = CreateHostSurfaceFromPGRAPH(zetaOffset, zetaFmt, rtWidth, rtHeight, true);

			if (pHostDS) {
				// D3D11 requires RTV and DSV dimensions to match.
				// If the new DS has different dimensions from the current color RT
				// (e.g. depth-only shadow pass with 512x512 DS vs 640x480 backbuffer),
				// unbind the color RT for depth-only rendering.  We keep the color
				// offset in g_LastBoundSurfaceState so the color section is correctly
				// skipped on subsequent depth-only draws.
				bool unboundColorForDepthOnly = false;
				if (g_pD3DCurrentHostRenderTarget) {
					D3D11_TEXTURE2D_DESC dsDesc = {}, rtDesc = {};
					pHostDS->GetDesc(&dsDesc);
					g_pD3DCurrentHostRenderTarget->GetDesc(&rtDesc);
					if (dsDesc.Width != rtDesc.Width || dsDesc.Height != rtDesc.Height) {
						// Unbind color RT — depth-only rendering.
						// Do NOT release g_pD3DCurrentRTV: it's owned by g_RTVCache.
						// Releasing here would leave a dangling pointer in the cache,
						// causing use-after-free when the same RT is rebound later.
						g_pD3DCurrentRTV = nullptr;
						g_pD3DCurrentHostRenderTarget = nullptr;
						unboundColorForDepthOnly = true;
					}
				} else if (g_pD3DCurrentRTV == nullptr && g_LastBoundSurfaceState.colorOffset != 0) {
					// Color was already unbound from a previous depth-only pass.
					// Check if the NEW DS matches the color surface dimensions,
					// which means we're transitioning out of depth-only mode.
					D3D11_TEXTURE2D_DESC dsDesc = {};
					pHostDS->GetDesc(&dsDesc);
					UINT expectedW = rtWidth * g_RenderUpscaleFactor;
					UINT expectedH = rtHeight * g_RenderUpscaleFactor;
					if (dsDesc.Width == expectedW && dsDesc.Height == expectedH) {
						// DS now matches color dimensions — force color rebind
						g_LastBoundSurfaceState.colorOffset = 0;
					}
				}
				CxbxSetDepthStencilSurface(pHostDS);
				UpdateDepthStencilFlags(pHostDS);
			}
		} else {
			CxbxSetDepthStencilSurface(nullptr);
		}
	}

	// Commit: record new surface state so subsequent calls see "nothing changed"
	g_LastBoundSurfaceState = surf;
}

// RTV cache: maps (texture pointer, mip slice) to its render target view, avoiding
// redundant CreateRenderTargetView calls for the same texture+mip combination.
std::unordered_map<RTVCacheKey, ID3D11RenderTargetView*, RTVCacheKeyHash> g_RTVCache;

void ClearRTVCache()
{
	for (auto &pair : g_RTVCache) {
		if (pair.second) pair.second->Release();
	}
	g_RTVCache.clear();
	// Reset the current RTV pointer if it was referencing a cached entry
	// (CxbxSetRenderTarget skips Release for cached RTVs, so the cache owns them)
	if (g_pD3DCurrentRTV != nullptr && g_pD3DCurrentRTV != g_pD3DBackBufferView) {
		g_pD3DCurrentRTV = nullptr;
	}
}

// ******************************************************************
// * Thick line GS bind/unbind helpers
// ******************************************************************
static bool CxbxIsLinePrimitive(uint32_t primitiveMode)
{
	return primitiveMode == NV097_SET_BEGIN_END_OP_LINES
	   	|| primitiveMode == NV097_SET_BEGIN_END_OP_LINE_STRIP
	   	|| primitiveMode == NV097_SET_BEGIN_END_OP_LINE_LOOP;
}

void CxbxBindThickLineGS(uint32_t primitiveMode)
{
	if (g_fLineWidth > 1.0f && CxbxIsLinePrimitive(primitiveMode) && g_pD3D11ThickLineGS) {
		g_pD3DDeviceContext->GSSetShader(g_pD3D11ThickLineGS, nullptr, 0);
	}
}

void CxbxUnbindThickLineGS(uint32_t primitiveMode)
{
	if (g_fLineWidth > 1.0f && CxbxIsLinePrimitive(primitiveMode) && g_pD3D11ThickLineGS) {
		// Restore point sprite GS or null
		if (g_bPointSpriteEnabled && g_pD3D11PointSpriteGS) {
			g_pD3DDeviceContext->GSSetShader(g_pD3D11PointSpriteGS, nullptr, 0);
		} else {
			g_pD3DDeviceContext->GSSetShader(nullptr, nullptr, 0);
		}
	}
}

// ******************************************************************
// * Render target / depth-stencil / viewport / scissor binding
// ******************************************************************

HRESULT CxbxSetRenderTarget(ID3D11Texture2D* pHostRenderTarget, UINT mipSlice, UINT arraySlice)
{
	LOG_INIT;
	HRESULT hRet;

	// D3D11 automatically unbinds any SRV referencing the new render target
	// (resource hazard prevention). Mark SRVs dirty so the next draw call
	// rebinds PS SRVs, without forcing expensive full texture re-upload.
	extern void CxbxMarkTextureSRVsDirty();
	CxbxMarkTextureSRVsDirty();
	if (pHostRenderTarget == nullptr) {
		g_pD3DCurrentHostRenderTarget = g_pD3DBackBufferSurface;
		// Don't release g_pD3DCurrentRTV — it's either the default
		// backbuffer view or a cached RTV (owned by g_RTVCache).
		// Releasing here would leave a dangling pointer in the cache,
		// causing use-after-free when the swap chain buffer rotates back.
		g_pD3DCurrentRTV = g_pD3DBackBufferView;
		g_pD3DDeviceContext->OMSetRenderTargets(1, &g_pD3DBackBufferView, g_pD3DDepthStencilView);
		hRet = S_OK;
	} else {
		g_pD3DCurrentHostRenderTarget = pHostRenderTarget;

		RTVCacheKey cacheKey = { pHostRenderTarget, mipSlice, arraySlice };

		// Check RTV cache first
		auto it = g_RTVCache.find(cacheKey);
		if (it != g_RTVCache.end()) {
			if (g_pD3DCurrentRTV != nullptr && g_pD3DCurrentRTV != g_pD3DBackBufferView) {
				// Don't release — it's in the cache
			}
			g_pD3DCurrentRTV = it->second;

			// If DS dimensions don't match the RT, unbind DS to avoid
			// D3D11 silently discarding the RT binding.
			ID3D11DepthStencilView* pDSV = g_pD3DDepthStencilView;
			if (pDSV != nullptr) {
				D3D11_TEXTURE2D_DESC textureDesc = {};
				pHostRenderTarget->GetDesc(&textureDesc);
				ID3D11Resource* dsRes = nullptr;
				pDSV->GetResource(&dsRes);
				if (dsRes) {
					D3D11_TEXTURE2D_DESC dsTexDesc = {};
					((ID3D11Texture2D*)dsRes)->GetDesc(&dsTexDesc);
					dsRes->Release();
					if (dsTexDesc.Width != textureDesc.Width || dsTexDesc.Height != textureDesc.Height) {
						pDSV = nullptr;
					}
				}
			}

			g_pD3DDeviceContext->OMSetRenderTargets(1, &g_pD3DCurrentRTV, pDSV);
			hRet = S_OK;
		} else {
			D3D11_TEXTURE2D_DESC textureDesc = {};
			pHostRenderTarget->GetDesc(&textureDesc);

			D3D11_RENDER_TARGET_VIEW_DESC renderTargetViewDesc{};
			renderTargetViewDesc.Format = textureDesc.Format;
			if (textureDesc.ArraySize > 1) {
				// Cubemap face or texture array — use TEXTURE2DARRAY view
				renderTargetViewDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2DARRAY;
				renderTargetViewDesc.Texture2DArray.MipSlice = mipSlice;
				renderTargetViewDesc.Texture2DArray.FirstArraySlice = arraySlice;
				renderTargetViewDesc.Texture2DArray.ArraySize = 1;
			} else {
				renderTargetViewDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
				renderTargetViewDesc.Texture2D.MipSlice = mipSlice;
			}

			ID3D11RenderTargetView* renderTargetView = nullptr;
			hRet = g_pD3DDevice->CreateRenderTargetView((ID3D11Resource*)pHostRenderTarget, &renderTargetViewDesc, &renderTargetView);
			DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateRenderTargetView");

			if (SUCCEEDED(hRet)) {
				g_RTVCache[cacheKey] = renderTargetView;
				if (g_pD3DCurrentRTV != nullptr && g_pD3DCurrentRTV != g_pD3DBackBufferView) {
					// Don't release — it's in the cache
				}
				g_pD3DCurrentRTV = renderTargetView;

				// If DS dimensions don't match the RT, unbind DS to avoid
				// D3D11 silently discarding the RT binding.
				ID3D11DepthStencilView* pDSV = g_pD3DDepthStencilView;
				if (pDSV != nullptr) {
					D3D11_DEPTH_STENCIL_VIEW_DESC dsDesc;
					pDSV->GetDesc(&dsDesc);
					ID3D11Resource* dsRes = nullptr;
					pDSV->GetResource(&dsRes);
					if (dsRes) {
						D3D11_TEXTURE2D_DESC dsTexDesc = {};
						((ID3D11Texture2D*)dsRes)->GetDesc(&dsTexDesc);
						dsRes->Release();
						if (dsTexDesc.Width != textureDesc.Width || dsTexDesc.Height != textureDesc.Height) {
							pDSV = nullptr; // Unbind mismatched DS
						}
					}
				}

				g_pD3DDeviceContext->OMSetRenderTargets(1, &renderTargetView, pDSV);
			}
		}
	}
	return hRet;
}

void CxbxSetDepthStencilSurface(ID3D11Texture2D* pHostDepthStencil)
{
	// D3D11 unbinds SRVs that conflict with the new DSV resource
	extern void CxbxMarkTextureSRVsDirty();
	CxbxMarkTextureSRVsDirty();

	ID3D11DepthStencilView* pDSV = nullptr;
	if (pHostDepthStencil != nullptr) {
		D3D11_TEXTURE2D_DESC texDesc = {};
		pHostDepthStencil->GetDesc(&texDesc);
		D3D11_DEPTH_STENCIL_VIEW_DESC dsvDesc = {};
		dsvDesc.Format = GetDepthDSVFormat(texDesc.Format);
		dsvDesc.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2D;
		dsvDesc.Texture2D.MipSlice = 0;
		g_pD3DDevice->CreateDepthStencilView(pHostDepthStencil, &dsvDesc, &pDSV);

		// D3D11 requires RTV and DSV dimensions to match, otherwise
		// OMSetRenderTargets silently unbinds both.  If the current RT
		// is a different size (e.g. 256x256 cubemap face vs 640x480 DS),
		// unbind the DSV rather than breaking the RT binding.
		if (pDSV != nullptr && g_pD3DCurrentHostRenderTarget != nullptr) {
			D3D11_TEXTURE2D_DESC rtDesc = {};
			g_pD3DCurrentHostRenderTarget->GetDesc(&rtDesc);
			if (texDesc.Width != rtDesc.Width || texDesc.Height != rtDesc.Height) {
				pDSV->Release();
				pDSV = nullptr;
			}
		}
	}
	if (g_pD3DDepthStencilView) { g_pD3DDepthStencilView->Release(); }
	g_pD3DDepthStencilView = pDSV;
	g_pD3DDeviceContext->OMSetRenderTargets(1, &g_pD3DCurrentRTV, g_pD3DDepthStencilView);
}

ID3D11Texture2D* CxbxGetCurrentRenderTarget()
{
	return g_pD3DCurrentHostRenderTarget;
}

void CxbxSetViewport(D3D11_VIEWPORT *pHostViewport)
{
	g_CurrentViewportWidth = pHostViewport->Width;
	g_CurrentViewportHeight = pHostViewport->Height;
	g_pD3DDeviceContext->RSSetViewports(1, pHostViewport);
}

void CxbxSetScissorRect(CONST RECT *pHostViewportRect)
{
	g_pD3DDeviceContext->RSSetScissorRects(1, pHostViewportRect);
}

