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
#include "devices\video\nv2a.h"        // PGRAPHState, nv2a_regs.h, GET_MASK, RI
#include <algorithm>                    // std::min

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
// ******************************************************************
void CxbxD3D11UpdatePipelineStateFromPGRAPH(PGRAPHState *pg)
{
	if (!pg) return;

	// ---- Blend state from NV_PGRAPH_BLEND (0x1804) ----
	{
		uint32_t blend = pg->regs[RI(NV_PGRAPH_BLEND)];

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

		g_bD3D11BlendStateDirty = true;
	}

	// ---- Blend color from NV_PGRAPH_BLENDCOLOR (0x1808) ----
	{
		uint32_t bc = pg->regs[RI(NV_PGRAPH_BLENDCOLOR)];
		// NV2A BLENDCOLOR is ARGB packed
		g_D3D11BlendFactor[0] = ((bc >> 16) & 0xFF) / 255.0f; // R
		g_D3D11BlendFactor[1] = ((bc >> 8) & 0xFF) / 255.0f;  // G
		g_D3D11BlendFactor[2] = (bc & 0xFF) / 255.0f;         // B
		g_D3D11BlendFactor[3] = ((bc >> 24) & 0xFF) / 255.0f;  // A
	}

	// ---- Color write mask from NV_PGRAPH_CONTROL_0 bits 26-29 ----
	{
		uint32_t ctrl0 = pg->regs[RI(NV_PGRAPH_CONTROL_0)];
		UINT8 writeMask = 0;
		if (ctrl0 & NV_PGRAPH_CONTROL_0_RED_WRITE_ENABLE)   writeMask |= D3D11_COLOR_WRITE_ENABLE_RED;
		if (ctrl0 & NV_PGRAPH_CONTROL_0_GREEN_WRITE_ENABLE) writeMask |= D3D11_COLOR_WRITE_ENABLE_GREEN;
		if (ctrl0 & NV_PGRAPH_CONTROL_0_BLUE_WRITE_ENABLE)  writeMask |= D3D11_COLOR_WRITE_ENABLE_BLUE;
		if (ctrl0 & NV_PGRAPH_CONTROL_0_ALPHA_WRITE_ENABLE) writeMask |= D3D11_COLOR_WRITE_ENABLE_ALPHA;
		g_D3D11BlendDesc.RenderTarget[0].RenderTargetWriteMask = writeMask;
	}

	// ---- Depth state from NV_PGRAPH_CONTROL_0 (0x194C) ----
	{
		uint32_t ctrl0 = pg->regs[RI(NV_PGRAPH_CONTROL_0)];

		g_D3D11DepthStencilDesc.DepthEnable = (ctrl0 & NV_PGRAPH_CONTROL_0_ZENABLE) ? TRUE : FALSE;
		g_D3D11DepthStencilDesc.DepthWriteMask = (ctrl0 & NV_PGRAPH_CONTROL_0_ZWRITEENABLE)
			? D3D11_DEPTH_WRITE_MASK_ALL : D3D11_DEPTH_WRITE_MASK_ZERO;

		// NV2A comparison func values (0-7) = D3D11 comparison func values (1-8) minus 1
		unsigned int zfunc = GET_MASK(ctrl0, NV_PGRAPH_CONTROL_0_ZFUNC);
		g_D3D11DepthStencilDesc.DepthFunc = (D3D11_COMPARISON_FUNC)(zfunc + 1);

		g_bD3D11DepthStencilStateDirty = true;
	}

	// ---- Stencil state from NV_PGRAPH_CONTROL_1 (0x1950) ----
	{
		uint32_t ctrl1 = pg->regs[RI(NV_PGRAPH_CONTROL_1)];

		g_D3D11DepthStencilDesc.StencilEnable = (ctrl1 & NV_PGRAPH_CONTROL_1_STENCIL_TEST_ENABLE) ? TRUE : FALSE;

		unsigned int sfunc = GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_FUNC);
		D3D11_COMPARISON_FUNC stencilFunc = (D3D11_COMPARISON_FUNC)(sfunc + 1);
		g_D3D11DepthStencilDesc.FrontFace.StencilFunc = stencilFunc;
		g_D3D11DepthStencilDesc.BackFace.StencilFunc  = stencilFunc;

		g_D3D11StencilRef = GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_REF);
		g_D3D11DepthStencilDesc.StencilReadMask  = (UINT8)GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_MASK_READ);
		g_D3D11DepthStencilDesc.StencilWriteMask = (UINT8)GET_MASK(ctrl1, NV_PGRAPH_CONTROL_1_STENCIL_MASK_WRITE);
	}

	// ---- Stencil ops from NV_PGRAPH_CONTROL_2 (0x1954) ----
	// NV2A stencil op values (1-8) match D3D11_STENCIL_OP (1-8) directly
	{
		uint32_t ctrl2 = pg->regs[RI(NV_PGRAPH_CONTROL_2)];

		D3D11_STENCIL_OP failOp  = (D3D11_STENCIL_OP)GET_MASK(ctrl2, NV_PGRAPH_CONTROL_2_STENCIL_OP_FAIL);
		D3D11_STENCIL_OP zfailOp = (D3D11_STENCIL_OP)GET_MASK(ctrl2, NV_PGRAPH_CONTROL_2_STENCIL_OP_ZFAIL);
		D3D11_STENCIL_OP zpassOp = (D3D11_STENCIL_OP)GET_MASK(ctrl2, NV_PGRAPH_CONTROL_2_STENCIL_OP_ZPASS);

		g_D3D11DepthStencilDesc.FrontFace.StencilFailOp      = failOp;
		g_D3D11DepthStencilDesc.FrontFace.StencilDepthFailOp = zfailOp;
		g_D3D11DepthStencilDesc.FrontFace.StencilPassOp      = zpassOp;
		g_D3D11DepthStencilDesc.BackFace.StencilFailOp       = failOp;
		g_D3D11DepthStencilDesc.BackFace.StencilDepthFailOp  = zfailOp;
		g_D3D11DepthStencilDesc.BackFace.StencilPassOp       = zpassOp;
	}

	// ---- Rasterizer state from NV_PGRAPH_SETUPRASTER (0x1990) ----
	{
		uint32_t setup = pg->regs[RI(NV_PGRAPH_SETUPRASTER)];

		// Fill mode: PGRAPH FRONTFACEMODE 0=FILL, 1=POINT, 2=LINE
		unsigned int fillMode = GET_MASK(setup, NV_PGRAPH_SETUPRASTER_FRONTFACEMODE);
		switch (fillMode) {
		case NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_FILL:  g_D3D11RasterizerDesc.FillMode = D3D11_FILL_SOLID; break;
		case NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_LINE:  g_D3D11RasterizerDesc.FillMode = D3D11_FILL_WIREFRAME; break;
		case NV_PGRAPH_SETUPRASTER_FRONTFACEMODE_POINT: g_D3D11RasterizerDesc.FillMode = D3D11_FILL_WIREFRAME; break; // no point fill in D3D11
		}

		// Cull mode
		if (!(setup & NV_PGRAPH_SETUPRASTER_CULLENABLE)) {
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

		// Front face winding: PGRAPH bit 23 — 0=CW, 1=CCW
		g_D3D11RasterizerDesc.FrontCounterClockwise = (setup & NV_PGRAPH_SETUPRASTER_FRONTFACE) ? TRUE : FALSE;

		// Line antialiasing
		g_D3D11RasterizerDesc.AntialiasedLineEnable = (setup & NV_PGRAPH_SETUPRASTER_LINESMOOTHENABLE) ? TRUE : FALSE;

		g_bD3D11RasterizerStateDirty = true;
	}

	// ---- Depth bias from NV_PGRAPH_ZOFFSETBIAS / ZOFFSETFACTOR ----
	{
		float zBias; std::memcpy(&zBias, &pg->regs[RI(NV_PGRAPH_ZOFFSETBIAS)], sizeof(float));
		float zFactor; std::memcpy(&zFactor, &pg->regs[RI(NV_PGRAPH_ZOFFSETFACTOR)], sizeof(float));
		// NV2A stores float bias directly; D3D11 DepthBias is an integer
		// scaled by the depth buffer's minimum representable value.
		// For D24: DepthBias * (1 / 2^24).
		g_D3D11RasterizerDesc.DepthBias = static_cast<INT>(zBias * (float)(1 << 24));
		g_D3D11RasterizerDesc.SlopeScaledDepthBias = zFactor;
		g_D3D11RasterizerDesc.DepthBiasClamp = 0.0f;
	}

	// ---- Point sprite enable from NV_PGRAPH_CONTROL_3 ----
	{
		uint32_t ctl3 = pg->regs[RI(NV_PGRAPH_CONTROL_3)];
		g_bPointSpriteEnabled = (ctl3 & NV_PGRAPH_CONTROL_3_POINTPARAMSENABLE) != 0;
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

	// D3D11 filter encoding: bit 4=minLinear, bit 2=magLinear, bit 0=mipLinear
	return (D3D11_FILTER)((minLinear ? 0x10 : 0) | (magLinear ? 0x04 : 0) | (mipLinear ? 0x01 : 0));
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

		// Release old sampler
		if (s_CachedSamplers[stage]) {
			s_CachedSamplers[stage]->Release();
			s_CachedSamplers[stage] = nullptr;
		}

		// Decode address modes
		unsigned int addrU = GET_MASK(texAddr, NV_PGRAPH_TEXADDRESS0_ADDRU);
		unsigned int addrV = GET_MASK(texAddr, NV_PGRAPH_TEXADDRESS0_ADDRV);
		unsigned int addrP = GET_MASK(texAddr, NV_PGRAPH_TEXADDRESS0_ADDRP);

		// Decode filter modes
		unsigned int minFilter = GET_MASK(texFilter, NV_PGRAPH_TEXFILTER0_MIN);
		unsigned int magFilter = GET_MASK(texFilter, NV_PGRAPH_TEXFILTER0_MAG);

		// LOD bias: 13-bit signed fixed-point (8.5 format)
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

		HRESULT hr = g_pD3DDevice->CreateSamplerState(&desc, &s_CachedSamplers[stage]);
		if (SUCCEEDED(hr)) {
			// Bind to all 3 slot groups: base (0-3), 3D (4-7), cube (8-11)
			// All pixel shaders use separate Texture2D/3D/Cube declarations
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
// * Window clip (scissor): NV_PGRAPH_WINDOWCLIPX0/Y0
// * Depth clip: NV_PGRAPH_ZCLIPMIN / NV_PGRAPH_ZCLIPMAX
// ******************************************************************
void CxbxD3D11UpdateViewportFromPGRAPH(PGRAPHState *pg)
{
	if (!pg) return;

	// Read viewport offset and scale from XFCTX constants
	float vpoff[4], vpscl[4];
	for (int i = 0; i < 4; i++) {
		std::memcpy(&vpoff[i], &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPOFF][i], sizeof(float));
		std::memcpy(&vpscl[i], &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL][i], sizeof(float));
	}

	// If the viewport scale constants are zero, PGRAPH hasn't been programmed
	// yet (the Xbox D3D runtime hasn't issued SET_VIEWPORT_OFFSET/SCALE).
	if (vpscl[0] == 0.0f && vpscl[1] == 0.0f) {
		return;
	}

	// Derive Xbox-style viewport rect from NV2A transform constants
	float xboxWidth  = vpscl[0] * 2.0f;
	float xboxHeight = fabsf(vpscl[1]) * 2.0f;
	float xboxX      = vpoff[0] - vpscl[0];
	float xboxY      = vpoff[1] + vpscl[1]; // vpscl[1] is negative

	// Read depth clip range
	float minZ, maxZ;
	std::memcpy(&minZ, &pg->regs[RI(NV_PGRAPH_ZCLIPMIN)], sizeof(float));
	std::memcpy(&maxZ, &pg->regs[RI(NV_PGRAPH_ZCLIPMAX)], sizeof(float));

	// Get host scaling factors (AA + render upscale)
	float aaScaleX, aaScaleY;
	GetMultiSampleScaleRaw(aaScaleX, aaScaleY);
	float Xscale = aaScaleX * g_RenderUpscaleFactor;
	float Yscale = aaScaleY * g_RenderUpscaleFactor;

	DWORD HostRenderTarget_Width, HostRenderTarget_Height;
	if (!GetHostRenderTargetDimensions(&HostRenderTarget_Width, &HostRenderTarget_Height)) {
		return; // can't set viewport without RT dimensions
	}

	// For passthrough mode (XYZRHW/pre-transformed vertices), CMAT is identity
	// because the game provides screen-space positions directly.  For normal FF,
	// CMAT = World*View*Proj*Viewport with large values.  Detect passthrough by
	// checking CMAT ≈ identity rather than VPSCL/VPOFF sign (which is always
	// negative in Y due to the Y-flip, even for normal FF viewports).
	{
		float cmat[4][4];
		for (int row = 0; row < 4; row++)
			std::memcpy(&cmat[row][0], &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_CMAT0 + row][0], 16);
		bool isPassthrough = true;
		for (int r = 0; r < 4 && isPassthrough; r++) {
			for (int c = 0; c < 4 && isPassthrough; c++) {
				float expected = (r == c) ? 1.0f : 0.0f;
				if (fabsf(cmat[r][c] - expected) > 0.01f)
					isPassthrough = false;
			}
		}
		if (isPassthrough) {
			D3D11_VIEWPORT hostViewport;
			hostViewport.TopLeftX = 0;
			hostViewport.TopLeftY = 0;
			hostViewport.Width    = static_cast<float>(HostRenderTarget_Width);
			hostViewport.Height   = static_cast<float>(HostRenderTarget_Height);
			hostViewport.MinDepth = 0.0f;
			hostViewport.MaxDepth = 1.0f;
			CxbxSetViewport(&hostViewport);

			RECT viewportRect = { 0, 0, (LONG)HostRenderTarget_Width, (LONG)HostRenderTarget_Height };
			CxbxSetScissorRect(&viewportRect);
			return;
		}
	}

	// Determine vertex shader mode from PGRAPH CSV0_D register.
	uint32_t pgraphVSMode = GET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_MODE);

	// For FF mode (FIXED and not passthrough), the viewport is already set by
	// UpdateFixedFunctionVertexShaderState() which derives it from CMAT.
	// VPSCL/VPOFF are NOT meaningful for FF mode (the viewport transform is
	// baked into CMAT), so we must not overwrite the FF viewport here.
	if (pgraphVSMode != NV097_SET_TRANSFORM_EXECUTION_MODE_MODE_PROGRAM) {
		return;
	}

	// Programmable VS: full-screen viewport, scissor clips to viewport bounds
	{
		D3D11_VIEWPORT hostViewport;
		hostViewport.TopLeftX = 0;
		hostViewport.TopLeftY = 0;
		hostViewport.Width    = static_cast<float>(HostRenderTarget_Width);
		hostViewport.Height   = static_cast<float>(HostRenderTarget_Height);
		hostViewport.MinDepth = 0.0f;
		hostViewport.MaxDepth = 1.0f;
		CxbxSetViewport(&hostViewport);

		g_D3D11RasterizerDesc.ScissorEnable = TRUE;
		g_bD3D11RasterizerStateDirty = true;

		RECT viewportRect;
		viewportRect.left   = static_cast<LONG>(xboxX * Xscale);
		viewportRect.top    = static_cast<LONG>(xboxY * Yscale);
		viewportRect.right  = std::min(static_cast<LONG>(viewportRect.left + (xboxWidth * Xscale)), (LONG)HostRenderTarget_Width);
		viewportRect.bottom = std::min(static_cast<LONG>(viewportRect.top + (xboxHeight * Yscale)), (LONG)HostRenderTarget_Height);
		CxbxSetScissorRect(&viewportRect);
	}
}

// ******************************************************************
// * Apply dirty states
// ******************************************************************
void CxbxD3D11ApplyDirtyStates()
{
	LOG_INIT;

	if (g_bD3D11RasterizerStateDirty) {
		HRESULT hr = g_pD3DDevice->CreateRasterizerState(&g_D3D11RasterizerDesc, g_pD3DRasterizerState.ReleaseAndGetAddressOf());
		DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateRasterizerState");
		if (SUCCEEDED(hr)) {
			g_pD3DDeviceContext->RSSetState(g_pD3DRasterizerState.Get());
		}
		g_bD3D11RasterizerStateDirty = false;
	}

	if (g_bD3D11DepthStencilStateDirty) {
		HRESULT hr = g_pD3DDevice->CreateDepthStencilState(&g_D3D11DepthStencilDesc, g_pD3DDepthStencilState.ReleaseAndGetAddressOf());
		DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateDepthStencilState");
		if (SUCCEEDED(hr)) {
			g_pD3DDeviceContext->OMSetDepthStencilState(g_pD3DDepthStencilState.Get(), g_D3D11StencilRef);
		}
		g_bD3D11DepthStencilStateDirty = false;
	}

	if (g_bD3D11BlendStateDirty) {
		HRESULT hr = g_pD3DDevice->CreateBlendState(&g_D3D11BlendDesc, g_pD3DBlendState.ReleaseAndGetAddressOf());
		DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateBlendState");
		if (SUCCEEDED(hr)) {
			g_pD3DDeviceContext->OMSetBlendState(g_pD3DBlendState.Get(), g_D3D11BlendFactor, g_D3D11SampleMask);
		}
		g_bD3D11BlendStateDirty = false;
	}

	CxbxD3D11FlushVertexShaderConstants();

	// Update GS constant buffer (shared by point sprite and thick line GS)
	// xy = inverse viewport dimensions, z = line width, w = unused
	{
		D3D11_VIEWPORT vp = {};
		UINT numVP = 1;
		g_pD3DDeviceContext->RSGetViewports(&numVP, &vp);
		if (vp.Width > 0 && vp.Height > 0 && g_pD3D11GSConstantBuffer) {
			float gsConstants[4] = { 1.0f / vp.Width, 1.0f / vp.Height, g_fLineWidth, 0.0f };
			CxbxD3D11UpdateDynamicBuffer(g_pD3D11GSConstantBuffer, gsConstants, sizeof(gsConstants));
			g_pD3DDeviceContext->GSSetConstantBuffers(0, 1, &g_pD3D11GSConstantBuffer);
		}
	}

	// Bind or unbind the point sprite geometry shader
	// (Thick line GS is bound at draw time since it depends on primitive type)
	if (g_bPointSpriteEnabled && g_pD3D11PointSpriteGS) {
		g_pD3DDeviceContext->GSSetShader(g_pD3D11PointSpriteGS, nullptr, 0);
	} else {
		g_pD3DDeviceContext->GSSetShader(nullptr, nullptr, 0);
	}
}

// ******************************************************************
// * Render target update from PGRAPH surface state
// ******************************************************************

// Track the last PGRAPH surface offsets we bound, so we only rebind on change
static xbox::addr_xt g_LastBoundColorOffset = ~0u;
static xbox::addr_xt g_LastBoundZetaOffset  = ~0u;

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

// Cache key for PGRAPH-created render targets / depth stencils
struct PgraphRTKey {
	xbox::addr_xt offset;
	DXGI_FORMAT format;
	UINT width;
	UINT height;

	bool operator==(const PgraphRTKey& other) const {
		return offset == other.offset && format == other.format
			&& width == other.width && height == other.height;
	}
};
struct PgraphRTKeyHash {
	size_t operator()(const PgraphRTKey& k) const {
		size_t h = std::hash<uint32_t>()(k.offset);
		h ^= std::hash<uint32_t>()(k.format) + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<uint32_t>()(k.width)  + 0x9e3779b9 + (h << 6) + (h >> 2);
		h ^= std::hash<uint32_t>()(k.height) + 0x9e3779b9 + (h << 6) + (h >> 2);
		return h;
	}
};
static std::unordered_map<PgraphRTKey, Microsoft::WRL::ComPtr<ID3D11Texture2D>, PgraphRTKeyHash> g_PgraphRTCache;

void CxbxResetPgraphSurfaceTracking()
{
	g_LastBoundColorOffset = ~0u;
	g_LastBoundZetaOffset = ~0u;
	g_PgraphBackBufferOffset = 0;
	g_pHostPgraphBackBuffer = nullptr;
	g_PgraphBackBufferWidth = 0;
	g_PgraphBackBufferHeight = 0;
	g_PgraphRTCache.clear();
}

ID3D11Texture2D* CxbxLookupPgraphRTByOffset(xbox::addr_xt offset)
{
	for (auto& entry : g_PgraphRTCache) {
		if (entry.first.offset == offset)
			return entry.second.Get();
	}
	return nullptr;
}

void CxbxInvalidatePgraphRTBinding()
{
	// Force CxbxD3D11UpdateRenderTargetFromPGRAPH to rebind on the next draw.
	// Must be called after binding a PGRAPH RT as a texture (SRV), because
	// D3D11 automatically unbinds the RTV when the same resource is bound as SRV.
	g_LastBoundColorOffset = ~0u;
	g_LastBoundZetaOffset = ~0u;
}

// Create a D3D11 render target or depth stencil directly from PGRAPH surface state
static ID3D11Texture2D* CreateHostSurfaceFromPGRAPH(
	xbox::addr_xt offset, DXGI_FORMAT format, UINT width, UINT height, bool isDepthStencil)
{
	UINT hostWidth = width * g_RenderUpscaleFactor;
	UINT hostHeight = height * g_RenderUpscaleFactor;

	PgraphRTKey key = { offset, format, hostWidth, hostHeight };
	auto it = g_PgraphRTCache.find(key);
	if (it != g_PgraphRTCache.end())
		return it->second.Get();

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
	g_PgraphRTCache[key] = std::move(pTexture);

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
	xbox::addr_xt colorOffset = pg->surface_color.offset;
	xbox::addr_xt zetaOffset  = pg->surface_zeta.offset;

	// Skip if nothing changed
	if (colorOffset == g_LastBoundColorOffset && zetaOffset == g_LastBoundZetaOffset)
		return;

	UINT rtWidth = pg->surface_shape.clip_width;
	UINT rtHeight = pg->surface_shape.clip_height;

	// Color render target
	if (colorOffset != g_LastBoundColorOffset && colorOffset != 0) {
		ID3D11Texture2D *pHostRT = nullptr;
		UINT mipSlice = 0;
		UINT faceIndex = 0;

		// Try the side-map first (populated by CreateDevice_End for the backbuffer)
		xbox::X_D3DSurface *pXboxRT = CxbxLookupSurfaceByDataAddr(colorOffset);
		if (pXboxRT) {
			pHostRT = GetHostSurface(pXboxRT, D3DUSAGE_RENDERTARGET);

			// Determine mip level and cubemap face for surfaces that are children of a texture
			xbox::X_D3DBaseTexture* pParent = pXboxRT->Parent;
			if (pParent != xbox::zeroptr && pXboxRT->Format == pParent->Format) {
				int face = 0;
				GetSurfaceFaceAndLevelWithinTexture(pXboxRT, pParent, mipSlice, face);
				faceIndex = static_cast<UINT>(face);
				if (GetXboxD3DResourceType(pParent) == xbox::X_D3DRTYPE_CUBETEXTURE) {
					auto pParentHost = (ID3D11Texture2D*)GetHostBaseTexture(pParent, D3DUSAGE_RENDERTARGET);
					if (pParentHost) {
						pHostRT = pParentHost;
					}
				}
			}
		} else {
			// No Xbox surface registered — create host RT directly from PGRAPH state
			DXGI_FORMAT colorFmt = NV097ColorFormatToDXGI(pg->surface_shape.color_format);
			pHostRT = CreateHostSurfaceFromPGRAPH(colorOffset, colorFmt, rtWidth, rtHeight, false);
		}

		if (pHostRT) {
			CxbxSetRenderTarget(pHostRT, mipSlice, faceIndex);
		}

		// Track the first color offset as the backbuffer; update pointer on re-bind
		if (g_PgraphBackBufferOffset == 0 && pHostRT) {
			g_PgraphBackBufferOffset = colorOffset;
		}
		if (colorOffset == g_PgraphBackBufferOffset) {
			g_pHostPgraphBackBuffer = pHostRT;
			g_PgraphBackBufferWidth = rtWidth;
			g_PgraphBackBufferHeight = rtHeight;
		}

		g_LastBoundColorOffset = colorOffset;
	}

	// Depth/stencil target
	if (zetaOffset != g_LastBoundZetaOffset) {
		if (zetaOffset != 0) {
			ID3D11Texture2D *pHostDS = nullptr;

			xbox::X_D3DSurface *pXboxDS = CxbxLookupSurfaceByDataAddr(zetaOffset);
			if (pXboxDS) {
				pHostDS = GetHostSurface(pXboxDS, D3DUSAGE_DEPTHSTENCIL);
			} else {
				// No Xbox surface registered — create host DS directly from PGRAPH state
				DXGI_FORMAT zetaFmt = NV097ZetaFormatToDXGI(pg->surface_shape.zeta_format);
				pHostDS = CreateHostSurfaceFromPGRAPH(zetaOffset, zetaFmt, rtWidth, rtHeight, true);
			}

			if (pHostDS) {
				// D3D11 requires RTV and DSV dimensions to match.
				// If the new DS has different dimensions from the current color RT
				// (e.g. depth-only shadow pass with 512x512 DS vs 640x480 backbuffer),
				// unbind the color RT and invalidate tracking so it gets rebound
				// when the color offset changes on the next pass.
				if (g_pD3DCurrentHostRenderTarget) {
					D3D11_TEXTURE2D_DESC dsDesc = {}, rtDesc = {};
					pHostDS->GetDesc(&dsDesc);
					g_pD3DCurrentHostRenderTarget->GetDesc(&rtDesc);
					if (dsDesc.Width != rtDesc.Width || dsDesc.Height != rtDesc.Height) {
						// Unbind color RT — depth-only rendering
						if (g_pD3DCurrentRTV && g_pD3DCurrentRTV != g_pD3DBackBufferView) {
							g_pD3DCurrentRTV->Release();
						}
						g_pD3DCurrentRTV = nullptr;
						g_pD3DCurrentHostRenderTarget = nullptr;
						g_LastBoundColorOffset = ~0u; // Force rebind on next color change
					}
				}
				CxbxSetDepthStencilSurface(pHostDS);
				UpdateDepthStencilFlags(pHostDS);
			}
		} else {
			CxbxSetDepthStencilSurface(nullptr);
		}
		g_LastBoundZetaOffset = zetaOffset;
	}
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
static bool CxbxIsLinePrimitive(xbox::X_D3DPRIMITIVETYPE type)
{
	return type == xbox::X_D3DPT_LINELIST
	   	|| type == xbox::X_D3DPT_LINESTRIP
	   	|| type == xbox::X_D3DPT_LINELOOP;
}

void CxbxBindThickLineGS(xbox::X_D3DPRIMITIVETYPE type)
{
	if (g_fLineWidth > 1.0f && CxbxIsLinePrimitive(type) && g_pD3D11ThickLineGS) {
		g_pD3DDeviceContext->GSSetShader(g_pD3D11ThickLineGS, nullptr, 0);
	}
}

void CxbxUnbindThickLineGS(xbox::X_D3DPRIMITIVETYPE type)
{
	if (g_fLineWidth > 1.0f && CxbxIsLinePrimitive(type) && g_pD3D11ThickLineGS) {
		// Restore point sprite GS or null
		if (g_bPointSpriteEnabled && g_pD3D11PointSpriteGS) {
			g_pD3DDeviceContext->GSSetShader(g_pD3D11PointSpriteGS, nullptr, 0);
		} else {
			g_pD3DDeviceContext->GSSetShader(nullptr, nullptr, 0);
		}
	}
}

