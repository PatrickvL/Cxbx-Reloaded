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
// *  (c) 2002-2003 kingofc <kingofc@freenet.de>
// *  2020 PatrickvL
// *
// *  All rights reserved
// *
// ******************************************************************
#define LOG_PREFIX CXBXR_MODULE::D3D8

#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\Rendering\RenderGlobals.h"
#include "core\hle\D3D8\Rendering\Backend\Shading\Shader.h"
#include "core\hle\D3D8\XbPixelShader.h"
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbD3D8Logging.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\kernel\init\CxbxKrnl.h"
// Texture format fixup constants (must match ApplyTexFmtFixup() in CxbxPixelShaderFunctions.hlsli)
static constexpr float TEXFMTFIXUP_IDENTITY = 0.0f;
static constexpr float TEXFMTFIXUP_GBAR     = 1.0f; // B8G8R8A8 uploaded as R8G8B8A8
static constexpr float TEXFMTFIXUP_ABGR     = 2.0f; // A8B8G8R8 uploaded as R8G8B8A8
static constexpr float TEXFMTFIXUP_LUM      = 3.0f; // Luminance: R8→(R,R,R,1)
static constexpr float TEXFMTFIXUP_ALUM     = 4.0f; // Alpha-luminance: R8G8→(R,R,R,G)
static constexpr float TEXFMTFIXUP_OPAQUEA  = 5.0f; // X8R8G8B8/X1R5G5B5: force alpha to 1.0
#include "devices\Xbox.h"              // For extern NV2ADevice* g_NV2A
#include "devices\video\nv2a.h"        // For NV2ADevice::GetDeviceState(), NV2AState, PGRAPHState, nv2a_regs.h
#include <assert.h>
#include <process.h>
#include <unordered_map>
#include "Rendering\RenderStates.h"
#include "Rendering\TextureStates.h"
#include <wrl/client.h>
#include <cstring> // For std::memcpy
#include "Rendering\Backend\Backend_D3D11.h"
#include "Rendering\Backend\Backend_D3D11_Internal.h"
#include "Rendering\Backend\Shading\PixelShaderCache.h"
#include "Rendering\Backend\Backend_D3D11_Profiler.h"

float AsFloat(uint32_t value)
{
	float f; std::memcpy(&f, &value, sizeof(f)); return f;
}

// Determines the Cxbx ColorSign requirement, as handled in the HLSL shaders by PerformColorSign()
float CxbxComponentColorSignFromXboxAndHost(bool XboxMarksComponentSigned, bool HostComponentIsSigned)
{
	// Equal "signedness" between Xbox and host implies we must not convert the component scale :
	if (XboxMarksComponentSigned == HostComponentIsSigned)
		return 0.0f;

	// Xbox wants the components to be signed (even though host has them unsigned)
	if (XboxMarksComponentSigned)
		return 1.0f; // Mark the component for scaling from unsigned_to_signed

	// Xbox doesn't want signed values, but host has them signed :
	return -1.0f; // Mark the component for scaling from signed_to_unsigned
}

float CxbxGetTexFmtFixup(int stage_nr)
{

	// Resolve texture via PGRAPH offset → side-map to avoid racing g_pXbox_SetTexture[].
	xbox::X_D3DBaseTexture *pXboxTex = xbox::zeroptr;
	{
		auto pg_ff = &(g_NV2A->GetDeviceState()->pgraph);
		uint32_t texCtl = pg_ff->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage_nr * 4)];
		if (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) {
			uint32_t texOffset = pg_ff->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage_nr * 4)];
			if (texOffset != 0)
				pXboxTex = CxbxLookupTextureByDataAddr(texOffset);
		}
	}
	if (pXboxTex == xbox::zeroptr)
		pXboxTex = g_pXbox_SetTexture[stage_nr]; // fallback
	if (pXboxTex == xbox::zeroptr)
		return TEXFMTFIXUP_IDENTITY;

	xbox::X_D3DFORMAT xboxFmt = GetXboxPixelContainerFormat((xbox::X_D3DPixelContainer*)pXboxTex);
	switch (xboxFmt) {
	case xbox::X_D3DFMT_X8R8G8B8:
	case xbox::X_D3DFMT_LIN_X8R8G8B8:
	case xbox::X_D3DFMT_X1R5G5B5:
	case xbox::X_D3DFMT_LIN_X1R5G5B5:
		return TEXFMTFIXUP_OPAQUEA;
	case xbox::X_D3DFMT_L8:
	case xbox::X_D3DFMT_LIN_L8:
	case xbox::X_D3DFMT_L16:
	case xbox::X_D3DFMT_LIN_L16:
		return TEXFMTFIXUP_LUM;
	case xbox::X_D3DFMT_A8L8:
	case xbox::X_D3DFMT_LIN_A8L8:
		return TEXFMTFIXUP_ALUM;
	// B8G8R8A8 and R8G8B8A8 use GBAR/ABGR swizzles when uploaded raw
	// (requires corresponding skip of CPU conversion in HostResourceCreate.cpp)
	// Exception: render targets use B8G8R8A8_UNORM and don't need a swizzle.
	case xbox::X_D3DFMT_B8G8R8A8:
	case xbox::X_D3DFMT_LIN_B8G8R8A8:
	{
		// Check if the host texture is B8G8R8A8_UNORM (render target) — data is already correct
		auto key = GetHostResourceKey(pXboxTex, stage_nr);
		auto& cache = GetResourceCache(key);
		auto it = cache.find(key);
		if (it != cache.end() && it->second.HostFormat == EMUFMT_A8R8G8B8)
			return TEXFMTFIXUP_IDENTITY;
		return TEXFMTFIXUP_GBAR;
	}
	case xbox::X_D3DFMT_R8G8B8A8:
	case xbox::X_D3DFMT_LIN_R8G8B8A8:
	{
		auto key = GetHostResourceKey(pXboxTex, stage_nr);
		auto& cache = GetResourceCache(key);
		auto it = cache.find(key);
		if (it != cache.end() && it->second.HostFormat == EMUFMT_A8R8G8B8)
			return TEXFMTFIXUP_IDENTITY;
		return TEXFMTFIXUP_ABGR;
	}
	default:
		break;
	}
	return TEXFMTFIXUP_IDENTITY;
}

D3DXCOLOR CxbxCalcColorSign(int stage_nr)
{
	// Read COLORSIGN from PGRAPH TEXFILTER register (bits 28-31: ASIGNED, RSIGNED, GSIGNED, BSIGNED).
	// The Xbox D3D runtime writes X_D3DTSS_COLORSIGN bits directly into the TEXFILTER register,
	// and the bit positions match: ASIGNED=bit28, RSIGNED=bit29, GSIGNED=bit30, BSIGNED=bit31.
	auto pg = &(g_NV2A->GetDeviceState()->pgraph);
	uint32_t texFilter = pg->regs[RI(NV_PGRAPH_TEXFILTER0 + stage_nr * 4)];
	DWORD XboxColorSign = texFilter & 0xF0000000; // Extract sign bits (matches X_D3DTSIGN layout)

	{ // This mimics behaviour of XDK LazySetShaderStageProgram, which we bypass due to our drawing patches without trampolines.
		// When bump environment mapping is enabled, check the shader stage program from PGRAPH.
		// COLOROP >= BUMPENVMAP maps to shader stage mode BUMPENVMAP(6) or BUMPENVMAP_LUMINANCE(7).
		static const uint32_t stageMasks[4] = {
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE0,
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE1,
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE2,
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE3
		};
		uint32_t shaderProg = pg->regs[RI(NV_PGRAPH_SHADERPROG)];
		uint32_t stageMode = GET_MASK(shaderProg, stageMasks[stage_nr]);
		// BUMPENVMAP=6, BUMPENVMAP_LUMINANCE=7 (same across all stages)
		if (stageMode == 6 || stageMode == 7)
			// Always mark the blue (alias for U) and green (alias for V) color channels as signed:
			XboxColorSign |= xbox::X_D3DTSIGN_GSIGNED | xbox::X_D3DTSIGN_BSIGNED;
	}

#if 0 // When this block is enabled, XDK samples BumpEarth and BumpLens turn red-ish, so keep this off for now...
	// Check if the pixel shader specifies bump mapping for this stage (TODO : How to handle this with the fixed function shader?)
	DWORD PSTextureModes = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_PSTEXTUREMODES);
	PS_TEXTUREMODES StageTextureMode = (PS_TEXTUREMODES)((PSTextureModes >> (stage_nr * 5)) & PS_TEXTUREMODES_MASK);
	if (StageTextureMode == PS_TEXTUREMODES_BUMPENVMAP || StageTextureMode == PS_TEXTUREMODES_BUMPENVMAP_LUM)
		XboxColorSign |= xbox::X_D3DTSIGN_GSIGNED | xbox::X_D3DTSIGN_BSIGNED;

#endif
	// Host D3DFMT's with one or more signed components : D3DFMT_V8U8, D3DFMT_Q8W8V8U8, D3DFMT_V16U16, D3DFMT_Q16W16V16U16, D3DFMT_CxV8U8
	DXGI_FORMAT H/*ostTextureFormat*/ = g_HostTextureFormats[stage_nr];
	// Guard: if the host format is unknown (stage not yet populated), skip all signed checks.
	// In D3D11, EMUFMT_L6V5U5 == DXGI_FORMAT_NOT_AVAILABLE == DXGI_FORMAT_UNKNOWN == 0,
	// so without this guard, uninitialized stages falsely match L6V5U5 signed detection.
	if (H == EMUFMT_UNKNOWN) {
		D3DXCOLOR zero(0, 0, 0, 0);
		return zero; // No signed conversion for unknown/uninitialized stages
	}
	// See https://docs.microsoft.com/en-us/windows/win32/direct3d9/bump-map-pixel-formats
	// No need to check for unused formats : D3DFMT_Q16W16V16U16, D3DFMT_CxV8U8, D3DFMT_A2W10V10U10
#if 0 // Original signed-ness checking code gave effectively this :
	// Host format     | Signed components
	// ----------------+------------------
	// D3DFMT_Q8W8V8U8 | A,R,G,B
	// D3DFMT_L6V5U5   |     G,B
	// D3DFMT_V8U8     |   R,G
	// D3DFMT_V16U16   |   R,G
	// D3DFMT_X8L8V8U8 |   R,G
	bool HostTextureFormatIsSignedForA = (H == EMUFMT_Q8W8V8U8);
	bool HostTextureFormatIsSignedForR = (H == EMUFMT_Q8W8V8U8)                         || (H == EMUFMT_V8U8) || (H == EMUFMT_V16U16) || (H == EMUFMT_X8L8V8U8);
	bool HostTextureFormatIsSignedForG = (H == EMUFMT_Q8W8V8U8) || (H == EMUFMT_L6V5U5) || (H == EMUFMT_V8U8) || (H == EMUFMT_V16U16) || (H == EMUFMT_X8L8V8U8);
	bool HostTextureFormatIsSignedForB = (H == EMUFMT_Q8W8V8U8) || (H == EMUFMT_L6V5U5);
#else // New, as experimentally discovered by medievil :
	// Host format     | Signed components
	// ----------------+------------------
	// D3DFMT_Q8W8V8U8 | A,R,G,B
	// D3DFMT_L6V5U5   | A,R
	// D3DFMT_V8U8     |   R,G
	// D3DFMT_V16U16   |   R,G
	// D3DFMT_X8L8V8U8 |   R,G
	// TODO : Verify D3DFMT_L6V5U5 indeed maps to A,R (instead of G,B).
	// If not, research why this (then incorret) change *does* improve both BumpEarth samples
	// (while keeping BumpLens and JSFR boost dash effect working). Perhaps duplicate signed range conversion in the shader?
	bool HostTextureFormatIsSignedForA = (H == EMUFMT_Q8W8V8U8) || (H == EMUFMT_L6V5U5);
	bool HostTextureFormatIsSignedForR = (H == EMUFMT_Q8W8V8U8) || (H == EMUFMT_L6V5U5) || (H == EMUFMT_V8U8) || (H == EMUFMT_V16U16) || (H == EMUFMT_X8L8V8U8);
	bool HostTextureFormatIsSignedForG = (H == EMUFMT_Q8W8V8U8)                         || (H == EMUFMT_V8U8) || (H == EMUFMT_V16U16) || (H == EMUFMT_X8L8V8U8);
	bool HostTextureFormatIsSignedForB = (H == EMUFMT_Q8W8V8U8);
#endif
	D3DXCOLOR CxbxColorSign;
	CxbxColorSign.r = CxbxComponentColorSignFromXboxAndHost(XboxColorSign & xbox::X_D3DTSIGN_RSIGNED, HostTextureFormatIsSignedForR); // Maps to COLORSIGN.r
	CxbxColorSign.g = CxbxComponentColorSignFromXboxAndHost(XboxColorSign & xbox::X_D3DTSIGN_GSIGNED, HostTextureFormatIsSignedForG); // Maps to COLORSIGN.g
	CxbxColorSign.b = CxbxComponentColorSignFromXboxAndHost(XboxColorSign & xbox::X_D3DTSIGN_BSIGNED, HostTextureFormatIsSignedForB); // Maps to COLORSIGN.b
	CxbxColorSign.a = CxbxComponentColorSignFromXboxAndHost(XboxColorSign & xbox::X_D3DTSIGN_ASIGNED, HostTextureFormatIsSignedForA); // Maps to COLORSIGN.a
	return CxbxColorSign;
}

static ID3D11PixelShader* g_pActivePixelShader = nullptr; // TODO : Reset when device resets!

void CxbxInvalidateActivePixelShader()
{
	// Called after the blit/present path which bypasses CxbxSetPixelShader
	// and binds its own PS directly. Without this, the next CxbxSetPixelShader
	// call would skip the bind because g_pActivePixelShader still holds the
	// pre-blit pointer even though the device now has the blit PS bound.
	g_pActivePixelShader = nullptr;
}

void CxbxSetPixelShader(ID3D11PixelShader* pPixelShader)
{
	// Here no call to (PS)GetPixelShader, but our own state tracking; See https://gamedev.stackexchange.com/a/88117
	if (g_pActivePixelShader == pPixelShader)
		return;

	// Switch to the converted pixel shader (if it's any different from our currently active
	// pixel shader, to avoid many unnecessary state changes on the local side).
	CxbxRawSetPixelShader(pPixelShader);
	g_pActivePixelShader = pPixelShader;
}

// Global copy of last-built aux CB, readable by the PS JIT for state hashing
PSAuxCBLayout g_LastPSAuxCB = {};

// Upload PGRAPH register combiner state to GPU buffers.
// PGRAPH is always authoritative — Xbox native D3D code pushes all combiner,
// texture, and fog state through PFIFO → PGRAPH before each draw.
void CxbxD3D11UploadRCInterpreterState()
{
	if (!g_pD3D11RCInterpreterAuxCB || !g_pD3D11PGRegsBuf)
		return;

	// PGRAPH source (populated by the puller thread via pushbuffer methods)
	PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;

	// --- Upload raw PGRAPH regs[] to the StructuredBuffer<uint> SRV ---
	// Only re-upload when regs actually changed (generation counter bumped
	// by nv097_dispatch_method on any register write).
	// NOTE: Both JIT and interpreter shaders read dynamic constants (C0/C1,
	// fog color, bump matrices) from this SRV at runtime, so upload is required.
	{
		static uint32_t s_LastRegsGeneration = ~0u;
		if (pg->regs_generation != s_LastRegsGeneration) {
			s_LastRegsGeneration = pg->regs_generation;
			CxbxD3D11UpdateDynamicBuffer(g_pD3D11PGRegsBuf, pg->regs, sizeof(pg->regs));
		}
	}
	// Bind the regs SRV to PS t12 (skip if already bound — pointer never changes)
	{
		static bool s_RegsSRVBound = false;
		if (!s_RegsSRVBound) {
			g_pD3DDeviceContext->PSSetShaderResources(CXBX_D3D11_PS_PGREGS_SRV_SLOT, 1, &g_pD3D11PGRegsSRV);
			s_RegsSRVBound = true;
		}
	}

	// --- Build the auxiliary cbuffer (software-computed fields only) ---
	// Skip rebuild if regs_generation hasn't changed (aux depends only on regs[])
	{
		static uint32_t s_LastAuxGeneration = ~0u;
		if (pg->regs_generation == s_LastAuxGeneration)
			return; // Aux CB and regs SRV are still current
		s_LastAuxGeneration = pg->regs_generation;
	}
	PSAuxCBLayout aux = {};

	// PSTextureModes: always from PGRAPH SHADERPROG
	DWORD psTextureModes = pg->regs[RI(NV_PGRAPH_SHADERPROG)];

	// --- AdjustTextureModes: fixup cubemap/volume texture modes ---
	// PGRAPH is always authoritative — SHADERPROG already contains correctly
	// adjusted texture modes from the Xbox D3D runtime push buffer.
	// Texture type is derived from PGRAPH TEXFMT0 (CUBEMAPENABLE +
	// DIMENSIONALITY) to avoid racing the game thread.
	{
		for (int i = 0; i < xbox::X_D3DTS_STAGECOUNT; i++) {
			uint32_t mode = (psTextureModes >> (i * 5)) & 0x1Fu;
			uint32_t clearMask = ~(0x1Fu << (i * 5));

			// Derive texture type from PGRAPH registers
			xbox::X_D3DRESOURCETYPE texType = xbox::X_D3DRTYPE_NONE;
			{
				uint32_t texCtl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
				if (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) {
					uint32_t texFmt = pg->regs[RI(NV_PGRAPH_TEXFMT0 + i * 4)];
					if (texFmt & NV_PGRAPH_TEXFMT0_CUBEMAPENABLE)
						texType = xbox::X_D3DRTYPE_CUBETEXTURE;
					else if (((texFmt & NV_PGRAPH_TEXFMT0_DIMENSIONALITY) >> 6) > 2)
						texType = xbox::X_D3DRTYPE_VOLUMETEXTURE;
					else
						texType = xbox::X_D3DRTYPE_TEXTURE;
				}
			}

			if (texType == xbox::X_D3DRTYPE_CUBETEXTURE && mode == PS_TEXTUREMODES_PROJECT2D) {
				psTextureModes = (psTextureModes & clearMask) | ((uint32_t)PS_TEXTUREMODES_CUBEMAP << (i * 5));
			}
			else if (texType == xbox::X_D3DRTYPE_CUBETEXTURE && mode == PS_TEXTUREMODES_DOT_STR_3D) {
				psTextureModes = (psTextureModes & clearMask) | ((uint32_t)PS_TEXTUREMODES_DOT_STR_CUBE << (i * 5));
			}
		}
	}
	aux.PSTextureModes.value = psTextureModes;

	// --- AdjustFinalCombiner: synthesize final combiner when not explicitly defined ---
	{
		uint32_t fcABCD = pg->regs[RI(NV_PGRAPH_COMBINESPECFOG0)];
		uint32_t fcEFG  = pg->regs[RI(NV_PGRAPH_COMBINESPECFOG1)];

		bool hasFinalCombiner = (fcABCD != 0) || (fcEFG != 0);
		if (!hasFinalCombiner) {
			bool fogEnable = (pg->regs[RI(NV_PGRAPH_CONTROL_3)] & NV_PGRAPH_CONTROL_3_FOGENABLE) != 0;
			bool specularEnable = (pg->regs[RI(NV_PGRAPH_CSV0_C)] & NV_PGRAPH_CSV0_C_SPECULAR_ENABLE) != 0;

			uint32_t regA = PS_REGISTER_FOG | PS_CHANNEL_ALPHA;
			uint32_t regB = PS_REGISTER_R0;
			uint32_t regC = fogEnable ? PS_REGISTER_FOG : PS_REGISTER_R0;
			uint32_t regD = specularEnable ? PS_REGISTER_V1 : PS_REGISTER_ZERO;
			fcABCD = (regA << 24) | (regB << 16) | (regC << 8) | regD;

			uint32_t regE = PS_REGISTER_ZERO;
			uint32_t regF = PS_REGISTER_ZERO;
			uint32_t regG = PS_REGISTER_R0 | PS_CHANNEL_ALPHA;
			fcEFG = (regE << 24) | (regF << 16) | (regG << 8);
		}

		aux.PSFinalCombinerInputsABCD.value = fcABCD;
		aux.PSFinalCombinerInputsEFG.value  = fcEFG;
	}

	// Color sign conversion — per-stage
	for (int stage = 0; stage < 4; stage++) {
		D3DXCOLOR cs = CxbxCalcColorSign(stage);
		aux.ColorSign[stage] = { cs.r, cs.g, cs.b, cs.a };
	}

	// Texture format channel fixup per stage
	aux.TexFmtFixup = { CxbxGetTexFmtFixup(0), CxbxGetTexFmtFixup(1),
	                     CxbxGetTexFmtFixup(2), CxbxGetTexFmtFixup(3) };

	// Color key per stage — read from PGRAPH (authoritative, no HLE dependency)
	for (int i = 0; i < 4; i++) {
		// COLORKEYOP is stored in TEXCTL0 bits 0-1 (COLORKEYMODE)
		uint32_t texCtl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
		uint32_t colorKeyMode = texCtl & NV_PGRAPH_TEXCTL0_0_COLORKEYMODE;
		aux.ColorKeyOp[i] = { static_cast<float>(colorKeyMode), 0.0f, 0.0f, 0.0f };

		// COLORKEYCOLOR is stored in NV_PGRAPH_COLORKEYCOLOR0..3
		D3DXCOLOR ckc(pg->regs[RI(NV_PGRAPH_COLORKEYCOLOR0 + i * 4)]);
		aux.ColorKeyColor[i] = { ckc.r, ckc.g, ckc.b, ckc.a };
	}

	// Alpha kill per stage — read from PGRAPH TEXCTL0 ALPHAKILLEN bit
	aux.AlphaKill = {
		static_cast<float>((pg->regs[RI(NV_PGRAPH_TEXCTL0_0)] & NV_PGRAPH_TEXCTL0_0_ALPHAKILLEN) ? 1 : 0),
		static_cast<float>((pg->regs[RI(NV_PGRAPH_TEXCTL0_1)] & NV_PGRAPH_TEXCTL0_0_ALPHAKILLEN) ? 1 : 0),
		static_cast<float>((pg->regs[RI(NV_PGRAPH_TEXCTL0_2)] & NV_PGRAPH_TEXCTL0_0_ALPHAKILLEN) ? 1 : 0),
		static_cast<float>((pg->regs[RI(NV_PGRAPH_TEXCTL0_3)] & NV_PGRAPH_TEXCTL0_0_ALPHAKILLEN) ? 1 : 0)
	};

	// Fog info: x=tableMode (from PGRAPH FOG_MODE), y/z/w unused by RC interpreter.
	// FogColor is read directly from g_PGRegs[] in the shader.
	{
		unsigned int fogMode = GET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_3)], NV_PGRAPH_CONTROL_3_FOG_MODE);
		aux.FogInfo = { static_cast<float>(fogMode), 0.0f, 0.0f, 0.0f };
		aux.FogEnable.value = (pg->regs[RI(NV_PGRAPH_CONTROL_3)] & NV_PGRAPH_CONTROL_3_FOGENABLE) ? 1u : 0u;
	}

	// Front-face factor for two-sided lighting — sourced from PGRAPH
	{
		float ff = 0.0f;
		uint32_t csv0c = pg->regs[RI(NV_PGRAPH_CSV0_C)];
		// NV2A LIGHT_MODEL_TWO_SIDE_ENABLE: bit 29 of CSV0_C
		bool twoSided = (csv0c & 0x20000000u) != 0;
		if (twoSided) {
			// NV_PGRAPH_SETUPRASTER_FRONTFACE: bit 23 — 0=CW, 1=CCW
			uint32_t setup = pg->regs[RI(NV_PGRAPH_SETUPRASTER)];
			bool ccwFront = (setup & NV_PGRAPH_SETUPRASTER_FRONTFACE) != 0;
			ff = ccwFront ? -1.0f : 1.0f;
		}
		aux.FrontFaceInfo = { ff, 0.0f, 0.0f, 0.0f };
	}

	// Shadow compare: per-stage flag indicating a depth texture is bound.
	// When active, the pixel shader compares the R texcoord against the
	// sampled depth value using NV_PGRAPH_SHADOWCTL as the comparison function.
	{
		float sc[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
		for (int i = 0; i < 4; i++) {
			uint32_t texCtl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
			if (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) {
				uint32_t texFmt = pg->regs[RI(NV_PGRAPH_TEXFMT0 + i * 4)];
				xbox::X_D3DFORMAT xboxFmt = GetXboxPixelContainerFormat(texFmt);
				if (EmuXBFormatIsDepthBuffer(xboxFmt)) {
					sc[i] = 1.0f;
				}
			}
		}
		aux.ShadowCompare = { sc[0], sc[1], sc[2], sc[3] };
	}

	// Upload aux cbuffer and bind to b0 (bind only once — buffer pointer is stable)
	CxbxD3D11UpdateDynamicBuffer(g_pD3D11RCInterpreterAuxCB, &aux, sizeof(aux));
	{
		static bool s_AuxCBBound = false;
		if (!s_AuxCBBound) {
			g_pD3DDeviceContext->PSSetConstantBuffers(CXBX_D3D11_PS_CB_SLOT, 1, &g_pD3D11RCInterpreterAuxCB);
			s_AuxCBBound = true;
		}
	}

	// Store a copy for the PS JIT to use for hashing
	g_LastPSAuxCB = aux;
}

void CxbxUpdateActivePixelShader() // NOPATCH
{
  // Always use the RC interpreter ubershader — PGRAPH combiners are authoritative.
  // Even when COMBINECTL == 0 (no combiner stages), the RC interpreter handles
  // this correctly as a passthrough (final combiner only).

  if (!g_pD3D11RCInterpreterPS) {
	if (!CxbxD3D11InitRCInterpreter()) {
		EmuLog(LOG_LEVEL::ERROR2, "RC Interpreter init failed");
		return;
	}
  }

  // Upload combiner state (aux CB + regs SRV). The JIT needs g_LastPSAuxCB for
  // key building, and the interpreter needs both the aux CB and regs SRV.
  CxbxD3D11UploadRCInterpreterState();

  // Try JIT-compiled pixel shader first
  try {
      ID3D11PixelShader* pJIT = g_PixelShaderCache.GetShader(g_pD3DDevice);
      if (pJIT) {
          InterlockedIncrement(&g_ProfilePSJITHits);
          CxbxSetPixelShader(pJIT);
          return;
      }
      // JIT returned nullptr — will use interpreter fallback
  } catch (const std::exception& e) {
      static int s_ExcCount = 0;
      if (s_ExcCount++ < 5)
          EmuLog(LOG_LEVEL::WARNING, "PS JIT exception: %s", e.what());
  } catch (...) {
      static int s_ExcCount2 = 0;
      if (s_ExcCount2++ < 5)
          EmuLog(LOG_LEVEL::WARNING, "PS JIT unknown exception");
  }

  // Fall back to the interpreter ubershader
  InterlockedIncrement(&g_ProfilePSInterpreterHits);
  CxbxSetPixelShader(g_pD3D11RCInterpreterPS);
}
