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
#include "core\hle\D3D8\Rendering\Shaders\Shader.h"
#include "core\hle\D3D8\XbPixelShader.h"
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbD3D8Logging.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\kernel\init\CxbxKrnl.h"
#include "core\hle\D3D8\Rendering\Shaders\CxbxFixedFunctionPixelShader.hlsli"
#include "common/FilePaths.hpp"
#include "devices\Xbox.h"              // For extern NV2ADevice* g_NV2A
#include "devices\video\nv2a.h"        // For NV2ADevice::GetDeviceState(), NV2AState, PGRAPHState, nv2a_regs.h
#include <assert.h>
#include <process.h>
#include <locale.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "Rendering\RenderStates.h"
#include "Rendering\TextureStates.h"
#include <wrl/client.h>
#include <cstring> // For std::memcpy
#include "Rendering\Backend\Backend_D3D11.h"
#include "Rendering\Backend\Backend_D3D11_Internal.h"

// The Xbox kernel's SetRenderState_FogColor swaps R↔B before calling
// SetRenderState_Simple, so D3D__RenderState stores the fog color in NV2A
// ABGR format (0xAABBGGRR) rather than D3DCOLOR ARGB (0xAARRGGBB).
// This helper reverses that swap when reading the stored value.
static inline DWORD FogColor_ABGR_to_ARGB(DWORD color)
{
	return (color & 0xFF00FF00) | ((color & 0x00FF0000) >> 16) | ((color & 0x000000FF) << 16);
}

std::string_view GetD3DTOPString(int d3dtop) {
	static constexpr std::string_view opToString[] = {
		"X_D3DTOP_DISABLE", // 0 (initialized for disabled stages)
		"X_D3DTOP_DISABLE", // 1
		"X_D3DTOP_SELECTARG1", // 2
		"X_D3DTOP_SELECTARG2", // 3
		"X_D3DTOP_MODULATE", // 4
		"X_D3DTOP_MODULATE2X", // 5
		"X_D3DTOP_MODULATE4X", // 6
		"X_D3DTOP_ADD", // 7
		"X_D3DTOP_ADDSIGNED", // 8
		"X_D3DTOP_ADDSIGNED2X", // 9
		"X_D3DTOP_SUBTRACT", // 10
		"X_D3DTOP_ADDSMOOTH", // 11
		"X_D3DTOP_BLENDDIFFUSEALPHA", // 12
		"X_D3DTOP_BLENDCURRENTALPHA", // 13
		"X_D3DTOP_BLENDTEXTUREALPHA", // 14
		"X_D3DTOP_BLENDFACTORALPHA", // 15
		"X_D3DTOP_BLENDTEXTUREALPHAPM", // 16
		"X_D3DTOP_PREMODULATE", // 17
		"X_D3DTOP_MODULATEALPHA_ADDCOLOR", // 18
		"X_D3DTOP_MODULATECOLOR_ADDALPHA", // 19
		"X_D3DTOP_MODULATEINVALPHA_ADDCOLOR", // 20
		"X_D3DTOP_MODULATEINVCOLOR_ADDALPHA", // 21
		"X_D3DTOP_DOTPRODUCT3", // 22
		"X_D3DTOP_MULTIPLYADD", // 23
		"X_D3DTOP_LERP", // 24
		"X_D3DTOP_BUMPENVMAP", // 25
		"X_D3DTOP_BUMPENVMAPLUMINANCE", // 26
	};

	if (d3dtop < 0 || d3dtop > 26) {
		EmuLog(LOG_LEVEL::ERROR2, "Unmapped texture operation %d", d3dtop);
		d3dtop = 0; // undefined
	}

	return opToString[d3dtop];
}

// Get a string equivalent of '<Texture Argument> + <Modifier>'
std::string GetD3DTASumString(int d3dta, bool allowModifier = true) {
	using namespace FixedFunctionPixelShader;

	static const std::string argToString[] = {
		"X_D3DTA_DIFFUSE", // 0
		"X_D3DTA_CURRENT", // 1
		"X_D3DTA_TEXTURE", // 2
		"X_D3DTA_TFACTOR", // 3
		"X_D3DTA_SPECULAR", // 4
		"X_D3DTA_TEMP", // 5
		"X_D3DTA_CONSTANT", // 6
		"UNDEFINED", // 7
	};

	// Write a texture argument
	const int flagMask = 0x30;
	int iFlags = d3dta & flagMask;
	int i = d3dta & ~flagMask;

	if (i < 0 || i > 6) {
		EmuLog(LOG_LEVEL::ERROR2, "Unmapped texture argument %d on texture arg", i);
		i = 7; // undefined
	}

	auto str = argToString[i];
	if (iFlags) {
		if (!allowModifier) {
			EmuLog(LOG_LEVEL::ERROR2, "Modifier not expected on texture argument");
		}

		if (iFlags == X_D3DTA_COMPLEMENT)
			str += " + X_D3DTA_COMPLEMENT";
		else if (iFlags == X_D3DTA_ALPHAREPLICATE)
			str += " + X_D3DTA_ALPHAREPLICATE";
		else {
			EmuLog(LOG_LEVEL::ERROR2, "Unmapped texture modifier %d", iFlags);
			str += " /* + UNKNOWN MODIFIER */";
		}
	}

	return str;
}

// TODO we have to create and cache shaders over and over and over and over
// Deduplicate this resource management
ID3D11PixelShader* GetFixedFunctionShader()
{
	using namespace FixedFunctionPixelShader;

	// TODO move this cache elsewhere - and flush it when the device is released!
	static std::unordered_map<uint64_t, ID3D11PixelShader*> ffPsCache = {};

	// Support hotloading hlsl
	static int pixelShaderVersion = -1;
	int shaderVersion = g_ShaderSources.Update();
	if (pixelShaderVersion != shaderVersion) {
		pixelShaderVersion = shaderVersion;
		CxbxRawSetPixelShader(nullptr);

		for (auto& hostShader : ffPsCache) {
			if (hostShader.second)
				hostShader.second->Release();
		}

		ffPsCache.clear();
	}

	// Create a key from state that will be baked in to the shader
	PsTextureHardcodedState states[4] = {};
	int sampleType[4] = { SAMPLE_NONE, SAMPLE_NONE, SAMPLE_NONE, SAMPLE_NONE };
	bool pointSpriteEnable = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_POINTSPRITEENABLE);

	bool previousStageDisabled = false;
	for (int i = 0; i < 4; i++) {
		// Determine COLOROP
		// This controls both the texture operation for the colour of the stage
		// and when to stop processing
		// Under certain circumstances we force it to be DISABLE
		auto colorOp = XboxTextureStates.Get(i, xbox::X_D3DTSS_COLOROP);

		// Usually we execute stages up to the first disabled stage
		// However, if point sprites are enabled, we just execute stage 3
		bool forceDisable =
			(!pointSpriteEnable && previousStageDisabled) ||
			(pointSpriteEnable && i < 3);

		// When a texture stage has D3DTSS_COLORARG1 equal to D3DTA_TEXTURE
		// and the texture is not enabled, this stage and all stages after
		// it are not processed.
		// Test cases: Red Dead Revolver, JSRF
		// https://docs.microsoft.com/en-us/windows/win32/direct3d9/texture-blending
		// Don't follow the D3D9 docs if SELECTARG2 is in use (PC D3D9 behaviour, nvidia quirk?)
		// Test case: Crash Nitro Kart (engine speed UI)
		// Use PGRAPH TEXCTL0 enable bit when available to avoid racing
		// g_pXbox_SetTexture[] which is written by the game thread.
		bool texturePresent = false;
		{
			auto pg_ff = &(g_NV2A->GetDeviceState()->pgraph);
			uint32_t texCtl = pg_ff->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
			texturePresent = (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) != 0;
		}
		if (!texturePresent
			&& (XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORARG1) & 0x7) == X_D3DTA_TEXTURE
			&& colorOp != xbox::X_D3DTOP_SELECTARG2)
		{
			forceDisable = true;
		}

		// Set the final COLOROP value
		states[i].COLOROP = forceDisable ? X_D3DTOP_DISABLE : colorOp;

		// If the stage is disabled we don't want its configuration to affect the key
		// Move on to the next stage
		if (colorOp == X_D3DTOP_DISABLE) {
			previousStageDisabled = true;
			continue;
		}

		// Get sample type from PGRAPH TEXFMT0 when available (avoids racing
		// g_pXbox_SetTexture[] which is written by the game thread).
		{
			auto pg_ff = &(g_NV2A->GetDeviceState()->pgraph);
			uint32_t texCtl = pg_ff->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
			if (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) {
				uint32_t texFmt = pg_ff->regs[RI(NV_PGRAPH_TEXFMT0 + i * 4)];
				if (texFmt & NV_PGRAPH_TEXFMT0_CUBEMAPENABLE)
					sampleType[i] = SAMPLE_CUBE;
				else if (((texFmt & NV_PGRAPH_TEXFMT0_DIMENSIONALITY) >> 6) > 2)
					sampleType[i] = SAMPLE_3D;
				else
					sampleType[i] = SAMPLE_2D;
			}
		}

		states[i].COLORARG0 = XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORARG0);
		states[i].COLORARG1 = XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORARG1);
		states[i].COLORARG2 = XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORARG2);

		auto alphaOp = XboxTextureStates.Get(i, xbox::X_D3DTSS_ALPHAOP);
		if (alphaOp == X_D3DTOP_DISABLE) LOG_TEST_CASE("Alpha stage disabled when colour stage is enabled");

		states[i].ALPHAOP = alphaOp;
		states[i].ALPHAARG0 = XboxTextureStates.Get(i, xbox::X_D3DTSS_ALPHAARG0);
		states[i].ALPHAARG1 = XboxTextureStates.Get(i, xbox::X_D3DTSS_ALPHAARG1);
		states[i].ALPHAARG2 = XboxTextureStates.Get(i, xbox::X_D3DTSS_ALPHAARG2);

		states[i].RESULTARG = XboxTextureStates.Get(i, xbox::X_D3DTSS_RESULTARG);
	}

	// Create a key from the shader state
	// Note currently this is padded since it's what we send to the GPU
	auto key = 3 * ComputeHash(states, sizeof(states))
		+ ComputeHash(sampleType, sizeof(sampleType));

	auto got = ffPsCache.find(key);
	if (got != ffPsCache.end()) {
		// We have a shader. Great!
		return got->second;
	}

	// Build and compile a new shader
	std::string hlslTemplate = g_ShaderSources.fixedFunctionPixelShaderHlsl;

	// In D3D9 it seems we need to know hardcode if we're doing a 2D or 3D lookup
	const std::string sampleTypePattern = "TEXTURE_SAMPLE_TYPE;";
	auto sampleTypeReplace = hlslTemplate.find(sampleTypePattern);
	std::string finalShader = hlslTemplate;

	if (sampleTypeReplace != std::string::npos) {
		static constexpr std::string_view typeToString[] = {
			"SAMPLE_NONE",
			"SAMPLE_2D",
			"SAMPLE_3D",
			"SAMPLE_CUBE"
		};

		std::stringstream sampleTypeString;
		sampleTypeString << "{"
			<< typeToString[sampleType[0]] << ", "
			<< typeToString[sampleType[1]] << ", "
			<< typeToString[sampleType[2]] << ", "
			<< typeToString[sampleType[3]] << "};";

		finalShader = hlslTemplate.replace(sampleTypeReplace, sampleTypePattern.size(), sampleTypeString.str());
	}

	// Hardcode the texture stage operations and arguments
	// So the shader handles exactly one combination of values
	const std::string stageDef = "// STAGE DEFINITIONS";
	auto stageDefInsert = finalShader.find(stageDef);
	if (stageDefInsert != std::string::npos) {
		stageDefInsert += stageDef.size();

		std::stringstream stageSetup;
		stageSetup << '\n';

		for (int i = 0; i < 4; i++) {
			// Even when a stage is disabled, we still have to fully initialize it's values, to prevent
			// "error X4000: variable 'stages' used without having been completely initialized"
			std::string target = "stages[" + std::to_string(i) + "].";

			auto s = states[i];
			stageSetup << target << "COLOROP = " << GetD3DTOPString(static_cast<int>(s.COLOROP)) << ";\n";

			stageSetup << target << "COLORARG0 = " << GetD3DTASumString(static_cast<int>(s.COLORARG0)) << ";\n";
			stageSetup << target << "COLORARG1 = " << GetD3DTASumString(static_cast<int>(s.COLORARG1)) << ";\n";
			stageSetup << target << "COLORARG2 = " << GetD3DTASumString(static_cast<int>(s.COLORARG2)) << ";\n";

			stageSetup << target << "ALPHAOP = " << GetD3DTOPString(static_cast<int>(s.ALPHAOP)) << ";\n";

			stageSetup << target << "ALPHAARG0 = " << GetD3DTASumString(static_cast<int>(s.ALPHAARG0)) << ";\n";
			stageSetup << target << "ALPHAARG1 = " << GetD3DTASumString(static_cast<int>(s.ALPHAARG1)) << ";\n";
			stageSetup << target << "ALPHAARG2 = " << GetD3DTASumString(static_cast<int>(s.ALPHAARG2)) << ";\n";

			stageSetup << target << "RESULTARG = " << GetD3DTASumString(static_cast<int>(s.RESULTARG), false) << ";\n";
			stageSetup << '\n';
		}

		finalShader = finalShader.insert(stageDefInsert, stageSetup.str());
	}

	// Compile the shader
	ID3DBlob* pShaderBlob = nullptr;

	auto hlslDir = std::filesystem::path(szFilePath_CxbxReloaded_Exe)
		.parent_path()
		.append("hlsl");

	auto pseudoFileName = "FixedFunctionPixelShader-" + std::to_string(key) + ".hlsl";
	auto pseudoSourceFile = hlslDir.append(pseudoFileName).string();
	EmuCompileShader(finalShader, "ps_5_0", &pShaderBlob, pseudoSourceFile.c_str(),
		/*asyncAllowed=*/false, /*useSharedCache=*/true);

	ID3D11PixelShader* pShader = nullptr;
	if (pShaderBlob) {
		// Create shader object for the device
		auto hRet = CxbxCreatePixelShader(pShaderBlob->GetBufferPointer(), pShaderBlob->GetBufferSize(), &pShader);
		if (hRet != S_OK) {
			EmuLog(LOG_LEVEL::ERROR2, "Failed to compile fixed function pixel shader");
		}
		pShaderBlob->Release();
	}

	// Insert the shader into the cache
	ffPsCache[key] = pShader;

	return pShader;
};

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
	using namespace FixedFunctionPixelShader;

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
	// Initially use what the running executable put in COLORSIGN :
	DWORD XboxColorSign = XboxTextureStates.Get(stage_nr, xbox::X_D3DTSS_COLORSIGN);

	{ // This mimicks behaviour of XDK LazySetShaderStageProgram, which we bypass due to our drawing patches without trampolines.
		// When bump environment mapping is enabled
		if (XboxTextureStates.Get(stage_nr, xbox::X_D3DTSS_COLOROP) >= xbox::X_D3DTOP_BUMPENVMAP)
			// Always mark the blue (alias for U) and green (alias for  V) color channels as signed :
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

// Set constant state for the fixed function pixel shader
void UpdateFixedFunctionPixelShaderState()
{
	using namespace FixedFunctionPixelShader;

	FixedFunctionPixelShaderState ffPsState;
	{ D3DXCOLOR c(XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_TEXTUREFACTOR)); ffPsState.TextureFactor = D3DXVECTOR4(c.r, c.g, c.b, c.a); }
	ffPsState.SpecularEnable = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_SPECULARENABLE) ? 1 : 0;
	ffPsState.FogEnable = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_FOGENABLE) ? 1 : 0;
	{
		DWORD raw = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_FOGCOLOR);
		DWORD converted = FogColor_ABGR_to_ARGB(raw);
		D3DXCOLOR c(converted);
		ffPsState.FogColor = D3DXVECTOR3(c.r, c.g, c.b);
	}
	ffPsState.FogTableMode = XboxRenderStates.GetXboxRenderState(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGTABLEMODE);
	ffPsState.FogDensity = XboxRenderStates.GetXboxRenderStateAsFloat(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGDENSITY);
	ffPsState.FogStart = XboxRenderStates.GetXboxRenderStateAsFloat(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGSTART);
	ffPsState.FogEnd = XboxRenderStates.GetXboxRenderStateAsFloat(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGEND);
	// Alpha test state (D3D11 has no fixed-function alpha test)
	ffPsState.AlphaTest.x = static_cast<float>(XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_ALPHATESTENABLE));
	ffPsState.AlphaTest.y = static_cast<float>(XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_ALPHAREF)) / 255.0f;
	ffPsState.AlphaTest.z = static_cast<float>(EmuXB2PC_D3DCMPFUNC((xbox::X_D3DCMPFUNC)XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_ALPHAFUNC)));
	ffPsState.AlphaTest.w = 0;
	// Texture state
	for (int i = 0; i < xbox::X_D3DTS_STAGECOUNT; i++) {
		auto stage = &ffPsState.stages[i];
		stage->COLORKEYOP = XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORKEYOP);
		auto CxbxColorSign = CxbxCalcColorSign(i);
		stage->COLORSIGN.x = CxbxColorSign.r;
		stage->COLORSIGN.y = CxbxColorSign.g;
		stage->COLORSIGN.z = CxbxColorSign.b;
		stage->COLORSIGN.w = CxbxColorSign.a;
		stage->ALPHAKILL = XboxTextureStates.Get(i, xbox::X_D3DTSS_ALPHAKILL);
		stage->BUMPENVMAT00 = AsFloat(XboxTextureStates.Get(i, xbox::X_D3DTSS_BUMPENVMAT00));
		stage->BUMPENVMAT01 = AsFloat(XboxTextureStates.Get(i, xbox::X_D3DTSS_BUMPENVMAT01));
		stage->BUMPENVMAT10 = AsFloat(XboxTextureStates.Get(i, xbox::X_D3DTSS_BUMPENVMAT10));
		stage->BUMPENVMAT11 = AsFloat(XboxTextureStates.Get(i, xbox::X_D3DTSS_BUMPENVMAT11));
		stage->BUMPENVLSCALE = AsFloat(XboxTextureStates.Get(i, xbox::X_D3DTSS_BUMPENVLSCALE));
		stage->BUMPENVLOFFSET = AsFloat(XboxTextureStates.Get(i, xbox::X_D3DTSS_BUMPENVLOFFSET));
		{ D3DXCOLOR c(XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORKEYCOLOR)); stage->COLORKEYCOLOR = D3DXVECTOR4(c.r, c.g, c.b, c.a); }
		stage->TEXFMTFIXUP = CxbxGetTexFmtFixup(i);
	}

	const int size = (sizeof(FixedFunctionPixelShaderState) + 16 - 1) / 16;
	CxbxSetPixelShaderConstantF(0, (float*)&ffPsState, size);
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
	CxbxD3D11UpdateDynamicBuffer(g_pD3D11PGRegsBuf, pg->regs, sizeof(pg->regs));
	// Bind the regs SRV to PS t12
	g_pD3DDeviceContext->PSSetShaderResources(CXBX_D3D11_PS_PGREGS_SRV_SLOT, 1, &g_pD3D11PGRegsSRV);

	// --- Build the auxiliary cbuffer (software-computed fields only) ---
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

	// Color key per stage
	for (int i = 0; i < 4; i++) {
		aux.ColorKeyOp[i] = { static_cast<float>(XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORKEYOP)), 0.0f, 0.0f, 0.0f };
		D3DXCOLOR ckc(XboxTextureStates.Get(i, xbox::X_D3DTSS_COLORKEYCOLOR));
		aux.ColorKeyColor[i] = { ckc.r, ckc.g, ckc.b, ckc.a };
	}

	// Alpha kill per stage (D3DTALPHAKILL_ENABLE = 4)
	aux.AlphaKill = {
		static_cast<float>(XboxTextureStates.Get(0, xbox::X_D3DTSS_ALPHAKILL) & 4 ? 1 : 0),
		static_cast<float>(XboxTextureStates.Get(1, xbox::X_D3DTSS_ALPHAKILL) & 4 ? 1 : 0),
		static_cast<float>(XboxTextureStates.Get(2, xbox::X_D3DTSS_ALPHAKILL) & 4 ? 1 : 0),
		static_cast<float>(XboxTextureStates.Get(3, xbox::X_D3DTSS_ALPHAKILL) & 4 ? 1 : 0)
	};

	// Fog info: x=tableMode (from PGRAPH FOG_MODE), y/z/w unused by RC interpreter.
	// FogColor is read directly from g_PGRegs[] in the shader.
	{
		unsigned int fogMode = GET_MASK(pg->regs[RI(NV_PGRAPH_CONTROL_3)], NV_PGRAPH_CONTROL_3_FOG_MODE);
		aux.FogInfo = { static_cast<float>(fogMode), 0.0f, 0.0f, 0.0f };
		aux.FogEnable.value = (pg->regs[RI(NV_PGRAPH_CONTROL_3)] & NV_PGRAPH_CONTROL_3_FOGENABLE) ? 1u : 0u;
	}

	// Front-face factor for two-sided lighting
	{
		float ff = 0.0f;
		if (XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_TWOSIDEDLIGHTING)) {
			bool cwFrontface = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_FRONTFACE) == 0x900;
			ff = cwFrontface ? 1.0f : -1.0f;
		}
		aux.FrontFaceInfo = { ff, 0.0f, 0.0f, 0.0f };
	}

	// Upload aux cbuffer and bind to b0
	CxbxD3D11UpdateDynamicBuffer(g_pD3D11RCInterpreterAuxCB, &aux, sizeof(aux));
	g_pD3DDeviceContext->PSSetConstantBuffers(CXBX_D3D11_PS_CB_SLOT, 1, &g_pD3D11RCInterpreterAuxCB);
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

  // Bind the ubershader
  CxbxSetPixelShader(g_pD3D11RCInterpreterPS);

  // Upload combiner state as cbuffer
  CxbxD3D11UploadRCInterpreterState();
  g_bRCInterpreterCBActive = true;
}
