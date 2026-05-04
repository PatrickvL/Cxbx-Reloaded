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
// *  (c) 2002-2004 Aaron Robinson <caustik@caustik.com>
// *                Kingofc <kingofc@freenet.de>
// *
// *  All rights reserved
// *
// ******************************************************************
#define LOG_PREFIX CXBXR_MODULE::VTXSH

//#define _DEBUG_TRACK_VS

#include "common\util\hasher.h" // For ComputeHash
#include "core\kernel\init\CxbxKrnl.h"
#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\Rendering\RenderGlobals.h" // For g_Xbox_VertexShader_Handle
#include "core\hle\D3D8\Rendering\RenderStates.h" // For XboxRenderStateConverter
#include "core\hle\D3D8\Rendering\Shaders\Shader.h" // For LoadPrecompiledCSO

#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbPushBuffer.h" // For g_NV2A
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11.h"
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11_Internal.h" // For g_pD3D11XFPRBuf, g_pD3D11PGRegsSRV
#include "core\hle\D3D8\XbD3D8Logging.h" // For DEBUG_D3DRESULT
#include "devices\xbox.h"
#include "core\hle\D3D8\XbConvert.h" // For NV2A_VP_UPLOAD_INST, NV2A_VP_UPLOAD_CONST_ID, NV2A_VP_UPLOAD_CONST
#include "devices\video\nv2a.h" // For D3DPUSH_DECODE
#include "common\Logging.h" // For LOG_INIT
#include "common\Settings.hpp" // for g_LibVersion_D3D8

#include "XbD3D8Types.h" // For X_D3DVSDE_*
#include <sstream>
#include <unordered_map>

#include <array>
#include <bitset>
#include <filesystem>

#include "nv2a_vsh_emulator.h"

// Variables set by [D3DDevice|CxbxImpl]_SetVertexShaderInput() :
                      unsigned g_Xbox_SetVertexShaderInput_Count = 0; // Read by GetXboxVertexAttributes
           xbox::X_STREAMINPUT g_Xbox_SetVertexShaderInput_Data[X_VSH_MAX_STREAMS] = { 0 }; // Active when g_Xbox_SetVertexShaderInput_Count > 0
 xbox::X_VERTEXATTRIBUTEFORMAT g_Xbox_SetVertexShaderInput_Attributes = { 0 }; // Read by GetXboxVertexAttributes when g_Xbox_SetVertexShaderInput_Count > 0

VertexShaderMode g_Xbox_VertexShaderMode = VertexShaderMode::FixedFunction;

                xbox::dword_xt g_Xbox_VertexShader_Handle = 0;

static xbox::X_D3DVertexShader g_Xbox_VertexShader_ForFVF = {};

// Retained bytecode for FixedFunction and Passthrough vertex shaders (needed for input layout creation)
static ID3DBlob* g_pD3D11FixedFunctionBytecode = nullptr;
static ID3DBlob* g_pD3D11PassthroughBytecode = nullptr;

extern bool g_bUsePassthroughHLSL; // defined in HostDevice.cpp

// Converts an Xbox FVF shader handle to X_D3DVertexShader
// Note : Temporary, until we reliably locate the Xbox internal state for this
// See D3DXDeclaratorFromFVF docs https://docs.microsoft.com/en-us/windows/win32/direct3d9/d3dxdeclaratorfromfvf
// and https://github.com/reactos/wine/blob/2e8dfbb1ad71f24c41e8485a39df01bb9304127f/dlls/d3dx9_36/mesh.c#L2041
static xbox::X_D3DVertexShader* XboxVertexShaderFromFVF(DWORD xboxFvf) // TODO : Rename CxbxFVFToXboxVertexAttributeFormat?
{
	using namespace xbox;

	// Note : FVFs don't tessellate, all slots read from stream zero, therefore
	// the following zero-initialization of StreamIndex (like all other fields)
	// is never updated below.
	g_Xbox_VertexShader_ForFVF = { 0 };

	// Shorthand, glueing earlier implementation below to global g_Xbox_VertexShader_ForFVF variable :
	X_VERTEXATTRIBUTEFORMAT &declaration = g_Xbox_VertexShader_ForFVF.VertexAttribute;

	static DWORD X_D3DVSDT_FLOAT[] = { 0, X_D3DVSDT_FLOAT1, X_D3DVSDT_FLOAT2, X_D3DVSDT_FLOAT3, X_D3DVSDT_FLOAT4 };

	static const DWORD InvalidXboxFVFBits = X_D3DFVF_RESERVED0 | X_D3DFVF_RESERVED1 /* probably D3DFVF_PSIZE if detected */
		| 0x0000F000 // Bits between texture count and the texture formats
		| 0xFF000000; // All bits above the four alllowed texture formats

	if (xboxFvf & InvalidXboxFVFBits) {
		// Test-case : JSRF (after "now Loading...") TODO : Figure out what's going on
		// LOG_TEST_CASE("Invalid Xbox FVF bits detected!");
	}

	// Position & Blendweights
	int nrPositionFloats = 3;
	int nrBlendWeights = 0;
	unsigned offset = 0;
	DWORD position = (xboxFvf & X_D3DFVF_POSITION_MASK);
	switch (position) {
		case 0: nrPositionFloats = 0; LOG_TEST_CASE("FVF without position"); break; // Note : Remove logging if this occurs often
		case X_D3DFVF_XYZ: /*nrPositionFloats is set to 3 by default*/ break;
		case X_D3DFVF_XYZRHW: nrPositionFloats = 4; g_Xbox_VertexShader_ForFVF.Flags |= X_VERTEXSHADER_FLAG_PASSTHROUGH; break;
		case X_D3DFVF_XYZB1: nrBlendWeights = 1; break;
		case X_D3DFVF_XYZB2: nrBlendWeights = 2; break;
		case X_D3DFVF_XYZB3: nrBlendWeights = 3; break;
		case X_D3DFVF_XYZB4: nrBlendWeights = 4; break;
		case X_D3DFVF_POSITION_MASK: /*Keep nrPositionFloats set to 3*/ LOG_TEST_CASE("FVF invalid (5th blendweight?)"); break;
		DEFAULT_UNREACHABLE;
	}

	// Assign vertex element (attribute) slots
	X_VERTEXSHADERINPUT* pSlot;

	// Write Position
	if (nrPositionFloats > 0) {
		pSlot = &declaration.Slots[X_D3DVSDE_POSITION];
		pSlot->Format = X_D3DVSDT_FLOAT[nrPositionFloats];
		pSlot->Offset = offset;
		offset += sizeof(float) * nrPositionFloats;
		// Write Blend Weights
		if (nrBlendWeights > 0) {
			pSlot = &declaration.Slots[X_D3DVSDE_BLENDWEIGHT];
			pSlot->Format = X_D3DVSDT_FLOAT[nrBlendWeights];
			pSlot->Offset = offset;
			offset += sizeof(float) * nrBlendWeights;
		}
	}
	else if (nrBlendWeights > 0) LOG_TEST_CASE("BlendWeights given without position?");

	// Write Normal, Diffuse, and Specular
	if (xboxFvf & X_D3DFVF_NORMAL) {
		if (position == X_D3DFVF_XYZRHW) {
			LOG_TEST_CASE("X_D3DFVF_NORMAL shouldn't use X_D3DFVF_XYZRHW");
		}

		pSlot = &declaration.Slots[X_D3DVSDE_NORMAL];
		pSlot->Format = X_D3DVSDT_FLOAT[3];
		pSlot->Offset = offset;
		offset += sizeof(float) * 3;
	}

	if (xboxFvf & X_D3DFVF_DIFFUSE) {
		g_Xbox_VertexShader_ForFVF.Flags |= X_VERTEXSHADER_FLAG_HASDIFFUSE;
		pSlot = &declaration.Slots[X_D3DVSDE_DIFFUSE];
		pSlot->Format = X_D3DVSDT_D3DCOLOR;
		pSlot->Offset = offset;
		offset += sizeof(DWORD) * 1;
	}

	if (xboxFvf & X_D3DFVF_SPECULAR) {
		g_Xbox_VertexShader_ForFVF.Flags |= X_VERTEXSHADER_FLAG_HASSPECULAR; 
		pSlot = &declaration.Slots[X_D3DVSDE_SPECULAR];
		pSlot->Format = X_D3DVSDT_D3DCOLOR;
		pSlot->Offset = offset;
		offset += sizeof(DWORD) * 1;
	}

	// Write Texture Coordinates
	int textureCount = (xboxFvf & X_D3DFVF_TEXCOUNT_MASK) >> X_D3DFVF_TEXCOUNT_SHIFT;
	if (textureCount > 4) {
		LOG_TEST_CASE("Limiting FVF to 4 textures");
		textureCount = 4; // Safeguard, since the X_D3DFVF_TEXCOUNT bitfield could contain invalid values (5 up to 15)
	}

	for (int i = 0; i < textureCount; i++) {
		auto FVFTextureFormat = (xboxFvf >> X_D3DFVF_TEXCOORDSIZE_SHIFT(i)) & 0x003;
#if 1
		int numberOfCoordinates = ((FVFTextureFormat + 1) & 3) + 1;
#else
		int numberOfCoordinates = 0;
		switch (FVFTextureFormat) { // Note : Below enums are not ordered; In a math expression mapped as :
			case X_D3DFVF_TEXTUREFORMAT1: numberOfCoordinates = 1; break; // input = 3 -> 4 -> 0 -> 1 = output
			case X_D3DFVF_TEXTUREFORMAT2: numberOfCoordinates = 2; break; // input = 0 -> 1 -> 1 -> 2 = output
			case X_D3DFVF_TEXTUREFORMAT3: numberOfCoordinates = 3; break; // input = 1 -> 2 -> 2 -> 3 = output
			case X_D3DFVF_TEXTUREFORMAT4: numberOfCoordinates = 4; break; // input = 2 -> 3 -> 3 -> 4 = output
			DEFAULT_UNREACHABLE;                                          // ((input   +1 ) &3 ) +1 ) = output
		}

		assert(numberOfCoordinates > 0);
#endif
		pSlot = &declaration.Slots[X_D3DVSDE_TEXCOORD0 + i];
		pSlot->Format = X_D3DVSDT_FLOAT[numberOfCoordinates];
		pSlot->Offset = offset;
		offset += sizeof(float) * numberOfCoordinates;
		// Update the VertexShader texture Dimensionality field here as well
		g_Xbox_VertexShader_ForFVF.Dimensionality[i] = numberOfCoordinates;
	}

	// Make sure all unused slots have a X_D3DVSDT_NONE format
	// TODO : Actually, maybe not, since this could avoid VshConvertToken_STREAMDATA_REG() calls!
	for (unsigned i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
		pSlot = &declaration.Slots[i];
		if (pSlot->Format == 0) {
			pSlot->Format = X_D3DVSDT_NONE;
		}
	}

	// Return the global g_Xbox_VertexShader_ForFVF variable 
	return &g_Xbox_VertexShader_ForFVF;
}

static xbox::X_D3DVertexShader* CxbxGetXboxVertexShaderForHandle(DWORD Handle)
{
	if (VshHandleIsVertexShader(Handle)) {
		return VshHandleToXboxVertexShader(Handle);
	} else {
		return XboxVertexShaderFromFVF(Handle);
	}
}

// TODO : Start using this function everywhere g_Xbox_VertexShader_Handle is accessed currently!
xbox::X_D3DVertexShader* GetXboxVertexShader()
{
	// LOG_INIT; // Allows use of DEBUG_D3DRESULT

	using namespace xbox;

	X_D3DVertexShader* pXboxVertexShader = xbox::zeroptr;

		// We use what we've last stored in the g_Xbox_VertexShader_Handle
		// variable via our D3DDevice_SetVertexShader and
		// D3DDevice_SelectVertexShader* patches.

		// Now, to convert, we do need to have a valid vertex shader :
		if (g_Xbox_VertexShader_Handle == 0) {
			// Expected when SetVertexShader patches are disabled — PGRAPH
			// vertex_attributes[] is the authoritative source instead.
			return nullptr;
		}

		pXboxVertexShader = CxbxGetXboxVertexShaderForHandle(g_Xbox_VertexShader_Handle);

	return pXboxVertexShader;
}

static bool UseXboxD3DVertexShaderTypeForVersion3948(const xbox::X_D3DVertexShader* pXboxVertexShader)
{
	// Don't check XDK version for our internal FVF vertex shader
	// because g_Xbox_VertexShader_ForFVF is an internal variable
	// that's compiled in as a xbox::X_D3DVertexShader
	if (pXboxVertexShader == &g_Xbox_VertexShader_ForFVF) {
		return false;
	}

	return g_LibVersion_D3D8 <= 3948;
}

static xbox::X_VERTEXATTRIBUTEFORMAT* CxbxGetVertexShaderAttributes(xbox::X_D3DVertexShader* pXboxVertexShader)
{
	if (UseXboxD3DVertexShaderTypeForVersion3948(pXboxVertexShader)) {
		auto pXboxVertexShader3948 = reinterpret_cast<xbox::X_D3DVertexShader3948*>(pXboxVertexShader);
		return &(pXboxVertexShader3948->VertexAttribute);
	}

	return &(pXboxVertexShader->VertexAttribute);
}

static DWORD* CxbxGetVertexShaderTokens(xbox::X_D3DVertexShader* pXboxVertexShader, DWORD* pNrTokens)
{
	if (UseXboxD3DVertexShaderTypeForVersion3948(pXboxVertexShader)) {
		auto pXboxVertexShader3948 = reinterpret_cast<xbox::X_D3DVertexShader3948*>(pXboxVertexShader);
		*pNrTokens = pXboxVertexShader3948->ProgramAndConstantsDwords;
		return &pXboxVertexShader3948->ProgramAndConstants[0];
	}

	*pNrTokens = pXboxVertexShader->ProgramAndConstantsDwords;
	return &pXboxVertexShader->ProgramAndConstants[0];
}

int GetXboxVertexDataComponentCount(int d3dvsdt) {
	using namespace xbox;
	switch (d3dvsdt) {
	case X_D3DVSDT_NORMPACKED3:
		return 3;
	case X_D3DVSDT_FLOAT2H:
		LOG_TEST_CASE("Attempting to use component count for X_D3DVSDT_FLOAT2H, which uses an odd (value, value, 0, value) layout");
		// This is a bit of an odd case. Will call it 4 since it writes a value to the 4th component...
		return 4;
	default:
		// Most data types have a representation consistent with the number of components
		const int countMask = 0x7;
		const int countShift = 4;
		return (d3dvsdt >> countShift) & countMask;
	}
}

xbox::X_VERTEXATTRIBUTEFORMAT* GetXboxVertexAttributeFormat()
{
	xbox::X_D3DVertexShader* pXboxVertexShader = GetXboxVertexShader();
	if (pXboxVertexShader == xbox::zeroptr) {
		// With SetVertexShader patches disabled, g_Xbox_VertexShader_Handle is
		// never set. Return nullptr so callers can fall back to PGRAPH state.
		return nullptr;
	}

	// If SetVertexShaderInput is active, its arguments overrule those of the active vertex shader
	if (g_Xbox_SetVertexShaderInput_Count > 0) {
		// Take overrides (on declarations and streaminputs, as optionally set by SetVertexShaderInput) into account :
		// Test-case : Crazy taxi 3
		LOG_TEST_CASE("SetVertexShaderInput_Attributes override in effect!");
		return &g_Xbox_SetVertexShaderInput_Attributes;
	}

	return CxbxGetVertexShaderAttributes(pXboxVertexShader);
}

// Reads the active Xbox stream input values (containing VertexBuffer, Offset and Stride) for the given stream index.
// (These values are set through SetStreamSource and can be overridden by SetVertexShaderInput.)
xbox::X_STREAMINPUT& GetXboxVertexStreamInput(unsigned XboxStreamNumber)
{
	// If SetVertexShaderInput is active, its arguments overrule those of SetStreamSource
	if (g_Xbox_SetVertexShaderInput_Count > 0) {
		return g_Xbox_SetVertexShaderInput_Data[XboxStreamNumber];
	}

	return g_Xbox_SetStreamSource[XboxStreamNumber];
}

#define DbgVshPrintf \
	LOG_CHECK_ENABLED(LOG_LEVEL::DEBUG) \
		if(g_bPrintfOn) printf


// Defined in XbVertexShaderDecoder.cpp
extern D3D11_INPUT_ELEMENT_DESC *EmuRecompileVshDeclaration(
	xbox::X_VERTEXATTRIBUTEFORMAT* pXboxDeclaration,
	bool bIsFixedFunction,
	CxbxVertexDeclaration *pCxbxVertexDeclaration
);

static bool FreeCxbxVertexDeclaration(CxbxVertexDeclaration *pCxbxVertexDeclaration)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	if (pCxbxVertexDeclaration) {
		if (pCxbxVertexDeclaration->pHostVertexDeclaration) {
			HRESULT hRet = pCxbxVertexDeclaration->pHostVertexDeclaration->Release();
			DEBUG_D3DRESULT(hRet, "pHostVertexDeclaration->Release()");
		}
		if (pCxbxVertexDeclaration->pD3D11InputElements) {
			free(pCxbxVertexDeclaration->pD3D11InputElements);
			pCxbxVertexDeclaration->pD3D11InputElements = nullptr;
		}
		free(pCxbxVertexDeclaration);
		return true;
	}

	return false;
}

VertexDeclarationKey GetXboxVertexAttributesKey(xbox::X_VERTEXATTRIBUTEFORMAT* pXboxVertexAttributeFormat)
{
	auto attributeHash = ComputeHash((void*)pXboxVertexAttributeFormat, sizeof(xbox::X_VERTEXATTRIBUTEFORMAT));
	// For now, we use different declarations depending on if the fixed function pipeline
	// is in use, even if the attributes are the same
	return g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction
		? attributeHash
		: attributeHash ^ 1;
}

std::unordered_map<VertexDeclarationKey, CxbxVertexDeclaration*> g_CxbxVertexDeclarations;

void RegisterCxbxVertexDeclaration(VertexDeclarationKey CacheKey, CxbxVertexDeclaration* pCxbxVertexDeclaration)
{
	auto it = g_CxbxVertexDeclarations.find(CacheKey);
	if (it != g_CxbxVertexDeclarations.end() && it->second != nullptr) {
		LOG_TEST_CASE("Overwriting existing Vertex Declaration");
		FreeCxbxVertexDeclaration(it->second); // Avoid memory leak
	}

	g_CxbxVertexDeclarations[CacheKey] = pCxbxVertexDeclaration;
}

CxbxVertexDeclaration* FetchCachedCxbxVertexDeclaration(VertexDeclarationKey CacheKey)
{
	auto it = g_CxbxVertexDeclarations.find(CacheKey);
	if (it != g_CxbxVertexDeclarations.end()) {
		return it->second;
	}

	return nullptr;
}

extern ID3D11VertexShader* CxbxCreateVertexShader(ID3DBlob* pCompiledShader, const char *shader_category)
{
	ID3D11VertexShader* pHostVertexShader = nullptr;

	if (g_pD3DDevice == nullptr) {
		EmuLog(LOG_LEVEL::WARNING, "Can't create %s vertex shader - no D3D device is set!", shader_category);
	}
	else {
		assert(pCompiledShader);

		HRESULT hRet;
		hRet = g_pD3DDevice->CreateVertexShader(
			(const void*)pCompiledShader->GetBufferPointer(),
			pCompiledShader->GetBufferSize(),
			nullptr,
			&pHostVertexShader
		);
		if (FAILED(hRet)) CxbxrAbort("Failed to create %s vertex shader", shader_category);
	}

	return pHostVertexShader;
}

ID3D11VertexShader* InitShader(const char* csoName, const char* label, ID3DBlob** ppRetainedBytecode = nullptr) {
	ID3D11VertexShader* shader = nullptr;

	ID3DBlob* pBlob = nullptr;
	LoadPrecompiledCSO(csoName, &pBlob);
	if (pBlob) {
		shader = CxbxCreateVertexShader(pBlob, label);
		if (ppRetainedBytecode) {
			*ppRetainedBytecode = pBlob; // Caller takes ownership
		} else {
			pBlob->Release();
		}
	}

	return shader;
}

// Upload NV2A XFPR (Transform Program RAM) and bind SRVs for the VS interpreter.
//
// Two StructuredBuffers feed the interpreter shader:
//   - g_PGRegs (t12): shared PGRAPH register array — already uploaded by
//     CxbxD3D11UploadRCInterpreterState(). We just bind it to the VS stage.
//   - g_XFPR (t5): pg->program_data[136][4] — the XFPR RAM mirror,
//     uploaded here.  On real NV2A hardware this is on-chip XF SRAM
//     behind the RDI interface, uploaded via NV097_SET_TRANSFORM_PROGRAM
//     with an auto-incrementing write pointer (CHEOPS_OFFSET.PROG_LD_PTR).
//
// The shader reads CHEOPS_PROGRAM_START from g_PGRegs to find the first
// active instruction slot and loops until FLD_FINAL.
void CxbxD3D11UploadVSInterpreterState(const xbox::dword_xt* /*pXboxMicrocode*/)
{
	if (!g_pD3D11XFPRBuf || !g_pD3D11PGRegsSRV)
		return;

	// PGRAPH source: upload the entire program_data[] array (XFPR mirror).
	// The shader selects the active program via CHEOPS_PROGRAM_START.
	PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;

	// Skip XFPR upload if program data hasn't changed (dirty flag set by NV097_SET_TRANSFORM_PROGRAM)
	static bool s_XFPRUploaded = false;
	if (!s_XFPRUploaded || pg->program_data_dirty) {
		CxbxD3D11UpdateDynamicBuffer(g_pD3D11XFPRBuf,
			pg->program_data, sizeof(pg->program_data));
		s_XFPRUploaded = true;
		// Note: program_data_dirty is cleared by RunVertexStateShader cache logic
	}

	// Bind VS interpreter SRVs once — pointers are stable for device lifetime
	static bool s_VSInterpreterSRVsBound = false;
	if (!s_VSInterpreterSRVsBound) {
		g_pD3DDeviceContext->VSSetShaderResources(CXBX_D3D11_VS_PGREGS_SRV_SLOT, 1, &g_pD3D11PGRegsSRV);
		g_pD3DDeviceContext->VSSetShaderResources(CXBX_D3D11_VS_XFPR_SRV_SLOT, 1, &g_pD3D11XFPRSRV);
		s_VSInterpreterSRVsBound = true;
	}
}

void CxbxUpdateHostVertexShader()
{
	// Vertex shaders are loaded once from embedded precompiled blobs (CSOs).
	// They persist for the lifetime of the D3D11 device; teardown is handled
	// by CxbxD3D11ReleaseBackendResources() on device release.
	static ID3D11VertexShader* fixedFunctionShader = nullptr;
	static ID3D11VertexShader* passthroughShader = nullptr;
	static bool shadersLoaded = false;

	if (!shadersLoaded) {
		shadersLoaded = true;
		CxbxSetVertexShader(nullptr);

		EmuLog(LOG_LEVEL::INFO, "Loading vertex shaders...");
		fixedFunctionShader = InitShader("CxbxFixedFunctionVS", "Fixed Function Vertex Shader", &g_pD3D11FixedFunctionBytecode);
		passthroughShader = InitShader("CxbxVSPassthroughVS", "Passthrough Vertex Shader", &g_pD3D11PassthroughBytecode);
		// VS interpreter is initialized lazily on first ShaderProgram draw via CxbxD3D11InitVSInterpreter()
	}

	// Select the active vertex shader based on current PGRAPH mode.
	// Called every draw; the actual shader objects are already loaded above.

	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction) {
		HRESULT hRet = CxbxSetVertexShader(fixedFunctionShader);
		if (FAILED(hRet)) CxbxrAbort("Failed to set fixed-function shader");
	}
	else if (g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough && g_bUsePassthroughHLSL) {
		HRESULT hRet = CxbxSetVertexShader(passthroughShader);
		if (FAILED(hRet)) CxbxrAbort("Failed to set passthrough shader");
	}
	else {
		// Read program tokens from PGRAPH program_data (the authoritative
		// source, written by the PFIFO puller from NV2A_VP_UPLOAD_INST).
		// The start address comes from CSV0_C CHEOPS_PROGRAM_START, which
		// the puller sets from NV097_SET_TRANSFORM_PROGRAM_START.
		xbox::dword_xt *pTokens = nullptr;
		{
			PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;
			uint32_t startAddr = GET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)],
				NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START);
			if (startAddr < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH) {
				pTokens = (xbox::dword_xt*)&pg->program_data[startAddr][0];
			}
		}
		if (!pTokens) {
			LOG_TEST_CASE("PGRAPH program_data not available");
			return;
		}

		if (g_bUseVSInterpreter && CxbxD3D11InitVSInterpreter()) {
			// Upload the raw NV2A microcode to the interpreter constant buffer
			CxbxD3D11UploadVSInterpreterState(pTokens);
			HRESULT hRet = CxbxSetVertexShader(g_pD3D11VSInterpreterVS);
			DEBUG_D3DRESULT(hRet, "CxbxSetVertexShader(VSInterpreter)");
		}
	}
}

CxbxVertexDeclaration* CxbxGetVertexDeclaration()
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	xbox::X_VERTEXATTRIBUTEFORMAT *pXboxVertexAttributeFormat = GetXboxVertexAttributeFormat();
	if (pXboxVertexAttributeFormat == nullptr) {
		// With SetVertexShader patches disabled, HLE attribute format is
		// unavailable. The vertex pull draw path reads PGRAPH directly instead.
		return nullptr;
	}

	auto XboxVertexAttributesKey = GetXboxVertexAttributesKey(pXboxVertexAttributeFormat);
	CxbxVertexDeclaration* pCxbxVertexDeclaration = FetchCachedCxbxVertexDeclaration(XboxVertexAttributesKey);
	if (pCxbxVertexDeclaration == nullptr) {
		pCxbxVertexDeclaration = (CxbxVertexDeclaration*)calloc(1, sizeof(CxbxVertexDeclaration));
		// calloc zero-initializes, but tessellation registers use -1 as "not present"
		pCxbxVertexDeclaration->autoNormalRegister = -1;
		pCxbxVertexDeclaration->autoNormalSourceRegister = -1;
		pCxbxVertexDeclaration->autoTexcoordRegister = -1;

		// Convert Xbox vertex attributes towards host Direct3D vertex element
		D3D11_INPUT_ELEMENT_DESC* pRecompiledVertexElements = EmuRecompileVshDeclaration(
			pXboxVertexAttributeFormat,
			g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction,
			pCxbxVertexDeclaration);

		// Create the vertex declaration
		pCxbxVertexDeclaration->pHostVertexDeclaration = CxbxCreateHostVertexDeclaration(pRecompiledVertexElements);

		// For D3D11, store a copy of the vertex elements for lazy input layout creation
		// Count the elements (terminated by an element with SemanticName==nullptr in D3D11)
		UINT elementCount = 0;
		if (pRecompiledVertexElements != nullptr) {
			while (pRecompiledVertexElements[elementCount].SemanticName != nullptr) {
				elementCount++;
			}
		}
		if (elementCount > 0) {
			pCxbxVertexDeclaration->pD3D11InputElements = (D3D11_INPUT_ELEMENT_DESC*)malloc(elementCount * sizeof(D3D11_INPUT_ELEMENT_DESC));
			memcpy(pCxbxVertexDeclaration->pD3D11InputElements, pRecompiledVertexElements, elementCount * sizeof(D3D11_INPUT_ELEMENT_DESC));
		}
		pCxbxVertexDeclaration->D3D11InputElementCount = elementCount;

		free(pRecompiledVertexElements);

		// Cache resulting declarations from given inputs
		pCxbxVertexDeclaration->Key = XboxVertexAttributesKey;
		RegisterCxbxVertexDeclaration(XboxVertexAttributesKey, pCxbxVertexDeclaration);
	}

	return pCxbxVertexDeclaration;
}

ID3DBlob* CxbxGetActiveVertexShaderBytecode()
{
	if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction)
		return g_pD3D11FixedFunctionBytecode;
	// Return passthrough bytecode when the passthrough HLSL shader is active
	if (g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough && g_bUsePassthroughHLSL)
		return g_pD3D11PassthroughBytecode;
	// VS interpreter provides its own bytecode for input layout creation
	if (g_bUseVSInterpreter && g_pD3D11VSInterpreterBytecode &&
		(g_Xbox_VertexShaderMode == VertexShaderMode::ShaderProgram ||
		 g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough))
		return g_pD3D11VSInterpreterBytecode;
	if (g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough)
		return g_pD3D11PassthroughBytecode;
	return nullptr;
}

ID3DBlob* CxbxGetFixedFunctionVertexShaderBytecode()
{
	return g_pD3D11FixedFunctionBytecode;
}

void CxbxUpdateHostVertexDeclaration()
{
	// Titles can specify default values for registers via calls like SetVertexData4f
	// HLSL shaders need to know whether to use vertex data or default vertex shader values
	// Any register not in the vertex declaration should be set to the default value
	float vertexDefaultFlags[X_VSH_MAX_ATTRIBUTES];

	// When g_Xbox_VertexShader_Handle is set, the HLE vertex declaration
	// tells us which registers are present in the vertex buffer.
	// When it's not set (patches disabled — Xbox native code drives NV2A),
	// derive the flags from PGRAPH vertex_attributes[]: count > 0 means
	// the attribute has stream data, count == 0 means use sticky default.
	if (g_Xbox_VertexShader_Handle != 0) {
		CxbxVertexDeclaration* pCxbxVertexDeclaration = CxbxGetVertexDeclaration();
		CxbxSetHostVertexDeclaration(pCxbxVertexDeclaration);

		for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
			vertexDefaultFlags[i] = pCxbxVertexDeclaration->vRegisterInDeclaration[i] ? 0.0f : 1.0f;
		}
	} else {
		// PGRAPH-driven: read which attributes are active from NV2A state
		CxbxSetHostVertexDeclaration(nullptr);
		PGRAPHState* pg = (g_NV2A != nullptr) ? &g_NV2A->GetDeviceState()->pgraph : nullptr;
		for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
			bool active = pg && (pg->vertex_attributes[i].count > 0);
			vertexDefaultFlags[i] = active ? 0.0f : 1.0f;
		}
	}

	// Only upload if the flags changed since last draw
	static float s_CachedVertexDefaultFlags[X_VSH_MAX_ATTRIBUTES] = {};
	static bool s_FirstDefaultFlagsCall = true;
	if (s_FirstDefaultFlagsCall || std::memcmp(vertexDefaultFlags, s_CachedVertexDefaultFlags, sizeof(vertexDefaultFlags)) != 0) {
		std::memcpy(s_CachedVertexDefaultFlags, vertexDefaultFlags, sizeof(vertexDefaultFlags));
		s_FirstDefaultFlagsCall = false;
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_BASE, vertexDefaultFlags, CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_SIZE);
	}
}

void CxbxrImpl_RunVertexStateShader(DWORD Address, CONST FLOAT *pData)
{
	// If pData is assigned, pData[0..3] is pushed towards nv2a transform data registers
	// then sends the nv2a a command to launch the vertex shader function located at Address
	if (Address >= NV2A_MAX_TRANSFORM_PROGRAM_LENGTH) {
		LOG_TEST_CASE("Address out of bounds");
		return;
	}

	NV2AState* dev = g_NV2A->GetDeviceState();
	PGRAPHState* pg = &(dev->pgraph);

	// Cache the parsed program globally; only re-parse when program_data changes
	static Nv2aVshProgram s_CachedProgram = {};
	static bool s_CachedProgramValid = false;

	if (pg->program_data_dirty || !s_CachedProgramValid) {
		if (s_CachedProgramValid) {
			nv2a_vsh_program_destroy(&s_CachedProgram);
		}
		s_CachedProgram = {};
		Nv2aVshParseResult result = nv2a_vsh_parse_program(
			&s_CachedProgram,
			pg->program_data[0],
			NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
		if (result != NV2AVPR_SUCCESS) {
			LOG_TEST_CASE("nv2a_vsh_parse_program failed (cached full parse)");
			s_CachedProgramValid = false;
			return;
		}
		// Guard against buffer overflow: force is_final on the last slot so
		// the executor always terminates within the 136-entry allocation,
		// even if no instruction in the program sets the final bit.
		s_CachedProgram.steps[NV2A_MAX_TRANSFORM_PROGRAM_LENGTH - 1].is_final = true;
		s_CachedProgramValid = true;
		pg->program_data_dirty = false;
	}

	// Create a view into the cached program starting at Address
	// Execution stops naturally at the step with is_final==true
	Nv2aVshProgram program;
	program.steps = s_CachedProgram.steps + Address;

	Nv2aVshCPUXVSSExecutionState state_linkage;
	Nv2aVshExecutionState state = nv2a_vsh_emu_initialize_xss_execution_state(
		&state_linkage, (float*)pg->vsh_constants); // Note : This wil memset(state_linkage, 0)
	if (pData != nullptr)
		//if pData != nullptr, then it contains v0.xyzw, we shall copy the binary content directly.
		memcpy(state_linkage.input_regs, pData, sizeof(state_linkage.input_regs));

	nv2a_vsh_emu_execute_track_context_writes(&state, &program, pg->vsh_constants_dirty);
	// Note: Above emulation's primary purpose is to update pg->vsh_constants and pg->vsh_constants_dirty
	// Do NOT call nv2a_vsh_program_destroy here — program.steps is a borrowed pointer
}
