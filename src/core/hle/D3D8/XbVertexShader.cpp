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
#include "core\hle\D3D8\Rendering\Shaders\VertexShaderCache.h" // For g_VertexShaderCache
#include "core\hle\D3D8\Rendering\Shaders\Shader.h" // For g_ShaderSources
#include "core\hle\D3D8\XbVertexBuffer.h" // For CxbxImpl_SetVertexData4f
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbPushBuffer.h" // For g_NV2A, HLE_get_NV2A_vertex_constant_float4_ptr
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
                bool g_bRecordingPushBuffer = false;

// Variable set by [D3DDevice|CxbxImpl]_LoadVertexShader() / [D3DDevice|CxbxImpl]_LoadVertexShaderProgram() (both through CxbxCopyVertexShaderFunctionSlots):
                xbox::dword_xt g_Xbox_VertexShader_FunctionSlots[(X_VSH_MAX_INSTRUCTION_COUNT + 1) * X_VSH_INSTRUCTION_SIZE] = { 0 }; // One extra for FLD_FINAL terminator


static xbox::X_D3DVertexShader g_Xbox_VertexShader_ForFVF = {};

static uint32_t                g_X_VERTEXSHADER_FLAG_PROGRAM; // X_VERTEXSHADER_FLAG_PROGRAM flag varies per XDK, so it is set on runtime
static uint32_t                g_X_VERTEXSHADER_FLAG_VALID_MASK; // For a test case

// Track the current active vertex shader key (used to retrieve bytecode for input layout creation)
static ShaderKey g_D3D11ActiveVertexShaderKey = 0;
static bool g_D3D11HasActiveShaderKey = false;
// Retained bytecode for FixedFunction and Passthrough vertex shaders (needed for input layout creation)
static ID3DBlob* g_pD3D11FixedFunctionBytecode = nullptr;
static ID3DBlob* g_pD3D11PassthroughBytecode = nullptr;

void CxbxVertexShaderSetFlags()
{
	// Set an appropriate X_VERTEXSHADER_FLAG_PROGRAM version and mask off the "wrong" one
	// to allow the test case to spot it
	if (g_LibVersion_D3D8 <= 3948) {
		g_X_VERTEXSHADER_FLAG_PROGRAM = X_VERTEXSHADER3948_FLAG_PROGRAM;
		g_X_VERTEXSHADER_FLAG_VALID_MASK = ~X_VERTEXSHADER_FLAG_PROGRAM;
	}
	else {
		g_X_VERTEXSHADER_FLAG_PROGRAM = X_VERTEXSHADER_FLAG_PROGRAM;
		g_X_VERTEXSHADER_FLAG_VALID_MASK = ~X_VERTEXSHADER3948_FLAG_PROGRAM;
	}
}

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
			LOG_TEST_CASE("Unassigned Xbox vertex shader!");
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

extern bool g_InlineVertexBuffer_DeclarationOverride; // TMP glue
extern xbox::X_VERTEXATTRIBUTEFORMAT g_InlineVertexBuffer_AttributeFormat; // TMP glue

xbox::X_VERTEXATTRIBUTEFORMAT* GetXboxVertexAttributeFormat()
{
	// Special case for CxbxImpl_End() based drawing
	if (g_InlineVertexBuffer_DeclarationOverride) {
		return &g_InlineVertexBuffer_AttributeFormat;
	}

	xbox::X_D3DVertexShader* pXboxVertexShader = GetXboxVertexShader();
	if (pXboxVertexShader == xbox::zeroptr) {
		// Despite possibly not being used, the pXboxVertexShader argument must always be assigned
		LOG_TEST_CASE("Xbox should always have a VertexShader set (even for FVF's)");
		return &g_Xbox_SetVertexShaderInput_Attributes; // WRONG result, but it's already strange this happens
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

xbox::dword_xt* GetCxbxVertexShaderSlotPtr(const DWORD SlotIndexAddress)
{
	if (SlotIndexAddress < X_VSH_MAX_INSTRUCTION_COUNT) {
		return &g_Xbox_VertexShader_FunctionSlots[SlotIndexAddress * X_VSH_INSTRUCTION_SIZE];
	} else {
		LOG_TEST_CASE("SlotIndexAddress out of range"); // FIXME : extend with value (once supported by LOG_TEST_CASE)
		return nullptr;
	}
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

extern ID3D11VertexShader* CxbxCreateVertexShader(ID3DBlob* pCompiledShader, const char *shader_category); // Implemented in VertexShaderCache.cpp

ID3D11VertexShader* InitShader(void (*compileFunc)(ID3DBlob**), const char* label, ID3DBlob** ppRetainedBytecode = nullptr) {
	ID3D11VertexShader* shader = nullptr;

	ID3DBlob* pBlob = nullptr;
	compileFunc(&pBlob);
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
	PGRAPHState *pg = nullptr;
	if (g_NV2A) {
		NV2AState *nv2a = g_NV2A->GetDeviceState();
		pg = &nv2a->pgraph;
	}

	if (pg) {
		CxbxD3D11UpdateDynamicBuffer(g_pD3D11XFPRBuf,
			pg->program_data, sizeof(pg->program_data));
	}

	// Bind the shared PGRAPH regs SRV to VS t12 (same buffer, different stage)
	g_pD3DDeviceContext->VSSetShaderResources(CXBX_D3D11_VS_PGREGS_SRV_SLOT, 1, &g_pD3D11PGRegsSRV);

	// Bind the XFPR SRV to VS t5
	g_pD3DDeviceContext->VSSetShaderResources(CXBX_D3D11_VS_XFPR_SRV_SLOT, 1, &g_pD3D11XFPRSRV);
}

void CxbxUpdateHostVertexShader()
{
	extern bool g_bUsePassthroughHLSL; // TMP glue
	// TODO: move render state to VertexShader.cpp
	static ID3D11VertexShader* fixedFunctionShader = nullptr; // TODO: move to shader cache
	static ID3D11VertexShader* passthroughShader = nullptr;
	static int vertexShaderVersion = -1;

	int shaderVersion = g_ShaderSources.Update();
	if (vertexShaderVersion != shaderVersion) {
		vertexShaderVersion = shaderVersion;
		CxbxSetVertexShader(nullptr);

		EmuLog(LOG_LEVEL::INFO, "Loading vertex shaders...");

		g_VertexShaderCache.Clear();

		if (fixedFunctionShader) {
			fixedFunctionShader->Release();
			fixedFunctionShader = nullptr;
		}
		if (g_pD3D11FixedFunctionBytecode) { g_pD3D11FixedFunctionBytecode->Release(); g_pD3D11FixedFunctionBytecode = nullptr; }
		fixedFunctionShader = InitShader(EmuCompileFixedFunction, "Fixed Function Vertex Shader", &g_pD3D11FixedFunctionBytecode);

		if (passthroughShader) {
			passthroughShader->Release();
			passthroughShader = nullptr;
		}
		if (g_pD3D11PassthroughBytecode) { g_pD3D11PassthroughBytecode->Release(); g_pD3D11PassthroughBytecode = nullptr; }
		passthroughShader = InitShader(EmuCompileXboxPassthrough, "Passthrough Vertex Shader", &g_pD3D11PassthroughBytecode);

		// Invalidate the VS interpreter so it recompiles from updated sources
		if (g_pD3D11VSInterpreterVS) { g_pD3D11VSInterpreterVS->Release(); g_pD3D11VSInterpreterVS = nullptr; }
		if (g_pD3D11VSInterpreterBytecode) { g_pD3D11VSInterpreterBytecode->Release(); g_pD3D11VSInterpreterBytecode = nullptr; }
		if (g_pD3D11XFPRSRV) { g_pD3D11XFPRSRV->Release(); g_pD3D11XFPRSRV = nullptr; }
		if (g_pD3D11XFPRBuf) { g_pD3D11XFPRBuf->Release(); g_pD3D11XFPRBuf = nullptr; }
	}

	// TODO Call this when state is dirty
	// Rather than every time state changes

	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction) {
		HRESULT hRet = CxbxSetVertexShader(fixedFunctionShader);
		if (FAILED(hRet)) CxbxrAbort("Failed to set fixed-function shader");
		g_D3D11HasActiveShaderKey = false; // Prevent stale programmable shader key from being used for input layout
	}
	else if (g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough && g_bUsePassthroughHLSL
		&& false // D3D11: passthrough goes through the shader cache (NV2A binary → template)
	) {
		HRESULT hRet = CxbxSetVertexShader(passthroughShader);
		if (FAILED(hRet)) CxbxrAbort("Failed to set passthrough shader");
		g_D3D11HasActiveShaderKey = false; // Prevent stale programmable shader key from being used for input layout
	}
	else {
		// Read program tokens from PGRAPH program_data (the authoritative
		// source, written by the PFIFO puller from NV2A_VP_UPLOAD_INST).
		// The start address comes from CSV0_C CHEOPS_PROGRAM_START, which
		// the puller sets from NV097_SET_TRANSFORM_PROGRAM_START.
		xbox::dword_xt *pTokens = nullptr;
		if (g_NV2A) {
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
			g_D3D11HasActiveShaderKey = false;
			DEBUG_D3DRESULT(hRet, "CxbxSetVertexShader(VSInterpreter)");
		} else {
			// Fallback: compile pipeline
			DWORD shaderSize;
			auto VertexShaderKey = g_VertexShaderCache.CreateShader(pTokens, &shaderSize);
			ID3D11VertexShader* pHostVertexShader = g_VertexShaderCache.GetShader(VertexShaderKey);
			// Track the active shader key so CxbxUpdateHostVertexDeclaration can create the input layout
			g_D3D11ActiveVertexShaderKey = VertexShaderKey;
			g_D3D11HasActiveShaderKey = true;
			HRESULT hRet = CxbxSetVertexShader(pHostVertexShader);
			DEBUG_D3DRESULT(hRet, "CxbxSetVertexShader(pHostVertexShader)");
		}
	}
}

CxbxVertexDeclaration* CxbxGetVertexDeclaration()
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	xbox::X_VERTEXATTRIBUTEFORMAT *pXboxVertexAttributeFormat = GetXboxVertexAttributeFormat();

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
	if (g_D3D11HasActiveShaderKey)
		return g_VertexShaderCache.GetShaderBytecode(g_D3D11ActiveVertexShaderKey);
	// VS interpreter provides its own bytecode for input layout creation
	if (g_bUseVSInterpreter && g_pD3D11VSInterpreterBytecode &&
		(g_Xbox_VertexShaderMode == VertexShaderMode::ShaderProgram ||
		 g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough))
		return g_pD3D11VSInterpreterBytecode;
	if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction)
		return g_pD3D11FixedFunctionBytecode;
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
	CxbxVertexDeclaration* pCxbxVertexDeclaration = CxbxGetVertexDeclaration();
	CxbxSetHostVertexDeclaration(pCxbxVertexDeclaration);

	// Titles can specify default values for registers via calls like SetVertexData4f
	// HLSL shaders need to know whether to use vertex data or default vertex shader values
	// Any register not in the vertex declaration should be set to the default value
	float vertexDefaultFlags[X_VSH_MAX_ATTRIBUTES];
	for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
		vertexDefaultFlags[i] = pCxbxVertexDeclaration->vRegisterInDeclaration[i] ? 0.0f : 1.0f;
	}
	CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_BASE, vertexDefaultFlags, CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_SIZE);
}

// Note : SetVertexShaderInputDirect needs no EMUPATCH CxbxImpl_..., since it just calls SetVertexShaderInput

void CxbxImpl_SetVertexShaderInput(DWORD Handle, UINT StreamCount, xbox::X_STREAMINPUT* pStreamInputs)
{
	using namespace xbox;

	// If Handle is NULL, all VertexShader input state is cleared.
	// Otherwise, Handle is the address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// (Thus, a FVF handle is an invalid argument.)

	if (Handle == NULL)
	{
		// Xbox doesn't remember a null-handle - this may be an XDK bug!
		// (Although, if that's skipped intentionally, we'd need to be very carefull about that!)
		// StreamCount and pStreamInputs arguments are ignored
		g_Xbox_SetVertexShaderInput_Count = 0;
	}
	else
	{
		assert(VshHandleIsVertexShader(Handle));
		assert(StreamCount > 0);
		assert(StreamCount <= X_VSH_MAX_STREAMS);
		assert(pStreamInputs != xbox::zeroptr);

		X_D3DVertexShader* pXboxVertexShader = VshHandleToXboxVertexShader(Handle);
		assert(pXboxVertexShader);

		// Xbox DOES store the Handle, but since it merely returns this through (unpatched) D3DDevice_GetVertexShaderInput, we don't have to.

		g_Xbox_SetVertexShaderInput_Count = StreamCount; // This > 0 indicates g_Xbox_SetVertexShaderInput_Data has to be used
		memcpy(g_Xbox_SetVertexShaderInput_Data, pStreamInputs, StreamCount * sizeof(xbox::X_STREAMINPUT)); // Make a copy of the supplied StreamInputs array

		g_Xbox_SetVertexShaderInput_Attributes = *CxbxGetVertexShaderAttributes(pXboxVertexShader); // Copy this vertex shaders's attribute slots
	}
}

// Note : SelectVertexShaderDirect needs no EMUPATCH CxbxImpl_..., since it just calls SelectVertexShader

void CxbxImpl_SelectVertexShader(DWORD Handle, DWORD Address)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	// Address always indicates a previously loaded vertex shader slot (from where the program is used).
	// Handle can be null if the current Xbox VertexShader is assigned
	// Handle can be an address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// If Handle is assigned, it becomes the new current Xbox VertexShader,
	// which resets a bit of state (nv2a execution mode, viewport, ?)
	// Either way, the given address slot is selected as the start of the current vertex shader program
	// g_Xbox_VertexShader_FunctionSlots_StartAddress is no longer written here;
	// the render thread reads the program start from NV_PGRAPH_CSV0_C
	// CHEOPS_PROGRAM_START, which the PFIFO puller sets from
	// NV097_SET_TRANSFORM_PROGRAM_START in the push buffer.

	// NOTE: Do NOT mirror the start address to pg->regs[CSV0_C] here.
	// The Xbox trampoline writes SET_TRANSFORM_PROGRAM_START to the push
	// buffer, and the PFIFO puller sets CSV0_C sequentially before the draw.
	// Writing from the game thread races with the puller reading CSV0_C at
	// draw time, causing the wrong VS program to be selected.

	// g_Xbox_VertexShaderMode is derived from PGRAPH CSV0_D by the
	// render thread in CxbxUpdateNativeD3DResources — don't race it.

	if (Handle) {
		if (!VshHandleIsVertexShader(Handle))
			LOG_TEST_CASE("Non-zero handle must be a VertexShader!");

		g_Xbox_VertexShader_Handle = Handle;
	}
}

// Set default values for attributes missing from vertex declaration
void SetFixedFunctionDefaultVertexAttributes(DWORD vshFlags) {
	// Test case: Mechassault (skybox)
	// Test case: KOTOR (overlay)
	auto decl = CxbxGetVertexDeclaration();
	for (int i = 0; i < xbox::X_D3DVSDE_TEXCOORD3; i++) {
		if (decl->vRegisterInDeclaration[i]) {
			continue; // only reset missing attributes
		}

		const float white[4] = { 1, 1, 1, 1 };
		const float black[4] = { 0, 0, 0, 0 };
		const float unset[4] = { 0, 0, 0, 1 };
		const float* value = unset;

		// Account for flags that override this reset behaviour
		if (i == xbox::X_D3DVSDE_DIFFUSE && !(vshFlags & X_VERTEXSHADER_FLAG_HASDIFFUSE) ||
			i == xbox::X_D3DVSDE_BACKDIFFUSE && !(vshFlags & X_VERTEXSHADER_FLAG_HASBACKDIFFUSE)) {
			value = white;
		}
		else if (i == xbox::X_D3DVSDE_SPECULAR && !(vshFlags & X_VERTEXSHADER_FLAG_HASSPECULAR) ||
			i == xbox::X_D3DVSDE_BACKSPECULAR && !(vshFlags & X_VERTEXSHADER_FLAG_HASBACKSPECULAR)) {
			value = black;
		}

		// Note : We avoid calling CxbxImpl_SetVertexData4f here, as that would
		// start populating g_InlineVertexBuffer_Table, which is not our intent here.
		CxbxSetVertexAttribute(i, value[0], value[1], value[2], value[3]);
	}
}

void CxbxImpl_SetVertexShader(DWORD Handle)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	CxbxD3D11IABypassInvalidateLayout();

	// Checks if the Handle has bit 0 set - if not, it's a FVF
	// which is converted to a global Xbox Vertex Shader struct
	// Otherwise bit 0 is cleared and the resulting address is
	// validated to be a valid Xbox Vertex Shader
	// D3D state fields are updated.
	// If the shader contains a program, the handle is passed to
	// D3DDevice_LoadVertexShader and D3DDevice_SelectVertexShader.
	// Otherwise the shader is send using push buffer commands.

	HRESULT hRet = S_OK;

	xbox::X_D3DVertexShader* pXboxVertexShader = CxbxGetXboxVertexShaderForHandle(Handle);

	if ((pXboxVertexShader->Flags & g_X_VERTEXSHADER_FLAG_VALID_MASK) != pXboxVertexShader->Flags) {
		LOG_TEST_CASE("Unknown vertex shader flag");
	}

	if (pXboxVertexShader->Flags & g_X_VERTEXSHADER_FLAG_PROGRAM) { // Global variable set from CxbxVertexShaderSetFlags
		// The Xbox SetVertexShader trampoline (called above) internally calls
		// LoadVertexShader and SelectVertexShader, which are intercepted by
		// our EMUPATCH stubs. Those call XB_TRMP(LoadVertexShader) and
		// XB_TRMP(SelectVertexShader), writing SET_TRANSFORM_PROGRAM and
		// SET_TRANSFORM_PROGRAM_START to the push buffer.  The PFIFO puller
		// processes these sequentially into PGRAPH before each draw.
		//
		// Do NOT call CxbxImpl_Load/Select here — they write directly to
		// pg->program_data and pg->regs[CSV0_C] on the game thread, racing
		// with the puller which reads that data at draw time.
		// Update only the HLE-side globals that the EMUPATCH stubs set:
		// g_Xbox_VertexShader_FunctionSlots_StartAddress is read from
		// PGRAPH CSV0_C by the render thread — don't race it.
		// g_Xbox_VertexShaderMode is derived from PGRAPH CSV0_D by the
		// render thread in CxbxUpdateNativeD3DResources — don't race it.
		g_Xbox_VertexShader_Handle = Handle;
	} else {
		// A shader without a program won't call LoadVertexShader nor SelectVertexShader
		g_Xbox_VertexShader_Handle = Handle;
		// g_Xbox_VertexShader_FunctionSlots_StartAddress is read from
		// PGRAPH CSV0_C by the render thread — don't race it.

		// NOTE: Do NOT write to pg->regs[CSV0_C] here.  The Xbox
		// SetVertexShader trampoline writes SET_TRANSFORM_EXECUTION_MODE
		// (FIXED) and SET_TRANSFORM_PROGRAM_START to the push buffer.
		// The puller processes these sequentially before the draw.

		SetFixedFunctionDefaultVertexAttributes(pXboxVertexShader->Flags);

		// Passthrough and fixed-function programs are pushed by the Xbox
		// SetVertexShader trampoline through the push buffer → PFIFO →
		// pg->program_data[].  No HLE-side upload needed.

		// g_Xbox_VertexShaderMode is derived from PGRAPH CSV0_D + VPSCL/VPOFF
		// by the render thread — don't write it here on the game thread.
	}
}

void CxbxImpl_DeleteVertexShader(DWORD Handle)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	// Handle is always address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// It's reference count is lowered. If it reaches zero (0), the struct is freed.

	xbox::X_D3DVertexShader* pXboxVertexShader = VshHandleToXboxVertexShader(Handle);
	if (pXboxVertexShader == nullptr) {
		return;
	}

	if (pXboxVertexShader->RefCount > 1) {
		return;
	}

	// TODO : Decide and implement what parts to free
	// RegisterCxbxVertexDeclaration(pCxbxVertexDeclaration->Key, nullptr);
	// g_VertexShaderCache.ReleaseShader(pCxbxVertexShader->Key);
}

// TODO : Remove SetVertexShaderConstant implementation and the patch once
// CxbxUpdateHostVertexShaderConstants is reliable (ie. : when we're able to flush the NV2A push buffer)
void CxbxImpl_SetVertexShaderConstant(INT Register, PVOID pConstantData, DWORD ConstantCount)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	// Xbox vertex shader constants range from -96 to 95
	// The host does not support negative, so we adjust to 0..191
	Register += X_D3DSCM_CORRECTION;

	if (Register < 0) LOG_TEST_CASE("Register < 0");
	if (Register + ConstantCount > X_D3DVS_CONSTREG_COUNT) LOG_TEST_CASE("Register + ConstantCount > X_D3DVS_CONSTREG_COUNT");

	// Write Vertex Shader constants in nv2a
	float* constant_floats = HLE_get_NV2A_vertex_constant_float4_ptr(Register);
	memcpy(constant_floats, pConstantData, ConstantCount * sizeof(float) * 4);

	// Mark the constant as dirty, so that CxbxUpdateHostVertexShaderConstants will pick it up
	auto nv2a = g_NV2A->GetDeviceState();
	for (DWORD i = 0; i < ConstantCount; i++) {
		nv2a->pgraph.vsh_constants_dirty[Register + i] = true;
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

	Nv2aVshProgram program = {}; // Note: This nulls program.steps
	// TODO : Retain program globally and perform nv2a_vsh_parse_program only when
	//        the address-range we're about to emulate was modified since last parse.
	// TODO : As a suggestion for this, parse all NV2A_MAX_TRANSFORM_PROGRAM_LENGTH slots,
	//        and here just point program.steps to global vsh_program_steps[Address].
	Nv2aVshParseResult result = nv2a_vsh_parse_program(
		&program, // Note : program.steps will be malloc'ed
		GetCxbxVertexShaderSlotPtr(Address), // TODO : At some point, use pg->program_data[Address] here instead
		NV2A_MAX_TRANSFORM_PROGRAM_LENGTH - Address);
	if (result != NV2AVPR_SUCCESS) {
		LOG_TEST_CASE("nv2a_vsh_parse_program failed");
		// TODO : Dump Nv2aVshParseResult as string and program for debugging purposes
		return;
	}

	Nv2aVshCPUXVSSExecutionState state_linkage;
	Nv2aVshExecutionState state = nv2a_vsh_emu_initialize_xss_execution_state(
		&state_linkage, (float*)pg->vsh_constants); // Note : This wil memset(state_linkage, 0)
	if (pData != nullptr)
		//if pData != nullptr, then it contains v0.xyzw, we shall copy the binary content directly.
		memcpy(state_linkage.input_regs, pData, sizeof(state_linkage.input_regs));

	nv2a_vsh_emu_execute_track_context_writes(&state, &program, pg->vsh_constants_dirty);
	// Note: Above emulation's primary purpose is to update pg->vsh_constants and pg->vsh_constants_dirty
	// therefor, nothing else needs to be done here, other than to cleanup

	nv2a_vsh_program_destroy(&program); // Note: program.steps will be free'ed
}
