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

#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbVertexShader.h"
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

// ****************************************************************************
// * Vertex shader function recompiler
// ****************************************************************************

namespace XboxVertexShaderDecoder
{
	// Xbox Vertex SHader microcode types

	enum VSH_OUTPUT_TYPE {
		OUTPUT_C = 0,
		OUTPUT_O
	};

	enum VSH_OUTPUT_MUX {
		OMUX_MAC = 0,
		OMUX_ILU
	};

	// Host intermediate vertex shader types

	enum VSH_FIELD_NAME {
		FLD_ILU = 0,
		FLD_MAC,
		FLD_CONST,
		FLD_V,
		// Input A
		FLD_A_NEG,
		FLD_A_SWZ_X,
		FLD_A_SWZ_Y,
		FLD_A_SWZ_Z,
		FLD_A_SWZ_W,
		FLD_A_R,
		FLD_A_MUX,
		// Input B
		FLD_B_NEG,
		FLD_B_SWZ_X,
		FLD_B_SWZ_Y,
		FLD_B_SWZ_Z,
		FLD_B_SWZ_W,
		FLD_B_R,
		FLD_B_MUX,
		// Input C
		FLD_C_NEG,
		FLD_C_SWZ_X,
		FLD_C_SWZ_Y,
		FLD_C_SWZ_Z,
		FLD_C_SWZ_W,
		FLD_C_R_HIGH,
		FLD_C_R_LOW,
		FLD_C_MUX,
		// Output
		FLD_OUT_MAC_MASK,
		FLD_OUT_R,
		FLD_OUT_ILU_MASK,
		FLD_OUT_O_MASK,
		FLD_OUT_ORB,
		FLD_OUT_ADDRESS,
		FLD_OUT_MUX,
		// Relative addressing
		FLD_A0X,
		// Final instruction
		FLD_FINAL
	};

	// Retrieves a number of bits in the instruction token
	static inline uint32_t VshGetFromToken(
		uint32_t* pShaderToken,
		uint8_t SubToken,
		uint8_t StartBit,
		uint8_t BitLength)
	{
		return (pShaderToken[SubToken] >> StartBit) & ~(0xFFFFFFFF << BitLength);
	}

	static uint8_t VshGetField(
		uint32_t* pShaderToken,
		VSH_FIELD_NAME FieldName)
	{
		// Used for xvu spec definition
		static const struct {
			uint8_t          SubToken;
			uint8_t          StartBit;
			uint8_t          BitLength;
		} FieldMapping[/*VSH_FIELD_NAME*/] = {
			// SubToken BitPos  BitSize
			{  1,   25,     3 }, // FLD_ILU,              
			{  1,   21,     4 }, // FLD_MAC,              
			{  1,   13,     8 }, // FLD_CONST,            
			{  1,    9,     4 }, // FLD_V,                
			// Input A
			{  1,    8,     1 }, // FLD_A_NEG,            
			{  1,    6,     2 }, // FLD_A_SWZ_X,          
			{  1,    4,     2 }, // FLD_A_SWZ_Y,          
			{  1,    2,     2 }, // FLD_A_SWZ_Z,          
			{  1,    0,     2 }, // FLD_A_SWZ_W,          
			{  2,   28,     4 }, // FLD_A_R,              
			{  2,   26,     2 }, // FLD_A_MUX,            
			// Input B
			{  2,   25,     1 }, // FLD_B_NEG,            
			{  2,   23,     2 }, // FLD_B_SWZ_X,          
			{  2,   21,     2 }, // FLD_B_SWZ_Y,          
			{  2,   19,     2 }, // FLD_B_SWZ_Z,          
			{  2,   17,     2 }, // FLD_B_SWZ_W,          
			{  2,   13,     4 }, // FLD_B_R,              
			{  2,   11,     2 }, // FLD_B_MUX,            
			// Input C
			{  2,   10,     1 }, // FLD_C_NEG,            
			{  2,    8,     2 }, // FLD_C_SWZ_X,          
			{  2,    6,     2 }, // FLD_C_SWZ_Y,          
			{  2,    4,     2 }, // FLD_C_SWZ_Z,          
			{  2,    2,     2 }, // FLD_C_SWZ_W,          
			{  2,    0,     2 }, // FLD_C_R_HIGH,         
			{  3,   30,     2 }, // FLD_C_R_LOW,          
			{  3,   28,     2 }, // FLD_C_MUX,            
			// Output
			{  3,   24,     4 }, // FLD_OUT_MAC_MASK,   
			{  3,   20,     4 }, // FLD_OUT_R,            
			{  3,   16,     4 }, // FLD_OUT_ILU_MASK,
			{  3,   12,     4 }, // FLD_OUT_O_MASK,
			{  3,   11,     1 }, // FLD_OUT_ORB,          
			{  3,    3,     8 }, // FLD_OUT_ADDRESS,      
			{  3,    2,     1 }, // FLD_OUT_MUX,          
			// Relative addressing
			{  3,    1,     1 }, // FLD_A0X,              
			// Final instruction
			{  3,    0,     1 }  // FLD_FINAL,            
		};

		return (uint8_t)(VshGetFromToken(pShaderToken,
			FieldMapping[FieldName].SubToken,
			FieldMapping[FieldName].StartBit,
			FieldMapping[FieldName].BitLength));
	}

	// Converts the C register address to disassembly format
	static inline int16_t ConvertCRegister(const int16_t CReg)
	{
		return ((((CReg >> 5) & 7) - 3) * 32) + (CReg & 31);
	}

	static VSH_IMD_PARAMETER VshGetIntermediateParam(
		uint32_t* pShaderToken,
		VSH_FIELD_NAME FLD_MUX,
		VSH_FIELD_NAME FLD_NEG,
		uint16_t R,
		uint16_t V,
		uint16_t C)
	{
		VSH_IMD_PARAMETER param{};
		param.Type = (VSH_IMD_PARAMETER_TYPE)VshGetField(pShaderToken, FLD_MUX);
		switch (param.Type) {
		case PARAM_R:
			param.Address = R;
			break;
		case PARAM_V:
			param.Address = V;
			break;
		case PARAM_C:
			param.Address = C;
			break;
		default:
			LOG_TEST_CASE("parameter type unknown");
		}

		int d = FLD_NEG - FLD_A_NEG;
		param.Neg = VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_NEG)) > 0;
		param.Swizzle[0] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_X));
		param.Swizzle[1] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_Y));
		param.Swizzle[2] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_Z));
		param.Swizzle[3] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_W));

		return param;
	}

	static bool VshConvertToIntermediate(uint32_t* pShaderToken, IntermediateVertexShader* pShader)
	{
		// First get the instruction(s).
		VSH_ILU ILU = (VSH_ILU)VshGetField(pShaderToken, FLD_ILU);
		VSH_MAC MAC = (VSH_MAC)VshGetField(pShaderToken, FLD_MAC);
		if (MAC > MAC_ARL) LOG_TEST_CASE("Unknown MAC");

		// MAC,ILU output R register
		int16_t RAddress = VshGetField(pShaderToken, FLD_OUT_R);

		// Test for paired opcodes
		bool bIsPaired = (MAC != MAC_NOP) && (ILU != ILU_NOP);

		VSH_IMD_MAC_OP MacOp{};
		VSH_IMD_ILU_OP IluOp{};

		// Set up input registers
		int16_t AR = VshGetField(pShaderToken, FLD_A_R);
		int16_t BR = VshGetField(pShaderToken, FLD_B_R);
		int16_t CR = VshGetField(pShaderToken, FLD_C_R_HIGH) << 2 | VshGetField(pShaderToken, FLD_C_R_LOW);
		int16_t V = VshGetField(pShaderToken, FLD_V);
		int16_t C = ConvertCRegister(VshGetField(pShaderToken, FLD_CONST));

		// Check if there's a MAC opcode
		if (MAC != MAC_NOP && MAC <= MAC_ARL) {
			MacOp.Opcode = MAC;

			if (bIsPaired && RAddress == 1) {
				// Ignore paired MAC opcodes that write to R1
			}
			else if (MAC == MAC_ARL) {
				MacOp.Dest.Type = IMD_DEST_A0X;
				MacOp.Dest.Mask = MASK_X;
			}
			else {
				MacOp.Dest.Type = IMD_DEST_R;
				MacOp.Dest.Address = RAddress;
				MacOp.Dest.Mask = VshGetField(pShaderToken, FLD_OUT_MAC_MASK);
			}

			if (MAC >= MAC_MOV) {
				MacOp.Parameters[MacOp.ParamCount++] = VshGetIntermediateParam(pShaderToken, FLD_A_MUX, FLD_A_NEG, AR, V, C);
			}

			if (MAC == MAC_MUL || (MAC >= MAC_MAD && MAC <= MAC_SGE)) {
				MacOp.Parameters[MacOp.ParamCount++] = VshGetIntermediateParam(pShaderToken, FLD_B_MUX, FLD_B_NEG, BR, V, C);
			}

			if (MAC == MAC_ADD || MAC == MAC_MAD) {
				MacOp.Parameters[MacOp.ParamCount++] = VshGetIntermediateParam(pShaderToken, FLD_C_MUX, FLD_C_NEG, CR, V, C);
			}
		}

		// Check if there's an ILU opcode
		if (ILU != ILU_NOP) {
			// Paired ILU opcodes will only write to R1
			IluOp.Opcode = ILU;
			IluOp.Dest.Type = IMD_DEST_R;
			IluOp.Dest.Address = bIsPaired ? 1 : RAddress;
			IluOp.Dest.Mask = VshGetField(pShaderToken, FLD_OUT_ILU_MASK);
			IluOp.Parameter = VshGetIntermediateParam(pShaderToken, FLD_C_MUX, FLD_C_NEG, CR, V, C);
		}

		// Output register
		VSH_OUTPUT_MUX OutputMux = (VSH_OUTPUT_MUX)VshGetField(pShaderToken, FLD_OUT_MUX);
		int16_t OutputAddress = VshGetField(pShaderToken, FLD_OUT_ADDRESS);
		VSH_IMD_DEST_TYPE OutputType;
		if ((VSH_OUTPUT_TYPE)VshGetField(pShaderToken, FLD_OUT_ORB) == OUTPUT_C) {
			OutputType = IMD_DEST_C;
			OutputAddress = ConvertCRegister(OutputAddress);
		}
		else { // OUTPUT_O:
			OutputType = IMD_DEST_O;
			OutputAddress = OutputAddress & 0xF;
		}

		VSH_IMD_INSTR imd{};
		imd.MAC = MacOp;
		imd.ILU = IluOp;
		imd.IndexesWithA0_X = VshGetField(pShaderToken, FLD_A0X) > 0;
		imd.ORegSource = OutputMux == OMUX_MAC ? SRC_MAC : SRC_ILU;
		imd.ORegDest.Type = OutputType;
		imd.ORegDest.Address = OutputAddress;
		imd.ORegDest.Mask = VshGetField(pShaderToken, FLD_OUT_O_MASK);

		pShader->Instructions.push_back(imd);

		return VshGetField(pShaderToken, FLD_FINAL) == 0;
	}
};

// Get the function size excluding the final field
size_t GetVshFunctionSize(const xbox::dword_xt* pXboxFunction) {
	auto curToken = (uint32_t*)pXboxFunction;

	while (!XboxVertexShaderDecoder::VshGetField(curToken, XboxVertexShaderDecoder::FLD_FINAL)) {
		curToken += X_VSH_INSTRUCTION_SIZE; // TODO use a struct to represent these instructions
	}

	curToken += X_VSH_INSTRUCTION_SIZE; // For the final instruction

	return (curToken - pXboxFunction) * sizeof(xbox::dword_xt);
}

// ****************************************************************************
// * Vertex shader declaration recompiler
// ****************************************************************************

class XboxVertexDeclarationConverter
{
protected:
	// Internal variables
	CxbxVertexDeclaration* pCurrentVertexDeclaration;
	CxbxVertexShaderStreamInfo* pCurrentVertexShaderStreamInfo = nullptr;
	bool IsFixedFunction;
	D3D11_INPUT_ELEMENT_DESC* pCurrentHostVertexElement;
	std::array<bool, 16> RegVIsPresentInDeclaration;

private:

	// VERTEX SHADER

	bool VshConvertToken_STREAMDATA_REG(DWORD VertexRegister, xbox::X_VERTEXSHADERINPUT &slot)
	{
		DWORD XboxVertexElementDataType = slot.Format;

		// Does this attribute use no storage present the vertex (check this as early as possible to avoid needless processing) ?
		if (XboxVertexElementDataType == xbox::X_D3DVSDT_NONE) {
			// Handle tessellating attributes
			// NOTE for D3D11: This code is currently unreachable because the
			// caller's `slot.Format > X_D3DVSDT_NONE` condition filters out
			// these slots before they reach here. The semantics below are kept
			// correct (TEXCOORD + register index) in case the code path is
			// activated in the future, but see the warning at the call site
			// about why it must NOT be activated for D3D11 normal draws.
			switch (slot.TessellationType) {
			case 0: return false; // AUTONONE
			case 1: // AUTONORMAL
				// Note : .Stream, .Offset and .Type are copied from pAttributeSlot->TessellationSource in a post-processing step below,
				// because these could all go through an Xbox to host conversion step, so must be copied over afterwards.
				// D3D11: Use TEXCOORD + register index to match the flat v[16] : TEXCOORD VS_INPUT
				pCurrentHostVertexElement->SemanticName = "TEXCOORD";
				pCurrentHostVertexElement->SemanticIndex = (UINT)VertexRegister;
				pCurrentHostVertexElement->InputSlotClass = D3D11_INPUT_PER_VERTEX_DATA;
				return true;
			case 2: // AUTOTEXCOORD
				// pCurrentHostVertexElement->Stream = 0; // The input stream is unused (but must be set to 0), which is the current default value
				// pCurrentHostVertexElement->AlignedByteOffset = 0; // The input offset is unused (but must be set to 0), which is the current default value
				pCurrentHostVertexElement->Format = DXGI_FORMAT_UNKNOWN; // The input type for D3DDECLMETHOD_UV must be D3DDECLTYPE_UNUSED (the output type implied by D3DDECLMETHOD_UV is D3DDECLTYPE_FLOAT2)
				// D3D11: Use TEXCOORD + register index to match the flat v[16] : TEXCOORD VS_INPUT
				pCurrentHostVertexElement->SemanticName = "TEXCOORD";
				pCurrentHostVertexElement->SemanticIndex = (UINT)VertexRegister;
				pCurrentHostVertexElement->InputSlotClass = D3D11_INPUT_PER_VERTEX_DATA;
				return true;
			default:
				LOG_TEST_CASE("invalid TessellationType");
				return false;
			}
		}

		WORD XboxVertexElementByteSize = 0; // When set above 0, implies NeedPatching
		DXGI_FORMAT HostVertexElementDataType = DXGI_FORMAT_UNKNOWN;
		WORD HostVertexElementByteSize = 0;

		switch (XboxVertexElementDataType) {
		case xbox::X_D3DVSDT_FLOAT1: // 0x12:
			HostVertexElementDataType = DXGI_FORMAT_R32_FLOAT;
			HostVertexElementByteSize = 1 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_FLOAT2: // 0x22:
			HostVertexElementDataType = DXGI_FORMAT_R32G32_FLOAT;
			HostVertexElementByteSize = 2 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_FLOAT3: // 0x32:
			HostVertexElementDataType = DXGI_FORMAT_R32G32B32_FLOAT;
			HostVertexElementByteSize = 3 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_FLOAT4: // 0x42:
			HostVertexElementDataType = DXGI_FORMAT_R32G32B32A32_FLOAT;
			HostVertexElementByteSize = 4 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_D3DCOLOR: // 0x40:
			// D3DCOLOR memory layout is [B,G,R,A]. DXGI_FORMAT_B8G8R8A8_UNORM is not universally
			// supported for vertex input (optional D3D11 feature), so use R8G8B8A8_UNORM and
			// swap R/B channels in the vertex buffer patching step.
			HostVertexElementDataType = DXGI_FORMAT_R8G8B8A8_UNORM;
			HostVertexElementByteSize = 1 * sizeof(D3DCOLOR);
			XboxVertexElementByteSize = 1 * sizeof(D3DCOLOR); // trigger NeedPatch for BGRA→RGBA swizzle
			break;
		case xbox::X_D3DVSDT_SHORT2: // 0x25:
			HostVertexElementDataType = DXGI_FORMAT_R16G16_SINT;
			HostVertexElementByteSize = 2 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_SHORT4: // 0x45:
			HostVertexElementDataType = DXGI_FORMAT_R16G16B16A16_SINT;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_NORMSHORT1: // 0x11:
			// Test-cases : Halo - Combat Evolved
			HostVertexElementDataType = DXGI_FORMAT_R16_SNORM;
			HostVertexElementByteSize = 1 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_NORMSHORT2: // 0x21:
			// Test-cases : Baldur's Gate: Dark Alliance 2, F1 2002, Gun, Halo - Combat Evolved, Scrapland
			HostVertexElementDataType = DXGI_FORMAT_R16G16_SNORM;
			HostVertexElementByteSize = 2 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_NORMSHORT3: // 0x31:
			// Test-cases : Cel Damage, Constantine, Destroy All Humans!
			HostVertexElementDataType = DXGI_FORMAT_R16G16B16A16_SNORM;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			XboxVertexElementByteSize = 3 * sizeof(xbox::short_xt);
			break;
		case xbox::X_D3DVSDT_NORMSHORT4: // 0x41:
			// Test-cases : Judge Dredd: Dredd vs Death, NHL Hitz 2002, Silent Hill 2, Sneakers, Tony Hawk Pro Skater 4
			HostVertexElementDataType = DXGI_FORMAT_R16G16B16A16_SNORM;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_NORMPACKED3: // 0x16:
			// Test-cases : Dashboard
			HostVertexElementDataType = DXGI_FORMAT_R32G32B32_FLOAT; // D3D11 Note : DXGI_FORMAT_R11G11B10_FLOAT can't be used (it contains floats, while X_D3DVSDT_NORMPACKED3 contains ints) and DXGI_FORMAT_R11G11B10_SNORM doesn't exist
			HostVertexElementByteSize = 3 * sizeof(FLOAT);
			XboxVertexElementByteSize = 1 * sizeof(xbox::dword_xt);
			break;
		case xbox::X_D3DVSDT_SHORT1: // 0x15:
			HostVertexElementDataType = DXGI_FORMAT_R16_SINT;
			HostVertexElementByteSize = 1 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_SHORT3: // 0x35:
			// Test-cases : Turok
			HostVertexElementDataType = DXGI_FORMAT_R16G16B16A16_SINT;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			XboxVertexElementByteSize = 3 * sizeof(xbox::short_xt);
			break;
		case xbox::X_D3DVSDT_PBYTE1: // 0x14:
			HostVertexElementDataType = DXGI_FORMAT_R8_UNORM;
			HostVertexElementByteSize = 1 * sizeof(BYTE);
			break;
		case xbox::X_D3DVSDT_PBYTE2: // 0x24:
			HostVertexElementDataType = DXGI_FORMAT_R8G8_UNORM;
			HostVertexElementByteSize = 2 * sizeof(BYTE);
			break;
		case xbox::X_D3DVSDT_PBYTE3: // 0x34:
			// Test-cases : Turok
			HostVertexElementDataType = DXGI_FORMAT_R8G8B8A8_UNORM; // D3D11 note : DXGI_FORMAT_R8G8B8_UNORM doesn't exist, so expand to 4 bytes
			HostVertexElementByteSize = 4 * sizeof(BYTE);
			XboxVertexElementByteSize = 3 * sizeof(xbox::byte_xt);
			break;
		case xbox::X_D3DVSDT_PBYTE4: // 0x44:
			// Test-case : Jet Set Radio Future, Panzer
			HostVertexElementDataType = DXGI_FORMAT_R8G8B8A8_UNORM;
			HostVertexElementByteSize = 4 * sizeof(BYTE);
			break;
		case xbox::X_D3DVSDT_FLOAT2H: // 0x72:
			HostVertexElementDataType = DXGI_FORMAT_R32G32B32A32_FLOAT;
			HostVertexElementByteSize = 4 * sizeof(FLOAT);
			XboxVertexElementByteSize = 3 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_NONE: // 0x02:
			assert(false); // Should already be handled above
			break;
		default:
			LOG_TEST_CASE("Unknown data type for D3DVSD_REG"); // TODO : Add : 0x % 02X\n", XboxVertexElementDataType);
			return false;
		}

		assert(HostVertexElementDataType != DXGI_FORMAT_UNKNOWN);
		assert(HostVertexElementByteSize > 0);

		// Select new stream, if needed
		if ((pCurrentVertexShaderStreamInfo == nullptr)
		 || (pCurrentVertexShaderStreamInfo->XboxStreamIndex != slot.StreamIndex)) {
			assert(slot.StreamIndex < X_VSH_MAX_STREAMS);
			assert(pCurrentVertexDeclaration->NumberOfVertexStreams < X_VSH_MAX_STREAMS);

			pCurrentVertexShaderStreamInfo =
				&(pCurrentVertexDeclaration->VertexStreams[
					pCurrentVertexDeclaration->NumberOfVertexStreams++]);
			pCurrentVertexShaderStreamInfo->NeedPatch = FALSE;
			pCurrentVertexShaderStreamInfo->XboxStreamIndex = (WORD)slot.StreamIndex;
			pCurrentVertexShaderStreamInfo->HostVertexStride = (WORD)slot.Offset;
			pCurrentVertexShaderStreamInfo->NumberOfVertexElements = 0;
			// Dxbx note : Use Dophin(s), FieldRender, MatrixPaletteSkinning and PersistDisplay as a testcase
		}

		// save patching information
		assert(pCurrentVertexShaderStreamInfo->NumberOfVertexElements < X_VSH_MAX_ATTRIBUTES);

		CxbxVertexShaderStreamElement* pCurrentVertexShaderStreamElementInfo =
			&(pCurrentVertexShaderStreamInfo->VertexElements[
				pCurrentVertexShaderStreamInfo->NumberOfVertexElements++]);

		if (XboxVertexElementByteSize == 0) {
			XboxVertexElementByteSize = HostVertexElementByteSize;
		} else {
			pCurrentVertexShaderStreamInfo->NeedPatch |= TRUE;
		}

		pCurrentVertexShaderStreamElementInfo->XboxType = XboxVertexElementDataType;
		pCurrentVertexShaderStreamElementInfo->XboxByteSize = XboxVertexElementByteSize;
		pCurrentVertexShaderStreamElementInfo->HostDataType = HostVertexElementDataType;
		pCurrentVertexShaderStreamElementInfo->HostByteSize = HostVertexElementByteSize;

		// Convert to host vertex element
		pCurrentHostVertexElement->InputSlot = pCurrentVertexShaderStreamInfo->XboxStreamIndex; // Use Xbox stream index on host
		// FIXME Don't assume vertex elements are contiguous!
		pCurrentHostVertexElement->AlignedByteOffset = pCurrentVertexShaderStreamInfo->HostVertexStride;
		pCurrentHostVertexElement->Format = pCurrentVertexShaderStreamElementInfo->HostDataType;
		// All VS types (FixedFunction, ShaderProgram) use a flat
		// float4 v[16] : TEXCOORD array as VS_INPUT, matching the NV2A's 16
		// generic vertex attribute registers (identity mapping).
		pCurrentHostVertexElement->SemanticName = "TEXCOORD";
		pCurrentHostVertexElement->SemanticIndex = (UINT)VertexRegister;
		pCurrentHostVertexElement->InputSlotClass = D3D11_INPUT_PER_VERTEX_DATA;
		pCurrentHostVertexElement->InstanceDataStepRate = 0;

		pCurrentVertexShaderStreamInfo->HostVertexStride += HostVertexElementByteSize;

		return true;
	}

public:
	D3D11_INPUT_ELEMENT_DESC* Convert(xbox::X_VERTEXATTRIBUTEFORMAT* pXboxDeclaration, bool bIsFixedFunction, CxbxVertexDeclaration* pCxbxVertexDeclaration)
	{
		// Get a preprocessed copy of the original Xbox Vertex Declaration
		pCurrentVertexDeclaration = pCxbxVertexDeclaration;

		IsFixedFunction = bIsFixedFunction;

		RegVIsPresentInDeclaration.fill(false);

		// Mapping between Xbox register and the resulting host vertex element
		D3D11_INPUT_ELEMENT_DESC* HostVertexElementPerRegister[X_VSH_MAX_ATTRIBUTES] = { 0 };

		// Reserve the maximum number of elements, plus one for the null terminator :
		unsigned HostDeclarationSize = (X_VSH_MAX_ATTRIBUTES + 1) * sizeof(D3D11_INPUT_ELEMENT_DESC);

		D3D11_INPUT_ELEMENT_DESC* HostVertexElements = (D3D11_INPUT_ELEMENT_DESC*)calloc(1, HostDeclarationSize);
		pCurrentHostVertexElement = HostVertexElements;

		std::array<byte, X_VSH_MAX_ATTRIBUTES> orderedRegisterIndices;
		for (byte i = 0; i < orderedRegisterIndices.size(); i++)
			orderedRegisterIndices[i] = i;

		// Make sure we convert registers in order of offset, per stream
		// TODO fix elements with identical positions?
		std::sort(orderedRegisterIndices.begin(), orderedRegisterIndices.end(),
			[pXboxDeclaration](const auto& x, const auto& y)
			{
				auto regX = pXboxDeclaration->Slots[x];
				auto regY = pXboxDeclaration->Slots[y];
				return std::tie(regX.StreamIndex, regX.Offset)
					 < std::tie(regY.StreamIndex, regY.Offset);
			});

		EmuLog(LOG_LEVEL::DEBUG, "Parsing vertex declaration");
		for (size_t i = 0; i < orderedRegisterIndices.size(); i++) {
			auto regIndex = orderedRegisterIndices[i];
			auto &slot = pXboxDeclaration->Slots[regIndex];
			if (slot.Format > xbox::X_D3DVSDT_NONE) {
				// NOTE: The > (not >=) intentionally excludes X_D3DVSDT_NONE (0x02),
				// which is the Format value for AUTONORMAL/AUTOTEXCOORD tessellation
				// attributes. Do NOT change this to >= for D3D11: these attributes
				// have no data in the vertex buffer, so adding them to the input
				// layout would cause reads from non-existent VB data on normal draws.
				// PatchDraw handles auto-generated attributes separately by building
				// its own input layout and populating the data via CPU tessellation.
				// The tessellation metadata (autoNormalRegister, autoTexcoordRegister)
				// is recorded in CxbxVertexDeclaration from a separate loop below.
				// (D3D9 used D3DDECLMETHOD_CROSSUV/UV to tell the runtime to
				// generate data during DrawRectPatch.)
				//
				// FUTURE: If this moves to HS/DS hardware tessellation, the DS would
				// generate these outputs directly, eliminating this concern entirely.
				// See PatchDraw.cpp header comment for the full HS/DS design.
				if (VshConvertToken_STREAMDATA_REG(regIndex, slot)) {
					// Add this register to the list of declared registers
					RegVIsPresentInDeclaration[regIndex] = true;
					// Remember a pointer to this register
					HostVertexElementPerRegister[regIndex] = pCurrentHostVertexElement;

					EmuLog(LOG_LEVEL::DEBUG, "\tXbox Stream %d, Offset %d, Format %d, Slot %d",
						slot.StreamIndex, slot.Offset, slot.Format, regIndex);
					EmuLog(LOG_LEVEL::DEBUG, "\tHost InputSlot %d, Offset %d, Format %d, Semantic %s-%d",
						pCurrentHostVertexElement->InputSlot,
						pCurrentHostVertexElement->AlignedByteOffset,
						pCurrentHostVertexElement->Format,
						pCurrentHostVertexElement->SemanticName ? pCurrentHostVertexElement->SemanticName : "(null)",
						pCurrentHostVertexElement->SemanticIndex);
					pCurrentHostVertexElement++;
				}
			}
		}

		// Terminate with a zero-initialized element (SemanticName will be nullptr)
		memset(pCurrentHostVertexElement, 0, sizeof(D3D11_INPUT_ELEMENT_DESC));

		// Post-process tessellation elements (identified by TessellationType)
		for (int AttributeIndex = 0; AttributeIndex < X_VSH_MAX_ATTRIBUTES; AttributeIndex++) {
			auto pHostElement = HostVertexElementPerRegister[AttributeIndex];
			if (pHostElement == nullptr) continue;
			auto &slot = pXboxDeclaration->Slots[AttributeIndex];
			if (slot.TessellationType == 1 /* AUTONORMAL */) {
				int TessellationSource = slot.TessellationSource;
				auto pSourceElement = HostVertexElementPerRegister[TessellationSource];
				if (pSourceElement != nullptr) {
					// Copy over the InputSlot, AlignedByteOffset and Format of the source element :
					pHostElement->InputSlot = pSourceElement->InputSlot;
					pHostElement->AlignedByteOffset = pSourceElement->AlignedByteOffset;
					pHostElement->Format = pSourceElement->Format;
				}
			}
		}

		// Ensure valid ordering of the vertex declaration
		// In particular "All vertex elements for a stream must be consecutive and sorted by offset"
		// Test case: King Kong (due to register redefinition)
		// Note : Xbox slots might use non-ordered stream indices, so we can't rely on the output ordering of our converted elements!
		std::sort(/*First=*/HostVertexElements, /*Last=*/pCurrentHostVertexElement, /*Pred=*/[] (const auto& x, const auto& y)
			{ return std::tie(x.InputSlot, x.AlignedByteOffset) < std::tie(y.InputSlot, y.AlignedByteOffset); });

		// Record which registers are in the vertex declaration
		for (size_t i = 0; i < RegVIsPresentInDeclaration.size(); i++) {
			pCxbxVertexDeclaration->vRegisterInDeclaration[i] = RegVIsPresentInDeclaration[i];
		}

		// Record tessellation auto-generation metadata for PatchDraw.
		// These attributes have Format == X_D3DVSDT_NONE so they are NOT in the host vertex
		// declaration / input layout (the source VB has no data for them). PatchDraw uses this
		// metadata to compute and inject the auto-generated data at draw time.
		for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
			auto &slot = pXboxDeclaration->Slots[i];
			if (slot.Format == xbox::X_D3DVSDT_NONE) {
				if (slot.TessellationType == 1 /* AUTONORMAL */) {
					pCxbxVertexDeclaration->autoNormalRegister = i;
					pCxbxVertexDeclaration->autoNormalSourceRegister = slot.TessellationSource;
				} else if (slot.TessellationType == 2 /* AUTOTEXCOORD */) {
					pCxbxVertexDeclaration->autoTexcoordRegister = i;
				}
			}
		}

		return HostVertexElements;
	}
};

D3D11_INPUT_ELEMENT_DESC *EmuRecompileVshDeclaration
(
	xbox::X_VERTEXATTRIBUTEFORMAT* pXboxDeclaration,
   	bool                  bIsFixedFunction,
   	CxbxVertexDeclaration *pCxbxVertexDeclaration
)
{
	XboxVertexDeclarationConverter Converter;

	D3D11_INPUT_ELEMENT_DESC* pHostVertexElements = Converter.Convert(pXboxDeclaration, bIsFixedFunction, pCxbxVertexDeclaration);

   	return pHostVertexElements;
}


// parse xbox vertex shader function into an intermediate format
extern void EmuParseVshFunction
(
	// Pointer to raw Xbox vertex shader instruction slots
	DWORD* pXboxFunction,
	IntermediateVertexShader* pShader
)
{
	// Decode the vertex shader program tokens into an intermediate representation
	auto pCurToken = (uint32_t*)pXboxFunction;

	// Decode until we hit a token marked final
	// Note : CxbxSetVertexShaderSlots makes sure this always stops
	// after X_VSH_MAX_INSTRUCTION_COUNT, by setting FLD_FINAL in there.
	while (XboxVertexShaderDecoder::VshConvertToIntermediate(pCurToken, pShader)) {
		pCurToken += X_VSH_INSTRUCTION_SIZE;
	}
}
