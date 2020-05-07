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

#include "core\kernel\init\CxbxKrnl.h"
#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\Direct3D9\Direct3D9.h" // For g_Xbox_VertexShader_Handle
#include "core\hle\D3D8\Direct3D9\VertexShaderSource.h" // For g_VertexShaderSource
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbD3D8Logging.h" // For DEBUG_D3DRESULT
#include "core\hle\D3D8\XbConvert.h" // For NV2A_VP_UPLOAD_INST, NV2A_VP_UPLOAD_CONST_ID, NV2A_VP_UPLOAD_CONST
#include "devices\video\nv2a.h" // For D3DPUSH_DECODE
#include "common\Logging.h" // For LOG_INIT

#include "XbD3D8Types.h" // For X_D3DVSDE_*
#include <sstream>
#include <unordered_map>
#include <array>
#include <bitset>

// External symbols :
extern xbox::X_STREAMINPUT g_Xbox_SetStreamSource[X_VSH_MAX_STREAMS]; // Declared in XbVertexBuffer.cpp

// Variables set by [D3DDevice|CxbxImpl]_SetVertexShaderInput() :
                    unsigned g_Xbox_SetVertexShaderInput_Count = 0; // Read by GetXboxVertexAttributes
          xbox::X_STREAMINPUT g_Xbox_SetVertexShaderInput_Data[X_VSH_MAX_STREAMS] = { 0 }; // Active when g_Xbox_SetVertexShaderInput_Count > 0
xbox::X_VERTEXATTRIBUTEFORMAT g_Xbox_SetVertexShaderInput_Attributes = { 0 }; // Read by GetXboxVertexAttributes when g_Xbox_SetVertexShaderInput_Count > 0

// Reads the active Xbox stream input values (containing VertexBuffer, Offset and Stride) for the given stream number.
// (These values are set through SetStreamSource and can be overridden by SetVertexShaderInput.)
// TODO : Start using this function everywhere g_Xbox_SetStreamSource is accessed currently!
xbox::X_STREAMINPUT& GetXboxVertexStreamInput(unsigned StreamNumber)
{
	// If SetVertexShaderInput is active, it's arguments overrule those of SetStreamSource
	if (g_Xbox_SetVertexShaderInput_Count > 0) {
		return g_Xbox_SetVertexShaderInput_Data[StreamNumber];
	}

	return g_Xbox_SetStreamSource[StreamNumber];
}

#define DbgVshPrintf \
	LOG_CHECK_ENABLED(LOG_LEVEL::DEBUG) \
		if(g_bPrintfOn) printf

// ****************************************************************************
// * Vertex shader function recompiler
// ****************************************************************************

class XboxVertexShaderDecoder
{
private:
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

	static void VshConvertIntermediateParam(VSH_IMD_PARAMETER& Param,
		uint32_t* pShaderToken,
		VSH_FIELD_NAME FLD_MUX,
		VSH_FIELD_NAME FLD_NEG,
		uint16_t R,
		uint16_t V,
		uint16_t C)
	{
		Param.ParameterType = (VSH_PARAMETER_TYPE)VshGetField(pShaderToken, FLD_MUX);
		switch (Param.ParameterType) {
		case PARAM_R:
			Param.Address = R;
			break;
		case PARAM_V:
			Param.Address = V;
			break;
		case PARAM_C:
			Param.Address = C;
			break;
		default:
			LOG_TEST_CASE("parameter type unknown");
		}

		int d = FLD_NEG - FLD_A_NEG;
		Param.Neg = VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_NEG)) > 0;
		Param.Swizzle[0] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_X));
		Param.Swizzle[1] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_Y));
		Param.Swizzle[2] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_Z));
		Param.Swizzle[3] = (VSH_SWIZZLE)VshGetField(pShaderToken, (VSH_FIELD_NAME)(d + FLD_A_SWZ_W));
	}

	void VshAddIntermediateInstruction(
		uint32_t* pShaderToken,
		IntermediateVertexShader* pShader,
		VSH_MAC MAC,
		VSH_ILU ILU,
		VSH_IMD_OUTPUT_TYPE output_type,
		int16_t output_address,
		int8_t output_mask)
	{
		// Is the output mask set?
		if (output_mask == 0) {
			return;
		}

		if (pShader->Instructions.size() >= VSH_MAX_INTERMEDIATE_COUNT) {
			CxbxKrnlCleanup("Shader exceeds conversion buffer!");
		}

		VSH_INTERMEDIATE_FORMAT intermediate;
		intermediate.MAC = MAC;
		intermediate.ILU = ILU;
		intermediate.Output.Type = output_type;
		intermediate.Output.Address = output_address;
		intermediate.Output.Mask = output_mask;
		// Get a0.x indirect constant addressing
		intermediate.IndexesWithA0_X = VshGetField(pShaderToken, FLD_A0X) > 0; // Applies to PARAM_C parameter reads

		int16_t R;
		int16_t V = VshGetField(pShaderToken, FLD_V);
		int16_t C = ConvertCRegister(VshGetField(pShaderToken, FLD_CONST));
		intermediate.ParamCount = 0;
		if (MAC >= MAC_MOV) {
			// Get parameter A
			R = VshGetField(pShaderToken, FLD_A_R);
			VshConvertIntermediateParam(intermediate.Parameters[intermediate.ParamCount++], pShaderToken, FLD_A_MUX, FLD_A_NEG, R, V, C);
		}

		if ((MAC == MAC_MUL) || ((MAC >= MAC_MAD) && (MAC <= MAC_SGE))) {
			// Get parameter B
			R = VshGetField(pShaderToken, FLD_B_R);
			VshConvertIntermediateParam(intermediate.Parameters[intermediate.ParamCount++], pShaderToken, FLD_B_MUX, FLD_B_NEG, R, V, C);
		}

		if ((ILU >= ILU_MOV) || (MAC == MAC_ADD) || (MAC == MAC_MAD)) {
			// Get parameter C
			R = VshGetField(pShaderToken, FLD_C_R_HIGH) << 2 | VshGetField(pShaderToken, FLD_C_R_LOW);
			VshConvertIntermediateParam(intermediate.Parameters[intermediate.ParamCount++], pShaderToken, FLD_C_MUX, FLD_C_NEG, R, V, C);
		}

		// Add the instruction to the shader
		pShader->Instructions.push_back(intermediate);
	}

public:
	bool VshConvertToIntermediate(uint32_t* pShaderToken, IntermediateVertexShader* pShader)
	{
		// First get the instruction(s).
		VSH_ILU ILU = (VSH_ILU)VshGetField(pShaderToken, FLD_ILU);
		VSH_MAC MAC = (VSH_MAC)VshGetField(pShaderToken, FLD_MAC);
		if (MAC > MAC_ARL) LOG_TEST_CASE("Unknown MAC");

		// Output register
		VSH_OUTPUT_MUX OutputMux = (VSH_OUTPUT_MUX)VshGetField(pShaderToken, FLD_OUT_MUX);
		int16_t OutputAddress = VshGetField(pShaderToken, FLD_OUT_ADDRESS);
		VSH_IMD_OUTPUT_TYPE OutputType;
		if ((VSH_OUTPUT_TYPE)VshGetField(pShaderToken, FLD_OUT_ORB) == OUTPUT_C) {
			OutputType = IMD_OUTPUT_C;
			OutputAddress = ConvertCRegister(OutputAddress);
		} else { // OUTPUT_O:
			OutputType = IMD_OUTPUT_O;
			OutputAddress = OutputAddress & 0xF;
		}

		// MAC,ILU output R register
		int16_t RAddress = VshGetField(pShaderToken, FLD_OUT_R);

		// Test for paired opcodes
		bool bIsPaired = (MAC != MAC_NOP) && (ILU != ILU_NOP);

		// Check if there's a MAC opcode
		if (MAC > MAC_NOP && MAC <= MAC_ARL) {
			if (bIsPaired && RAddress == 1) {
				// Ignore paired MAC opcodes that write to R1
			} else {
				if (MAC == MAC_ARL) {
					VshAddIntermediateInstruction(pShaderToken, pShader, MAC, ILU_NOP, IMD_OUTPUT_A0X, 0, MASK_X);
				} else {
					VshAddIntermediateInstruction(pShaderToken, pShader, MAC, ILU_NOP, IMD_OUTPUT_R, RAddress, VshGetField(pShaderToken, FLD_OUT_MAC_MASK));
				}
			}

			// Check if we must add a muxed MAC opcode as well
			if (OutputMux == OMUX_MAC) {
				VshAddIntermediateInstruction(pShaderToken, pShader, MAC, ILU_NOP, OutputType, OutputAddress, VshGetField(pShaderToken, FLD_OUT_O_MASK));
			}
		}

		// Check if there's an ILU opcode
		if (ILU != ILU_NOP) {
			// Paired ILU opcodes will only write to R1
			VshAddIntermediateInstruction(pShaderToken, pShader, MAC_NOP, ILU, IMD_OUTPUT_R, bIsPaired ? 1 : RAddress, VshGetField(pShaderToken, FLD_OUT_ILU_MASK));
			// Check if we must add a muxed ILU opcode as well
			if (OutputMux == OMUX_ILU) {
				VshAddIntermediateInstruction(pShaderToken, pShader, MAC_NOP, ILU, OutputType, OutputAddress, VshGetField(pShaderToken, FLD_OUT_O_MASK));
			}
		}

		return VshGetField(pShaderToken, FLD_FINAL) == 0;
	}

};

// ****************************************************************************
// * Vertex shader declaration recompiler
// ****************************************************************************

extern D3DCAPS g_D3DCaps;

class XboxVertexDeclarationConverter
{
protected:
	// Internal variables
	CxbxVertexDeclaration* pVertexDeclarationToSet;
	CxbxVertexShaderStreamInfo* pCurrentVertexShaderStreamInfo = nullptr;
	bool IsFixedFunction;
	D3DVERTEXELEMENT* pRecompiled;
	std::array<bool, 16> RegVIsPresentInDeclaration;

private:
	#define D3DDECLUSAGE_UNSUPPORTED ((D3DDECLUSAGE)-1)

	D3DDECLUSAGE Xb2PCRegisterType(DWORD VertexRegister, BYTE &UsageIndex)
	{
		UsageIndex = 0;
		switch (VertexRegister) {
			case xbox::X_D3DVSDE_POSITION    /*= 0*/: return (pRecompiled->Type == D3DDECLTYPE_FLOAT4) ? D3DDECLUSAGE_POSITIONT : D3DDECLUSAGE_POSITION;
			case xbox::X_D3DVSDE_BLENDWEIGHT /*= 1*/: return D3DDECLUSAGE_BLENDWEIGHT;
			case xbox::X_D3DVSDE_NORMAL      /*= 2*/: return D3DDECLUSAGE_NORMAL;
			case xbox::X_D3DVSDE_DIFFUSE     /*= 3*/: return D3DDECLUSAGE_COLOR;
			case xbox::X_D3DVSDE_SPECULAR    /*= 4*/: UsageIndex = 1; return D3DDECLUSAGE_COLOR;
			case xbox::X_D3DVSDE_FOG         /*= 5*/: return D3DDECLUSAGE_FOG;
			case xbox::X_D3DVSDE_POINTSIZE   /*= 6*/: return D3DDECLUSAGE_PSIZE;
			case xbox::X_D3DVSDE_BACKDIFFUSE /*= 7*/: UsageIndex = 2; return D3DDECLUSAGE_COLOR;
			case xbox::X_D3DVSDE_BACKSPECULAR/*= 8*/: UsageIndex = 3; return D3DDECLUSAGE_COLOR;
			case xbox::X_D3DVSDE_TEXCOORD0   /*= 9*/: return D3DDECLUSAGE_TEXCOORD;
			case xbox::X_D3DVSDE_TEXCOORD1   /*=10*/: UsageIndex = 1; return D3DDECLUSAGE_TEXCOORD;
			case xbox::X_D3DVSDE_TEXCOORD2   /*=11*/: UsageIndex = 2; return D3DDECLUSAGE_TEXCOORD;
			case xbox::X_D3DVSDE_TEXCOORD3   /*=12*/: UsageIndex = 3; return D3DDECLUSAGE_TEXCOORD;
			default /*13-15*/ : return D3DDECLUSAGE_UNSUPPORTED;
		}
	}

	// VERTEX SHADER

#if 0 // TODO : Support tesselation for X_VERTEXSHADERINPUT's
	void VshConvertToken_TESSELATOR(DWORD *pXboxToken)
	{
		BYTE Index;

		if(*pXboxToken & X_D3DVSD_MASK_TESSUV)
		{
			DWORD VertexRegister    = VshGetVertexRegister(*pXboxToken);
			DWORD NewVertexRegister = VertexRegister;

			NewVertexRegister = Xb2PCRegisterType(VertexRegister, Index);
			// TODO : Expand on the setting of this TESSUV register element :
			pRecompiled->Usage = D3DDECLUSAGE(NewVertexRegister);
			pRecompiled->UsageIndex = Index;
		}
		else // D3DVSD_TESSNORMAL
		{
			DWORD VertexRegisterIn  = VshGetVertexRegisterIn(*pXboxToken);
			DWORD VertexRegisterOut = VshGetVertexRegister(*pXboxToken);

			DWORD NewVertexRegisterIn  = VertexRegisterIn;
			DWORD NewVertexRegisterOut = VertexRegisterOut;

			NewVertexRegisterIn = Xb2PCRegisterType(VertexRegisterIn, Index);
			// TODO : Expand on the setting of this TESSNORMAL input register element :
			pRecompiled->Usage = D3DDECLUSAGE(NewVertexRegisterIn);
			pRecompiled->UsageIndex = Index;

			NewVertexRegisterOut = Xb2PCRegisterType(VertexRegisterOut, Index);
			// TODO : Expand on the setting of this TESSNORMAL output register element :
			pRecompiled++;
			pRecompiled->Usage = D3DDECLUSAGE(NewVertexRegisterOut);
			pRecompiled->UsageIndex = Index;
		}
	}
#endif
	void VshConvertToken_STREAM(DWORD StreamNumber)
	{
		// new stream
		pCurrentVertexShaderStreamInfo = &(pVertexDeclarationToSet->VertexStreams[StreamNumber]);
		pCurrentVertexShaderStreamInfo->NeedPatch = FALSE;
		pCurrentVertexShaderStreamInfo->DeclPosition = FALSE;
		pCurrentVertexShaderStreamInfo->CurrentStreamNumber = 0;
		pCurrentVertexShaderStreamInfo->HostVertexStride = 0;
		pCurrentVertexShaderStreamInfo->NumberOfVertexElements = 0;

		// Dxbx note : Use Dophin(s), FieldRender, MatrixPaletteSkinning and PersistDisplay as a testcase

		pCurrentVertexShaderStreamInfo->CurrentStreamNumber = (WORD)StreamNumber;
		pVertexDeclarationToSet->NumberOfVertexStreams++;
		// TODO : Keep a bitmask for all StreamNumber's seen?
	}

	void VshConvert_RegisterVertexElement(
		UINT XboxVertexElementDataType,
		UINT XboxVertexElementByteSize,
		UINT HostVertexElementByteSize,
		BOOL NeedPatching)
	{
		CxbxVertexShaderStreamElement* pCurrentElement = &(pCurrentVertexShaderStreamInfo->VertexElements[pCurrentVertexShaderStreamInfo->NumberOfVertexElements]);
		pCurrentElement->XboxType = XboxVertexElementDataType;
		pCurrentElement->XboxByteSize = XboxVertexElementByteSize;
		pCurrentElement->HostByteSize = HostVertexElementByteSize;
		pCurrentVertexShaderStreamInfo->NumberOfVertexElements++;
		pCurrentVertexShaderStreamInfo->NeedPatch |= NeedPatching;
	}

	void VshConvertToken_STREAMDATA_REG(DWORD VertexRegister, xbox::X_VERTEXSHADERINPUT &slot)
	{
		BOOL NeedPatching = FALSE;

		// Add this register to the list of declared registers
		RegVIsPresentInDeclaration[VertexRegister] = true;

		DWORD XboxVertexElementDataType = slot.Format;
		WORD XboxVertexElementByteSize = 0;
		BYTE HostVertexElementDataType = 0;
		WORD HostVertexElementByteSize = 0;

		switch (XboxVertexElementDataType)
		{
		case xbox::X_D3DVSDT_FLOAT1: // 0x12:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT1;
			HostVertexElementByteSize = 1 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_FLOAT2: // 0x22:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT2;
			HostVertexElementByteSize = 2 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_FLOAT3: // 0x32:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
			HostVertexElementByteSize = 3 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_FLOAT4: // 0x42:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
			HostVertexElementByteSize = 4 * sizeof(FLOAT);
			break;
		case xbox::X_D3DVSDT_D3DCOLOR: // 0x40:
			HostVertexElementDataType = D3DDECLTYPE_D3DCOLOR;
			HostVertexElementByteSize = 1 * sizeof(D3DCOLOR);
			break;
		case xbox::X_D3DVSDT_SHORT2: // 0x25:
			HostVertexElementDataType = D3DDECLTYPE_SHORT2;
			HostVertexElementByteSize = 2 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_SHORT4: // 0x45:
			HostVertexElementDataType = D3DDECLTYPE_SHORT4;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			break;
		case xbox::X_D3DVSDT_NORMSHORT1: // 0x11:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT2N;
				HostVertexElementByteSize = 2 * sizeof(SHORT);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT1;
				HostVertexElementByteSize = 1 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 1 * sizeof(xbox::SHORT);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_NORMSHORT2: // 0x21:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT2N;
				HostVertexElementByteSize = 2 * sizeof(SHORT);
				// No need for patching in D3D9
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT2;
				HostVertexElementByteSize = 2 * sizeof(FLOAT);
				XboxVertexElementByteSize = 2 * sizeof(xbox::SHORT);
				NeedPatching = TRUE;
			}
			break;
		case xbox::X_D3DVSDT_NORMSHORT3: // 0x31:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT4N;
				HostVertexElementByteSize = 4 * sizeof(SHORT);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
				HostVertexElementByteSize = 3 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 3 * sizeof(xbox::SHORT);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_NORMSHORT4: // 0x41:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT4N;
				HostVertexElementByteSize = 4 * sizeof(SHORT);
				// No need for patching in D3D9
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
				HostVertexElementByteSize = 4 * sizeof(FLOAT);
				XboxVertexElementByteSize = 4 * sizeof(xbox::SHORT);
				NeedPatching = TRUE;
			}
			break;
		case xbox::X_D3DVSDT_NORMPACKED3: // 0x16:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
			HostVertexElementByteSize = 3 * sizeof(FLOAT);
			XboxVertexElementByteSize = 1 * sizeof(xbox::DWORD);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_SHORT1: // 0x15:
			HostVertexElementDataType = D3DDECLTYPE_SHORT2;
			HostVertexElementByteSize = 2 * sizeof(SHORT);
			XboxVertexElementByteSize = 1 * sizeof(xbox::SHORT);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_SHORT3: // 0x35:
			HostVertexElementDataType = D3DDECLTYPE_SHORT4;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			XboxVertexElementByteSize = 3 * sizeof(xbox::SHORT);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_PBYTE1: // 0x14:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT1;
				HostVertexElementByteSize = 1 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 1 * sizeof(xbox::BYTE);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_PBYTE2: // 0x24:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT2;
				HostVertexElementByteSize = 2 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 2 * sizeof(xbox::BYTE);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_PBYTE3: // 0x34:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
				HostVertexElementByteSize = 3 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 3 * sizeof(xbox::BYTE);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_PBYTE4: // 0x44:
			// Test-case : Panzer
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
				// No need for patching when D3D9 supports D3DDECLTYPE_UBYTE4N
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
				HostVertexElementByteSize = 4 * sizeof(FLOAT);
				XboxVertexElementByteSize = 4 * sizeof(xbox::BYTE);
				NeedPatching = TRUE;
			}
			break;
		case xbox::X_D3DVSDT_FLOAT2H: // 0x72:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
			HostVertexElementByteSize = 4 * sizeof(FLOAT);
			XboxVertexElementByteSize = 3 * sizeof(FLOAT);
			NeedPatching = TRUE;
			break;
		case xbox::X_D3DVSDT_NONE: // 0x02:
			// No host element data, so no patching
			break;
		default:
			//LOG_TEST_CASE("Unknown data type for D3DVSD_REG: 0x%02X\n", XboxVertexElementDataType);
			break;
		}

		// On X_D3DVSDT_NONE skip this token
		if (XboxVertexElementDataType == xbox::X_D3DVSDT_NONE)
		{
			// Xbox elements with X_D3DVSDT_NONE have size zero, so there's no need to register those.
			return;
		}

		// Select new stream, if needed
		if ((pCurrentVertexShaderStreamInfo == nullptr)
		 || (pCurrentVertexShaderStreamInfo->CurrentStreamNumber != slot.IndexOfStream)) {
			VshConvertToken_STREAM(slot.IndexOfStream);
		}

		// save patching information
		VshConvert_RegisterVertexElement(
			XboxVertexElementDataType,
			NeedPatching ? XboxVertexElementByteSize : HostVertexElementByteSize,
			HostVertexElementByteSize,
			NeedPatching);

/*
	typedef struct _D3DVERTEXELEMENT9
	{
		WORD    Stream;     // Stream index
		WORD    Offset;     // Offset in the stream in bytes
		BYTE    Type;       // Data type
		BYTE    Method;     // Processing method
		BYTE    Usage;      // Semantics
		BYTE    UsageIndex; // Semantic index
	} D3DVERTEXELEMENT9, * LPD3DVERTEXELEMENT9;
*/

		pRecompiled->Stream = pCurrentVertexShaderStreamInfo->CurrentStreamNumber;
		pRecompiled->Offset = pCurrentVertexShaderStreamInfo->HostVertexStride;
		pRecompiled->Type = HostVertexElementDataType;
		pRecompiled->Method = D3DDECLMETHOD_DEFAULT;
		if (IsFixedFunction) {
			pRecompiled->Usage = Xb2PCRegisterType(VertexRegister, /*&*/pRecompiled->UsageIndex);
		}
		else {
			// D3DDECLUSAGE_TEXCOORD can be useds for any user-defined data
			// We need this because there is no reliable way to detect the real usage
			// Xbox has no concept of 'usage types', it only requires a list of attribute register numbers.
			// So we treat them all as 'user-defined' with an Index of the Vertex Register Index
			// this prevents information loss in shaders due to non-matching dcl types!
			pRecompiled->Usage = D3DDECLUSAGE_TEXCOORD;
			pRecompiled->UsageIndex = (BYTE)VertexRegister;
		}

		pRecompiled++;

		pCurrentVertexShaderStreamInfo->HostVertexStride += HostVertexElementByteSize;
	}

public:
	D3DVERTEXELEMENT* Convert(xbox::X_VERTEXATTRIBUTEFORMAT* pXboxDeclaration, bool bIsFixedFunction, CxbxVertexDeclaration* pCxbxVertexDeclaration)
	{
		// Get a preprocessed copy of the original Xbox Vertex Declaration
		pVertexDeclarationToSet = pCxbxVertexDeclaration;

		IsFixedFunction = bIsFixedFunction;

		RegVIsPresentInDeclaration.fill(false);

		// For Direct3D9, we need to reserve at least twice the number of elements, as one token can generate two registers (in and out) :
		unsigned HostDeclarationSize = X_VSH_MAX_ATTRIBUTES * sizeof(D3DVERTEXELEMENT) * 2;

		D3DVERTEXELEMENT* Result = (D3DVERTEXELEMENT*)calloc(1, HostDeclarationSize);
		pRecompiled = Result;
		uint8_t* pRecompiledBufferOverflow = ((uint8_t*)pRecompiled) + HostDeclarationSize;


		for (size_t i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
			auto &slot = pXboxDeclaration->Slots[i];
			if (slot.Format > 0) {
				// Set Direct3D9 vertex element (declaration) members :
				VshConvertToken_STREAMDATA_REG(i, slot);
			}
		}

		*pRecompiled = D3DDECL_END();

		// Ensure valid ordering of the vertex declaration (http://doc.51windows.net/Directx9_SDK/graphics/programmingguide/gettingstarted/vertexdeclaration/vertexdeclaration.htm)
		// In particular "All vertex elements for a stream must be consecutive and sorted by offset"
		// Test case: King Kong (due to register redefinition)
		std::sort(Result, pRecompiled, [] (const auto& x, const auto& y)
			{ return std::tie(x.Stream, x.Method, x.Offset) < std::tie(y.Stream, y.Method, y.Offset); });

		// Record which registers are in the vertex declaration
		for (size_t i = 0; i < RegVIsPresentInDeclaration.size(); i++) {
			pCxbxVertexDeclaration->vRegisterInDeclaration[i] = RegVIsPresentInDeclaration[i];
		}

		return Result;
	}
};

D3DVERTEXELEMENT *EmuRecompileVshDeclaration
(
	xbox::X_VERTEXATTRIBUTEFORMAT* pXboxDeclaration,
    bool                  bIsFixedFunction,
    CxbxVertexDeclaration *pCxbxVertexDeclaration
)
{
	XboxVertexDeclarationConverter Converter;

	D3DVERTEXELEMENT* pHostVertexElements = Converter.Convert(pXboxDeclaration, bIsFixedFunction, pCxbxVertexDeclaration);

    return pHostVertexElements;
}

extern void FreeVertexDynamicPatch(CxbxVertexShader *pVertexShader)
{
    pVertexShader->Declaration.NumberOfVertexStreams = 0;
}

// Checks for failed vertex shaders, and shaders that would need patching
boolean VshHandleIsValidShader(DWORD XboxVertexShaderHandle)
{
#if 0
	//printf( "VS = 0x%.08X\n", XboxVertexShaderHandle );

    CxbxVertexShader *pCxbxVertexShader = GetCxbxVertexShader(XboxVertexShaderHandle);
    if (pCxbxVertexShader) {
        if (pCxbxVertexShader->XboxStatus != 0)
        {
            return FALSE;
        }
        /*
        for (uint32 i = 0; i < pCxbxVertexShader->VertexShaderInfo.NumberOfVertexStreams; i++)
        {
            if (pCxbxVertexShader->VertexShaderInfo.VertexStreams[i].NeedPatch)
            {
                // Just for caching purposes
                pCxbxVertexShader->XboxStatus = 0x80000001;
                return FALSE;
            }
        }
        */
    }
#endif
    return TRUE;
}

extern boolean IsValidCurrentShader(void)
{
	// Dxbx addition : There's no need to call
	// XTL_EmuIDirect3DDevice_GetVertexShader, just check g_Xbox_VertexShader_Handle :
	return VshHandleIsValidShader(g_Xbox_VertexShader_Handle);
}

// Vertex shader state
static DWORD g_CxbxVertexShaderSlotAddress = 0;
static DWORD g_CxbxVertexShaderSlots[X_VSH_MAX_INSTRUCTION_COUNT * X_VSH_INSTRUCTION_SIZE] = { 0 };

DWORD* GetCxbxVertexShaderSlotPtr(const DWORD SlotIndexAddress)
{
	if (SlotIndexAddress < X_VSH_MAX_INSTRUCTION_COUNT) {
		return &g_CxbxVertexShaderSlots[SlotIndexAddress * X_VSH_INSTRUCTION_SIZE];
	} else {
		LOG_TEST_CASE("SlotIndexAddress out of range"); // FIXME : extend with value (once supported by LOG_TEST_CASE)
		return nullptr;
	}
}

CxbxVertexDeclaration *GetCxbxVertexDeclaration(DWORD XboxVertexShaderHandle)
{
    CxbxVertexShader *pCxbxVertexShader = GetCxbxVertexShader(XboxVertexShaderHandle);

    for (uint32_t i = 0; i < pCxbxVertexShader->Declaration.NumberOfVertexStreams; i++)
    {
        if (pCxbxVertexShader->Declaration.VertexStreams[i].NeedPatch)
        {
            return &pCxbxVertexShader->Declaration;
        }
    }
    return nullptr;
}

std::unordered_map<DWORD, CxbxVertexShader*> g_CxbxVertexShaders;

CxbxVertexShader* GetCxbxVertexShader(DWORD XboxVertexShaderHandle)
{
	if (VshHandleIsVertexShader(XboxVertexShaderHandle)) {
		auto it = g_CxbxVertexShaders.find(XboxVertexShaderHandle);
		if (it != g_CxbxVertexShaders.end()) {
			return it->second;
		}
	}

	return nullptr;
}

void RegisterCxbxVertexShader(DWORD XboxVertexShaderHandle, CxbxVertexShader* shader)
{
	auto it = g_CxbxVertexShaders.find(XboxVertexShaderHandle);
	if (it != g_CxbxVertexShaders.end() && it->second != nullptr && shader != nullptr) {
		LOG_TEST_CASE("Overwriting existing Vertex Shader");
	}

	g_CxbxVertexShaders[XboxVertexShaderHandle] = shader;
}

IDirect3DVertexDeclaration* CxbxCreateVertexDeclaration(D3DVERTEXELEMENT *pDeclaration)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	IDirect3DVertexDeclaration* pHostVertexDeclaration = nullptr;
	HRESULT hRet = g_pD3DDevice->CreateVertexDeclaration(pDeclaration, &pHostVertexDeclaration);
	DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateVertexDeclaration");
	return pHostVertexDeclaration;
}

// Converts an Xbox FVF shader handle to X_VERTEXATTRIBUTEFORMAT
// Note : Temporary, until we reliably locate the Xbox internal state for this
// See D3DXDeclaratorFromFVF docs https://docs.microsoft.com/en-us/windows/win32/direct3d9/d3dxdeclaratorfromfvf
// and https://github.com/reactos/wine/blob/2e8dfbb1ad71f24c41e8485a39df01bb9304127f/dlls/d3dx9_36/mesh.c#L2041
xbox::X_VERTEXATTRIBUTEFORMAT XboxFVFToXboxVertexAttributeFormat(DWORD xboxFvf)
{
	using namespace xbox;

	X_VERTEXATTRIBUTEFORMAT declaration = { 0 }; // FVFs don't tesselate, all slots read from stream zero

	static DWORD X_D3DVSDT_FLOAT[] = { 0, X_D3DVSDT_FLOAT1, X_D3DVSDT_FLOAT2, X_D3DVSDT_FLOAT3, X_D3DVSDT_FLOAT4 };

	static const DWORD InvalidXboxFVFBits = X_D3DFVF_RESERVED0 | X_D3DFVF_RESERVED1 /* probably D3DFVF_PSIZE if detected */
		| 0x0000F000 // Bits between texture count and the texture formats
		| 0xFF000000; // All bits above the four alllowed texture formats

	if (xboxFvf & InvalidXboxFVFBits) {
		LOG_TEST_CASE("Invalid Xbox FVF bits detected!");
	}

	// Position & Blendweights
	int nrPositionFloats = 3;
	int nrBlendWeights = 0;
	unsigned offset = 0;
	DWORD position = (xboxFvf & X_D3DFVF_POSITION_MASK);
	switch (position) {
		case 0: nrPositionFloats = 0; LOG_TEST_CASE("FVF without position"); break; // Note : Remove logging if this occurs often
		case X_D3DFVF_XYZ: /*nrPositionFloats is set to 3 by default*/ break;
		case X_D3DFVF_XYZRHW: nrPositionFloats = 4; break;
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
		pSlot = &declaration.Slots[X_D3DVSDE_DIFFUSE];
		pSlot->Format = X_D3DVSDT_D3DCOLOR;
		pSlot->Offset = offset;
		offset += sizeof(DWORD) * 1;
	}

	if (xboxFvf & X_D3DFVF_SPECULAR) {
		pSlot = &declaration.Slots[X_D3DVSDE_SPECULAR];
		pSlot->Format = X_D3DVSDT_D3DCOLOR;
		pSlot->Offset = offset;
		offset += sizeof(DWORD) * 1;
	}

	// Write Texture Coordinates
	int textureCount = (xboxFvf & X_D3DFVF_TEXCOUNT_MASK) >> X_D3DFVF_TEXCOUNT_SHIFT;
	assert(textureCount <= 4); // Safeguard, since the X_D3DFVF_TEXCOUNT bitfield could contain invalid values (5 up to 15)
	for (int i = 0; i < textureCount; i++) {
		int numberOfCoordinates = 0;
		auto FVFTextureFormat = (xboxFvf >> X_D3DFVF_TEXCOORDSIZE_SHIFT(i)) & 0x003;
		switch (FVFTextureFormat) {
			case X_D3DFVF_TEXTUREFORMAT1: numberOfCoordinates = 1; break;
			case X_D3DFVF_TEXTUREFORMAT2: numberOfCoordinates = 2; break;
			case X_D3DFVF_TEXTUREFORMAT3: numberOfCoordinates = 3; break;
			case X_D3DFVF_TEXTUREFORMAT4: numberOfCoordinates = 4; break;
			DEFAULT_UNREACHABLE;
		}

		assert(numberOfCoordinates > 0);
		pSlot = &declaration.Slots[X_D3DVSDE_TEXCOORD0 + i];
		pSlot->Format = X_D3DVSDT_FLOAT[numberOfCoordinates];
		pSlot->Offset = offset;
		offset += sizeof(float) * numberOfCoordinates;
	}

	return declaration;
}

void SetCxbxVertexDeclaration(CxbxVertexDeclaration& pCxbxVertexDeclaration) {
	LOG_INIT

	HRESULT hRet;

	// Set vertex declaration
	hRet = g_pD3DDevice->SetVertexDeclaration(pCxbxVertexDeclaration.pHostVertexDeclaration);
	DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexDeclaration");

	// Titles can specify default values for registers via calls like SetVertexData4f
	// HLSL shaders need to know whether to use vertex data or default vertex shader values
	// Any register not in the vertex declaration should be set to the default value
	float vertexDefaultFlags[X_VSH_MAX_ATTRIBUTES];
	for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
		vertexDefaultFlags[i] = pCxbxVertexDeclaration.vRegisterInDeclaration[i] ? 0.0f : 1.0f;
	}
	g_pD3DDevice->SetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_BASE, vertexDefaultFlags, 4);
}

// TODO Call this when state is dirty in UpdateNativeD3DResources
// Rather than every time state changes
void SetVertexShaderFromSlots()
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	auto pTokens = GetCxbxVertexShaderSlotPtr(g_CxbxVertexShaderSlotAddress);
	if (pTokens) {
		// Create a vertex shader from the tokens
		DWORD shaderSize;
		auto shaderKey = g_VertexShaderSource.CreateShader(pTokens, &shaderSize);
		IDirect3DVertexShader* pHostVertexShader = g_VertexShaderSource.GetShader(shaderKey);
		HRESULT hRet = g_pD3DDevice->SetVertexShader(pHostVertexShader);
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShader");
	}
}

void SetCxbxVertexShaderHandle(CxbxVertexShader* pCxbxVertexShader)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	HRESULT hRet;

	// Get vertex shader if we have a key
	auto pHostShader = pCxbxVertexShader->VertexShaderKey
		? g_VertexShaderSource.GetShader(pCxbxVertexShader->VertexShaderKey)
		: nullptr;

	// Set vertex shader
	hRet = g_pD3DDevice->SetVertexShader(pHostShader);
	DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShader");

	SetCxbxVertexDeclaration(pCxbxVertexShader->Declaration);
}

void CxbxSetVertexShaderSlots(DWORD* pTokens, DWORD Address, DWORD NrInstructions)
{
	int upToSlot = Address + NrInstructions;
	if (upToSlot > X_VSH_MAX_INSTRUCTION_COUNT) {
		LOG_TEST_CASE("Shader does not fit in vertex shader slots");
		return;
	}

	auto CxbxVertexShaderSlotPtr = GetCxbxVertexShaderSlotPtr(Address);
	if (CxbxVertexShaderSlotPtr == nullptr) {
		return;
	}

	memcpy(CxbxVertexShaderSlotPtr, pTokens, NrInstructions * X_VSH_INSTRUCTION_SIZE_BYTES);
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

		X_D3DVertexShader* pXboxVertexShader = VshHandleToXboxVertexShader(Handle);
		assert(pXboxVertexShader);

		// Xbox DOES store the Handle, but since it merely returns this through (unpatched) D3DDevice_GetVertexShaderInput, we don't have to.

		g_Xbox_SetVertexShaderInput_Count = StreamCount; // This > 0 indicates g_Xbox_SetVertexShaderInput_Data has to be used
		if (StreamCount > 0) {
			assert(StreamCount <= X_VSH_MAX_STREAMS);
			assert(pStreamInputs != xbnullptr);
			memcpy(g_Xbox_SetVertexShaderInput_Data, pStreamInputs, StreamCount * sizeof(xbox::X_STREAMINPUT)); // Make a copy of the supplied StreamInputs array
		}

		g_Xbox_SetVertexShaderInput_Attributes = pXboxVertexShader->VertexAttribute; // Copy this vertex shaders's attribute slots
	}
}

// Note : SelectVertexShaderDirect needs no EMUPATCH CxbxImpl_..., since it just calls SelectVertexShader

IDirect3DVertexDeclaration* CxbxGetVertexDeclaration(DWORD Handle)
{
	LOG_INIT;

	xbox::X_VERTEXATTRIBUTEFORMAT XboxVertexAttributeFormat;

	// TODO : Cache resulting declarations from given inputs

	bool bIsFixedFunction = VshHandleIsFVF(Handle);
	if (bIsFixedFunction) {
		XboxVertexAttributeFormat = XboxFVFToXboxVertexAttributeFormat(Handle);
	}
	else {
		xbox::X_D3DVertexShader* pXboxVertexShader = VshHandleToXboxVertexShader(Handle);
		XboxVertexAttributeFormat = pXboxVertexShader->VertexAttribute;
	}

	CxbxVertexShader* pCxbxVertexShader = (CxbxVertexShader*)calloc(1, sizeof(CxbxVertexShader));

	D3DVERTEXELEMENT* pRecompiledDeclaration = nullptr;
	pRecompiledDeclaration = EmuRecompileVshDeclaration(
		&XboxVertexAttributeFormat,
		bIsFixedFunction,
		&pCxbxVertexShader->Declaration);

	// Create the vertex declaration
	pCxbxVertexShader->Declaration.pHostVertexDeclaration = CxbxCreateVertexDeclaration(pRecompiledDeclaration);

	free(pRecompiledDeclaration);

	return pCxbxVertexShader->Declaration.pHostVertexDeclaration;
}

void CxbxImpl_SelectVertexShader(DWORD Handle, DWORD Address)
{
	LOG_INIT; // 
	// Address always indicates a previously loaded vertex shader slot (from where the program is used).
	// Handle can be null if the current Xbox VertexShader is assigned
	// Handle can be an address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// If Handle is assigned, it becomes the new current Xbox VertexShader,
	// which resets a bit of state (nv2a execution mode, viewport, ?)
	// Either way, the given address slot is selected as the start of the current vertex shader program
	g_CxbxVertexShaderSlotAddress = Address;

	if (Handle) {
		if (!VshHandleIsVertexShader(Handle))
			LOG_TEST_CASE("Non-zero handle must be a VertexShader!");

		g_Xbox_VertexShader_Handle = Handle;

		// Set the shader handle declaration
	}

	IDirect3DVertexDeclaration * pHostVertexDeclaration = CxbxGetVertexDeclaration(Handle);
	HRESULT hRet = g_pD3DDevice->SetVertexDeclaration(pHostVertexDeclaration);

#if 0
	if (pFunction != xbnullptr)
	{
		// Parse and compile the shader
		DWORD xboxFunctionSize = 0;
		pCxbxVertexShader->VertexShaderKey = g_VertexShaderSource.CreateShader(pFunction, &xboxFunctionSize);
	}

	// Save the status, to remove things later
	// pCxbxVertexShader->XboxStatus = hRet; // Not even used by VshHandleIsValidShader()

	RegisterCxbxVertexShader(*pHandle, pCxbxVertexShader);
#endif

	SetVertexShaderFromSlots();
}

void CxbxImpl_LoadVertexShaderProgram(CONST DWORD* pFunction, DWORD Address)
{
	// pFunction is a X_VSH_SHADER_HEADER pointer
	// D3DDevice_LoadVertexShaderProgram splits the given function buffer into batch-wise pushes to the NV2A
	// However, we can suffice by copying the program into our slots (and make sure these slots get converted into a vertex shader)

	// Copy shader instructions to shader slots
	auto shaderHeader = *((xbox::X_VSH_SHADER_HEADER*) pFunction);
	if (shaderHeader.Version != VERSION_XVS)
		LOG_TEST_CASE("Non-regular (state or read/write) shader detected!");

	auto tokens = (DWORD*)&pFunction[1];
	CxbxSetVertexShaderSlots(tokens, Address, shaderHeader.NumInst);

	SetVertexShaderFromSlots();
}

void CxbxImpl_LoadVertexShader(DWORD Handle, DWORD Address)
{
	// Handle is always address of an X_D3DVertexShader struct, thus always or-ed with 1 (X_D3DFVF_RESERVED0)
	// Address is the slot (offset) from which the program must be written onwards (as whole DWORDS)
	// D3DDevice_LoadVertexShader pushes the program contained in the Xbox VertexShader struct to the NV2A

	xbox::X_D3DVertexShader* pXboxVertexShader = VshHandleToXboxVertexShader(Handle);

	auto pNV2ATokens = &pXboxVertexShader->FunctionData[0];
	DWORD NrTokens = pXboxVertexShader->TotalSize;

#if 1 // TODO : Remove dirty hack (?once CreateVertexShader trampolines to Xbox code that sets TotalSize correctly?) :
	if (NrTokens == 0)
		NrTokens = 10000;
#endif

	unsigned ConstantAddress = 0;
	DWORD* pEnd = pNV2ATokens + NrTokens;
	while (pNV2ATokens < pEnd) {
		DWORD dwMethod, dwSubChannel, nrDWORDS;
		D3DPUSH_DECODE(*pNV2ATokens++, dwMethod, dwSubChannel, nrDWORDS);
		if (nrDWORDS == 0) { LOG_TEST_CASE("Zero-length NV2A method detected!"); break; }
		switch (dwMethod) {
		case NV2A_VP_UPLOAD_INST(0): { // = 0x00000B00
			if ((nrDWORDS & 3) != 0) LOG_TEST_CASE("NV2A_VP_UPLOAD_INST arguments should be a multiple of 4!");
			unsigned nrSlots = nrDWORDS / X_VSH_INSTRUCTION_SIZE;
			CxbxSetVertexShaderSlots(pNV2ATokens, Address, nrSlots);
			Address += nrSlots;
			break;
		}
		case NV2A_VP_UPLOAD_CONST_ID: // = 0x00001EA4
			if (nrDWORDS != 1) LOG_TEST_CASE("NV2A_VP_UPLOAD_CONST_ID should have one argument!");
			ConstantAddress = *pNV2ATokens;
			break;
		case NV2A_VP_UPLOAD_CONST(0): { // = 0x00000B80
			if ((nrDWORDS & 3) != 0) LOG_TEST_CASE("NV2A_VP_UPLOAD_CONST arguments should be a multiple of 4!");
			unsigned nrConstants = nrDWORDS / X_VSH_INSTRUCTION_SIZE;
			// TODO : FIXME : Implement and call SetVertexShaderConstants(pNV2ATokens, ConstantAddress, nrConstants);
			ConstantAddress += nrConstants;
			break;
		}
		default:
			// TODO : Remove this break-out hack once NrTokens is reliable and instead have: DEFAULT_UNREACHABLE;
			pEnd = pNV2ATokens;
			break;
		}

		pNV2ATokens += nrDWORDS;
	}

	SetVertexShaderFromSlots();
}

void CxbxImpl_SetVertexShader(DWORD Handle)
{
	LOG_INIT // Allows use of DEBUG_D3DRESULT

	// Checks if the Handle has bit 0 set - if not, it's a FVF
	// which is converted to a global Xbox Vertex Shader struct
	// Otherwise bit 0 is cleared and the resulting address is
	// validated to be a valid Xbox Vertex Shader
	// D3D state fields are updated.
	// If the shader contains a program, the handle is passed to
	// D3DDevice_LoadVertexShader and D3DDevice_SelectVertexShader.
	// Otherwise the shader is send using push buffer commands.

	HRESULT hRet = D3D_OK;

	g_Xbox_VertexShader_Handle = Handle;

	xbox::X_D3DVertexShader* pXboxVertexShader;
	if (VshHandleIsVertexShader(Handle)) {
		pXboxVertexShader = VshHandleToXboxVertexShader(Handle);
		if (pXboxVertexShader->Flags & X_VERTEXSHADER_FLAG_PROGRAM) {
#if 0 // Since the D3DDevice_SetVertexShader patch already called it's trampoline, these calls have already been executed :
			CxbxImpl_LoadVertexShader(Handle, 0);
			CxbxImpl_SelectVertexShader(Handle, 0);
#endif
		}
		else {
			if (pXboxVertexShader->Flags & X_VERTEXSHADER_FLAG_PASSTHROUGH) {
				LOG_TEST_CASE("TODO : Select Pass-through program HLSL Shader");
			}
			else {
				LOG_TEST_CASE("TODO : Select Fixed Function HLSL Shader");
			}
		}
	}
	else {
		hRet = g_pD3DDevice->SetVertexShader(nullptr);
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShader");
		// TODO : Avoid SetFVF by using XboxFVFToXboxVertexAttributeFormat(Handle);
		hRet = g_pD3DDevice->SetFVF(Handle);
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetFVF");
	}
}

void CxbxImpl_DeleteVertexShader(DWORD Handle)
{
	LOG_INIT // Allows use of DEBUG_D3DRESULT

	// Handle is always address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// It's reference count is lowered. If it reaches zero (0), the struct is freed.

	if (VshHandleIsVertexShader(Handle))
	{
		CxbxVertexShader* pCxbxVertexShader = GetCxbxVertexShader(Handle); // Fetch from cache
		if (pCxbxVertexShader == nullptr) {
			return; // Avoid crash if no shader was cached yet
		}

		RegisterCxbxVertexShader(Handle, nullptr); // Remove from cache

		if (pCxbxVertexShader->Declaration.pHostVertexDeclaration) {
			HRESULT hRet = pCxbxVertexShader->Declaration.pHostVertexDeclaration->Release();
			DEBUG_D3DRESULT(hRet, "g_pD3DDevice->DeleteVertexShader(pHostVertexDeclaration)");
			pCxbxVertexShader->Declaration.pHostVertexDeclaration = nullptr;
		}

		// Release the host vertex shader
		g_VertexShaderSource.ReleaseShader(pCxbxVertexShader->VertexShaderKey);

		FreeVertexDynamicPatch(pCxbxVertexShader);

		free(pCxbxVertexShader);
	}
}

void CxbxImpl_SetVertexShaderConstant(INT Register, PVOID pConstantData, DWORD ConstantCount)
{
	LOG_INIT // Allows use of DEBUG_D3DRESULT

/*#ifdef _DEBUG_TRACK_VS_CONST
	for (uint32_t i = 0; i < ConstantCount; i++)
	{
		printf("SetVertexShaderConstant, c%d  = { %f, %f, %f, %f }\n",
			   Register + i,
			   *((float*)pConstantData + 4 * i),
			   *((float*)pConstantData + 4 * i + 1),
			   *((float*)pConstantData + 4 * i + 2),
			   *((float*)pConstantData + 4 * i + 3));
	}
#endif*/ // _DEBUG_TRACK_VS_CONST

// Xbox vertex shader constants range from -96 to 95
// The host does not support negative, so we adjust to 0..191
	Register += X_D3DSCM_CORRECTION;

	if (Register < 0) LOG_TEST_CASE("Register < 0");
	if (Register + ConstantCount > X_D3DVS_CONSTREG_COUNT) LOG_TEST_CASE("Register + ConstantCount > X_D3DVS_CONSTREG_COUNT");
	HRESULT hRet;
	hRet = g_pD3DDevice->SetVertexShaderConstantF(
		Register,
		(float*)pConstantData,
		ConstantCount
	);
	DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShaderConstant");

	if (FAILED(hRet))
	{
		EmuLog(LOG_LEVEL::WARNING, "We're lying about setting a vertex shader constant!");
		hRet = D3D_OK;
	}
}

// parse xbox vertex shader function into an intermediate format
extern void EmuParseVshFunction
(
	DWORD* pXboxFunction,
	DWORD* pXboxFunctionSize,
	IntermediateVertexShader* pShader
)
{
	auto VshDecoder = XboxVertexShaderDecoder();

	*pXboxFunctionSize = 0;

	// FIXME tidy handling of the header vs headerless cases
	// Normally, pXboxFunction has a shader header before the shader tokens
	// But we can also load shader tokens directly from the Xbox vertex shader slots too

	bool headerless = pXboxFunction[0] == 0; // if its a token instead of a header, first DWORD is unused
	auto headerSize = headerless ? 0 : sizeof(xbox::X_VSH_SHADER_HEADER);

	// Decode the vertex shader program tokens into an intermediate representation
	uint32_t* pCurToken = (uint32_t*)((uintptr_t)pXboxFunction + headerSize);

	if (headerless) {
		// We've been fed shader slots. Make up a header...
		pShader->Header.Version = VERSION_XVS;
		pShader->Header.NumInst = (uint16_t)pShader->Instructions.size();

		// Decode until we hit a token marked final
		while (VshDecoder.VshConvertToIntermediate(pCurToken, pShader)) {
			pCurToken += X_VSH_INSTRUCTION_SIZE;
		}
	}
	else {
		pShader->Header = *(xbox::X_VSH_SHADER_HEADER*)pXboxFunction;
		// Decode only up to the number of instructions in the header
		// The last instruction may not be marked final:
		// Test case: Multiple Vertex Shaders sample
		for (int i = 0; i < pShader->Header.NumInst; i++) {
			if (!VshDecoder.VshConvertToIntermediate(pCurToken, pShader)) {
				if (i < pShader->Header.NumInst - 1) {
					LOG_TEST_CASE("Shader instructions after final instruction");
				}
				break;
			}
			pCurToken += X_VSH_INSTRUCTION_SIZE;
		}
	}

	// The size of the shader is
	pCurToken += X_VSH_INSTRUCTION_SIZE; // always at least one token
	*pXboxFunctionSize = (intptr_t)pCurToken - (intptr_t)pXboxFunction;
}
