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

#define _DEBUG_TRACK_VS

// prevent name collisions
namespace xboxkrnl
{
#include <xboxkrnl/xboxkrnl.h>
};

#include "core\kernel\init\CxbxKrnl.h"
#include "core\kernel\memory-manager\VMManager.h" // for g_VMManager
#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\Direct3D9\Direct3D9.h" // For g_Xbox_VertexShader_Handle
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbPushBuffer.h" // For PUSH_METHOD, PUSH_COUNT
#include "core\hle\Intercept.hpp" // For XB_trampoline
#include "devices\video\nv2a_regs.h" // For NV097_SET_TRANSFORM_PROGRAM et al

#include <sstream>
#include <unordered_map>
#include <array>
#include <bitset>

// External symbols :
extern XTL::X_STREAMINPUT g_Xbox_SetStreamSource[X_VSH_MAX_STREAMS]; // Declared in XbVertexBuffer.cpp
extern XTL::X_VERTEXSHADERCONSTANTMODE g_Xbox_VertexShaderConstantMode; // Declared in Direct3D9.cpp
void* GetDataFromXboxResource(XTL::X_D3DResource* pXboxResource); // Implemented in Direct3D9.cpp

XboxVertexShaderConverter XboxVertexShaders;

// Variables set by [D3DDevice|CxbxImpl]_SetVertexShaderInput() :
                    unsigned g_Xbox_SetVertexShaderInput_Count = 0; // Read by GetXboxVertexAttributes
          XTL::X_STREAMINPUT g_Xbox_SetVertexShaderInput_Data[X_VSH_MAX_STREAMS] = { 0 }; // Active when g_Xbox_SetVertexShaderInput_Count > 0
XTL::X_VERTEXATTRIBUTEFORMAT g_Xbox_SetVertexShaderInput_Attributes = { 0 }; // Read by GetXboxVertexAttributes when g_Xbox_SetVertexShaderInput_Count > 0

// Variables set by [D3DDevice|CxbxImpl]_SetVertexShader() and [D3DDevice|CxbxImpl]_SelectVertexShader() :
                  XTL::DWORD g_Xbox_VertexShader_Handle = 0;
                  XTL::DWORD g_Xbox_VertexShader_FunctionSlots_StartAddress = 0;

// Variable set by [D3DDevice|CxbxImpl]_LoadVertexShader() / [D3DDevice|CxbxImpl]_LoadVertexShaderProgram() (both through CxbxCopyVertexShaderFunctionSlots):
				  XTL::DWORD g_Xbox_VertexShader_FunctionSlots[X_VSH_MAX_INSTRUCTION_COUNT * X_VSH_INSTRUCTION_SIZE]; // Each slot takes either 4 DWORDS (for instructions) or 4 floats (for constants)

typedef uint16_t binary16_t; // Quick and dirty way to indicate IEEE754-2008 'half-precision floats'

bool XboxVertexShaderConverter::Init()
{
	for (int i = 0; i < X_VSH_MAX_INSTRUCTION_COUNT; i++)
		g_Xbox_VertexShader_FunctionSlots[(i * X_VSH_INSTRUCTION_SIZE) + 3] = 1; // Set FLD_FINAL bit on all slots

	// Symbol IDs MUST correspond to how they're registered in XbSymbolDatabase!
	static const std::string D3DDeviceStr = "D3DDEVICE";
	static const std::string VertexShaderStr = "D3DDevice__m_VertexShader_OFFSET";

	// Set g_Xbox_D3DDevice to point to the Xbox D3D Device
	auto it = g_SymbolAddresses.find(D3DDeviceStr);
	if (it != g_SymbolAddresses.end()) {
		g_Xbox_D3DDevice = (DWORD*)it->second;
	} else {
		LOG_TEST_CASE("Couldn't locate D3DDEVICE!");
		return false;
	}

	// Set g_Xbox_D3DDevice to point to the Xbox D3D Device
	it = g_SymbolAddresses.find(VertexShaderStr);
	if (it != g_SymbolAddresses.end()) {
		xbaddr XREF_OFFSET_D3DDEVICE_M_VERTEXSHADER = it->second;
		g_XboxAddr_pVertexShader = (DWORD*)((intptr_t)*g_Xbox_D3DDevice + XREF_OFFSET_D3DDEVICE_M_VERTEXSHADER);
	} else {
		LOG_TEST_CASE("Couldn't locate VertexShader!");
		return false;
	}

	return true;
}

// Reads the active Xbox stream input values (containing VertexBuffer, Offset and Stride) for the given stream number.
// (These values are set through SetStreamSource and can be overridden by SetVertexShaderInput.)
XTL::X_STREAMINPUT& GetXboxVertexStreamInput(unsigned StreamNumber)
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

	enum VSH_SWIZZLE {
		SWIZZLE_X = 0,
		SWIZZLE_Y,
		SWIZZLE_Z,
		SWIZZLE_W
	};

	#define MASK_X 0x008
	#define MASK_Y 0x004
	#define MASK_Z 0x002
	#define MASK_W 0x001

	enum VSH_OREG_NAME {
		OREG_OPOS,    //  0
		OREG_UNUSED1, //  1
		OREG_UNUSED2, //  2
		OREG_OD0,     //  3
		OREG_OD1,     //  4
		OREG_OFOG,    //  5
		OREG_OPTS,    //  6
		OREG_OB0,     //  7
		OREG_OB1,     //  8
		OREG_OT0,     //  9
		OREG_OT1,     // 10
		OREG_OT2,     // 11
		OREG_OT3,     // 12
		OREG_UNUSED3, // 13
		OREG_UNUSED4, // 14
		OREG_A0X      // 15 - all values of the 4 bits are used
	};

	enum VSH_OUTPUT_TYPE {
		OUTPUT_C = 0,
		OUTPUT_O
	};

	enum VSH_PARAMETER_TYPE {
		PARAM_UNKNOWN = 0,
		PARAM_R,          // Temporary (scRatch) registers
		PARAM_V,          // Vertex registers
		PARAM_C,          // Constant registers, set by SetVertexShaderConstant
		PARAM_O // = 0??
	};

	enum VSH_OUTPUT_MUX {
		OMUX_MAC = 0,
		OMUX_ILU
	};

	enum VSH_ILU { // Dxbx note : ILU stands for 'Inverse Logic Unit' opcodes
		ILU_NOP = 0,
		ILU_MOV,
		ILU_RCP,
		ILU_RCC,
		ILU_RSQ,
		ILU_EXP,
		ILU_LOG,
		ILU_LIT // = 7 - all values of the 3 bits are used
	};

	enum VSH_MAC { // Dxbx note : MAC stands for 'Multiply And Accumulate' opcodes
		MAC_NOP = 0,
		MAC_MOV,
		MAC_MUL,
		MAC_ADD,
		MAC_MAD,
		MAC_DP3,
		MAC_DPH,
		MAC_DP4,
		MAC_DST,
		MAC_MIN,
		MAC_MAX,
		MAC_SLT,
		MAC_SGE,
		MAC_ARL
		// ??? 14
		// ??? 15 - 2 values of the 4 bits are undefined
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

	enum VSH_IMD_OUTPUT_TYPE {
		IMD_OUTPUT_C,
		IMD_OUTPUT_R,
		IMD_OUTPUT_O,
		IMD_OUTPUT_A0X
	} ;

	typedef struct _VSH_IMD_OUTPUT {
		VSH_IMD_OUTPUT_TYPE Type;
		int16_t             Address;
		int8_t              Mask;
	} VSH_IMD_OUTPUT;

	typedef struct _VSH_IMD_PARAMETER {
		VSH_PARAMETER_TYPE  ParameterType;   // Parameter type, R, V or C
		bool                Neg;             // true if negated, false if not
		VSH_SWIZZLE         Swizzle[4];      // The four swizzles
		int16_t             Address;         // Register address
	} VSH_IMD_PARAMETER;

	typedef struct _VSH_INTERMEDIATE_FORMAT {
		VSH_MAC                  MAC;
		VSH_ILU                  ILU;
		VSH_IMD_OUTPUT           Output;
		unsigned                 ParamCount;
		VSH_IMD_PARAMETER        Parameters[3];
		// There is only a single address register in Microsoft DirectX 8.0.
		// The address register, designated as a0.x, may be used as signed
		// integer offset in relative addressing into the constant register file.
		//     c[a0.x + n]
		bool                     IndexesWithA0_X;
	} VSH_INTERMEDIATE_FORMAT;

	// State variables :

	uint16_t                 IntermediateCount = 0;
	VSH_INTERMEDIATE_FORMAT  Intermediate[VSH_MAX_INTERMEDIATE_COUNT] = {};

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

		if (IntermediateCount >= VSH_MAX_INTERMEDIATE_COUNT) {
			CxbxKrnlCleanup("Shader exceeds conversion buffer!");
		}

		VSH_INTERMEDIATE_FORMAT* pIntermediate = &(Intermediate[IntermediateCount++]);
		pIntermediate->MAC = MAC;
		pIntermediate->ILU = ILU;
		pIntermediate->Output.Type = output_type;
		pIntermediate->Output.Address = output_address;
		pIntermediate->Output.Mask = output_mask;
		// Get a0.x indirect constant addressing
		pIntermediate->IndexesWithA0_X = VshGetField(pShaderToken, FLD_A0X) > 0; // Applies to PARAM_C parameter reads

		int16_t R;
		int16_t V = VshGetField(pShaderToken, FLD_V);
		int16_t C = ConvertCRegister(VshGetField(pShaderToken, FLD_CONST));
		pIntermediate->ParamCount = 0;
		if (MAC >= MAC_MOV) {
			// Get parameter A
			R = VshGetField(pShaderToken, FLD_A_R);
			VshConvertIntermediateParam(pIntermediate->Parameters[pIntermediate->ParamCount++], pShaderToken, FLD_A_MUX, FLD_A_NEG, R, V, C);
		}

		if ((MAC == MAC_MUL) || ((MAC >= MAC_MAD) && (MAC <= MAC_SGE))) {
			// Get parameter B
			R = VshGetField(pShaderToken, FLD_B_R);
			VshConvertIntermediateParam(pIntermediate->Parameters[pIntermediate->ParamCount++], pShaderToken, FLD_B_MUX, FLD_B_NEG, R, V, C);
		}

		if ((ILU >= ILU_MOV) || (MAC == MAC_ADD) || (MAC == MAC_MAD)) {
			// Get parameter C
			R = VshGetField(pShaderToken, FLD_C_R_HIGH) << 2 | VshGetField(pShaderToken, FLD_C_R_LOW);
			VshConvertIntermediateParam(pIntermediate->Parameters[pIntermediate->ParamCount++], pShaderToken, FLD_C_MUX, FLD_C_NEG, R, V, C);
		}
	}

public:
	bool VshIsEndToken(uint32_t* pShaderToken)
	{
		return VshGetField(pShaderToken, FLD_FINAL) > 0;
	}

	bool VshConvertToIntermediate(uint32_t* pShaderToken)
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
					VshAddIntermediateInstruction(pShaderToken, MAC, ILU_NOP, IMD_OUTPUT_A0X, 0, MASK_X);
				} else {
					VshAddIntermediateInstruction(pShaderToken, MAC, ILU_NOP, IMD_OUTPUT_R, RAddress, VshGetField(pShaderToken, FLD_OUT_MAC_MASK));
				}
			}

			// Check if we must add a muxed MAC opcode as well
			if (OutputMux == OMUX_MAC) {
				VshAddIntermediateInstruction(pShaderToken, MAC, ILU_NOP, OutputType, OutputAddress, VshGetField(pShaderToken, FLD_OUT_O_MASK));
			}
		}

		// Check if there's an ILU opcode
		if (ILU != ILU_NOP) {
			// Paired ILU opcodes will only write to R1
			VshAddIntermediateInstruction(pShaderToken, MAC_NOP, ILU, IMD_OUTPUT_R, bIsPaired ? 1 : RAddress, VshGetField(pShaderToken, FLD_OUT_ILU_MASK));
			// Check if we must add a muxed ILU opcode as well
			if (OutputMux == OMUX_ILU) {
				VshAddIntermediateInstruction(pShaderToken, MAC_NOP, ILU, OutputType, OutputAddress, VshGetField(pShaderToken, FLD_OUT_O_MASK));
			}
		}

		return !VshIsEndToken(pShaderToken);
	}

	// HLSL generation - TODO : Move this to another (friend) class??
private:
	static void OutputHlsl(std::stringstream& hlsl, VSH_IMD_OUTPUT& dest)
	{
		static const char* OReg_Name[/*VSH_OREG_NAME*/] = {
			"oPos",
			"???",
			"???",
			"oD0",
			"oD1",
			"oFog",
			"oPts",
			"oB0",
			"oB1",
			"oT0",
			"oT1",
			"oT2",
			"oT3",
			"???",
			"???",
			"a0.x"
		};

		switch (dest.Type) {
		case IMD_OUTPUT_C:
			// Access the HLSL capital C[] constants array, with the index bias applied :
			// TODO : Avoid out-of-bound writes (perhaps writing to a reserved index?)
			hlsl << "C[" << dest.Address + X_D3DSCM_CORRECTION << "]";
			LOG_TEST_CASE("Vertex shader writes to constant table");
			break;
		case IMD_OUTPUT_R:
			hlsl << "r" << dest.Address;
			break;
		case IMD_OUTPUT_O:
			assert(dest.Address < OREG_A0X);
			hlsl << OReg_Name[dest.Address];
			break;
		case IMD_OUTPUT_A0X:
			hlsl << "a0";
			break;
		default:
			assert(false);
			break;
		}

		// Write the mask as a separate argument to the opcode defines
		// (No space, so that "dest,mask, ..." looks close to "dest.mask, ...")
		hlsl << ",";
		if (dest.Mask & MASK_X) hlsl << "x";
		if (dest.Mask & MASK_Y) hlsl << "y";
		if (dest.Mask & MASK_Z) hlsl << "z";
		if (dest.Mask & MASK_W) hlsl << "w";
	}

	static void ParameterHlsl(std::stringstream& hlsl, VSH_IMD_PARAMETER& param, bool IndexesWithA0_X)
	{
		static char* RegisterName[/*VSH_PARAMETER_TYPE*/] = {
			"?", // PARAM_UNKNOWN = 0,
			"r", // PARAM_R,          // Temporary (scRatch) registers
			"v", // PARAM_V,          // Vertex registers
			"c", // PARAM_C,          // Constant registers, set by SetVertexShaderConstant
			"oPos" // PARAM_O // = 0??
		};

		if (param.Neg) {
			hlsl << "-";
		}

		if (param.ParameterType == PARAM_C) {
			// Access constant registers through our HLSL c() function,
			// which allows dumping negative indices (like Xbox shaders),
			// and which returns zero when out-of-bounds indices are passed in:
			if (IndexesWithA0_X) {
				if (param.Address == 0) {
					hlsl << "c(a0.x)"; // Hide the offset if it's 0
				} else if (param.Address < 0) {
					hlsl << "c(a0.x" << param.Address << ")"; // minus is part of the offset
				} else {
					hlsl << "c(a0.x+" << param.Address << ")"; // show addition character
				}
			} else {
				hlsl << "c(" << param.Address << ")";
			}
		} else {
			hlsl << RegisterName[param.ParameterType] << param.Address;
		}

		// Write the swizzle if we need to
		// Only bother printing the swizzle if it is not the default .xyzw
		if (!(param.Swizzle[0] == SWIZZLE_X &&
			param.Swizzle[1] == SWIZZLE_Y &&
			param.Swizzle[2] == SWIZZLE_Z &&
			param.Swizzle[3] == SWIZZLE_W)) {
			// We'll try to simplify swizzles if we can
			// If all swizzles are the same, we only need to write one out
			unsigned swizzles = 1;

			// Otherwise, we need to use the full swizzle
			if (param.Swizzle[0] != param.Swizzle[1] ||
				param.Swizzle[0] != param.Swizzle[2] ||
				param.Swizzle[0] != param.Swizzle[3]) {
				// Note, we can't remove trailing repeats, like in VS asm,
				// as it may change the type from float4 to float3, float2 or float1!
				swizzles = 4;
			}

			hlsl << ".";
			for (unsigned i = 0; i < swizzles; i++) {
				hlsl << "xyzw"[param.Swizzle[i]];
			}
		}
	}

public:
	bool BuildShader(std::stringstream& hlsl)
	{
		// HLSL strings for all MAC opcodes, indexed with VSH_MAC
		static std::string VSH_MAC_HLSL[/*VSH_MAC*/] = {
			/*MAC_NOP:*/"",
			/*MAC_MOV:*/"x_mov",
			/*MAC_MUL:*/"x_mul",
			/*MAC_ADD:*/"x_add",
			/*MAC_MAD:*/"x_mad",
			/*MAC_DP3:*/"x_dp3",
			/*MAC_DPH:*/"x_dph",
			/*MAC_DP4:*/"x_dp4",
			/*MAC_DST:*/"x_dst",
			/*MAC_MIN:*/"x_min",
			/*MAC_MAX:*/"x_max",
			/*MAC_SLT:*/"x_slt",
			/*MAC_SGE:*/"x_sge",
			/*MAC_ARL:*/"x_arl",
						"",
						"" // VSH_MAC 2 final values of the 4 bits are undefined/unknown  TODO : Investigate their effect (if any) and emulate that as well
		};

		// HLSL strings for all ILU opcodes, indexed with VSH_ILU
		static std::string VSH_ILU_HLSL[/*VSH_ILU*/] = {
			/*ILU_NOP:*/"",
			/*ILU_MOV:*/"x_mov",
			/*ILU_RCP:*/"x_rcp",
			/*ILU_RCC:*/"x_rcc",
			/*ILU_RSQ:*/"x_rsq",
			/*ILU_EXP:*/"x_expp",
			/*ILU_LOG:*/"x_logp",
			/*ILU_LIT:*/"x_lit" // = 7 - all values of the 3 bits are used
		};

		for (int i = 0; i < IntermediateCount; i++) {
			VSH_INTERMEDIATE_FORMAT& IntermediateInstruction = Intermediate[i];

			std::string str;
			if (IntermediateInstruction.MAC > MAC_NOP) {
				str = VSH_MAC_HLSL[IntermediateInstruction.MAC];
			} else {
				str = VSH_ILU_HLSL[IntermediateInstruction.ILU];
			}

			hlsl << "\n  " << str << "("; // opcode
			OutputHlsl(hlsl, IntermediateInstruction.Output);
			for (unsigned i = 0; i < IntermediateInstruction.ParamCount; i++) {
				hlsl << ", ";
				ParameterHlsl(hlsl, IntermediateInstruction.Parameters[i], IntermediateInstruction.IndexesWithA0_X);
			}

			hlsl << ");";
		}

		return IntermediateCount > 0;
	}
};

// ****************************************************************************
// * Vertex shader declaration recompiler
// ****************************************************************************

extern D3DCAPS g_D3DCaps;

std::array<bool, 16> RegVIsPresentInDeclaration; // TODO : Move inside XboxVertexDeclarationConverter again

class XboxVertexDeclarationConverter
{
protected:
	// Internal variables
	CxbxVertexShaderInfo* pVertexShaderInfoToSet;
	CxbxVertexShaderStreamInfo* pCurrentVertexShaderStreamInfo = nullptr;
	bool IsFixedFunction;
	D3DVERTEXELEMENT* pRecompiled;
	//std::array<bool, 16> RegVIsPresentInDeclaration; // TODO : Re-enable again

public:
	// Output
	DWORD XboxDeclarationCount;

private:
	#define D3DDECLUSAGE_UNSUPPORTED ((D3DDECLUSAGE)-1)

	static D3DDECLUSAGE Xb2PCRegisterType
	(
		DWORD VertexRegister,
		BYTE& PCUsageIndex
	)
	{
		D3DDECLUSAGE PCRegisterType;
		PCUsageIndex = 0;

		switch (VertexRegister)
		{
		case XTL::X_D3DVSDE_VERTEX: // -1
			PCRegisterType = D3DDECLUSAGE_UNSUPPORTED;
			break;
		case XTL::X_D3DVSDE_POSITION: // 0
			PCRegisterType = D3DDECLUSAGE_POSITION;
			break;
		case XTL::X_D3DVSDE_BLENDWEIGHT: // 1
			PCRegisterType = D3DDECLUSAGE_BLENDWEIGHT;
			break;
		case XTL::X_D3DVSDE_NORMAL: // 2
			PCRegisterType = D3DDECLUSAGE_NORMAL;
			break;
		case XTL::X_D3DVSDE_DIFFUSE: // 3
			PCRegisterType = D3DDECLUSAGE_COLOR; PCUsageIndex = 0;
			break;
		case XTL::X_D3DVSDE_SPECULAR: // 4
			PCRegisterType = D3DDECLUSAGE_COLOR; PCUsageIndex = 1;
			break;
		case XTL::X_D3DVSDE_FOG: // 5
			PCRegisterType = D3DDECLUSAGE_FOG;
			break;
		case XTL::X_D3DVSDE_POINTSIZE: // 6
			PCRegisterType = D3DDECLUSAGE_PSIZE;
			break;
		case XTL::X_D3DVSDE_BACKDIFFUSE: // 7
			PCRegisterType = D3DDECLUSAGE_COLOR; PCUsageIndex = 2;
			break;
		case XTL::X_D3DVSDE_BACKSPECULAR: // 8
			PCRegisterType = D3DDECLUSAGE_COLOR; PCUsageIndex = 3;
			break;
		case XTL::X_D3DVSDE_TEXCOORD0: // 9
			PCRegisterType = D3DDECLUSAGE_TEXCOORD; PCUsageIndex = 0;
			break;
		case XTL::X_D3DVSDE_TEXCOORD1: // 10
			PCRegisterType = D3DDECLUSAGE_TEXCOORD; PCUsageIndex = 1;
			break;
		case XTL::X_D3DVSDE_TEXCOORD2: // 11
			PCRegisterType = D3DDECLUSAGE_TEXCOORD; PCUsageIndex = 2;
			break;
		case XTL::X_D3DVSDE_TEXCOORD3: // 12
			PCRegisterType = D3DDECLUSAGE_TEXCOORD; PCUsageIndex = 3;
			break;
		default:
			PCRegisterType = D3DDECLUSAGE_UNSUPPORTED;
			break;
		}

		return PCRegisterType;
	}

	static char* XboxVertexRegisterAsString(DWORD VertexRegister)
	{
		switch (VertexRegister)
		{
		case XTL::X_D3DVSDE_VERTEX: // -1
			return "D3DVSDE_VERTEX /* xbox ext. */";
		case XTL::X_D3DVSDE_POSITION: // 0
			return "D3DVSDE_POSITION";
		case XTL::X_D3DVSDE_BLENDWEIGHT: // 1
			return "D3DVSDE_BLENDWEIGHT";
		case XTL::X_D3DVSDE_NORMAL: // 2
			return "D3DVSDE_NORMAL";
		case XTL::X_D3DVSDE_DIFFUSE: // 3
			return "D3DVSDE_DIFFUSE";
		case XTL::X_D3DVSDE_SPECULAR: // 4
			return "D3DVSDE_SPECULAR";
		case XTL::X_D3DVSDE_FOG: // 5
			return "D3DVSDE_FOG";
		case XTL::X_D3DVSDE_POINTSIZE: // 6
			return "D3DVDSE_POINTSIZE";
		case XTL::X_D3DVSDE_BACKDIFFUSE: // 7
			return "D3DVSDE_BACKDIFFUSE /* xbox ext. */";
		case XTL::X_D3DVSDE_BACKSPECULAR: // 8
			return "D3DVSDE_BACKSPECULAR /* xbox ext. */";
		case XTL::X_D3DVSDE_TEXCOORD0: // 9
			return "D3DVSDE_TEXCOORD0";
		case XTL::X_D3DVSDE_TEXCOORD1: // 10
			return "D3DVSDE_TEXCOORD1";
		case XTL::X_D3DVSDE_TEXCOORD2: // 11
			return "D3DVSDE_TEXCOORD2";
		case XTL::X_D3DVSDE_TEXCOORD3: // 12
			return "D3DVSDE_TEXCOORD3";
		case 13:
			return "13 /* unknown register */";
		case 14:
			return "14 /* unknown register */";
		case 15:
			return "15 /* unknown register */";
		default:
			return "16 /* or higher, unknown register */";
		}
	}

	// VERTEX SHADER

	static DWORD VshGetDeclarationCount(DWORD *pXboxDeclaration)
	{
		DWORD Pos = 0;
		while (pXboxDeclaration[Pos] != X_D3DVSD_END())
		{
			Pos++;
		}

		return Pos + 1;
	}

	static inline DWORD VshGetTokenType(DWORD XboxToken)
	{
		return (XboxToken & X_D3DVSD_TOKENTYPEMASK) >> X_D3DVSD_TOKENTYPESHIFT;
	}

	static inline WORD VshGetVertexStream(DWORD XboxToken)
	{
		return (XboxToken & X_D3DVSD_STREAMNUMBERMASK) >> X_D3DVSD_STREAMNUMBERSHIFT;
	}

	static inline DWORD VshGetVertexRegister(DWORD XboxToken)
	{
		DWORD regNum = (XboxToken & X_D3DVSD_VERTEXREGMASK) >> X_D3DVSD_VERTEXREGSHIFT;
		return regNum;
	}

	static inline DWORD VshGetVertexRegisterIn(DWORD XboxToken)
	{
		DWORD regNum = (XboxToken & X_D3DVSD_VERTEXREGINMASK) >> X_D3DVSD_VERTEXREGINSHIFT;
		return regNum;
	}

	void VshDumpXboxDeclaration(DWORD* pXboxDeclaration)
	{
		DbgVshPrintf("DWORD dwVSHDecl[] =\n{\n");
		unsigned iNumberOfVertexStreams = 0;
		bool bStreamNeedsPatching = false;
		auto pXboxToken = pXboxDeclaration;
		while (*pXboxToken != X_D3DVSD_END()) // X_D3DVSD_TOKEN_END 
		{
			DWORD Step = 1;

			switch (VshGetTokenType(*pXboxToken)) {
			case XTL::X_D3DVSD_TOKEN_NOP: {
				DbgVshPrintf("\tD3DVSD_NOP(),\n");
				break;
			}
			case XTL::X_D3DVSD_TOKEN_STREAM: {
				if (*pXboxToken & X_D3DVSD_STREAMTESSMASK) {
					DbgVshPrintf("\tD3DVSD_STREAM_TESS(),\n");
				} else {
					if (iNumberOfVertexStreams > 0) {
						DbgVshPrintf("\t// NeedPatching: %d\n", bStreamNeedsPatching);
					}
					DWORD StreamNumber = VshGetVertexStream(*pXboxToken);
					DbgVshPrintf("\tD3DVSD_STREAM(%u),\n", StreamNumber);
					iNumberOfVertexStreams++;
					bStreamNeedsPatching = false;
				}
				break;
			}
			case XTL::X_D3DVSD_TOKEN_STREAMDATA: {
				if (*pXboxToken & X_D3DVSD_MASK_SKIP) {
 					WORD SkipCount = (*pXboxToken & X_D3DVSD_SKIPCOUNTMASK) >> X_D3DVSD_SKIPCOUNTSHIFT;
					if (*pXboxToken & X_D3DVSD_MASK_SKIPBYTES) {
						DbgVshPrintf("\tD3DVSD_SKIPBYTES(%d), /* xbox ext. */\n", SkipCount);
					} else {
						DbgVshPrintf("\tD3DVSD_SKIP(%d),\n", SkipCount);
					}
				} else {
					DWORD VertexRegister = VshGetVertexRegister(*pXboxToken);
					if (IsFixedFunction) {
						DbgVshPrintf("\t\tD3DVSD_REG(%s, ", XboxVertexRegisterAsString(VertexRegister));
					} else {
						DbgVshPrintf("\t\tD3DVSD_REG(%d, ", (BYTE)VertexRegister);
					}

					DWORD XboxVertexElementDataType = (*pXboxToken & X_D3DVSD_DATATYPEMASK) >> X_D3DVSD_DATATYPESHIFT;
					switch (XboxVertexElementDataType) {
					case XTL::X_D3DVSDT_FLOAT1: // 0x12:
						DbgVshPrintf("D3DVSDT_FLOAT1");
						break;
					case XTL::X_D3DVSDT_FLOAT2: // 0x22:
						DbgVshPrintf("D3DVSDT_FLOAT2");
						break;
					case XTL::X_D3DVSDT_FLOAT3: // 0x32:
						DbgVshPrintf("D3DVSDT_FLOAT3");
						break;
					case XTL::X_D3DVSDT_FLOAT4: // 0x42:
						DbgVshPrintf("D3DVSDT_FLOAT4");
						break;
					case XTL::X_D3DVSDT_D3DCOLOR: // 0x40:
						DbgVshPrintf("D3DVSDT_D3DCOLOR");
						break;
					case XTL::X_D3DVSDT_SHORT2: // 0x25:
						DbgVshPrintf("D3DVSDT_SHORT2");
						break;
					case XTL::X_D3DVSDT_SHORT4: // 0x45:
						DbgVshPrintf("D3DVSDT_SHORT4");
						break;
					case XTL::X_D3DVSDT_NORMSHORT1: // 0x11:
						DbgVshPrintf("D3DVSDT_NORMSHORT1 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_NORMSHORT2: // 0x21:
						if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
							DbgVshPrintf("D3DVSDT_NORMSHORT2");
						} else {
							DbgVshPrintf("D3DVSDT_NORMSHORT2 /* xbox ext. */");
							bStreamNeedsPatching = true;
						}
						break;
					case XTL::X_D3DVSDT_NORMSHORT3: // 0x31:
						DbgVshPrintf("D3DVSDT_NORMSHORT3 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_NORMSHORT4: // 0x41:
						if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
							DbgVshPrintf("D3DVSDT_NORMSHORT4");
							// No need for patching in D3D9
						} else {
							DbgVshPrintf("D3DVSDT_NORMSHORT4 /* xbox ext. */");
							bStreamNeedsPatching = true;
						}
						break;
					case XTL::X_D3DVSDT_NORMPACKED3: // 0x16:
						DbgVshPrintf("D3DVSDT_NORMPACKED3 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_SHORT1: // 0x15:
						DbgVshPrintf("D3DVSDT_SHORT1 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_SHORT3: // 0x35:
						DbgVshPrintf("D3DVSDT_SHORT3 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_PBYTE1: // 0x14:
						DbgVshPrintf("D3DVSDT_PBYTE1 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_PBYTE2: // 0x24:
						DbgVshPrintf("D3DVSDT_PBYTE2 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_PBYTE3: // 0x34:
						DbgVshPrintf("D3DVSDT_PBYTE3 /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_PBYTE4: // 0x44:
						if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
							DbgVshPrintf("D3DVSDT_PBYTE4");
						} else {
							DbgVshPrintf("D3DVSDT_PBYTE4 /* xbox ext. */");
							bStreamNeedsPatching = true;
						}
						break;
					case XTL::X_D3DVSDT_FLOAT2H: // 0x72:
						DbgVshPrintf("D3DVSDT_FLOAT2H /* xbox ext. */");
						bStreamNeedsPatching = true;
						break;
					case XTL::X_D3DVSDT_NONE: // 0x02:
						DbgVshPrintf("D3DVSDT_NONE /* xbox ext. */");
						break;
					default:
						DbgVshPrintf("Unknown data type for D3DVSD_REG: 0x%02X\n", XboxVertexElementDataType);
						break;
					}

					DbgVshPrintf("),\n");
				};
				break;
			}
			case XTL::X_D3DVSD_TOKEN_TESSELLATOR: {
				DWORD VertexRegisterOut = VshGetVertexRegister(*pXboxToken);
				if (*pXboxToken & X_D3DVSD_MASK_TESSUV) {
					DbgVshPrintf("\tD3DVSD_TESSUV(%s),\n", XboxVertexRegisterAsString(VertexRegisterOut));
				} else { // D3DVSD_TESSNORMAL
					DWORD VertexRegisterIn = VshGetVertexRegisterIn(*pXboxToken);
					DbgVshPrintf("\tD3DVSD_TESSNORMAL(%s, %s),\n",
						XboxVertexRegisterAsString(VertexRegisterIn),
						XboxVertexRegisterAsString(VertexRegisterOut));
				}
				break;
			}
			case XTL::X_D3DVSD_TOKEN_CONSTMEM: {
				DWORD ConstantAddress = (*pXboxToken & X_D3DVSD_CONSTADDRESSMASK) >> X_D3DVSD_CONSTADDRESSSHIFT;
				DWORD Count = (*pXboxToken & X_D3DVSD_CONSTCOUNTMASK) >> X_D3DVSD_CONSTCOUNTSHIFT;
				DbgVshPrintf("\tD3DVSD_CONST(%d, %d),\n", ConstantAddress, Count);
				LOG_TEST_CASE("X_D3DVSD_TOKEN_CONSTMEM");
				Step = Count * 4 + 1;
				break;
			}
			case XTL::X_D3DVSD_TOKEN_EXT: {
				DWORD ExtInfo = (*pXboxToken & X_D3DVSD_EXTINFOMASK) >> X_D3DVSD_EXTINFOSHIFT;
				DWORD Count = (*pXboxToken & X_D3DVSD_EXTCOUNTMASK) >> X_D3DVSD_EXTCOUNTSHIFT;
				DbgVshPrintf("\tD3DVSD_EXT(%d, %d),\n", ExtInfo, Count);
				LOG_TEST_CASE("X_D3DVSD_TOKEN_EXT");
				Step = Count * 4 + 1; // TODO : Is this correct?
				break;
			}
			default:
				DbgVshPrintf("Unknown token type: %d\n", VshGetTokenType(*pXboxToken));
				break;
			}

			pXboxToken += Step;
		}

		if (iNumberOfVertexStreams > 0) {
			DbgVshPrintf("\t// NeedPatching: %d\n", bStreamNeedsPatching);
		}

		DbgVshPrintf("\tD3DVSD_END()\n};\n");

		DbgVshPrintf("// NbrStreams: %d\n", iNumberOfVertexStreams);
	}

	static void VshConvertToken_NOP(DWORD *pXboxToken)
	{
		if(*pXboxToken != X_D3DVSD_NOP())
		{
			LOG_TEST_CASE("Token NOP found, but extra parameters are given!");
		}
	}

	static DWORD VshConvertToken_CONSTMEM(DWORD *pXboxToken)
	{
		// DWORD ConstantAddress = (*pXboxToken & X_D3DVSD_CONSTADDRESSMASK) >> X_D3DVSD_CONSTADDRESSSHIFT;
		DWORD Count           = (*pXboxToken & X_D3DVSD_CONSTCOUNTMASK) >> X_D3DVSD_CONSTCOUNTSHIFT;
		LOG_TEST_CASE("CONST"); // TODO : Implement
		return Count * 4 + 1;
	}

	void VshConvertToken_TESSELATOR(DWORD *pXboxToken)
	{
		BYTE Index;

		if(*pXboxToken & X_D3DVSD_MASK_TESSUV)
		{
			DWORD VertexRegister    = VshGetVertexRegister(*pXboxToken);
			DWORD NewVertexRegister = VertexRegister;

			NewVertexRegister = Xb2PCRegisterType(VertexRegister, Index);
			// TODO : Expand on the setting of this TESSUV register element :
			pRecompiled->Method = D3DDECLMETHOD_UV; // TODO : Is this correct?
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
			pRecompiled->Method = 0; // TODO ?
			pRecompiled->Usage = D3DDECLUSAGE(NewVertexRegisterIn);
			pRecompiled->UsageIndex = Index;

			NewVertexRegisterOut = Xb2PCRegisterType(VertexRegisterOut, Index);
			// TODO : Expand on the setting of this TESSNORMAL output register element :
			pRecompiled++;
			pRecompiled->Method = D3DDECLMETHOD_CROSSUV; // TODO : Is this correct?
			pRecompiled->Usage = D3DDECLUSAGE(NewVertexRegisterOut);
			pRecompiled->UsageIndex = Index;
		}
	}

	void VshConvertToken_STREAM(DWORD *pXboxToken)
	{
		// D3DVSD_STREAM_TESS
		if(*pXboxToken & X_D3DVSD_STREAMTESSMASK)
		{
			// TODO
		}
		else // D3DVSD_STREAM
		{
			DWORD StreamNumber = VshGetVertexStream(*pXboxToken);

			// new stream
			pCurrentVertexShaderStreamInfo = &(pVertexShaderInfoToSet->VertexStreams[StreamNumber]);
			pCurrentVertexShaderStreamInfo->NeedPatch = FALSE;
			pCurrentVertexShaderStreamInfo->DeclPosition = FALSE;
			pCurrentVertexShaderStreamInfo->CurrentStreamNumber = 0;
			pCurrentVertexShaderStreamInfo->HostVertexStride = 0;
			pCurrentVertexShaderStreamInfo->NumberOfVertexElements = 0;

			// Dxbx note : Use Dophin(s), FieldRender, MatrixPaletteSkinning and PersistDisplay as a testcase

			pCurrentVertexShaderStreamInfo->CurrentStreamNumber = VshGetVertexStream(*pXboxToken);
			pVertexShaderInfoToSet->NumberOfVertexStreams++;
			// TODO : Keep a bitmask for all StreamNumber's seen?
		}
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

	void VshConvert_SkipBytes(int SkipBytesCount)
	{
		if (SkipBytesCount % sizeof(DWORD)) {
			LOG_TEST_CASE("D3DVSD_SKIPBYTES not divisble by 4!");
		}
#if 0 // Potential optimization, for now disabled for simplicity :
		else {
			// Skip size is a whole multiple of 4 bytes;
			// Is stream patching not needed up until this element?
			if (!pCurrentVertexShaderStreamInfo->NeedPatch) {
				// Then we can get away with increasing the host stride,
				// which avoids otherwise needless vertex buffer patching :
				pCurrentVertexShaderStreamInfo->HostVertexStride += SkipBytesCount;
				return;
			}
		}
#endif

		// Register a 'skip' element, so that Xbox data will be skipped
		// without increasing host stride - this does require patching :
		VshConvert_RegisterVertexElement(XTL::X_D3DVSDT_NONE, SkipBytesCount, /*HostSize=*/0, /*NeedPatching=*/TRUE);
	}

	void VshConvertToken_STREAMDATA_SKIP(DWORD *pXboxToken)
	{
		WORD SkipCount = (*pXboxToken & X_D3DVSD_SKIPCOUNTMASK) >> X_D3DVSD_SKIPCOUNTSHIFT;
		VshConvert_SkipBytes(SkipCount * sizeof(DWORD));
	}

	void VshConvertToken_STREAMDATA_SKIPBYTES(DWORD* pXboxToken)
	{
		WORD SkipBytesCount = (*pXboxToken & X_D3DVSD_SKIPCOUNTMASK) >> X_D3DVSD_SKIPCOUNTSHIFT;
		VshConvert_SkipBytes(SkipBytesCount);
	}

	void VshConvertToken_STREAMDATA_REG(DWORD *pXboxToken)
	{
		DWORD VertexRegister = VshGetVertexRegister(*pXboxToken);
		BOOL NeedPatching = FALSE;
		BYTE Index;
		BYTE HostVertexRegisterType;

		if (IsFixedFunction) {
			HostVertexRegisterType = Xb2PCRegisterType(VertexRegister, Index);
		} else {
			// D3DDECLUSAGE_TEXCOORD can be useds for any user-defined data
			// We need this because there is no reliable way to detect the real usage
			// Xbox has no concept of 'usage types', it only requires a list of attribute register numbers.
			// So we treat them all as 'user-defined' with an Index of the Vertex Register Index
			// this prevents information loss in shaders due to non-matching dcl types!
			HostVertexRegisterType = D3DDECLUSAGE_TEXCOORD;
			Index = (BYTE)VertexRegister;
		}

		// Add this register to the list of declared registers
		RegVIsPresentInDeclaration[VertexRegister] = true;

		DWORD XboxVertexElementDataType = (*pXboxToken & X_D3DVSD_DATATYPEMASK) >> X_D3DVSD_DATATYPESHIFT;
		WORD XboxVertexElementByteSize = 0;
		BYTE HostVertexElementDataType = 0;
		WORD HostVertexElementByteSize = 0;

		switch (XboxVertexElementDataType)
		{
		case XTL::X_D3DVSDT_FLOAT1: // 0x12:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT1;
			HostVertexElementByteSize = 1 * sizeof(FLOAT);
			break;
		case XTL::X_D3DVSDT_FLOAT2: // 0x22:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT2;
			HostVertexElementByteSize = 2 * sizeof(FLOAT);
			break;
		case XTL::X_D3DVSDT_FLOAT3: // 0x32:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
			HostVertexElementByteSize = 3 * sizeof(FLOAT);
			break;
		case XTL::X_D3DVSDT_FLOAT4: // 0x42:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
			HostVertexElementByteSize = 4 * sizeof(FLOAT);
			break;
		case XTL::X_D3DVSDT_D3DCOLOR: // 0x40:
			HostVertexElementDataType = D3DDECLTYPE_D3DCOLOR;
			HostVertexElementByteSize = 1 * sizeof(D3DCOLOR);
			break;
		case XTL::X_D3DVSDT_SHORT2: // 0x25:
			HostVertexElementDataType = D3DDECLTYPE_SHORT2;
			HostVertexElementByteSize = 2 * sizeof(SHORT);
			break;
		case XTL::X_D3DVSDT_SHORT4: // 0x45:
			HostVertexElementDataType = D3DDECLTYPE_SHORT4;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			break;
		case XTL::X_D3DVSDT_NORMSHORT1: // 0x11:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT2N;
				HostVertexElementByteSize = 2 * sizeof(SHORT);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT1;
				HostVertexElementByteSize = 1 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 1 * sizeof(XTL::SHORT);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_NORMSHORT2: // 0x21:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT2N;
				HostVertexElementByteSize = 2 * sizeof(SHORT);
				// No need for patching in D3D9
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT2;
				HostVertexElementByteSize = 2 * sizeof(FLOAT);
				XboxVertexElementByteSize = 2 * sizeof(XTL::SHORT);
				NeedPatching = TRUE;
			}
			break;
		case XTL::X_D3DVSDT_NORMSHORT3: // 0x31:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT4N;
				HostVertexElementByteSize = 4 * sizeof(SHORT);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
				HostVertexElementByteSize = 3 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 3 * sizeof(XTL::SHORT);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_NORMSHORT4: // 0x41:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
				HostVertexElementDataType = D3DDECLTYPE_SHORT4N;
				HostVertexElementByteSize = 4 * sizeof(SHORT);
				// No need for patching in D3D9
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
				HostVertexElementByteSize = 4 * sizeof(FLOAT);
				XboxVertexElementByteSize = 4 * sizeof(XTL::SHORT);
				NeedPatching = TRUE;
			}
			break;
		case XTL::X_D3DVSDT_NORMPACKED3: // 0x16:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
			HostVertexElementByteSize = 3 * sizeof(FLOAT);
			XboxVertexElementByteSize = 1 * sizeof(XTL::DWORD);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_SHORT1: // 0x15:
			HostVertexElementDataType = D3DDECLTYPE_SHORT2;
			HostVertexElementByteSize = 2 * sizeof(SHORT);
			XboxVertexElementByteSize = 1 * sizeof(XTL::SHORT);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_SHORT3: // 0x35:
			HostVertexElementDataType = D3DDECLTYPE_SHORT4;
			HostVertexElementByteSize = 4 * sizeof(SHORT);
			XboxVertexElementByteSize = 3 * sizeof(XTL::SHORT);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_PBYTE1: // 0x14:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT1;
				HostVertexElementByteSize = 1 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 1 * sizeof(XTL::BYTE);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_PBYTE2: // 0x24:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT2;
				HostVertexElementByteSize = 2 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 2 * sizeof(XTL::BYTE);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_PBYTE3: // 0x34:
			if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
				HostVertexElementDataType = D3DDECLTYPE_UBYTE4N;
				HostVertexElementByteSize = 4 * sizeof(BYTE);
			}
			else
			{
				HostVertexElementDataType = D3DDECLTYPE_FLOAT3;
				HostVertexElementByteSize = 3 * sizeof(FLOAT);
			}
			XboxVertexElementByteSize = 3 * sizeof(XTL::BYTE);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_PBYTE4: // 0x44:
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
				XboxVertexElementByteSize = 4 * sizeof(XTL::BYTE);
				NeedPatching = TRUE;
			}
			break;
		case XTL::X_D3DVSDT_FLOAT2H: // 0x72:
			HostVertexElementDataType = D3DDECLTYPE_FLOAT4;
			HostVertexElementByteSize = 4 * sizeof(FLOAT);
			XboxVertexElementByteSize = 3 * sizeof(FLOAT);
			NeedPatching = TRUE;
			break;
		case XTL::X_D3DVSDT_NONE: // 0x02:
			// No host element data, so no patching
			break;
		default:
			//LOG_TEST_CASE("Unknown data type for D3DVSD_REG: 0x%02X\n", XboxVertexElementDataType);
			break;
		}

		// On X_D3DVSDT_NONE skip this token
		if (XboxVertexElementDataType == XTL::X_D3DVSDT_NONE)
		{
			// Xbox elements with X_D3DVSDT_NONE have size zero, so there's no need to register those.
			// Note, that for skip tokens, we DO call VshConvert_RegisterVertexElement with a X_D3DVSDT_NONE!
			return;
		}

		// save patching information
		VshConvert_RegisterVertexElement(
			XboxVertexElementDataType,
			NeedPatching ? XboxVertexElementByteSize : HostVertexElementByteSize,
			HostVertexElementByteSize,
			NeedPatching);

		pRecompiled->Stream = pCurrentVertexShaderStreamInfo->CurrentStreamNumber;
		pRecompiled->Offset = pCurrentVertexShaderStreamInfo->HostVertexStride;
		pRecompiled->Type = HostVertexElementDataType;
		pRecompiled->Method = D3DDECLMETHOD_DEFAULT;
		pRecompiled->Usage = HostVertexRegisterType;
		pRecompiled->UsageIndex = Index;

		pRecompiled++;

		pCurrentVertexShaderStreamInfo->HostVertexStride += HostVertexElementByteSize;
	}

	void VshConvertToken_STREAMDATA(DWORD *pXboxToken)
	{
		if (*pXboxToken & X_D3DVSD_MASK_SKIP)
		{
			// For D3D9, use D3DDECLTYPE_UNUSED ?
			if (*pXboxToken & X_D3DVSD_MASK_SKIPBYTES) {
				VshConvertToken_STREAMDATA_SKIPBYTES(pXboxToken);
			} else {
				VshConvertToken_STREAMDATA_SKIP(pXboxToken);
			}
		}
		else // D3DVSD_REG
		{
			VshConvertToken_STREAMDATA_REG(pXboxToken);
		}
	}

	DWORD VshRecompileToken(DWORD *pXboxToken)
	{
		DWORD Step = 1;

		switch(VshGetTokenType(*pXboxToken))
		{
		case XTL::X_D3DVSD_TOKEN_NOP:
			VshConvertToken_NOP(pXboxToken);
			break;
		case XTL::X_D3DVSD_TOKEN_STREAM:
		{
			VshConvertToken_STREAM(pXboxToken);
			break;
		}
		case XTL::X_D3DVSD_TOKEN_STREAMDATA:
		{
			VshConvertToken_STREAMDATA(pXboxToken);
			break;
		}
		case XTL::X_D3DVSD_TOKEN_TESSELLATOR:
		{
			VshConvertToken_TESSELATOR(pXboxToken);
			break;
		}
		case XTL::X_D3DVSD_TOKEN_CONSTMEM:
		{
			Step = VshConvertToken_CONSTMEM(pXboxToken);
			break;
		}
		default:
			//LOG_TEST_CASE("Unknown token type: %d\n", VshGetTokenType(*pXboxToken));
			break;
		}

		return Step;
	}

	static DWORD* RemoveXboxDeclarationRedefinition(DWORD* pXboxDeclaration)
	{
		// Detect and remove register redefinitions by preprocessing the Xbox Vertex Declaration
		// Test Case: King Kong

		// Find the last token
		DWORD* pXboxToken = pXboxDeclaration;
		while (*pXboxToken != X_D3DVSD_END()){
			pXboxToken++;
		}

		// Operate on a copy of the Xbox declaration, rather than messing with the Xbox's memory
		auto declarationBytes = sizeof(DWORD) * (pXboxToken - pXboxDeclaration + 1);
		auto pXboxDeclarationCopy = (DWORD*)malloc(declarationBytes);
		memcpy(pXboxDeclarationCopy, pXboxDeclaration, declarationBytes);
		pXboxToken = pXboxDeclarationCopy + (pXboxToken - pXboxDeclaration); // Move to end of the copy

		// Remember if we've seen a given output register
		std::bitset<16> seen;

		// We want to keep later definitions, and remove earlier ones
		// Scan back from the end of the declaration, and replace redefinitions with nops
		while (pXboxToken > pXboxDeclarationCopy) {
			auto type = VshGetTokenType(*pXboxToken);
			if (type == XTL::X_D3DVSD_TOKEN_STREAMDATA && !(*pXboxToken & X_D3DVSD_MASK_SKIP) ||
				type == XTL::X_D3DVSD_TOKEN_TESSELLATOR)
			{
				auto outputRegister = VshGetVertexRegister(*pXboxToken);
				if (seen[outputRegister])
				{
					// Blank out tokens for mapped registers
					*pXboxToken = X_D3DVSD_NOP();
					EmuLog(LOG_LEVEL::DEBUG, "Replacing duplicate definition of register %d with D3DVSD_NOP", outputRegister);
				}
				else
				{
					// Mark register as seen
					seen[outputRegister] = true;
				}
			}

			pXboxToken--;
		}

		return pXboxDeclarationCopy;
	}

public:
	D3DVERTEXELEMENT *Convert(DWORD* pXboxDeclaration, bool bIsFixedFunction, CxbxVertexShaderInfo* pCxbxVertexShaderInfo)
	{
		// Get a preprocessed copy of the original Xbox Vertex Declaration
		auto pXboxVertexDeclarationCopy = RemoveXboxDeclarationRedefinition(pXboxDeclaration);

		pVertexShaderInfoToSet = pCxbxVertexShaderInfo;

		IsFixedFunction = bIsFixedFunction;

		RegVIsPresentInDeclaration.fill(false);

		// First of all some info:
		// We have to figure out which flags are set and then
		// we have to patch their params

		// some token values
		// 0xFFFFFFFF - end of the declaration
		// 0x00000000 - nop (means that this value is ignored)

		// Calculate size of declaration
		unsigned XboxDeclarationCount = VshGetDeclarationCount(pXboxVertexDeclarationCopy);
		// For Direct3D9, we need to reserve at least twice the number of elements, as one token can generate two registers (in and out) :
		unsigned HostDeclarationSize = XboxDeclarationCount * sizeof(D3DVERTEXELEMENT) * 2;
	
		D3DVERTEXELEMENT *Result = (D3DVERTEXELEMENT *)calloc(1, HostDeclarationSize);
		pRecompiled = Result;
		uint8_t *pRecompiledBufferOverflow = ((uint8_t*)pRecompiled) + HostDeclarationSize;

		VshDumpXboxDeclaration(pXboxDeclaration);

		auto pXboxToken = pXboxVertexDeclarationCopy;
		while (*pXboxToken != X_D3DVSD_END())
		{
			if ((uint8_t*)pRecompiled >= pRecompiledBufferOverflow) {
				DbgVshPrintf("Detected buffer-overflow, breaking out...\n");
				break;
			}

			DWORD Step = VshRecompileToken(pXboxToken);
			pXboxToken += Step;
		}

		*pRecompiled = D3DDECL_END();

		// Ensure valid ordering of the vertex declaration (http://doc.51windows.net/Directx9_SDK/graphics/programmingguide/gettingstarted/vertexdeclaration/vertexdeclaration.htm)
		// In particular "All vertex elements for a stream must be consecutive and sorted by offset"
		// Test case: King Kong (due to register redefinition)
		std::sort(Result, pRecompiled, [] (const auto& x, const auto& y)
			{ return std::tie(x.Stream, x.Method, x.Offset) < std::tie(y.Stream, y.Method, y.Offset); });

		// Free the preprocessed declaration copy
		free(pXboxVertexDeclarationCopy);

		for (size_t i = 0; i < RegVIsPresentInDeclaration.size(); i++) {
			pCxbxVertexShaderInfo->vRegisterInDeclaration[i] = RegVIsPresentInDeclaration[i];
			EmuLog(LOG_LEVEL::DEBUG, "Vertex regs used: v%d %d", i, pCxbxVertexShaderInfo->vRegisterInDeclaration[i]);
		}

		return Result;
	}
};

D3DVERTEXELEMENT *EmuRecompileVshDeclaration
(
    DWORD                *pXboxDeclaration,
    bool                  bIsFixedFunction,
//    DWORD                *pXboxDeclarationCount,
    CxbxVertexShaderInfo *pCxbxVertexShaderInfo
)
{
	XboxVertexDeclarationConverter Converter;

	D3DVERTEXELEMENT* pHostVertexElements = Converter.Convert(pXboxDeclaration, bIsFixedFunction, pCxbxVertexShaderInfo);

//	*pXboxDeclarationCount = Converter.XboxDeclarationCount;

    return pHostVertexElements;
}

extern void FreeVertexDynamicPatch(CxbxVertexShader *pVertexShader)
{
    pVertexShader->VertexShaderInfo.NumberOfVertexStreams = 0;
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

CxbxVertexShaderInfo *GetCxbxVertexShaderInfo(DWORD XboxVertexShaderHandle)
{
    CxbxVertexShader *pCxbxVertexShader = GetCxbxVertexShader(XboxVertexShaderHandle);
if (pCxbxVertexShader) // TMP : Avoid access violation
    for (uint32_t i = 0; i < pCxbxVertexShader->VertexShaderInfo.NumberOfVertexStreams; i++)
    {
        if (pCxbxVertexShader->VertexShaderInfo.VertexStreams[i].NeedPatch)
        {
            return &pCxbxVertexShader->VertexShaderInfo;
        }
    }
    return nullptr;
}

typedef DWORD CxbxVertexShaderKey_t; // TODO : Hash X_VERTEXATTRIBUTEFORMAT (vertex declaration), X_STREAMINPUT (vertex buffer mapping) and actual vertex buffer contents

// Maps (non-FVF?) Xbox vertex shader handle to CxbxVertexShader (a struct containing a host Xbox vertex shader handle and the original members)
std::unordered_map<CxbxVertexShaderKey_t, CxbxVertexShader*> g_CxbxVertexShaders;

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

void SetCxbxVertexShader(DWORD XboxVertexShaderHandle, CxbxVertexShader* shader)
{
	auto it = g_CxbxVertexShaders.find(XboxVertexShaderHandle);
	if (it != g_CxbxVertexShaders.end() && it->second != nullptr && shader != nullptr) {
		LOG_TEST_CASE("Overwriting existing Vertex Shader");
	}

	g_CxbxVertexShaders[XboxVertexShaderHandle] = shader;
}

XTL::X_D3DVertexShader* GetXboxVertexShader()
{
	using namespace XTL;

	X_D3DVertexShader* pXboxVertexShader = xbnullptr;

	// Only when we're sure of the location of the Xbox Device.m_pVertexShader variable
	if (XboxVertexShaders.g_XboxAddr_pVertexShader) {
		// read that (so that we get access to internal vertex shaders, like those generated
		// to contain the attribute-information for FVF shaders) :
		pXboxVertexShader = (X_D3DVertexShader*)(*XboxVertexShaders.g_XboxAddr_pVertexShader);
	}
	else
	{
		LOG_TEST_CASE("Unknown pVertexShader symbol location!");
		// Otherwise, we have no choice but to use what we've last stored in the
		// g_Xbox_VertexShader_Handle variable via our D3DDevice_SetVertexShader
		// and D3DDevice_SelectVertexShader* patches.

		// Note, that once we have a fail-safe way to determine the location of the
		// Xbox Device.m_pVertexShader symbol, the FVF and the accompanying Address,
		// we no longer need this statement block, nor patches on D3DDevice_SetVertexShader
		// nor D3DDevice_SelectVertexShader* !

		// Now, to convert, we do need to have a valid vertex shader :
		if (g_Xbox_VertexShader_Handle == 0) {
			LOG_TEST_CASE("Unassigned Xbox vertex shader!");
			return nullptr;
		}

		if (!VshHandleIsVertexShader(g_Xbox_VertexShader_Handle)) {
			LOG_TEST_CASE("Xbox vertex shader lacks X_D3DFVF_RESERVED0 bit!");
			return nullptr;
		}

		pXboxVertexShader = VshHandleToXboxVertexShader(g_Xbox_VertexShader_Handle);
	}

	return pXboxVertexShader;
}

XTL::X_VERTEXATTRIBUTEFORMAT* GetXboxVertexAttributes()
{
	XTL::X_D3DVertexShader* pXboxVertexShader = GetXboxVertexShader();
	if (pXboxVertexShader == xbnullptr)
	{
		// Despite possibly not being used, the pXboxVertexShader argument must always be assigned
		LOG_TEST_CASE("Xbox should always have a VertexShader set (even for FVF's)");
		return &g_Xbox_SetVertexShaderInput_Attributes; // WRONG result, but it's already strange this happens
	}

	// If SetVertexShaderInput is active, it's arguments overrule those of the active vertex shader
	if (g_Xbox_SetVertexShaderInput_Count > 0) {
		// Take overrides (on declarations and streaminputs, as optionally set by SetVertexShaderInput) into account :
		return &g_Xbox_SetVertexShaderInput_Attributes;
	}

	return &pXboxVertexShader->VertexAttribute;
}

typedef struct {
	int XboxDeclarationFormat; // X_D3DVSDT_* key
	int XboxNrOfUnits;
	int XboxSizeInBytes;
	//int NrOfDimensions;
	D3DDECLTYPE HostDeclType;
	bool IsCompatible;
	int HostNrOfUnits;
	int HostSizeInBytes;
} XboxVertexAttributeDeclarationDecoded_t;

// VSDT validity check (hand-optimized to avoid the memory accesses that compiled switch-statements would incur)
bool X_D3DVSDT_IsValid(const uint32_t VSDT)
{
	if (VSDT & 0xFFFFFF88) return false; // Bit 31-7 or 3 should never be set

	using namespace XTL;
	// With that out of the way, compact the valid bits (6,5,4,2,1 and 0) efficiently into the least significant 6 bits :
#define _COMPACT_VALID_BITS(x) (((x & 1) << 2) | (x >> 1))
// Define a 64-bit mask containing a set bit for each valid (compacted) VSDT value :
#define _BITMASK(x) ((uint64_t)1 << _COMPACT_VALID_BITS(x))
	static constexpr uint64_t AllValidTypesBitmask = 0
		| _BITMASK(X_D3DVSDT_D3DCOLOR)
		| _BITMASK(X_D3DVSDT_PBYTE1)
		| _BITMASK(X_D3DVSDT_PBYTE2)
		| _BITMASK(X_D3DVSDT_PBYTE3)
		| _BITMASK(X_D3DVSDT_PBYTE4)
		| _BITMASK(X_D3DVSDT_NORMSHORT1)
		| _BITMASK(X_D3DVSDT_NORMSHORT2)
		| _BITMASK(X_D3DVSDT_NORMSHORT3)
		| _BITMASK(X_D3DVSDT_NORMSHORT4)
		| _BITMASK(X_D3DVSDT_SHORT1)
		| _BITMASK(X_D3DVSDT_SHORT2)
		| _BITMASK(X_D3DVSDT_SHORT3)
		| _BITMASK(X_D3DVSDT_SHORT4)
		| _BITMASK(X_D3DVSDT_NONE)
		| _BITMASK(X_D3DVSDT_FLOAT1)
		| _BITMASK(X_D3DVSDT_FLOAT2)
		| _BITMASK(X_D3DVSDT_FLOAT3)
		| _BITMASK(X_D3DVSDT_FLOAT4)
		| _BITMASK(X_D3DVSDT_FLOAT2H)
		| _BITMASK(X_D3DVSDT_NORMPACKED3);

	// Now, return validity of the given VSDT value by checking it's appearance in this bit mask :
	// (Note, shifting the bitmask and returning bit 0 is faster than checking the VSDT'th bit)
	return (AllValidTypesBitmask >> _COMPACT_VALID_BITS(VSDT)) & 1;
#undef _COMPACT_VALID_BITS
#undef _MASK
}

/* Determine size for the following X_D3DVSDT_* types (value) :
	X_D3DVSDT_D3DCOLOR    (0x40) =  4 : sizeof(uint32_t) * 1 (equal to sizeof(uint8_t) * 4)
	X_D3DVSDT_PBYTE1      (0x14) =  1 : sizeof(uint8_t ) * 1
	X_D3DVSDT_PBYTE2      (0x24) =  2 : sizeof(uint8_t ) * 2
	X_D3DVSDT_PBYTE3      (0x34) =  3 : sizeof(uint8_t ) * 3
	X_D3DVSDT_PBYTE4      (0x44) =  4 : sizeof(uint8_t ) * 4

	X_D3DVSDT_NORMSHORT1  (0x11) =  2 : sizeof(uint16_t) * 1
	X_D3DVSDT_NORMSHORT2  (0x21) =  4 : sizeof(uint16_t) * 2
	X_D3DVSDT_NORMSHORT3  (0x31) =  6 : sizeof(uint16_t) * 3
	X_D3DVSDT_NORMSHORT4  (0x41) =  8 : sizeof(uint16_t) * 4
	X_D3DVSDT_SHORT1      (0x15) =  2 : sizeof(uint16_t) * 1
	X_D3DVSDT_SHORT2      (0x25) =  4 : sizeof(uint16_t) * 2
	X_D3DVSDT_SHORT3      (0x35) =  6 : sizeof(uint16_t) * 3
	X_D3DVSDT_SHORT4      (0x45) =  8 : sizeof(uint16_t) * 4

	X_D3DVSDT_NONE        (0x02) =  0 : sizeof(float   ) * 0
	X_D3DVSDT_FLOAT1      (0x12) =  4 : sizeof(float   ) * 1
	X_D3DVSDT_FLOAT2      (0x22) =  8 : sizeof(float   ) * 2
	X_D3DVSDT_FLOAT3      (0x32) = 12 : sizeof(float   ) * 3
	X_D3DVSDT_FLOAT4      (0x42) = 16 : sizeof(float   ) * 4
	X_D3DVSDT_FLOAT2H     (0x72) = 12 : sizeof(float   ) * 3 (NOT 7!)

	X_D3DVSDT_NORMPACKED3 (0x16) =  4 : sizeof(uint32_t) * 1

No other inputs should be given (they will produce unreliable sizes).
*/
unsigned X_D3DVSDT_SizeInBytes(const uint8_t VSDT) // VSDT = VSDGetDataType(dwDecl) a.k.a. SizeAndType
{
	// First handle the case that doesn't fit the code below :
	if (VSDT == XTL::X_D3DVSDT_FLOAT2H) // == 0x72
		return 12; // == 3 * sizeof(float)

	unsigned Shift = VSDT & 0x3;
	// Shift becomes 0 for VSDT 0x?4 (X_D3DVSDT_PBYTE1, etc) and 0x?0 (X_D3DVSDT_D3DCOLOR),
	// Shift becomes 1 for VSDT 0x?5 (X_D3DVSDT_SHORT1, etc) and 0x?1 (X_D3DVSDT_NORMSHORT1, etc),
	// Shift becomes 2 for VSDT 0x?6 (X_D3DVSDT_NORMPACKED3) and 0x?2 (X_D3DVSDT_NONE, X_D3DVSDT_FLOAT1, etc).
	// A Shift of 0 equals factor (1 << 0) == 1 == sizeof(uint8_t)
	// A Shift of 1 equals factor (1 << 1) == 2 == sizeof(uint16_t)
	// A Shift of 2 equals factor (1 << 2) == 4 == sizeof(float) == sizeof(uint32_t)

	unsigned Size = VSDT >> 4; // No need to AND with 0xF, since VSDT (as an uint8_t) has only 8 bits anyway

	// Note, X_D3DVSDT_D3DCOLOR    (0x40) has Size 4, resulting in 4 << 0 == 4, which still equals sizeof(uint32_t)!
	// Note, X_D3DVSDT_NONE        (0x02) has Size 0, resulting in 0 << 2 == 0, which still equals 0!
	// Note, X_D3DVSDT_NORMPACKED3 (0x16) has Size 1, resulting in 1 << 2 == 4, which still equals sizeof(uint32_t)!
	return Size << Shift;
}

void CxbxDecodeVertexAttributeFormat(DWORD AttributeFormat, XboxVertexAttributeDeclarationDecoded_t* pDecoded)
{
	assert(X_D3DVSDT_IsValid(AttributeFormat));

	int XboxSizeInBytes = X_D3DVSDT_SizeInBytes((uint8_t)AttributeFormat);

	int NV2AType = (AttributeFormat & 0x0F) >> 0;
	int NV2ASize = (AttributeFormat & 0xF0) >> 4;

	assert(NV2AType <= 6 && NV2AType != 3); // only 0,1,2 and 4,5,6 are valid types (0:X_D3DVSDT_D3DCOLOR,1:X_D3DVSDT_NORMSHORT*,2:X_D3DVSDT_FLOAT* and 4:X_D3DVSDT_PBYTE*,5:X_D3DVSDT_SHORT*,6:X_D3DVSDT_NORMPACKED3)
	assert(NV2ASize <= 4 || NV2ASize == 7); // only 0..4 and 7 are valid sizes

	int XboxBytesPerUnit = 1 << (NV2AType & 3); // 0 and 4 use byte units, 1 and 5 use short units, 2 and 6 use float/dword units.
	int XboxNrOfUnits = (AttributeFormat == XTL::X_D3DVSDT_FLOAT2H) ? 3 : NV2ASize; // Identity-map sizes to units, except X_D3DVSDT_FLOAT2H becomes 3 units (the only one that has NV2ASize == 7)
	assert(XboxSizeInBytes == XboxNrOfUnits * XboxBytesPerUnit); // Temporary check if X_D3DVSDT_SizeInBytes() output matches the previous determination method
	//int NrOfDimensions = (AttributeFormat == X_D3DVSDT_NORMPACKED3) ? 3 : (AttributeFormat == X_D3DVSDT_FLOAT2H) ? 4 : XboxNrOfUnits;

#define X ((D3DDECLTYPE)-1) // Marks an invalid entry
	static const D3DDECLTYPE c_HostDeclTypeMap[7/*NV2AType*/][8/*NV2ASize*/] = {
		// 0x?0 : X_D3DVSDT_D3DCOLOR=0x40
		{ X, X, X, X, D3DDECLTYPE_D3DCOLOR, X, X, X },
		// 0x?1 : X_D3DVSDT_NORMSHORT1=0x11, X_D3DVSDT_NORMSHORT2=0x21, X_D3DVSDT_NORMSHORT3=0x31, X_D3DVSDT_NORMSHORT4=0x41
		{ X, D3DDECLTYPE_SHORT2N, D3DDECLTYPE_SHORT2N, D3DDECLTYPE_SHORT4N, D3DDECLTYPE_SHORT4N, X, X, X },
		// 0x?2 : X_D3DVSDT_NONE=0x02, X_D3DVSDT_FLOAT1=0x12, X_D3DVSDT_FLOAT2=0x22,X_D3DVSDT_FLOAT3=0x32,X_D3DVSDT_FLOAT4=0x42, X_D3DVSDT_FLOAT2H=0x72
		{ D3DDECLTYPE_UNUSED, D3DDECLTYPE_FLOAT1, D3DDECLTYPE_FLOAT2, D3DDECLTYPE_FLOAT3, D3DDECLTYPE_FLOAT4, X, X, D3DDECLTYPE_FLOAT4 },
		// 0x03 : unused
		{ X, X, X, X, X, X, X, X },
		// 0x?4 : X_D3DVSDT_PBYTE1=0x14, X_D3DVSDT_PBYTE2=0x24, X_D3DVSDT_PBYTE3=0x34, X_D3DVSDT_PBYTE4=0x44
		{ X, D3DDECLTYPE_UBYTE4N, D3DDECLTYPE_UBYTE4N, D3DDECLTYPE_UBYTE4N, D3DDECLTYPE_UBYTE4N, X, X, X },
		// 0x?5 : X_D3DVSDT_SHORT1=0x15, X_D3DVSDT_SHORT2=0x25, X_D3DVSDT_SHORT3=0x35, X_D3DVSDT_SHORT4=0x45
		{ X, D3DDECLTYPE_SHORT2, D3DDECLTYPE_SHORT2, D3DDECLTYPE_SHORT4, D3DDECLTYPE_SHORT4, X, X, X },
		// 0x?6 : X_D3DVSDT_NORMPACKED3 = 0x16
		{ X, D3DDECLTYPE_FLOAT3, X, X, X, X, X, X }
	};
#undef X

	// Lookup optimal host D3DDECLTYPE (might have a different number of units and/or unit size) :
	D3DDECLTYPE HostDeclType = c_HostDeclTypeMap[NV2AType][NV2ASize];

	// Definition of all known host D3DDECLTYPE enums (number of units, size of each unit and optional capability to check)
	static const struct {
		int NrOfUnits; int SizeOfUnit; DWORD D3DDTCAPS;
	} c_HostDeclTypeInfo[(int)D3DDECLTYPE_UNUSED + 1] = {
		/*D3DDECLTYPE_FLOAT1     =  0 */ { 1, sizeof(float) },
		/*D3DDECLTYPE_FLOAT2     =  1 */ { 2, sizeof(float) },
		/*D3DDECLTYPE_FLOAT3     =  2 */ { 3, sizeof(float) },
		/*D3DDECLTYPE_FLOAT4     =  3 */ { 4, sizeof(float) },
		/*D3DDECLTYPE_D3DCOLOR   =  4 */ { 4, sizeof(uint8_t) }, // Makes sizeof(D3DCOLOR)
		/*D3DDECLTYPE_UBYTE4     =  5 */ { 4, sizeof(uint8_t), D3DDTCAPS_UBYTE4 },
		/*D3DDECLTYPE_SHORT2     =  6 */ { 2, sizeof(int16_t) },
		/*D3DDECLTYPE_SHORT4     =  7 */ { 4, sizeof(int16_t) },
		/*D3DDECLTYPE_UBYTE4N    =  8 */ { 4, sizeof(uint8_t), D3DDTCAPS_UBYTE4N },
		/*D3DDECLTYPE_SHORT2N    =  9 */ { 2, sizeof(int16_t), D3DDTCAPS_SHORT2N },
		/*D3DDECLTYPE_SHORT4N    = 10 */ { 4, sizeof(int16_t), D3DDTCAPS_SHORT4N },
		// For completeness, these are not used in HostDeclTypeMap :
		/*D3DDECLTYPE_USHORT2N   = 11 */ { 2, sizeof(int16_t), D3DDTCAPS_USHORT2N },
		/*D3DDECLTYPE_USHORT4N   = 12 */ { 4, sizeof(int16_t), D3DDTCAPS_USHORT4N },
		/*D3DDECLTYPE_UDEC3      = 13 */ { 1, sizeof(DWORD), D3DDTCAPS_UDEC3 },
		/*D3DDECLTYPE_DEC3N      = 14 */ { 1, sizeof(DWORD), D3DDTCAPS_DEC3N },
		/*D3DDECLTYPE_FLOAT16_2  = 15 */ { 2, sizeof(binary16_t), D3DDTCAPS_FLOAT16_2 },
		/*D3DDECLTYPE_FLOAT16_4  = 16 */ { 4, sizeof(binary16_t), D3DDTCAPS_FLOAT16_4 },
		/*D3DDECLTYPE_UNUSED     = 17 */ { 0 } // dynamic size
	};

	// Is there a capability flag to check?
	if (c_HostDeclTypeInfo[HostDeclType].D3DDTCAPS) {
		// and is it unsupported on this driver ?
		if ((g_D3DCaps.DeclTypes & c_HostDeclTypeInfo[HostDeclType].D3DDTCAPS) == 0) {
			// Then fallback to the required number of floats :
			HostDeclType = (D3DDECLTYPE)((int)D3DDECLTYPE_FLOAT1 + XboxNrOfUnits - 1);
		}
	}

	// Gather info on the host declaration (potentially a fallback), and check if it is Xbox compatible :
	int HostNrOfUnits = c_HostDeclTypeInfo[HostDeclType].NrOfUnits;
	int HostSizeInBytes = HostNrOfUnits * c_HostDeclTypeInfo[HostDeclType].SizeOfUnit;
	bool IsCompatible = (XboxNrOfUnits == HostNrOfUnits) && (XboxSizeInBytes == HostSizeInBytes);

	// Return the results :
	pDecoded->XboxDeclarationFormat = AttributeFormat; // Key during conversion
	pDecoded->XboxNrOfUnits = XboxNrOfUnits; // Needed for conversion
	pDecoded->XboxSizeInBytes = XboxSizeInBytes; // Amount of space taken in vertex
	//pDecoded->NrOfDimensions = NrOfDimensions; //
	pDecoded->HostDeclType = HostDeclType;
	pDecoded->IsCompatible = IsCompatible;
	pDecoded->HostNrOfUnits = HostNrOfUnits;
	pDecoded->HostSizeInBytes = HostSizeInBytes;
}

IDirect3DVertexBuffer* CxbxConvertXboxVertexBufferCompletely(
	XTL::X_D3DVertexBuffer* pXboxVertexBuffer,
	unsigned VertexDataSizeInBytes
)
{
	void* pXboxVertexData = GetDataFromXboxResource(pXboxVertexBuffer);

	IDirect3DVertexBuffer* pNewHostVertexBuffer = nullptr;

	HRESULT hRet = g_pD3DDevice->CreateVertexBuffer(
		VertexDataSizeInBytes,
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		0,
		D3DPOOL_DEFAULT,
		&pNewHostVertexBuffer,
		nullptr
	);

	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Failed to create vertex buffer");
	}

	uint8_t* pHostVertexData = nullptr;
	hRet = pNewHostVertexBuffer->Lock(0, 0, (D3DLockData * *)& pHostVertexData, D3DLOCK_DISCARD);
	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Couldn't lock vertex buffer");
	}

	memcpy(pHostVertexData, pXboxVertexData, VertexDataSizeInBytes);

	hRet = pNewHostVertexBuffer->Unlock();
	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Couldn't unlock vertex buffer");
	}

	return pNewHostVertexBuffer;
}

IDirect3DVertexBuffer* CxbxConvertXboxVertexBufferSingleAttribute(
	XTL::X_D3DVertexBuffer* pXboxVertexBuffer,
	unsigned VerticesInBuffer,
	XboxVertexAttributeDeclarationDecoded_t* pDecodedAttribute,
	unsigned AttributeOffsetInStream
)
{
	uint8_t* pXboxVertexData = (uint8_t*)GetDataFromXboxResource(pXboxVertexBuffer);

	unsigned VertexDataSizeInBytes = VerticesInBuffer * pDecodedAttribute->HostSizeInBytes;

	IDirect3DVertexBuffer* pNewHostVertexBuffer = nullptr;

	HRESULT hRet = g_pD3DDevice->CreateVertexBuffer(
		VertexDataSizeInBytes,
		D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
		0,
		D3DPOOL_DEFAULT,
		&pNewHostVertexBuffer,
		nullptr
	);

	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Failed to create vertex buffer");
	}

	uint8_t* pHostVertexData = nullptr;
	hRet = pNewHostVertexBuffer->Lock(0, 0, (D3DLockData * *)& pHostVertexData, D3DLOCK_DISCARD);
	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Couldn't lock vertex buffer");
	}

	pXboxVertexData += AttributeOffsetInStream;

	unsigned XboxStreamStride = 0; // TODO : Get this from somewhere
	for (unsigned i = 0; i < VerticesInBuffer; i++) {
		// TODO : Copy-convert Xbox vertex attribute to host
		pXboxVertexData += XboxStreamStride;
		pHostVertexData += pDecodedAttribute->HostSizeInBytes;
	}

	hRet = pNewHostVertexBuffer->Unlock();
	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Couldn't unlock vertex buffer");
	}

	return pNewHostVertexBuffer;
}

void DumpXboxVertexAttributes(
	XTL::X_VERTEXATTRIBUTEFORMAT* pVertexAttributes
)
{
	using namespace XTL;

	DbgVshPrintf("VERTEXSHADERINPUT Slots[] =\n{\n");

	// Here, parse the Xbox declaration, dumping it like it was an Xbox declaration
	for (int AttributeIndex = 0; AttributeIndex < X_VSH_MAX_ATTRIBUTES; AttributeIndex++) {
		// Fetch the Xbox stream from the pVertexAttributes slots array (this pointer honors overrides) :
		X_VERTEXSHADERINPUT* pAttributeSlot = &(pVertexAttributes->Slots[AttributeIndex]);

		unsigned AttributeFormat = pAttributeSlot->Format;
		unsigned Dump_Type = (AttributeFormat & 0x0F) >> 0;
		unsigned Dump_Size = (AttributeFormat & 0xF0) >> 4;

		if (Dump_Size >= 5) Dump_Type = 3; // Dump all Size's above 4 as an illegal type (X_D3DVSDT_FLOAT2H dumping is fixed up two lines down)
		if (Dump_Type >= 7) Dump_Type = 3; // Dump all Type's above 6 as an illegal type
		if (Dump_Type == 3) Dump_Size = 0; // Dump all illegal types without a size suffix
		if (AttributeFormat == X_D3DVSDT_FLOAT2H)     Dump_Type = 7; // 0x72 = Size 7, Type 2 : Dump c_FormatTypeStr[7] (Dump_Size is already reset to zero)
		if (AttributeFormat == X_D3DVSDT_NONE)        Dump_Type = 8; // 0x02 = Size 0, Type 2 : Dump c_FormatTypeStr[8] (Dump_Size is already zero)
		if (AttributeFormat == X_D3DVSDT_NORMPACKED3) Dump_Size = 0; // 0x16 = Size 1, Type 6 : Don't dump Size suffix (since c_FormatTypeStr[6] already has a "3" suffix)
		if (AttributeFormat == X_D3DVSDT_D3DCOLOR)    Dump_Size = 0; // 0x40 = Size 4, Type 0 : Don't dump Size suffix (since X_D3DVSDT_D3DCOLOR has no size suffix)

		static constexpr char* c_FormatTypeStr[9] = { // Indexed with Dump_Type
			/*[0]*/"X_D3DVSDT_D3DCOLOR",
			/*[1]*/"X_D3DVSDT_NORMSHORT",
			/*[2]*/"X_D3DVSDT_FLOAT",
			/*[3]*/"Illegal!",
			/*[4]*/"X_D3DVSDT_PBYTE",
			/*[5]*/"X_D3DVSDT_SHORT",
			/*[6]*/"X_D3DVSDT_NORMPACKED3", // Longest, 19 characters
			/*[7]*/"X_D3DVSDT_FLOAT2H",
			/*[8]*/"X_D3DVSDT_NONE"
		};

		static constexpr char* c_FormatSizeStr[5] = { // Indexed with Dump_Size
			/*[0]*/" ", // no suffix
			/*[1]*/"1",
			/*[2]*/"2",
			/*[3]*/"3",
			/*[4]*/"4"
		};

		unsigned Dump_Tess = pAttributeSlot->TesselationType;
		if (Dump_Tess > 3) Dump_Tess = 3; // Illegal!

		static constexpr char* c_TessTypeStr[4] = { // Indexed with Dump_Tess
			/*[0]*/"AUTONONE",
			/*[1]*/"AUTONORMAL",
			/*[2]*/"AUTOTEXCOORD", // Longest, 12 characters
			/*[3]*/"Illegal!"
		};

		DbgVshPrintf("\t/*%2d*/{Stream:%2d, Offset:%2d, Format:0x%.02X/*=%19s%s*/, TessType:%1d/*=%-12s*/, TessSource:%2d},",
			AttributeIndex,
			pAttributeSlot->IndexOfStream,
			pAttributeSlot->Offset,
			pAttributeSlot->Format, c_FormatTypeStr[Dump_Type], c_FormatSizeStr[Dump_Size],
			pAttributeSlot->TesselationType, c_TessTypeStr[Dump_Tess],
			pAttributeSlot->TesselationSource
		);

		// Decode the Xbox NV2A attribute format, including mapping to host :
		XboxVertexAttributeDeclarationDecoded_t DecodedAttribute;
		CxbxDecodeVertexAttributeFormat(AttributeFormat, &DecodedAttribute);

		static constexpr char* c_HostDeclTypeStr[1 + (int)D3DDECLTYPE::D3DDECLTYPE_UNUSED] = {
			"D3DDECLTYPE_FLOAT1",
			"D3DDECLTYPE_FLOAT2",
			"D3DDECLTYPE_FLOAT3",
			"D3DDECLTYPE_FLOAT4",
			"D3DDECLTYPE_D3DCOLOR", // Longest, 20 characters
			"D3DDECLTYPE_UBYTE4",
			"D3DDECLTYPE_SHORT2",
			"D3DDECLTYPE_SHORT4",
			"D3DDECLTYPE_UBYTE4N",
			"D3DDECLTYPE_SHORT2N",
			"D3DDECLTYPE_SHORT4N",
			// For completeness, these are not used in HostDeclType :
			"D3DDECLTYPE_USHORT2N",
			"D3DDECLTYPE_USHORT4N",
			"D3DDECLTYPE_UDEC3",
			"D3DDECLTYPE_DEC3N",
			"D3DDECLTYPE_FLOAT16_2",
			"D3DDECLTYPE_FLOAT16_4",
			"D3DDECLTYPE_UNUSED",
		};

		DbgVshPrintf(" // Xbox Units:%d, Bytes:%2d  Host Units:%d, Bytes:%2d, DeclType:%20s%s\n",
			DecodedAttribute.XboxNrOfUnits,
			DecodedAttribute.XboxSizeInBytes,
			//DecodedAttribute.NrOfDimensions,
			DecodedAttribute.HostNrOfUnits,
			DecodedAttribute.HostSizeInBytes,
			c_HostDeclTypeStr[DecodedAttribute.HostDeclType],
			DecodedAttribute.IsCompatible ? "" : " /* xbox ext. */"
		);
	}

	DbgVshPrintf("};\n");
}

// Converts Xbox vertex declaration towards a host vertex declaration,
// splitting up each attribute into a separate stream per host element,
// so that we can use an unmodified host clone of the Xbox vertex buffer
// for all compatible data types, and separate host streams for conversions.
void CxbxUpdateActiveVertexDeclaration(XTL::X_VERTEXATTRIBUTEFORMAT* pVertexAttributes)
{
	LOG_INIT;

	RegVIsPresentInDeclaration.fill(false); // Glue to VshWriteShader / VshConvertShader

	HRESULT hRet = D3D_OK;
	std::array<D3DVERTEXELEMENT, X_VSH_MAX_ATTRIBUTES + 1> HostVertexElements; // +1 is for the closing D3DDECL_END() when all possible attributes are declared

	// Initialize all elements to D3DDECL_END()
	// Later we'll need to create a valid declaration from only the elements that we set
	// D3DDECL_END() has a stream number out of the valid range, so all we need to do is sort the array
	HostVertexElements.fill(D3DDECL_END());
	assert(HostVertexElements[0].Stream > X_VSH_MAX_ATTRIBUTES);

	// Here, parse the Xbox declaration, converting it to a host declaration
	UINT StreamOffsets[16] = { 0 };
	for (int AttributeIndex = 0; AttributeIndex < X_VSH_MAX_ATTRIBUTES; AttributeIndex++) {
		// Fetch the Xbox stream from the pVertexAttributes slots array (this pointer honors overrides) :
		XTL::X_VERTEXSHADERINPUT* pAttributeSlot = &(pVertexAttributes->Slots[AttributeIndex]);
		// We'll call g_pD3DDevice->SetStreamSource for each attribute with these (initially empty) arguments :
		IDirect3DVertexBuffer* pHostVertexBuffer = nullptr;
		// Does this attribute use no storage present the vertex (check this as early as possible to avoid needless processing) ?
		if (pAttributeSlot->Format == XTL::X_D3DVSDT_NONE)
		{
			// Handle tesselating attributes
			switch (pAttributeSlot->TesselationType) {
			case 0: break; // AUTONONE
			case 1: // AUTONORMAL
				// Note : .Stream, .Offset and .Type are copied from pAttributeSlot->TesselationSource in a post-processing step below,
				// because these could all go through an Xbox to host conversion step, so must be copied over afterwards.
				HostVertexElements[AttributeIndex].Method = D3DDECLMETHOD_CROSSUV; // for D3DVSD_TESSNORMAL
				HostVertexElements[AttributeIndex].Usage = D3DDECLUSAGE_NORMAL; // TODO : Is this correct?
				HostVertexElements[AttributeIndex].UsageIndex = 1; // TODO : Is this correct?
				break;
			case 2: // AUTOTEXCOORD
				// HostVertexElements[AttributeIndex].Stream = 0; // The input stream is unused (but must be set to 0), which is already done above
				// HostVertexElements[AttributeIndex].Offset = 0; // The input offset is unused (but must be set to 0), which is already done above
				HostVertexElements[AttributeIndex].Type = D3DDECLTYPE_UNUSED; // The input type for D3DDECLMETHOD_UV must be D3DDECLTYPE_UNUSED (the output type implied by D3DDECLMETHOD_UV is D3DDECLTYPE_FLOAT2)
				HostVertexElements[AttributeIndex].Method = D3DDECLMETHOD_UV; // For X_D3DVSD_MASK_TESSUV
				HostVertexElements[AttributeIndex].Usage = D3DDECLUSAGE_NORMAL; // Note : In Fixed Function Vertex Pipeline, D3DDECLMETHOD_UV must specify usage D3DDECLUSAGE_TEXCOORD or D3DDECLUSAGE_BLENDWEIGHT. TODO : So, what to do?
				HostVertexElements[AttributeIndex].UsageIndex = 1; // TODO ; Is this correct?
				break;
			default:
				assert(false); // invalid TesselationType
			}
		}
		else
		{
			// Each attribute specifies the Xbox stream number it must be read from :
			unsigned XboxStreamIndex = pAttributeSlot->IndexOfStream;
			// When input streams are overriden, the stream indices referenced by attribute slots should stay in bounds!
			assert((g_Xbox_SetVertexShaderInput_Count == 0) || (XboxStreamIndex < g_Xbox_SetVertexShaderInput_Count));

			// Decode the Xbox NV2A attribute format, including mapping to host :
			XboxVertexAttributeDeclarationDecoded_t DecodedAttribute;
			CxbxDecodeVertexAttributeFormat(pAttributeSlot->Format, &DecodedAttribute);
			// All non-X_D3DVSDT_NONE formats should have valid data :
			assert(DecodedAttribute.XboxSizeInBytes > 0);

			int HostAttributeOffset;

			// Can we use the Xbox attribute as-is?
			if (DecodedAttribute.IsCompatible) {
				// All we need to do, is use the host counterpart with unmodified Xbox contents :
				HostAttributeOffset = StreamOffsets[XboxStreamIndex];
				// TODO : What about pXboxStreamInput->Offset; // Can only become non-zero when g_Xbox_SetVertexShaderInput_Data is active
			}
			else
			{
				// TODO : Should we map StreamOffset?
				HostAttributeOffset = 0; // The dedicated stream contains only this attribute, so it's offset in the 'vertex' stream must be zero
			}

			StreamOffsets[XboxStreamIndex] += DecodedAttribute.XboxSizeInBytes; // Note : This implies stream offsets are listed in increasing order!

			static const struct {
				D3DDECLUSAGE Usage;
				BYTE UsageIndex = 0;
			} c_XboxAtrributeInfo[X_VSH_MAX_ATTRIBUTES] = { // TODO : Review this - perhaps everything
				/*[ 0]:*/{ D3DDECLUSAGE_POSITION }, // for X_D3DVSDE_POSITION -- TODO : Use D3DDECLUSAGE_POSITIONT instead?
				/*[ 1]:*/{ D3DDECLUSAGE_BLENDWEIGHT }, // for X_D3DVSDE_BLENDWEIGHT
				/*[ 2]:*/{ D3DDECLUSAGE_NORMAL, 0 } , // for X_D3DVSDE_NORMAL -- TODO : Is this correct?
				/*[ 3]:*/{ D3DDECLUSAGE_COLOR, 0 } , // for X_D3DVSDE_DIFFUSE
				/*[ 4]:*/{ D3DDECLUSAGE_COLOR, 1 } , // for X_D3DVSDE_SPECULAR
				/*[ 5]:*/{ D3DDECLUSAGE_FOG } , // for X_D3DVSDE_FOG
				/*[ 6]:*/{ D3DDECLUSAGE_PSIZE } , // for X_D3DVSDE_POINTSIZE
				/*[ 7]:*/{ D3DDECLUSAGE_COLOR, 2 } , // for X_D3DVSDE_BACKDIFFUSE
				/*[ 8]:*/{ D3DDECLUSAGE_COLOR, 3 } , // for X_D3DVSDE_BACKSPECULAR
				/*[ 9]:*/{ D3DDECLUSAGE_TEXCOORD, 0 } , // for X_D3DVSDE_TEXCOORD0
				/*[10]:*/{ D3DDECLUSAGE_TEXCOORD, 1 } , // for X_D3DVSDE_TEXCOORD1
				/*[11]:*/{ D3DDECLUSAGE_TEXCOORD, 2 } , // for X_D3DVSDE_TEXCOORD2
				/*[12]:*/{ D3DDECLUSAGE_TEXCOORD, 3 } , // for X_D3DVSDE_TEXCOORD3
				/*[13]:*/{ D3DDECLUSAGE_TEXCOORD, 4 } , // TODO : Is this a generic attribute?
				/*[14]:*/{ D3DDECLUSAGE_TEXCOORD, 5 } , // TODO : Is this a generic attribute?
				/*[15]:*/{ D3DDECLUSAGE_TEXCOORD, 6 } , // TODO : Is this a generic attribute?
			};

			HostVertexElements[AttributeIndex].Stream = AttributeIndex; // Stream index matches attribute index (each one has it's own stream)
			HostVertexElements[AttributeIndex].Offset = HostAttributeOffset; // Offset of this attribute in it's stream, measured in bytes
			HostVertexElements[AttributeIndex].Type = DecodedAttribute.HostDeclType; // Host vertex data type
			HostVertexElements[AttributeIndex].Method = D3DDECLMETHOD_DEFAULT; // Processing method. The input type for D3DDECLMETHOD_DEFAULT can be anything. The output type is the same as the input type.
			if (false/*TODO : Honour IsFixedFunction?*/) {
				HostVertexElements[AttributeIndex].Usage = c_XboxAtrributeInfo[AttributeIndex].Usage; // Attribute semantics
				HostVertexElements[AttributeIndex].UsageIndex = c_XboxAtrributeInfo[AttributeIndex].UsageIndex; // Attribute semantic index
			}
			else {
				// D3DDECLUSAGE_TEXCOORD can be useds for any user-defined data
				// We need this because there is no reliable way to detect the real usage
				// Xbox has no concept of 'usage types', it only requires a list of attribute register numbers.
				// So we treat them all as 'user-defined' with an Index of the Vertex Register Index
				// this prevents information loss in shaders due to non-matching dcl types!
				HostVertexElements[AttributeIndex].Usage = D3DDECLUSAGE_TEXCOORD; // Attribute semantics
				HostVertexElements[AttributeIndex].UsageIndex = AttributeIndex; // Attribute semantic index
			}

			// Add this register to the list of declared registers
			RegVIsPresentInDeclaration[AttributeIndex] = true; // glue
		}
	}

	// Post-process host vertex elements that have a D3DDECLMETHOD_CROSSUV method :
	for (int AttributeIndex = 0; AttributeIndex < X_VSH_MAX_ATTRIBUTES; AttributeIndex++) {
		if (HostVertexElements[AttributeIndex].Method == D3DDECLMETHOD_CROSSUV)
		{
			int TesselationSource = pVertexAttributes->Slots[AttributeIndex].TesselationSource;
			// Copy over the Stream, Offset and Type of the host vertex element that serves as 'TesselationSource' :
			HostVertexElements[AttributeIndex].Stream = HostVertexElements[TesselationSource].Stream;
			HostVertexElements[AttributeIndex].Offset = HostVertexElements[TesselationSource].Offset;
			HostVertexElements[AttributeIndex].Type = HostVertexElements[TesselationSource].Type;
			// Note, the input type for D3DDECLMETHOD_CROSSUV can be D3DDECLTYPE_FLOAT[43], D3DDECLTYPE_D3DCOLOR, D3DDECLTYPE_UBYTE4, or D3DDECLTYPE_SHORT4
			// (the output type implied by D3DDECLMETHOD_CROSSUV is D3DDECLTYPE_FLOAT3).
			// TODO : Should we assert this?
		}
	}

	// Note, that host doesn't allow D3DDECLTYPE_UNUSED (apart from D3DDECLMETHOD_UV),
	// so we might need to condense the declarations here afterwards (and hope that the elements'
	// Usage+UsageIndex are correctly connected to the corresponding vertex shader registers) !?!

	// Sort elements according to restrictions (see http://doc.51windows.net/Directx9_SDK/graphics/programmingguide/gettingstarted/vertexdeclaration/vertexdeclaration.htm
	// All unmodified D3DDECL_END() elements (there's at least one) are moved to the end, giving a valid D3D9 VertexDeclaration :
	std::sort(HostVertexElements.begin(), HostVertexElements.end(), [](D3DVERTEXELEMENT& x, D3DVERTEXELEMENT& y)
		{ return std::tie(x.Stream, x.Method, x.Offset, x.Usage, x.UsageIndex) < std::tie(y.Stream, y.Method, y.Offset, y.Usage, y.UsageIndex); });

	// Set the vertex declaration we've prepare above :
	IDirect3DVertexDeclaration* pHostVertexDeclaration = nullptr;
	hRet = g_pD3DDevice->CreateVertexDeclaration(HostVertexElements.data(), &pHostVertexDeclaration);
	//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateVertexDeclaration"); // TODO : Why does this fail?!?
	hRet = g_pD3DDevice->SetVertexDeclaration(pHostVertexDeclaration);
	//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexDeclaration");
	// Now that the declaration is set (or not, when conversion failed?)
	if (pHostVertexDeclaration)
	{
		// Throw away the temporary object reference already (host still holds a reference anyway):
		hRet = pHostVertexDeclaration->Release();
		//DEBUG_D3DRESULT(hRet, "pHostVertexDeclaration->Release");
	}
}

// Converts Xbox vertex stream buffers towards a host vertex buffers,
// using an unmodified host clone of the Xbox vertex buffer
// for all compatible data types, and separate host streams for incompatible ones.
void CxbxConvertActiveVertexStreams(
	XTL::X_VERTEXATTRIBUTEFORMAT* pVertexAttributes,
	unsigned VerticesInBuffer
)
{
	LOG_INIT;

	// Here, parse the Xbox declaration, converting it to a host declaration
	for (int AttributeIndex = 0; AttributeIndex < X_VSH_MAX_ATTRIBUTES; AttributeIndex++) {
		// Fetch the Xbox stream from the pVertexAttributes slots array (this pointer honors overrides) :
		XTL::X_VERTEXSHADERINPUT* pAttributeSlot = &(pVertexAttributes->Slots[AttributeIndex]);
		// We'll call g_pD3DDevice->SetStreamSource for each attribute with these (initially empty) arguments :
		IDirect3DVertexBuffer* pHostVertexBuffer = nullptr;
		UINT StreamStride = 0;
		UINT StreamOffset = 0;
		// Does this attribute use no storage present the vertex (check this as early as possible to avoid needless processing) ?
		if (pAttributeSlot->Format != XTL::X_D3DVSDT_NONE) {
			// Each attribute specifies the Xbox stream number it must be read from :
			unsigned XboxStreamIndex = pAttributeSlot->IndexOfStream;
			// When input streams are overriden, the stream indices referenced by attribute slots should stay in bounds!
			assert((g_Xbox_SetVertexShaderInput_Count == 0) || (XboxStreamIndex < g_Xbox_SetVertexShaderInput_Count));

			// Fetch the input stream from the pStreamInputs array (this pointer honors overrides) :
			XTL::X_STREAMINPUT* pXboxStreamInput = &GetXboxVertexStreamInput(XboxStreamIndex);
			// Fetch the Xbox vertex buffer for this input stream :
			XTL::X_D3DVertexBuffer* pXboxVertexBuffer = pXboxStreamInput->VertexBuffer;
			// TODO : What if pXboxVertexBuffer is null? (This implies that the Xbox vertex declaration mentions a stream that's not yet been set using D3DDevice_SetVertexShaderInput or D3DDevice_SetStreamSource*
			// For now, just skip the vertex element :
			if (pXboxVertexBuffer == xbnullptr)
				continue;

			// Decode the Xbox NV2A attribute format, including mapping to host :
			XboxVertexAttributeDeclarationDecoded_t DecodedAttribute;
			CxbxDecodeVertexAttributeFormat(pAttributeSlot->Format, &DecodedAttribute);
			// All non-X_D3DVSDT_NONE formats should have valid data :
			assert(pXboxVertexBuffer);
			assert(DecodedAttribute.XboxSizeInBytes > 0);

			// Can we use the Xbox attribute as-is?
			if (DecodedAttribute.IsCompatible) {
				// All we need to do, is use the host counterpart with unmodified Xbox contents :
				StreamStride = pXboxStreamInput->Stride;
				pHostVertexBuffer = CxbxConvertXboxVertexBufferCompletely(pXboxVertexBuffer, VerticesInBuffer * StreamStride);
				StreamOffset = pXboxStreamInput->Offset; // Can only become non-zero when g_Xbox_SetVertexShaderInput_Data is active
			}
			else
			{
				pHostVertexBuffer = CxbxConvertXboxVertexBufferSingleAttribute(pXboxVertexBuffer, VerticesInBuffer, &DecodedAttribute, pXboxStreamInput->Offset);
				// TODO : Should we map StreamOffset?
				assert(StreamStride == 0); // This dedicated VertexBuffer has no other contents, so it has no stride. TODO : Or must we set DecodedAttribute.HostSizeInBytes?
			}
		}

		// Give each attribute (vertex element in host terms) it's own vertex stream :
		HRESULT hRet = g_pD3DDevice->SetStreamSource(
			/*StreamNumber=*/AttributeIndex,
			/*pStreamData=*/pHostVertexBuffer,
			/*OffsetInBytes=*/StreamOffset + pAttributeSlot->Offset,
			/*Stride=*/StreamStride);
		//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetStreamSource");
	}
}

std::string DebugPrependLineNumbers(std::string shaderString) {
	std::stringstream shader(shaderString);
	auto debugShader = std::stringstream();

	int i = 1;
	for (std::string line; std::getline(shader, line); ) {
		auto lineNumber = std::to_string(i++);
		auto paddedLineNumber = lineNumber.insert(0, 3 - lineNumber.size(), ' ');
		debugShader << "/* " << paddedLineNumber << " */ " << line << "\n";
	}

	return debugShader.str();
}

// Recompile and set the result on host via SetVertexShader :
// recompile xbox vertex shader function
void CxbxRecompileVertexProgram()
{
	using namespace XTL;


	//CxbxVertexShader* pCxbxVertexShader = GetCxbxVertexShader(XboxVertexShaderHandle);
	//if (pCxbxVertexShader = nullptr) {
	//	SetCxbxVertexShader(, pCxbxVertexShader);
	//}

	DWORD* pXboxFunction = &g_Xbox_VertexShader_FunctionSlots[g_Xbox_VertexShader_FunctionSlots_StartAddress];

	XboxVertexShaderDecoder VshDecoder;

	// Before we begin, detect the size of the function :
	uint32_t* pToken = (uint32_t*)pXboxFunction;
	while (!VshDecoder.VshIsEndToken(pToken)) {
		pToken += X_VSH_INSTRUCTION_SIZE;
	}
	pToken += X_VSH_INSTRUCTION_SIZE;

	unsigned XboxFunctionSizeInBytes = (intptr_t)pToken - (intptr_t)pXboxFunction;
	unsigned XboxFunctionSize = XboxFunctionSizeInBytes / X_VSH_INSTRUCTION_SIZE / sizeof(DWORD);

	if (XboxFunctionSize > X_VSH_MAX_INSTRUCTION_COUNT) {
		LOG_TEST_CASE("Parsing didn't result in a shader program?");
	}

// bool* pbUseDeclarationOnly,

	// Include HLSL header and footer as raw strings :
	static std::string hlsl_template[2] = {
		#include "core\hle\D3D8\Direct3D9\CxbxVertexShaderTemplate.hlsl"
	};

	// Decode the vertex shader program tokens into an intermediate representation
	while (VshDecoder.VshConvertToIntermediate(pToken)) {
		pToken += X_VSH_INSTRUCTION_SIZE;
	}

	auto hlsl_stream = std::stringstream();
	hlsl_stream << hlsl_template[0]; // Start with the HLSL template header
	if (!VshDecoder.BuildShader(hlsl_stream)) {
		// Do not attempt to compile empty shaders
		// This is a declaration only shader, so there is no function to recompile
//		*pbUseDeclarationOnly = true;
		return;// D3D_OK;
	}

	hlsl_stream << hlsl_template[1]; // Finish with the HLSL template footer
	std::string hlsl_str = hlsl_stream.str();

	DbgVshPrintf("--- HLSL conversion ---\n");
	DbgVshPrintf(DebugPrependLineNumbers(hlsl_str).c_str());
	DbgVshPrintf("-----------------------\n");

//	bool bNoReservedConstants;

	// D3DVERTEXELEMENT * pRecompiledDeclaration,
	ID3DBlob* pRecompiledShader = nullptr;
	ID3DBlob* pErrors = nullptr;
	HRESULT hRet = 0;

	// Level 0 for fastest runtime compilation
	// TODO Can we recompile an optimized shader in the background?
	UINT flags1 = D3DCOMPILE_OPTIMIZATION_LEVEL0;

	for (int i = 0; i < 2; i++) {
		hRet = D3DCompile(
			hlsl_str.c_str(),
			hlsl_str.length(),
			nullptr, // pSourceName
			nullptr, // pDefines
			nullptr, // pInclude // TODO precompile x_* HLSL functions?
			"main", // shader entry poiint
			"vs_3_0", // shader profile
			flags1, // flags1
			0, // flags2
			&pRecompiledShader, // out
			&pErrors // ppErrorMsgs out
		);
		if (FAILED(hRet)) {
			// Attempt to retry in compatibility mode, this allows some vertex-state shaders to compile
			// Test Case: Spy vs Spy
			flags1 |= D3DCOMPILE_ENABLE_BACKWARDS_COMPATIBILITY;
		} else
			break;
	}

	if (FAILED(hRet)) {
		LOG_TEST_CASE("Couldn't assemble recompiled vertex shader");
		//EmuLog(LOG_LEVEL::WARNING, "Couldn't assemble recompiled vertex shader");
	}

	// Determine the log level
	auto hlslErrorLogLevel = FAILED(hRet) ? LOG_LEVEL::ERROR2 : LOG_LEVEL::DEBUG;
	if (pErrors) {
		// Log HLSL compiler errors
		EmuLog(hlslErrorLogLevel, "%s", (char*)(pErrors->GetBufferPointer()));
		pErrors->Release();
		pErrors = nullptr;
	}

	LOG_CHECK_ENABLED(LOG_LEVEL::DEBUG)
		if (g_bPrintfOn)
			if (!FAILED(hRet)) {
				// Log disassembly
				hRet = D3DDisassemble(
					pRecompiledShader->GetBufferPointer(),
					pRecompiledShader->GetBufferSize(),
					D3D_DISASM_ENABLE_DEFAULT_VALUE_PRINTS | D3D_DISASM_ENABLE_INSTRUCTION_NUMBERING,
					NULL,
					&pErrors
				);
				if (pErrors) {
					EmuLog(hlslErrorLogLevel, "%s", (char*)(pErrors->GetBufferPointer()));
					pErrors->Release();
				}
			}

	IDirect3DVertexShader* pHostVertexShader = nullptr;

	if (SUCCEEDED(hRet) && pRecompiledShader != nullptr) {
		hRet = g_pD3DDevice->CreateVertexShader((DWORD*)(pRecompiledShader->GetBufferPointer()), &pHostVertexShader);
		//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateVertexShader");
	}

	if (pRecompiledShader != nullptr) {
		pRecompiledShader->Release();
		pRecompiledShader = nullptr;
	}
/*
	//	pCxbxVertexShader->pXboxDeclarationCopy = (DWORD*)malloc(XboxDeclarationCount * sizeof(DWORD));
	//	memcpy(pCxbxVertexShader->pXboxDeclarationCopy, pDeclaration, XboxDeclarationCount * sizeof(DWORD));
	pCxbxVertexShader->XboxFunctionSize = 0;
	pCxbxVertexShader->pXboxFunctionCopy = NULL;
	pCxbxVertexShader->XboxVertexShaderType = X_VST_NORMAL; // TODO : This can vary
	pCxbxVertexShader->XboxNrAddressSlots = (XboxFunctionSize - sizeof(X_VSH_SHADER_HEADER)) / X_VSH_INSTRUCTION_SIZE_BYTES;
	pCxbxVertexShader->HostFVF = 0;
	pCxbxVertexShader->pHostVertexShader = nullptr;
	//	pCxbxVertexShader->XboxDeclarationCount = XboxDeclarationCount;
	//	pCxbxVertexShader->HostDeclarationSize = HostDeclarationSize;
		// Save the status, to remove things later
		// pCxbxVertexShader->XboxStatus = hRet; // Not even used by VshHandleIsValidShader()

	if (SUCCEEDED(hRet))
	{
		if (pFunction != NULL)
		{
			pCxbxVertexShader->XboxFunctionSize = XboxFunctionSize;
			pCxbxVertexShader->pXboxFunctionCopy = (DWORD*)malloc(XboxFunctionSize);
			memcpy(pCxbxVertexShader->pXboxFunctionCopy, pFunction, XboxFunctionSize);
		}

		pCxbxVertexShader->pHostVertexShader = pHostVertexShader;
	}
	else
	{
		LOG_TEST_CASE("Falling back to FVF shader");
		pCxbxVertexShader->HostFVF = D3DFVF_XYZ | D3DFVF_TEX0;
	}

	// Register the host Vertex Shader
	SetCxbxVertexShader(Handle, pCxbxVertexShader);


	if (pHostVertexShader)
	{
		HRESULT hRet = g_pD3DDevice->SetVertexShader(pHostVertexShader);
		//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShader()");
		pHostVertexShader->Release();
		//DEBUG_D3DRESULT(hRet, "pHostVertexShader->Release()");
	}
*/
	// TODO : Pick parts from Cxbx_CreateVertexShader()
	// like UpdateViewPortOffsetAndScaleConstants();
}

void CxbxUpdateActiveVertexShader(unsigned VerticesInBuffer)
{
	using namespace XTL;

	LOG_INIT;

	assert(g_pD3DDevice);

	// Take overrides on declarations (as optionally set by SetVertexShaderInput) into account :
	X_VERTEXATTRIBUTEFORMAT* pVertexAttributes = GetXboxVertexAttributes();

	DumpXboxVertexAttributes(pVertexAttributes);

	CxbxUpdateActiveVertexDeclaration(pVertexAttributes);
	CxbxRecompileVertexProgram();
	CxbxConvertActiveVertexStreams(pVertexAttributes, VerticesInBuffer);
}

void CxbxCopyVertexShaderFunctionSlots(unsigned Address, unsigned NumberOfSlots, void* FunctionData)
{
	if (Address >= X_VSH_MAX_INSTRUCTION_COUNT) {
		LOG_TEST_CASE("Address would start too high!");
		return;
	}

	if (NumberOfSlots == 0) {
		LOG_TEST_CASE("Received zero slots?");
		return;
	}

	if (Address + NumberOfSlots > X_VSH_MAX_INSTRUCTION_COUNT) {
		LOG_TEST_CASE("Copy would overflow available instruction slots!");
		return;
	}

	// Each slot is either 4 DWORDS or 4 floats
	DWORD SizeInBytes = NumberOfSlots * X_VSH_INSTRUCTION_SIZE_BYTES;
	// Copy this data towards our storage for this, at the indicated address offset :
	memcpy(&(g_Xbox_VertexShader_FunctionSlots[Address * X_VSH_INSTRUCTION_SIZE]), FunctionData, SizeInBytes);
}

// Our vertex shader related EMU_PATCH implementations :

// D3DDevice_LoadVertexShaderProgram splits the given function buffer into batch-wise pushes to the NV2A
// This implies that the given pFunction input (after the header DWORD) is a RAW vertex shader program!
// This is different from LoadVertexShader, which contains NV2A opcodes in it's function DWORD buffer
void CxbxImpl_LoadVertexShaderProgram(DWORD* pFunction, DWORD Address)
{
	using namespace XTL;

	assert(pFunction);

	X_VSH_SHADER_HEADER* ShaderHeader = (X_VSH_SHADER_HEADER*)pFunction;

	// Verify the vertex shader header
	{
		switch (ShaderHeader->Version) {
		case VERSION_XVS: // 0x2078 // 'x ' Xbox vertex shader
			break;
		case VERSION_XVSS: // 0x7378 // 'xs' Xbox vertex state shader
			LOG_TEST_CASE("Might not support vertex state shaders?");
			break;
		case VERSION_XVSW: // 0x7778 // 'xw' Xbox vertex read/write shader
			EmuLog(LOG_LEVEL::WARNING, "Might not support vertex read/write shaders?");
			//hRet = E_FAIL;
			break;
		default:
			LOG_TEST_CASE("Unexpected vertex shader header version?");
			EmuLog(LOG_LEVEL::WARNING, "Unknown vertex shader version 0x%02X", ShaderHeader->Version);
			//hRet = E_FAIL;
			return;
		}

		if (ShaderHeader->NumInst == 0)
			LOG_TEST_CASE("Unexpected vertex shader header zero length?");

		if (ShaderHeader->NumInst > X_VSH_MAX_INSTRUCTION_COUNT)
			LOG_TEST_CASE("Too large vertex shader header length?");
	}

	// The second word in the function indicates it's size, expressed in 4-DWORD-sized slots :
	DWORD NumberOfSlots = ShaderHeader->NumInst;
	// Copy this data towards our storage for this, at the indicated address offset :
	CxbxCopyVertexShaderFunctionSlots(Address, NumberOfSlots, &pFunction[1]);
}

void CxbxImpl_LoadVertexShader(DWORD Handle, DWORD Address)
{
	using namespace XTL;

	if (!VshHandleIsVertexShader(Handle)) {
		LOG_TEST_CASE("D3DDevice_LoadVertexShader called without a shader!");
		return;
	}

	X_D3DVertexShader* XboxVertexShader = VshHandleToXboxVertexShader(Handle);
	if (Address + XboxVertexShader->FunctionSize > X_VSH_MAX_INSTRUCTION_COUNT) {
		LOG_TEST_CASE("D3DDevice_LoadVertexShader would exceed too far!");
		return;
	}


	static UINT ConstantStartRegister = 0;
	unsigned CurrentProgramIndex = Address; // Note : Make this static, if below code would be extended to parse NV097_SET_TRANSFORM_PROGRAM_LOAD
	HRESULT hRet;
	DWORD* ProgramData = XboxVertexShader->FunctionData;

	// Parse the Xbox Vertex shader program, and transfer the constants present therein to host.
	while (*ProgramData != 0) { // Each program ends with a zero DWORD
		// TODO : Limit to XboxVertexShader->FunctionSize
		// Parse the limited set of Vertex program NV2A Push buffer commands :
		NV2AMETHOD NV2ACommand = *ProgramData++;
		switch (PUSH_METHOD(NV2ACommand)) {
		case NV097_SET_TRANSFORM_PROGRAM: { // == 0x00000B00
			// Copy a batch of instructions :
			unsigned NumberOfDWORDs = PUSH_COUNT(NV2ACommand); // Fetch the number of instruction dwords
			assert(NumberOfDWORDs > 0); // There should at least be one instruction
			assert((NumberOfDWORDs & 3) == 0); // Instructions are expected to come in 4 DWORD pairs
			unsigned NumberOfTokens = NumberOfDWORDs / X_VSH_INSTRUCTION_SIZE;
			// Copy the function data into the indicated address slots :
			CxbxCopyVertexShaderFunctionSlots(CurrentProgramIndex, NumberOfTokens, ProgramData);
			// Skip this number of DWORD's for the next batch :
			CurrentProgramIndex += NumberOfDWORDs;
			ProgramData += NumberOfDWORDs;
			continue;
		}
		case NV097_SET_TRANSFORM_CONSTANT_LOAD: // == 0x00001EA4
			ConstantStartRegister = *ProgramData++; // Remember the starting constant register
			assert(ConstantStartRegister < X_D3DVS_CONSTREG_COUNT);
			continue;
		case NV097_SET_TRANSFORM_CONSTANT: { // == 0x00000B80
			UINT Vector4fCount = PUSH_COUNT(NV2ACommand); // Fetch the number of constants
			// TODO : Are constants counted in DWORD's too?
			assert(Vector4fCount > 0 && Vector4fCount < 8); // There can be at most 8 constants per batch
			// Set this batch of constant data on host :
			hRet = g_pD3DDevice->SetVertexShaderConstantF(ConstantStartRegister, (float*)ProgramData, Vector4fCount);
			//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShaderConstantF");
			// Each constant takes 4 floats, so skip that number of DWORD's for the next batch :
			unsigned NumberOfDWORDs = Vector4fCount * X_VSH_INSTRUCTION_SIZE;
			ProgramData += NumberOfDWORDs;
			continue;
		}
		default:
			assert(false); // Unexpected NV2AMethod in vertex program!
		}

		break; // 
	}
}

// Note : SelectVertexShaderDirect needs no EMUPATCH CxbxImpl_..., since it just calls SelectVertexShader

void CxbxImpl_SelectVertexShader(DWORD Handle, DWORD Address)
{
	using namespace XTL;

	// Handle can be null if the current Xbox VertexShader is assigned
	if (Handle)
	{
		assert(!VshHandleIsFVF(Handle));
		// Handle can be an address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
		// If Handle is assigned, it becomes the new current Xbox VertexShader,
		// which resets a bit of state (nv2a execution mode, viewport, ?)
		g_Xbox_VertexShader_Handle = Handle; // TODO : Remove bit here, or on use??
	}

	// Either way, the given address slot is selected as the start of the current vertex shader program
	// Address always indicates a previously loaded vertex shader slot (from where the program is used).
	g_Xbox_VertexShader_FunctionSlots_StartAddress = Address;
}

void CxbxImpl_SetVertexShader(DWORD Handle)
{
	using namespace XTL;

	// Is the user-supplied handle a FVF?
	if (VshHandleIsFVF(Handle))
	{
		// Then we're more interested in what the D3DDevice_SetVertexShader trampoline
		// stored in the Xbox D3Device.m_pVertexShader field :
		// Note : This requires CxbxImpl_SetVertexShader to be called _AFTER_ the trampoline in D3DDevice_SetVertexShader!!
		// Note : XboxVertexShaders.g_XboxAddr_pVertexShader is located at the start of CreateDevice (which calls D3DDevice_SetVertexShader)
		if (XboxVertexShaders.g_XboxAddr_pVertexShader)
			g_Xbox_VertexShader_Handle = *XboxVertexShaders.g_XboxAddr_pVertexShader | X_D3DFVF_RESERVED0; // Note : Mark Xbox vertex-shader as non-FVF (we should put the pointer in another variable and deprecate this handle)
		else
			g_Xbox_VertexShader_Handle = Handle;
	}
	else
		g_Xbox_VertexShader_Handle = Handle;

	// g_Xbox_VertexShader_FVF = VshHandleIsFVF(Handle) ? Handle : 0; // enable if we ever need the Xbox FVF

	g_Xbox_VertexShader_FunctionSlots_StartAddress = 0;
}

void CxbxImpl_SetVertexShaderInput(DWORD Handle, UINT StreamCount, XTL::X_STREAMINPUT* pStreamInputs)
{
	using namespace XTL;

	// If Handle is NULL, all VertexShader input state is cleared.
	// Otherwise, Handle is the address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// (Thus, a FVF handle is an invalid argument.)

	if (Handle == NULL)
	{
		// Xbox doesn't remember a null-handle - this may be an XDK bug!
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
			memcpy(g_Xbox_SetVertexShaderInput_Data, pStreamInputs, StreamCount * sizeof(XTL::X_STREAMINPUT)); // Make a copy of the supplied StreamInputs array
		}

		g_Xbox_SetVertexShaderInput_Attributes = pXboxVertexShader->VertexAttribute; // Copy this vertex shaders's attribute slots
	}
}

#if 0
// Old code from D3DDevice_CreateVertexShader patch
// (kept enabled to make sure this at least compiles)
HRESULT Cxbx_CreateVertexShader
(
	DWORD* pDeclaration,
	DWORD* pFunction,
	DWORD Handle
)
{
	using namespace XTL;

	HRESULT hRet;
	// The Xbox CreateVertexShader function does the following:
	// Allocates an Xbox VertexShader struct
	// Sets reference count to 1
	// Puts Usage in VertexShader->Flags
	// If pFunction is not null, it points to DWORDS shader type, length and a binary compiled xbox vertex shader
	// If pDeclaration is not null, it's parsed, resulting in a number of constants
	// Parse results are pushed to the push buffer
	// Sets other fields
	// Handle recieved the addres of the new shader, or-ed with 1 (D3DFVF_RESERVED0)

	if (g_pD3DDevice == nullptr) {
		LOG_TEST_CASE("D3DDevice_CreateVertexShader called before Direct3D_CreateDevice");
		// We lie to allow the game to continue for now, but it probably won't work well
		return STATUS_SUCCESS;
	}

	// HACK: TODO: support this situation
	if (pDeclaration == NULL) {
		LOG_TEST_CASE("Vertex shader without declaration");
		return D3D_OK;
	}

	// Now, we can create the host vertex shader
//	DWORD             XboxDeclarationCount = 0;
//	DWORD             HostDeclarationSize = 0;
	CxbxVertexShader* pCxbxVertexShader = (CxbxVertexShader*)calloc(1, sizeof(CxbxVertexShader));
	D3DVERTEXELEMENT* pRecompiledDeclaration = nullptr;

	pRecompiledDeclaration = EmuRecompileVshDeclaration((DWORD*)pDeclaration,
		/*bIsFixedFunction=*/pFunction == NULL,
		//		&XboxDeclarationCount,
		//		&HostDeclarationSize,
		&pCxbxVertexShader->VertexShaderInfo);

	// Create the vertex declaration
	hRet = g_pD3DDevice->CreateVertexDeclaration(pRecompiledDeclaration, &pCxbxVertexShader->pHostVertexDeclaration);
	//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateVertexDeclaration");

	if (FAILED(hRet)) {
		// NOTE: This is a fatal error because it ALWAYS triggers a crash within DrawVertices if not set
		CxbxKrnlCleanup("Failed to create Vertex Declaration");
	}
	g_pD3DDevice->SetVertexDeclaration(pCxbxVertexShader->pHostVertexDeclaration);
	//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexDeclaration");

	ID3DBlob*  pRecompiledBuffer = nullptr;
	DWORD         XboxFunctionSize = 0;
	DWORD* pRecompiledFunction = NULL;
	if (SUCCEEDED(hRet) && pFunction)
	{
		bool bUseDeclarationOnly = false;

		hRet = EmuRecompileVshFunction((DWORD*)pFunction,
			/*bNoReservedConstants=*/g_Xbox_VertexShaderConstantMode == X_D3DSCM_NORESERVEDCONSTANTS,
			pRecompiledDeclaration,
			&bUseDeclarationOnly,
			&XboxFunctionSize,
			&pRecompiledBuffer);
		if (SUCCEEDED(hRet))
		{
			if (!bUseDeclarationOnly)
				pRecompiledFunction = (DWORD*)pRecompiledBuffer->GetBufferPointer();
			else
				pRecompiledFunction = NULL;
		}
		else
		{
			pRecompiledFunction = NULL;
			EmuLog(LOG_LEVEL::WARNING, "Couldn't recompile vertex shader function.");
		}
	}

	//EmuLog(LOG_LEVEL::DEBUG, "MaxVertexShaderConst = %d", g_D3DCaps.MaxVertexShaderConst);

	IDirect3DVertexShader* pHostVertexShader = nullptr;
	if (SUCCEEDED(hRet) && pRecompiledFunction != nullptr)
	{
		hRet = g_pD3DDevice->CreateVertexShader
		(
			pRecompiledFunction,
			&pHostVertexShader
		);
		//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateVertexShader");
	}

	if (pRecompiledBuffer != nullptr)
	{
		pRecompiledBuffer->Release();
		pRecompiledBuffer = nullptr;
	}

	free(pRecompiledDeclaration);

	//	pCxbxVertexShader->pXboxDeclarationCopy = (DWORD*)malloc(XboxDeclarationCount * sizeof(DWORD));
	//	memcpy(pCxbxVertexShader->pXboxDeclarationCopy, pDeclaration, XboxDeclarationCount * sizeof(DWORD));
	pCxbxVertexShader->XboxFunctionSize = 0;
	pCxbxVertexShader->pXboxFunctionCopy = NULL;
	pCxbxVertexShader->XboxVertexShaderType = X_VST_NORMAL; // TODO : This can vary
	pCxbxVertexShader->XboxNrAddressSlots = (XboxFunctionSize - sizeof(X_VSH_SHADER_HEADER)) / X_VSH_INSTRUCTION_SIZE_BYTES;
	pCxbxVertexShader->HostFVF = 0;
	pCxbxVertexShader->pHostVertexShader = nullptr;
	//	pCxbxVertexShader->XboxDeclarationCount = XboxDeclarationCount;
	//	pCxbxVertexShader->HostDeclarationSize = HostDeclarationSize;
		// Save the status, to remove things later
		// pCxbxVertexShader->XboxStatus = hRet; // Not even used by VshHandleIsValidShader()

	if (SUCCEEDED(hRet))
	{
		if (pFunction != NULL)
		{
			pCxbxVertexShader->XboxFunctionSize = XboxFunctionSize;
			pCxbxVertexShader->pXboxFunctionCopy = (DWORD*)malloc(XboxFunctionSize);
			memcpy(pCxbxVertexShader->pXboxFunctionCopy, pFunction, XboxFunctionSize);
		}

		pCxbxVertexShader->pHostVertexShader = pHostVertexShader;
	}
	else
	{
		LOG_TEST_CASE("Falling back to FVF shader");
		pCxbxVertexShader->HostFVF = D3DFVF_XYZ | D3DFVF_TEX0;
	}

	// Register the host Vertex Shader
	SetCxbxVertexShader(Handle, pCxbxVertexShader);

	if (FAILED(hRet))
	{
#ifdef _DEBUG_TRACK_VS
		if (pFunction)
		{
			char pFileName[30];
			static int FailedShaderCount = 0;
			X_VSH_SHADER_HEADER* pHeader = (X_VSH_SHADER_HEADER*)pFunction;
			EmuLog(LOG_LEVEL::WARNING, "Couldn't create vertex shader!");
			sprintf(pFileName, "failed%05d.xvu", FailedShaderCount);
			FILE* f = fopen(pFileName, "wb");
			if (f)
			{
				fwrite(pFunction, sizeof(X_VSH_SHADER_HEADER) + pHeader->NumInst * 16, 1, f);
				fclose(f);
			}
			FailedShaderCount++;
		}
#endif // _DEBUG_TRACK_VS
		//hRet = D3D_OK;
	}

	return hRet;
}
#endif

// recompile xbox vertex shader function
extern HRESULT EmuRecompileVshFunction
(
	DWORD* pXboxFunction,
	bool          bNoReservedConstants,
	D3DVERTEXELEMENT* pRecompiledDeclaration,
	bool* pbUseDeclarationOnly,
	DWORD* pXboxFunctionSize,
	ID3DBlob** ppRecompiledShader
)
{
	XTL::X_VSH_SHADER_HEADER* pXboxVertexShaderHeader = (XTL::X_VSH_SHADER_HEADER*)pXboxFunction;
	uint32_t* pToken;
	XboxVertexShaderDecoder VshDecoder;
	ID3DBlob* pErrors = nullptr;
	HRESULT             hRet = 0;

	// TODO: support this situation..
	if (pXboxFunction == xbnullptr) {
		return E_FAIL;
	}

	// Initialize output arguments to zero
	*pbUseDeclarationOnly = false;
	*pXboxFunctionSize = 0;
	*ppRecompiledShader = nullptr;

	switch (pXboxVertexShaderHeader->Version) {
	case VERSION_XVS:
		break;
	case VERSION_XVSS:
		LOG_TEST_CASE("Might not support vertex state shaders?");
		break;
	case VERSION_XVSW:
		EmuLog(LOG_LEVEL::WARNING, "Might not support vertex read/write shaders?");
		hRet = E_FAIL;
		break;
	default:
		EmuLog(LOG_LEVEL::WARNING, "Unknown vertex shader version 0x%02X", pXboxVertexShaderHeader->Version);
		hRet = E_FAIL;
		break;
	}

	if (!SUCCEEDED(hRet)) return hRet;

	// Include HLSL header and footer as raw strings :
	static std::string hlsl_template[2] = {
		#include "core\hle\D3D8\Direct3D9\CxbxVertexShaderTemplate.hlsl"
	};

	// Decode the vertex shader program tokens into an intermediate representation
	pToken = (uint32_t*)((uintptr_t)pXboxFunction + sizeof(XTL::X_VSH_SHADER_HEADER));
	while (VshDecoder.VshConvertToIntermediate(pToken)) {
		pToken += X_VSH_INSTRUCTION_SIZE;
	}

	// The size of the shader is
	pToken += X_VSH_INSTRUCTION_SIZE; // always at least one token
	*pXboxFunctionSize = (intptr_t)pToken - (intptr_t)pXboxFunction;

	auto hlsl_stream = std::stringstream();
	hlsl_stream << hlsl_template[0]; // Start with the HLSL template header
	if (!VshDecoder.BuildShader(hlsl_stream)) {
		// Do not attempt to compile empty shaders
		// This is a declaration only shader, so there is no function to recompile
		*pbUseDeclarationOnly = true;
		return D3D_OK;
	}

	hlsl_stream << hlsl_template[1]; // Finish with the HLSL template footer
	std::string hlsl_str = hlsl_stream.str();

	DbgVshPrintf("--- HLSL conversion ---\n");
	DbgVshPrintf(DebugPrependLineNumbers(hlsl_str).c_str());
	DbgVshPrintf("-----------------------\n");

	// Level 0 for fastest runtime compilation
	// TODO Can we recompile an optimized shader in the background?
	UINT flags1 = D3DCOMPILE_OPTIMIZATION_LEVEL0;

	hRet = D3DCompile(
		hlsl_str.c_str(),
		hlsl_str.length(),
		nullptr, // pSourceName
		nullptr, // pDefines
		nullptr, // pInclude // TODO precompile x_* HLSL functions?
		"main", // shader entry poiint
		"vs_3_0", // shader profile
		flags1, // flags1
		0, // flags2
		ppRecompiledShader, // out
		&pErrors // ppErrorMsgs out
	);
	if (FAILED(hRet)) {
		// Attempt to retry in compatibility mode, this allows some vertex-state shaders to compile
		// Test Case: Spy vs Spy
		flags1 |= D3DCOMPILE_ENABLE_BACKWARDS_COMPATIBILITY;
		hRet = D3DCompile(
			hlsl_str.c_str(),
			hlsl_str.length(),
			nullptr, // pSourceName
			nullptr, // pDefines
			nullptr, // pInclude // TODO precompile x_* HLSL functions?
			"main", // shader entry poiint
			"vs_3_0", // shader profile
			flags1, // flags1
			0, // flags2
			ppRecompiledShader, // out
			&pErrors // ppErrorMsgs out
		);

		if (FAILED(hRet)) {
			LOG_TEST_CASE("Couldn't assemble recompiled vertex shader");
			//EmuLog(LOG_LEVEL::WARNING, "Couldn't assemble recompiled vertex shader");
		}
	}

	// Determine the log level
	auto hlslErrorLogLevel = FAILED(hRet) ? LOG_LEVEL::ERROR2 : LOG_LEVEL::DEBUG;
	if (pErrors) {
		// Log HLSL compiler errors
		EmuLog(hlslErrorLogLevel, "%s", (char*)(pErrors->GetBufferPointer()));
		pErrors->Release();
		pErrors = nullptr;
	}

	LOG_CHECK_ENABLED(LOG_LEVEL::DEBUG)
	if (g_bPrintfOn)
	if (!FAILED(hRet)) {
		// Log disassembly
		hRet = D3DDisassemble(
			(*ppRecompiledShader)->GetBufferPointer(),
			(*ppRecompiledShader)->GetBufferSize(),
			D3D_DISASM_ENABLE_DEFAULT_VALUE_PRINTS | D3D_DISASM_ENABLE_INSTRUCTION_NUMBERING,
			NULL,
			&pErrors
		);
		if (pErrors) {
			EmuLog(hlslErrorLogLevel, "%s", (char*)(pErrors->GetBufferPointer()));
			pErrors->Release();
		}
	}

	return hRet;
}
