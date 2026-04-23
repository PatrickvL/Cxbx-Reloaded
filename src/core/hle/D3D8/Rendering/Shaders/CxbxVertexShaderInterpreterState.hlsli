// CxbxVertexShaderInterpreterState.hlsli — shared C++ / HLSL header
//
// Defines instruction field constants for the NV2A vertex shader interpreter.
// Included by both the HLSL vertex shader and the C++ backend.
//
// Data flow (no cbuffer — all data via SRVs):
//   - g_PGRegs (StructuredBuffer<uint>, t12): shared PGRAPH register array,
//     provides CHEOPS_PROGRAM_START and any future register-backed VS fields.
//     Same GPU buffer bound to PS and VS stages.
//   - g_XFPR (StructuredBuffer<uint4>, t5): NV2A XFPR (Transform Program
//     RAM) — 136 instruction slots × 128-bit containers (92 bits used per
//     instruction in Kelvin ISA encoding).  On real hardware this is on-chip
//     XF SRAM accessed via the RDI (Register Direct Interface); the CPU
//     reaches it through the NV097_SET_TRANSFORM_PROGRAM method range with
//     the write pointer in NV_PGRAPH_CHEOPS_OFFSET.PROG_LD_PTR.
//     The shader reads CHEOPS_PROGRAM_START from g_PGRegs to find the
//     first active slot and loops until FLD_FINAL.

#ifdef __cplusplus
#pragma once
#include <cstdint>
#else
// HLSL side: XFPR (Transform Program RAM) SRV, uploaded from pg->program_data[136][4]
StructuredBuffer<uint4> g_XFPR : register(t5);
#endif

// ============================================================
// NV2A vertex shader instruction field layout
// (matches XbVertexShaderDecoder.cpp FieldMapping)
// ============================================================

// XFPR capacity: 136 instruction slots (indices 0-0x87), shared across all VPEs.
// On Kelvin each slot is a 128-bit container with 92 bits of actual instruction data.
#ifdef __cplusplus
static constexpr uint32_t XFPR_LENGTH = 136;
#else
static const uint XFPR_LENGTH = 136;
#endif

// ============================================================
// Instruction field bit positions and sizes
// Organized by DWORD (SubToken) within the 128-bit instruction
// ============================================================

// DWORD 1 (SubToken 1)
#ifdef __cplusplus
#define VSI_CONST static constexpr uint32_t
#else
#define VSI_CONST static const uint
#endif

// ILU opcode: SubToken 1, bits [27:25] (3 bits)
VSI_CONST VSI_FLD_ILU_SHIFT       = 25;
VSI_CONST VSI_FLD_ILU_MASK        = 0x7;

// MAC opcode: SubToken 1, bits [24:21] (4 bits)
VSI_CONST VSI_FLD_MAC_SHIFT       = 21;
VSI_CONST VSI_FLD_MAC_MASK        = 0xF;

// Constant register index: SubToken 1, bits [20:13] (8 bits)
VSI_CONST VSI_FLD_CONST_SHIFT     = 13;
VSI_CONST VSI_FLD_CONST_MASK      = 0xFF;

// Vertex register index: SubToken 1, bits [12:9] (4 bits)
VSI_CONST VSI_FLD_V_SHIFT         = 9;
VSI_CONST VSI_FLD_V_MASK          = 0xF;

// Input A: SubToken 1 bits [8:0] + SubToken 2 bits [31:26]
VSI_CONST VSI_FLD_A_NEG_BIT1      = 8;   // SubToken 1, bit 8
VSI_CONST VSI_FLD_A_SWZ_X_SHIFT1  = 6;   // SubToken 1, bits [7:6]
VSI_CONST VSI_FLD_A_SWZ_Y_SHIFT1  = 4;   // SubToken 1, bits [5:4]
VSI_CONST VSI_FLD_A_SWZ_Z_SHIFT1  = 2;   // SubToken 1, bits [3:2]
VSI_CONST VSI_FLD_A_SWZ_W_SHIFT1  = 0;   // SubToken 1, bits [1:0]

// DWORD 2 (SubToken 2)
VSI_CONST VSI_FLD_A_R_SHIFT       = 28;  // bits [31:28] (4 bits)
VSI_CONST VSI_FLD_A_R_MASK        = 0xF;
VSI_CONST VSI_FLD_A_MUX_SHIFT     = 26;  // bits [27:26] (2 bits)
VSI_CONST VSI_FLD_A_MUX_MASK      = 0x3;

// Input B
VSI_CONST VSI_FLD_B_NEG_BIT2      = 25;  // SubToken 2, bit 25
VSI_CONST VSI_FLD_B_SWZ_X_SHIFT2  = 23;  // bits [24:23]
VSI_CONST VSI_FLD_B_SWZ_Y_SHIFT2  = 21;  // bits [22:21]
VSI_CONST VSI_FLD_B_SWZ_Z_SHIFT2  = 19;  // bits [20:19]
VSI_CONST VSI_FLD_B_SWZ_W_SHIFT2  = 17;  // bits [18:17]
VSI_CONST VSI_FLD_B_R_SHIFT       = 13;  // bits [16:13] (4 bits)
VSI_CONST VSI_FLD_B_R_MASK        = 0xF;
VSI_CONST VSI_FLD_B_MUX_SHIFT     = 11;  // bits [12:11] (2 bits)
VSI_CONST VSI_FLD_B_MUX_MASK      = 0x3;

// Input C
VSI_CONST VSI_FLD_C_NEG_BIT2      = 10;  // SubToken 2, bit 10
VSI_CONST VSI_FLD_C_SWZ_X_SHIFT2  = 8;   // bits [9:8]
VSI_CONST VSI_FLD_C_SWZ_Y_SHIFT2  = 6;   // bits [7:6]
VSI_CONST VSI_FLD_C_SWZ_Z_SHIFT2  = 4;   // bits [5:4]
VSI_CONST VSI_FLD_C_SWZ_W_SHIFT2  = 2;   // bits [3:2]
VSI_CONST VSI_FLD_C_R_HIGH_SHIFT2 = 0;   // bits [1:0] (2 bits, high part)
VSI_CONST VSI_FLD_C_R_HIGH_MASK   = 0x3;

// DWORD 3 (SubToken 3)
VSI_CONST VSI_FLD_C_R_LOW_SHIFT3  = 30;  // bits [31:30] (2 bits, low part)
VSI_CONST VSI_FLD_C_R_LOW_MASK    = 0x3;
VSI_CONST VSI_FLD_C_MUX_SHIFT3    = 28;  // bits [29:28] (2 bits)
VSI_CONST VSI_FLD_C_MUX_MASK      = 0x3;

// Output fields
VSI_CONST VSI_FLD_OUT_MAC_MASK_SHIFT = 24; // bits [27:24] (4 bits)
VSI_CONST VSI_FLD_OUT_MAC_MASK_MASK  = 0xF;
VSI_CONST VSI_FLD_OUT_R_SHIFT       = 20; // bits [23:20] (4 bits)
VSI_CONST VSI_FLD_OUT_R_MASK        = 0xF;
VSI_CONST VSI_FLD_OUT_ILU_MASK_SHIFT = 16; // bits [19:16] (4 bits)
VSI_CONST VSI_FLD_OUT_ILU_MASK_MASK  = 0xF;
VSI_CONST VSI_FLD_OUT_O_MASK_SHIFT  = 12; // bits [15:12] (4 bits)
VSI_CONST VSI_FLD_OUT_O_MASK_MASK   = 0xF;
VSI_CONST VSI_FLD_OUT_ORB_BIT3     = 11;  // bit 11 (0=context, 1=output)
VSI_CONST VSI_FLD_OUT_ADDRESS_SHIFT = 3;   // bits [10:3] (8 bits)
VSI_CONST VSI_FLD_OUT_ADDRESS_MASK  = 0xFF;
VSI_CONST VSI_FLD_OUT_MUX_BIT3     = 2;   // bit 2 (0=MAC, 1=ILU)
VSI_CONST VSI_FLD_A0X_BIT3         = 1;   // bit 1
VSI_CONST VSI_FLD_FINAL_BIT3       = 0;   // bit 0

// MUX values for input source selection
VSI_CONST VSI_MUX_R = 1; // Temporary register (matches PARAM_R in hardware encoding)
VSI_CONST VSI_MUX_V = 2; // Vertex input register (matches PARAM_V)
VSI_CONST VSI_MUX_C = 3; // Constant register (matches PARAM_C)

// MAC opcode values
VSI_CONST VSI_MAC_NOP = 0;
VSI_CONST VSI_MAC_MOV = 1;
VSI_CONST VSI_MAC_MUL = 2;
VSI_CONST VSI_MAC_ADD = 3;
VSI_CONST VSI_MAC_MAD = 4;
VSI_CONST VSI_MAC_DP3 = 5;
VSI_CONST VSI_MAC_DPH = 6;
VSI_CONST VSI_MAC_DP4 = 7;
VSI_CONST VSI_MAC_DST = 8;
VSI_CONST VSI_MAC_MIN = 9;
VSI_CONST VSI_MAC_MAX = 10;
VSI_CONST VSI_MAC_SLT = 11;
VSI_CONST VSI_MAC_SGE = 12;
VSI_CONST VSI_MAC_ARL = 13;

// ILU opcode values
VSI_CONST VSI_ILU_NOP = 0;
VSI_CONST VSI_ILU_MOV = 1;
VSI_CONST VSI_ILU_RCP = 2;
VSI_CONST VSI_ILU_RCC = 3;
VSI_CONST VSI_ILU_RSQ = 4;
VSI_CONST VSI_ILU_EXP = 5;
VSI_CONST VSI_ILU_LOG = 6;
VSI_CONST VSI_ILU_LIT = 7;

// Writemask bits (same encoding as NV2A)
VSI_CONST VSI_MASK_X = 0x8;
VSI_CONST VSI_MASK_Y = 0x4;
VSI_CONST VSI_MASK_Z = 0x2;
VSI_CONST VSI_MASK_W = 0x1;

#undef VSI_CONST
