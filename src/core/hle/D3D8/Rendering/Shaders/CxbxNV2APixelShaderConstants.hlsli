// CxbxNV2APixelShaderConstants.hlsli — shared C++ / HLSL header
//
// Defines NV2A register combiner constants that are used by both the
// HLSL CxbxRegisterCombinerInterpreter ubershader and the C++ backend.
// Include from either language; the preprocessor adapts types automatically.
//
// In C++ this header is included via:
//   #include "CxbxNV2APixelShaderConstants.hlsli"
// In HLSL this header is included via:
//   #include "CxbxNV2APixelShaderConstants.hlsli"

#ifdef __cplusplus
#pragma once
#include <cstdint>
// C++: use uint32_t for constants; HLSL uses `uint` natively.
#define NV2A_CONST static constexpr uint32_t
// If XbPixelShader.h was already included, its enums define the same names.
// Skip our definitions to avoid redefinition errors.
#ifdef XBPIXELSHADER_H
#define NV2A_PS_CONSTANTS_ALREADY_DEFINED
#endif
#else
// HLSL: `static const uint` is the idiomatic constant form.
#define NV2A_CONST static const uint
#endif

#ifndef NV2A_PS_CONSTANTS_ALREADY_DEFINED

// ============================================================
// PS_REGISTER — 4-bit register index (bits [3:0] of an input byte)
// ============================================================
NV2A_CONST PS_REGISTER_ZERO         = 0x00;
NV2A_CONST PS_REGISTER_DISCARD      = 0x00; // write alias for ZERO
NV2A_CONST PS_REGISTER_C0           = 0x01;
NV2A_CONST PS_REGISTER_C1           = 0x02;
NV2A_CONST PS_REGISTER_FOG          = 0x03;
NV2A_CONST PS_REGISTER_V0           = 0x04;
NV2A_CONST PS_REGISTER_V1           = 0x05;
// 0x06, 0x07 — reserved
NV2A_CONST PS_REGISTER_T0           = 0x08;
NV2A_CONST PS_REGISTER_T1           = 0x09;
NV2A_CONST PS_REGISTER_T2           = 0x0a;
NV2A_CONST PS_REGISTER_T3           = 0x0b;
NV2A_CONST PS_REGISTER_R0           = 0x0c;
NV2A_CONST PS_REGISTER_R1           = 0x0d;
NV2A_CONST PS_REGISTER_V1R0_SUM     = 0x0e;
NV2A_CONST PS_REGISTER_EF_PROD      = 0x0f;
NV2A_CONST PS_REGISTER_MASK         = 0x0f;

// ============================================================
// PS_CHANNEL — bit 4 of an input byte
// ============================================================
NV2A_CONST PS_CHANNEL_RGB           = 0x00;
NV2A_CONST PS_CHANNEL_BLUE          = 0x00; // alias when used as alpha source
NV2A_CONST PS_CHANNEL_ALPHA         = 0x10;
NV2A_CONST PS_CHANNEL_MASK          = 0x10;

// ============================================================
// PS_INPUTMAPPING — bits [7:5] of an input byte (pre-shifted)
// ============================================================
NV2A_CONST PS_INPUTMAPPING_UNSIGNED_IDENTITY = 0x00;
NV2A_CONST PS_INPUTMAPPING_UNSIGNED_INVERT   = 0x20;
NV2A_CONST PS_INPUTMAPPING_EXPAND_NORMAL     = 0x40;
NV2A_CONST PS_INPUTMAPPING_EXPAND_NEGATE     = 0x60;
NV2A_CONST PS_INPUTMAPPING_HALFBIAS_NORMAL   = 0x80;
NV2A_CONST PS_INPUTMAPPING_HALFBIAS_NEGATE   = 0xa0;
NV2A_CONST PS_INPUTMAPPING_SIGNED_IDENTITY   = 0xc0;
NV2A_CONST PS_INPUTMAPPING_SIGNED_NEGATE     = 0xe0;
NV2A_CONST PS_INPUTMAPPING_MASK              = 0xe0;

// ============================================================
// PS_TEXTUREMODES — 5-bit value per stage
// ============================================================
NV2A_CONST PS_TEXTUREMODES_NONE                  = 0x00;
NV2A_CONST PS_TEXTUREMODES_PROJECT2D             = 0x01;
NV2A_CONST PS_TEXTUREMODES_PROJECT3D             = 0x02;
NV2A_CONST PS_TEXTUREMODES_CUBEMAP               = 0x03;
NV2A_CONST PS_TEXTUREMODES_PASSTHRU              = 0x04;
NV2A_CONST PS_TEXTUREMODES_CLIPPLANE             = 0x05;
NV2A_CONST PS_TEXTUREMODES_BUMPENVMAP            = 0x06;
NV2A_CONST PS_TEXTUREMODES_BUMPENVMAP_LUM        = 0x07;
NV2A_CONST PS_TEXTUREMODES_BRDF                  = 0x08;
NV2A_CONST PS_TEXTUREMODES_DOT_ST                = 0x09;
NV2A_CONST PS_TEXTUREMODES_DOT_ZW                = 0x0a;
NV2A_CONST PS_TEXTUREMODES_DOT_RFLCT_DIFF        = 0x0b;
NV2A_CONST PS_TEXTUREMODES_DOT_RFLCT_SPEC        = 0x0c;
NV2A_CONST PS_TEXTUREMODES_DOT_STR_3D            = 0x0d;
NV2A_CONST PS_TEXTUREMODES_DOT_STR_CUBE          = 0x0e;
NV2A_CONST PS_TEXTUREMODES_DPNDNT_AR             = 0x0f;
NV2A_CONST PS_TEXTUREMODES_DPNDNT_GB             = 0x10;
NV2A_CONST PS_TEXTUREMODES_DOTPRODUCT            = 0x11;
NV2A_CONST PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST  = 0x12;
NV2A_CONST PS_TEXTUREMODES_MASK                  = 0x1f;

// ============================================================
// PS_FINALCOMBINERSETTING — settings byte of PSFinalCombinerInputsEFG
// ============================================================
NV2A_CONST PS_FINALCOMBINERSETTING_CLAMP_SUM     = 0x80;
NV2A_CONST PS_FINALCOMBINERSETTING_COMPLEMENT_V1 = 0x40;
NV2A_CONST PS_FINALCOMBINERSETTING_COMPLEMENT_R0 = 0x20;

// ============================================================
// PS_COMBINERCOUNTFLAGS — packed into PSCombinerCount bits [23:8]
//   PS_COMBINERCOUNT(count, flags) = (flags << 8) | count
// These are flag values *before* the <<8 shift.
// ============================================================
NV2A_CONST PS_COMBINERCOUNT_MUX_LSB   = 0x0000;
NV2A_CONST PS_COMBINERCOUNT_MUX_MSB   = 0x0001;
NV2A_CONST PS_COMBINERCOUNT_SAME_C0   = 0x0000;
NV2A_CONST PS_COMBINERCOUNT_UNIQUE_C0 = 0x0010;
NV2A_CONST PS_COMBINERCOUNT_SAME_C1   = 0x0000;
NV2A_CONST PS_COMBINERCOUNT_UNIQUE_C1 = 0x0100;

// ============================================================
// PS_COMBINEROUTPUT — flags field of PSRGBOutputs / PSAlphaOutputs
//   PS_COMBINEROUTPUTS(ab, cd, mux_sum, flags) = (flags<<12)|(mux_sum<<8)|(ab<<4)|cd
// These are flag values *before* the <<12 shift.
// ============================================================
NV2A_CONST PS_COMBINEROUTPUT_CD_DOT_PRODUCT    = 0x01; // RGB only
NV2A_CONST PS_COMBINEROUTPUT_AB_DOT_PRODUCT    = 0x02; // RGB only
NV2A_CONST PS_COMBINEROUTPUT_AB_CD_MUX         = 0x04; // 3rd output = MUX(AB,CD) on R0.a
NV2A_CONST PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS = 0x08; // subtract 0.5 before scaling
NV2A_CONST PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA  = 0x40; // RGB only
NV2A_CONST PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA  = 0x80; // RGB only

#endif // NV2A_PS_CONSTANTS_ALREADY_DEFINED

// ============================================================
// Packing layout: PS_COMBINERINPUTS(a, b, c, d)
//   = (a << 24) | (b << 16) | (c << 8) | d
// Shift amounts for extracting each 8-bit input byte.
// ============================================================
NV2A_CONST PS_COMBINERINPUTS_A_SHIFT = 24;
NV2A_CONST PS_COMBINERINPUTS_B_SHIFT = 16;
NV2A_CONST PS_COMBINERINPUTS_C_SHIFT =  8;
NV2A_CONST PS_COMBINERINPUTS_D_SHIFT =  0;

// ============================================================
// Packing layout: PS_COMBINEROUTPUTS(ab, cd, mux_sum, flags)
//   = (flags << 12) | (mux_sum << 8) | (ab << 4) | cd
// Shift amounts for extracting each 4-bit output register nibble.
// ============================================================
NV2A_CONST PS_COMBINEROUTPUTS_CD_SHIFT      =  0;
NV2A_CONST PS_COMBINEROUTPUTS_AB_SHIFT      =  4;
NV2A_CONST PS_COMBINEROUTPUTS_MUX_SUM_SHIFT =  8;
NV2A_CONST PS_COMBINEROUTPUTS_FLAGS_SHIFT   = 12;

#undef NV2A_CONST
