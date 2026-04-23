// CxbxRegisterCombinerInterpreterState.hlsli — shared C++ / HLSL header
//
// Defines the AUXILIARY constant buffer layout for the register combiner
// interpreter ubershader.  Contains only software-computed fields that
// have no direct PGRAPH register equivalent.
//
// Register-backed combiner state (PSAlphaInputs, PSRGBOutputs, PSConstant0,
// PSCombinerCount, etc.) is now read directly from the raw PGRAPH regs[]
// StructuredBuffer<uint> (see CxbxPGRAPHRegs.hlsli).
//
// This aux cbuffer also carries a few fields that ARE derived from PGRAPH
// registers but require C++-side adjustment before the shader sees them
// (PSTextureModes after AdjustTextureModes, PSFinalCombinerInputs after
// AdjustFinalCombiner).

#ifdef __cplusplus
#pragma once
#include <cstdint>

// C++: 16-byte-aligned scalar uint occupying one SM5 constant register.
struct alignas(16) RCI_UintReg { uint32_t value; uint32_t _pad[3]; };

// C++: 16-byte float4 matching HLSL's native float4.
struct alignas(16) RCI_Float4 { float x, y, z, w; };

#define RCI_BEGIN struct PSAuxCBLayout {
#define RCI_END   };
#define RCI_UINT(name)       RCI_UintReg name
#define RCI_UINT_ARRAY(name, n)  RCI_UintReg name[n]
#define RCI_FLOAT4(name)     RCI_Float4 name
#define RCI_FLOAT4_ARRAY(name, n) RCI_Float4 name[n]

#else
// HLSL cbuffer definition.  Single uint fields are padded to 16 bytes to
// match the C++ alignas(16) layout.
#define RCI_BEGIN cbuffer PSAuxCBLayout : register(b0) {
#define RCI_END   };
#define RCI_UINT(name)       uint name; uint3 _pad_##name
#define RCI_UINT_ARRAY(name, n)  uint name[n]; uint3 _pad_##name
#define RCI_FLOAT4(name)     float4 name
#define RCI_FLOAT4_ARRAY(name, n) float4 name[n]

#endif

// ============================================================
// Auxiliary cbuffer — software-computed fields only.
// Register-backed fields are accessed via g_PGRegs[].
//
// Field order MUST match between HLSL and C++ — do not reorder.
// ============================================================
RCI_BEGIN
    // --- Adjusted register values (C++ modifies before upload) ---
    RCI_UINT(PSTextureModes);                   // After AdjustTextureModes
    RCI_UINT(PSFinalCombinerInputsABCD);        // After AdjustFinalCombiner
    RCI_UINT(PSFinalCombinerInputsEFG);         // After AdjustFinalCombiner

    // --- Software-computed per-stage state ---
    RCI_FLOAT4_ARRAY(ColorSign, 4);             // Per-stage sign conversion
    RCI_FLOAT4(TexFmtFixup);                    // Per-stage format fixup
    RCI_FLOAT4_ARRAY(ColorKeyOp, 4);            // Per-stage color key operation
    RCI_FLOAT4_ARRAY(ColorKeyColor, 4);         // Per-stage color key color
    RCI_FLOAT4(AlphaKill);                      // Per-stage alpha kill

    // --- Fog state (not yet migrated to PGRAPH) ---
    RCI_FLOAT4(FogInfo);                        // x=tableMode, y=density, z=start, w=end
    RCI_UINT(FogEnable);                        // Fog enable flag

    // --- Misc runtime state ---
    RCI_FLOAT4(FrontFaceInfo);                  // x=FrontFaceFactor
RCI_END

// Clean up macros
#undef RCI_BEGIN
#undef RCI_END
#undef RCI_UINT
#undef RCI_UINT_ARRAY
#undef RCI_FLOAT4
#undef RCI_FLOAT4_ARRAY

#ifdef __cplusplus
static_assert(sizeof(PSAuxCBLayout) == 320, "PSAuxCBLayout size mismatch");
#endif
