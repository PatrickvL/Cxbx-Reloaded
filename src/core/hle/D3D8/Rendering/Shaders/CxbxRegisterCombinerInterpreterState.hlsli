// CxbxRegisterCombinerInterpreterState.hlsli — shared C++ / HLSL header
//
// Defines the constant buffer layout for the register combiner interpreter
// ubershader.  Included by both the HLSL pixel shader and the C++ backend
// to ensure the cbuffer layout stays in sync.
//
// Field ordering matches X_D3DPIXELSHADERDEF, with the following differences:
//   - PSConstant0/1 are float4 RGBA [0..1] (Xbox stores packed ARGB DWORDs;
//     the C++ upload code converts at upload time).
//   - PSFinalCombinerConstant is a float4[2] array instead of two separate DWORDs.
//   - ColorSign[4] and FogColor are appended (not part of X_D3DPIXELSHADERDEF).
//   - Software-only fields (PSC0Mapping, PSC1Mapping, PSFinalCombinerConstants)
//     are omitted.

#ifdef __cplusplus
#pragma once
#include <cstdint>

// C++: 16-byte-aligned scalar uint occupying one SM5 constant register.
// SM5 cbuffer packing places each uint array element in its own 16-byte slot.
struct alignas(16) RCI_UintReg { uint32_t value; uint32_t _pad[3]; };

// C++: 16-byte float4 matching HLSL's native float4.
struct alignas(16) RCI_Float4 { float x, y, z, w; };

// Begin the struct definition for C++
#define RCI_BEGIN struct RCInterpreterCBLayout {
#define RCI_END   };
#define RCI_UINT(name)       RCI_UintReg name
#define RCI_UINT_ARRAY(name, n)  RCI_UintReg name[n]
#define RCI_FLOAT4(name)     RCI_Float4 name
#define RCI_FLOAT4_ARRAY(name, n) RCI_Float4 name[n]

#else
// HLSL: the cbuffer keyword defines the layout directly.
// Single uint fields must be padded to 16 bytes (one constant register) to
// match the C++ alignas(16) layout.  Without this, SM5 cbuffer packing rules
// pack consecutive scalars together, shifting all subsequent field offsets.
#define RCI_BEGIN cbuffer RCInterpreterCBLayout : register(b0) {
#define RCI_END   };
#define RCI_UINT(name)       uint name; uint3 _pad_##name
#define RCI_UINT_ARRAY(name, n)  uint name[n]
#define RCI_FLOAT4(name)     float4 name
#define RCI_FLOAT4_ARRAY(name, n) float4 name[n]

#endif

// ============================================================
// Shared cbuffer / struct layout
//
// Field order MUST match between HLSL and C++ — do not reorder.
// ============================================================
RCI_BEGIN
    RCI_UINT_ARRAY(PSAlphaInputs, 8);          // Alpha combiner A..D input specs
    RCI_UINT(PSFinalCombinerInputsABCD);        // (A<<24)|(B<<16)|(C<<8)|D
    RCI_UINT(PSFinalCombinerInputsEFG);         // (E<<24)|(F<<16)|(G<<8)|settings
    RCI_FLOAT4_ARRAY(PSConstant0, 8);           // Per-stage C0 color constant [0..1]
    RCI_FLOAT4_ARRAY(PSConstant1, 8);           // Per-stage C1 color constant [0..1]
    RCI_UINT_ARRAY(PSAlphaOutputs, 8);
    RCI_UINT_ARRAY(PSRGBInputs, 8);             // RGB combiner A..D input specs
    RCI_UINT(PSCompareMode);                    // Clip-plane comparison mode
    RCI_FLOAT4_ARRAY(PSFinalCombinerConstant, 2); // FC0, FC1
    RCI_UINT_ARRAY(PSRGBOutputs, 8);
    RCI_UINT(PSCombinerCount);                  // (flags<<8)|numStages
    RCI_UINT(PSTextureModes);                   // 4 x 5-bit modes
    RCI_UINT(PSDotMapping);                     // Dot-product normal mapping
    RCI_UINT(PSInputTexture);                   // Input-texture for dependent modes
    RCI_FLOAT4_ARRAY(ColorSign, 4);             // Per-stage: 0=keep, >0=u->s, <0=s->u
    RCI_FLOAT4(FogColor);                       // rgb=fog color constant; a=unused
    // --- Post-processing state (matches compiled PS c23..c43) ---
    RCI_FLOAT4(TexFmtFixup);                    // Per-stage fixup: 0=id,1=.gbar,2=.abgr,3=lum,4=alum
    RCI_FLOAT4(AlphaTest);                      // x=enable, y=ref [0..1], z=func [D3DCMPFUNC]
    RCI_FLOAT4_ARRAY(ColorKeyOp, 4);            // Per-stage color key operation
    RCI_FLOAT4_ARRAY(ColorKeyColor, 4);         // Per-stage color key color
    RCI_FLOAT4_ARRAY(BEM, 4);                   // Per-stage bump env material matrix
    RCI_FLOAT4_ARRAY(LUM, 4);                   // Per-stage bump luminance (scale, offset)
    RCI_FLOAT4(FogInfo);                        // x=tableMode, y=density, z=start, w=end
    RCI_UINT(FogEnable);                        // Fog enable flag
RCI_END

// Clean up macros to avoid polluting the global namespace
#undef RCI_BEGIN
#undef RCI_END
#undef RCI_UINT
#undef RCI_UINT_ARRAY
#undef RCI_FLOAT4
#undef RCI_FLOAT4_ARRAY

#ifdef __cplusplus
static_assert(sizeof(RCInterpreterCBLayout) == 1312, "RC cbuffer layout size mismatch");

// Convert a packed DWORD ARGB color (0xAARRGGBB) to RCI_Float4 RGBA [0..1]
inline RCI_Float4 DwordColorToFloat4(uint32_t color)
{
    RCI_Float4 f;
    f.x = ((color >> 16) & 0xFF) / 255.0f; // R
    f.y = ((color >> 8)  & 0xFF) / 255.0f; // G
    f.z = ( color        & 0xFF) / 255.0f; // B
    f.w = ((color >> 24) & 0xFF) / 255.0f; // A
    return f;
}
#endif
