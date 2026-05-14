// CxbxRegisterCombinerInterpreterState.hlsli — shared C++ / HLSL header
//
// Defines the AUXILIARY constant buffer layout for the register combiner
// interpreter ubershader.  Contains ONLY fields that require host-side
// information not available in NV2A PGRAPH registers.
//
// All other PS auxiliary state (PSTextureModes, FinalCombinerInputs,
// ColorKeyOp, ColorKeyColor, AlphaKill, FogInfo, FogEnable, FrontFaceInfo,
// ShadowCompare) is now derived directly from the PGRAPH SRV in-shader.
// See CxbxPSAuxFromPGRAPH.hlsli.

#ifdef __cplusplus
#pragma once
#include <cstdint>

// C++: 16-byte float4 matching HLSL's native float4.
struct alignas(16) RCI_Float4 { float x, y, z, w; };

#define RCI_BEGIN struct PSAuxCBLayout {
#define RCI_END   };
#define RCI_FLOAT4(name)     RCI_Float4 name
#define RCI_FLOAT4_ARRAY(name, n) RCI_Float4 name[n]

#else
// HLSL cbuffer definition.
#define RCI_BEGIN cbuffer PSAuxCBLayout : register(b0) {
#define RCI_END   };
#define RCI_FLOAT4(name)     float4 name
#define RCI_FLOAT4_ARRAY(name, n) float4 name[n]

#endif

// ============================================================
// Auxiliary cbuffer — host-dependent fields only.
// These require information not in NV2A PGRAPH registers:
//   - ColorSign: needs host DXGI texture format signedness
//   - TexFmtFixup: needs host resource cache (channel swizzle)
//   - DepthScale: needs xfctx (viewport Z scale, not in regs[])
//   - DepthTexAlias: needs host RT lookup (RT-as-texture detection)
//
// Field order MUST match between HLSL and C++ — do not reorder.
// ============================================================
RCI_BEGIN
    RCI_FLOAT4_ARRAY(ColorSign, 4);             // Per-stage per-channel sign conversion (0/1/-1)
    RCI_FLOAT4(TexFmtFixup);                    // Per-stage format fixup code
    RCI_FLOAT4(DepthScale);                     // x=VPSCL.z (viewport Z scale for DOT_ZW)
    RCI_FLOAT4(DepthTexAlias);                  // Per-stage: 1.0=D24S8, 2.0=D16 (host RT alias)
RCI_END

// Clean up macros
#undef RCI_BEGIN
#undef RCI_END
#undef RCI_FLOAT4
#undef RCI_FLOAT4_ARRAY

#ifdef __cplusplus
static_assert(sizeof(PSAuxCBLayout) == 112, "PSAuxCBLayout size mismatch");
#endif
