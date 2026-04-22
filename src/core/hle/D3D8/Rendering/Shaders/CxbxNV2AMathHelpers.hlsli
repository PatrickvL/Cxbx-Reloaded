// CxbxNV2AMathHelpers.hlsli
//
// NV2A-accurate math helpers shared by vertex and pixel shaders.
//
// NV2A hardware treats 0 * anything = 0, even when the other operand is
// inf or NaN. Standard GPU float math produces NaN for 0 * inf, which
// breaks vertex transforms and register combiner blending. These helpers
// enforce the NV2A zero-propagation rule per component.

#ifndef CXBX_NV2A_MATH_HELPERS_HLSLI
#define CXBX_NV2A_MATH_HELPERS_HLSLI

// ============================================================
// NV2A-accurate multiply: 0 * anything = 0, even 0 * inf
// ============================================================
float4 nv2a_mul(float4 a, float4 b)
{
    bool4 isZero = (a == 0.0f) | (b == 0.0f);
    return isZero ? (float4)0.0f : a * b;
}

float3 nv2a_mul3(float3 a, float3 b)
{
    bool3 isZero = (a == 0.0f) | (b == 0.0f);
    return isZero ? (float3)0.0f : a * b;
}

float nv2a_mul1(float a, float b)
{
    bool isZero = (a == 0.0f) | (b == 0.0f);
    return isZero ? 0.0f : a * b;
}

// ============================================================
// NV2A-accurate dot products using nv2a_mul per component
// ============================================================
float nv2a_dot3(float3 a, float3 b)
{
    float3 m = nv2a_mul3(a, b);
    return m.x + m.y + m.z;
}

float nv2a_dot4(float4 a, float4 b)
{
    float4 m = nv2a_mul(a, b);
    return m.x + m.y + m.z + m.w;
}

#endif // CXBX_NV2A_MATH_HELPERS_HLSLI
