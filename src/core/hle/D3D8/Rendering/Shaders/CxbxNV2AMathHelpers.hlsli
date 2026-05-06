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

// IEEE 754 infinity constants as raw uint bit patterns.
// Used by LOG(0).  The reference C emulator (nv2a_vsh_cpu) returns
// -INFINITY for this case, and xemu does the same.  However the reference
// code carries a "TODO: Validate this on HW" comment — real NV2A silicon
// (Kelvin-class, 2001) may clamp to a large finite value instead of
// producing a true IEEE infinity.  We match the existing emulator consensus
// for now; if hardware tests reveal different behaviour, replacing this
// with a large negative float (e.g. -FLT_MAX / -3.4e38) would be the fix.
// Note: FXC rejects literal division by zero (-1.0f/0.0f), so we use
// asfloat() on the raw IEEE 754 bit patterns instead.
static const float CXBX_POS_INF = asfloat(0x7F800000u);
static const float CXBX_NEG_INF = asfloat(0xFF800000u);

// ============================================================
// NV2A-accurate multiply: 0 * anything = 0, even 0 * inf
// Additionally, non-inf * non-inf cannot produce inf (clamp to ±FLT_MAX).
// Matches fix_inf_mult() in nv2a_vsh_cpu.c.
//
// Implementation: fully branchless.  The zero rule fires far more often
// than the overflow rule in real shaders (every W-component transform
// hits 0*something).  The overflow path is nearly free since GPUs execute
// both sides of a ternary anyway (movc).
// ============================================================
float4 nv2a_mul(float4 a, float4 b)
{
    float4 result = a * b;
    // Rule 1: 0 * anything = 0 (overrides NaN from 0*inf)
    // select() compiles to movc — no divergence.
    bool4 zeroMask = (a == 0.0f) | (b == 0.0f);
    // Rule 2: finite * finite cannot overflow to inf
    // Combine with zero mask: if zero, result is 0 regardless of overflow.
    uint4 bits = asuint(result);
    float4 clamped = asfloat((bits & 0xFF000000u) | 0x7FFFFFu);
    bool4 overflow = isinf(result) & !isinf(a) & !isinf(b);
    // Single select chain: zero wins over overflow wins over raw.
    result = overflow ? clamped : result;
    result = zeroMask ? (float4)0.0f : result;
    return result;
}

float3 nv2a_mul3(float3 a, float3 b)
{
    float3 result = a * b;
    bool3 zeroMask = (a == 0.0f) | (b == 0.0f);
    uint3 bits = asuint(result);
    float3 clamped = asfloat((bits & 0xFF000000u) | 0x7FFFFFu);
    bool3 overflow = isinf(result) & !isinf(a) & !isinf(b);
    result = overflow ? clamped : result;
    result = zeroMask ? (float3)0.0f : result;
    return result;
}

float nv2a_mul1(float a, float b)
{
    float result = a * b;
    if ((a == 0.0f) || (b == 0.0f)) return 0.0f;
    uint bits = asuint(result);
    return (isinf(result) && !isinf(a) && !isinf(b))
         ? asfloat((bits & 0xFF000000u) | 0x7FFFFFu)
         : result;
}

// ============================================================
// NV2A-accurate dot products using nv2a_mul per component.
// The sum itself also cannot overflow to inf (fix_inf in reference).
// Branchless: compute clamped value unconditionally, select via movc.
// ============================================================
float nv2a_dot3(float3 a, float3 b)
{
    float3 m = nv2a_mul3(a, b);
    float result = m.x + m.y + m.z;
    uint bits = asuint(result);
    float clamped = asfloat((bits & 0xFF000000u) | 0x7FFFFFu);
    return isinf(result) ? clamped : result;
}

float nv2a_dot4(float4 a, float4 b)
{
    float4 m = nv2a_mul(a, b);
    float result = m.x + m.y + m.z + m.w;
    uint bits = asuint(result);
    float clamped = asfloat((bits & 0xFF000000u) | 0x7FFFFFu);
    return isinf(result) ? clamped : result;
}

#endif // CXBX_NV2A_MATH_HELPERS_HLSLI
