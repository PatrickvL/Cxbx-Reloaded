// CxbxNV2AVshOps.hlsli
//
// NV2A vertex shader instruction implementations shared by JIT and interpreter.
// Both paths #include this file to guarantee identical behavior.
//
// Naming: mac_* for MAC unit ops, ilu_* for ILU unit ops.
// Internally these use the generic nv2a_* math helpers from CxbxNV2AMathHelpers.hlsli.
//
// Requires: CxbxNV2AMathHelpers.hlsli (nv2a_mul, nv2a_mul1, nv2a_dot3, nv2a_dot4, CXBX_NEG_INF)

#ifndef CXBX_NV2A_VSH_OPS_HLSLI
#define CXBX_NV2A_VSH_OPS_HLSLI

// ============================================================
// MAC (vector) unit operations
// ============================================================

// ARL: floor(src.x) with small bias to counter GPU float under-rounding
// of byte-normalised vertex attributes.
//
// When the Xbox CPU uploads a byte vertex attribute like 17, the NV2A
// normalises it to 17/255 via its fixed-function input unit (well-defined
// rounding).  GPU shader floats may represent this as slightly less than
// the true value (e.g. 16.999... instead of 17.0) after the shader
// multiplies back by 255, so a naive floor() would yield 16.
// Adding a small bias before floor() compensates for this.
//
// Origin: xqemu PR #79 "Add ARL-bias to work around OpenGL float behaviour"
//   https://github.com/xqemu/xqemu/pull/79
// Background: xqemu issue #78 "GLSL floats are not suitable for VS emulation"
//   https://github.com/xqemu/xqemu/issues/78
//
// Per the NV_vertex_program spec (section 2.14.1.11), the floor operations
// in ARL and EXP "must operate identically".  xqemu issue #105 notes that
// applying the bias to EXP's floor too would be the correct approach, but
// doing so risks breaking EXP's result.y fractional guarantee (expected in
// [0,1)).  We therefore apply the bias only to ARL (matching xemu) and
// leave EXP using exact floor() -- see ilu_exp() below.
//
// Known limitation: the bias can cause floor(N - epsilon) -> N when the
// true mathematical result should be N-1 (e.g. 16.999 -> 17).  This is
// considered less common than the byte-normalisation under-rounding it fixes.
int mac_arl(float4 a) { return (int)floor(a.x + 0.001); }

float4 mac_mov(float4 a) { return a; }
float4 mac_mul(float4 a, float4 b) { return nv2a_mul(a, b); }
float4 mac_add(float4 a, float4 c) { return a + c; }
float4 mac_mad(float4 a, float4 b, float4 c) { return nv2a_mul(a, b) + c; }
float4 mac_dp3(float4 a, float4 b) { return nv2a_dot3(a.xyz, b.xyz).xxxx; }
float4 mac_dph(float4 a, float4 b) { return (nv2a_dot3(a.xyz, b.xyz) + b.w).xxxx; }
float4 mac_dp4(float4 a, float4 b) { return nv2a_dot4(a, b).xxxx; }
float4 mac_dst(float4 a, float4 b) { return float4(1.0, nv2a_mul1(a.y, b.y), a.z, b.w); }
float4 mac_min(float4 a, float4 b) { return min(a, b); }
float4 mac_max(float4 a, float4 b) { return max(a, b); }

// NV2A SLT treats -0 < +0 as true (bit pattern comparison).
// Matches nv2a_less_than() in nv2a_vsh_cpu.c.
// Fully branchless: GPU `<` already handles normal cases; only the
// -0 vs +0 edge case needs a fixup via bitwise comparison.
float4 mac_slt(float4 a, float4 b)
{
    float4 result = float4(a < b);
    // IEEE says -0 == +0, but NV2A says -0 < +0.
    // Fixup: OR in 1.0 where a is -0 and b is +0.
    uint4 fixup = uint4(asuint(a) == 0x80000000u) & uint4(asuint(b) == 0u);
    return max(result, asfloat(fixup & asuint(1.0f)));
}

// SGE is defined as !SLT (1 - less_than) per nv2a_vsh_cpu.c.
float4 mac_sge(float4 a, float4 b)
{
    return 1.0f - mac_slt(a, b);
}

// ============================================================
// ILU (scalar) unit operations — operate on src.x
// ============================================================

// ILU MOV preserves all 4 components (the disassembler does NOT force
// .xxxx swizzle for MOV/LIT — only RCP/RCC/RSQ/EXP/LOG get that treatment).
float4 ilu_mov(float4 src) { return src; }

// RCP: 1/src.x, replicated.  HLSL division by zero produces ±inf per
// IEEE 754 (well-defined, unlike C).  The sign of zero is preserved via
// the sign bit of the division result.  No branches needed.
float4 ilu_rcp(float4 src)
{
    return (1.0 / src.x).xxxx;
}

// RCC: reciprocal with NV2A range clamping.
// After 1/x, clamp magnitude to [2^-64, 2^64] preserving sign.
// Matches xemu's clampAwayZeroInf and nv2a_vsh_cpu_rcc.
// Branchless: clamp abs, restore sign bit.
float4 ilu_rcc(float4 src)
{
    float rv = 1.0 / src.x;
    float av = abs(rv);
    // Clamp magnitude to [2^-64, 2^64] using exact IEEE bit patterns
    // 0x1F800000 = 2^-64, 0x5F800000 = 2^64
    av = clamp(av, asfloat(0x1F800000u), asfloat(0x5F800000u));
    // Restore original sign
    float result = asfloat(asuint(av) | (asuint(rv) & 0x80000000u));
    return result.xxxx;
}

// RSQ: 1/sqrt(|src.x|), replicated.  HLSL rsqrt(0) = +inf per spec.
float4 ilu_rsq(float4 src)
{
    return rsqrt(abs(src.x)).xxxx;
}

// EXP: {2^floor(src.x), frac(src.x), 2^src.x, 1.0}
// Uses exact floor() — NOT the ARL bias floor (see mac_arl comment above).
// Per NV_vertex_program spec (section 2.14.1.11), the floor operations in
// ARL and EXP "must operate identically".  xqemu issue #105 notes that
// applying the bias here too would be correct, but would break EXP's
// result.y fractional guarantee (expected in [0,1)).  We match xemu:
// ARL gets the bias; EXP does not.
float4 ilu_exp(float4 src)
{
    float s = src.x;
    float fl = floor(s);
    return float4(exp2(fl), s - fl, exp2(s), 1.0);
}

// LOG: {exponent, mantissa, log2(|src.x|), 1.0}
// Matches nv2a_vsh_cpu_log — uses IEEE bit extraction for .x and .y rather
// than log2-based decomposition (which loses precision for denormals).
//   .x = biased exponent (floor(log2(|src|)) for normal floats)
//   .y = mantissa in [1.0, 2.0) — raw IEEE mantissa bits | 0x3F800000
//   .z = log2(|src|) — full-precision intrinsic
//   .w = 1.0
// Special case: LOG(0) = (-inf, 1, -inf, 1).  See CXBX_NEG_INF for HW notes.
// Branchless: compute normal path (NaN when t==0 is harmless, selected away
// by the ternary → movc).  Edge cases selected via single ternary chain.
float4 ilu_log(float4 src)
{
    float t = abs(src.x);
    // Normal path: extract IEEE exponent and mantissa via bit ops
    uint t_bits = asuint(t);
    float exponent = (float)((int)((t_bits >> 23u) & 0xFFu) - 127);
    float mantissa = asfloat((t_bits & 0x7FFFFFu) | 0x3F800000u);
    float4 normal_result = float4(exponent, mantissa, log2(t), 1.0);
    // Edge cases: t==0 → (-inf, 1, -inf, 1); t==inf → (+inf, 1, +inf, 1)
    float4 zero_result = float4(CXBX_NEG_INF, 1.0f, CXBX_NEG_INF, 1.0f);
    float4 inf_result  = float4(CXBX_POS_INF, 1.0f, CXBX_POS_INF, 1.0f);
    return (t == 0.0f) ? zero_result : (isinf(t) ? inf_result : normal_result);
}

// LIT: matches nv2a_vsh_cpu_lit — kMax = 127.9961f.
// Branchless: conditional multiply pattern avoids divergence.
// Result: {1, diffuse>0 ? diffuse : 0, (diffuse>0 && blinn>0) ? pow(blinn, power) : 0, 1}
float4 ilu_lit(float4 src)
{
    static const float kMax = 127.9961f;
    float diffuse = src.x;
    float blinn = src.y;
    float specPower = clamp(src.w, -kMax, kMax);
    float dGt0 = (diffuse > 0.0f) ? 1.0f : 0.0f;
    float bGt0 = (blinn > 0.0f) ? 1.0f : 0.0f;
    // pow(blinn, specPower) — guard blinn with max to avoid pow(0/neg, x) domain issues
    float litZ = dGt0 * bGt0 * pow(max(blinn, 0.000001f), specPower);
    return float4(1.0f, dGt0 * diffuse, litZ, 1.0f);
}

#endif // CXBX_NV2A_VSH_OPS_HLSLI
