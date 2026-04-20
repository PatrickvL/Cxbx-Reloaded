// RegisterCombinerInterpreter.hlsl
//
// DX11 / SM5 Xbox NV2A register combiner interpreter ubershader.
//
// Optimizations applied (relative to the SM3.0 original):
//   Step 2 – Register file is a flat float4[16] array; all switch-based
//             register dispatch is replaced by direct indexed reads/writes.
//   Step 3 – NUM_STAGES compile-time specialization: define NUM_STAGES=N
//             (1..8) when compiling to produce a variant whose combiner loop
//             unrolls to exactly N iterations, letting the compiler eliminate
//             dead stage code.  Default (NUM_STAGES=8) is the safe fallback.
//   Step 4 – RGB and alpha combiner sub-stages are merged into one function
//             call per stage, halving call overhead and letting the compiler
//             schedule input fetches for both paths together.
//   All fmod/floor bit-extraction replaced with native bitwise operators.
//
// Bug fixes applied (relative to the SM3.0 original):
//   – Texture register written to the correct T[stage] slot.
//     (Was always writing T3 due to max() used instead of direct index.)
//   – R0.rgb initialised from T0.rgb, not white.
//     (NV2A spec: R0 starts as a copy of T0.)
//   – PS_CHANNEL_RGB input passes full rgba through intact.
//     (Was replacing .a with .b via .rgbb swizzle.)
//   – CLIPPLANE discards without sampling the texture.
//     (Was calling Sample2D after the discard test.)
//   – DOT_RFLCT_DIFF uses the constructed normal directly as the cubemap
//     direction; no erroneous reflect() computation.
//
// Host-side packing change from the SM3.0 version:
//   All 32-bit Xbox DWORDs (PSRGBInputs, PSAlphaOutputs, etc.) are now
//   passed as raw uint in the cbuffer rather than byte-split float4.
//   The C++ side must be updated to write these as uint.
//   Color constants (PSConstant0/1, FogColor, …) remain float4 [0..1].

// ------------------------------------------------------------
// Step 3: compile-time stage count specialization
// ------------------------------------------------------------
#ifndef NUM_STAGES
#define NUM_STAGES 8        // Safe default; override with 1..8 per variant
#endif

// ============================================================
// Shared constants and cbuffer layout (defined once in .hlsli headers,
// shared with C++ backend code)
// ============================================================
#include "NV2APixelShaderConstants.hlsli"
#include "RegisterCombinerInterpreterState.hlsli"

// ============================================================
// Textures and samplers
// Individual declarations per stage avoid DX11 error X4539
// (sampler array elements cannot be used with different intrinsics)
// ============================================================

Texture2D    Tex2D_0   : register(t0);
Texture2D    Tex2D_1   : register(t1);
Texture2D    Tex2D_2   : register(t2);
Texture2D    Tex2D_3   : register(t3);
Texture3D    Tex3D_0   : register(t4);
Texture3D    Tex3D_1   : register(t5);
Texture3D    Tex3D_2   : register(t6);
Texture3D    Tex3D_3   : register(t7);
TextureCube  TexCube_0 : register(t8);
TextureCube  TexCube_1 : register(t9);
TextureCube  TexCube_2 : register(t10);
TextureCube  TexCube_3 : register(t11);
SamplerState Samp0     : register(s0);
SamplerState Samp1     : register(s1);
SamplerState Samp2     : register(s2);
SamplerState Samp3     : register(s3);

// ============================================================
// Pixel shader input
// Note: CxbxPixelShaderHelpers.hlsli has a cross-API version of this struct,
// but we can't include it here because it also declares D3D9-style samplers
// that conflict with our D3D11-style Texture2D/SamplerState declarations.
// ============================================================
struct PS_INPUT
{
    float4 iPos : SV_Position;
    float4 iD0  : COLOR0;       // Diffuse  (front-facing)
    float4 iD1  : COLOR1;       // Specular (front-facing)
    float  iFog : FOG;
    float  iPts : PSIZE;
    float4 iB0  : TEXCOORD4;    // Diffuse  (back-facing)
    float4 iB1  : TEXCOORD5;    // Specular (back-facing)
    float4 iT0  : TEXCOORD0;
    float4 iT1  : TEXCOORD1;
    float4 iT2  : TEXCOORD2;
    float4 iT3  : TEXCOORD3;
    bool   iFF  : SV_IsFrontFace;
};

// ============================================================
// Step 2: register file
//
// All 16 PS_REGISTER_* indices mapped to float4 Regs[16].
// Direct indexing replaces the original switch chains entirely.
//
// Slot ownership:
//   0  ZERO/DISCARD  read-only zero; writes silently dropped
//   1  C0            resolved from cbuffer; not stored here
//   2  C1            resolved from cbuffer; not stored here
//   3  FOG           read/write
//   4  V0            read/write
//   5  V1            read/write
//   6, 7             reserved; writes silently dropped
//   8..11  T0..T3    written by FetchTexture, then read/write
//   12..13 R0..R1    read/write
//   14 V1R0_SUM      written in final combiner EFG phase, read in ABCD
//   15 EF_PROD       written in final combiner EFG phase, read in ABCD
// ============================================================

float4 RegRead(float4 Regs[16], uint idx)
{
    return Regs[idx & 0xFu];
}

void RegWrite(inout float4 Regs[16], uint idx, float4 val)
{
    // Silently drop writes to read-only and reserved slots
    if (idx == PS_REGISTER_ZERO ||
        idx == PS_REGISTER_C0   ||
        idx == PS_REGISTER_C1   ||
        idx == 6u               ||
        idx == 7u)
        return;
    Regs[idx & 0xFu] = val;
}

// Modify only the .rgb channels of a register; preserve .a
void RegWriteRGB(inout float4 Regs[16], uint idx, float3 rgb)
{
    if (idx == PS_REGISTER_DISCARD) return;
    RegWrite(Regs, idx, float4(rgb, RegRead(Regs, idx).a));
}

// Modify only the .a channel of a register; preserve .rgb
void RegWriteA(inout float4 Regs[16], uint idx, float a)
{
    if (idx == PS_REGISTER_DISCARD) return;
    RegWrite(Regs, idx, float4(RegRead(Regs, idx).rgb, a));
}

// ============================================================
// Input mapping
// ============================================================

float4 ApplyInputMapping(uint mapping, float4 v)
{
    switch (mapping & 0xE0u) // isolate the 3 mapping bits
    {
        case PS_INPUTMAPPING_UNSIGNED_IDENTITY: return max(0.0f, v);
        case PS_INPUTMAPPING_UNSIGNED_INVERT:   return 1.0f - clamp(v, 0.0f, 1.0f);
        case PS_INPUTMAPPING_EXPAND_NORMAL:     return  2.0f * max(0.0f, v) - 1.0f;
        case PS_INPUTMAPPING_EXPAND_NEGATE:     return -2.0f * max(0.0f, v) + 1.0f;
        case PS_INPUTMAPPING_HALFBIAS_NORMAL:   return        max(0.0f, v) - 0.5f;
        case PS_INPUTMAPPING_HALFBIAS_NEGATE:   return  0.5f - max(0.0f, v);
        case PS_INPUTMAPPING_SIGNED_IDENTITY:   return v;
        default: /* PS_INPUTMAPPING_SIGNED_NEGATE */ return -v;
    }
}

// ============================================================
// Output mapping
//
// flags = PS_COMBINEROUTPUT flags field (rgbOut >> 12 or aOut >> 12).
// Bits [5:3] encode the output mapping:
//   bit 3 = PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS: subtract 0.5 before scaling
//   bits [5:4] = scale: 00=x1  01=x2  10=x4  11=/2
// ============================================================

float4 ApplyOutputMapping(uint flags, float4 v)
{
    float bias  = (flags & PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS) ? -0.5f : 0.0f;
    float scale;
    switch ((flags >> 4u) & 3u) // scale bits are [5:4] of flags
    {
        case 1u:  scale = 2.0f;  break;
        case 2u:  scale = 4.0f;  break;
        case 3u:  scale = 0.5f;  break;
        default:  scale = 1.0f;  break;
    }
    return clamp((v + bias) * scale, -1.0f, 1.0f);
}

// ============================================================
// Input register resolution
//
// stageIdx:  0..7 = color combiner stage
//            8    = final combiner EFG inputs
//            9    = final combiner ABCD inputs
// ============================================================

float4 ResolveInput(float4 Regs[16], uint regByte, bool isAlpha, uint stageIdx,
                    bool flagUniqueC0, bool flagUniqueC1)
{
    uint regIdx    = regByte & 0x0Fu;
    uint mapping   = regByte & 0xE0u;
    bool useAlphaC = (regByte & PS_CHANNEL_ALPHA) != 0u;
    bool isFinal   = (stageIdx >= 8u);
    bool isFinalAB = (stageIdx == 9u);

    // Clamp stage index to [0,7] for PSConstant0/1 array access.
    // The isFinal guard prevents out-of-range reads at runtime, but the
    // HLSL compiler cannot prove this statically (X3504).
    uint safeStage = min(stageIdx, 7u);

    float4 val;
    switch (regIdx)
    {
        case PS_REGISTER_C0:
            val = isFinal    ? PSFinalCombinerConstant[0]
                : flagUniqueC0 ? PSConstant0[safeStage]
                               : PSConstant0[0];
            break;

        case PS_REGISTER_C1:
            val = isFinal    ? PSFinalCombinerConstant[1]
                : flagUniqueC1 ? PSConstant1[safeStage]
                               : PSConstant1[0];
            break;

        case PS_REGISTER_FOG:
        {
            float4 fog = RegRead(Regs, PS_REGISTER_FOG);
            // Color stages may read FOG.rgb only; final combiner reads FOG.a only
            val = isFinal ? float4(0.0f, 0.0f, 0.0f, fog.a)
                          : float4(fog.rgb, 1.0f);
            break;
        }

        case PS_REGISTER_V1R0_SUM:
        case PS_REGISTER_EF_PROD:
            // These are only valid as ABCD inputs to the final combiner;
            // reading them elsewhere returns zero
            val = isFinalAB ? RegRead(Regs, regIdx) : (float4)0.0f;
            break;

        default:
            val = RegRead(Regs, regIdx);
            break;
    }

    // Channel select
    // Bug fix: PS_CHANNEL_RGB must pass the full float4 through intact.
    // Combiner operations then use .rgb for RGB work and .a for alpha work.
    // The SM3.0 version incorrectly used .rgbb, clobbering alpha with blue.
    if (useAlphaC)
        val = val.aaaa;

    // Remap final-combiner-invalid mappings to their nearest valid equivalents
    if (isFinal)
    {
        switch (mapping)
        {
            case PS_INPUTMAPPING_EXPAND_NORMAL:
            case PS_INPUTMAPPING_HALFBIAS_NORMAL:
            case PS_INPUTMAPPING_SIGNED_IDENTITY:
                mapping = PS_INPUTMAPPING_UNSIGNED_IDENTITY;
                break;
            case PS_INPUTMAPPING_EXPAND_NEGATE:
            case PS_INPUTMAPPING_HALFBIAS_NEGATE:
            case PS_INPUTMAPPING_SIGNED_NEGATE:
                mapping = PS_INPUTMAPPING_UNSIGNED_INVERT;
                break;
        }
    }

    return ApplyInputMapping(mapping, val);
}

// ============================================================
// Color sign conversion (Xbox X_D3DTSS_COLORSIGN extension)
// ============================================================

float4 ApplyColorSign(float4 sign, float4 t)
{
    if (sign.r > 0.0f) t.r = t.r * 2.0f - 1.0f; else if (sign.r < 0.0f) t.r = t.r * 0.5f + 0.5f;
    if (sign.g > 0.0f) t.g = t.g * 2.0f - 1.0f; else if (sign.g < 0.0f) t.g = t.g * 0.5f + 0.5f;
    if (sign.b > 0.0f) t.b = t.b * 2.0f - 1.0f; else if (sign.b < 0.0f) t.b = t.b * 0.5f + 0.5f;
    if (sign.a > 0.0f) t.a = t.a * 2.0f - 1.0f; else if (sign.a < 0.0f) t.a = t.a * 0.5f + 0.5f;
    return t;
}

// ============================================================
// Texture sampling helpers (individual to avoid X4539)
// ============================================================

float4 Sample2D  (uint s, float2 uv)
{
    switch (s) {
        case 0:  return Tex2D_0.Sample(Samp0, uv);
        case 1:  return Tex2D_1.Sample(Samp1, uv);
        case 2:  return Tex2D_2.Sample(Samp2, uv);
        default: return Tex2D_3.Sample(Samp3, uv);
    }
}

float4 Sample3D  (uint s, float3 uvw)
{
    switch (s) {
        case 0:  return Tex3D_0.Sample(Samp0, uvw);
        case 1:  return Tex3D_1.Sample(Samp1, uvw);
        case 2:  return Tex3D_2.Sample(Samp2, uvw);
        default: return Tex3D_3.Sample(Samp3, uvw);
    }
}

float4 SampleCube(uint s, float3 dir)
{
    switch (s) {
        case 0:  return TexCube_0.Sample(Samp0, dir);
        case 1:  return TexCube_1.Sample(Samp1, dir);
        case 2:  return TexCube_2.Sample(Samp2, dir);
        default: return TexCube_3.Sample(Samp3, dir);
    }
}

// ============================================================
// Texture stage fetch
// ============================================================

void FetchTexture(inout float4 Regs[16], uint stage, float4 coords, uint mode)
{
    float4 val = float4(0.0f, 0.0f, 0.0f, 1.0f); // default: opaque black

    switch (mode)
    {
    case PS_TEXTUREMODES_NONE:
        // Leave T[stage] at its initialised value (0,0,0,1); no write
        return;

    case PS_TEXTUREMODES_PROJECT2D:
        val = Sample2D(stage, coords.xy / coords.w);
        break;

    case PS_TEXTUREMODES_PROJECT3D:
        val = Sample3D(stage, coords.xyz / coords.w);
        break;

    case PS_TEXTUREMODES_CUBEMAP:
        val = SampleCube(stage, coords.xyz);
        break;

    case PS_TEXTUREMODES_PASSTHRU:
        val = saturate(coords);
        break;

    case PS_TEXTUREMODES_CLIPPLANE:
        // Discard pixel if any coordinate is negative.
        // Bug fix: T[stage] stays (0,0,0,1) — no sampling after the test.
        if (coords.x < 0.0f || coords.y < 0.0f ||
            coords.z < 0.0f || coords.w < 0.0f)
            discard;
        // val stays (0,0,0,1)
        break;

    case PS_TEXTUREMODES_BUMPENVMAP:
    case PS_TEXTUREMODES_BUMPENVMAP_LUM:
        // Sample the bump source texture; perturbation of the next stage's
        // coordinates is handled by vertex shader / texture coordinate routing
        // (TODO: full perturbation pass)
        val = Sample2D(stage, coords.xy);
        break;

    case PS_TEXTUREMODES_DPNDNT_AR:
    {
        // Dependent lookup: use .a and .r from T[stage-1] as (u,v)
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        val = Sample2D(stage, prev.ar);
        break;
    }

    case PS_TEXTUREMODES_DPNDNT_GB:
    {
        // Dependent lookup: use .g and .b from T[stage-1] as (u,v)
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        val = Sample2D(stage, prev.gb);
        break;
    }

    case PS_TEXTUREMODES_DOTPRODUCT:
    {
        // Compute dot(tex_coords, T[stage-1].rgb); store scalar in .x
        // Subsequent DOT_ST / DOT_ZW / DOT_STR modes read T[stage].x
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  d    = dot(coords.xyz, prev.xyz);
        val = float4(d, 0.0f, 0.0f, 1.0f);
        break;
    }

    case PS_TEXTUREMODES_DOT_ST:
    {
        // Use dot results from T[stage-2].x and T[stage-1].x as (s,t)
        float s = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float t = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u)).x;
        val = Sample2D(stage, float2(s, t));
        break;
    }

    case PS_TEXTUREMODES_DOT_ZW:
    {
        // (stage-2 dot result, current dot result) stored as (z, w)
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  s    = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float  t    = dot(coords.xyz, prev.xyz);
        val = float4(0.0f, 0.0f, s, t);
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_DIFF:
    {
        // Build a normal from three sequential DOTPRODUCT results; cubemap lookup.
        // Bug fix: normal is used directly as the cubemap direction.
        // There is no reflection computation for the DIFF mode (that is SPEC).
        float  nx   = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float  ny   = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u)).x;
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  nz   = dot(coords.xyz, prev.xyz);
        val = SampleCube(stage, float3(nx, ny, nz));
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_SPEC:
    {
        // Build normal, compute specular reflection vector, cubemap lookup.
        // Eye vector is assembled from the .w component of T1, T2, T3
        // (NV2A hardcodes these indices).
        float  nx   = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float  ny   = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u)).x;
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  nz   = dot(coords.xyz, prev.xyz);
        float3 N    = normalize(float3(nx, ny, nz));
        float3 E    = normalize(float3(RegRead(Regs, PS_REGISTER_T1).w,
                                       RegRead(Regs, PS_REGISTER_T2).w,
                                       RegRead(Regs, PS_REGISTER_T3).w));
        float3 R    = 2.0f * dot(N, E) * N - E;
        val = SampleCube(stage, R);
        break;
    }

    case PS_TEXTUREMODES_DOT_STR_3D:
    {
        float  s    = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float  t    = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u)).x;
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  r    = dot(coords.xyz, prev.xyz);
        val = Sample3D(stage, float3(s, t, r));
        break;
    }

    case PS_TEXTUREMODES_DOT_STR_CUBE:
    {
        float  s    = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float  t    = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u)).x;
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  r    = dot(coords.xyz, prev.xyz);
        val = SampleCube(stage, float3(s, t, r));
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST:
    {
        // Like DOT_RFLCT_SPEC but eye vector is a constant.
        // TODO: wire SetEyeVector() / D3DRS_PSINPUTTEXTURE to a cbuffer entry
        float  nx   = RegRead(Regs, PS_REGISTER_T0 + (stage - 2u)).x;
        float  ny   = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u)).x;
        float4 prev = RegRead(Regs, PS_REGISTER_T0 + (stage - 1u));
        float  nz   = dot(coords.xyz, prev.xyz);
        float3 N    = normalize(float3(nx, ny, nz));
        float3 E    = float3(0.0f, 0.0f, 1.0f); // placeholder
        float3 R    = 2.0f * dot(N, E) * N - E;
        val = SampleCube(stage, R);
        break;
    }

    case PS_TEXTUREMODES_BRDF:
        // TODO: proper BRDF requires eye / light sigma inputs
        val = Sample2D(stage, coords.xy);
        break;

    default:
        val = Sample2D(stage, coords.xy);
        break;
    }

    // Apply COLORSIGN fixup, then write to the correct T register.
    // Bug fix: direct indexed write replaces the broken max() clamp
    // that always directed all four stages to write T3.
    val = ApplyColorSign(ColorSign[stage], val);
    RegWrite(Regs, PS_REGISTER_T0 + stage, val);
}

// ============================================================
// Step 4: combined RGB + alpha combiner stage
//
// The SM3.0 version called do_color_combiner_stage() twice per stage
// (once with is_alpha=false, once with is_alpha=true), decoding output
// bytes redundantly in each call.
//
// This version decodes all per-stage state once, fetches all eight inputs
// in a single block (letting the compiler schedule register reads together),
// and computes both RGB and alpha results before any register writes.
//
// Write ordering: RGB writes precede alpha writes.  When the same register
// appears in both paths (the common case, e.g. AB → R0 and alphaAB → R0),
// the alpha write reads the updated RGB and overlays only the .a component,
// which is the correct NV2A behaviour.
// ============================================================

void DoCombinerStage(inout float4 Regs[16], uint stage,
                     bool flagMuxMsb, bool flagUniqueC0, bool flagUniqueC1)
{
    uint rgbIn  = PSRGBInputs[stage];
    uint aIn    = PSAlphaInputs[stage];
    uint rgbOut = PSRGBOutputs[stage];
    uint aOut   = PSAlphaOutputs[stage];

    // --- Decode output control bits ---
    uint rgbFlags = rgbOut >> PS_COMBINEROUTPUTS_FLAGS_SHIFT;
    uint aFlags   = aOut   >> PS_COMBINEROUTPUTS_FLAGS_SHIFT;
    bool flagCDDot  = (rgbFlags & PS_COMBINEROUTPUT_CD_DOT_PRODUCT)   != 0u;
    bool flagABDot  = (rgbFlags & PS_COMBINEROUTPUT_AB_DOT_PRODUCT)   != 0u;
    bool flagRGBMux = (rgbFlags & PS_COMBINEROUTPUT_AB_CD_MUX)        != 0u;
    bool flagAMux   = (aFlags   & PS_COMBINEROUTPUT_AB_CD_MUX)        != 0u;
    bool cdBlue2A   = (rgbFlags & PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA) != 0u;
    bool abBlue2A   = (rgbFlags & PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA) != 0u;

    // Output destination registers (DISCARD == 0 == no-op write)
    uint rgbRegCD  = (rgbOut >> PS_COMBINEROUTPUTS_CD_SHIFT)      & 0xFu;
    uint rgbRegAB  = (rgbOut >> PS_COMBINEROUTPUTS_AB_SHIFT)      & 0xFu;
    uint rgbRegSum = (rgbOut >> PS_COMBINEROUTPUTS_MUX_SUM_SHIFT) & 0xFu;
    uint aRegCD    = (aOut   >> PS_COMBINEROUTPUTS_CD_SHIFT)      & 0xFu;
    uint aRegAB    = (aOut   >> PS_COMBINEROUTPUTS_AB_SHIFT)      & 0xFu;
    uint aRegSum   = (aOut   >> PS_COMBINEROUTPUTS_MUX_SUM_SHIFT) & 0xFu;

    // --- Fetch all eight inputs in one block ---
    // PS_COMBINERINPUTS(A, B, C, D) packing; grouping reads lets the compiler schedule together
    float4 rgbA = ResolveInput(Regs, (rgbIn >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu, false, stage, flagUniqueC0, flagUniqueC1);
    float4 rgbB = ResolveInput(Regs, (rgbIn >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu, false, stage, flagUniqueC0, flagUniqueC1);
    float4 rgbC = ResolveInput(Regs, (rgbIn >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu, false, stage, flagUniqueC0, flagUniqueC1);
    float4 rgbD = ResolveInput(Regs, (rgbIn >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu, false, stage, flagUniqueC0, flagUniqueC1);
    float   aA  = ResolveInput(Regs, (aIn   >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu, true,  stage, flagUniqueC0, flagUniqueC1).a;
    float   aB  = ResolveInput(Regs, (aIn   >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu, true,  stage, flagUniqueC0, flagUniqueC1).a;
    float   aC  = ResolveInput(Regs, (aIn   >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu, true,  stage, flagUniqueC0, flagUniqueC1).a;
    float   aD  = ResolveInput(Regs, (aIn   >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu, true,  stage, flagUniqueC0, flagUniqueC1).a;

    // --- Compute AB and CD products ---
    // Dot product applies to RGB only; alpha always multiplies scalars
    float3 rgbAB = flagABDot ? (float3)dot(rgbA.rgb, rgbB.rgb) : rgbA.rgb * rgbB.rgb;
    float3 rgbCD = flagCDDot ? (float3)dot(rgbC.rgb, rgbD.rgb) : rgbC.rgb * rgbD.rgb;
    float   aAB  = aA * aB;
    float   aCD  = aC * aD;

    // --- SUM or MUX ---
    // MUX reads R0.a; MSB mode thresholds at 0.5, LSB mode checks the integer bit
    float r0a = Regs[PS_REGISTER_R0].a;
    bool muxSel = flagMuxMsb
        ? (r0a >= 0.5f)
        : (((uint)(r0a * 255.0f + 0.5f) & 1u) != 0u);

    float3 rgbABCD = flagRGBMux ? (muxSel ? rgbCD : rgbAB) : (rgbAB + rgbCD);
    float   aABCD  = flagAMux   ? (muxSel ? aCD   : aAB  ) : (aAB   + aCD  );

    // --- Apply output mapping to all six results ---
    float3 outRGB_AB  = ApplyOutputMapping(rgbFlags, float4(rgbAB,   0.0f)).rgb;
    float3 outRGB_CD  = ApplyOutputMapping(rgbFlags, float4(rgbCD,   0.0f)).rgb;
    float3 outRGB_Sum = ApplyOutputMapping(rgbFlags, float4(rgbABCD, 0.0f)).rgb;
    float  outA_AB    = ApplyOutputMapping(aFlags,   float4(0.0f, 0.0f, 0.0f, aAB  )).a;
    float  outA_CD    = ApplyOutputMapping(aFlags,   float4(0.0f, 0.0f, 0.0f, aCD  )).a;
    float  outA_Sum   = ApplyOutputMapping(aFlags,   float4(0.0f, 0.0f, 0.0f, aABCD)).a;

    // --- RGB writes ---
    // AB: write RGB; optionally propagate .b to .a (BlueToAlpha)
    if (rgbRegAB != PS_REGISTER_DISCARD) {
        float a = abBlue2A ? outRGB_AB.b : RegRead(Regs, rgbRegAB).a;
        RegWrite(Regs, rgbRegAB, float4(outRGB_AB, a));
    }
    // CD: same pattern
    if (rgbRegCD != PS_REGISTER_DISCARD) {
        float a = cdBlue2A ? outRGB_CD.b : RegRead(Regs, rgbRegCD).a;
        RegWrite(Regs, rgbRegCD, float4(outRGB_CD, a));
    }
    // AB+CD sum/mux: write RGB only; spec requires DISCARD when any DOT flag is set
    if (rgbRegSum != PS_REGISTER_DISCARD && !flagABDot && !flagCDDot)
        RegWriteRGB(Regs, rgbRegSum, outRGB_Sum);

    // --- Alpha writes (after RGB writes; see ordering note above) ---
    RegWriteA(Regs, aRegAB,  outA_AB);
    RegWriteA(Regs, aRegCD,  outA_CD);
    RegWriteA(Regs, aRegSum, outA_Sum);
}

// ============================================================
// Final combiner
// ============================================================

float4 DoFinalCombiner(inout float4 Regs[16], bool flagUniqueC0, bool flagUniqueC1)
{
    float4 R0 = RegRead(Regs, PS_REGISTER_R0);

    // If both ABCD and EFG are zero the final combiner is unused
    if (PSFinalCombinerInputsABCD == 0u && PSFinalCombinerInputsEFG == 0u)
        return R0;

    // PSFinalCombinerInputsEFG = PS_COMBINERINPUTS(E, F, G, settings)
    //   = (E<<24) | (F<<16) | (G<<8) | settings
    uint efg      = PSFinalCombinerInputsEFG;
    // PS_COMBINERINPUTS(E, F, G, settings) for the final combiner
    uint settings = (efg >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu;
    uint eReg     = (efg >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu;
    uint fReg     = (efg >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu;
    uint gReg     = (efg >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu;

    // --- Resolve E, F (RGB) and G (alpha) — stageIdx 8 = final combiner EFG ---
    float3 E = ResolveInput(Regs, eReg, false, 8u, flagUniqueC0, flagUniqueC1).rgb;
    float3 F = ResolveInput(Regs, fReg, false, 8u, flagUniqueC0, flagUniqueC1).rgb;
    float  G = ResolveInput(Regs, gReg, true,  8u, flagUniqueC0, flagUniqueC1).a;

    // Compute E*F and store in EF_PROD for potential use by ABCD inputs
    RegWrite(Regs, PS_REGISTER_EF_PROD, float4(E * F, 1.0f));

    // --- Optional complement and clamp on V1 and R0 ---
    // These modify the sum inputs, not the stored register values
    float3 v1 = RegRead(Regs, PS_REGISTER_V1).rgb;
    float3 r0 = R0.rgb;
    if (settings & PS_FINALCOMBINERSETTING_COMPLEMENT_V1) v1 = 1.0f - v1;
    if (settings & PS_FINALCOMBINERSETTING_COMPLEMENT_R0) r0 = 1.0f - r0;

    float3 v1r0sum = v1 + r0;
    if (settings & PS_FINALCOMBINERSETTING_CLAMP_SUM) v1r0sum = saturate(v1r0sum);

    // Store V1+R0 sum for potential use by ABCD inputs
    RegWrite(Regs, PS_REGISTER_V1R0_SUM, float4(v1r0sum, 1.0f));

    // --- Resolve A, B, C, D — stageIdx 9 = final combiner ABCD ---
    // V1R0_SUM and EF_PROD are now valid for reading at this stage
    uint abcd = PSFinalCombinerInputsABCD;
    float4 A = ResolveInput(Regs, (abcd >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu, false, 9u, flagUniqueC0, flagUniqueC1);
    float4 B = ResolveInput(Regs, (abcd >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu, false, 9u, flagUniqueC0, flagUniqueC1);
    float4 C = ResolveInput(Regs, (abcd >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu, false, 9u, flagUniqueC0, flagUniqueC1);
    float4 D = ResolveInput(Regs, (abcd >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu, false, 9u, flagUniqueC0, flagUniqueC1);

    // Final RGB = A*B + (1-A)*C + D, clamped to [0,1]
    // Final alpha = G
    float4 result;
    result.rgb = min(lerp(C.rgb, B.rgb, A.rgb) + D.rgb, 1.0f);
    result.a   = G;
    return result;
}

// ============================================================
// Pixel shader entry point
// ============================================================

float4 main(PS_INPUT input) : SV_Target
{
    // --- Decode PSCombinerCount ---
    // PSCombinerCount = PS_COMBINERCOUNT(count, flags) = (flags<<8) | count
    uint numStages    = clamp(PSCombinerCount & 0xFFu, 1u, 8u);
    uint ccFlags      = PSCombinerCount >> 8u; // extract flags portion
    bool flagMuxMsb   = (ccFlags & PS_COMBINERCOUNT_MUX_MSB)   != 0u;
    bool flagUniqueC0 = (ccFlags & PS_COMBINERCOUNT_UNIQUE_C0) != 0u;
    bool flagUniqueC1 = (ccFlags & PS_COMBINERCOUNT_UNIQUE_C1) != 0u;

    // --- Decode PSTextureModes ---
    // Four 5-bit fields packed sequentially; stage N occupies bits[N*5+4 : N*5]
    uint texMode[4];
    texMode[0] = (PSTextureModes      ) & 0x1Fu;
    texMode[1] = (PSTextureModes >>  5) & 0x1Fu;
    texMode[2] = (PSTextureModes >> 10) & 0x1Fu;
    texMode[3] = (PSTextureModes >> 15) & 0x1Fu;

    // --- Initialise the register file ---
    // Zero-init sets ZERO (index 0) permanently to 0.
    // T registers are explicitly set to (0,0,0,1) as the NV2A default.
    float4 Regs[16];
    [unroll] for (uint i = 0u; i < 16u; i++) Regs[i] = 0.0f;
    Regs[PS_REGISTER_T0] = float4(0.0f, 0.0f, 0.0f, 1.0f);
    Regs[PS_REGISTER_T1] = float4(0.0f, 0.0f, 0.0f, 1.0f);
    Regs[PS_REGISTER_T2] = float4(0.0f, 0.0f, 0.0f, 1.0f);
    Regs[PS_REGISTER_T3] = float4(0.0f, 0.0f, 0.0f, 1.0f);

    // --- Texture stages (must run sequentially; later stages can read earlier T regs) ---
    float4 texCoords[4] = { input.iT0, input.iT1, input.iT2, input.iT3 };
    [unroll] for (uint ts = 0u; ts < 4u; ts++)
        FetchTexture(Regs, ts, texCoords[ts], texMode[ts]);

    // --- Set vertex-derived registers ---
    bool  isFront = input.iFF;
    float4 diffuse  = isFront ? input.iD0 : input.iB0;
    float4 specular = isFront ? input.iD1 : input.iB1;
    RegWrite(Regs, PS_REGISTER_V0, diffuse);
    RegWrite(Regs, PS_REGISTER_V1, specular);
    // FOG: rgb from the fog color constant, alpha from the vertex fog factor
    RegWrite(Regs, PS_REGISTER_FOG, float4(FogColor.rgb, saturate(input.iFog)));

    // R0: Bug fix – initialise from T0 (not white).
    // NV2A spec: R0 starts as a full copy of T0; R0.a is T0.a.
    RegWrite(Regs, PS_REGISTER_R0, RegRead(Regs, PS_REGISTER_T0));
    // R1 remains (0,0,0,0) from the zero-init above

    // --- Step 3: combiner stage loop ---
    // [unroll] with a compile-time NUM_STAGES constant tells the compiler the
    // exact iteration count, which for specialized variants (NUM_STAGES < 8)
    // eliminates the unreachable tail stages entirely.
    // The runtime guard `stage < numStages` provides safety for the default
    // NUM_STAGES=8 variant when fewer stages are actually active.
    [unroll]
    for (uint stage = 0u; stage < (uint)NUM_STAGES; stage++)
    {
        if (stage < numStages)
            DoCombinerStage(Regs, stage, flagMuxMsb, flagUniqueC0, flagUniqueC1);
    }

    return DoFinalCombiner(Regs, flagUniqueC0, flagUniqueC1);
}
