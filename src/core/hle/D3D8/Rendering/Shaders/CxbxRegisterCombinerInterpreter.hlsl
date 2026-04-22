// CxbxRegisterCombinerInterpreter.hlsl
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
//   Step 3b– NUM_TEXTURE_STAGES compile-time specialization (1..4).
//             Single-texture shaders avoid instantiating FetchTexture 4×.
//   Step 4 – RGB and alpha combiner sub-stages are merged into one function
//             call per stage, halving call overhead and letting the compiler
//             schedule input fetches for both paths together.
//   Step 5 – ResolveInput split into ResolveStageInput / ResolveFinalInput,
//             eliminating isFinal branches from the hot stage-loop path.
//   Step 6 – Output mapping bias/scale hoisted: decoded once per stage,
//             applied inline to all six results (removes 6 switches/stage).
//   Step 7 – Vectorized ApplyColorSign (movc pipeline, no per-channel branches).
//   Step 8 – MUX selector guarded: R0.a decode skipped when no MUX flag set.
//             Sum computation skipped when both DOT flags set.
//   Step 9 – ResolveStageInputAlpha: scalar alpha resolver reads only .a
//             from the register file and calls ApplyInputMappingScalar,
//             eliminating ~12 wasted float ops per combiner stage.
//   Step 10– ApplyInputMappingScalar: scalar float counterpart to
//             ApplyInputMapping, used exclusively by ResolveStageInputAlpha.
//   Step 11– round() replaces floor(x + 0.5f) in ApplyDotMapping and
//             the MUX R0.a decode — single GPU instruction, same semantics.
//   Step 12– tBase in FetchTexture: PS_REGISTER_T0 + stage extracted once,
//             simplifying multi-stage dot-mode offsets (tBase - 1u, etc.).
//   Step 13– ResolveStageInput/Alpha switch eliminated: Regs[14]/[15] are
//             zero-initialized and never written until DoFinalCombiner, so
//             V1R0_SUM/EF_PROD already read as zero; FOG falls through to
//             the default case.  Both functions collapse to a direct array
//             read, removing a branch from the hottest path (8×/stage).
//   Step 14– DecodeOutputScale via exp2: scale values are powers of 2
//             (1, 2, 4, 0.5); exponent = ((m+1)&3)-1.  Single SFU op
//             replaces a 4-arm switch / 3-movc chain.
//   All fmod/floor bit-extraction replaced with native bitwise operators.
//
// Bug fixes applied (relative to the SM3.0 original):
//   – Texture register written to the correct T[stage] slot.
//     (Was always writing T3 due to max() used instead of direct index.)
//   – R0.a initialised from T0.a; R0.rgb starts at zero.
//     (NV2A spec: only R0.a copies from T0, rgb is zero.)
//   – PS_CHANNEL_RGB input passes full rgba through intact.
//     (Was replacing .a with .b via .rgbb swizzle.)
//   – CLIPPLANE discards without sampling the texture.
//     (Was calling Sample2D after the discard test.)
//   – DOT_RFLCT_DIFF uses the constructed normal directly as the cubemap
//     direction; no erroneous reflect() computation.
//   – FOG register in color stages preserves real fog.a for alpha channel
//     reads (was fabricating 1.0 via float4(fog.rgb, 1.0f) + .aaaa).
//   – nv2a_mul3 per-component zero propagation: || is scalar in SM5,
//     so (a == 0 || b == 0) collapsed bool3 to scalar via any(), zeroing
//     the entire result when any one component was zero.  Fixed with
//     bitwise | for per-component bool3 OR.
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

#ifndef NUM_TEXTURE_STAGES
#define NUM_TEXTURE_STAGES 4 // Safe default; override with 1..4 per variant
#endif

// ============================================================
// Shared constants and cbuffer layout (defined once in .hlsli headers,
// shared with C++ backend code)
// ============================================================
#include "CxbxNV2APixelShaderConstants.hlsli"
#include "CxbxRegisterCombinerInterpreterState.hlsli"
#include "CxbxNV2AMathHelpers.hlsli"

// Shared pure-math pixel shader helpers (ApplyTexFmtFixup, PerformColorSign,
// PerformColorKeyOp, PerformAlphaTest, CalculateFogFactor, etc.)
#include "CxbxPixelShaderFunctions.hlsli"

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

#include "CxbxPixelShaderInput.hlsli"

// ============================================================
// Step 2: register file
//
// All 16 PS_REGISTER_* indices mapped to float4 Regs[16].
// Direct indexing replaces the original switch chains entirely.
//
// Slot ownership:
//   0  ZERO/DISCARD  read-only zero; writes silently dropped
//   1  C0            initialised from cbuffer per-stage; read/write
//   2  C1            initialised from cbuffer per-stage; read/write
//   3  FOG           read/write
//   4  V0            read/write
//   5  V1            read/write
//   6, 7             reserved; writes silently dropped
//   8..11  T0..T3    written by FetchTexture, then read/write
//   12..13 R0..R1    read/write
//   14 V1R0_SUM      written in final combiner EFG phase, read in ABCD
//   15 EF_PROD       written in final combiner EFG phase, read in ABCD
// ============================================================

#define RegRead(Regs, idx) ((Regs)[(idx) & 0xFu])

// For callsites where the index is a literal constant known != 0
#define RegWriteDirect(Regs, idx, val) ((Regs)[(idx) & 0xFu] = (val))

// Silently drop writes to the read-only slot.
// The reserved slots (6,7) are also effectively read-only 
// since they are never initialised or read by the combiner logic;
// silently dropping writes to them matches NV2A behavior and
// avoids the need for explicit checks in the combiner code.
// C0/C1 ARE writable on NV2A hardware (confirmed by xemu).
// They are initialised from the cbuffer per-stage, but combiner
// output stages can overwrite them for subsequent reads.
void RegWrite(inout float4 Regs[16], uint idx, float4 val)
{
    uint i = idx & 0xFu;
    Regs[i] = (i == PS_REGISTER_ZERO) ? Regs[i] : val;
}

// Modify only the .rgb channels of a register; preserve .a
void RegWriteRGB(inout float4 Regs[16], uint idx, float3 rgb)
{
    uint i = idx & 0xFu;
    Regs[i].xyz = (i == PS_REGISTER_ZERO) ? Regs[i].xyz : rgb;
}

// Modify only the .a channel of a register; preserve .rgb
void RegWriteA(inout float4 Regs[16], uint idx, float a)
{
    uint i = idx & 0xFu;
    Regs[i].a = (i == PS_REGISTER_ZERO) ? Regs[i].a : a;
}

// ============================================================
// Input mapping
//
// PS_INPUTMAPPING values (bits [7:5] of the register byte):
//   0x00 = UNSIGNED_IDENTITY  max(0, v)
//   0x20 = UNSIGNED_INVERT    1 - sat(v)
//   0x40 = EXPAND_NORMAL      2 * max(0, v) - 1
//   0x60 = EXPAND_NEGATE     -2 * max(0, v) + 1
//   0x80 = HALFBIAS_NORMAL    max(0, v) - 0.5
//   0xA0 = HALFBIAS_NEGATE    0.5 - max(0, v)
//   0xC0 = SIGNED_IDENTITY    v
//   0xE0 = SIGNED_NEGATE     -v
//
// All 8 cases reduce to: decode a 2-bit mode from bits [7:6],
// compute a base expression per mode, then conditionally negate
// (bit [5]).  Implemented as branchless movc chains — no switch,
// no jump table, uniform execution across the wavefront.
// ============================================================

float4 ApplyInputMapping(uint mapping, float4 v)
{
    float4 clamped = max(0.0f, v);
    uint   mode    = (mapping >> 6u) & 3u;
    bool   negate  = (mapping & 0x20u) != 0u;

    float4 mUnsigned =        clamped;          // mode 0: max(0, v)
    float4 mExpand   = 2.0f * clamped - 1.0f;   // mode 1: 2*max(0,v)-1
    float4 mHalfbias =        clamped - 0.5f;   // mode 2: max(0,v)-0.5
    float4 mSigned   =        v;                // mode 3: v

    // Select base via movc chain
    float4 base = mUnsigned;
    base = (mode == 1u) ? mExpand   : base;
    base = (mode == 2u) ? mHalfbias : base;
    base = (mode == 3u) ? mSigned   : base;

    // Negate: mode 0 special-cases to 1-sat(v); all others are -base
    float4 neg = (mode == 0u) ? (1.0f - saturate(v)) : -base;

    return negate ? neg : base;
}

float ApplyInputMappingScalar(uint mapping, float v)
{
    float  clamped = max(0.0f, v);
    uint   mode    = (mapping >> 6u) & 3u;
    bool   negate  = (mapping & 0x20u) != 0u;

    float mUnsigned =        clamped;
    float mExpand   = 2.0f * clamped - 1.0f;
    float mHalfbias =        clamped - 0.5f;
    float mSigned   =        v;

    float base = mUnsigned;
    base = (mode == 1u) ? mExpand   : base;
    base = (mode == 2u) ? mHalfbias : base;
    base = (mode == 3u) ? mSigned   : base;

    float neg = (mode == 0u) ? (1.0f - saturate(v)) : -base;

    return negate ? neg : base;
}

// ============================================================
// Output mapping helpers
//
// flags = PS_COMBINEROUTPUT flags field (rgbOut >> 12 or aOut >> 12).
// Bits [5:3] encode the output mapping:
//   bit 3 = PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS: subtract 0.5 before scaling
//   bits [5:4] = scale: 00=x1  01=x2  10=x4  11=/2
//
// DecodeOutputBias/Scale extract per-stage-constant values once;
// the combiner stage applies them inline to all six results.
// ============================================================

float DecodeOutputBias(uint flags)
{
    return (flags & PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS) ? -0.5f : 0.0f;
}

float DecodeOutputScale(uint flags)
{
    // Scale values are powers of 2: m=0→1, m=1→2, m=2→4, m=3→0.5
    // Exponent sequence (0, 1, 2, -1) = ((m + 1) & 3) - 1
    uint m = (flags >> 4u) & 3u;
    return exp2((float)((m + 1u) & 3u) - 1.0f);
}

// ============================================================
// Input register resolution — split into stage vs final paths
//
// ResolveStageInput: called from DoCombinerStage (stageIdx 0..7)
//   Regs[V1R0_SUM] and Regs[EF_PROD] are zero-initialized and never
//   written until DoFinalCombiner, so all register indices (including
//   FOG, V1R0_SUM, EF_PROD) resolve correctly via a direct array read.
//
// ResolveFinalInput: called from DoFinalCombiner (EFG stageIdx=8,
//   ABCD stageIdx=9).  Remaps invalid mappings; V1R0_SUM/EF_PROD
//   are valid only at stageIdx=9.
// ============================================================

float4 ResolveStageInput(float4 Regs[16], uint regByte)
{
    uint regIdx    = regByte & 0x0Fu;
    uint mapping   = regByte & 0xE0u;
    bool useAlphaC = (regByte & PS_CHANNEL_ALPHA) != 0u;

    float4 val = Regs[regIdx];

    if (useAlphaC)
        val = val.aaaa;

    return ApplyInputMapping(mapping, val);
}

float ResolveStageInputAlpha(float4 Regs[16], uint regByte)
{
    uint regIdx  = regByte & 0x0Fu;
    uint mapping = regByte & 0xE0u;

    return ApplyInputMappingScalar(mapping, Regs[regIdx].a);
}

float4 ResolveFinalInput(float4 Regs[16], uint regByte, bool isFinalAB,
                         bool flagUniqueC0, bool flagUniqueC1)
{
    uint regIdx    = regByte & 0x0Fu;
    uint mapping   = regByte & 0xE0u;
    bool useAlphaC = (regByte & PS_CHANNEL_ALPHA) != 0u;

    float4 val;
    switch (regIdx)
    {
        // C0/C1 read from Regs[] (pre-loaded from PSFinalCombinerConstant
        // before DoFinalCombiner, matching NV2A/xemu behavior).

        case PS_REGISTER_FOG:
        {
            // Final combiner sees only FOG.a; rgb reads as zero
            float4 fog = Regs[PS_REGISTER_FOG];
            val = float4(0.0f, 0.0f, 0.0f, fog.a);
            break;
        }

        case PS_REGISTER_V1R0_SUM:
        case PS_REGISTER_EF_PROD:
            val = isFinalAB ? Regs[regIdx & 0xFu] : (float4)0.0f;
            break;

        default:
            val = Regs[regIdx & 0xFu];
            break;
    }

    if (useAlphaC)
        val = val.aaaa;

    // Remap final-combiner-invalid mappings to nearest valid equivalents
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

    return ApplyInputMapping(mapping, val);
}

// ============================================================
// PSInputTexture source-stage decoder
//
// PSInputTexture is a packed uint:
//   Stage 0: no input
//   Stage 1: always 0
//   Stage 2: bit 16 (1 bit: 0 or 1)
//   Stage 3: bits 20-21 (2 bits: 0, 1, or 2)
// ============================================================

uint GetSourceStage(uint stage)
{
    if (stage <= 1u) return 0u; // stage 0 has no pred; stage 1 always reads 0
    if (stage == 2u) return (PSInputTexture >> 16u) & 0x1u;
    /* stage 3 */   return (PSInputTexture >> 20u) & 0x3u;
}

// ============================================================
// PSCompareMode decoder
//
// 4 bits per stage (RSTQ), each bit selects LT (1) or GE (0).
// ============================================================

void ApplyCompareMode(uint stage, float4 coords)
{
    uint bits = (PSCompareMode >> (stage * 4u)) & 0xFu;
    // Each bit: 0 = GE (clip if >= 0), 1 = LT (clip if < 0)
    // Per NV2A: bit set means "discard if coord < 0"
    bool killR = (bits & 1u) != 0u ? (coords.x < 0.0f) : (coords.x >= 0.0f);
    bool killS = (bits & 2u) != 0u ? (coords.y < 0.0f) : (coords.y >= 0.0f);
    bool killT = (bits & 4u) != 0u ? (coords.z < 0.0f) : (coords.z >= 0.0f);
    bool killQ = (bits & 8u) != 0u ? (coords.w < 0.0f) : (coords.w >= 0.0f);
    if (killR || killS || killT || killQ)
        discard;
}

// ============================================================
// PSDotMapping decoder
//
// 3 bits per stage (stage 1 in bits 0-2, stage 2 in bits 4-6, stage 3 in bits 8-10).
// Returns the remapped float3 value of the source texture register for dot product.
// ============================================================

float3 ApplyDotMapping(uint stage, float4 src)
{
    uint mapping = (PSDotMapping >> ((stage - 1u) * 4u)) & 0x7u;

    // PS_DOTMAPPING values (from NV2A docs):
    // 0 = ZERO_TO_ONE        — identity
    // 1 = MINUS1_TO_1_D3D    — (v - 128) / 127
    // 2 = MINUS1_TO_1_GL     — two's complement
    // 3 = MINUS1_TO_1        — (v < 128 ? v : v-256) / 127
    // 4 = HILO_1             — 16-bit unsigned
    // 5 = HILO_HEMISPHERE_D3D
    // 6 = HILO_HEMISPHERE_GL
    // 7 = HILO_HEMISPHERE
    //
    // Most games use mapping 0 or 1. Implement the common cases;
    // HILO modes are extremely rare and would need 16-bit decode.
    if (mapping == 0u)
        return src.rgb; // ZERO_TO_ONE: identity [0,1]

    if (mapping == 1u) {
        // MINUS1_TO_1_D3D: (byte - 128) / 127
        float3 b = round(saturate(src.rgb) * 255.0f);
        return (b - 128.0f) / 127.0f;
    }
    // Fallback: treat as identity for unimplemented HILO/GL modes
    return src.rgb;
}

// ============================================================
// Post-process a sampled texel (matches compiled PS pipeline)
// ============================================================

float4 PostProcessTexel(uint stage, float4 t)
{
    // 1. Channel swizzle / luminance fixup
    t = ApplyTexFmtFixup(t, (int)TexFmtFixup[stage]);

    // 2. Color sign conversion
    [branch] if (any(ColorSign[stage] != 0.0f))
        t = PerformColorSign(ColorSign[stage], t);

    // 3. Color key
    [branch] if (ColorKeyOp[stage].x != 0.0f)
        t = PerformColorKeyOp((int)ColorKeyOp[stage].x, ColorKeyColor[stage], t);

    // 4. Alpha kill
    PerformAlphaKill((int)AlphaKill[stage], t);

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

void FetchTexture(inout float4 Regs[16], uint stage, uint mode)
{
    float4 val = float4(0.0f, 0.0f, 0.0f, 1.0f); // default: opaque black
    uint   tBase = PS_REGISTER_T0 + stage;

    // Read texture coordinates from T register (pre-loaded from VS output;
    // may have been perturbed by a prior BUMPENVMAP stage).
    float4 coords = Regs[tBase];

    switch (mode)
    {
    case PS_TEXTUREMODES_NONE:
        // Leave T[stage] at its initialised value (0,0,0,1); no write
        return;

    case PS_TEXTUREMODES_PROJECT2D:
        // No explicit projective divide — the VS output already accounts for it.
        // (The compiled PS also uses tex2D without /w.)
        val = Sample2D(stage, coords.xy);
        break;

    case PS_TEXTUREMODES_PROJECT3D:
        val = Sample3D(stage, coords.xyz);
        break;

    case PS_TEXTUREMODES_CUBEMAP:
        val = SampleCube(stage, coords.xyz);
        break;

    case PS_TEXTUREMODES_PASSTHRU:
        val = saturate(coords);
        break;

    case PS_TEXTUREMODES_CLIPPLANE:
        // Per-component compare mode from PSCompareMode (4 bits per stage)
        ApplyCompareMode(stage, coords);
        // val stays (0,0,0,1)
        break;

    case PS_TEXTUREMODES_BUMPENVMAP:
    case PS_TEXTUREMODES_BUMPENVMAP_LUM:
    {
        // Bump source = source stage's already-computed texel (matches compiled PS
        // BumpEnv() macro which reads src(ts), not the current stage's sample).
        float4 bumpSrc = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float4 bem = BEM[stage];
        // Perturb THIS stage's own texcoords, then sample with the perturbed coords
        float u = coords.x + bem.x * bumpSrc.r + bem.z * bumpSrc.g;
        float v = coords.y + bem.y * bumpSrc.r + bem.w * bumpSrc.g;
        val = Sample2D(stage, float2(u, v));
        // Apply luminance scaling for BUMPENVMAP_LUM
        if (mode == PS_TEXTUREMODES_BUMPENVMAP_LUM) {
            // LUM[stage].x = BumpEnvLScale, .y = BumpEnvLOffset
            float lumFactor = LUM[stage].x * bumpSrc.b + LUM[stage].y;
            val.rgb *= lumFactor;
        }
        break;
    }

    case PS_TEXTUREMODES_DPNDNT_AR:
    {
        // Dependent lookup: use .a and .r from source stage as (u,v)
        float4 src = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        val = Sample2D(stage, src.ar);
        break;
    }

    case PS_TEXTUREMODES_DPNDNT_GB:
    {
        // Dependent lookup: use .g and .b from source stage as (u,v)
        float4 src = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        val = Sample2D(stage, src.gb);
        break;
    }

    case PS_TEXTUREMODES_DOTPRODUCT:
    {
        // Apply dot mapping to the source stage's texture register, then dot with coords
        float4 src = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm  = ApplyDotMapping(stage, src);
        float  d   = dot(coords.xyz, dm);
        val = float4(d, 0.0f, 0.0f, 1.0f);
        break;
    }

    case PS_TEXTUREMODES_DOT_ST:
    {
        // texm3x2tex: Normal2(ts) = (dot_[ts-1], dot_[ts])
        // Compute current dot, then use (previous dot, current dot) as (s,t)
        float4 src = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm  = ApplyDotMapping(stage, src);
        float  d   = dot(coords.xyz, dm);
        float  s   = Regs[tBase - 1u].x; // dot_[ts-1]
        float  t   = d;                   // dot_[ts]
        val = Sample2D(stage, float2(s, t));
        break;
    }

    case PS_TEXTUREMODES_DOT_ZW:
    {
        // texm3x2depth: compute n = (prev_dot, current_dot), depth = n.x / n.y
        float4 src = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm  = ApplyDotMapping(stage, src);
        float  d   = dot(coords.xyz, dm);
        float  prevDot = Regs[tBase - 1u].x;
        // Avoid division by near-zero (matches compiled PS guard)
        float  depth = (abs(d) < 0.00001f) ? 1.0f : (prevDot / d);
        val = float4(depth, depth, depth, depth);
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_DIFF:
    {
        // texm3x3diff: restricted to stage 2. Normal2(ts) = (dot_[ts-1], dot_[ts], 0).
        // Normal is used directly as the cubemap direction (no reflection for DIFF).
        float4 src  = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm   = ApplyDotMapping(stage, src);
        float  nx   = Regs[tBase - 1u].x; // dot_[ts-1]
        float  ny   = dot(coords.xyz, dm);  // dot_[ts]
        val = SampleCube(stage, float3(nx, ny, 0.0f));
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_SPEC:
    {
        // Build normal, compute specular reflection vector, cubemap lookup.
        // Eye vector is assembled from the .w component of T1, T2, T3
        // (NV2A hardcodes these indices).
        float  nx   = Regs[tBase - 2u].x;
        float  ny   = Regs[tBase - 1u].x;
        float4 src  = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm   = ApplyDotMapping(stage, src);
        float  nz   = dot(coords.xyz, dm);
        float3 N    = normalize(float3(nx, ny, nz));
        float3 E    = normalize(float3(Regs[PS_REGISTER_T1].w,
                                       Regs[PS_REGISTER_T2].w,
                                       Regs[PS_REGISTER_T3].w));
        float3 R    = 2.0f * dot(N, E) * N - E;
        val = SampleCube(stage, R);
        break;
    }

    case PS_TEXTUREMODES_DOT_STR_3D:
    {
        float  s    = Regs[tBase - 2u].x;
        float  t    = Regs[tBase - 1u].x;
        float4 src  = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm   = ApplyDotMapping(stage, src);
        float  r    = dot(coords.xyz, dm);
        val = Sample3D(stage, float3(s, t, r));
        break;
    }

    case PS_TEXTUREMODES_DOT_STR_CUBE:
    {
        float  s    = Regs[tBase - 2u].x;
        float  t    = Regs[tBase - 1u].x;
        float4 src  = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm   = ApplyDotMapping(stage, src);
        float  r    = dot(coords.xyz, dm);
        val = SampleCube(stage, float3(s, t, r));
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST:
    {
        // Like DOT_RFLCT_SPEC but eye vector is a constant.
        // TODO: wire SetEyeVector() / D3DRS_PSINPUTTEXTURE to a cbuffer entry
        float  nx   = Regs[tBase - 2u].x;
        float  ny   = Regs[tBase - 1u].x;
        float4 src  = Regs[PS_REGISTER_T0 + GetSourceStage(stage)];
        float3 dm   = ApplyDotMapping(stage, src);
        float  nz   = dot(coords.xyz, dm);
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

    // Post-process: format fixup, color sign, color key (matches compiled PS pipeline)
    val = PostProcessTexel(stage, val);
    RegWriteDirect(Regs, tBase, val);
}

// NV2A-accurate multiply helpers are in CxbxNV2AMathHelpers.hlsli

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
    float4 rgbA = ResolveStageInput(Regs, (rgbIn >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu);
    float4 rgbB = ResolveStageInput(Regs, (rgbIn >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu);
    float4 rgbC = ResolveStageInput(Regs, (rgbIn >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu);
    float4 rgbD = ResolveStageInput(Regs, (rgbIn >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu);
    float   aA  = ResolveStageInputAlpha(Regs, (aIn   >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu);
    float   aB  = ResolveStageInputAlpha(Regs, (aIn   >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu);
    float   aC  = ResolveStageInputAlpha(Regs, (aIn   >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu);
    float   aD  = ResolveStageInputAlpha(Regs, (aIn   >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu);

    // --- Compute AB and CD products ---
    // Dot product applies to RGB only; alpha always multiplies scalars
    // NV2A enforces 0*anything=0 (even 0*inf); use nv2a_mul to match hardware
    float3 rgbAB = flagABDot ? (float3)dot(rgbA.rgb, rgbB.rgb) : nv2a_mul3(rgbA.rgb, rgbB.rgb);
    float3 rgbCD = flagCDDot ? (float3)dot(rgbC.rgb, rgbD.rgb) : nv2a_mul3(rgbC.rgb, rgbD.rgb);
    float   aAB  = nv2a_mul1(aA, aB);
    float   aCD  = nv2a_mul1(aC, aD);

    // --- SUM or MUX ---
    // Only decode R0.a when a MUX flag is actually set
    bool muxSel = false;
    if (flagRGBMux || flagAMux) {
        float r0a = Regs[PS_REGISTER_R0].a;
        uint r0aBits = (uint)round(saturate(r0a) * 255.0f);
        muxSel = flagMuxMsb
            ? (r0a >= 0.5f)
            : ((r0aBits & 1u) != 0u);
    }

    // Skip sum computation when the result won't be written
    bool writeSumRGB = !flagABDot && !flagCDDot;
    float3 rgbABCD = flagRGBMux ? (muxSel ? rgbCD : rgbAB)
                   : writeSumRGB ? (rgbAB + rgbCD) : (float3)0.0f;
    float   aABCD  = flagAMux   ? (muxSel ? aCD   : aAB  ) : (aAB   + aCD  );

    // --- Hoisted output mapping: decode bias/scale once per stage ---
    float rgbBias  = DecodeOutputBias(rgbFlags);
    float rgbScale = DecodeOutputScale(rgbFlags);
    float aBias    = DecodeOutputBias(aFlags);
    float aScale   = DecodeOutputScale(aFlags);

    float3 outRGB_AB  = clamp((rgbAB   + rgbBias) * rgbScale, -1.0f, 1.0f);
    float3 outRGB_CD  = clamp((rgbCD   + rgbBias) * rgbScale, -1.0f, 1.0f);
    float3 outRGB_Sum = clamp((rgbABCD + rgbBias) * rgbScale, -1.0f, 1.0f);
    float  outA_AB    = clamp((aAB     + aBias)   * aScale,   -1.0f, 1.0f);
    float  outA_CD    = clamp((aCD     + aBias)   * aScale,   -1.0f, 1.0f);
    float  outA_Sum   = clamp((aABCD   + aBias)   * aScale,   -1.0f, 1.0f);

    // --- RGB writes ---
    // AB: write RGB; optionally propagate .b to .a (BlueToAlpha)
    // Skip BlueToAlpha read when alpha write will overwrite immediately
    {
        bool abAlphaOverwrite = (aRegAB == rgbRegAB) && (aRegAB != PS_REGISTER_DISCARD);
        float a = (abBlue2A && !abAlphaOverwrite) ? outRGB_AB.b : Regs[rgbRegAB & 0xFu].a;
        RegWrite(Regs, rgbRegAB, float4(outRGB_AB, a));
    }
    // CD: same pattern
    {
        bool cdAlphaOverwrite = (aRegCD == rgbRegCD) && (aRegCD != PS_REGISTER_DISCARD);
        float a = (cdBlue2A && !cdAlphaOverwrite) ? outRGB_CD.b : Regs[rgbRegCD & 0xFu].a;
        RegWrite(Regs, rgbRegCD, float4(outRGB_CD, a));
    }
    // AB+CD sum/mux: write RGB only; spec requires DISCARD when any DOT flag is set
    if (writeSumRGB)
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
    [branch] if (PSFinalCombinerInputsABCD == 0u && PSFinalCombinerInputsEFG == 0u)
        return R0;

    // PSFinalCombinerInputsEFG = PS_COMBINERINPUTS(E, F, G, settings)
    //   = (E<<24) | (F<<16) | (G<<8) | settings
    uint efg      = PSFinalCombinerInputsEFG;
    // PS_COMBINERINPUTS(E, F, G, settings) for the final combiner
    uint settings = (efg >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu;
    uint eReg     = (efg >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu;
    uint fReg     = (efg >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu;
    uint gReg     = (efg >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu;

    // --- Resolve E, F (RGB) and G (alpha) — EFG phase (not ABCD) ---
    float3 E = ResolveFinalInput(Regs, eReg, false, flagUniqueC0, flagUniqueC1).rgb;
    float3 F = ResolveFinalInput(Regs, fReg, false, flagUniqueC0, flagUniqueC1).rgb;
    float  G = ResolveFinalInput(Regs, gReg, false, flagUniqueC0, flagUniqueC1).a;

    // Compute E*F and store in EF_PROD for potential use by ABCD inputs
    RegWriteDirect(Regs, PS_REGISTER_EF_PROD, float4(E * F, 1.0f));

    // --- Optional complement and clamp on V1 and R0 ---
    // These modify the sum inputs, not the stored register values
    float3 v1 = RegRead(Regs, PS_REGISTER_V1).rgb;
    float3 r0 = R0.rgb;
    if (settings & PS_FINALCOMBINERSETTING_COMPLEMENT_V1) v1 = 1.0f - v1;
    if (settings & PS_FINALCOMBINERSETTING_COMPLEMENT_R0) r0 = 1.0f - r0;

    float3 v1r0sum = v1 + r0;
    if (settings & PS_FINALCOMBINERSETTING_CLAMP_SUM) v1r0sum = saturate(v1r0sum);

    // Store V1+R0 sum for potential use by ABCD inputs
    RegWriteDirect(Regs, PS_REGISTER_V1R0_SUM, float4(v1r0sum, 1.0f));

    // --- Resolve A, B, C, D — ABCD phase (V1R0_SUM / EF_PROD now valid) ---
    uint abcd = PSFinalCombinerInputsABCD;
    float4 A = ResolveFinalInput(Regs, (abcd >> PS_COMBINERINPUTS_A_SHIFT) & 0xFFu, true, flagUniqueC0, flagUniqueC1);
    float4 B = ResolveFinalInput(Regs, (abcd >> PS_COMBINERINPUTS_B_SHIFT) & 0xFFu, true, flagUniqueC0, flagUniqueC1);
    float4 C = ResolveFinalInput(Regs, (abcd >> PS_COMBINERINPUTS_C_SHIFT) & 0xFFu, true, flagUniqueC0, flagUniqueC1);
    float4 D = ResolveFinalInput(Regs, (abcd >> PS_COMBINERINPUTS_D_SHIFT) & 0xFFu, true, flagUniqueC0, flagUniqueC1);

    // Final RGB = A*B + (1-A)*C + D, clamped to [0,1]
    // Final alpha = G
    float4 result;
    result.rgb = saturate(lerp(C.rgb, B.rgb, A.rgb) + D.rgb);
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
    uint4 texMode = uint4(
        (PSTextureModes      ) & 0x1Fu,
        (PSTextureModes >>  5) & 0x1Fu,
        (PSTextureModes >> 10) & 0x1Fu,
        (PSTextureModes >> 15) & 0x1Fu
    );

    // --- Initialise the register file ---
    // Zero-init sets ZERO (index 0) permanently to 0.
    float4 Regs[16];
    [unroll] for (uint i = 0u; i < 16u; i++) Regs[i] = 0.0f;

    // --- Texture stages (must run sequentially; later stages can read earlier T regs) ---
    // Pre-load T registers with VS output tex coords. BUMPENVMAP stages
    // may perturb a later T register before that stage samples.
    Regs[PS_REGISTER_T0] = input.iT0;
#if NUM_TEXTURE_STAGES >= 2
    Regs[PS_REGISTER_T1] = input.iT1;
#endif
#if NUM_TEXTURE_STAGES >= 3
    Regs[PS_REGISTER_T2] = input.iT2;
#endif
#if NUM_TEXTURE_STAGES >= 4
    Regs[PS_REGISTER_T3] = input.iT3;
#endif

    FetchTexture(Regs, 0u, texMode.x);
#if NUM_TEXTURE_STAGES >= 2
    FetchTexture(Regs, 1u, texMode.y);
#endif
#if NUM_TEXTURE_STAGES >= 3
    FetchTexture(Regs, 2u, texMode.z);
#endif
#if NUM_TEXTURE_STAGES >= 4
    FetchTexture(Regs, 3u, texMode.w);
#endif

    // --- Set vertex-derived registers ---
    // Use FRONTFACE_FACTOR to match compiled PS winding-order correction:
    // 0 = always front, +/-1 = two-sided with CW/CCW convention
    bool isFront = (input.iFF ? 1.0f : -1.0f) * FrontFaceInfo.x >= 0.0f;
    float4 diffuse  = isFront ? input.iD0 : input.iB0;
    float4 specular = isFront ? input.iD1 : input.iB1;
    RegWriteDirect(Regs, PS_REGISTER_V0, diffuse);
    RegWriteDirect(Regs, PS_REGISTER_V1, specular);
    // FOG: rgb from the fog color constant, alpha from the vertex fog factor
    RegWriteDirect(Regs, PS_REGISTER_FOG, float4(FogColor.rgb, saturate(input.iFog)));

    // R0 initialization: NV2A spec says R0.rgb starts at 0, R0.a starts from T0.a.
    // xemu matches this: "r0 = vec4(0); r0.a = t0.a;".
    // Previous code copied all of T0 into R0 which was incorrect.
    RegWriteA(Regs, PS_REGISTER_R0, RegRead(Regs, PS_REGISTER_T0).a);
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
        {
            // Initialise C0/C1 from cbuffer constants before each stage.
            // UNIQUE_C0/C1: each stage gets its own constant; otherwise all
            // stages share constant[0].  Combiner outputs CAN write to C0/C1
            // (confirmed by xemu / NV2A hardware), so this must happen before
            // DoCombinerStage, not inside ResolveStageInput.
            Regs[PS_REGISTER_C0] = flagUniqueC0 ? PSConstant0[stage] : PSConstant0[0];
            Regs[PS_REGISTER_C1] = flagUniqueC1 ? PSConstant1[stage] : PSConstant1[0];
            DoCombinerStage(Regs, stage, flagMuxMsb, flagUniqueC0, flagUniqueC1);
        }
    }

    // Initialise C0/C1 for the final combiner from FinalCombinerConstants
    Regs[PS_REGISTER_C0] = PSFinalCombinerConstant[0];
    Regs[PS_REGISTER_C1] = PSFinalCombinerConstant[1];
    float4 result = DoFinalCombiner(Regs, flagUniqueC0, flagUniqueC1);

    // --- Alpha test ---
    PerformAlphaTest(AlphaTest.xyz, result.a);

    // --- Fog blending ---
    // FogInfo: x=tableMode, y=density, z=start, w=end
    [branch] if (FogEnable != 0u) {
        float fogFactor = CalculateFogFactor((int)FogInfo.x, FogInfo.y,
                                             FogInfo.z, FogInfo.w, input.iFog);
        result.rgb = lerp(FogColor.rgb, result.rgb, fogFactor);
    }

    return result;
}
