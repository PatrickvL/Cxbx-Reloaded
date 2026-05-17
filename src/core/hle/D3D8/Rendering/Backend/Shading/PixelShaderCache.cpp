// PixelShaderCache.cpp — Runtime NV2A register combiner → HLSL JIT compiler + cache
//
// Reads the current PGRAPH register combiner topology (input routing,
// output destinations, texture modes, final combiner config) and generates
// a self-contained HLSL pixel shader with straight-line math — no loops,
// no register file array, no runtime input decode.
//
// Dynamic per-frame data (C0/C1 colors, fog color, bump matrices) is still
// read from g_PGRegs StructuredBuffer at runtime; only the STRUCTURE is baked.
//
// Compile with D3DCompile (ps_5_0, O3) and cache by hash of topology.

#define LOG_PREFIX CXBXR_MODULE::PXSH

#include "PixelShaderCache.h"
#include "CxbxJITIncludeHandler.h"
#include "ShaderDiskCache.h"
#include "core/kernel/init/CxbxKrnl.h"
#include "common/Logging.h"
#include "core/hle/D3D8/Rendering/Shaders/CxbxNV2APixelShaderConstants.hlsli"
#include "core/hle/D3D8/Rendering/Shaders/CxbxRegisterCombinerInterpreterState.hlsli"
#include "devices/Xbox.h"
#include "devices/video/nv2a.h"
#include "../Backend_D3D11_Profiler.h"

#include <unordered_map>
#include <string>
#include <sstream>
#include <mutex>
#include <cstring>
#include <d3dcompiler.h>
#include "common/util/hasher.h"

// Externally stored last aux CB (built by CxbxD3D11UploadRCInterpreterState)
extern PSAuxCBLayout g_LastPSAuxCB;

// ============================================================
// State capture — everything that affects generated HLSL
// ============================================================

// Bump this whenever the JIT HLSL infrastructure changes (e.g. resource
// binding type changes, new includes, register slot moves). This ensures
// stale disk-cached .cso files are not reused after incompatible changes.
static constexpr uint32_t PS_JIT_CACHE_VERSION = 3; // v3: PSAuxCB fields derived from PGRAPH SRV

struct PSJITKey {
    uint32_t cacheVersion;      // invalidates disk cache on HLSL infra changes
    uint32_t numStages;
    uint32_t combinectl;        // full COMBINECTL (includes flags)
    uint32_t rgbInputs[8];
    uint32_t alphaInputs[8];
    uint32_t rgbOutputs[8];
    uint32_t alphaOutputs[8];
    uint32_t fcABCD;            // adjusted final combiner
    uint32_t fcEFG;
    uint32_t textureModes;      // adjusted from aux CB
    uint32_t shaderCtl;
    uint32_t shaderClipMode;
    uint32_t fogEnable;
    float    frontFaceFactor;
    float    alphaKill[4];
    float    colorKeyOp[4];
    float    shadowCompare[4];
    float    colorSign[16];     // 4 stages × float4 (r,g,b,a)
    float    texFmtFixup[4];
};

static uint64_t HashKey(const PSJITKey& key)
{
    return ComputeHash(&key, sizeof(key));
}

// ============================================================
// Cache
// ============================================================
struct CachedPSShader {
    ID3D11PixelShader* pPS;
};

static std::unordered_map<uint64_t, CachedPSShader> g_PSJITCache;
static std::mutex g_PSJITMutex;

PixelShaderCache g_PixelShaderCache;

void PixelShaderCache::Clear()
{
    std::lock_guard<std::mutex> lock(g_PSJITMutex);
    for (auto& pair : g_PSJITCache) {
        if (pair.second.pPS) pair.second.pPS->Release();
    }
    g_PSJITCache.clear();
}

// ============================================================
// Dead-code analysis — which registers / channels are actually read
// ============================================================

// Channel flags (2-bit per register slot)
static constexpr uint32_t CH_NONE = 0;
static constexpr uint32_t CH_A    = 1;  // alpha only
static constexpr uint32_t CH_RGB  = 2;  // rgb only
static constexpr uint32_t CH_RGBA = 3;  // both

struct ReadMasks {
    uint32_t regChans;   // 16 slots × 2 bits = 32 bits (regIdx 0..15)
    uint16_t stageC0;    // one bit per combiner stage (0..7), bit 8 = final combiner
    uint16_t stageC1;    // one bit per combiner stage (0..7), bit 8 = final combiner
    uint16_t texNeeded;  // one bit per texture stage (0..3) — needed by combiners or cross-stage deps
};

static inline void SetChan(uint32_t& regChans, uint32_t regIdx, uint32_t ch)
{
    uint32_t shift = (regIdx & PS_REGISTER_MASK) * 2;
    regChans |= (ch << shift);
}

static inline uint32_t GetChan(uint32_t regChans, uint32_t regIdx)
{
    uint32_t shift = (regIdx & PS_REGISTER_MASK) * 2;
    return (regChans >> shift) & 3;
}

// Extract one of four input bytes (A/B/C/D) packed into a 32-bit word.
// Slot 0 = bits 31:24, slot 1 = bits 23:16, slot 2 = bits 15:8, slot 3 = bits 7:0.
static inline uint32_t InputByte(uint32_t word, int slot)
{
    return (word >> ((3 - slot) * 8)) & 0xFF;
}

// Compute the source stage for a texture dependency at a given stage.
static inline uint32_t GetTexSrcStage(uint32_t stage, uint32_t shaderCtl)
{
    return (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
}

// Unpack the 4 texture modes from the packed PSTextureModes word.
static inline void UnpackTexModes(uint32_t packed, uint32_t texModes[4])
{
    texModes[0] = (packed      ) & PS_TEXTUREMODES_MASK;
    texModes[1] = (packed >>  5) & PS_TEXTUREMODES_MASK;
    texModes[2] = (packed >> 10) & PS_TEXTUREMODES_MASK;
    texModes[3] = (packed >> 15) & PS_TEXTUREMODES_MASK;
}

// ============================================================
// HLSL code generation helpers
// ============================================================

static const char* RegName(uint32_t regIdx)
{
    switch (regIdx) {
    case PS_REGISTER_ZERO:     return nullptr;  // ZERO — handled specially
    case PS_REGISTER_C0:       return "C0";
    case PS_REGISTER_C1:       return "C1";
    case PS_REGISTER_FOG:      return "FOG";
    case PS_REGISTER_V0:       return "V0";
    case PS_REGISTER_V1:       return "V1";
    case PS_REGISTER_T0:       return "T0";
    case PS_REGISTER_T1:       return "T1";
    case PS_REGISTER_T2:       return "T2";
    case PS_REGISTER_T3:       return "T3";
    case PS_REGISTER_R0:       return "R0";
    case PS_REGISTER_R1:       return "R1";
    case PS_REGISTER_V1R0_SUM: return "V1R0_SUM";
    case PS_REGISTER_EF_PROD:  return "EF_PROD";
    default:                   return "R0"; // reserved → R0 fallback
    }
}

// Emit a combiner input expression with channel select + input mapping applied.
// isAlpha=false → float4 (RGB path, .aaaa swizzle when bit 4 set)
// isAlpha=true  → scalar (alpha path, .a or .b selection)
static std::string EmitInput(uint32_t inputByte, bool isAlpha)
{
    uint32_t regIdx  = inputByte & PS_REGISTER_MASK;
    uint32_t mapping = inputByte & PS_INPUTMAPPING_MASK;
    bool useAlphaChan = (inputByte & 0x10) != 0;

    const char* maxPfx = isAlpha ? "max(0.0, " : "max((float4)0.0, ";

    // ZERO register constant folding
    if (regIdx == PS_REGISTER_ZERO) {
        switch (mapping) {
        case PS_INPUTMAPPING_UNSIGNED_IDENTITY: return isAlpha ? "0.0"    : "(float4)0.0";     // 0x00
        case PS_INPUTMAPPING_UNSIGNED_INVERT:   return isAlpha ? "1.0"    : "(float4)1.0";     // 0x20
        case PS_INPUTMAPPING_EXPAND_NORMAL:     return isAlpha ? "(-1.0)" : "(float4)(-1.0)";  // 0x40
        case PS_INPUTMAPPING_EXPAND_NEGATE:     return isAlpha ? "1.0"    : "(float4)1.0";     // 0x60
        case PS_INPUTMAPPING_HALFBIAS_NORMAL:   return isAlpha ? "(-0.5)" : "(float4)(-0.5)";  // 0x80
        case PS_INPUTMAPPING_HALFBIAS_NEGATE:   return isAlpha ? "0.5"    : "(float4)0.5";     // 0xa0
        case PS_INPUTMAPPING_SIGNED_IDENTITY:   return isAlpha ? "0.0"    : "(float4)0.0";     // 0xc0
        case PS_INPUTMAPPING_SIGNED_NEGATE:     return isAlpha ? "0.0"    : "(float4)0.0";     // 0xe0
        default:                                return isAlpha ? "0.0"    : "(float4)0.0";
        }
    }

    // Value expression
    std::string val = RegName(regIdx);
    if (isAlpha)
        val += useAlphaChan ? ".a" : ".b";
    else if (useAlphaChan)
        val += ".aaaa";

    // Apply input mapping
    switch (mapping) {
    case PS_INPUTMAPPING_UNSIGNED_IDENTITY: return std::string(maxPfx) + val + ")";                       // 0x00
    case PS_INPUTMAPPING_UNSIGNED_INVERT:   return "(1.0 - saturate(" + val + "))";                       // 0x20
    case PS_INPUTMAPPING_EXPAND_NORMAL:     return "(2.0 * " + std::string(maxPfx) + val + ") - 1.0)";    // 0x40
    case PS_INPUTMAPPING_EXPAND_NEGATE:     return "(-(2.0 * " + std::string(maxPfx) + val + ") - 1.0))"; // 0x60
    case PS_INPUTMAPPING_HALFBIAS_NORMAL:   return "(" + std::string(maxPfx) + val + ") - 0.5)";          // 0x80
    case PS_INPUTMAPPING_HALFBIAS_NEGATE:   return "(-(" + std::string(maxPfx) + val + ") - 0.5))";       // 0xa0
    case PS_INPUTMAPPING_SIGNED_IDENTITY:   return val;                                                   // 0xc0
    case PS_INPUTMAPPING_SIGNED_NEGATE:     return "(-" + val + ")";                                      // 0xe0
    default:                                return val;
    }
}

// Emit a final combiner input (restricted mappings)
static std::string EmitFinalInput(uint32_t inputByte, bool isFinalABCD)
{
    uint32_t regIdx  = inputByte & PS_REGISTER_MASK;
    uint32_t mapping = inputByte & PS_INPUTMAPPING_MASK;
    bool useAlpha    = (inputByte & 0x10) != 0;

    // Restrict mapping for final combiner: anything >= EXPAND_NORMAL → just bit 5
    if (mapping >= PS_INPUTMAPPING_EXPAND_NORMAL)
        mapping = mapping & PS_INPUTMAPPING_UNSIGNED_INVERT;

    std::string val;
    if (regIdx == PS_REGISTER_ZERO) {
        if (mapping == PS_INPUTMAPPING_UNSIGNED_INVERT) return "(float4)1.0"; // 0x20
        return "(float4)0.0";
    }

    val = RegName(regIdx);

    // FOG special handling: rgb→0, alpha passthrough
    if (regIdx == PS_REGISTER_FOG) {
        val = "float4(0.0, 0.0, 0.0, FOG.a)";
    }

    // V1R0_SUM / EF_PROD only valid in ABCD phase
    if ((regIdx == PS_REGISTER_V1R0_SUM || regIdx == PS_REGISTER_EF_PROD) && !isFinalABCD) {
        val = "(float4)0.0";
    }

    if (useAlpha)
        val += ".aaaa";

    // Apply the restricted mapping (only UNSIGNED_IDENTITY or UNSIGNED_INVERT)
    if (mapping == PS_INPUTMAPPING_UNSIGNED_INVERT) // 0x20
        return "(1.0 - saturate(" + val + "))";
    else
        return "max((float4)0.0, " + val + ")";
}

// ============================================================
// Texture mode code generation
// ============================================================

static void EmitTextureFetch(std::ostringstream& ss, uint32_t stage, uint32_t mode,
                             uint32_t shaderCtl, float shadowCompare,
                             float colorSign4[4], float texFmtFixup,
                             float colorKeyOp, float alphaKill)
{
    std::string sIdx = std::to_string(stage);
    std::string tReg = "T" + sIdx;
    std::string coords = "input.iT" + sIdx;

    if (mode == PS_TEXTUREMODES_NONE) {
        // NONE — T register keeps VS texcoord (already initialized)
        return;
    }

    // Emit the sample
    switch (mode) {
    case PS_TEXTUREMODES_PROJECT2D:
        ss << "    { float3 proj = " << coords << ".xyz / " << coords << ".w;\n";
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", proj.xy);\n";
        if (shadowCompare != 0.0f) {
            ss << "      " << tReg << " = ApplyShadowCompare(PG_UINT(0x19A4) & 7u, " << tReg << ", proj.z);\n";
        }
        ss << "    }\n";
        break;

    case PS_TEXTUREMODES_PROJECT3D:
        ss << "    { float3 proj = " << coords << ".xyz / " << coords << ".w;\n";
        if (shadowCompare != 0.0f) {
            ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", proj.xy);\n";
            ss << "      " << tReg << " = ApplyShadowCompare(PG_UINT(0x19A4) & 7u, " << tReg << ", proj.z);\n";
        } else {
            ss << "      " << tReg << " = Tex3D_" << sIdx << ".Sample(Samp" << sIdx << ", proj);\n";
        }
        ss << "    }\n";
        break;

    case PS_TEXTUREMODES_CUBEMAP:
        ss << "    " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx << ", " << coords << ".xyz);\n";
        break;

    case PS_TEXTUREMODES_PASSTHRU:
        ss << "    " << tReg << " = saturate(" << coords << ");\n";
        break;

    case PS_TEXTUREMODES_CLIPPLANE:
        ss << "    ApplyCompareMode((PG_UINT(0x1994) >> " << (stage * 4) << "u) & 0xFu, " << coords << ");\n";
        return; // no post-process

    case PS_TEXTUREMODES_BUMPENVMAP:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        uint32_t bumOff = (stage - 1) * 4;
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float4 bem = float4(PG_FLOAT(0x" << std::hex << (0x181C + bumOff) << "), PG_FLOAT(0x" << (0x1828 + bumOff) << "), PG_FLOAT(0x" << (0x1834 + bumOff) << "), PG_FLOAT(0x" << (0x1840 + bumOff) << "));\n" << std::dec;
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", float2("
           << coords << ".x + bem.x * src.r + bem.z * src.g, "
           << coords << ".y + bem.y * src.r + bem.w * src.g));\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_BUMPENVMAP_LUM:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        uint32_t bumOff = (stage - 1) * 4;
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float4 bem = float4(PG_FLOAT(0x" << std::hex << (0x181C + bumOff) << "), PG_FLOAT(0x" << (0x1828 + bumOff) << "), PG_FLOAT(0x" << (0x1834 + bumOff) << "), PG_FLOAT(0x" << (0x1840 + bumOff) << "));\n";
        ss << "      float lumS = PG_FLOAT(0x" << (0x1858 + (stage - 1) * 4) << ");\n";
        ss << "      float lumO = PG_FLOAT(0x" << (0x184C + (stage - 1) * 4) << ");\n" << std::dec;
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", float2("
           << coords << ".x + bem.x * src.r + bem.z * src.g, "
           << coords << ".y + bem.y * src.r + bem.w * src.g));\n";
        ss << "      " << tReg << ".rgb *= lumS * src.b + lumO;\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_BRDF:
        ss << "    " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", " << coords << ".xy);\n";
        break;

    case PS_TEXTUREMODES_DPNDNT_AR:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", T" << srcStage << ".ar);\n";
        break;
    }

    case PS_TEXTUREMODES_DPNDNT_GB:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", T" << srcStage << ".gb);\n";
        break;
    }

    case PS_TEXTUREMODES_DOTPRODUCT:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = float4(dot(" << coords << ".xyz, dm), 0.0, 0.0, 1.0);\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_ST:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx
           << ", float2(T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_ZW:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      float d = dot(" << coords << ".xyz, dm);\n";
        // WARNING: Do NOT add an epsilon guard here (e.g. "abs(d) < 1e-5 ? 1.0 : ...").
        // For D24S8, the VS passes c2 = (0, 0, 1/16777215) so d ≈ 5.96e-8 which is
        // a perfectly valid denominator. An epsilon guard would force depth to 1.0,
        // making the z-sprite always render at near-plane and breaking occlusion.
        // HLSL handles d=0 gracefully: INF is clamped by saturate() in the output.
        ss << "      float depth = T" << (stage - 1) << ".x / d;\n";
        ss << "      " << tReg << " = depth.xxxx;\n";
        ss << "      fragDepth = depth;\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_DIFF:
    {
        // xemu reference: normal = (dot[stage-1], dot[stage], dot_next)
        // where dot_next peeks at the next stage's dot mapping and source.
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        uint32_t dotMapping = (shaderCtl >> ((stage - 1) * 4)) & 7;
        uint32_t nextStage = stage + 1;
        uint32_t nextSrcStage = GetTexSrcStage(nextStage, shaderCtl);
        uint32_t nextDotMapping = (shaderCtl >> ((nextStage - 1) * 4)) & 7;
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << dotMapping << ", src);\n";
        ss << "      float curDot = dot(" << coords << ".xyz, dm);\n";
        ss << "      float4 nextSrc = T" << nextSrcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << nextSrcStage << "] > 0.5)\n";
        ss << "        nextSrc = (DepthTexAlias[" << nextSrcStage << "] >= 1.5) ? RemapD16ToColor(nextSrc.r)\n";
        ss << "            : RemapD24S8ToColor(nextSrc.r, TexStencil_" << nextSrcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 nextDm = ApplyDotMapping(" << nextDotMapping << ", nextSrc);\n";
        ss << "      float nextDot = dot(T" << nextStage << ".xyz, nextDm);\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(T" << (stage - 1) << ".x, curDot, nextDot));\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_SPEC:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      float3 N = normalize(float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "      float3 E = normalize(eyeVec);\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx << ", 2.0 * dot(N, E) * N - E);\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_STR_3D:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = Tex3D_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_STR_CUBE:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "    }\n";
        break;
    }

    case PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST:
    {
        uint32_t srcStage = GetTexSrcStage(stage, shaderCtl);
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      [flatten] if (DepthTexAlias[" << srcStage << "] > 0.5)\n";
        ss << "        src = (DepthTexAlias[" << srcStage << "] >= 1.5) ? RemapD16ToColor(src.r)\n";
        ss << "            : RemapD24S8ToColor(src.r, TexStencil_" << srcStage << ".Load(int3((int2)input.iPos.xy, 0)).g);\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      float3 N = normalize(float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(2.0*N.z*N.x, 2.0*N.z*N.y, 2.0*N.z*N.z - 1.0));\n";
        ss << "    }\n";
        break;
    }

    default:
        // Unknown mode — leave T register as texcoord passthrough
        return;
    }

    // Post-process the sampled texel (for modes that produce a texel)
    // PASSTHRU (0x04) runs post-process for hardware parity — per-stage
    // fixup/colorSign/alphaKill values should be identity when no texture
    // is bound, but we honor whatever the game sets.
    if (mode == PS_TEXTUREMODES_CLIPPLANE) return; // CLIPPLANE has no texel

    // Texture format fixup
    if (texFmtFixup != 0.0f) {
        int fixup = (int)texFmtFixup;
        switch (fixup) {
        case 1: ss << "    " << tReg << " = " << tReg << ".gbar;\n"; break;
        case 2: ss << "    " << tReg << " = " << tReg << ".abgr;\n"; break;
        case 3: ss << "    " << tReg << " = float4(" << tReg << ".r, " << tReg << ".r, " << tReg << ".r, " << tReg << ".a);\n"; break;
        case 4: ss << "    " << tReg << " = float4(" << tReg << ".r, " << tReg << ".r, " << tReg << ".r, " << tReg << ".g);\n"; break;
        case 5: ss << "    " << tReg << " = float4(" << tReg << ".rgb, 1.0);\n"; break;
        }
    }

    // Color sign conversion
    bool hasColorSign = (colorSign4[0] != 0.0f || colorSign4[1] != 0.0f ||
                         colorSign4[2] != 0.0f || colorSign4[3] != 0.0f);
    if (hasColorSign) {
        ss << "    " << tReg << " = PerformColorSign(ColorSign[" << stage << "], " << tReg << ");\n";
    }

    // Color key
    if (colorKeyOp != 0.0f) {
        const char* swizzle[] = { "x", "y", "z", "w" };
        ss << "    " << tReg << " = PerformColorKeyOp((int)DeriveColorKeyOp()." << swizzle[stage] << ", DeriveColorKeyColor(" << stage << "), " << tReg << ");\n";
    }

    // Alpha kill
    if (alphaKill != 0.0f) {
        const char* swizzle[] = { "x", "y", "z", "w" };
        ss << "    PerformAlphaKill((int)DeriveAlphaKill()." << swizzle[stage] << ", " << tReg << ");\n";
    }
}

// ============================================================
// Dead-code collection — scan all combiner + final inputs
// ============================================================

// Scan one combiner or final-combiner input byte and mark the register + channel it reads.
// Also sets the per-stage C0/C1 bits when the input references those regs.
// For combiner stages: isFinal=false, stage=0..7.
// For final combiner:  isFinal=true, isFinalABCD distinguishes EFG vs ABCD phase, stage=8.
static void ScanInput(ReadMasks& m, uint32_t inputByte, uint32_t stage,
                      bool isFinal = false, bool isFinalABCD = false)
{
    uint32_t regIdx = inputByte & PS_REGISTER_MASK;
    if (regIdx == PS_REGISTER_ZERO) return;

    // V1R0_SUM / EF_PROD only valid in final combiner ABCD phase
    if (isFinal && !isFinalABCD && (regIdx == PS_REGISTER_V1R0_SUM || regIdx == PS_REGISTER_EF_PROD))
        return;

    // Bit 4: for RGB inputs selects .aaaa swizzle; for alpha inputs selects .a vs .b
    // .b is the blue channel → part of RGB, not alpha
    bool useAlpha = (inputByte & 0x10) != 0;
    uint32_t ch = isFinal ? (useAlpha ? CH_A : CH_RGBA)
                          : (useAlpha ? CH_A : CH_RGB);

    SetChan(m.regChans, regIdx, ch);

    // Track per-stage C0/C1 usage
    if (regIdx == PS_REGISTER_C0) m.stageC0 |= (1u << stage);
    if (regIdx == PS_REGISTER_C1) m.stageC1 |= (1u << stage);
}

// Compute which T registers are needed by texture-mode cross-stage dependencies.
// Walk stages in reverse so later-stage dependencies propagate to earlier stages.
static uint16_t ScanTextureDeps(const PSJITKey& key, const uint32_t texModes[4],
                                uint16_t combinerNeeded)
{
    uint16_t needed = combinerNeeded; // start with what combiners need

    // Reverse pass — if stage i is needed, mark its source dependencies
    for (int i = 3; i >= 0; i--) {
        if (!(needed & (1u << i))) continue;

        uint32_t mode = texModes[i];
        uint32_t srcStage;
        switch (mode) {
        case PS_TEXTUREMODES_BUMPENVMAP: case PS_TEXTUREMODES_BUMPENVMAP_LUM:
        case PS_TEXTUREMODES_DPNDNT_AR: case PS_TEXTUREMODES_DPNDNT_GB:
        case PS_TEXTUREMODES_DOTPRODUCT:
        case PS_TEXTUREMODES_DOT_ST: case PS_TEXTUREMODES_DOT_ZW:
            srcStage = GetTexSrcStage(i, key.shaderCtl);
            needed |= (1u << srcStage);
            // DOT_ST, DOT_ZW also read T[stage-1]
            if ((mode == PS_TEXTUREMODES_DOT_ST || mode == PS_TEXTUREMODES_DOT_ZW) && i >= 1)
                needed |= (1u << (i - 1));
            break;

        case PS_TEXTUREMODES_DOT_RFLCT_DIFF:
        {
            srcStage = GetTexSrcStage(i, key.shaderCtl);
            needed |= (1u << srcStage);
            if (i >= 1) needed |= (1u << (i - 1));
            // Peeks at next stage's source
            uint32_t nextStage = i + 1;
            if (nextStage < 4) {
                needed |= (1u << GetTexSrcStage(nextStage, key.shaderCtl));
                needed |= (1u << nextStage);
            }
            break;
        }

        case PS_TEXTUREMODES_DOT_RFLCT_SPEC: case PS_TEXTUREMODES_DOT_STR_3D: case PS_TEXTUREMODES_DOT_STR_CUBE: case PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST:
            srcStage = GetTexSrcStage(i, key.shaderCtl);
            needed |= (1u << srcStage);
            if (i >= 1) needed |= (1u << (i - 1));
            if (i >= 2) needed |= (1u << (i - 2));
            break;

        default:
            break;
        }
    }
    return needed;
}

// Collect all read masks from the full combiner + final combiner configuration.
static ReadMasks CollectReadRegs(const PSJITKey& key)
{
    ReadMasks m = {};
    uint32_t numStages = key.numStages;
    uint32_t ccFlags = key.combinectl >> 8;

    // --- Combiner stage inputs ---
    for (uint32_t stage = 0; stage < numStages; stage++) {
        uint32_t rgbIn = key.rgbInputs[stage];
        uint32_t aIn   = key.alphaInputs[stage];

        for (int slot = 0; slot < 4; slot++) {
            ScanInput(m, InputByte(rgbIn, slot), stage);
            ScanInput(m, InputByte(aIn, slot), stage);
        }

        // MUX implicitly reads R0.a
        uint32_t rgbFlags = key.rgbOutputs[stage] >> 12;
        uint32_t aFlags   = key.alphaOutputs[stage] >> 12;
        if ((rgbFlags & PS_COMBINEROUTPUT_AB_CD_MUX) || (aFlags & PS_COMBINEROUTPUT_AB_CD_MUX))
            SetChan(m.regChans, PS_REGISTER_R0, CH_A); // R0.a
    }

    // --- Final combiner inputs ---
    uint32_t fcABCD = key.fcABCD;
    uint32_t fcEFG  = key.fcEFG;

    if (fcABCD != 0 || fcEFG != 0) {
        // EFG phase (not ABCD — V1R0_SUM/EF_PROD not yet valid)
        for (int slot = 0; slot < 3; slot++)
            ScanInput(m, InputByte(fcEFG, slot), 8, true, false);

        // ABCD phase
        for (int slot = 0; slot < 4; slot++)
            ScanInput(m, InputByte(fcABCD, slot), 8, true, true);

        // V1R0_SUM transitively reads V1.rgb and R0.rgb
        if (GetChan(m.regChans, PS_REGISTER_V1R0_SUM)) {
            SetChan(m.regChans, PS_REGISTER_V1, CH_RGB);
            SetChan(m.regChans, PS_REGISTER_R0, CH_RGB);
        }
    } else {
        // No final combiner → output is R0
        SetChan(m.regChans, PS_REGISTER_R0, CH_RGBA);
    }

    // --- Fog reads FOG ---
    if (key.fogEnable)
        SetChan(m.regChans, PS_REGISTER_FOG, CH_A);

    // --- Texture cross-stage dependencies ---
    // Start with T registers referenced by combiners
    uint16_t combinerTexNeeded = 0;
    for (uint32_t i = 0; i < 4; i++) {
        if (GetChan(m.regChans, PS_REGISTER_T0 + i))
            combinerTexNeeded |= (1u << i);
    }

    uint32_t texModes[4];
    UnpackTexModes(key.textureModes, texModes);

    // Also mark stages with side effects (CLIPPLANE, DOT_ZW) as needed
    for (uint32_t i = 0; i < 4; i++) {
        if (texModes[i] == PS_TEXTUREMODES_CLIPPLANE || // CLIPPLANE — side effect
            texModes[i] == PS_TEXTUREMODES_DOT_ZW)      // DOT_ZW — writes fragment depth
            combinerTexNeeded |= (1u << i);
    }

    m.texNeeded = ScanTextureDeps(key, texModes, combinerTexNeeded);

    return m;
}

// ============================================================
// Generate complete HLSL source
// ============================================================

static std::string GenerateHLSL(const PSJITKey& key)
{
    std::ostringstream ss;

    uint32_t numStages = key.numStages;
    uint32_t ccFlags = key.combinectl >> 8;
    bool flagMuxMsb   = (ccFlags & PS_COMBINERCOUNT_MUX_MSB) != 0;
    bool flagUniqueC0 = (ccFlags & PS_COMBINERCOUNT_UNIQUE_C0) != 0;
    bool flagUniqueC1 = (ccFlags & PS_COMBINERCOUNT_UNIQUE_C1) != 0;

    // Determine which texture types are needed per stage
    uint32_t texModes[4];
    UnpackTexModes(key.textureModes, texModes);

    // Dead-code analysis
    ReadMasks masks = CollectReadRegs(key);
    auto isRead     = [&](uint32_t regIdx) { return GetChan(masks.regChans, regIdx) != CH_NONE; };
    auto isReadChan = [&](uint32_t regIdx, uint32_t ch) { return (GetChan(masks.regChans, regIdx) & ch) != 0; };

    // DOT_ZW (z-sprite) produces per-pixel depth → needs SV_Depth output
    bool hasDotZW = (texModes[0] == 0x0A || texModes[1] == 0x0A ||
                     texModes[2] == 0x0A || texModes[3] == 0x0A);

    // ---- Header / resource declarations ----
    ss << "// Auto-generated by CxbxPixelShaderJIT\n";
    ss << "// Stages: " << numStages << " TexModes: "
       << texModes[0] << "/" << texModes[1] << "/" << texModes[2] << "/" << texModes[3] << "\n\n";

    // PGRAPH register accessors and color unpacking — shared with RC interpreter
    ss << "#include \"CxbxPGRAPHRegs.hlsli\"\n\n";

    // Texture/sampler declarations — only what's needed
    for (uint32_t i = 0; i < 4; i++) {
        bool need2D = false, need3D = false, needCube = false;
        switch (texModes[i]) {
        case PS_TEXTUREMODES_PROJECT2D: case PS_TEXTUREMODES_BUMPENVMAP: case PS_TEXTUREMODES_BUMPENVMAP_LUM: case PS_TEXTUREMODES_BRDF: case PS_TEXTUREMODES_DOT_ST:
        case PS_TEXTUREMODES_DPNDNT_AR: case PS_TEXTUREMODES_DPNDNT_GB:
            need2D = true; break;
        case PS_TEXTUREMODES_PROJECT3D: case PS_TEXTUREMODES_DOT_STR_3D:
            need3D = true; break;
        case PS_TEXTUREMODES_CUBEMAP: case PS_TEXTUREMODES_DOT_RFLCT_DIFF: case PS_TEXTUREMODES_DOT_RFLCT_SPEC: case PS_TEXTUREMODES_DOT_STR_CUBE: case PS_TEXTUREMODES_DOT_RFLCT_SPEC_CONST:
            needCube = true; break;
        }
        // PROJECT3D with shadow uses 2D
        if (texModes[i] == PS_TEXTUREMODES_PROJECT3D && key.shadowCompare[i] != 0.0f)
            need2D = true, need3D = false;

        if (need2D)   ss << "Texture2D<float4>   Tex2D_" << i << "   : register(t" << i << ");\n";
        if (need3D)   ss << "Texture3D<float4>   Tex3D_" << i << "   : register(t" << (i+4) << ");\n";
        if (needCube) ss << "TextureCube<float4> TexCube_" << i << " : register(t" << (i+8) << ");\n";
        if (texModes[i] != PS_TEXTUREMODES_NONE)
            ss << "SamplerState Samp" << i << " : register(s" << i << ");\n";
    }
    // Stencil SRVs for depth-as-color remapping (X24_TYPELESS_G8_UINT, slots t16-t19)
    for (uint32_t i = 0; i < 4; i++)
        ss << "Texture2D<uint2> TexStencil_" << i << " : register(t" << (16 + i) << ");\n";
    ss << "\n";

    // Aux cbuffer — shared layout with RC interpreter
    ss << "#include \"CxbxRegisterCombinerInterpreterState.hlsli\"\n\n";

    // PS_INPUT — shared with VS output and RC interpreter
    ss << "#include \"CxbxPixelShaderInput.hlsli\"\n\n";

    // Shared helper functions (nv2a_mul, PerformColorSign, ApplyShadowCompare,
    // ApplyCompareMode, RemapD24S8ToColor, RemapD16ToColor, etc.)
    ss << "#include \"CxbxNV2AMathHelpers.hlsli\"\n";
    ss << "#include \"CxbxPixelShaderFunctions.hlsli\"\n";
    ss << "#include \"CxbxPSAuxFromPGRAPH.hlsli\"\n\n";

    // ApplyShadowCompare, ApplyCompareMode, and ApplyDotMapping are all
    // provided by CxbxPixelShaderFunctions.hlsli (included above)

    // ---- main() ----
    if (hasDotZW) {
        ss << "struct PS_OUTPUT {\n";
        ss << "    float4 color : SV_Target;\n";
        ss << "    float depth : SV_Depth;\n";
        ss << "};\n\n";
        ss << "PS_OUTPUT main(PS_INPUT input)\n{\n";
    } else {
        ss << "float4 main(PS_INPUT input) : SV_Target\n{\n";
    }

    // Declare register variables
    ss << "    float4 T0 = input.iT0;\n";
    ss << "    float4 T1 = input.iT1;\n";
    ss << "    float4 T2 = input.iT2;\n";
    ss << "    float4 T3 = input.iT3;\n";

    // Per-pixel depth variable for DOT_ZW (z-sprite)
    if (hasDotZW) {
        ss << "    float fragDepth = 0.0;\n";
    }

    // Eye vector for DOT_RFLCT_SPEC
    bool needsEyeVec = false;
    for (int i = 0; i < 4; i++)
        if (texModes[i] == PS_TEXTUREMODES_DOT_RFLCT_SPEC) needsEyeVec = true;
    if (needsEyeVec)
        ss << "    float3 eyeVec = float3(input.iT1.w, input.iT2.w, input.iT3.w);\n";

    // Texture fetches (sequential — later stages can depend on earlier)
    ss << "\n    // --- Texture fetches ---\n";
    for (uint32_t i = 0; i < 4; i++) {
        if (texModes[i] != PS_TEXTUREMODES_NONE && (masks.texNeeded & (1u << i))) {
            float cs4[4] = { key.colorSign[i*4+0], key.colorSign[i*4+1],
                             key.colorSign[i*4+2], key.colorSign[i*4+3] };
            EmitTextureFetch(ss, i, texModes[i], key.shaderCtl,
                             key.shadowCompare[i], cs4, key.texFmtFixup[i],
                             key.colorKeyOp[i], key.alphaKill[i]);
        }
    }

    // Vertex colors with front/back face selection
    ss << "\n    // --- Vertex colors ---\n";
    if (key.frontFaceFactor != 0.0f && (isRead(PS_REGISTER_V0) || isRead(PS_REGISTER_V1))) {
        ss << "    float faceSign = input.iFF ? 1.0 : -1.0;\n";
        ss << "    bool isFront = (faceSign * DeriveFrontFaceFactor()) >= 0.0;\n";
        if (isRead(PS_REGISTER_V0))
            ss << "    float4 V0 = isFront ? input.iD0 : input.iB0;\n";
        if (isRead(PS_REGISTER_V1))
            ss << "    float4 V1 = isFront ? input.iD1 : input.iB1;\n";
    } else {
        if (isRead(PS_REGISTER_V0))
            ss << "    float4 V0 = input.iD0;\n";
        if (isRead(PS_REGISTER_V1))
            ss << "    float4 V1 = input.iD1;\n";
    }
    if (isRead(PS_REGISTER_FOG))
        ss << "    float4 FOG = float4(PG_COLOR_ARGB(0x1980).rgb, saturate(input.iFog));\n";
    ss << "    float4 R0 = float4(0.0, 0.0, 0.0, T0.a);\n";
    if (isRead(PS_REGISTER_R1))
        ss << "    float4 R1 = (float4)0.0;\n";
    if (isRead(PS_REGISTER_V1R0_SUM))
        ss << "    float4 V1R0_SUM = (float4)0.0;\n";
    if (isRead(PS_REGISTER_EF_PROD))
        ss << "    float4 EF_PROD = (float4)0.0;\n";
    // Declare C0/C1 only if any stage or the final combiner references them
    bool anyC0 = (masks.stageC0 != 0);
    bool anyC1 = (masks.stageC1 != 0);
    if (anyC0) ss << "    float4 C0;\n";
    if (anyC1) ss << "    float4 C1;\n";

    // Hoist shared (non-unique) C0/C1 loads before the stage loop
    if (!flagUniqueC0 && anyC0 && (masks.stageC0 & 0xFF))
        ss << "    C0 = PG_COLOR_ARGB(0x" << std::hex << 0x1880 << ");\n" << std::dec;
    if (!flagUniqueC1 && anyC1 && (masks.stageC1 & 0xFF))
        ss << "    C1 = PG_COLOR_ARGB(0x" << std::hex << 0x18A0 << ");\n" << std::dec;

    // Combiner stages
    ss << "\n    // --- Combiner stages ---\n";
    for (uint32_t stage = 0; stage < numStages; stage++) {
        ss << "    // Stage " << stage << "\n";

        // Load per-stage C0/C1 from PGRAPH (only when unique per stage)
        if (flagUniqueC0 && (masks.stageC0 & (1u << stage))) {
            uint32_t c0Off = 0x1880 + stage * 4;
            ss << "    C0 = PG_COLOR_ARGB(0x" << std::hex << c0Off << ");\n" << std::dec;
        }
        if (flagUniqueC1 && (masks.stageC1 & (1u << stage))) {
            uint32_t c1Off = 0x18A0 + stage * 4;
            ss << "    C1 = PG_COLOR_ARGB(0x" << std::hex << c1Off << ");\n" << std::dec;
        }

        // Decode output control
        uint32_t rgbOut = key.rgbOutputs[stage];
        uint32_t aOut   = key.alphaOutputs[stage];
        uint32_t rgbFlags = rgbOut >> 12;
        uint32_t aFlags   = aOut >> 12;
        bool flagCDDot  = (rgbFlags & PS_COMBINEROUTPUT_CD_DOT_PRODUCT)   != 0;
        bool flagABDot  = (rgbFlags & PS_COMBINEROUTPUT_AB_DOT_PRODUCT)   != 0;
        bool flagRGBMux = (rgbFlags & PS_COMBINEROUTPUT_AB_CD_MUX)        != 0;
        bool flagAMux   = (aFlags   & PS_COMBINEROUTPUT_AB_CD_MUX)        != 0;
        bool cdBlue2A   = (rgbFlags & PS_COMBINEROUTPUT_CD_BLUE_TO_ALPHA) != 0;
        bool abBlue2A   = (rgbFlags & PS_COMBINEROUTPUT_AB_BLUE_TO_ALPHA) != 0;

        // Output destinations
        uint32_t rgbRegAB  = (rgbOut >> 4) & PS_REGISTER_MASK;
        uint32_t rgbRegCD  = (rgbOut     ) & PS_REGISTER_MASK;
        uint32_t rgbRegSum = (rgbOut >> 8) & PS_REGISTER_MASK;
        uint32_t aRegAB    = (aOut   >> 4) & PS_REGISTER_MASK;
        uint32_t aRegCD    = (aOut       ) & PS_REGISTER_MASK;
        uint32_t aRegSum   = (aOut   >> 8) & PS_REGISTER_MASK;

        // SUM output only exists when neither AB nor CD use dot product
        bool writeSumRGB = !flagABDot && !flagCDDot;

        // Whole-stage skip: if none of the output destinations are read, skip
        bool anyOutRead = isReadChan(rgbRegAB, CH_RGB) || isReadChan(rgbRegCD, CH_RGB)
                       || (writeSumRGB && isReadChan(rgbRegSum, CH_RGB))
                       || isReadChan(aRegAB, CH_A)    || isReadChan(aRegCD, CH_A)
                       || isReadChan(aRegSum, CH_A);
        if (!anyOutRead) {
            ss << "    // (stage " << stage << " skipped — outputs not read)\n";
            continue;
        }

        // Determine which halves of the stage are needed
        bool needSum_rgb = writeSumRGB && isReadChan(rgbRegSum, CH_RGB);
        bool needSum_a   = isReadChan(aRegSum, CH_A);
        bool needAB_rgb  = isReadChan(rgbRegAB, CH_RGB) || needSum_rgb;
        bool needCD_rgb  = isReadChan(rgbRegCD, CH_RGB) || needSum_rgb;
        bool needAB_a    = isReadChan(aRegAB, CH_A) || needSum_a;
        bool needCD_a    = isReadChan(aRegCD, CH_A) || needSum_a;
        bool needRGB     = needAB_rgb || needCD_rgb;
        bool needA       = needAB_a || needCD_a;
        bool needMux     = (flagRGBMux && needSum_rgb) || (flagAMux && needSum_a);

        // Decode inputs
        uint32_t rgbIn = key.rgbInputs[stage];
        uint32_t aIn   = key.alphaInputs[stage];

        // Compute AB, CD
        ss << "    {\n";

        if (needRGB) {
            std::string rgbA = EmitInput(InputByte(rgbIn, 0), false);
            std::string rgbB = EmitInput(InputByte(rgbIn, 1), false);
            std::string rgbC = EmitInput(InputByte(rgbIn, 2), false);
            std::string rgbD = EmitInput(InputByte(rgbIn, 3), false);

            if (needAB_rgb) {
                ss << "        float4 rgbA = " << rgbA << ";\n";
                ss << "        float4 rgbB = " << rgbB << ";\n";
                if (flagABDot)
                    ss << "        float3 rgbAB = (float3)dot(rgbA.rgb, rgbB.rgb);\n";
                else
                    ss << "        float3 rgbAB = nv2a_mul3(rgbA.rgb, rgbB.rgb);\n";
            }
            if (needCD_rgb) {
                ss << "        float4 rgbC = " << rgbC << ";\n";
                ss << "        float4 rgbD = " << rgbD << ";\n";
                if (flagCDDot)
                    ss << "        float3 rgbCD = (float3)dot(rgbC.rgb, rgbD.rgb);\n";
                else
                    ss << "        float3 rgbCD = nv2a_mul3(rgbC.rgb, rgbD.rgb);\n";
            }
        }

        if (needA) {
            std::string aA = EmitInput(InputByte(aIn, 0), true);
            std::string aB = EmitInput(InputByte(aIn, 1), true);
            std::string aC = EmitInput(InputByte(aIn, 2), true);
            std::string aD = EmitInput(InputByte(aIn, 3), true);

            if (needAB_a) {
                ss << "        float aA = " << aA << ";\n";
                ss << "        float aB = " << aB << ";\n";
                ss << "        float aAB = nv2a_mul1(aA, aB);\n";
            }
            if (needCD_a) {
                ss << "        float aC = " << aC << ";\n";
                ss << "        float aD = " << aD << ";\n";
                ss << "        float aCD = nv2a_mul1(aC, aD);\n";
            }
        }

        // MUX selector
        if (needMux) {
            if (flagMuxMsb)
                ss << "        bool muxSel = (R0.a >= 0.5);\n";
            else
                ss << "        bool muxSel = (((uint)(saturate(R0.a)*255.0+0.5) & 1u) != 0u);\n";
        }

        // SUM/MUX results
        if (needSum_rgb) {
            if (flagRGBMux)
                ss << "        float3 rgbSUM = muxSel ? rgbCD : rgbAB;\n";
            else
                ss << "        float3 rgbSUM = rgbAB + rgbCD;\n";
        }

        if (needSum_a) {
            if (flagAMux)
                ss << "        float aSUM = muxSel ? aCD : aAB;\n";
            else
                ss << "        float aSUM = aAB + aCD;\n";
        }

        // Output mapping
        float rgbBias = 0.0f, rgbScale = 1.0f;
        float aBias = 0.0f, aScale = 1.0f;
        if (needRGB) {
            if (rgbFlags & PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS) rgbBias = -0.5f;
            switch ((rgbFlags >> 4) & 3) {
            case 1: rgbScale = 2.0f; break;
            case 2: rgbScale = 4.0f; break;
            case 3: rgbScale = 0.5f; break;
            }
        }
        if (needA) {
            if (aFlags & PS_COMBINEROUTPUT_OUTPUTMAPPING_BIAS) aBias = -0.5f;
            switch ((aFlags >> 4) & 3) {
            case 1: aScale = 2.0f; break;
            case 2: aScale = 4.0f; break;
            case 3: aScale = 0.5f; break;
            }
        }

        // Emit output expressions
        auto emitClamp = [&](const std::string& expr, float bias, float scale) -> std::string {
            std::ostringstream out;
            if (bias == 0.0f && scale == 1.0f) {
                out << "clamp(" << expr << ", -1.0, 1.0)";
            } else if (bias == 0.0f) {
                out << "clamp((" << expr << ") * " << scale << ", -1.0, 1.0)";
            } else {
                out << "clamp((" << expr << " + " << bias << ") * " << scale << ", -1.0, 1.0)";
            }
            return out.str();
        };

        // RGB writes
        auto emitRGBWrite = [&](uint32_t dest, const std::string& rgbExpr,
                                bool blue2Alpha, uint32_t aRegSame) {
            if (dest == PS_REGISTER_DISCARD) return;
            const char* dName = RegName(dest);
            if (!dName) return;
            bool aOverwrite = (aRegSame == dest) && (aRegSame != PS_REGISTER_DISCARD);
            if (blue2Alpha && !aOverwrite) {
                ss << "        { float3 _rgb = " << rgbExpr << ";\n";
                ss << "          " << dName << " = float4(_rgb, _rgb.b); }\n";
            } else {
                ss << "        " << dName << ".rgb = " << rgbExpr << ";\n";
            }
        };

        if (needRGB) {
            if (needAB_rgb && isReadChan(rgbRegAB, CH_RGB))
                emitRGBWrite(rgbRegAB, emitClamp("rgbAB", rgbBias, rgbScale), abBlue2A, aRegAB);
            if (needCD_rgb && isReadChan(rgbRegCD, CH_RGB))
                emitRGBWrite(rgbRegCD, emitClamp("rgbCD", rgbBias, rgbScale), cdBlue2A, aRegCD);
            if (needSum_rgb && rgbRegSum != PS_REGISTER_DISCARD && isReadChan(rgbRegSum, CH_RGB)) {
                const char* dName = RegName(rgbRegSum);
                if (dName)
                    ss << "        " << dName << ".rgb = " << emitClamp("rgbSUM", rgbBias, rgbScale) << ";\n";
            }
        }

        // Alpha writes
        if (needA) {
            auto emitAlphaWrite = [&](uint32_t dest, const std::string& aExpr) {
                if (dest == PS_REGISTER_DISCARD) return;
                if (!isReadChan(dest, CH_A)) return;
                const char* dName = RegName(dest);
                if (!dName) return;
                ss << "        " << dName << ".a = " << aExpr << ";\n";
            };

            if (needAB_a)  emitAlphaWrite(aRegAB,  emitClamp("aAB", aBias, aScale));
            if (needCD_a)  emitAlphaWrite(aRegCD,  emitClamp("aCD", aBias, aScale));
            if (needSum_a) emitAlphaWrite(aRegSum, emitClamp("aSUM", aBias, aScale));
        }

        ss << "    }\n";
    }

    // ---- Final Combiner ----
    ss << "\n    // --- Final Combiner ---\n";
    uint32_t fcABCD = key.fcABCD;
    uint32_t fcEFG  = key.fcEFG;

    if (fcABCD == 0 && fcEFG == 0) {
        ss << "    float4 result = R0;\n";
    } else {
        // Load final combiner C0/C1
        if (masks.stageC0 & (1u << 8))
            ss << "    C0 = PG_COLOR_ARGB(0x19AC);\n";
        if (masks.stageC1 & (1u << 8))
            ss << "    C1 = PG_COLOR_ARGB(0x19B0);\n";

        // EFG phase
        uint32_t settings = (fcEFG) & 0xFF;
        uint32_t eReg = InputByte(fcEFG, 0);
        uint32_t fReg = InputByte(fcEFG, 1);
        uint32_t gReg = InputByte(fcEFG, 2);

        ss << "    float  fcG = (" << EmitFinalInput(gReg, false) << ").a;\n";

        if (isRead(PS_REGISTER_EF_PROD)) {
            ss << "    float3 fcE = (" << EmitFinalInput(eReg, false) << ").rgb;\n";
            ss << "    float3 fcF = (" << EmitFinalInput(fReg, false) << ").rgb;\n";
            ss << "    EF_PROD = float4(fcE * fcF, 1.0);\n";
        }

        // V1R0_SUM — only if referenced by ABCD inputs
        if (isRead(PS_REGISTER_V1R0_SUM)) {
            ss << "    {\n";
            ss << "        float3 v1s = V1.rgb;\n";
            ss << "        float3 r0s = R0.rgb;\n";
            if (settings & 0x40)
                ss << "        v1s = 1.0 - v1s;\n";
            if (settings & 0x20)
                ss << "        r0s = 1.0 - r0s;\n";
            ss << "        float3 sum = v1s + r0s;\n";
            if (settings & 0x80)
                ss << "        sum = saturate(sum);\n";
            ss << "        V1R0_SUM = float4(sum, 1.0);\n";
            ss << "    }\n";
        }

        // ABCD phase
        uint32_t aReg = InputByte(fcABCD, 0);
        uint32_t bReg = InputByte(fcABCD, 1);
        uint32_t cReg = InputByte(fcABCD, 2);
        uint32_t dReg = InputByte(fcABCD, 3);

        ss << "    float4 fcA = " << EmitFinalInput(aReg, true) << ";\n";
        ss << "    float4 fcB = " << EmitFinalInput(bReg, true) << ";\n";
        ss << "    float4 fcC = " << EmitFinalInput(cReg, true) << ";\n";
        ss << "    float4 fcD = " << EmitFinalInput(dReg, true) << ";\n";

        ss << "    float4 result;\n";
        ss << "    result.rgb = saturate(lerp(fcC.rgb, fcB.rgb, fcA.rgb) + fcD.rgb);\n";
        ss << "    result.a = fcG;\n";
    }

    // Alpha test
    ss << "\n    // --- Alpha test ---\n";
    ss << "    { uint c0 = PG_UINT(0x194C);\n";
    ss << "      float ae = (c0 & 0x1000u) ? 1.0 : 0.0;\n";
    ss << "      float ar = float(c0 & 0xFFu) / 255.0;\n";
    ss << "      float af = float((c0 & 0xF00u) >> 8);\n";
    ss << "      PerformAlphaTest(float3(ae, ar, af), result.a); }\n";

    // Fog blending
    if (key.fogEnable) {
        ss << "\n    // --- Fog ---\n";
        ss << "    result.rgb = lerp(PG_COLOR_ARGB(0x1980).rgb, result.rgb, saturate(input.iFog));\n";
    }

    if (hasDotZW) {
        ss << "\n    PS_OUTPUT psOut;\n";
        ss << "    psOut.color = result;\n";
        ss << "    psOut.depth = saturate(fragDepth / DepthScale.x);\n";
        ss << "    return psOut;\n}\n";
    } else {
        ss << "\n    return result;\n}\n";
    }

    return ss.str();
}

// ============================================================
// JIT entry point
// ============================================================

ID3D11PixelShader* PixelShaderCache::GetShader(ID3D11Device* pDevice)
{
    extern NV2ADevice* g_NV2A;
    if (!g_NV2A) return nullptr;

    PGRAPHState* pg = &g_NV2A->GetDeviceState()->pgraph;

    // Fast path: if none of the relevant dirty groups changed since last
    // call, the combiner topology is identical — return cached result
    // directly without rebuilding the key or hashing.
    // PS JIT reads SHADER, TEXTURE, BLEND, and RASTERIZER registers.
    static uint32_t s_LastPSRegsGen = ~0u;
    static ID3D11PixelShader* s_LastPSResult = nullptr;
    static PSJITKey s_LastKey = {};
    uint32_t psRegsGen = pg->dirty[NV2A_DIRTY_SHADER] + pg->dirty[NV2A_DIRTY_TEXTURE]
                       + pg->dirty[NV2A_DIRTY_BLEND] + pg->dirty[NV2A_DIRTY_RASTERIZER];
    if (psRegsGen == s_LastPSRegsGen)
        return s_LastPSResult;

    // Capture current state — read directly from PGRAPH registers
    // (no dependency on g_LastPSAuxCB or CxbxD3D11UploadRCInterpreterState)
    PSJITKey key = {};
    key.cacheVersion = PS_JIT_CACHE_VERSION;
    key.combinectl = pg->regs[0x1940 >> 2];
    key.numStages = key.combinectl & 0xFF;
    if (key.numStages == 0) key.numStages = 1;
    if (key.numStages > 8)  key.numStages = 8;

    for (uint32_t i = 0; i < 8; i++) {
        key.rgbInputs[i]   = pg->regs[(0x1900 + i * 4) >> 2];
        key.alphaInputs[i] = pg->regs[(0x18C0 + i * 4) >> 2];
        key.rgbOutputs[i]  = pg->regs[(0x1920 + i * 4) >> 2];
        key.alphaOutputs[i]= pg->regs[(0x18E0 + i * 4) >> 2];
    }

    key.shaderCtl      = pg->regs[0x1998 >> 2];
    key.shaderClipMode = pg->regs[0x1994 >> 2];

    // Read from the last-built aux CB (already computed by CxbxD3D11UploadRCInterpreterState)
    const PSAuxCBLayout& aux = g_LastPSAuxCB;

    // Read textureModes directly from PGRAPH (NV_PGRAPH_SHADERPROG = 0x199C)
    // to avoid stale g_LastPSAuxCB issues.
    key.textureModes = pg->regs[0x199C >> 2]; // NV_PGRAPH_SHADERPROG

    // Final combiner inputs — derive from PGRAPH (synthesize default if not set)
    {
        uint32_t fcABCD = pg->regs[0x1944 >> 2]; // COMBINESPECFOG0
        uint32_t fcEFG  = pg->regs[0x1948 >> 2]; // COMBINESPECFOG1
        if (fcABCD == 0 && fcEFG == 0) {
            bool fogEnable = (pg->regs[0x1958 >> 2] & 0x100) != 0;
            bool specEnable = (pg->regs[0x0FB8 >> 2] & 0x10000) != 0;
            // NV2A register encoding: FOG=0x03, R0=0x0C, V1=0x05, ZERO=0x00, ALPHA=0x10
            uint32_t regA = 0x13; // FOG | CHANNEL_ALPHA
            uint32_t regB = 0x0C; // R0
            uint32_t regC = fogEnable ? 0x03u : 0x0Cu; // FOG or R0
            uint32_t regD = specEnable ? 0x05u : 0x00u; // V1 or ZERO
            fcABCD = (regA << 24) | (regB << 16) | (regC << 8) | regD;
            fcEFG = (0x00u << 24) | (0x00u << 16) | (0x1Cu << 8); // ZERO, ZERO, R0|ALPHA
        }
        key.fcABCD = fcABCD;
        key.fcEFG  = fcEFG;
    }

    // Fog enable — from PGRAPH CONTROL_3
    key.fogEnable = (pg->regs[0x1958 >> 2] & 0x100) ? 1u : 0u;

    // Front face factor — from PGRAPH CSV0_C + SETUPRASTER
    {
        uint32_t csv0c = pg->regs[0x0FB8 >> 2];
        bool twoSided = (csv0c & 0x20000000u) != 0;
        if (twoSided) {
            bool ccwFront = (pg->regs[0x1990 >> 2] & 0x00800000u) != 0;
            key.frontFaceFactor = ccwFront ? -1.0f : 1.0f;
        } else {
            key.frontFaceFactor = 0.0f;
        }
    }

    // Color key, alpha kill, shadow compare — from PGRAPH TEXCTL0/TEXFMT0
    for (int i = 0; i < 4; i++) {
        uint32_t texCtl = pg->regs[(0x19CC + i * 4) >> 2];
        key.colorKeyOp[i]    = (float)(texCtl & 0x03u);
        key.alphaKill[i]     = (texCtl & 0x04u) ? 1.0f : 0.0f;

        // Shadow compare: check if TEXFMT COLOR field is a depth format
        float sc = 0.0f;
        if (texCtl & 0x40000000u) { // ENABLE
            uint32_t texFmt = pg->regs[(0x1A04 + i * 4) >> 2];
            uint32_t colorCode = (texFmt >> 8) & 0x7Fu;
            if ((colorCode >= 0x2A && colorCode <= 0x2D) || (colorCode >= 0x2E && colorCode <= 0x31))
                sc = 1.0f;
        }
        key.shadowCompare[i] = sc;

        key.texFmtFixup[i]   = (&aux.TexFmtFixup.x)[i];
        key.colorSign[i*4+0] = aux.ColorSign[i].x;
        key.colorSign[i*4+1] = aux.ColorSign[i].y;
        key.colorSign[i*4+2] = aux.ColorSign[i].z;
        key.colorSign[i*4+3] = aux.ColorSign[i].w;
    }

    // Second fast path: if the key matches the last one (combiner state unchanged
    // despite dirty groups bumping from non-combiner register writes),
    // skip the expensive hash + mutex + map lookup.
    if (memcmp(&key, &s_LastKey, sizeof(PSJITKey)) == 0) {
        s_LastPSRegsGen = psRegsGen;
        return s_LastPSResult;
    }

    // Hash and cache lookup
    uint64_t hash = HashKey(key);
    {
        std::lock_guard<std::mutex> lock(g_PSJITMutex);
        auto it = g_PSJITCache.find(hash);
        if (it != g_PSJITCache.end()) {
            s_LastKey = key;
            s_LastPSRegsGen = psRegsGen;
            s_LastPSResult = it->second.pPS;
            return s_LastPSResult; // nullptr = known failure
        }
    }

    EmuLog(LOG_LEVEL::DEBUG, "PS JIT: new key hash=%016llX stages=%u texModes=0x%08X fcABCD=0x%08X fcEFG=0x%08X combinectl=0x%08X",
           hash, key.numStages, key.textureModes, key.fcABCD, key.fcEFG, key.combinectl);

    // Generate HLSL
    InterlockedIncrement(&g_ProfilePSJITCompiles);
    CXBX_PROFILE_SCOPE(PROF_PS_JIT_COMPILE);

    // Try disk cache first
    ID3DBlob* pCode = ShaderDiskCache::TryLoad(hash);
    if (pCode) {
        ID3D11PixelShader* pPS = nullptr;
        HRESULT hr = pDevice->CreatePixelShader(pCode->GetBufferPointer(), pCode->GetBufferSize(), nullptr, &pPS);
        pCode->Release();
        if (SUCCEEDED(hr)) {
            std::lock_guard<std::mutex> lock(g_PSJITMutex);
            g_PSJITCache[hash] = { pPS };
            s_LastKey = key;
            s_LastPSRegsGen = psRegsGen;
            s_LastPSResult = pPS;
            return pPS;
        }
        // Fall through to recompile if cached blob is invalid
    }

    std::string hlsl = GenerateHLSL(key);

    // Compile
    pCode = nullptr;
    ID3DBlob* pErrors = nullptr;
    CxbxJITIncludeHandler includeHandler;
    HRESULT hr = D3DCompile(
        hlsl.c_str(), hlsl.size(),
        "CxbxPixelShaderJIT", nullptr, &includeHandler,
        "main", "ps_5_0",
        D3DCOMPILE_OPTIMIZATION_LEVEL3, 0,
        &pCode, &pErrors);

    if (FAILED(hr)) {
        EmuLog(LOG_LEVEL::WARNING, "PS JIT compile failed (0x%08X)", hr);
        if (pErrors) {
            EmuLog(LOG_LEVEL::WARNING, "PS JIT errors: %s", (const char*)pErrors->GetBufferPointer());
            pErrors->Release();
        }
        // Cache the failure
        std::lock_guard<std::mutex> lock(g_PSJITMutex);
        g_PSJITCache[hash] = { nullptr };
        s_LastKey = key;
        s_LastPSRegsGen = psRegsGen;
        s_LastPSResult = nullptr;
        return nullptr;
    }
    if (pErrors) pErrors->Release();

    // Create pixel shader
    ID3D11PixelShader* pPS = nullptr;
    hr = pDevice->CreatePixelShader(pCode->GetBufferPointer(), pCode->GetBufferSize(), nullptr, &pPS);

    if (FAILED(hr)) {
        pCode->Release();
        EmuLog(LOG_LEVEL::WARNING, "PS JIT CreatePixelShader failed (0x%08X)", hr);
        std::lock_guard<std::mutex> lock(g_PSJITMutex);
        g_PSJITCache[hash] = { nullptr };
        s_LastKey = key;
        s_LastPSRegsGen = psRegsGen;
        s_LastPSResult = nullptr;
        return nullptr;
    }

    // Save to disk cache for next session
    ShaderDiskCache::Save(hash, pCode);
    pCode->Release();

    EmuLog(LOG_LEVEL::DEBUG, "PS JIT: compiled new shader (hash=%016llX, stages=%u)", hash, key.numStages);

    // Cache
    std::lock_guard<std::mutex> lock(g_PSJITMutex);
    g_PSJITCache[hash] = { pPS };
    s_LastKey = key;
    s_LastPSRegsGen = psRegsGen;
    s_LastPSResult = pPS;
    return pPS;
}
