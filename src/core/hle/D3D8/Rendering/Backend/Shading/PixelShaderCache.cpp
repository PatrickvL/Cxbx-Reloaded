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
// Compile with D3DCompile (ps_5_0, O3) and cache by FNV-1a hash of topology.

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

// Externally stored last aux CB (built by CxbxD3D11UploadRCInterpreterState)
extern PSAuxCBLayout g_LastPSAuxCB;

// ============================================================
// State capture — everything that affects generated HLSL
// ============================================================
struct PSJITKey {
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
    uint64_t hash = 0xcbf29ce484222325ULL;
    const uint8_t* data = reinterpret_cast<const uint8_t*>(&key);
    for (size_t i = 0; i < sizeof(key); i++) {
        hash ^= data[i];
        hash *= 0x100000001b3ULL;
    }
    return hash;
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
// HLSL code generation helpers
// ============================================================

static const char* RegName(uint32_t regIdx)
{
    switch (regIdx) {
    case 0x00: return nullptr;  // ZERO — handled specially
    case 0x01: return "C0";
    case 0x02: return "C1";
    case 0x03: return "FOG";
    case 0x04: return "V0";
    case 0x05: return "V1";
    case 0x08: return "T0";
    case 0x09: return "T1";
    case 0x0a: return "T2";
    case 0x0b: return "T3";
    case 0x0c: return "R0";
    case 0x0d: return "R1";
    case 0x0e: return "V1R0_SUM";
    case 0x0f: return "EF_PROD";
    default:   return "R0"; // reserved → R0 fallback
    }
}

// Emit an RGB input expression with channel select + input mapping applied
static std::string EmitRGBInput(uint32_t inputByte)
{
    uint32_t regIdx  = inputByte & 0x0F;
    uint32_t mapping = inputByte & 0xE0;
    bool useAlpha    = (inputByte & 0x10) != 0;

    // Value expression
    std::string val;
    if (regIdx == 0) {
        // ZERO register
        val = "(float4)0.0";
    } else {
        val = RegName(regIdx);
        if (useAlpha)
            val += ".aaaa";
    }

    // Apply input mapping
    // For ZERO register, some mappings produce known constants
    if (regIdx == 0) {
        switch (mapping) {
        case 0x00: return "(float4)0.0";      // UNSIGNED_IDENTITY(0) = 0
        case 0x20: return "(float4)1.0";      // UNSIGNED_INVERT(0) = 1-0 = 1
        case 0x40: return "(float4)(-1.0)";   // EXPAND_NORMAL(0) = -1
        case 0x60: return "(float4)1.0";      // EXPAND_NEGATE(0) = 1
        case 0x80: return "(float4)(-0.5)";   // HALFBIAS_NORMAL(0) = -0.5
        case 0xa0: return "(float4)0.5";      // HALFBIAS_NEGATE(0) = 0.5
        case 0xc0: return "(float4)0.0";      // SIGNED_IDENTITY(0) = 0
        case 0xe0: return "(float4)0.0";      // SIGNED_NEGATE(0) = -0 = 0
        default:   return "(float4)0.0";
        }
    }

    switch (mapping) {
    case 0x00: return "max((float4)0.0, " + val + ")";                   // UNSIGNED_IDENTITY
    case 0x20: return "(1.0 - saturate(" + val + "))";                   // UNSIGNED_INVERT
    case 0x40: return "(2.0 * max((float4)0.0, " + val + ") - 1.0)";     // EXPAND_NORMAL
    case 0x60: return "(-(2.0 * max((float4)0.0, " + val + ") - 1.0))";  // EXPAND_NEGATE
    case 0x80: return "(max((float4)0.0, " + val + ") - 0.5)";           // HALFBIAS_NORMAL
    case 0xa0: return "(-(max((float4)0.0, " + val + ") - 0.5))";        // HALFBIAS_NEGATE
    case 0xc0: return val;                                               // SIGNED_IDENTITY
    case 0xe0: return "(-" + val + ")";                                  // SIGNED_NEGATE
    default:   return val;
    }
}

// Emit an alpha input expression (scalar)
static std::string EmitAlphaInput(uint32_t inputByte)
{
    uint32_t regIdx  = inputByte & 0x0F;
    uint32_t mapping = inputByte & 0xE0;
    bool useAlpha    = (inputByte & 0x10) != 0;

    // Value expression — scalar (.a or .b)
    std::string val;
    if (regIdx == 0) {
        val = "0.0";
    } else {
        val = std::string(RegName(regIdx)) + (useAlpha ? ".a" : ".b");
    }

    // ZERO register constant folding
    if (regIdx == 0) {
        switch (mapping) {
        case 0x00: return "0.0";
        case 0x20: return "1.0";
        case 0x40: return "(-1.0)";
        case 0x60: return "1.0";
        case 0x80: return "(-0.5)";
        case 0xa0: return "0.5";
        case 0xc0: return "0.0";
        case 0xe0: return "0.0";
        default:   return "0.0";
        }
    }

    switch (mapping) {
    case 0x00: return "max(0.0, " + val + ")";                   // UNSIGNED_IDENTITY
    case 0x20: return "(1.0 - saturate(" + val + "))";           // UNSIGNED_INVERT
    case 0x40: return "(2.0 * max(0.0, " + val + ") - 1.0)";     // EXPAND_NORMAL
    case 0x60: return "(-(2.0 * max(0.0, " + val + ") - 1.0))";  // EXPAND_NEGATE
    case 0x80: return "(max(0.0, " + val + ") - 0.5)";           // HALFBIAS_NORMAL
    case 0xa0: return "(-(max(0.0, " + val + ") - 0.5))";        // HALFBIAS_NEGATE
    case 0xc0: return val;                                       // SIGNED_IDENTITY
    case 0xe0: return "(-" + val + ")";                          // SIGNED_NEGATE
    default:   return val;
    }
}

// Emit a final combiner input (restricted mappings)
static std::string EmitFinalInput(uint32_t inputByte, bool isFinalABCD)
{
    uint32_t regIdx  = inputByte & 0x0F;
    uint32_t mapping = inputByte & 0xE0;
    bool useAlpha    = (inputByte & 0x10) != 0;

    // Restrict mapping for final combiner: anything >= 0x40 → just bit 5
    if (mapping >= 0x40)
        mapping = mapping & 0x20;

    std::string val;
    if (regIdx == 0) {
        // ZERO
        if (mapping == 0x20) return "(float4)1.0";
        return "(float4)0.0";
    }

    val = RegName(regIdx);

    // FOG special handling: rgb→0, alpha passthrough
    if (regIdx == 0x03) {
        val = "float4(0.0, 0.0, 0.0, FOG.a)";
    }

    // V1R0_SUM / EF_PROD only valid in ABCD phase
    if ((regIdx == 0x0e || regIdx == 0x0f) && !isFinalABCD) {
        val = "(float4)0.0";
    }

    if (useAlpha)
        val += ".aaaa";

    // Apply the restricted mapping (only UNSIGNED_IDENTITY or UNSIGNED_INVERT)
    if (mapping == 0x20)
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

    if (mode == 0) {
        // NONE — T register keeps VS texcoord (already initialized)
        return;
    }

    // Emit the sample
    switch (mode) {
    case 0x01: // PROJECT2D
        ss << "    { float3 proj = " << coords << ".xyz / " << coords << ".w;\n";
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", proj.xy);\n";
        if (shadowCompare != 0.0f) {
            ss << "      " << tReg << " = ApplyShadowCompare(" << stage << ", " << tReg << ", proj);\n";
        }
        ss << "    }\n";
        break;

    case 0x02: // PROJECT3D
        ss << "    { float3 proj = " << coords << ".xyz / " << coords << ".w;\n";
        if (shadowCompare != 0.0f) {
            ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", proj.xy);\n";
            ss << "      " << tReg << " = ApplyShadowCompare(" << stage << ", " << tReg << ", proj);\n";
        } else {
            ss << "      " << tReg << " = Tex3D_" << sIdx << ".Sample(Samp" << sIdx << ", proj);\n";
        }
        ss << "    }\n";
        break;

    case 0x03: // CUBEMAP
        ss << "    " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx << ", " << coords << ".xyz);\n";
        break;

    case 0x04: // PASSTHRU
        ss << "    " << tReg << " = saturate(" << coords << ");\n";
        break;

    case 0x05: // CLIPPLANE
        ss << "    ApplyCompareMode(" << stage << ", " << coords << ");\n";
        return; // no post-process

    case 0x06: // BUMPENVMAP
    {
        uint32_t bumOff = (stage - 1) * 4;
        ss << "    { float4 src = T" << ((shaderCtl >> (stage == 3 ? 20 : 16)) & (stage == 3 ? 3 : 1)) << ";\n";
        ss << "      float4 bem = float4(PG_FLOAT(0x" << std::hex << (0x181C + bumOff) << "), PG_FLOAT(0x" << (0x1828 + bumOff) << "), PG_FLOAT(0x" << (0x1834 + bumOff) << "), PG_FLOAT(0x" << (0x1840 + bumOff) << "));\n" << std::dec;
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", float2("
           << coords << ".x + bem.x * src.r + bem.z * src.g, "
           << coords << ".y + bem.y * src.r + bem.w * src.g));\n";
        ss << "    }\n";
        break;
    }

    case 0x07: // BUMPENVMAP_LUM
    {
        uint32_t bumOff = (stage - 1) * 4;
        ss << "    { float4 src = T" << ((shaderCtl >> (stage == 3 ? 20 : 16)) & (stage == 3 ? 3 : 1)) << ";\n";
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

    case 0x08: // BRDF
        ss << "    " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", " << coords << ".xy);\n";
        break;

    case 0x0F: // DPNDNT_AR
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", T" << srcStage << ".ar);\n";
        break;
    }

    case 0x10: // DPNDNT_GB
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx << ", T" << srcStage << ".gb);\n";
        break;
    }

    case 0x11: // DOTPRODUCT
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = float4(dot(" << coords << ".xyz, dm), 0.0, 0.0, 1.0);\n";
        ss << "    }\n";
        break;
    }

    case 0x09: // DOT_ST
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = Tex2D_" << sIdx << ".Sample(Samp" << sIdx
           << ", float2(T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "    }\n";
        break;
    }

    case 0x0A: // DOT_ZW
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      float d = dot(" << coords << ".xyz, dm);\n";
        ss << "      float depth = (abs(d) < 0.00001) ? 1.0 : (T" << (stage - 1) << ".x / d);\n";
        ss << "      " << tReg << " = depth.xxxx;\n";
        ss << "    }\n";
        break;
    }

    case 0x0B: // DOT_RFLCT_DIFF
    {
        // xemu reference: normal = (dot[stage-1], dot[stage], dot_next)
        // where dot_next peeks at the next stage's dot mapping and source.
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        uint32_t dotMapping = (shaderCtl >> ((stage - 1) * 4)) & 7;
        uint32_t nextStage = stage + 1;
        uint32_t nextSrcStage = (nextStage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1);
        uint32_t nextDotMapping = (shaderCtl >> ((nextStage - 1) * 4)) & 7;
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << dotMapping << ", src);\n";
        ss << "      float curDot = dot(" << coords << ".xyz, dm);\n";
        ss << "      float4 nextSrc = T" << nextSrcStage << ";\n";
        ss << "      float3 nextDm = ApplyDotMapping(" << nextDotMapping << ", nextSrc);\n";
        ss << "      float nextDot = dot(T" << nextStage << ".xyz, nextDm);\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(T" << (stage - 1) << ".x, curDot, nextDot));\n";
        ss << "    }\n";
        break;
    }

    case 0x0C: // DOT_RFLCT_SPEC
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      float3 N = normalize(float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "      float3 E = normalize(eyeVec);\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx << ", 2.0 * dot(N, E) * N - E);\n";
        ss << "    }\n";
        break;
    }

    case 0x0D: // DOT_STR_3D
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = Tex3D_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "    }\n";
        break;
    }

    case 0x0E: // DOT_STR_CUBE
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
        ss << "      float3 dm = ApplyDotMapping(" << ((shaderCtl >> ((stage - 1) * 4)) & 7) << ", src);\n";
        ss << "      " << tReg << " = TexCube_" << sIdx << ".Sample(Samp" << sIdx
           << ", float3(T" << (stage - 2) << ".x, T" << (stage - 1) << ".x, dot(" << coords << ".xyz, dm)));\n";
        ss << "    }\n";
        break;
    }

    case 0x12: // DOT_RFLCT_SPEC_CONST
    {
        uint32_t srcStage = (stage <= 1) ? 0 : ((stage == 3) ? ((shaderCtl >> 20) & 3) : ((shaderCtl >> 16) & 1));
        ss << "    { float4 src = T" << srcStage << ";\n";
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
    if (mode == 0x04) return; // PASSTHRU already saturated, no post-process
    if (mode == 0x05) return; // CLIPPLANE has no texel

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
        ss << "    " << tReg << " = PerformColorKeyOp((int)ColorKeyOp[" << stage << "].x, ColorKeyColor[" << stage << "], " << tReg << ");\n";
    }

    // Alpha kill
    if (alphaKill != 0.0f) {
        ss << "    PerformAlphaKill(1, " << tReg << ");\n";
    }
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
    texModes[0] = (key.textureModes      ) & 0x1F;
    texModes[1] = (key.textureModes >>  5) & 0x1F;
    texModes[2] = (key.textureModes >> 10) & 0x1F;
    texModes[3] = (key.textureModes >> 15) & 0x1F;

    // ---- Header / resource declarations ----
    ss << "// Auto-generated by CxbxPixelShaderJIT\n";
    ss << "// Stages: " << numStages << " TexModes: "
       << texModes[0] << "/" << texModes[1] << "/" << texModes[2] << "/" << texModes[3] << "\n\n";

    // StructuredBuffer for PGRAPH regs
    ss << "StructuredBuffer<uint> g_PGRegs : register(t12);\n";
    ss << "uint PG_UINT(uint byteOff) { return g_PGRegs[byteOff >> 2]; }\n";
    ss << "float PG_FLOAT(uint byteOff) { return asfloat(g_PGRegs[byteOff >> 2]); }\n";
    ss << "float4 UnpackABGR(uint c) {\n";
    ss << "    return float4(float(c & 0xFFu)/255.0, float((c>>8)&0xFFu)/255.0,\n";
    ss << "                  float((c>>16)&0xFFu)/255.0, float((c>>24)&0xFFu)/255.0);\n";
    ss << "}\n";
    ss << "float4 UnpackARGB(uint c) {\n";
    ss << "    return float4(float((c>>16)&0xFFu)/255.0, float((c>>8)&0xFFu)/255.0,\n";
    ss << "                  float(c&0xFFu)/255.0, float((c>>24)&0xFFu)/255.0);\n";
    ss << "}\n";
    ss << "float4 PG_COLOR(uint byteOff) { return UnpackABGR(PG_UINT(byteOff)); }\n";
    ss << "float4 PG_COLOR_ARGB(uint byteOff) { return UnpackARGB(PG_UINT(byteOff)); }\n\n";

    // Texture/sampler declarations — only what's needed
    for (uint32_t i = 0; i < 4; i++) {
        bool need2D = false, need3D = false, needCube = false;
        switch (texModes[i]) {
        case 0x01: case 0x06: case 0x07: case 0x08: case 0x09:
        case 0x0F: case 0x10:
            need2D = true; break;
        case 0x02: case 0x0D:
            need3D = true; break;
        case 0x03: case 0x0B: case 0x0C: case 0x0E: case 0x12:
            needCube = true; break;
        }
        // PROJECT3D with shadow uses 2D
        if (texModes[i] == 0x02 && key.shadowCompare[i] != 0.0f)
            need2D = true, need3D = false;

        if (need2D)   ss << "Texture2D<float4>   Tex2D_" << i << "   : register(t" << i << ");\n";
        if (need3D)   ss << "Texture3D<float4>   Tex3D_" << i << "   : register(t" << (i+4) << ");\n";
        if (needCube) ss << "TextureCube<float4> TexCube_" << i << " : register(t" << (i+8) << ");\n";
        if (texModes[i] != 0)
            ss << "SamplerState Samp" << i << " : register(s" << i << ");\n";
    }
    ss << "\n";

    // Aux cbuffer
    ss << "cbuffer PSAuxCBLayout : register(b0) {\n";
    ss << "    uint PSTextureModes; uint3 _p0;\n";
    ss << "    uint PSFinalCombinerInputsABCD; uint3 _p1;\n";
    ss << "    uint PSFinalCombinerInputsEFG; uint3 _p2;\n";
    ss << "    float4 ColorSign[4];\n";
    ss << "    float4 TexFmtFixup;\n";
    ss << "    float4 ColorKeyOp[4];\n";
    ss << "    float4 ColorKeyColor[4];\n";
    ss << "    float4 AlphaKill;\n";
    ss << "    float4 FogInfo;\n";
    ss << "    uint FogEnable; uint3 _p3;\n";
    ss << "    float4 FrontFaceInfo;\n";
    ss << "    float4 ShadowCompare;\n";
    ss << "};\n\n";

    // PS_INPUT — shared with VS output and RC interpreter
    ss << "#include \"CxbxPixelShaderInput.hlsli\"\n\n";

    // Shared helper functions (nv2a_mul, PerformColorSign, etc.)
    ss << "#include \"CxbxNV2AMathHelpers.hlsli\"\n";
    ss << "#include \"CxbxPixelShaderFunctions.hlsli\"\n\n";

    // Shadow compare helper
    bool anyShadow = (key.shadowCompare[0] != 0.0f || key.shadowCompare[1] != 0.0f ||
                      key.shadowCompare[2] != 0.0f || key.shadowCompare[3] != 0.0f);
    if (anyShadow) {
        ss << "float4 ApplyShadowCompare(uint stage, float4 sampled, float3 coords) {\n";
        ss << "    uint sf = PG_UINT(0x19A4) & 7;\n";
        ss << "    if (sf==0) return sampled;\n";
        ss << "    if (sf==7) return 1.0.xxxx;\n";
        ss << "    float depth = sampled.r, ref = coords.z;\n";
        ss << "    uint cmp = (depth<ref?1u:0u)|(depth==ref?2u:0u)|(depth>ref?4u:0u);\n";
        ss << "    return ((cmp&sf)!=0u ? 1.0 : 0.0).xxxx;\n";
        ss << "}\n\n";
    }

    // Clip plane helper
    bool anyClip = false;
    for (int i = 0; i < 4; i++) if (texModes[i] == 0x05) anyClip = true;
    if (anyClip) {
        ss << "void ApplyCompareMode(uint stage, float4 coords) {\n";
        ss << "    uint bits = (PG_UINT(0x1994) >> (stage*4)) & 0xF;\n";
        ss << "    bool k0 = (bits&1)!=0 ? (coords.x>=0.0) : (coords.x<0.0);\n";
        ss << "    bool k1 = (bits&2)!=0 ? (coords.y>=0.0) : (coords.y<0.0);\n";
        ss << "    bool k2 = (bits&4)!=0 ? (coords.z>=0.0) : (coords.z<0.0);\n";
        ss << "    bool k3 = (bits&8)!=0 ? (coords.w>=0.0) : (coords.w<0.0);\n";
        ss << "    if (k0||k1||k2||k3) discard;\n";
        ss << "}\n\n";
    }

    // ApplyDotMapping is provided by CxbxPixelShaderFunctions.hlsli (included above)

    // ---- main() ----
    ss << "float4 main(PS_INPUT input) : SV_Target\n{\n";

    // Declare register variables
    ss << "    float4 T0 = input.iT0;\n";
    ss << "    float4 T1 = input.iT1;\n";
    ss << "    float4 T2 = input.iT2;\n";
    ss << "    float4 T3 = input.iT3;\n";

    // Eye vector for DOT_RFLCT_SPEC
    bool needsEyeVec = false;
    for (int i = 0; i < 4; i++)
        if (texModes[i] == 0x0C) needsEyeVec = true;
    if (needsEyeVec)
        ss << "    float3 eyeVec = float3(input.iT1.w, input.iT2.w, input.iT3.w);\n";

    // Texture fetches (sequential — later stages can depend on earlier)
    ss << "\n    // --- Texture fetches ---\n";
    for (uint32_t i = 0; i < 4; i++) {
        if (texModes[i] != 0) {
            float cs4[4] = { key.colorSign[i*4+0], key.colorSign[i*4+1],
                             key.colorSign[i*4+2], key.colorSign[i*4+3] };
            EmitTextureFetch(ss, i, texModes[i], key.shaderCtl,
                             key.shadowCompare[i], cs4, key.texFmtFixup[i],
                             key.colorKeyOp[i], key.alphaKill[i]);
        }
    }

    // Vertex colors with front/back face selection
    ss << "\n    // --- Vertex colors ---\n";
    if (key.frontFaceFactor != 0.0f) {
        ss << "    float faceSign = input.iFF ? 1.0 : -1.0;\n";
        ss << "    bool isFront = (faceSign * FrontFaceInfo.x) >= 0.0;\n";
        ss << "    float4 V0 = isFront ? input.iD0 : input.iB0;\n";
        ss << "    float4 V1 = isFront ? input.iD1 : input.iB1;\n";
    } else {
        ss << "    float4 V0 = input.iD0;\n";
        ss << "    float4 V1 = input.iD1;\n";
    }
    ss << "    float4 FOG = float4(PG_COLOR_ARGB(0x1980).rgb, saturate(input.iFog));\n";
    ss << "    float4 R0 = float4(0.0, 0.0, 0.0, T0.a);\n";
    ss << "    float4 R1 = (float4)0.0;\n";
    ss << "    float4 V1R0_SUM = (float4)0.0;\n";
    ss << "    float4 EF_PROD = (float4)0.0;\n";
    ss << "    float4 C0, C1;\n";

    // Combiner stages
    ss << "\n    // --- Combiner stages ---\n";
    for (uint32_t stage = 0; stage < numStages; stage++) {
        ss << "    // Stage " << stage << "\n";

        // Load C0/C1 from PGRAPH
        uint32_t c0Off = flagUniqueC0 ? (0x1880 + stage * 4) : 0x1880;
        uint32_t c1Off = flagUniqueC1 ? (0x18A0 + stage * 4) : 0x18A0;
        ss << "    C0 = PG_COLOR(0x" << std::hex << c0Off << ");\n";
        ss << "    C1 = PG_COLOR(0x" << c1Off << ");\n" << std::dec;

        // Decode output control
        uint32_t rgbOut = key.rgbOutputs[stage];
        uint32_t aOut   = key.alphaOutputs[stage];
        uint32_t rgbFlags = rgbOut >> 12;
        uint32_t aFlags   = aOut >> 12;
        bool flagCDDot  = (rgbFlags & 0x01) != 0;
        bool flagABDot  = (rgbFlags & 0x02) != 0;
        bool flagRGBMux = (rgbFlags & 0x04) != 0;
        bool flagAMux   = (aFlags   & 0x04) != 0;
        bool cdBlue2A   = (rgbFlags & 0x40) != 0;
        bool abBlue2A   = (rgbFlags & 0x80) != 0;

        // Output destinations
        uint32_t rgbRegAB  = (rgbOut >> 4) & 0xF;
        uint32_t rgbRegCD  = (rgbOut     ) & 0xF;
        uint32_t rgbRegSum = (rgbOut >> 8) & 0xF;
        uint32_t aRegAB    = (aOut   >> 4) & 0xF;
        uint32_t aRegCD    = (aOut       ) & 0xF;
        uint32_t aRegSum   = (aOut   >> 8) & 0xF;

        // Decode inputs
        uint32_t rgbIn = key.rgbInputs[stage];
        uint32_t aIn   = key.alphaInputs[stage];

        std::string prefix = "s" + std::to_string(stage) + "_";

        // RGB inputs
        std::string rgbA = EmitRGBInput((rgbIn >> 24) & 0xFF);
        std::string rgbB = EmitRGBInput((rgbIn >> 16) & 0xFF);
        std::string rgbC = EmitRGBInput((rgbIn >>  8) & 0xFF);
        std::string rgbD = EmitRGBInput((rgbIn      ) & 0xFF);

        // Alpha inputs
        std::string aA = EmitAlphaInput((aIn >> 24) & 0xFF);
        std::string aB = EmitAlphaInput((aIn >> 16) & 0xFF);
        std::string aC = EmitAlphaInput((aIn >>  8) & 0xFF);
        std::string aD = EmitAlphaInput((aIn      ) & 0xFF);

        // Compute AB, CD
        ss << "    {\n";
        ss << "        float4 rgbA = " << rgbA << ";\n";
        ss << "        float4 rgbB = " << rgbB << ";\n";
        ss << "        float4 rgbC = " << rgbC << ";\n";
        ss << "        float4 rgbD = " << rgbD << ";\n";
        ss << "        float aA = " << aA << ";\n";
        ss << "        float aB = " << aB << ";\n";
        ss << "        float aC = " << aC << ";\n";
        ss << "        float aD = " << aD << ";\n";

        // RGB products
        if (flagABDot)
            ss << "        float3 rgbAB = (float3)dot(rgbA.rgb, rgbB.rgb);\n";
        else
            ss << "        float3 rgbAB = nv2a_mul3(rgbA.rgb, rgbB.rgb);\n";

        if (flagCDDot)
            ss << "        float3 rgbCD = (float3)dot(rgbC.rgb, rgbD.rgb);\n";
        else
            ss << "        float3 rgbCD = nv2a_mul3(rgbC.rgb, rgbD.rgb);\n";

        // Alpha products
        ss << "        float aAB = nv2a_mul1(aA, aB);\n";
        ss << "        float aCD = nv2a_mul1(aC, aD);\n";

        // SUM/MUX
        bool writeSumRGB = !flagABDot && !flagCDDot;
        if (flagRGBMux || flagAMux) {
            if (flagMuxMsb)
                ss << "        bool muxSel = (R0.a >= 0.5);\n";
            else
                ss << "        bool muxSel = (((uint)(saturate(R0.a)*255.0+0.5) & 1u) != 0u);\n";
        }

        if (writeSumRGB) {
            if (flagRGBMux)
                ss << "        float3 rgbSUM = muxSel ? rgbCD : rgbAB;\n";
            else
                ss << "        float3 rgbSUM = rgbAB + rgbCD;\n";
        }

        if (flagAMux)
            ss << "        float aSUM = muxSel ? aCD : aAB;\n";
        else
            ss << "        float aSUM = aAB + aCD;\n";

        // Output mapping
        bool hasBias = (rgbFlags & 0x08) != 0;
        uint32_t scaleMode = (rgbFlags >> 4) & 3;
        float rgbBias = hasBias ? -0.5f : 0.0f;
        float rgbScale = 1.0f;
        switch (scaleMode) {
        case 0: rgbScale = 1.0f; break;
        case 1: rgbScale = 2.0f; break;
        case 2: rgbScale = 4.0f; break;
        case 3: rgbScale = 0.5f; break;
        }

        bool aBiasBit = (aFlags & 0x08) != 0;
        uint32_t aScaleMode = (aFlags >> 4) & 3;
        float aBias = aBiasBit ? -0.5f : 0.0f;
        float aScale = 1.0f;
        switch (aScaleMode) {
        case 0: aScale = 1.0f; break;
        case 1: aScale = 2.0f; break;
        case 2: aScale = 4.0f; break;
        case 3: aScale = 0.5f; break;
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

        std::string outRGB_AB = emitClamp("rgbAB", rgbBias, rgbScale);
        std::string outRGB_CD = emitClamp("rgbCD", rgbBias, rgbScale);
        std::string outRGB_Sum = writeSumRGB ? emitClamp("rgbSUM", rgbBias, rgbScale) : "(float3)0.0";
        std::string outA_AB  = emitClamp("aAB", aBias, aScale);
        std::string outA_CD  = emitClamp("aCD", aBias, aScale);
        std::string outA_Sum = emitClamp("aSUM", aBias, aScale);

        // RGB writes
        auto emitRGBWrite = [&](uint32_t dest, const std::string& rgbExpr,
                                bool blue2Alpha, uint32_t aRegSame) {
            if (dest == 0) return; // DISCARD
            const char* dName = RegName(dest);
            if (!dName) return;
            bool aOverwrite = (aRegSame == dest) && (aRegSame != 0);
            if (blue2Alpha && !aOverwrite) {
                ss << "        { float3 _rgb = " << rgbExpr << ";\n";
                ss << "          " << dName << " = float4(_rgb, _rgb.b); }\n";
            } else {
                ss << "        " << dName << ".rgb = " << rgbExpr << ";\n";
            }
        };

        emitRGBWrite(rgbRegAB, outRGB_AB, abBlue2A, aRegAB);
        emitRGBWrite(rgbRegCD, outRGB_CD, cdBlue2A, aRegCD);
        if (writeSumRGB && rgbRegSum != 0) {
            const char* dName = RegName(rgbRegSum);
            if (dName)
                ss << "        " << dName << ".rgb = " << outRGB_Sum << ";\n";
        }

        // Alpha writes
        auto emitAlphaWrite = [&](uint32_t dest, const std::string& aExpr) {
            if (dest == 0) return; // DISCARD
            const char* dName = RegName(dest);
            if (!dName) return;
            ss << "        " << dName << ".a = " << aExpr << ";\n";
        };

        emitAlphaWrite(aRegAB, outA_AB);
        emitAlphaWrite(aRegCD, outA_CD);
        emitAlphaWrite(aRegSum, outA_Sum);

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
        ss << "    C0 = PG_COLOR(0x19AC);\n";
        ss << "    C1 = PG_COLOR(0x19B0);\n";

        // EFG phase
        uint32_t settings = (fcEFG) & 0xFF;
        uint32_t eReg = (fcEFG >> 24) & 0xFF;
        uint32_t fReg = (fcEFG >> 16) & 0xFF;
        uint32_t gReg = (fcEFG >>  8) & 0xFF;

        ss << "    float3 fcE = (" << EmitFinalInput(eReg, false) << ").rgb;\n";
        ss << "    float3 fcF = (" << EmitFinalInput(fReg, false) << ").rgb;\n";
        ss << "    float  fcG = (" << EmitFinalInput(gReg, false) << ").a;\n";

        ss << "    EF_PROD = float4(fcE * fcF, 1.0);\n";

        // V1R0_SUM
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

        // ABCD phase
        uint32_t aReg = (fcABCD >> 24) & 0xFF;
        uint32_t bReg = (fcABCD >> 16) & 0xFF;
        uint32_t cReg = (fcABCD >>  8) & 0xFF;
        uint32_t dReg = (fcABCD      ) & 0xFF;

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

    ss << "\n    return result;\n}\n";

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

    // Fast path: if PGRAPH registers haven't changed since last call,
    // the combiner topology is identical — return cached result directly
    // without rebuilding the key or hashing.
    static uint32_t s_LastPSRegsGen = ~0u;
    static ID3D11PixelShader* s_LastPSResult = nullptr;
    static PSJITKey s_LastKey = {};
    if (pg->regs_generation == s_LastPSRegsGen)
        return s_LastPSResult;

    // Capture current state — read directly from PGRAPH registers
    // (no dependency on g_LastPSAuxCB or CxbxD3D11UploadRCInterpreterState)
    PSJITKey key = {};
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

    key.textureModes = aux.PSTextureModes.value;
    key.fcABCD       = aux.PSFinalCombinerInputsABCD.value;
    key.fcEFG        = aux.PSFinalCombinerInputsEFG.value;
    key.fogEnable    = aux.FogEnable.value;
    key.frontFaceFactor = aux.FrontFaceInfo.x;

    for (int i = 0; i < 4; i++) {
        key.colorKeyOp[i]   = aux.ColorKeyOp[i].x;
        key.shadowCompare[i]= (&aux.ShadowCompare.x)[i];
        key.texFmtFixup[i]  = (&aux.TexFmtFixup.x)[i];
        key.colorSign[i*4+0]= aux.ColorSign[i].x;
        key.colorSign[i*4+1]= aux.ColorSign[i].y;
        key.colorSign[i*4+2]= aux.ColorSign[i].z;
        key.colorSign[i*4+3]= aux.ColorSign[i].w;
    }
    key.alphaKill[0] = aux.AlphaKill.x;
    key.alphaKill[1] = aux.AlphaKill.y;
    key.alphaKill[2] = aux.AlphaKill.z;
    key.alphaKill[3] = aux.AlphaKill.w;

    // Second fast path: if the key matches the last one (combiner state unchanged
    // despite regs_generation bumping from non-combiner register writes like VS
    // constants), skip the expensive hash + mutex + map lookup.
    if (memcmp(&key, &s_LastKey, sizeof(PSJITKey)) == 0) {
        s_LastPSRegsGen = pg->regs_generation;
        return s_LastPSResult;
    }

    // Hash and cache lookup
    uint64_t hash = HashKey(key);
    {
        std::lock_guard<std::mutex> lock(g_PSJITMutex);
        auto it = g_PSJITCache.find(hash);
        if (it != g_PSJITCache.end()) {
            s_LastKey = key;
            s_LastPSRegsGen = pg->regs_generation;
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
            s_LastPSRegsGen = pg->regs_generation;
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
        s_LastPSRegsGen = pg->regs_generation;
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
        s_LastPSRegsGen = pg->regs_generation;
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
    s_LastPSRegsGen = pg->regs_generation;
    s_LastPSResult = pPS;
    return pPS;
}
