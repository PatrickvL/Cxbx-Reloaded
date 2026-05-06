// VertexShaderCache.cpp — Runtime NV2A→HLSL JIT compiler + cache
//
// Translates NV2A vertex transform microcode into straight-line HLSL,
// compiles with D3DCompile, caches by FNV-1a hash of program tokens.
// A typical 20-instruction NV2A program becomes ~30 lines of HLSL
// with zero loops — orders of magnitude faster than interpreting.

#define LOG_PREFIX CXBXR_MODULE::VTXSH

#include "VertexShaderCache.h"
#include "ShaderDiskCache.h"
#include "../Backend_D3D11_Profiler.h"
#include "core/kernel/init/CxbxKrnl.h" // EmuLog
#include "common/Logging.h"

#include <unordered_map>
#include <string>
#include <sstream>
#include <mutex>
#include <cstring>
#include <d3dcompiler.h>

// ---------------------------------------------------------------------------
// CxbxJITIncludeHandler — ID3DInclude that resolves #include directives by
// loading files relative to the emulator executable. The JIT-compiled HLSL
// references shared headers (CxbxVertexShaderCommon.hlsli, etc.) that ship
// next to the exe in a "hlsl\" sub-directory (copied there by CMake POST_BUILD).
// ---------------------------------------------------------------------------
class CxbxJITIncludeHandler : public ID3DInclude
{
public:
    CxbxJITIncludeHandler()
    {
        char exePath[MAX_PATH] = {};
        GetModuleFileNameA(nullptr, exePath, MAX_PATH);
        char* lastSlash = strrchr(exePath, '\\');
        if (lastSlash) *(lastSlash + 1) = '\0';
        m_shaderDir = std::string(exePath) + "hlsl\\";
    }

    STDMETHOD(Open)(D3D_INCLUDE_TYPE /*IncludeType*/, LPCSTR pFileName,
        LPCVOID /*pParentData*/, LPCVOID* ppData, UINT* pBytes) override
    {
        std::string fullPath = m_shaderDir + pFileName;
        HANDLE hFile = CreateFileA(fullPath.c_str(), GENERIC_READ, FILE_SHARE_READ,
            nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (hFile == INVALID_HANDLE_VALUE)
            return E_FAIL;

        DWORD fileSize = GetFileSize(hFile, nullptr);
        if (fileSize == INVALID_FILE_SIZE) { CloseHandle(hFile); return E_FAIL; }

        char* pData = new (std::nothrow) char[fileSize];
        if (!pData) { CloseHandle(hFile); return E_OUTOFMEMORY; }

        DWORD bytesRead = 0;
        BOOL ok = ReadFile(hFile, pData, fileSize, &bytesRead, nullptr);
        CloseHandle(hFile);
        if (!ok || bytesRead != fileSize) { delete[] pData; return E_FAIL; }

        *ppData = pData;
        *pBytes = bytesRead;
        return S_OK;
    }

    STDMETHOD(Close)(LPCVOID pData) override
    {
        delete[] static_cast<const char*>(pData);
        return S_OK;
    }

private:
    std::string m_shaderDir;
};

// nv2a_vsh_cpu disassembler for instruction decoding
extern "C" {
#include "nv2a_vsh_disassembler.h"
}

// ============================================================
// Simple hash for program cache keying
// ============================================================
static uint64_t HashProgram(const uint32_t program_data[][4], uint32_t startAddr, uint32_t count)
{
    // FNV-1a 64-bit
    uint64_t hash = 0xcbf29ce484222325ULL;
    const uint8_t* data = reinterpret_cast<const uint8_t*>(&program_data[startAddr][0]);
    size_t len = count * 4 * sizeof(uint32_t);
    for (size_t i = 0; i < len; i++) {
        hash ^= data[i];
        hash *= 0x100000001b3ULL;
    }
    return hash;
}

// ============================================================
// Cache
// ============================================================
struct CachedShader {
    ID3D11VertexShader* pVS;
    ID3DBlob* pBytecode;
};

static std::unordered_map<uint64_t, CachedShader> g_JITCache;
static std::mutex g_JITMutex;

VertexShaderCache g_VertexShaderCache;

void VertexShaderCache::Clear()
{
    std::lock_guard<std::mutex> lock(g_JITMutex);
    for (auto& pair : g_JITCache) {
        if (pair.second.pVS) pair.second.pVS->Release();
        if (pair.second.pBytecode) pair.second.pBytecode->Release();
    }
    g_JITCache.clear();
}

// ============================================================
// HLSL code generation helpers
// ============================================================
static const char* kSwizzleChars = "xyzw";

static std::string EmitSwizzle(const uint8_t swz[4])
{
    // If identity swizzle (.xyzw), emit nothing
    if (swz[0] == NV2ASW_X && swz[1] == NV2ASW_Y && swz[2] == NV2ASW_Z && swz[3] == NV2ASW_W)
        return "";
    std::string s = ".";
    s += kSwizzleChars[swz[0]];
    s += kSwizzleChars[swz[1]];
    s += kSwizzleChars[swz[2]];
    s += kSwizzleChars[swz[3]];
    return s;
}

static std::string EmitWritemask(Nv2aVshWritemask wm)
{
    if (wm == NV2AWM_XYZW) return "";
    std::string s = ".";
    if (wm & NV2AWM_X) s += 'x'; // = 8
    if (wm & NV2AWM_Y) s += 'y'; // = 4
    if (wm & NV2AWM_Z) s += 'z'; // = 2
    if (wm & NV2AWM_W) s += 'w'; // = 1
    return s;
}

static std::string EmitInputReg(const Nv2aVshInput& input)
{
    std::string reg;
    switch (input.type) {
    case NV2ART_TEMPORARY:
        if (input.index == 12) reg = "oPos"; // r12 aliases oPos
        else reg = "r" + std::to_string(input.index);
        break;
    case NV2ART_INPUT:
        reg = "v[" + std::to_string(input.index) + "]";
        break;
    case NV2ART_CONTEXT:
        if (input.is_relative)
            reg = "C[" + std::to_string(input.index) + " + a0]";
        else
            reg = "C[" + std::to_string(input.index) + "]";
        break;
    default:
        reg = "float4(0,0,0,0)";
        break;
    }

    std::string swz = EmitSwizzle(input.swizzle);
    std::string expr = reg + swz;
    if (input.is_negated)
        expr = "(-" + expr + ")";
    return expr;
}

static std::string EmitOutputReg(const Nv2aVshOutput& output)
{
    switch (output.type) {
    case NV2ART_TEMPORARY:
        if (output.index == 12) return "oPos";
        return "r" + std::to_string(output.index);
    case NV2ART_OUTPUT:
        switch (output.index) {
        case 0: return "oPos";
        case 3: return "oD0";
        case 4: return "oD1";
        case 5: return "oFog";
        case 6: return "oPts";
        case 7: return "oB0";
        case 8: return "oB1";
        case 9: return "oT0";
        case 10: return "oT1";
        case 11: return "oT2";
        case 12: return "oT3";
        default: return "oPos"; // shouldn't happen
        }
    case NV2ART_CONTEXT:
        // Context writes — emit to a temp and handle specially
        return "ctx_" + std::to_string(output.index);
    case NV2ART_ADDRESS:
        return "a0_tmp";
    default:
        return "discard_tmp";
    }
}

static void EmitOperation(std::ostringstream& ss, const Nv2aVshOperation& op, bool isILU)
{
    if (op.opcode == NV2AOP_NOP) return;

    // Get output destination
    std::string dest;
    std::string wmask;
    bool hasOutput = false;
    for (int i = 0; i < 2; i++) {
        if (op.outputs[i].type != NV2ART_NONE) {
            dest = EmitOutputReg(op.outputs[i]);
            wmask = EmitWritemask(op.outputs[i].writemask);
            hasOutput = true;

            // Handle ARL specially
            if (op.opcode == NV2AOP_ARL) {
                std::string src = EmitInputReg(op.inputs[0]);
                ss << "    a0 = (int)floor(" << src << ".x + 0.001);\n";
                return;
            }

            // Get input expressions
            std::string a = EmitInputReg(op.inputs[0]);
            std::string b = (op.inputs[1].type != NV2ART_NONE) ? EmitInputReg(op.inputs[1]) : "float4(0,0,0,0)";
            std::string c = (op.inputs[2].type != NV2ART_NONE) ? EmitInputReg(op.inputs[2]) : "float4(0,0,0,0)";

            std::string expr;
            switch (op.opcode) {
            case NV2AOP_MOV: expr = a; break;
            case NV2AOP_MUL: expr = "nv2a_mul(" + a + ", " + b + ")"; break;
            case NV2AOP_ADD: expr = "(" + a + " + " + b + ")"; break;
            case NV2AOP_MAD: expr = "(nv2a_mul(" + a + ", " + b + ") + " + c + ")"; break;
            case NV2AOP_DP3: expr = "nv2a_dot3(" + a + ".xyz, " + b + ".xyz).xxxx"; break;
            case NV2AOP_DPH: expr = "(nv2a_dot3(" + a + ".xyz, " + b + ".xyz) + " + b + ".w).xxxx"; break;
            case NV2AOP_DP4: expr = "nv2a_dot4(" + a + ", " + b + ").xxxx"; break;
            case NV2AOP_DST: expr = "float4(1.0, nv2a_mul1(" + a + ".y, " + b + ".y), " + a + ".z, " + b + ".w)"; break;
            case NV2AOP_MIN: expr = "min(" + a + ", " + b + ")"; break;
            case NV2AOP_MAX: expr = "max(" + a + ", " + b + ")"; break;
            case NV2AOP_SLT: expr = "float4(" + a + " < " + b + ")"; break;
            case NV2AOP_SGE: expr = "float4(" + a + " >= " + b + ")"; break;
            case NV2AOP_RCP: expr = "(1.0 / " + a + ".x).xxxx"; break;
            case NV2AOP_RCC: {
                // Sign-preserving clamp: clamp(|1/x|) then restore original sign bit
                expr = "asfloat(asuint(clamp(abs(1.0 / " + a + ".x), 5.42101e-20, 1.84467e+19)) | (asuint(1.0 / " + a + ".x) & 0x80000000u)).xxxx";
                break;
            }
            case NV2AOP_RSQ: expr = "rsqrt(abs(" + a + ".x)).xxxx"; break;
            case NV2AOP_EXP: expr = "float4(exp2(floor(" + a + ".x)), " + a + ".x - floor(" + a + ".x), exp2(" + a + ".x), 1.0)"; break;
            case NV2AOP_LOG: {
                // LOG(0) = (-inf, 1, -inf, 1) per NV2A hardware
                std::string t = "abs(" + a + ".x)";
                expr = "(" + t + " == 0.0) ? float4(asfloat(0xFF800000u), 1.0, asfloat(0xFF800000u), 1.0) : float4(floor(log2(" + t + ")), " + t + " * exp2(-floor(log2(" + t + "))), log2(" + t + "), 1.0)";
                break;
            }
            case NV2AOP_LIT: {
                // LIT: result = {1.0, max(src.x, 0), (src.x>0 && src.y>0) ? pow(src.y, src.w) : 0, 1.0}
                // Simplified for common case
                expr = "float4(1.0, max(" + a + ".x, 0.0), (" + a + ".x > 0 && " + a + ".y > 0) ? exp2(clamp(" + a + ".w, -128, 128) * log2(max(" + a + ".y, 0.00001))) : 0.0, 1.0)";
                break;
            }
            default: expr = "float4(0,0,0,0)"; break;
            }

            // Emit assignment with writemask
            if (wmask.empty()) {
                ss << "    " << dest << " = " << expr << ";\n";
            } else {
                ss << "    " << dest << wmask << " = (" << expr << ")" << wmask << ";\n";
            }

            // Handle context writes (write-back to C[] array)
            if (op.outputs[i].type == NV2ART_CONTEXT) {
                // We emitted to a temp; now handle the context register write
                // This is rare — most programs don't use it
            }
        }
    }
}

// ============================================================
// Main translation function: NV2A program → HLSL source
// ============================================================
static std::string TranslateToHLSL(const uint32_t program_data[][4], uint32_t startAddr, uint32_t instrCount)
{
    std::ostringstream ss;

    // Preamble — include all shared headers
    ss << "// JIT-compiled NV2A vertex shader (" << instrCount << " instructions)\n";
    ss << "#include \"CxbxVertexShaderCommon.hlsli\"\n";
    ss << "#include \"CxbxVertexFetch.hlsli\"\n";
    ss << "#include \"CxbxScreenspaceTransform.hlsli\"\n";
    ss << "#include \"CxbxNV2AMathHelpers.hlsli\"\n\n";
    ss << "#define X_D3DVS_CONSTREG_COUNT 192\n";
    ss << "uniform float4 C[X_D3DVS_CONSTREG_COUNT] : register(c0);\n\n";
    ss << "VS_OUTPUT main(const VS_INPUT xIn)\n{\n";
    ss << "    float4 v[16];\n";
    ss << "    FetchAllAttributes(ResolveVertexIndex(xIn.vertexId), v);\n\n";

    // Declare temp registers — r0-r11 start at zero (NV2A hardware)
    ss << "    float4 r0=0, r1=0, r2=0, r3=0, r4=0, r5=0;\n";
    ss << "    float4 r6=0, r7=0, r8=0, r9=0, r10=0, r11=0;\n";
    // Output registers — must match NV2A hardware defaults:
    //   oPos/oT0-oT3: W=1 (homogeneous), oD0/oD1/oB0/oB1/oFog: all 1 (white/full)
    ss << "    float4 oPos = float4(0,0,0,1);\n";
    ss << "    float4 oD0 = float4(1,1,1,1), oD1 = float4(1,1,1,1);\n";
    ss << "    float4 oFog = float4(1,1,1,1), oPts = float4(0,0,0,0);\n";
    ss << "    float4 oB0 = float4(1,1,1,1), oB1 = float4(1,1,1,1);\n";
    ss << "    float4 oT0 = float4(0,0,0,1), oT1 = float4(0,0,0,1);\n";
    ss << "    float4 oT2 = float4(0,0,0,1), oT3 = float4(0,0,0,1);\n";
    ss << "    int a0 = 0;\n\n";

    // Translate each instruction
    for (uint32_t i = 0; i < instrCount; i++) {
        Nv2aVshStep step = {};
        Nv2aVshParseResult result = nv2a_vsh_parse_step(&step, program_data[startAddr + i]);
        if (result != NV2AVPR_SUCCESS) {
            return ""; // Can't translate — bail out
        }

        // For paired instructions, check if ILU reads a register that MAC writes.
        // On NV2A both execute in parallel (read-before-write), but in HLSL they're sequential.
        // Save the pre-MAC value of any register ILU reads that MAC also writes to.
        bool needSave = false;
        std::string saveReg;
        if (step.mac.opcode != NV2AOP_NOP && step.ilu.opcode != NV2AOP_NOP) {
            // Determine MAC write destination (temp register)
            int macWriteIdx = -1;
            if (step.mac.outputs[0].type == NV2ART_TEMPORARY)
                macWriteIdx = (int)step.mac.outputs[0].index;
            
            // Check if ILU reads from that same temp register
            if (macWriteIdx >= 0 && step.ilu.inputs[0].type == NV2ART_TEMPORARY
                && (int)step.ilu.inputs[0].index == macWriteIdx) {
                needSave = true;
                if (macWriteIdx == 12) saveReg = "oPos";
                else saveReg = "r" + std::to_string(macWriteIdx);
            }
        }

        ss << "    // Instruction " << i << "\n";
        if (needSave) {
            ss << "    { float4 _save = " << saveReg << ";\n";
            EmitOperation(ss, step.mac, false);
            // Patch: replace the ILU input source with the saved value
            // For simplicity, emit the ILU using _save directly
            // Override ILU input by emitting with temporary rename
            ss << "    " << saveReg << " = _save; // restore for ILU read\n";
            EmitOperation(ss, step.ilu, true);
            ss << "    }\n";
        } else {
            EmitOperation(ss, step.mac, false);
            EmitOperation(ss, step.ilu, true);
        }
    }

    // Output footer
    ss << "\n    VS_OUTPUT xOut;\n";
    ss << "#include \"CxbxVertexOutputFooter.hlsli\"\n";
    ss << "    return xOut;\n";
    ss << "}\n";

    return ss.str();
}

// ============================================================
// Validation: check if a program can be safely JIT'd
// Returns false if any instruction has missing operands or context writes
// ============================================================
static bool ValidateProgram(const uint32_t program_data[][4], uint32_t startAddr, uint32_t instrCount)
{
    for (uint32_t i = 0; i < instrCount; i++) {
        Nv2aVshStep step = {};
        Nv2aVshParseResult result = nv2a_vsh_parse_step(&step, program_data[startAddr + i]);
        if (result != NV2AVPR_SUCCESS) return false;

        // Check MAC operation inputs
        auto& mac = step.mac;
        if (mac.opcode != NV2AOP_NOP) {
            // All non-NOP operations need input A
            if (mac.inputs[0].type == NV2ART_NONE && mac.opcode != NV2AOP_NOP)
                return false;
            // Two-operand ops need input B
            switch (mac.opcode) {
            case NV2AOP_MUL: case NV2AOP_ADD: case NV2AOP_DP3: case NV2AOP_DP4:
            case NV2AOP_DPH: case NV2AOP_DST: case NV2AOP_MIN: case NV2AOP_MAX:
            case NV2AOP_SLT: case NV2AOP_SGE:
                if (mac.inputs[1].type == NV2ART_NONE) return false;
                break;
            case NV2AOP_MAD:
                if (mac.inputs[1].type == NV2ART_NONE) return false;
                if (mac.inputs[2].type == NV2ART_NONE) return false;
                break;
            default: break;
            }
            // Context writes not supported
            for (int o = 0; o < 2; o++) {
                if (mac.outputs[o].type == NV2ART_CONTEXT) return false;
            }
        }

        // Check ILU operation inputs
        auto& ilu = step.ilu;
        if (ilu.opcode != NV2AOP_NOP) {
            if (ilu.inputs[0].type == NV2ART_NONE) return false;
            for (int o = 0; o < 2; o++) {
                if (ilu.outputs[o].type == NV2ART_CONTEXT) return false;
            }
        }
    }
    return true;
}

// ============================================================
// Public API
// ============================================================
bool g_bEnableVSJIT = true; // Set false to disable JIT and use interpreter

ID3D11VertexShader* VertexShaderCache::GetShader(
    const uint32_t program_data[][4],
    uint32_t startAddr,
    ID3D11Device* pDevice,
    ID3DBlob** ppBytecode)
{
    *ppBytecode = nullptr;
    if (!g_bEnableVSJIT) return nullptr;

    // Determine program length (scan for final bit)
    uint32_t instrCount = 0;
    for (uint32_t i = 0; (startAddr + i) < 136; i++) {
        instrCount = i + 1;
        uint32_t dw3 = program_data[startAddr + i][3];
        if (dw3 & 1) // FINAL bit
            break;
    }

    if (instrCount == 0) return nullptr;

    // Validate: reject programs with missing operands or context writes
    if (!ValidateProgram(program_data, startAddr, instrCount))
        return nullptr;

    // Hash and check cache
    uint64_t hash = HashProgram(program_data, startAddr, instrCount);

    {
        std::lock_guard<std::mutex> lock(g_JITMutex);
        auto it = g_JITCache.find(hash);
        if (it != g_JITCache.end()) {
            *ppBytecode = it->second.pBytecode;
            if (*ppBytecode) (*ppBytecode)->AddRef();
            return it->second.pVS;
        }
    }

    // Cache miss — try disk cache first, then translate and compile
    InterlockedIncrement(&g_ProfileVSJITCompiles);
    CXBX_PROFILE_SCOPE(PROF_VS_JIT_COMPILE);

    ID3DBlob* pCode = ShaderDiskCache::TryLoad(hash);
    if (pCode) {
        // Disk cache hit — create shader from cached bytecode
        ID3D11VertexShader* pVS = nullptr;
        HRESULT hr = pDevice->CreateVertexShader(pCode->GetBufferPointer(), pCode->GetBufferSize(), nullptr, &pVS);
        if (SUCCEEDED(hr)) {
            std::lock_guard<std::mutex> lock(g_JITMutex);
            g_JITCache[hash] = { pVS, pCode };
            *ppBytecode = pCode;
            pCode->AddRef();
            return pVS;
        }
        pCode->Release();
        // Fall through to recompile if cached blob is invalid
    }

    std::string hlsl = TranslateToHLSL(program_data, startAddr, instrCount);
    if (hlsl.empty()) return nullptr; // Translation failed

    // Build include path for D3DCompile
    // CxbxJITIncludeHandler resolves includes relative to <exe_dir>\hlsl\
    // at runtime (CMake POST_BUILD copies all .hlsli headers there).
    pCode = nullptr;
    ID3DBlob* pErrors = nullptr;
    CxbxJITIncludeHandler includeHandler;

    HRESULT hr = D3DCompile(
        hlsl.c_str(), hlsl.size(),
        "JIT_VS.hlsl",    // source name (for error messages only)
        nullptr,           // defines
        &includeHandler,   // custom include handler
        "main",
        "vs_5_0",
        D3DCOMPILE_OPTIMIZATION_LEVEL3,
        0,
        &pCode,
        &pErrors);

    if (FAILED(hr)) {
        if (pErrors) {
            EmuLog(LOG_LEVEL::WARNING, "VS JIT compile failed: %s",
                (const char*)pErrors->GetBufferPointer());
            pErrors->Release();
        }
        std::lock_guard<std::mutex> lock(g_JITMutex);
        g_JITCache[hash] = { nullptr, nullptr };
        return nullptr;
    }
    if (pErrors) pErrors->Release();

    // Create the vertex shader
    ID3D11VertexShader* pVS = nullptr;
    hr = pDevice->CreateVertexShader(pCode->GetBufferPointer(), pCode->GetBufferSize(), nullptr, &pVS);
    if (FAILED(hr)) {
        EmuLog(LOG_LEVEL::WARNING, "VS JIT CreateVertexShader failed: 0x%08X", hr);
        pCode->Release();
        std::lock_guard<std::mutex> lock(g_JITMutex);
        g_JITCache[hash] = { nullptr, nullptr };
        return nullptr;
    }

    EmuLog(LOG_LEVEL::INFO, "VS JIT: compiled %u-instruction program (hash=%016llX)",
        instrCount, hash);

    // Save to disk cache for next session
    ShaderDiskCache::Save(hash, pCode);

    // Cache it
    {
        std::lock_guard<std::mutex> lock(g_JITMutex);
        g_JITCache[hash] = { pVS, pCode };
    }

    *ppBytecode = pCode;
    pCode->AddRef(); // caller gets a ref
    return pVS;
}
