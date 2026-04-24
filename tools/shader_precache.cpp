// shader_precache.cpp — Build-time tool to precompile interpreter shaders
// into the ShaderCache/00000000-All/ directory using the same hash key
// that EmuCompileShader() uses at runtime.
//
// Usage: shader_precache.exe <hlsl_dir> <output_cache_dir>
//
// Compiles the three static interpreter/passthrough shaders:
//   CxbxRegisterCombinerInterpreter.hlsl  (ps_5_0)
//   CxbxVertexShaderInterpreter.hlsl      (vs_5_0)
//   CxbxVertexShaderPassthrough.hlsl      (vs_5_0)
//
// Each .cso is named by XXH3_64bits(hlsl_source + "|" + profile), matching
// the key that EmuCompileShader() computes. At runtime the emulator will
// find these in the cache and skip compilation entirely.

#ifndef XXH_INLINE_ALL
#define XXH_INLINE_ALL
#endif
#include "xxhash.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <windows.h>

// D3DCompile types — loaded dynamically so we don't need to link d3dcompiler.lib
typedef HRESULT(WINAPI* pD3DCompile)(
    const void* pSrcData, SIZE_T SrcDataSize,
    const char* pSourceName,
    const void* pDefines,
    void* pInclude, // ID3DInclude*
    const char* pEntrypoint,
    const char* pTarget,
    UINT Flags1, UINT Flags2,
    void** ppCode,   // ID3DBlob**
    void** ppErrors); // ID3DBlob**

// Minimal ID3DBlob interface — we only need GetBufferPointer/Size/Release
struct ID3DBlob_Vtbl {
    // IUnknown
    HRESULT(__stdcall* QueryInterface)(void*, const void*, void**);
    ULONG(__stdcall* AddRef)(void*);
    ULONG(__stdcall* Release)(void*);
    // ID3DBlob
    void* (__stdcall* GetBufferPointer)(void*);
    SIZE_T(__stdcall* GetBufferSize)(void*);
};

static void* BlobGetPtr(void* blob) {
    auto vtbl = *(ID3DBlob_Vtbl**)blob;
    return vtbl->GetBufferPointer(blob);
}
static SIZE_T BlobGetSize(void* blob) {
    auto vtbl = *(ID3DBlob_Vtbl**)blob;
    return vtbl->GetBufferSize(blob);
}
static void BlobRelease(void* blob) {
    auto vtbl = *(ID3DBlob_Vtbl**)blob;
    vtbl->Release(blob);
}

// D3DCOMPILE flags (from d3dcompiler.h)
// Level encoding: bits 14 and 15 form a 2-bit field.
//   O0 = (1 << 14)            = 0x4000
//   O1 = 0                    = 0x0000  (default)
//   O2 = (1 << 14)|(1 << 15)  = 0xC000
//   O3 = (1 << 15)            = 0x8000
#define MY_D3DCOMPILE_OPTIMIZATION_LEVEL0 (1 << 14)
#define MY_D3DCOMPILE_OPTIMIZATION_LEVEL1 0
#define MY_D3DCOMPILE_OPTIMIZATION_LEVEL3 (1 << 15)
// D3D_COMPILE_STANDARD_FILE_INCLUDE = (ID3DInclude*)(UINT_PTR)1
#define MY_D3D_COMPILE_STANDARD_FILE_INCLUDE ((void*)(UINT_PTR)1)

struct ShaderJob {
    const char* filename;
    const char* profile;
};

static const ShaderJob g_Jobs[] = {
    { "CxbxRegisterCombinerInterpreter.hlsl", "ps_5_0" },
    { "CxbxVertexShaderInterpreter.hlsl",     "vs_5_0" },
    { "CxbxVertexShaderPassthrough.hlsl",     "vs_5_0" },
    { "CxbxFixedFunctionVertexShader.hlsl",   "vs_5_0" },
};

static std::string ReadFile(const std::string& path) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return {};
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    std::string buf(sz, '\0');
    fread(&buf[0], 1, sz, f);
    fclose(f);
    return buf;
}

static bool WriteFile_(const std::string& path, const void* data, size_t size) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    size_t written = fwrite(data, 1, size, f);
    fclose(f);
    return written == size;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <hlsl_dir> <output_cache_dir>\n", argv[0]);
        return 1;
    }

    std::string hlslDir = argv[1];
    std::string cacheDir = argv[2];

    // Ensure trailing backslash
    if (!hlslDir.empty() && hlslDir.back() != '\\' && hlslDir.back() != '/')
        hlslDir += '\\';
    if (!cacheDir.empty() && cacheDir.back() != '\\' && cacheDir.back() != '/')
        cacheDir += '\\';

    // Create output directory
    CreateDirectoryA(cacheDir.c_str(), nullptr);

    // Load d3dcompiler_47.dll
    HMODULE hD3DCompiler = LoadLibraryW(L"d3dcompiler_47.dll");
    if (!hD3DCompiler) {
        fprintf(stderr, "ERROR: Cannot load d3dcompiler_47.dll\n");
        return 1;
    }
    auto pfnD3DCompile = (pD3DCompile)GetProcAddress(hD3DCompiler, "D3DCompile");
    if (!pfnD3DCompile) {
        fprintf(stderr, "ERROR: Cannot find D3DCompile in d3dcompiler_47.dll\n");
        return 1;
    }

    int failures = 0;
    int successes = 0;

    for (const auto& job : g_Jobs) {
        std::string hlslPath = hlslDir + job.filename;
        std::string hlsl = ReadFile(hlslPath);
        if (hlsl.empty()) {
            fprintf(stderr, "WARNING: Cannot read %s, skipping\n", hlslPath.c_str());
            failures++;
            continue;
        }

        // Compute cache key: XXH3_64bits(hlsl + "|" + profile)
        // Must match EmuCompileShader() in Shader.cpp
        std::string cacheInput = hlsl + "|" + job.profile;
        uint64_t cacheHash = XXH3_64bits(cacheInput.c_str(), cacheInput.size());

        // Compile at maximum optimization (O3). At build time we can afford the
        // extra compile time; the resulting DXBC bytecode is a standardized
        // intermediate format that all D3D11 drivers (and vkd3d/DXVK) consume,
        // so O3 output is universally compatible.
        // Cascade: O3 -> O1 -> O0, matching the runtime fallback strategy.
        struct OptLevel { UINT flags; const char* name; };
        static const OptLevel levels[] = {
            { MY_D3DCOMPILE_OPTIMIZATION_LEVEL3, "O3" },
            { MY_D3DCOMPILE_OPTIMIZATION_LEVEL1, "O1" },
            { MY_D3DCOMPILE_OPTIMIZATION_LEVEL0, "O0" },
        };

        void* pBlob = nullptr;
        void* pErrors = nullptr;
        HRESULT hr = E_FAIL;
        const char* usedLevel = "?";

        for (const auto& level : levels) {
            pBlob = nullptr;
            pErrors = nullptr;
            hr = pfnD3DCompile(
                hlsl.c_str(), hlsl.size(),
                hlslPath.c_str(),
                nullptr,
                MY_D3D_COMPILE_STANDARD_FILE_INCLUDE,
                "main",
                job.profile,
                level.flags,
                0,
                &pBlob,
                &pErrors);

            if (SUCCEEDED(hr)) {
                usedLevel = level.name;
                if (pErrors) BlobRelease(pErrors);
                break;
            }
            // Log failure and try next level
            fprintf(stderr, "  %s compile failed for %s", level.name, job.filename);
            if (pErrors) {
                fprintf(stderr, ": %s", (const char*)BlobGetPtr(pErrors));
                BlobRelease(pErrors);
            }
            fprintf(stderr, "\n");
        }

        if (FAILED(hr)) {
            fprintf(stderr, "ERROR: All optimization levels failed for %s\n", job.filename);
            failures++;
            continue;
        }

        // Write <hash>.cso
        char csoFilename[32];
        snprintf(csoFilename, sizeof(csoFilename), "%016llx.cso", (unsigned long long)cacheHash);
        std::string csoPath = cacheDir + csoFilename;

        if (!WriteFile_(csoPath, BlobGetPtr(pBlob), BlobGetSize(pBlob))) {
            fprintf(stderr, "ERROR: Cannot write %s\n", csoPath.c_str());
            BlobRelease(pBlob);
            failures++;
            continue;
        }

        printf("OK: %s (%s, %s) -> %s (%.1f KB)\n",
            job.filename, job.profile, usedLevel, csoFilename,
            BlobGetSize(pBlob) / 1024.0);

        BlobRelease(pBlob);
        successes++;
    }

    FreeLibrary(hD3DCompiler);

    printf("\nPrecache complete: %d succeeded, %d failed\n", successes, failures);
    return failures > 0 ? 1 : 0;
}
