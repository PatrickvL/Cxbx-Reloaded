#pragma once
// CxbxVertexShaderJIT.h — Runtime NV2A vertex shader JIT compiler
//
// Translates NV2A transform program microcode into native HLSL, compiles
// via D3DCompile at runtime, and caches by program hash. This replaces
// the GPU interpreter loop for 10-100x faster vertex shading.

#include <d3d11.h>
#include <d3dcompiler.h>
#include <cstdint>

// Try to get a JIT-compiled VS for the current PGRAPH program.
// Returns the compiled shader (cached on second call), or nullptr if
// JIT fails (caller should fall back to interpreter).
// Also returns the bytecode blob (for input layout creation) via ppBytecode.
ID3D11VertexShader* CxbxJITVertexShader(
    const uint32_t program_data[][4],
    uint32_t startAddr,
    ID3D11Device* pDevice,
    ID3DBlob** ppBytecode);

// Release all cached JIT shaders (call on device reset)
void CxbxJITVertexShaderClearCache();
