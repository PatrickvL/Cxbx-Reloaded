#pragma once
// CxbxPixelShaderJIT.h — Runtime NV2A register combiner JIT compiler
//
// Translates NV2A PGRAPH register combiner topology (input routing,
// output destinations, texture modes) into straight-line HLSL, compiles
// via D3DCompile at runtime, and caches by state hash.
// The generated shader reads dynamic values (C0/C1 colors, fog, bump
// matrices) from g_PGRegs at runtime but eliminates all indirection,
// loops, register file indexing, and input mapping dispatch.

#include <d3d11.h>
#include <cstdint>

// Try to get a JIT-compiled PS for the current PGRAPH combiner state.
// Returns the compiled shader (cached on subsequent calls with same
// topology), or nullptr if JIT fails (caller falls back to interpreter).
ID3D11PixelShader* CxbxJITPixelShader(ID3D11Device* pDevice);

// Release all cached JIT pixel shaders (call on device reset)
void CxbxJITPixelShaderClearCache();
