#pragma once
// PixelShaderCache.h — Runtime NV2A register combiner JIT compiler + cache
//
// Translates NV2A PGRAPH register combiner topology (input routing,
// output destinations, texture modes) into straight-line HLSL, compiles
// via D3DCompile at runtime, and caches by state hash.

#include <d3d11.h>
#include <cstdint>

class PixelShaderCache {
public:
    // Get (or compile on miss) a JIT'd pixel shader for the current PGRAPH state.
    // Returns nullptr if JIT fails (caller falls back to interpreter).
    ID3D11PixelShader* GetShader(ID3D11Device* pDevice);

    // Release all cached JIT pixel shaders (call on device reset)
    void Clear();
};

extern PixelShaderCache g_PixelShaderCache;
