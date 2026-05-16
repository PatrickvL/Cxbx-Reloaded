# JIT Shader Architecture

## VS JIT (VertexShaderCache.cpp)

- Translates NV2A transform program microcode → straight-line HLSL → D3DCompile (vs_5_0, O3)
- Cache key: rapidhash of program tokens
- Entry: `g_VertexShaderCache.GetShader()` called from Backend_D3D11_VertexShader.cpp
- Fallback: VS interpreter when JIT misses or fails
- Source: `src/core/hle/D3D8/Rendering/Backend/Shading/VertexShaderCache.cpp`

## PS JIT (PixelShaderCache.cpp)

- Translates PGRAPH register combiner TOPOLOGY → straight-line HLSL → D3DCompile (ps_5_0, O3)
- Bakes: input routing, output destinations, texture modes, final combiner structure
- Dynamic (still from g_PGRegs): C0/C1 constants, FogColor, bump matrices
- Cache key: rapidhash of PSJITKey struct (topology state)
- Entry: `g_PixelShaderCache.GetShader()` called from `CxbxUpdateActivePixelShader()` in Backend_D3D11_PixelShader.cpp
- Bindings: g_PGRegs (t12), PSAuxCBLayout (b0), textures (t0-t11), samplers (s0-s3)
- Fallback: RC interpreter ubershader (g_pD3D11RCInterpreterPS)
- Source: `src/core/hle/D3D8/Rendering/Backend/Shading/PixelShaderCache.cpp`

## State Flow

1. `CxbxUpdateActivePixelShader()` calls `CxbxD3D11UploadRCInterpreterState()` (uploads aux CB + PGRegs)
2. Tries `g_PixelShaderCache.GetShader()` — reads g_LastPSAuxCB for topology hash
3. On JIT hit: sets JIT PS and returns
4. On JIT miss/fail: falls back to interpreter (g_pD3D11RCInterpreterPS)

## Variant Removal

- Removed: RC PS variants (S2T2/S2T4/S4T2/S4T4/S8T2/S8T4), VS NoCtx variant
- Kept: base RC interpreter (fallback), base VS interpreter (fallback)
- JIT replaces the variant optimization with per-program specialization

## Shader Disk Cache (ShaderDiskCache.cpp)

- Persists compiled D3D11 bytecode to disk to avoid recompilation across runs
- Path: `src/core/hle/D3D8/Rendering/Backend/Shading/ShaderDiskCache.cpp`
- Indexed by rapidhash of the JIT cache key
- Loads on startup; writes on new JIT compilations
- Avoids D3DCompile() calls for previously-seen shaders (major startup speedup)

## JIT ↔ Interpreter Remaining Differences

### VS
| Area | Difference | Impact |
|------|-----------|--------|
| Context writes | JIT blocks; Interp supports | Falls back gracefully |

Note: RCC sign preservation and LOG(0) bugs were fixed in commit 39fd1b808. Both JIT and interpreter now share identical implementations via CxbxNV2AVshOps.hlsli.

### PS
| Area | Difference | Impact |
|------|-----------|--------|
| nv2a_mul zero×inf | Both return 0 (xemu returns NaN) | Correct for NV2A hardware |
