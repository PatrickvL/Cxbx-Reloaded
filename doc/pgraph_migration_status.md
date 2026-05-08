# PGRAPH Migration Status

*Last updated: May 8, 2026*

---

## Current Architecture

### VertexShaderMode
- No explicit enum — mode determined at runtime by `NV2AIsFixedFunctionMode()` reading `NV_PGRAPH_CSV0_D_MODE` register (FIXED→FF, PROGRAM→SP)
- XYZRHW draws use MODE=PROGRAM with a trivial passthrough VS (CxbxVertexShaderPassthrough.hlsl) — JIT compiles it like any other program
- Removed: passthrough heuristics (VPSCL/CMAT), passthrough viewport bypass, g_bUsePassthroughHLSL

### Pixel Shader Path
- RC interpreter ubershader is the only pixel shader path (unconditional, no recompiler fallback)
- Deleted: PixelShader.cpp/h, CxbxPixelShaderTemplate.hlsl, XbPixelShader.cpp (recompiler)
- Retained: XbPixelShader.h (PS state enums), XbPixelShaderCompiler.cpp (compiler utilities, texture format fixup)
- PS selection: PGRAPH COMBINECTL != 0 → RC interpreter; COMBINECTL == 0 → fixed-function PS
- PS JIT cache accelerates common combiner topologies — see [jit_architecture.md](jit_architecture.md) for details

### State Authority
- `PGRAPHState::regs[]` is the ground truth for all GPU state
- All D3D EMUPATCH'es disabled — Xbox D3D runtime pushes state to PFIFO→PGRAPH natively
- No HLE bridges remain; renderer reads only from PGRAPH registers + small aux cbuffer

---

## Completed Migration Steps

| Step | Description | Status |
|------|-------------|--------|
| 1 | Remove OpenGL LLE backend | ✅ (commit f3c63999) |
| 2 | PFIFO flush primitive | ✅ (commit 181c53ed) |
| 3.1 | RC interpreter core combiner regs from PGRAPH | ✅ |
| 3.2 | RC interpreter post-processing from PGRAPH | ✅ (fog, alpha, BEM, LUM, final combiner) |
| 3.3 | Fog from PGRAPH CONTROL_3 | ✅ (commit d9311588) |
| 4.1 | VS microcode from PGRAPH program_data[] | ✅ (commit 55c680c3) |
| 4.2 | VS constants from PGRAPH vsh_constants[] | ✅ |
| 5.1 | Vertex fetch from PGRAPH vertex_attributes[] | ✅ (commit 203ba1e3) |
| 6.1 | Render target from PGRAPH surface_color/zeta | ✅ (commit 20e012c0) |
| 6.2 | Blend/depth-stencil/rasterizer from PGRAPH | ✅ (commit 0a2566a2) |
| 6.3 | Viewport/scissor from PGRAPH | ✅ (commit fa49e373) |
| 7.1 | Texture lookup from PGRAPH TEXOFFSET | ✅ (commit c09c2bd7) |
| 7.2 | Remove m_Textures device fallback | ✅ (commit 16c310e2) |
| 8.1 | Puller-driven draw_arrays | ✅ |
| 8.3 | Puller context flag (deadlock prevention) | ✅ |
| 9.2 | Shader recompiler removal | ✅ |
| PB | Push buffer migration (pfifo_submit_pushbuffer) | ✅ |
| 7.3 | Palette textures from PGRAPH DMA context | ✅ (commits 49cd110f5, 02a7982ae) |
| 9.3 | RunVertexStateShader → NV097_LAUNCH_TRANSFORM_PROGRAM | ✅ (commit d795c73fc) |
| 9.4 | BlockOnTime unpatched (native PFIFO flush) | ✅ (commit 13d460e35) |
| 10.1 | Backend reorganization (topic-oriented files) | ✅ (commit 9395cc3e9) |
| 10.2 | Shader JIT with disk cache (Backend/Shading/) | ✅ (commits 97fdc5eaa, f218b69b9, 09c536aa2) |
| 10.3 | SurfaceShape → NV2ASurfaceState (register-backed) | ✅ (commit 6010ff660) |
| 10.4 | DMA context A/B for texture+palette address resolution | ✅ (commit 49cd110f5) |
| 10.5 | Structured PGRAPH helper accessors (NV2A_PGRAPH_Helpers) | ✅ (commit fef2f05a5) |

### HLE Dependencies Fully Removed
- XboxRenderStates.Apply() — all 23 reads → PGRAPH registers
- XboxTextureStates.Apply() — sampler states from PGRAPH
- Fog (mode/density/start/end) → PGRAPH CONTROL_3 + FOGPARAM0/1
- FrontFace → PGRAPH CSV0_C + SETUPRASTER
- Point sprite → PGRAPH CONTROL_3 POINTPARAMSENABLE
- ALPHAKILL → PGRAPH TEXCTL0 per stage
- g_ZScale → removed entirely (uses NV2ASurfaceState.zetaFormat)
- g_pXbox_Palette_Data/Size → removed; PGRAPH reads palette via DMA context resolution
- 826 lines of dead Apply/SetDirty/Deferred code removed
- ~3600 lines of dead HLE vertex/patch infrastructure removed

---

## Remaining Work

### State Authority
All D3D EMUPATCHes are disabled. The renderer reads exclusively from PGRAPH registers
and PGRAPHState fields. No HLE patch is needed for rendering to function.

### Still HLE-sourced (no PGRAPH register equivalent)
- `ColorSign[4]` — host-side texture format compensation
- `TexFmtFixup` — host-side format compensation
- `ColorKeyOp/Color[4]` — Xbox D3D extension, no PGRAPH register
- `TEXCOORDINDEX` — FF/passthrough texcoord remapping

These are software-only concepts with no NV2A register backing; they remain in a small
auxiliary cbuffer (PSAuxCBLayout) uploaded alongside the PGRAPH register SRV.

---

## Sample Rendering Status (May 2026)

| Sample | Status | Notes |
|--------|--------|-------|
| Vertices | ✅ | RGB triangle, 60fps |
| Textures | ✅ | Textured cylinder, 35fps |
| Lights | ✅ | Yellow lit cylinder, 31fps |
| Matrices | ✅ | Rotating colored triangle |
| Meshes | ✅ | Textured tiger, 60fps |
| VertexShaders | ✅ | Color-interpolated cone, 45fps |
| BumpEarth | ✅ | Earth globe with bump+cloud textures |
| Dolphin | ✅ | Dolphin+ocean+sky |
| Ripple | ✅ | Animated rippling mesh (VS program) |
| Billboard | ✅ | Forest scene with billboarded trees |
| PixelShader | ✅ | RC interpreter, text+cylinder |
| PointSprites | ✅ | Fixed: blend factor + GS + sizing (commits cad5d46f8, a44400fd2) |
| CubeMap | ❌ | Black teapot (pre-existing) |
| TwoSidedLighting | ❌ | Text only, no geometry (pre-existing) |
| Fog | ❌ | No geometry (pre-existing) |
| MotionBlur | ❌ | Pre-existing |

---

## Known Issues

See [rc_unit_test_matrix.md](rc_unit_test_matrix.md) and [vs_unit_test_matrix.md](vs_unit_test_matrix.md) for the complete per-sample known issues list.

### Key Root Causes Fixed
1. **Gray Dolphin**: g_pXbox_PixelShader was NULL during puller draws → wrong PS path. Fix: COMBINECTL != 0 selects RC interpreter.
2. **Black screen**: SetHostResource missing D3DUsage → back buffer recreated blank every frame.
3. **CubeMap** (7 bugs): Stale PSTextureModes, missing RT bind flags, DS dimension mismatch, cubemap face routing.
4. **Shadow mapping**: D3D11 requires RTV/DSV dimension match → unbind RTV for depth-only rendering.

---

## Key Commits (dx11 branch)

| Commit | Description |
|--------|-------------|
| f3c63999 | Remove OpenGL LLE backend |
| 181c53ed | PFIFO flush primitive |
| 55c680c3 | VS interpreter program data from PGRAPH |
| 203ba1e3 | Vertex fetch from PGRAPH |
| 0a2566a2 | Blend/depth-stencil/rasterizer from PGRAPH |
| fa49e373 | Viewport/scissor from PGRAPH |
| c09c2bd7 | Texture lookup from PGRAPH TEXOFFSET |
| d9311588 | Fog from PGRAPH CONTROL_3 |
| a51fa2879 | Dead code removal (826 lines) |
| 5f39c0c3 | SetRenderState_Simple removed |
| c6d5221df | 31 state-caching patches disabled |
| 8fadaa438 | PGRAPH-driven clear |
| 70b30da8e | DOT_RFLCT_DIFF, CLIPPLANE, shadow compare fixes |
| 938ab3e3d | RT-as-texture SRV rebinding fix |
| cad5d46f8 | Fix point sprite rendering (blend, GS, textures, sizing) |
| a44400fd2 | Fix PointSpriteEnable register (SETUPRASTER not CONTROL_3) |
| 44bb61759 | Remove VertexShaderMode enum; use NV2AIsFixedFunctionMode() |
| 09e121dbc | Replace all hashing with rapidhash |
| d795c73fc | Implement NV097_LAUNCH_TRANSFORM_PROGRAM |
| 02a7982ae | Remove g_pXbox_Palette_Data globals, read PGRAPH directly |
| 49cd110f5 | DMA context A/B for texture and palette address resolution |
| 6010ff660 | Replace SurfaceShape with register-backed NV2ASurfaceState |
| fef2f05a5 | Add structured PGRAPH helper accessors (NV2A_PGRAPH_Helpers) |
| 17eb2426b | Remove ~3600 lines of dead HLE vertex/patch infrastructure |
| 09c536aa2 | Reorganize shader JIT into Backend/Shading with disk cache |
| 13d460e35 | Unpatch D3D_BlockOnTime: native PFIFO flush |
| 5bd6999d4 | GPU→CPU render target readback via VEH fault handler |

---

## PGRAPH Callback Architecture

Draw callbacks registered in XbPushBuffer.cpp `CxbxInitD3D11Renderer()`:
- `pgraph_draw` — Main draw dispatch (BEGIN_END with END)
- `pgraph_draw_state_update` — Pre-draw state sync (BEGIN_END with BEGIN)
- `pgraph_draw_clear` — NV097_CLEAR_SURFACE handler
- `pgraph_draw_patch` — Hardware tessellation (NV097_SET_BEGIN_PATCH)
- `pgraph_launch_transform_program` — Vertex state shader execution
- `pgraph_zpass_begin/end/collect` — Visibility test (occlusion query)

## FF Lighting from PGRAPH

- Light enable mask: `pg->regs[NV_PGRAPH_CSV0_D / 4] & NV_PGRAPH_CSV0_D_LIGHTS` (2 bits/light)
- NV2A type mapping: 1(INFINITE)→3(DIRECTIONAL), 2(LOCAL)→1(POINT), 3(SPOT)→2(SPOT)
- Colors from `pg->ltctxb[]` (pre-multiplied by material by Xbox D3D runtime)
- Material forced to white — ltctxb already contains light×material product
- Scene ambient from `pg->ltctxa[FR_AMB]` / `pg->ltctxa[BR_AMB]`
