# PGRAPH Migration Status

*Last updated: May 16, 2026*

---

## Current Architecture

### VertexShaderMode
- No explicit enum — mode determined at runtime by `NV2AIsFixedFunctionMode()` reading `NV_PGRAPH_CSV0_D_MODE` register (FIXED→FF, PROGRAM→SP)
- XYZRHW draws use MODE=PROGRAM with a trivial passthrough VS (CxbxVertexShaderPassthrough.hlsl) — JIT compiles it like any other program
- Removed: passthrough heuristics (VPSCL/CMAT), passthrough viewport bypass, g_bUsePassthroughHLSL

### Pixel Shader Path
- RC interpreter ubershader is the only pixel shader path (unconditional, no recompiler fallback)
- Deleted: PixelShader.cpp/h, CxbxPixelShaderTemplate.hlsl, XbPixelShader.cpp (recompiler)
- Retained: XbPixelShader.h (PS state enums), Backend_D3D11_PixelShader.cpp (compiler utilities, texture format fixup)
- PS selection: always RC interpreter/JIT — PGRAPH combiners are authoritative regardless of COMBINECTL
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
auxiliary cbuffer (PSAuxCBLayout, 112 bytes) uploaded alongside the PGRAPH register SRV.
PSAuxCBLayout also contains `DepthScale` and `DepthTexAlias` for depth texture handling.

---

## Rendering Status

See [rendering_test_status.md](rendering_test_status.md) for the consolidated XDK sample status,
known issues, fixed issues, key commits, and game compatibility tracking.

---

## PGRAPH Callback Architecture

Draw callbacks registered in XbPushBuffer.cpp `D3D11_init_pgraph_plugins()`:
- `draw` — Main draw dispatch (BEGIN_END with END)
- `draw_state_update` — Pre-draw state sync (BEGIN_END with BEGIN)
- `draw_clear` — NV097_CLEAR_SURFACE handler
- `draw_patch` — Hardware tessellation (NV097_SET_BEGIN_PATCH)
- `flip_stall` — PCRTC flip synchronization
- `launch_transform_program` — Vertex state shader execution
- `zpass_begin/end/collect` — Visibility test (occlusion query)

## FF Lighting from PGRAPH

- Light enable mask: `pg->regs[NV_PGRAPH_CSV0_D / 4] & NV_PGRAPH_CSV0_D_LIGHTS` (2 bits/light)
- NV2A type mapping: 1(INFINITE)→3(DIRECTIONAL), 2(LOCAL)→1(POINT), 3(SPOT)→2(SPOT)
- Colors from `pg->xf.ltctxb[]` (pre-multiplied by material by Xbox D3D runtime)
- Material forced to white — ltctxb already contains light×material product
- Scene ambient from `pg->xf.ltctxa[FR_AMB]` / `pg->xf.ltctxa[BR_AMB]`
