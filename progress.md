# Vertex Distortion Investigation — Progress Tracker

## Problem Statement
GTA San Andreas shows vertex coordinate distortion (spiked/exploded vertices on character models).
The `pfifo_submit_pushbuffer` path processes push buffer commands through NV2A PGRAPH.

## Critical Architecture Finding

**Push buffer draws use a SPLIT state model:**
- **Vertex attributes** (offset, stride, format): Read from PGRAPH `pg->vertex_attributes[]` ✓
- **Vertex defaults** (inline_value): Read from PGRAPH ✓
- **Vertex shader program**: Read from HLE `g_Xbox_VertexShader_FunctionSlots[]` ✗ (NOT PGRAPH `pg->program_data[][]`)
- **Render states** (blend, cull, skin, etc.): Read from HLE `XboxRenderStates` ✗
- **Textures**: Read from HLE state ✗
- **Viewport, RT, pipeline**: PGRAPH overrides are `#if 0` DISABLED — HLE state only ✗

`HLE_draw_state_update()` has `LOG_INCOMPLETE()` — only reads fog color from PGRAPH.
`CxbxUpdateNativeD3DResources()` syncs everything from HLE-tracked state, not PGRAPH.

**Exception**: The VS **interpreter** path reads `pg->program_data[][]` directly (bypasses HLE slots).
A TODO in the code notes: "use pg->program_data[Address] here instead" of HLE slots.

## Root Causes (Re-Ranked After Investigation)

### RC6: Vertex Shader / Render State Desync (NEW — HIGHEST PRIORITY)
**Status**: IDENTIFIED — needs fix  
**Description**: When push buffer commands set vertex shader program (NV097_SET_TRANSFORM_PROGRAM),
render states, textures, etc., these go into PGRAPH registers but are NOT read by the HLE sync path.
The compile path for VS reads from `g_Xbox_VertexShader_FunctionSlots[]` which is only set by HLE
patches (SetVertexShader/SelectVertexShader). Push buffer VS uploads go to `pg->program_data[][]`
which is only used by the VS interpreter path. If a push buffer draws with a different VS than
the last HLE-set VS, the wrong shader is used → wrong vertex transformation → vertex distortion.
**Impact**: Could cause vertex distortion if push buffer uploads/changes vertex shader programs
**XDK Samples**: PushBuffer sample works because it uses same VS set via D3D8 API before recording
**Fix Direction**: Either (a) sync VS program from PGRAPH before push buffer draws, or (b) always
use the VS interpreter path for push buffer draws, or (c) read pg->program_data in compile path

### RC1: Typed SRV Misalignment (Demoted — Unlikely for Common Layouts)
**Status**: INVESTIGATED — low risk  
**Description**: `uint elem = byteOff / 4u` in typed SRV path truncates for non-aligned access.  
**Finding**: Requires `(attr.offset + vtxIdx * stride) % 4 != 0`. With diagnostics enabled, the
PushBuffer sample shows all offsets and strides are 4-byte aligned. NV2A hardware typically uses
aligned vertex layouts. The bug is REAL in theory but unlikely to trigger for common game layouts.
Still need GTA SA diagnostics to confirm.
**Diagnostic added**: `C:\Temp\vtx_diag.txt` — flags `[MISALIGNED_OFFSET]`, `[MISALIGNED_STRIDE]`

### RC2: DMA Context Base Address Ignored (Demoted — By Design)
**Status**: INVESTIGATED — low risk by design  
**Description**: DMA context methods (NV097_SET_CONTEXT_DMA_VERTEX_A/B) are intentionally commented
out in EmuNV2A_PGRAPH.cpp with note: "vertex attribute offsets already store physical byte offsets
into contiguous memory, and both DMA contexts point to the same region."
`dma_select` is stored but never read. No `dma_vertex_a`/`dma_vertex_b` fields in PGRAPHState.
**Finding**: Both DMA contexts on Xbox map all of physical memory with base=0. By design, no fix needed.

### RC3: Stride=0 Constant Attribute (Possible — PGRAPH path affected)
**Status**: INVESTIGATED — possible but secondary  
**Description**: PGRAPH path passes stride=0 directly to shader. HLE path has fallback:
`if (stride == 0) stride = streamInfo.HostVertexStride`. Shader with stride=0 computes
`byteOff = streamBase + vtxIdx * 0 + elemOffset = streamBase + elemOffset` (same for all vertices).
This is correct if the single value is valid, but wrong if the game intends inline_value to be used.
**Diagnostic added**: `[STRIDE_ZERO]` flag in vtx_diag.txt

### RC4: Skinning Mode Desync (Subsumed by RC6)
**Status**: Subsumed by RC6 — skinning mode is a render state, which uses stale HLE state for
push buffer draws. If push buffer changes skin mode, HLE render state won't reflect it.

### RC5: NORMPACKED3 Decode (Ruled Out)
**Status**: TESTED — CompressedVertices sample renders correctly. Teapot with NORMPACKED3 normals
displays correct lighting. Not the cause of vertex spikes.

---

## XDK Sample Test Results

### PGRAPH Path (Push Buffer Draws)
Only 1 of 20 tested samples uses the PGRAPH attribute path:

| Sample | PGRAPH Attrs | Formats | Result |
|--------|-------------|---------|--------|
| PushBuffer | 2 (attr[0], attr[3]) | FLOAT3 + D3DCOLOR, stride 16, aligned | PASS ✓ |

### HLE Path (Standard D3D8 API Draws)
All other samples use the HLE fallback (PGRAPH vertex_attributes have count=0):

| Sample | Vertex Formats | Result | Notes |
|--------|---------------|--------|-------|
| CompressedVertices | FLOAT3 + NORMPACKED3 | PASS ✓ | Teapot with packed normals |
| MatrixPaletteSkinning | FLOAT3/2 + SHORT3 (multi-stream) | PASS ✓ | Skinned model |
| VertexBlend | FLOAT3 + FLOAT1 blend | PASS ✓ | MS logo with vertex blending |
| Dolphin | FLOAT3×3 streams (morphing) | PASS ✓ | Multi-stream vertex fetch |
| DisplacementMap | D3DCOLOR×3 streams | PASS ✓ | D3DCOLOR as indices/weights |
| Explosion | Standard formats | PASS ✓ | Complex particle effects |
| Trees | Standard formats | PASS ✓ | Billboard trees |
| ShadowBuffer | Standard formats | BLACK | 0 FPS — unrelated issue |
| BumpDemo | Standard formats | PASS ✓ | Bump-mapped gears |
| Cartoon | Standard formats | PASS ✓ | Toon-shaded character model |
| Water | Standard formats | PASS ✓ | Reflective water surface |
| Fur | Standard formats | PASS ✓ | Fur rendering (3 FPS) |
| Lights | Standard formats | PASS ✓ | Basic lit geometry |
| Meshes | Standard formats | PASS ✓ | Mesh rendering |
| Strip | Standard formats | PASS ✓ | Triangle strip rendering |
| PointSprites | (verified earlier) | PASS ✓ | Push buffer path |
| BeginPush | (verified earlier) | PASS ✓ | Inline push path |
| Fog | (verified earlier) | PASS ✓ | Fog modes |

### Format Coverage Gap
**No XDK sample exercises typed SRV formats (SHORT2N, SHORT4N, PBYTE1-4) via PGRAPH path.**
All HLE-path samples use standard D3D8 API which bypasses PGRAPH attribute setup.
GTA SA likely uses compact vertex formats (D3DCOLOR, SHORT2N for blend data) via push buffers.

---

## Diagnostic Instrumentation

**File**: `Backend_D3D11_IABypass.cpp` (PGRAPH attribute setup path)  
**Output**: `C:\Temp\vtx_diag.txt` (first 5000 draws)  
**Format**: Per-draw header + per-attribute details  
**Flags**: `[MISALIGNED_OFFSET]`, `[MISALIGNED_STRIDE]`, `[STRIDE_ZERO]`, `*** ISSUE ***`

To use: run any game, then check `C:\Temp\vtx_diag.txt` for flagged issues.

---

## Investigation Log

### Session 1 (Current)
1. Surveyed 97 compiled XDK samples, 72 with source code
2. Analyzed vertex format usage across all Graphics samples
3. Studied CxbxVertexFetch.hlsli and Backend_D3D11_IABypass.cpp in detail
4. **Key discovery**: CxbxUpdateNativeD3DResources syncs from HLE state, not PGRAPH
   - Vertex shader compile path reads HLE `g_Xbox_VertexShader_FunctionSlots[]`
   - Push buffer VS uploads go to PGRAPH `pg->program_data[][]`  
   - `HLE_draw_state_update()` has LOG_INCOMPLETE — only reads fog color from PGRAPH
   - PGRAPH viewport/RT/pipeline overrides are `#if 0` disabled
5. Added vertex attribute diagnostic to PGRAPH path
6. Tested 18 XDK samples — all pass except ShadowBuffer (black screen, separate issue)
7. Only PushBuffer sample exercises PGRAPH attribute path; all others use HLE fallback
8. CompressedVertices confirms NORMPACKED3 decode is correct (RC5 ruled out)

## Next Steps (Priority Order)

### 1. Run GTA SA with diagnostics
Run the game and check `C:\Temp\vtx_diag.txt` to see:
- What vertex formats GTA SA uses in push buffer draws
- Whether any misaligned offsets or strides are flagged
- Whether stride=0 is used for any attributes
- How many attributes are active (indicates vertex complexity)

### 2. Investigate RC6 — Vertex Shader Desync
This is the most likely cause of vertex distortion for push buffer draws:
- Check which VS mode is active during push buffer draws (compile vs interpreter)
- If compile mode: VS program from HLE slots may not match push buffer VS upload
- Options: (a) sync HLE slots from `pg->program_data[][]` in draw_state_update,
  (b) force interpreter mode for push buffer draws, (c) mirror program_data to HLE slots
  in pgraph_handle_method for NV097_SET_TRANSFORM_PROGRAM

### 3. Unpatch D3D functions incrementally (informed by reference branch)
Remove EMUPATCH'es that conflict with push buffer rendering, starting with the
highest-impact ones. The reference branch provides the order and approach.

### 4. Fix identified issues
Based on diagnostics from steps 1-3, implement fixes for confirmed root causes.

---

## Reference Branch Analysis: `jackchentwkh/pushbuffer_based_rendering`

**Branch**: `jackchentwkh/pushbuffer_based_rendering` (503 commits, D3D9 era, old LLE path)  
**Added as remote**: `git remote add jackchentwkh https://github.com/jackchentwkh/Cxbx-Reloaded.git`  
**Related branches**: `Pushbuffer_Restructure_HLE_Rebase_V3` (202 commits, cleaner rebase)

### Key Architecture Decisions

1. **Push buffer commands flow through NV2A PGRAPH** — same as our `pfifo_submit_pushbuffer`
2. **HLE draw callbacks** connected to PGRAPH draw triggers (same as our `pgraph_draw_arrays` etc.)
3. **Most D3D8 SetRenderState/SetTextureState patches REMOVED** — Xbox code writes push buffer
   commands that flow through PGRAPH, so state tracks itself via NV2A registers
4. **Resource management patches KEPT** — Create/Lock/Destroy still needed for host GPU mirroring
5. **Draw call patches KEPT** — DrawVertices, DrawIndexedVertices, DrawVerticesUP still patched
   (needed to capture draw parameters and trigger host rendering)
6. **Vertex shader patches partially removed** — SetVertexShaderConstant unpatched (flows through
   push buffer), but SetVertexShader/SelectVertexShader kept on the branch's final state
7. **Defined HLE API enum (`X_D3DAPI_ENUM`)** for push buffer token-based HLE method dispatch

### Patches REMOVED in Reference Branch (Moved to `Direct3D9.cpp.unused-patches`)

**Render State patches removed** (35+):
- ALL individual `SetRenderState_*` patches (CullMode, FillMode, ZEnable, StencilEnable, etc.)
- `SetRenderState_Deferred` — the batched deferred state path
- ALL `SetTextureState_*` patches (TexCoordIndex, BorderColor, BumpEnv, ColorKeyColor)
- `SetScissors`, `SetTile`, `SetViewport`

**Pixel Shader patches removed**:
- `SetPixelShader`, `SetPixelShaderConstant`, `SetPixelShaderProgram`
- `CreatePixelShader`, `DeletePixelShader`

**Resource management patches removed**:
- `CreateTexture`, `CreateCubeTexture`, `CreateVolumeTexture`, `CreateVertexBuffer`,
  `CreateIndexBuffer`, `CreatePalette`, `CreateImageSurface`, `CreateStateBlock`
- ALL Lock/GetDesc patches (Lock2DSurface, LockRect, LockBox, GetDesc, GetSurfaceLevel)
- `D3DResource_AddRef`, `D3DResource_IsBusy`, `D3DResource_GetType`
- `SetPalette` (replaced with LLE palette setup from NV2A)

**Push buffer management patches removed**:
- `KickOff`, `KickPushBuffer`, `MakeSpace`, `MakeRequestedSpace`, `XMETAL_StartPush`
- `BeginPushBuffer`, `EndPushBuffer` (push buffer recording let through to hardware)
- `GetPushBufferOffset`

**Vertex shader patches removed** (later re-added with modifications):
- `SetVertexShaderConstant` (all variants)
- `SelectVertexShaderDirect`, `SetVertexShaderInputDirect`

**Display/Misc patches removed**:
- `GetViewport`, `GetViewportOffsetAndScale`
- `SetBackMaterial`, `GetBackMaterial`, `GetLight`
- `GetDisplayMode`, `GetDeviceCaps`, `GetCreationParameters`
- `SetGammaRamp`, `SetFlickerFilter`, `SetSoftDisplayFilter`, `PersistDisplay`

### Patches KEPT in Reference Branch (Still in `Direct3D9.cpp`)

**Always needed — Host GPU resource management**:
- `Direct3D_CreateDevice` — initializes host D3D device
- `D3DDevice_Clear` — needs host RT translation
- `D3DDevice_CopyRects` — host surface blit
- `D3DDevice_Present`, `D3DDevice_Swap` — present to host window
- `D3DDevice_SetRenderTarget`, `D3DDevice_SetRenderTargetFast` — host RT binding
- `D3DDevice_EnableOverlay`, `D3DDevice_UpdateOverlay` — overlay display
- `D3DDevice_SetBackBufferScale` — host backbuffer sizing
- `D3DResource_BlockUntilNotBusy` — sync primitive
- `D3D_DestroyResource` — host resource cleanup
- `D3D_BlockOnTime` — timing sync
- `D3D_CommonSetRenderTarget` — shared RT setup
- `D3DDevice_BlockOnFence`, `D3DDevice_InsertFence`, `D3DDevice_IsFencePending` — GPU sync
- `D3DDevice_BlockUntilVerticalBlank`, `D3DDevice_SetSwapCallback`,
  `D3DDevice_SetVerticalBlankCallback` — VBlank sync

**Draw calls — still needed to trigger host rendering**:
- `D3DDevice_DrawVertices`, `D3DDevice_DrawVerticesUP`
- `D3DDevice_DrawIndexedVertices`, `D3DDevice_DrawIndexedVerticesUP`
- `D3DDevice_DrawRectPatch`, `D3DDevice_DrawTriPatch`
- `D3DDevice_Begin`, `D3DDevice_End`
- `D3DDevice_SetVertexData*` (2f, 2s, 4f, 4s, 4ub, Color) — inline vertex submission
- `D3DDevice_FlushVertexCache`, `D3DDevice_PrimeVertexCache`

**Push buffer recording — needed for RunPushBuffer**:
- `D3DDevice_BeginPush_4/8`, `D3DDevice_EndPush`
- `D3DDevice_BeginPushBuffer`, `D3DDevice_EndPushBuffer`
- `D3DDevice_RunPushBuffer`

**Vertex shader management — HLE bookkeeping**:
- `D3DDevice_CreateVertexShader`, `D3DDevice_DeleteVertexShader`
- `D3DDevice_SetVertexShader` (kept, but branch attempted to unpatch it)
- `D3DDevice_SelectVertexShader`
- `D3DDevice_LoadVertexShader`, `D3DDevice_LoadVertexShaderProgram`
- `D3DDevice_SetVertexShaderConstant*` (re-added after initial unpatch)
- `D3DDevice_SetVertexShaderInput`
- `D3DDevice_RunVertexStateShader`
- `D3DDevice_GetVertexShader*` (Get variants for readback)

**Stream/Index setup — HLE state tracking**:
- `D3DDevice_SetStreamSource` (all variants)
- `D3DDevice_SetIndices`
- `D3DDevice_SetTexture`, `D3DDevice_SwitchTexture`

**Fixed-function state — HLE tracking for FF VS**:
- `D3DDevice_SetTransform`, `D3DDevice_MultiplyTransform`
- `D3DDevice_SetLight`, `D3DDevice_LightEnable`
- `D3DDevice_SetMaterial`, `D3DDevice_GetMaterial`
- `D3DDevice_SetModelView`, `D3DDevice_GetModelView`
- `D3DDevice_SetShaderConstantMode`
- `D3DDevice_SetScreenSpaceOffset`
- `D3DDevice_SetDepthClipPlanes`
- `D3DDevice_SetStipple`

**Pixel shader — kept for HLE bridge**:
- `D3DDevice_SetPixelShader` (kept in active code)

**Other kept**:
- `D3DDevice_SetRenderState_Simple` — the main deferred RS write path
- `D3D_LazySetPointParams`, `D3D_SetCommonDebugRegisters`
- `D3DDevice_BeginVisibilityTest`, `D3DDevice_EndVisibilityTest`,
  `D3DDevice_GetVisibilityTestResult`
- `D3DDevice_ApplyStateBlock` (kept despite CreateStateBlock removed)
- `D3DResource_Register`, `D3DTexture_GetSurfaceLevel*`
- `CDevice_SetStateUP`, `CDevice_SetStateVB` — internal D3D state setup
- `Lock2DSurface`, `Lock3DSurface` (kept in active code despite being in unused-patches too)

### Key Unpatch Commits (Chronological Order)

1. `ff0d54005` — Unpatch `D3DDevice_SetVertexShaderConstant()` (VS constants flow through push buffer)
2. `e5ee08cef` — Unpatch VS-related APIs + SetVertexShaderConstant variants, add passthrough VS
3. `053fbcec4` — Unpatch fence-related APIs (InsertFence, BlockOnFence, IsFencePending)
4. `03eaa98a0` — Unpatch `XGSetTextureHeader`, `XGSetSurfaceHeader`, `XGSetVertexBufferHeader`,
   `SetTile`, `GetTile`
5. `d7650d13c` — Unpatch `GetPersistedSurface` (let kernel handle it)
6. `705f9034e` — Unpatch `Lock2DSurface` (trying to fix DynamicGamma)
7. `da8f2bcd2` — Unpatch `SetGammaRamp`, implement from PGRAPH
8. `66e10cd8a` — Unpatch `SetShaderConstantMode` (constant mode tracked via push buffer)
9. `d332e0b42` — Unpatch `SetStreamSource` and variants
10. `2730adb06` — Unpatch `LoadVertexShader`, `SelectVertexShader` (VS loads via push buffer)
11. `3ad509d21` — Unpatch `RunPushBuffer` (avoid reentrance of pgraph_handle_method)
12. `87796fc56` — Unpatch `SetMaterial`, `GetMaterial`
13. `e866a56a5` — Unpatch `SetVertexShader` (biggest single change — requires D3D globals init)
14. `1f0feed25` — Unpatch `SetTransform`, remove Xbox view matrix usage in FF VS
15. `969867d23` — Unpatch `SetViewport` (viewport composed from NV2A)
16. `7e597ee4e` — Unpatch `SetScreenSpaceOffset` (variables never used)
17. `a9c138dc7` — Unpatch `RunVertexStateShader`, implement NV097_LAUNCH_TRANSFORM_PROGRAM

### Lessons from the Reference Branch

1. **Incremental approach**: Patches were removed one at a time with testing between each
2. **SetVertexShaderConstant was first** — VS constants are the lowest-risk unpatch because
   they flow directly through push buffer → PGRAPH → already read from `pg->vsh_constants[]`
3. **SetRenderState patches were bulk-removed** — all individual RS patches removed at once,
   replaced by reading from Xbox kernel `D3D__RenderState[]` memory (which Xbox code writes to)
4. **SetVertexShader was the hardest** — required initializing D3D globals in CreateDevice_End
   and careful separation of NV2A vs Xbox globals
5. **Some unpatches were REVERTED** (e.g., `BlockUntilVerticalBlank` was unpatched then reverted)
6. **Resource patches stayed** — texture/VB/IB creation and locking still needed for host mirroring
7. **Draw patches stayed** — DrawVertices etc. still needed to trigger host GPU rendering

### Applicability to Our Codebase (dx11 branch)

**Our advantages over the reference branch:**
- D3D11 backend with interpreting shaders (VS interpreter, RC interpreter) that already read PGRAPH
- Page-tracked 64 MiB mirror for vertex/index data (no per-resource host VB/IB creation needed)
- VS constants already read from PGRAPH `pg->vsh_constants[]`
- Textures already have PGRAPH fallback via `NV_PGRAPH_TEXOFFSET`
- Pixel shader state already bridged HLE→PGRAPH→GPU via `pg->regs[]` upload
- Fog already partially reads from PGRAPH (color + _ABS mode)

**What we still need from the reference branch approach:**
- Remove individual `SetRenderState_*` patches (let Xbox code write to `D3D__RenderState[]`)
- Remove `SetTextureState_*` patches (let Xbox code write to `D3D__TextureState[]`)
- Unpatch `SetViewport` (read from PGRAPH or Xbox kernel memory)
- Unpatch `SetPixelShader` + `SetPixelShaderConstant` (already bridged to PGRAPH)
- Unpatch `SetStreamSource` (PGRAPH path already reads vertex_attributes directly)
- Consider unpatching `SetVertexShader`/`SelectVertexShader`/`LoadVertexShader` (most complex)
- Consider unpatching draw calls eventually (let push buffer draws handle everything)

**What we should NOT copy from the reference branch:**
- The D3D9 rendering code (we use D3D11)
- The OpenGL LLE rendering code (we use interpreting shaders)
- The HLE API enum / push buffer token injection system (we use `pfifo_submit_pushbuffer`)
- The vertex buffer converter code (we use shader-based vertex fetch)
