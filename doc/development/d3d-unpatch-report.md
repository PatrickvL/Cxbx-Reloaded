# D3D Patch Removal Report

## Session Summary

Starting from commit `aefa8756a`, we systematically disabled D3D EMUPATCH entries
that are unnecessary under the NV2A-first architecture. Each batch was built, tested
against 5 XDK samples (Dolphin, Meshes, PushBuffer, Textures, Ripple), and committed.

### Core Principle

Xbox D3D API functions emit NV2A push buffer commands via PFIFO → PGRAPH. Our NV2A
emulation processes those commands regardless of their origin. Host resources are
purely derived from NV2A state — NOT from Xbox D3D objects. Therefore, any patch that
merely intercepts the Xbox D3D → NV2A push buffer flow is unnecessary.

## Commits (6 total, 8 files changed, +101/-71)

| Commit | Description | Entries Disabled |
|--------|-------------|-----------------|
| `2b672d542` | Disable SetRenderTarget and Reset patches; PGRAPH backbuffer tracking | 5 |
| `388c0ad4e` | Disable PersistDisplay and GetBackBuffer patches | 5 |
| `2827283d8` | Disable DestroyResource patches | 2 |
| `4125f8e04` | Disable BeginPush/EndPush patches | 3 |
| `36279d0bf` | Disable DrawRectPatch/DrawTriPatch patches | 2 |
| `4f9fedec4` | Disable SetSwapCallback and SetBackBufferScale patches | 2 |
| **Total** | | **19** |

## Test Results (post all batches, commit `4f9fedec4`)

| Sample | FPS | Status |
|--------|-----|--------|
| Dolphin | ~10.77 | Dolphin + seafloor rendering correctly |
| Meshes | ~59.98 | Tiger visible, lit, textured |
| PushBuffer | ~56.03 | Colored triangle + text |
| Textures | ~41.53 | Textured cylinder |
| Ripple | ~5.13 | Ripple wave |

All samples render correctly. No regressions observed.

---

## Remaining Active D3D PATCH_ENTRYs (25)

### KEEP — Must Remain Patched (9 entries)

**Direct3D_CreateDevice** (4 entries)
- Xbox native crashes (CMiniport access, etc.). Host D3D11 device must be created here.
- Future: move host init to NV2A device init, but keep patch to prevent Xbox crash.

**D3DDevice_Swap / Present** (3 entries)
- This IS the host present path: PFIFO flush → blit PGRAPH backbuffer → host window.
- Xbox native Swap writes PCRTC flip registers; we intercept to blit to DXGI.
- Frame limiter, VBlank counters live here.

**D3DDevice_CopyRects** (1 entry)
- Xbox native uses NV2A 2D blit engine (NV062/NV09F) which operates on VRAM.
- Host GPU textures are authoritative; VRAM copies are stale/empty.
- Dual-residency problem: native blit would copy stale VRAM, not host GPU content.
- KEEP until render-to-VRAM readback is implemented.

**D3DDevice_SetGammaRamp / GetGammaRamp** (2 entries — likely permanent)
- Xbox native writes PRAMDAC gamma LUT via MMIO.
- We use DXGI gamma control. PRAMDAC LUT format differs from DXGI gamma.
- High effort to convert; likely always patched.

### UNPATCH-SOON — Minor NV2A Work Needed (8 entries)

**D3DDevice_RunVertexStateShader** (2 entries)
- Xbox native uses NV097_LAUNCH_TRANSFORM_PROGRAM (0x1E90) — currently unhandled.
- **Needed**: Implement NV097_LAUNCH_TRANSFORM_PROGRAM handler in PGRAPH.
  This should invoke the fixed-function or programmable vertex shader on specified
  data and write results to output registers or context DMA.

**D3DDevice_InsertCallback** (1 entry)
- Xbox native writes NV097_NO_OPERATION(param≠0) to push buffer.
- PGRAPH already fires software interrupt on NOP with non-zero param.
- **Needed**: Verify interrupt delivery → kernel ISR → D3D callback dispatch works
  end-to-end. May already work — just needs testing.

**D3D_BlockOnTime** (2 entries)
- Xbox native uses semaphore, wait-for-idle, and nop methods pushed into the
  command buffer. Completion signals fire when DMA_PUT inline processing
  (pfifo_run_pusher) executes the pushed methods.
- **Status**: Unpatched. Works correctly with inline command processing on DMA_PUT writes.

**D3DDevice_BeginVisibilityTest / EndVisibilityTest / GetVisibilityTestResult** (3 entries)
- Xbox native uses NV097_SET_ZPASS_PIXEL_COUNT_ENABLE + NV097_GET_REPORT.
- These write results to Xbox-visible memory (report semaphore).
- **Needed**: Implement NV097_SET_ZPASS_PIXEL_COUNT_ENABLE and NV097_GET_REPORT
  using D3D11 occlusion queries, triggered from PGRAPH method handlers.
  Write results back to Xbox report memory so native GetVisibilityTestResult reads them.

### UNPATCH-LATER — NV2A Feature Implementation Needed (6 entries)

**D3DDevice_EnableOverlay / UpdateOverlay** (4 entries)
- Xbox native writes PVIDEO registers.
- Our Swap reads g_OverlayProxy (HLE global).
- **Needed**: Read PVIDEO registers (NV_PVIDEO_BUFFER, NV_PVIDEO_OFFSET,
  NV_PVIDEO_SIZE_IN/OUT, NV_PVIDEO_POINT_IN/OUT) in Swap path instead of
  HLE overlay proxy. Then unpatch.

**D3DDevice_BlockUntilVerticalBlank** (1 entry)
- Xbox native waits on PCRTC VBlank interrupt.
- **Needed**: Implement PCRTC interrupt emulation — periodic timer → NV_PCRTC_INTR
  at display refresh rate. This also benefits many other titles that busy-wait on VBlank.

**D3DDevice_GetDisplayFieldStatus** (1 entry)
- Xbox native reads PCRTC raster/field registers.
- **Needed**: Implement NV_PCRTC_RASTER register and interlace field detection.

### Summary Table — What Remains Needed to Unpatch More

| Blocked By | Entries | NV2A Feature |
|------------|---------|-------------|
| PGRAPH: NV097_LAUNCH_TRANSFORM_PROGRAM | 2 | Vertex shader launch method |
| PGRAPH: NV097_SET_ZPASS_PIXEL_COUNT_ENABLE + GET_REPORT | 3 | Occlusion query → report memory |
| PFIFO: accurate DMA_GET (remove fast-path hack) | 2 | Proper DMA progress tracking |
| PGRAPH: NV097_NO_OPERATION interrupt → kernel ISR | 1 | Interrupt delivery verification |
| PVIDEO: register reads in Swap path | 4 | Video overlay registers |
| PCRTC: VBlank interrupt | 1 | VBlank timer interrupt |
| PCRTC: raster/field registers | 1 | Display field status |
| PRAMDAC: gamma LUT → DXGI conversion | 2 | Gamma ramp mapping |
| Host: render-to-VRAM readback | 1 | Resolve dual-residency for CopyRects |
| Host: device creation infrastructure | 4 | Can't run Xbox native CreateDevice |

### Priority Order for Future Work

1. **InsertCallback** — May already work. Just test and disable. (1 entry)
2. **Visibility tests** — D3D11 occlusion queries from PGRAPH. (3 entries)
3. **RunVertexStateShader** — NV097_LAUNCH_TRANSFORM_PROGRAM. (2 entries)
4. **Overlay** — Read PVIDEO registers in Swap. (4 entries)
5. **BlockUntilVerticalBlank** — PCRTC VBlank interrupt. (1 entry)
6. **GetDisplayFieldStatus** — PCRTC raster registers. (1 entry)
