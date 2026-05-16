# DX11 Backend — Remaining Work

> Working document. Updated: 2026-05-16  
> Branch: `dx11` (HEAD: a7017c704)  
> Last audit: 2026-05-16 — verified each item against code

---

## 1. D3D11 Backend Optimizations

### 1.1 Draw Batching
**Status:** Failed first attempt — regressed clear ops and label/text rendering.  
**Problem:** Consecutive draws with identical pipeline state should be coalesced into a single `Draw()` call, but the initial implementation suppressed clear-only draw calls.  
**Next step:** Investigate why batched draws skip clears. Likely the NOP-check (no color write + no depth test + no stencil test) incorrectly classifies clear operations when batched.

### ~~1.2 AVX2 Bulk Bitmap Scan~~
**Status:** ✅ DONE.  
`CxbxPageTrackerIsTextureDirty()` now uses AVX2 `_mm256_testz_si256` to test 256 pages (1 MiB) per iteration when the CPU supports it. Runtime CPUID detection falls back to the scalar DWORD loop on older hardware. Bitmap is `alignas(64)` so aligned 256-bit loads are safe.

### ~~1.3 RT Cache Eviction~~
**Status:** ✅ DONE.  
Elastic LRU eviction in `CxbxPgraphRTCacheEvict()` — called each frame from `D3D11_flip_stall`. Cache grows freely up to 64 entries (high watermark); when exceeded, the oldest entries by last-access frame are evicted down to 32 (low watermark). Current backbuffer is pinned from eviction.

### ~~1.4 RT Cache + HLE Resource Eviction~~
**Status:** ✅ DONE.  
RT cache: replaced O(n) linear scan with offset-keyed hash map, multi-entry per offset (format/size variants), elastic LRU eviction (256/128 watermarks).  
HLE cache: added PruneResourceCache() with watermark-based LRU (2048/1024) — previously grew unbounded.  
Note: the two caches remain separate structures (incompatible key domains: VRAM offset vs Xbox D3DResource pointer). Original "merge" aspiration was impractical.

### ~~1.5 Mip Tail Packing~~
**Status:** ✅ NOT A BUG — verified correct against xemu.  
The `dwMipRowPitch` clamp to `blockSize` minimum ensures sub-4x4 compressed mips advance by exactly one block per level (matching xemu's `physical_width/4 * physical_height/4 * block_size`). `numRows = (h+3)/4` correctly gives 1 for heights ≤ 4. D3D11 `UpdateSubresource` receives the correct `RowPitch = blockSize` for BC minimum subresources.

### 1.6 Tiled RT Readback Correctness
**Status:** Known limitation.  
**Problem:** VEH readback path uses `CopyLinearWithPitch` without tile untransform. CPU reads of tiled RTs via the Contiguous window see Morton-encoded data.  
**Workaround:** GPU→GPU fast path (RT-as-texture) avoids readback entirely for the common case. Only affects explicit CPU reads (screenshots, save-game thumbnails).  
**Fix:** Apply `TiledToLinear` per pixel in the readback memcpy.

### 1.7 PFB/PVIDEO Shader Upload
**Status:** Buffer space reserved, upload functions exist, call sites commented out.  
**Trigger:** Enable when shaders need to read tile configuration or overlay state directly.  
**Files:** `Backend_D3D11_PageTracker.cpp` (lines 1020/1033), shader includes

### ~~1.8 PVIDEO Overlay Compositing~~
**Status:** ✅ DONE.  
`D3D11_flip_stall` reads PVIDEO registers directly (`NV_PVIDEO_OFFSET`, `NV_PVIDEO_FORMAT`, etc.) and composites the overlay onto the host backbuffer. `EnableOverlay`/`UpdateOverlay` patches are disabled.

---

## 2. EMUPATCH Removal (Phase 3 of LLE Plan)

### ~~2.1 Batch 1 — Pure State Patches~~ ✅ ALL DISABLED
All PATCH_ENTRY lines are commented out in `Patches.cpp`. Puller populates PGRAPH regs[] and interpreters read from there. Verified:
- `SetRenderState_Simple` (L239), `SetTexture` (L262), `SetVertexShader` (L276)
- `SetStreamSource` (L254), `SetViewport` (L296), `SetPixelShader` (L232)
- `SetVertexShaderConstant`, `SetTransform`, `SetScissors` — all disabled

**Remaining cleanup:** ~~Move dead implementations to unused/dead code file (per project principle).~~ ✅ DONE — `RaiseException` moved to `Xapi.cpp.unused-patches`; `XAudioCreateAdpcmFormat` already `#if 0` guarded.

### 2.2 Batch 2 — Resource Management Patches
**Status:** Needs verification per-patch (some may already be disabled).
| Patch | Notes |
|-------|-------|
| `D3DDevice_SetPalette` | Palette resolved from PGRAPH TEXPALETTE regs |
| `D3DDevice_SetVertexData*` | Inline vertex data via NV097 methods |
| `D3DDevice_LoadVertexShader` / `DeleteVertexShader` | XFPR from puller |
| `D3DDevice_SetGammaRamp` | DAC state |

### ~~2.3 Batch 3 — Sync Patches~~ ✅ DONE
All sync patches disabled AND implementations removed. Confirmed in `d3d-unpatch-report.md`:
| Patch | Status |
|-------|--------|
| `D3DDevice_BlockOnFence` / `InsertFence` / `IsFencePending` | ✅ Disabled + removed |
| `D3DDevice_BlockOnTime` | ✅ Disabled + removed (native semaphore works) |
| `D3DDevice_BeginVisibilityTest` / `EndVisibilityTest` / `GetVisibilityTestResult` | ✅ Disabled + removed (PGRAPH handles NV097) |
| `D3DResource_BlockUntilNotBusy` | ✅ Disabled + removed |

### ~~2.4 Batch 4 — Presentation Patches~~ ✅ DISABLED
| Patch | Status |
|-------|--------|
| `D3DDevice_Present` / `Swap` / `PersistDisplay` | Disabled (L188/L298) |
| `D3DDevice_EnableOverlay` / `UpdateOverlay` | Disabled |

---

## ~~3. NV2A GPU Synchronization~~ ✅ DONE

### ~~3.1 Semaphore/Fence Support~~ ✅ DONE
`NV097_BACK_END_WRITE_SEMAPHORE_RELEASE` (EmuNV2A_PGRAPH.cpp L1614) resolves DMA context, writes `parameter` to `semaphore_offset` via `stl_le_p`. Full implementation.

### ~~3.2 ZPass / Visibility Queries~~ ✅ DONE
`D3D11_zpass_begin/end/collect` implemented in `XbPushBuffer.cpp` (L225/251/270). `NV097_GET_REPORT` (EmuNV2A_PGRAPH.cpp L1342) collects query, writes timestamp + pixel count to report DMA. Wired to `g_pgraph_backend` function pointers.

---

## 4. Display Path (PCRTC/PVIDEO)

### ~~4.1 PCRTC Scan-Out~~ ✅ DONE
`D3D11_flip_stall` reads `d->pcrtc.start` and resolves the corresponding host RT via `CxbxLookupPgraphRTByOffset()`. Falls back to `g_pHostPgraphBackBuffer` when `pcrtc.start` is zero or not yet in the cache (early boot). Overlay scaling derives Xbox framebuffer dimensions from the resolved display surface.

### ~~4.2 PVIDEO Overlay~~ ✅ DONE
`D3D11_flip_stall` reads PVIDEO registers (`NV_PVIDEO_OFFSET`, `NV_PVIDEO_FORMAT`, etc.) directly and composites the overlay. `EnableOverlay`/`UpdateOverlay` patches are disabled.

---

## 5. GPU Tessellation (Patches/Higher-Order Surfaces)

**Status:** ✅ CPU-side tessellation DONE. GPU compute shader WIP (stashed).  
Full CPU-side FD (forward differencing) tessellation in `PatchDraw.cpp` — multi-attribute evaluation triggered on `SET_END_PATCH`, renders via `CxbxD3D11VertexFetchDraw`.  
**Remaining:** The stashed GPU CS path (`git stash push --staged -m "GPU tessellation CS work in progress"`) would move evaluation to a compute shader for performance. Not blocking correctness.  
**Test title:** XDK Patch sample.

---

## 6. Correctness Issues Under Investigation

### 6.1 Gauntlet Dark Legacy
**Status:** Active testing (most recent terminal runs).  
**Issues:** Unknown — need to characterize rendering artifacts from diagnostic output.

### 6.2 Vertex Fetch Diagnostics
**Status:** Diagnostic logging was used during testing (output files cleaned up). No permanent diagnostic code in tree.  
**Context:** Attribute fetch correctness for titles with complex vertex formats.

---

## 7. Code Organization (Completed)

- [x] Move `XbPixelShaderCompiler.cpp` → `Backend_D3D11_PixelShader.cpp`
- [x] Move `XbVertexShader.cpp` → `Backend_D3D11_VertexShader.cpp`
- [x] Split `RenderGlobals.h` → `Backend_D3D11.h`
- [x] Move D3D11 files into `Backend/` folder (Step 11.2)
- [x] Vulkan SDK detection in CMake (Step 11.4)

---

## 8. Documentation Debt

### 8.1 NV2A_LLE_Migration_Plan.md Status Update
**Problem:** "Current Architecture" diagram still shows old HLE path as active. Phases 1–2 are largely complete (puller-driven, interpreters read PGRAPH). No status markers beyond one ✅.  
**Action:** Add completion markers to Phase 1 (1.1–1.5) and Phase 2 (2.1–2.5) sections.

---

## 9. Long-Term (Post-DX11)

| Item | Phase | Notes |
|------|-------|-------|
| Vulkan backend | 4 | HLSL→SPIR-V via DXC; replace D3D11 device/swapchain |
| GPU compute pushbuffer | 5 | Process NV2A commands on host GPU |
| PFIFO batch-read | 1.5 | memcpy entire GET→PUT segment |

---

## Priority Order

1. **§6.1 Gauntlet correctness** — active debugging
2. **§1.1 Draw batching** — highest perf gain remaining
3. ~~**§2.1 Dead code cleanup**~~ ✅ DONE
4. ~~**§1.4 RT cache + HLE resource eviction**~~ ✅ DONE
5. **§5 GPU tessellation CS** — performance (CPU path works)
