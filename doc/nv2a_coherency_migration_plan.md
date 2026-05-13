# NV2A Coherency: Migration Plan — Current → Proposed Design

This plan describes how to safely evolve the memory coherency system from the current
implementation to the proposed unified design, minimizing regression risk at each step.

---

## Current State Summary

| Component | Current | Target |
|-----------|---------|--------|
| GPU buffer | `s_pMirrorBuf` — DYNAMIC, SRV only, 64 MiB RAM | Combined `s_pGpuMem` — DEFAULT, SRV+UAV, 64 MiB + MMIO |
| Upload method | `Map/Unmap` (DISCARD/NO_OVERWRITE) | `UpdateSubresource` with `D3D11_BOX` |
| MMIO storage | C struct `NV2AState` with `pg->regs[]` | Flat 16 MiB `s_pNV2AMMIO` allocation |
| PGRAPH→GPU | Separate `StructuredBuffer<uint>` at t12 | Appended to combined buffer, same t0 binding |
| PFB→GPU | Not uploaded | Appended to combined buffer |
| Deswizzle CS source | Per-texture staging buffer | Reads directly from combined buffer |
| Bitmaps | `uint32_t[512]`, non-atomic | `uint8_t[2048]`, `alignas(64)`, interlocked |
| Staging textures | Created per-readback, released immediately | Cached per-RT |

---

## Migration Phases

### Phase 1: Cache Staging Textures Per-RT

**Goal:** Eliminate per-readback `CreateTexture2D` allocation in the VEH handler.

**Steps:**
1. Add an `ID3D11Texture2D* pStagingTex` field to `RegisteredRT`.
2. In `CxbxPageTrackerRegisterRT`, create the staging texture alongside the RT registration (same dimensions, `USAGE_STAGING`, `CPU_ACCESS_READ`).
3. In the VEH handler (`CxbxPageTrackerHandleFault`) and `CxbxPageTrackerFlushGPUDirtyToMirror`, use the cached staging texture instead of creating one on the fly.
4. Release the staging texture in `CxbxPageTrackerShutdown` and when an RT is evicted from the registered list.

**Validation:** Run 5+ titles with RT-as-VB usage (DisplacementMap, DOA3 shadows). Confirm no readback regressions.

**Risk:** Low — strictly performance optimization; same data flow.

---

### Phase 2: Upgrade Bitmaps to Atomic Byte-Granularity

**Goal:** Prepare for multi-threaded bitmap access (required before moving GPU-dirty marking off the puller thread).

**Steps:**
1. Replace `uint32_t s_GpuDirtyBitmap[512]` with `alignas(64) uint8_t s_GpuDirtyBitmap[2048]`.
2. Replace `uint32_t s_TextureDirtyBitmap[512]` with `alignas(64) uint8_t s_TextureDirtyBitmap[2048]`.
3. Change `SetBit`/`ClearBit`/`TestBit` to use `_InterlockedOr8`/`_InterlockedAnd8` for the coherency bitmaps.
4. Keep `s_TiledCommittedBitmap` and `s_IdentityAllocBitmap` as `uint32_t[512]` (single-threaded access).
5. Update `CxbxPageTrackerIsTextureDirty` and `CxbxPageTrackerClearTextureDirty` to work with byte arrays.

**Validation:** Bit-level unit tests. Full title regression pass — bitmap behavior must be identical.

**Risk:** Low — semantically equivalent with added thread safety.

---

### Phase 3: Switch to DEFAULT + UpdateSubresource

**Goal:** Enable UAV binding on the mirror buffer (prerequisite for CS writing directly to it).

**Steps:**
1. Change `s_pMirrorBuf` from `D3D11_USAGE_DYNAMIC` to `D3D11_USAGE_DEFAULT`.
2. Add `D3D11_BIND_UNORDERED_ACCESS` to bind flags.
3. Remove `D3D11_CPU_ACCESS_WRITE`.
4. Replace all `Map`/`Unmap` calls in `CxbxPageTrackerFlushToGPU` with `UpdateSubresource` + `D3D11_BOX`:
   - Bulk path: `UpdateSubresource(s_pMirrorBuf, 0, nullptr, CONTIG_BASE, CONTIG_SIZE, 0)` (full buffer, no box).
   - Incremental path: one `UpdateSubresource` per coalesced page run with a `D3D11_BOX{byteStart, 0, 0, byteEnd, 1, 1}`.
5. Replace `Map`/`Unmap` in `CxbxPageTrackerFlushGPUDirtyToMirror` with per-page `UpdateSubresource`.
6. Remove `s_bFirstFlushOfFrame` gating logic (no DISCARD/NO_OVERWRITE distinction needed with DEFAULT).
7. Create `RWByteAddressBuffer` UAV over `s_pMirrorBuf` for future CS use.

**Validation:**
- Performance benchmark: compare frame times before/after. UpdateSubresource may be slightly slower for the bulk path but enables future gains.
- Visual regression: run 10+ titles, screenshot comparison.
- Verify typed SRV views (SNORM16x2, UNORM8x4) still work with DEFAULT usage.

**Risk:** Medium — changes the GPU upload path entirely. The DISCARD optimization (buffer orphaning) is lost; if titles with heavy CPU→GPU traffic regress, consider a hybrid approach with a staging buffer intermediary.

**Rollback:** Keep the DYNAMIC path behind a `#ifdef` until confidence is high.

---

### Phase 4: Append PGRAPH Block to Mirror Buffer

**Goal:** Eliminate the separate `g_pD3D11PGRegsBuf` (t12) binding.

**Steps:**
1. Increase `s_pMirrorBuf` size from 64 MiB to 64 MiB + 8 KB (PGRAPH).
2. Define `GPU_PGRAPH_BASE = 0x04000000u` offset constant.
3. After each `pg->regs_generation` change, call `UpdateSubresource` with a box at `[GPU_PGRAPH_BASE, GPU_PGRAPH_BASE + 8KB)` to upload `pg->regs[]`.
4. Create HLSL header `CxbxGpuMemAccess.hlsli` with:
   ```hlsl
   ByteAddressBuffer g_GpuMem : register(t0);
   #define PG_UINT_NEW(reg) g_GpuMem.Load(0x04000000u + (reg))
   ```
5. Dual-path shaders: `#ifdef USE_COMBINED_BUFFER` selects between `g_GpuMem.Load(PGRAPH_BASE + reg)` and the old `g_PGRegs[reg >> 2]`. Both paths must produce identical results.
6. Once validated, remove `g_pD3D11PGRegsBuf`, its SRV, and the t12 binding.
7. Free the t12 slot.

**Validation:**
- Pixel-perfect comparison: RC interpreter output with old t12 path vs. new combined-buffer path.
- PS JIT output must be unchanged (it bakes register offsets as constants).

**Risk:** Medium — shader changes touch every draw. The dual-path `#ifdef` provides safe rollback.

---

### Phase 5: Append PFB Block

**Goal:** Make tile configuration accessible from shaders.

**Steps:**
1. Increase buffer to 64 MiB + 8 KB (PGRAPH) + 4 KB (PFB).
2. Define `GPU_PFB_BASE = 0x04002000u`.
3. Upload `d->pfb.regs[]` on tile register writes.
4. Add `PFBLoad(reg)` macro to `CxbxGpuMemAccess.hlsli`.
5. Initially, no shader code reads PFB — this phase just makes it available for future use (e.g., tiling-aware texture fetch, compute-based RT readback with untile).

**Validation:** Buffer size increase only; no shader behavior change.

**Risk:** Very low — additive only.

---

### Phase 6: Deswizzle CS Reads from Combined Buffer

**Goal:** Eliminate per-texture staging buffer allocation for deswizzle CS.

**Steps:**
1. Modify `CxbxUnswizzleCS.hlsl` to accept a `SurfaceBase` cbuffer parameter and read from `g_GpuMem` at `SurfaceBase + swizzledOffset` instead of a separate `g_SrcBuffer`.
2. Remove the per-texture staging buffer upload path in `HostResourceCreate.cpp`.
3. The deswizzle CS now reads directly from the mirror buffer's RAM region.
4. Ensure the mirror buffer is up-to-date before the CS dispatch (it already is, since `FlushToGPU` runs before texture bind).

**Validation:**
- Texture correctness: compare deswizzled output texel-by-texel against CPU reference.
- Performance: one fewer `CreateBuffer`+`UpdateSubresource` per dirty texture.

**Risk:** Medium — changes the deswizzle data source. If the mirror buffer's region is stale (race with flush), textures would show corruption. Gate behind a feature flag initially.

---

### Phase 7: Flat MMIO Allocation (Optional, Longer-Term)

**Goal:** Replace `NV2AState` struct with a flat 16 MiB `s_pNV2AMMIO` allocation.

**Steps:**
1. Allocate 16 MiB `VirtualAlloc` with `PAGE_READWRITE`.
2. Redirect all `pg->regs[off >> 2] = val` writes to `*(uint32_t*)(s_pNV2AMMIO + 0x400000 + off) = val`.
3. Redirect all `d->pfb.regs[off >> 2]` to `*(uint32_t*)(s_pNV2AMMIO + 0x100000 + off)`.
4. Upload entire blocks via `memcpy` from the flat allocation to the GPU buffer.
5. Remove the `NV2AState` struct `regs[]` arrays (keep `program_data[]` and `vsh_constants[]` as separate sub-state).

**Validation:** Full regression suite — every register read/write must match.

**Risk:** High — touches every NV2A engine handler. Requires extensive testing across many titles. Should be deferred until all other phases are stable.

**Note:** This phase is optional. The current struct approach works fine with `memcpy(mapped.pData, pg->regs, 8192)`. The flat allocation mainly benefits code simplicity and documentation alignment, not runtime performance.

---

## Phase Dependencies

```
Phase 1 (staging cache)        ─── independent
Phase 2 (atomic bitmaps)       ─── independent
Phase 3 (DEFAULT + UpdateSubresource) ── Phase 4 (PGRAPH append)
                                              │
                                              ├─ Phase 5 (PFB append)
                                              │
                                              └─ Phase 6 (CS from combined buf)

Phase 7 (flat MMIO) is independent, can happen at any point after Phase 4.
Phases 1 and 2 are independent of each other and can be done in parallel.
```

---

## Risk Mitigation Strategy

1. **Feature flags:** Each phase introduces changes behind `#ifdef CXBX_COMBINED_GPU_BUFFER` / `#ifdef CXBX_ATOMIC_BITMAPS` / etc. Both paths coexist until the new path is validated.

2. **Pixel-perfect screenshots:** Capture reference screenshots from 10+ titles before each phase. Compare after. Any pixel difference triggers investigation.

3. **Performance counters:** Track per-frame metrics:
   - Number of `UpdateSubresource` calls (Phase 3)
   - Number of dirty pages flushed
   - VEH readback count per frame
   - Staging texture allocations per frame (should drop to 0 after Phase 1)

4. **Incremental commits:** Each phase is a single logical commit (may be split into sub-commits for review). Bisection must be possible per-phase.

5. **No behavioral coupling:** Phases are designed so that reverting one does not break others. The combined buffer phases (3-6) have a hard dependency chain, but within that chain each step is independently revertable to the prior state.

---

## Success Criteria

- After Phase 3: Mirror buffer is DEFAULT+UAV, typed SRVs work, no visual regression.
- After Phase 4: `g_PGRegs` t12 binding removed, all combiner/VS state reads go through t0.
- After Phase 6: Per-texture staging buffers eliminated, deswizzle reads from t0.
- After all phases: Single `ByteAddressBuffer` at t0 provides all memory access (RAM + PGRAPH + PFB). Staging texture allocation is zero per frame. Bitmap operations are lock-free.
