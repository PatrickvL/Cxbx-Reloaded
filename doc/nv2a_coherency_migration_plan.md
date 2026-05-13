# NV2A Coherency: Migration Plan — Current → Proposed Design

This plan describes how to safely evolve the memory coherency system from the current
implementation to the proposed unified design, minimizing regression risk at each step.

---

## Current State Summary

| Component | Current | Target |
|-----------|---------|--------|
| GPU buffer | `s_pMirrorBuf` — DEFAULT, SRV+UAV, 64 MiB RAM | Combined `s_pGpuMem` — DEFAULT, SRV+UAV, 64 MiB + MMIO |
| Upload method | `UpdateSubresource` with `D3D11_BOX` | `UpdateSubresource` with `D3D11_BOX` |
| MMIO storage | C struct `NV2AState` with `pg->regs[]` | Flat 16 MiB `s_pNV2AMMIO` allocation |
| PGRAPH→GPU | Separate `StructuredBuffer<uint>` at t12 | Appended to combined buffer, same t0 binding |
| PFB→GPU | Not uploaded | Appended to combined buffer |
| Deswizzle CS source | Per-texture staging buffer | Reads directly from combined buffer |
| Bitmaps | `volatile uint32_t[512]`, `alignas(64)`, interlocked | `volatile uint32_t[512]`, `alignas(64)`, interlocked |
| Staging textures | Cached per-RT | Cached per-RT |

---

## Migration Phases

### Phase 1: Cache Staging Textures Per-RT ✅ COMPLETE

**Goal:** Eliminate per-readback `CreateTexture2D` allocation in the VEH handler.

**Steps:**
1. Add an `ID3D11Texture2D* pStagingTex` field to `RegisteredRT`.
2. In `CxbxPageTrackerRegisterRT`, create the staging texture alongside the RT registration (same dimensions, `USAGE_STAGING`, `CPU_ACCESS_READ`).
3. In the VEH handler (`CxbxPageTrackerHandleFault`) and `CxbxPageTrackerFlushGPUDirtyToMirror`, use the cached staging texture instead of creating one on the fly.
4. Release the staging texture in `CxbxPageTrackerShutdown` and when an RT is evicted from the registered list.

**Result:** +27% FPS on DisplacementMap (69 → 88 FPS). No regressions observed.

**Risk:** Low — strictly performance optimization; same data flow.

---

### Phase 2: Upgrade Bitmaps to Thread-Safe Atomic Access ✅ COMPLETE

**Goal:** Prepare for multi-threaded bitmap access (required before moving GPU-dirty marking off the puller thread).

**Steps:**
1. Replace `uint32_t s_GpuDirtyBitmap[512]` with `alignas(64) volatile uint32_t s_GpuDirtyBitmap[512]`.
2. Replace `uint32_t s_TextureDirtyBitmap[512]` with `alignas(64) volatile uint32_t s_TextureDirtyBitmap[512]`.
3. Add `SetBitAtomic`/`ClearBitAtomic`/`TestBitAtomic` helpers that cast to `volatile long*` at the `_InterlockedOr`/`_InterlockedAnd` call sites.
4. Keep `s_TiledCommittedBitmap` as non-volatile `uint32_t[512]` (single-threaded access).
5. Update `CxbxPageTrackerIsTextureDirty` and `CxbxPageTrackerClearTextureDirty` to use `uint32_t` masks with casts to `long` only at interlocked API calls.

**Result:** No performance regression (83 FPS on DisplacementMap, same as Phase 1). Thread-safe VEH access confirmed working.

**Design notes:**
- 32-bit granularity chosen over byte-granularity (`uint8_t[2048]`) because `_InterlockedOr8`/`_InterlockedAnd8` are not available on x86.
- Arrays are `volatile uint32_t` (natural unsigned type); casts to `volatile long*` are confined to the interlocked intrinsic call sites.
- `_BitScanForward` calls use explicit `(unsigned long)` casts on the `uint32_t` bits argument.

**Risk:** Low — semantically equivalent with added thread safety.

---

### Phase 3: Switch to DEFAULT + UpdateSubresource ✅ COMPLETE

**Goal:** Enable UAV binding on the mirror buffer (prerequisite for CS writing directly to it).

**Steps:**
1. Change `s_pMirrorBuf` from `D3D11_USAGE_DYNAMIC` to `D3D11_USAGE_DEFAULT`.
2. Add `D3D11_BIND_UNORDERED_ACCESS` to bind flags.
3. Remove `D3D11_CPU_ACCESS_WRITE`.
4. Replace all `Map`/`Unmap` calls in `CxbxPageTrackerFlushToGPU` with `UpdateSubresource` + `D3D11_BOX`:
   - Bulk path: `UpdateSubresource(s_pMirrorBuf, 0, nullptr, CONTIG_BASE, CONTIG_SIZE, 0)` (full buffer, no box).
   - Incremental path: one `UpdateSubresource` per coalesced page run with a `D3D11_BOX{byteStart, 0, 0, byteEnd, 1, 1}`.
5. Replace `Map`/`Unmap` in `CxbxPageTrackerFlushGPUDirtyToMirror` with a single `UpdateSubresource` + `D3D11_BOX` for the VB page range.
6. Keep `s_bFirstFlushOfFrame` gating (once-per-frame flush) but remove DISCARD/NO_OVERWRITE branching.
7. Create `RWByteAddressBuffer` UAV (`s_pMirrorUAV`) over `s_pMirrorBuf` for future CS use.

**Result:** DisplacementMap 79 FPS (vs 83 FPS Phase 2). Minor regression expected — UpdateSubresource has higher per-call overhead than Map+memcpy, but enables UAV binding and future CS path. No visual regressions.

**Design alignment:** Implementation now matches the proposed design in `nv2a_emulation_memory_coherency.md` §2.2: DEFAULT usage, SRV+UAV bind flags, UpdateSubresource with D3D11_BOX for incremental updates. The buffer is ready for Phase 4 (append PGRAPH block) without further buffer recreation.

**Risk:** Medium — DISCARD buffer-orphaning optimization is lost. Acceptable given the small FPS delta observed.

---

### Phase 4: Append PGRAPH Block to Mirror Buffer ✅ COMPLETE

**Goal:** Eliminate the separate `g_pD3D11PGRegsBuf` (t12) binding.

**Result:** PGRAPH registers (8 KB) are now appended to the combined mirror buffer at offset `GPU_PGRAPH_BASE = 0x04000000u`. The separate StructuredBuffer and its SRV have been removed. HLSL shaders use `ByteAddressBuffer g_PGRegs : register(t12)` loading from `GPU_PGRAPH_BASE + byteOff`. Upload uses `CxbxPageTrackerUploadPGRAPH()` via `UpdateSubresource` with a D3D11_BOX. PS JIT cache version bumped to v2 to invalidate stale disk-cached shaders compiled against the old StructuredBuffer declaration.

**Validated:** DisplacementMap, Trees, Fur, CubeMap, DolphinClassic, SphereMap — all render correctly.

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
