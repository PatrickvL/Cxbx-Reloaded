# NV2A Resource Management — Complete Lifecycle Report

## Overview

This document describes all GPU resource types managed by the Cxbx-Reloaded DX11 backend,
their lifecycle from Xbox RAM through PGRAPH registers to host GPU objects, synchronization
requirements, current implementation status, and recommended improvements.

The NV2A GPU has no dedicated VRAM — it shares the Xbox's 64 MiB unified memory.
Every GPU resource (textures, render targets, vertex buffers, palettes, programs) is
stored in this contiguous region at physical addresses `0x00000000–0x03FFFFFF`
(mapped to virtual `0x80000000–0x83FFFFFF` on Xbox, and the same on our host emulation).

---

## 1. Resource Types

### 1.1 Textures (Sampled Images)

**PGRAPH Registers (per stage, 4 stages):**
| Register | Offset | Content |
|----------|--------|---------|
| TEXOFFSET0-3 | 0x1A24+i*4 | Physical VRAM byte offset |
| TEXFMT0-3 | 0x1A04+i*4 | Color format, dimensionality, mip count, cubemap flag, border |
| TEXCTL0-3 | 0x19CC+i*4 | Enable (bit 30), alpha kill, max aniso, LOD clamp |
| TEXCTL1-3 | 0x19DC+i*4 | Image pitch (linear textures only) |
| TEXFILTER0-3 | 0x19F4+i*4 | Min/mag filter, LOD bias, convolution kernel, color sign |
| TEXIMAGERECT0-3 | 0x1A14+i*4 | Width/height in pixels (linear/NPOT textures) |
| TEXADDRESS0-3 | 0x19BC+i*4 | U/V/W wrap modes (clamp, repeat, mirror, border) |
| TEXPALETTE0-3 | 0x1A34+i*4 | Palette DMA context, length (256/128/64/32), offset |

**Data Layout in Xbox RAM:**
- **Swizzled**: Morton Z-curve order (interleaved X/Y/Z bits). Always power-of-two dimensions.
  Dimensions decoded from TEXFMT log2 USIZE/VSIZE/PSIZE fields.
- **Linear**: Row-major, arbitrary dimensions. Pitch from TEXCTL1; dimensions from TEXIMAGERECT.
- **Compressed (DXT1/3/5)**: 4×4 block tiles, row pitch = blocks_per_row × block_size.
- **Cubemap**: 6 faces × mip chain, packed sequentially. Face size aligned to `max(pitch, 256)`.
- **Volume (3D)**: Depth slices × mip chain. Each mip halves depth.
- **Mip Tail**: NV2A packs mip levels smaller than 1 block into a single "tail" region at
  hardware-specific offsets. Not the same as naive `offset += mipSlicePitch` advancement.

**Formats** (subset — 43 total in FormatInfos[]):
| Xbox Format | BPP | Type | Host DXGI | Notes |
|-------------|-----|------|-----------|-------|
| X_D3DFMT_A8R8G8B8 (0x06) | 32 | Swizzled | B8G8R8A8_UNORM | Direct, no conversion |
| X_D3DFMT_DXT1 (0x0C) | 4 | Compressed | BC1_UNORM | Block-compressed |
| X_D3DFMT_P8 (0x0B) | 8 | Palettized | R8G8B8A8_UNORM | Palette expand required |
| X_D3DFMT_YUY2 (0x24) | 16 | Linear | EMUFMT_YUY2 | **BUG: no conversion, SRV creation fails** |
| X_D3DFMT_LIN_A8R8G8B8 (0x12) | 32 | Linear | B8G8R8A8_UNORM | Direct |
| X_D3DFMT_R8G8B8A8 (0x3C) | 32 | Swizzled | R8G8B8A8_UNORM | Channel swap needed |

**When synced to GPU:**
- On first bind (texture stage enabled + TEXOFFSET non-zero). A stage is considered
  enabled if TEXCTL0 enable bit is set OR SHADERPROG specifies a non-NONE mode for
  that stage (commit cad5d46f8 — needed for point sprites using stage 3).
- When `CxbxPageTrackerIsTextureDirty(offset, size)` returns true (CPU wrote new data).
- When format/dimensions change (cache key mismatch → re-create).

**When synced back to CPU:** Never (textures are read-only from GPU perspective).

**Current Upload Path:**
```
CxbxUpdateHostTextures() → GetHostBaseTextureWithFormat() → [cache hit?]
  → YES: return cached (BUG: no dirty check on fast path)
  → NO: EmuVerifyResourceIsRegistered() → HostResourceRequiresUpdate()
    → CxbxPageTrackerIsTextureDirty() → CreateHostResource()
      → UploadPixelContainerMips() → [CS swizzle | CS palette | CPU fallback]
        → CxbxPageTrackerClearTextureDirty()
```

---

### 1.2 Render Targets (Color + Depth)

**PGRAPH State (via `NV2ASurfaceState`):**
| Field | Content |
|-------|---------|
| colorOffset | Physical VRAM offset of color RT |
| colorPitch | Row pitch in bytes |
| zetaOffset | Physical VRAM offset of depth/stencil RT |
| zetaPitch | Row pitch in bytes |
| clipWidth/clipHeight | Render region dimensions |
| colorFormat | NV097 color format enum |
| zetaFormat | NV097 zeta format enum |
| antiAliasing | 1x / 2xH / 4x |

`NV2ASurfaceState` is decoded directly from PGRAPH registers (no side-struct).
The old `SurfaceShape` struct was removed in commit 6010ff660.

**Data Layout:** Render targets are linear-pitch surfaces. The NV2A writes directly to VRAM
at `surface_color.offset`. Multiple RTs can exist at different offsets (shadow maps, cubemap
faces, reflections).

**When synced to GPU:** Never uploaded from CPU — GPU creates and renders to them directly.
The host D3D11 RT is created on first bind, matched by `{offset, format, width, height}` key.

**When synced back to CPU:**
- When Xbox CPU reads VRAM pages that an RT was rendered to (e.g., CPU-based post-processing,
  reading shadow map results, video capture).
- **IMPLEMENTED** (commit 5bd6999d4): VEH fault handler performs staged readback via
  `TryEnterCriticalSection` → staging texture → `CopyResource` → `Map` → row-by-row
  `memcpy` back to Xbox RAM. Entire RT is read back at once to amortize GPU stall.

**RT-as-Texture Path:** When TEXOFFSET matches a known RT offset, the host RT texture is
bound directly as SRV (no CPU round-trip needed — this is correct and efficient).
Implemented via `CxbxLookupPgraphRTByOffset()` in `HostSync.cpp:213-230`.

**RT Cache:**
```cpp
struct PgraphRTKey { xbox::addr_xt offset; DXGI_FORMAT format; UINT width, height; };
std::unordered_map<PgraphRTKey, ComPtr<ID3D11Texture2D>> g_PgraphRTCache;
```
- No eviction policy (grows unbounded).
- Cleared only on device reset (`CxbxResetPgraphSurfaceTracking()`).

---

### 1.3 Front/Back Buffers (Presentation)

**PGRAPH State:** The backbuffer is identified by comparing RT dimensions against
`g_EmuCDPD.HostPresentationParameters.BackBuffer{Width,Height}`.

**Data Flow:**
```
Game renders to RT at offset X (dimensions match presentation params)
  → CxbxD3D11UpdateRenderTargetFromPGRAPH() recognizes it as backbuffer
  → g_PgraphBackBufferOffsets.insert(X), g_pHostPgraphBackBuffer = pHostRT
  → NV097_FLIP_STALL triggers D3D11_flip_stall()
  → Blit g_pHostPgraphBackBuffer → swap chain backbuffer
  → g_pSwapChain->Present()
  → Re-acquire swap chain buffer (FLIP_DISCARD rotates buffers)
```

**Multi-buffer tracking** (PR #55 fix): Games with double/triple buffering cycle between
multiple RT offsets for the backbuffer. All matching offsets are accumulated in
`g_PgraphBackBufferOffsets` (unordered_set), and `g_pHostPgraphBackBuffer` always
points to the most recently rendered one.

**Auto-present fallback:** For raw pushbuffer games that never call FLIP_STALL,
the puller thread auto-presents on VBlank if `!g_pgraph_explicit_flip_stall_seen`
and `draw_dirty` is set.  During FMV (overlay active, no 3D draws), the DPC
thread wakes the puller so it composites and presents the PVIDEO overlay.

**When synced to CPU:** Same as RT readback — only if Xbox CPU reads the framebuffer
pages (extremely rare; typically only for screenshots or save-game thumbnails).

---

### 1.4 Video Overlay (PVIDEO)

**NV2A Registers:**
| Register | Content |
|----------|---------|
| NV_PVIDEO_OFFSET0/1 | VRAM offset of overlay plane (YUV surface) |
| NV_PVIDEO_SIZE_IN | Source dimensions |
| NV_PVIDEO_POINT_OUT | Destination position on screen |
| NV_PVIDEO_FORMAT | Format (YUY2, UYVY, etc.) |

**Status:** Implemented. `D3D11_flip_stall` reads PVIDEO registers to composite
overlay during the backbuffer blit. `CxbxSyncTiledRangeToContiguous()` is called
before reading frame data to handle games that write decoded video frames via the
tiled mapping (`0xF0000000`). Fixed in commit 57a058155.

**When synced:** The overlay surface lives in Xbox RAM and is typically updated by the
CPU (video decoder writes YUV frames). Tiled pages are synced to contiguous memory
every frame the overlay is active, before the compositor reads the data.

---

### 1.5 Vertex Data

**PGRAPH State:**
```cpp
VertexAttribute vertex_attributes[16];  // Per-attribute: format, size, count, offset, stride
uint32_t draw_arrays_start/count[];     // NV097_DRAW_ARRAYS batches
uint32_t inline_elements[];             // 16-bit index array
uint32_t inline_array[];                // Packed vertex data
float*   inline_buffer[];               // Per-attribute float4 arrays (Begin/End path)
```

**Three Draw Modes:**
| Mode | Source | Populated By |
|------|--------|-------------|
| Draw Arrays | Mirror buffer at `vertex_attributes[i].offset + vertex * stride` | SET_VERTEX_DATA_ARRAY_OFFSET/FORMAT |
| Inline Array | `pg->inline_array[]` (packed all-attrs-per-vertex) | NV097_INLINE_ARRAY method |
| Inline Buffer | `vertex_attributes[i].inline_buffer[]` (float4 per attr per vert) | SET_VERTEX_DATA2F/4F/4UB + ARRAY_ELEMENT |

**When synced to GPU:**
- **Draw Arrays**: Vertex data lives in the 64 MiB mirror buffer (`s_pMirrorBuf`).
  The page tracker flushes CPU-written pages to this GPU buffer before each draw.
  The VS samples it via ByteAddressBuffer SRV at byte offsets.
- **Inline Array/Buffer**: Uploaded to a DYNAMIC D3D11 buffer per-draw, then discarded.

**When synced back to CPU:** Never (vertex data is GPU-read-only).

---

### 1.6 Vertex Shader Programs

**PGRAPH State:**
```cpp
uint32_t program_data[NV2A_MAX_TRANSFORM_PROGRAM_LENGTH][VSH_TOKEN_SIZE];  // 136 × 4 uint32
uint32_t vsh_constants[NV2A_VERTEXSHADER_CONSTANTS][4];  // 192 × float4
```

**When synced to GPU:**
- Programs: JIT-compiled to native HLSL vertex shaders and cached (by hash, with disk
  cache). Falls back to interpreter (StructuredBuffer SRV) only when JIT validation
  fails (e.g., context writes). Added in commits f218b69b9 and 09c536aa2.
- Constants: Uploaded as cbuffer b0 every draw (192 × float4 = 3072 bytes).

**When synced back to CPU:** Never.

---

### 1.7 Palette Data

**PGRAPH Registers:** `TEXPALETTE0-3` contains DMA context select, palette length, and
byte offset within the DMA region.

**Current Source:** PGRAPH registers directly. Read `TEXPALETTE0 + stage*4`, extract offset
bits, resolve against DMA context A/B base address via `NV2AResolvePaletteAddress()` →
physical palette pointer in Xbox RAM. Contains 256/128/64/32 ARGB entries per stage.
The old HLE-based `g_pXbox_Palette_Data[]` globals were removed in commit 02a7982ae.

**When synced to GPU:** On texture upload (palette expanded CPU-side or via CS).
**When synced back to CPU:** Never.

---

### 1.8 Register Combiner (Pixel Shader) State

**PGRAPH Registers:** 68 registers in range 0x1880–0x1948 (combiner factors, alpha/color
input/output per stage, final combiner EFG).

**When synced:** Combiner topology is JIT-compiled to native HLSL pixel shaders and cached
by hash (commit 97fdc5eaa). The cache skips hash+mutex+map when combiner state is unchanged
(commit beec465c1). Dynamic per-frame data (C0/C1 colors, fog, bump matrices) is still
uploaded as part of the PGRAPH register SRV snapshot every draw. The interpreter path
remains as a fallback when JIT compilation fails.

---

## 2. Synchronization Architecture

### 2.1 CPU → GPU (Page Tracker)

```
┌─────────────────────────────────────────────────────────────┐
│                    Xbox RAM (64 MiB)                          │
│  0x80000000 ─ 0x83FFFFFF  (MEM_WRITE_WATCH enabled)         │
└──────────────────────────┬──────────────────────────────────┘
                           │ CPU writes (game logic, streaming, DMA)
                           ▼
┌──────────────────────────────────────────────────────────────┐
│              GetWriteWatch() on frame boundary                │
│  Returns list of dirty 4 KB pages since last reset           │
└──────────────────────────┬───────────────────────────────────┘
                           │
              ┌────────────┼────────────┐
              ▼            ▼            ▼
     s_TextureDirtyBitmap  │   s_pMirrorBuf (GPU)
     (texture staleness)   │   (vertex data source)
                           │
                           ▼
              CxbxPageTrackerFlushToGPU():
              - If >25% dirty + first-of-frame: MAP_WRITE_DISCARD + full memcpy
              - Else: MAP_WRITE_NO_OVERWRITE + coalesced page runs
```

### 2.2 GPU → CPU (Render Target Readback)

```
┌──────────────────────────────────────────────────────────┐
│  D3D11 renders to host RT at VRAM offset X               │
│  CxbxPageTrackerMarkGPUDirty(X, RT_size)                 │
│  Pages in s_GpuDirtyBitmap set, VirtualProtect(NOACCESS) │
└──────────────────────────┬───────────────────────────────┘
                           │ CPU reads address in RT range
                           ▼
┌──────────────────────────────────────────────────────────┐
│  Page fault → VEH → CxbxPageTrackerHandleFault()         │
│  Check s_GpuDirtyBitmap[pageIdx]:                        │
│    → TryEnterCriticalSection (serialize with puller)     │
│    → Find RegisteredRT covering the faulted page         │
│    → Restore PAGE_READWRITE on entire RT first           │
│    → CopyResource → staging → Map → memcpy to Xbox RAM  │
│    → Clear GPU-dirty bits for all RT pages               │
│    → LeaveCriticalSection                                │
└──────────────────────────────────────────────────────────┘
```

### 2.3 RT-as-Texture (GPU → GPU, no CPU round-trip)

```
┌──────────────────────────────────────────┐
│  Game renders to RT at offset X          │
│  g_PgraphRTCache[{X,fmt,w,h}] = pTex2D  │
└──────────────────────┬───────────────────┘
                       │ Later draw samples TEXOFFSET = X
                       ▼
┌──────────────────────────────────────────┐
│  CxbxLookupPgraphRTByOffset(X)          │
│  Returns existing ID3D11Texture2D*       │
│  CreateShaderResourceView → bind as SRV  │
│  CxbxInvalidatePgraphRTBinding() if      │
│  same texture was bound as RTV           │
└──────────────────────────────────────────┘
```

This is the correct fast path. No CPU round-trip needed because the same host
texture object serves as both RTV and SRV (just not simultaneously — the RTV
must be unbound before SRV binding, per D3D11 rules).

---

## 3. Confirmed Bugs (Ranked by Impact)

### CRITICAL

**#1 — Texture cache fast-path bypasses dirty check** ✅ FIXED (commit 7c8b275c5)
- **Where:** `GetHostBaseTextureWithFormat()` ([HostResource.cpp:441-447](src/core/hle/D3D8/Rendering/HostResource.cpp#L441))
- **What:** Cache hit returned immediately without calling `HostResourceRequiresUpdate()`.
  CPU writes detected by page tracker were never processed for cached textures.
- **Impact:** Any dynamically-updated texture (streamed content, animated textures, video
  frames, font glyphs) showed stale content after initial upload.
- **Resolution:** Added `else if` branch on cache-hit path that calls
  `CxbxPageTrackerIsTextureDirty()` and triggers `EmuVerifyResourceIsRegistered()` if
  dirty. Render targets are skipped (they're GPU-managed).
- **Risk:** Low — adds one bitmap scan (few DWORD reads) per texture per draw.

**#2 — GPU RT readback not implemented** ✅ FIXED (commit 5bd6999d4)
- **Where:** `CxbxPageTrackerHandleFault()` ([Backend_D3D11_PageTracker.cpp:191-270](src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_PageTracker.cpp#L191))
- **What:** When CPU reads a page previously rendered to by the GPU, the fault handler
  just cleared the dirty bit and restored access — but didn't copy GPU data back.
  The CPU read stale/zero data.
- **Impact:** Shadow maps read by CPU, environment capture, CPU-based post-processing,
  video capture/thumbnails all read garbage.
- **Resolution:** Implemented full staged readback in the VEH fault handler:
  - `RegisteredRT` metadata table tracks active RTs (offset, pitch, size, host texture)
  - On CPU read fault: find covering RT → create staging texture → CopyResource → Map →
    memcpy back to Xbox RAM row-by-row (handling host/Xbox pitch differences)
  - Entire RT is read back at once (amortizes GPU stall across all pages in the RT)
  - Pages are restored to PAGE_READWRITE BEFORE memcpy (prevents nested faults)
  - Thread safety: `TryEnterCriticalSection` serializes with puller thread — if puller
    holds the lock, readback is skipped (graceful degradation to stale data)
  - CPU writes to GPU-dirty pages skip readback (CPU is overwriting)
  - Upscaled RTs skip readback (would need a downscale resolve pass)
- **Design note:** "Once a render target lives in host GPU memory, it can theoretically
  be used as texture without having to download it back into Xbox RAM. Only when Xbox CPU
  accesses such regions, it would need a copy back." — The RT-as-texture SRV path already
  handles the GPU→GPU case correctly. Readback is only needed for genuine CPU reads.

### HIGH

**#3 — Page tracker once-per-frame guard** ❌ NOT A BUG (validated)
- **Where:** `CxbxPageTrackerFlushToGPU()` ([Backend_D3D11_PageTracker.cpp:388-397](src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_PageTracker.cpp#L388))
- **Analysis:** The once-per-frame guard was introduced in commit 7cfdadb07 to eliminate
  a 280× performance regression (GetWriteWatch ~74μs × 20k+ draws/frame = 1.5s overhead).
  Removing it was attempted and reverted — the guard is correct because:
  1. The Xbox CPU completes ALL memory writes before kicking the pushbuffer
  2. The first-of-frame flush populates `s_TextureDirtyBitmap` with all pages written
     since last frame
  3. Bug #1 fix (dirty check on cache hit) uses `IsTextureDirty()` — a cheap bitmap
     bit-test, not a syscall — which works correctly on every draw
  4. Mid-frame CPU writes are structurally impossible in the Xbox pushbuffer model
- **Conclusion:** The perceived "one-frame-late" issue was actually Bug #1 (cache bypass).
  With Bug #1 fixed, the once-per-frame flush is both correct and performant.

**#4 — Multi-buffer backbuffer tracking (was single-offset)**
- **Where:** `Backend_D3D11_State.cpp` backbuffer identification
- **What:** Previously used a single `g_PgraphBackBufferOffset` set on first match — any
  subsequent buffer from a double/triple-buffered game was not recognized.
- **Impact:** Alternate frames showed stale content (blitting wrong buffer).
- **Status:** Fixed in PR #55 (set-based tracking).

**#5 — Swap chain buffer not re-acquired after FLIP_DISCARD Present**
- **Where:** `CxbxPresent()` in `HostRender.cpp`
- **What:** `DXGI_SWAP_EFFECT_FLIP_DISCARD` rotates buffers — `GetBuffer(0)` returns a
  different texture after each Present. Stale RTV/surface pointers caused blits to
  already-queued-for-display surfaces.
- **Status:** Fixed in PR #55 (re-acquire after Present).

### MEDIUM

**#6 — Synthetic texture Common field cache key mismatch** ✅ FIXED (commit b7eac22f0)
- **Where:** `HostSync.cpp:222` — `synth.Common = X_D3DCOMMON_TYPE_TEXTURE | 1`
- **What:** Missing `X_D3DCOMMON_D3DCREATED` flag meant synthetic textures hashed to
  different cache keys than real textures at the same address.
- **Impact:** Duplicate host resources for same VRAM content; possible SRV mismatch.
- **Resolution:** ORed in `X_D3DCOMMON_D3DCREATED`:
  `synth.Common = X_D3DCOMMON_TYPE_TEXTURE | X_D3DCOMMON_D3DCREATED | 1;`

**#7 — Palette sourced from HLE patch, not PGRAPH** ✅ FIXED (commits 341ad00df, 02a7982ae)
- **Where:** `HostSync.cpp` — `CxbxUpdateHostTextures()` palette resolution path
- **What:** Palette pointer required `D3DDevice_SetPalette` HLE patch to be active.
  If disabled or called late, P8 textures zero-filled (black).
- **Resolution:** First added PGRAPH fallback (341ad00df), then fully removed the HLE
  `g_pXbox_Palette_Data[]` globals (02a7982ae). Palette address is now resolved entirely
  from PGRAPH `TEXPALETTE0 + stage*4` registers via DMA context A/B, matching xemu's
  hardware-accurate approach.

**#8 — YUY2/UYVY textures not converted** ✅ FIXED (commit a2cfd7374)
- **Where:** `XbConvert.cpp:115-116` — no warning string → `EmuXBFormatRequiresConversion` = false
- **What:** Mapped to custom EMUFMT_YUY2/UYVY which are not valid DXGI_FORMAT values for SRV.
  CreateShaderResourceView silently failed → texture unbound → black.
- **Impact:** FMV playback, video textures in any game using YUV surfaces.
- **Resolution:** Added warning string + `Texture` ResourceType to FormatInfos[] entries.
  This triggers `EmuXBFormatRequiresConversion()` → `ConvertD3DTextureToARGBBuffer()` which
  already had working row converters (`____YUY2ToARGBRow_C` / `____UYVYToARGBRow_C`).

**#9 — Linear texture pitch < 64 guard produces wrong dimensions** ✅ FIXED (commit 160e8b55c)
- **Where:** `HostSync.cpp:236` — `if (width > 0 && height > 0 && pitch >= 64)`
- **What:** If pitch < 64, `synth.Size` stayed 0, forcing swizzled dimension decode
  (log2 bits from TEXFMT) — which is wrong for a linear texture.
- **Impact:** Very narrow linear textures (1-wide L8, lookup tables) rendered as garbage.
- **Resolution:** Added `if (pitch < 64) pitch = 64;` clamp before Size encoding. The
  NV2A hardware enforces a minimum 64-byte pitch for linear textures.

**#10 — DXT mip pitch clamp after mip2dSize computation** ❌ NOT A BUG (validated)
- **Where:** `HostResourceUpload.cpp:290-294`
- **Analysis:** The clamp occurs at the END of the loop iteration (after mip2dSize is used),
  meaning it takes effect for the NEXT iteration. The mip2dSize calculation at the top of
  the loop uses `dwMipRowPitch` which was already clamped from the previous iteration's end.
  This is correct: the first mip level's pitch is always >= blockSize (set from the full
  texture row pitch), and subsequent mips inherit the clamped value.
- **Conclusion:** The pitch clamp ordering is correct as-is.

### LOW

**#11 — RT cache grows unbounded**
- No eviction. Long-running games with many RT switches leak host textures.
- Fix: LRU eviction when cache exceeds N entries (e.g., 64).

**#12 — No mip-tail packing awareness**
- NV2A packs sub-block mips into a "tail" with hardware-specific alignment.
- Our upload assumes naive `offset += mipSlicePitch` per level.
- Impact: Barely visible (only affects lowest mip levels of high-mipcount textures).

**#13 — RT usage flag not set for RT-capable textures**
- `GetHostBaseTextureWithFormat` passes `D3DUsage=0`. If a texture is also used as RT
  but not found in `g_PgraphRTCache`, it's created without RENDER_TARGET flag and can't
  serve as both RTV and SRV.

---

## 4. Ideal Resource Management Architecture

### Design Principles

1. **Single source of truth: PGRAPH registers + Xbox RAM.**
   No HLE side-maps, no `g_pXbox_*` globals. Everything derivable from PGRAPH state.

2. **Unified resource cache keyed by {VRAM offset, format, dimensions}.**
   No separate RT cache vs texture cache. One cache, one lookup.

3. **Lazy creation, eager invalidation.**
   Create host resources on first access. Invalidate when page tracker detects writes.

4. **GPU→GPU for RT-as-texture (zero-copy).**
   Same D3D11 texture serves as RTV and SRV. Only unbind RTV when binding SRV.

5. **GPU→CPU readback only for genuine CPU reads.**
   Page-fault driven. Amortize by reading back entire RT, not just one page.

6. **No per-draw hashing (unlike xemu).**
   Page tracker dirty bits are O(1) to check. Per-draw CRC hashing is too expensive
   for a real-time emulator targeting 60 fps.

### Proposed Unified Cache

```cpp
struct ResourceEntry {
    ComPtr<ID3D11Texture2D> pTexture;     // Host texture (can be RTV + SRV)
    ComPtr<ID3D11ShaderResourceView> pSRV;
    ComPtr<ID3D11RenderTargetView> pRTV;  // Non-null if RT-capable
    ComPtr<ID3D11DepthStencilView> pDSV;  // Non-null if depth
    DXGI_FORMAT format;
    UINT width, height, depth, mipLevels;
    uint32_t vramOffset;                  // Xbox physical offset
    uint32_t vramSize;                    // Total bytes in VRAM
    uint64_t lastUploadFrame;             // Frame counter of last CPU→GPU upload
    uint64_t lastGpuWriteFrame;           // Frame counter of last GPU render
    uint64_t lastAccessFrame;             // LRU tracking
    bool isRenderTarget;                  // GPU writes to this
    bool isCpuReadable;                   // Has been read back
};

// Single global: O(1) lookup by VRAM offset (most common case)
std::unordered_map<uint32_t /*vramOffset*/, ResourceEntry> g_ResourceCache;
```

### Lookup Flow (per draw, per texture stage)

```
uint32_t offset = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage*4)];
auto it = g_ResourceCache.find(offset);

if (it != g_ResourceCache.end()) {
    // Validate format/dimensions haven't changed
    if (FormatMatches(it->second, pg, stage)) {
        // Check dirty pages (skip for RTs — they're GPU-managed)
        if (!it->second.isRenderTarget &&
            CxbxPageTrackerIsTextureDirty(offset, it->second.vramSize)) {
            ReuploadTexture(it->second, pg, stage);
        }
        BindSRV(it->second.pSRV, stage);
        return;
    }
    // Format changed — destroy and re-create
    g_ResourceCache.erase(it);
}

// Create new
ResourceEntry entry = CreateResourceFromPGRAPH(pg, stage, offset);
g_ResourceCache[offset] = entry;
BindSRV(entry.pSRV, stage);
```

### RT Bind Flow

```
uint32_t colorOffset = pg->surface_color.offset;
auto it = g_ResourceCache.find(colorOffset);

if (it != g_ResourceCache.end() && it->second.isRenderTarget) {
    if (DimensionsMatch(it->second, pg)) {
        BindRTV(it->second.pRTV);
        return;
    }
    g_ResourceCache.erase(it);
}

ResourceEntry entry = CreateRTFromPGRAPH(pg, colorOffset);
entry.isRenderTarget = true;
g_ResourceCache[colorOffset] = entry;
BindRTV(entry.pRTV);
CxbxPageTrackerMarkGPUDirty(colorOffset, entry.vramSize);
```

### Benefits
- Eliminates the "texture at RT offset" special case — same entry serves both.
- Eliminates duplicate entries for same VRAM content (no Common-field key mismatch).
- Enables proper LRU eviction (single cache to manage).
- Makes dirty-check mandatory for all texture binds (fixes Bug #1).
- Simplifies RT-as-texture (just unbind RTV, bind SRV from same entry).

---

## 5. Phased Implementation Plan

### Phase 1: Fix Critical Bugs (no architecture change)

| Fix | File | Effort | Status |
|-----|------|--------|--------|
| Dirty check on texture cache fast path | HostResource.cpp | 10 lines | ✅ Done (7c8b275c5) |
| Once-per-frame flush guard | Backend_D3D11_PageTracker.cpp | — | ❌ Not a bug |
| YUY2/UYVY CPU conversion | XbConvert.cpp | 2 lines | ✅ Done (a2cfd7374) |
| Synthetic Common field fix | HostSync.cpp | 1 line | ✅ Done (b7eac22f0) |
| DXT mip pitch clamp ordering | HostResourceUpload.cpp | — | ❌ Not a bug |
| Linear pitch < 64 clamp | HostSync.cpp | 3 lines | ✅ Done (160e8b55c) |

### Phase 2: RT Readback + Palette Fallback

| Fix | File | Effort | Status |
|-----|------|--------|--------|
| GPU→CPU RT readback via VEH | Backend_D3D11_PageTracker.cpp + State.cpp | 180 lines | ✅ Done (5bd6999d4) |
| PGRAPH-based palette (HLE globals removed) | HostSync.cpp | 40 lines | ✅ Done (341ad00df, 02a7982ae) |
| RT usage flag from PGRAPH RT table | — | — | ❌ Not needed (two-cache arch handles it) |

### Phase 3: Unified Resource Cache (architecture change)

| Task | Files | Effort | Risk |
|------|-------|--------|------|
| Design ResourceEntry struct | New header | 50 lines | Low |
| Merge g_PgraphRTCache + ResourceCache | Backend_D3D11_State.cpp + HostResource.cpp | 200 lines | High |
| Single-lookup CxbxBindTexture() | HostSync.cpp | 100 lines | Medium |
| LRU eviction policy | New utility | 40 lines | Low |
| Remove dead code (g_TexturesByDataAddr, etc.) | Multiple | 150 lines deleted | Medium |

### Phase 4: Completeness

| Task | Files | Effort | Risk |
|------|-------|--------|------|
| Mip-tail packing awareness | HostResourceUpload.cpp | 30 lines | Medium |
| PVIDEO overlay composition | XbPushBuffer.cpp + new shader | 150 lines | Medium |
| Volume texture + 3D mip validation | HostResourceUpload.cpp | 20 lines | Low |

---

## 6. Comparison: Cxbx vs xemu

| Aspect | Cxbx-Reloaded (DX11) | xemu (OpenGL) |
|--------|----------------------|---------------|
| **Texture staleness** | Page tracker dirty bits (O(1) bitmap check, fixed) | Per-draw pixel data CRC hash (expensive but always correct) |
| **RT-as-texture** | Same host texture, SRV/RTV swap (correct, fast) | glFramebufferTexture2D → detach + bind as texture (similar) |
| **RT readback** | VEH fault → staging texture → CopyResource → memcpy (fixed) | GL FBO read (glReadPixels) |
| **Palette source** | PGRAPH TEXPALETTE register via DMA context (hardware-accurate, commit 02a7982ae) | PGRAPH TEXPALETTE register (hardware-accurate) |
| **YUV textures** | CPU YUY→ARGB conversion (fixed) | CPU YUV→RGB conversion |
| **Mip tail** | Naive offset advancement | Hardware-accurate packing |
| **Cache key** | {Common, Data, Format, Size, PaletteHash} | {TextureShape, data_hash, data_addr} |
| **Cache eviction** | None (grows forever) | GHashTable (auto via GLib) |
| **Vertex fetch** | GPU mirror buffer + VS ByteAddressBuffer sample | Traditional VBO upload per-draw |
| **Vertex shader** | JIT to HLSL (with interpreter fallback) | Static GLSL translation |
| **Surface format** | PGRAPH regs → DXGI_FORMAT mapping | PGRAPH regs → GL internal format |
| **DMA contexts** | Full RAMHT resolution (commits 49cd110f5, 513a2c2ea) | Full RAMHT resolution |
| **Performance model** | Zero-copy vertex fetch + page-tracked textures (fast) | Per-draw hash + VBO realloc (correct but slower) |

**Key insight:** Cxbx's architecture is fundamentally more performant than xemu's
(no per-draw hashing, zero-copy vertex fetch). The page-tracker model is sound —
with Bug #1 fixed, the dirty-check bitmap provides correct staleness detection on
every cache hit at negligible cost (a few DWORD reads per texture per draw).
Phases 1 and 2 are complete. Remaining work is architectural cleanup (Phase 3).

---

## 7. Summary

The NV2A DX11 backend manages 8 resource types through a page-tracker-based dirty
tracking system that is architecturally sound but has implementation gaps:

1. **Textures** — Working. Fast-path dirty check added (commit 7c8b275c5).
2. **Render Targets** — Working for GPU→GPU and GPU→CPU readback (commit 5bd6999d4).
3. **Front/Back Buffers** — Fixed (PR #55 multi-buffer + FLIP_DISCARD).
4. **Video Overlay** — Incomplete (patches disabled, composition stubbed).
5. **Vertex Data** — Working (mirror buffer + ByteAddressBuffer fetch).
6. **VS Programs/Constants** — Working (uploaded every draw).
7. **Palettes** — Working. HLE path + PGRAPH fallback (commit 341ad00df).
8. **RC/PS State** — Working (register SRV snapshot every draw).

Phases 1 and 2 are complete: Bugs #1, #2, #6, #7, #8, #9 are fixed; #3, #10, #13
were validated as non-issues. Remaining work is Phase 3 (unified resource cache
architecture) and Phase 4 (completeness).
