# NV2A Emulation: Full Memory Coherency and Rendering System

*Revised against Cxbx-Reloaded DX11 backend implementation, xemu, and Kelvin/NV20-family hardware documentation.*

## Table of Contents

1. [Physical Memory Architecture](#1-physical-memory-architecture)
2. [Host Memory Representation](#2-host-memory-representation)
   - 2.1 [s_pXboxRAM — Xbox RAM Host Mirror](#21-the-canonical-host-mirror--s_pxboxram)
   - 2.2 [s_pGpuMem — Combined XboxRAM + MMIO GPU Buffer](#22-the-gpu-side-mirror--s_pgpumem-combined-xboxram--mmio-blocks)
   - 2.3 [NV2A MMIO Backing Storage](#23-nv2a-mmio-backing-storage)
   - 2.4 [Uploading MMIO Backing to the Host GPU](#24-uploading-mmio-backing-to-the-host-gpu)
   - 2.5 [XFPR, XFCTX, and Derived Fields](#25-pgraph-sub-state--xfpr-xfctx-and-derived-fields)
   - 2.6 [Staging Buffer](#26-the-staging-buffer)
3. [Address Mirrors and Physical Resolution](#3-address-mirrors-and-physical-resolution)
4. [Dirty Tracking — Bit-per-Page Arrays](#4-dirty-tracking--bit-per-page-arrays)
5. [CPU-Side Write Detection](#5-cpu-side-write-detection)
6. [GPU-Side Write Detection — NOACCESS VEH](#6-gpu-side-write-detection--noaccess-veh)
7. [Coherency Synchronization Protocol](#7-coherency-synchronization-protocol)
8. [Thread Architecture and PFIFO Synchronization](#8-thread-architecture-and-pfifo-synchronization)
9. [Tiling, Swizzling, and Z Compression](#9-tiling-swizzling-and-z-compression)
10. [Host Resource Aliases — Texture Pool and RT Cache](#10-host-resource-aliases--texture-pool-and-rt-cache)
11. [NV2A Rendering Operations](#11-nv2a-rendering-operations)
12. [Vertex Fetch — IA Bypass Architecture](#12-vertex-fetch--ia-bypass-architecture)
13. [Vertex and Pixel Shader Emulation](#13-vertex-and-pixel-shader-emulation)
14. [HLSL SM5.0 Format Conversion Shaders](#14-hlsl-sm50-format-conversion-shaders)
15. [Tile Register Interception](#15-tile-register-interception)
16. [Complete Coherency State Machine](#16-complete-coherency-state-machine)

---

## 1. Physical Memory Architecture

The original Xbox contains 64 MiB of unified DDR SDRAM managed by the MCPX northbridge. All participants — the Pentium III CPU, the NV2A GPU, the audio DSP, and the DMA engines — share this single pool. There is no dedicated VRAM.

The NV2A is an NV20-family GPU (Kelvin class, comparable to a GeForce3/4 Ti). Key components relevant to emulation:

- A programmable vertex shader unit (VSH) with a 136-instruction RISC ISA, 16 input attribute registers (v[0]–v[15]), **192** constant vector registers (c[0]–c[191]), 12 temporaries (r[0]–r[11]), and 13 output registers (oPos, oD0/oD1, oT0–oT3, oB0/oB1, oFog, oPts).
- A fixed-function pixel pipeline built around up to 8 general register combiner stages plus a final combiner (NV_register_combiners / NV_register_combiners2 extension semantics).
- A hierarchical Z-compression unit (zcomp) with tag-based lossless Z compaction at 4×4 block granularity.
- A tiled surface aperture (PFB_TILE registers) that causes the NV2A rasterizer to pre-swizzle render target and depth buffer addresses in Morton/Z-order before they leave the chip, improving GPU cache locality.

**NV2A register state is MMIO, not RAM.** The NV2A exposes its internal engine registers through a 16 MiB MMIO window at 0xFD000000–0xFDFFFFFF on the Xbox, distinct from the 64 MiB unified RAM. The relevant register blocks within that window, their offsets and sizes, are listed in §2.3. PGRAPH alone (at 0xFD400000) is the primary source of 3D rendering state, but tile configuration (PFB), overlay (PVIDEO), timer (PTIMER), and GPU instance memory (PRAMIN) are also relevant to correct emulation. The host should allocate backing storage for this MMIO space, both so MMIO handlers can operate by direct memory write rather than dispatched struct-field assignment, and so contiguous slices of the backing can be uploaded to the host GPU as buffer SRVs and consumed directly by rendering shaders.

---

## 2. Host Memory Representation

### 2.1 The Canonical Host Mirror — s_pXboxRAM

A single host virtual allocation represents the 64 MiB of Xbox physical RAM. `MEM_WRITE_WATCH` is requested so CPU writes can be harvested in bulk via `GetWriteWatch` without paying a VEH fault on every store:

```cpp
static constexpr uint32_t XBOX_RAM_SIZE = 64 * 1024 * 1024;
static constexpr uint32_t PAGE_SIZE     = 4096;
static constexpr uint32_t PAGE_COUNT    = XBOX_RAM_SIZE / PAGE_SIZE; // 16384
static constexpr uint32_t BITMAP_BYTES  = PAGE_COUNT / 8;            // 2048

// Primary host mirror with MEM_WRITE_WATCH backing the Contiguous window.
// In Cxbx-Reloaded this is allocated at the Contiguous base (0x80000000).
uint8_t* s_pXboxRAM = (uint8_t*)VirtualAlloc(
    (void*)0x80000000,
    XBOX_RAM_SIZE,
    MEM_RESERVE | MEM_COMMIT | MEM_WRITE_WATCH,
    PAGE_READWRITE);
```

Xbox x86 code runs **directly on the host x86/x64 processor** — Cxbx-Reloaded does not JIT-recompile guest code. The Xbox Pentium III instruction set is a strict subset of x86-64, so native execution is possible once the Xbox address space is mapped into the host process. Guest stores land in this buffer at the appropriate physical offset without any translation layer. It is the ground truth for all CPU-visible memory.

### 2.2 The GPU-Side Mirror — s_pMirrorBuf (Xbox RAM Only)

**Current implementation:** A single `D3D11_USAGE_DYNAMIC` buffer (`s_pMirrorBuf`) on the host GPU holds the 64 MiB Xbox RAM image. It is bound as a raw `ByteAddressBuffer` SRV at t0 for vertex fetch (IA-bypass architecture, §12). NV2A MMIO register state is uploaded separately (§2.4).

```cpp
// Current implementation (Backend_D3D11_PageTracker.cpp)
static constexpr uint32_t CONTIG_SIZE = 64 * 1024 * 1024; // 64 MiB

D3D11_BUFFER_DESC desc = {};
desc.ByteWidth  = CONTIG_SIZE;
desc.Usage      = D3D11_USAGE_DYNAMIC;
desc.BindFlags  = D3D11_BIND_SHADER_RESOURCE;
desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
desc.MiscFlags  = D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS;
g_pD3DDevice->CreateBuffer(&desc, nullptr, &s_pMirrorBuf);

// Raw ByteAddressBuffer SRV — bound at t0
D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
srvDesc.Format               = DXGI_FORMAT_R32_TYPELESS;
srvDesc.ViewDimension        = D3D11_SRV_DIMENSION_BUFFEREX;
srvDesc.BufferEx.Flags       = D3D11_BUFFEREX_SRV_FLAG_RAW;
srvDesc.BufferEx.NumElements = CONTIG_SIZE / 4;
g_pD3DDevice->CreateShaderResourceView(s_pMirrorBuf, &srvDesc, &s_pMirrorSRV);
```

Additionally, typed SRV views are created over the same buffer for hardware-accelerated vertex attribute format decode:
- `R16G16_SNORM` at t2 — for S1 (SNORM16) vertex attributes
- `R8G8B8A8_UNORM` at t3 — for UB_OGL / UB_D3D vertex attributes

#### Proposed Design: Combined RAM + MMIO Buffer

The long-term goal is to replace `s_pMirrorBuf` with a combined `D3D11_USAGE_DEFAULT` buffer (`s_pGpuMem`) that holds both Xbox RAM and appended NV2A MMIO register blocks. This would give shaders a single SRV binding for all read access (XboxRAM + NV2A state) via offset macros, eliminate the separate PGRAPH StructuredBuffer (currently at t12), and enable UAV writes from compute shaders:

```
Proposed GPU buffer layout:
  [0x00000000 .. 0x03FFFFFF]  Xbox physical RAM        64 MiB
  [0x04000000 .. 0x04000FFF]  PFB block   (appended)   4 KB
  [0x04001000 .. 0x04002FFF]  PGRAPH block (appended)  8 KB
  [0x04003000 .. 0x04003FFF]  PVIDEO block (appended)  4 KB
  ─────────────────────────────────────────────────────────────
  Total:                                              ~64 MiB + 16 KB
```

The combined design would use `D3D11_USAGE_DEFAULT` with `D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS` and updates via `UpdateSubresource` with `D3D11_BOX` for per-page granularity.

#### Current: Per-Page Updates via Map/Unmap

The current implementation uses `Map`/`Unmap` with `WRITE_DISCARD` (first flush per frame, orphans the buffer) or `WRITE_NO_OVERWRITE` (subsequent flushes). Dirty pages from `GetWriteWatch` are coalesced into contiguous runs and `memcpy`'d into the mapped region:

```cpp
// Current implementation (Backend_D3D11_PageTracker.cpp — CxbxPageTrackerFlushToGPU)
if (count > PAGE_COUNT / 4 && s_bFirstFlushOfFrame) {
    // Bulk: DISCARD + full 64 MiB memcpy
    D3D11_MAPPED_SUBRESOURCE mapped = {};
    g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
    memcpy(mapped.pData, (void*)CONTIG_BASE, CONTIG_SIZE);
    g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
} else {
    // Incremental: NO_OVERWRITE, coalesced page-run memcpy
    D3D11_MAPPED_SUBRESOURCE mapped = {};
    g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_NO_OVERWRITE, 0, &mapped);
    // Walk GetWriteWatch results, merge consecutive pages, memcpy runs
    g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
}
```

For MMIO block updates (triggered after register writes, not on a page-fault schedule):

```cpp
// Proposed (not yet implemented):
void FlushPGRAPHToGPU() {
    D3D11_BOX box = { GPU_PGRAPH_BASE, 0, 0, GPU_PGRAPH_BASE + 0x2000u, 1, 1 };
    gD3DCtx->UpdateSubresource(s_pGpuMem, 0, &box,
                               s_pNV2AMMIO + 0x400000u, 0x2000u, 0);
}
void FlushPFBToGPU() {
    D3D11_BOX box = { GPU_PFB_BASE, 0, 0, GPU_PFB_BASE + 0x1000u, 1, 1 };
    gD3DCtx->UpdateSubresource(s_pGpuMem, 0, &box,
                               s_pNV2AMMIO + 0x100000u, 0x1000u, 0);
}
```

These are cheap (8 KB and 4 KB respectively) and can be called unconditionally before each draw, or gated on a dirty flag set by the MMIO handler when any PGRAPH or PFB register changes.

#### HLSL Accessor Macros

**Current implementation:** Shaders use two separate bindings:
- `ByteAddressBuffer g_XboxRAM : register(t0)` — the 64 MiB mirror buffer for vertex/index data
- `StructuredBuffer<uint> g_PGRegs : register(t12)` — PGRAPH register block (2048 × uint32), accessed via `PG_UINT(offset)` / `PG_FLOAT(offset)` helpers defined in `CxbxPGRAPHRegs.hlsli`

PFB registers are not accessible from shaders in the current implementation.

**Proposed design:** A single `ByteAddressBuffer` binding for all memory reads, with macros that translate documentation offsets to GPU-buffer offsets:

```hlsl
// Proposed: t0 — the combined GPU memory buffer
ByteAddressBuffer g_GpuMem : register(t0);

// Xbox physical RAM — identity-mapped within the buffer
#define GpuRamLoad(physAddr)        g_GpuMem.Load(physAddr)
#define GpuRamLoad4(physAddr)       g_GpuMem.Load4(physAddr)

// NV2A PGRAPH registers — offset by GPU_PGRAPH_BASE
#define PGRAPH_BASE  0x04001000u
#define PFB_BASE     0x04000000u
#define PVIDEO_BASE  0x04003000u

#define PGRAPHLoad(reg)    g_GpuMem.Load(PGRAPH_BASE + (reg))
#define PFBLoad(reg)       g_GpuMem.Load(PFB_BASE    + (reg))
#define PVIDEOLoad(reg)    g_GpuMem.Load(PVIDEO_BASE + (reg))

// Example usage:
uint numStages = (PGRAPHLoad(NV_PGRAPH_COMBINECTL) >> 0) & 0xFu;
uint tile0Base  = PFBLoad(NV_PFB_TILE(0));
float4 vConst3 = asfloat(GpuRamLoad4(VS_CONSTANTS_BASE + 3 * 16));
```

### 2.3 NV2A MMIO Backing Storage

The NV2A MMIO window spans 16 MiB at 0xFD000000 on Xbox.

**Current implementation:** NV2A state is stored in a C struct hierarchy. `NV2AState` (defined in `nv2a_int.h`) contains per-engine sub-structs, each with a `regs[]` array indexed by block-relative register offset divided by 4. Key structures:

```cpp
// Current implementation (nv2a_int.h)
struct PGRAPHState {
    uint32_t regs[2048];             // PGRAPH registers indexed by (offset/4)
    float    vsh_constants[192][4];  // Transform context RAM (XFCTX)
    uint32_t program_data[136][4];   // Transform program RAM (XFPR)
    // ... vertex attributes, surface state, etc.
};

struct NV2AState {
    struct { uint32_t regs[...]; } pmc;
    struct { uint32_t regs[...]; } pfifo;
    struct { uint32_t regs[...]; } pfb;
    PGRAPHState pgraph;
    struct { uint32_t regs[...]; } pcrtc;
    // ... other engine blocks
};
```

MMIO accesses from guest code trigger an access violation (the 0xFD000000 range has no backing memory). The unified VEH handler catches these faults and dispatches to per-engine read/write handlers that operate on the `NV2AState` struct fields.

**Proposed design:** Allocate a contiguous 16 MiB flat backing region so that MMIO handlers can operate by direct memory write at `backing + (guestAddr - 0xFD000000)` rather than by dispatching to struct fields. This simplifies register upload (contiguous `memcpy` of a block to the GPU buffer) and makes the memory layout match documentation offsets exactly:

```cpp
// Proposed (not yet implemented):
static constexpr uint32_t NV2A_MMIO_BASE = 0xFD000000;
static constexpr uint32_t NV2A_MMIO_SIZE = 0x01000000; // 16 MiB

uint8_t* s_pNV2AMMIO = (uint8_t*)VirtualAlloc(
    nullptr,
    NV2A_MMIO_SIZE,
    MEM_RESERVE | MEM_COMMIT,
    PAGE_READWRITE);
```

The relevant engine blocks and their offsets within this backing region:

```
Block      Offset      Size   Key contents
─────────  ──────────  ─────  ──────────────────────────────────────────────
PMC        0x000000    4 KB   Interrupt routing, engine enables
PBUS       0x001000    4 KB   Bus control, PCI config mirror
PFIFO      0x002000    8 KB   DMA GET/PUT, CACHE1, channel control
PVIDEO     0x008000    4 KB   Video overlay address, format, position
PTIMER     0x009000    4 KB   GPU clock counter (TIME_0 / TIME_1)
PFB        0x100000    4 KB   Tile regions (TILE[0-7], TLIMIT, TSIZE, ZCOMP)
PGRAPH     0x400000    8 KB   All 3D state: combiners, VS, texture, surface
PCRTC      0x600000    4 KB   Display scan-out base address, interrupt
PRAMDAC    0x680000    4 KB   DAC, PLL, cursor
PRAMIN     0x700000    1 MB   GPU instance memory (RAMHT, RAMFC, DMA objects)
USER       0x800000    4 MB   PFIFO DMA submission (MMIO kick registers)
UREMAP     0xC00000    4 MB   Mirror of USER
```

MMIO read/write handlers receive a block-relative offset and a pointer directly into the backing:

```cpp
uint32_t* NV2ARegPtr(uint32_t guestAddr) {
    uint32_t offset = guestAddr - NV2A_MMIO_BASE;
    assert(offset + 4 <= NV2A_MMIO_SIZE);
    return reinterpret_cast<uint32_t*>(s_pNV2AMMIO + offset);
}

// Side-effect-free register write (most PGRAPH state registers):
void NV2ARegWrite(uint32_t guestAddr, uint32_t value) {
    *NV2ARegPtr(guestAddr) = value;
}

// Side-effect-bearing writes (PFIFO kick, PMC interrupt ack, PCRTC scan-out
// base, PTIMER clock write) still go through the full dispatch table but
// also commit to the backing store so reads see consistent state:
void NV2ARegWriteWithSideEffect(uint32_t guestAddr, uint32_t value,
                                 NV2AState* d) {
    *NV2ARegPtr(guestAddr) = value;
    DispatchMMIOSideEffect(guestAddr, value, d);
}
```

**PRAMIN** (0x700000, 1 MiB) contains the GPU's instance memory: the RAMHT (register hash table for DMA object handles), RAMFC (FIFO context save areas), and the DMA object descriptors the PFIFO puller resolves to find surface base addresses. PRAMIN is backed by a region of Xbox physical RAM (the top of the 64 MiB pool), so reads/writes to the PRAMIN MMIO block alias the same physical storage as the corresponding offset in `s_pXboxRAM`. The emulator can implement this by mapping `s_pNV2AMMIO + 0x700000` as a view over `s_pXboxRAM + PRAMIN_PHYS_BASE` rather than a separate allocation, keeping both in sync without copies.

**USER / UREMAP** (0x800000–0xFFFFFF, 8 MiB combined) is the DMA submission aperture. Guest writes to USER channel slots (e.g. `NV_USER_DMA_PUT`) kick the PFIFO pusher. These are not general state registers and do not need to be uploaded to the GPU. The backing store records the last written values for read-back consistency; the actual side effect (pusher wake) is dispatched from the MMIO handler.

### 2.4 Uploading MMIO Backing to the Host GPU

**Current implementation:** PGRAPH registers are uploaded to the host GPU as a separate `StructuredBuffer<uint>` (`g_pD3D11PGRegsBuf`, 2048 × uint32 = 8 KB) bound at slot t12. The upload uses `D3D11_USAGE_DYNAMIC` + `Map/Unmap`, gated by a dirty generation counter that avoids redundant uploads when PGRAPH state is unchanged between draws:

```cpp
// Current implementation (XbPixelShaderCompiler.cpp — CxbxD3D11UploadRCInterpreterState)
if (s_lastPGRegsGeneration != pg->regs_generation) {
    D3D11_MAPPED_SUBRESOURCE mapped = {};
    g_pD3DDeviceContext->Map(g_pD3D11PGRegsBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
    memcpy(mapped.pData, pg->regs, 2048 * sizeof(uint32_t));
    g_pD3DDeviceContext->Unmap(g_pD3D11PGRegsBuf, 0);
    s_lastPGRegsGeneration = pg->regs_generation;
}
// Bound at t12 — shaders access via g_PGRegs[offset >> 2]
```

Shaders access PGRAPH state through `StructuredBuffer<uint> g_PGRegs : register(t12)` with helper macros `PG_UINT(reg)` and `PG_FLOAT(reg)` defined in `CxbxPGRAPHRegs.hlsli`.

**PFB registers are NOT uploaded to the GPU.** Tile region configuration is consumed CPU-side only (for surface invalidation and upload path decisions). Shaders cannot currently read PFB_TILE state.

**Proposed design (combined buffer):** With the combined buffer design (§2.2), MMIO register blocks would become part of `s_pGpuMem` — appended after the 64 MiB RAM region. Per-draw state would be refreshed via `UpdateSubresource` with a `D3D11_BOX` for the PGRAPH (8 KB) and PFB (4 KB) regions. Shaders would access both Xbox RAM and NV2A register state through the single `g_GpuMem` ByteAddressBuffer at t0, eliminating the separate t12 binding.

If the combined design is not used, the current approach (separate `StructuredBuffer<uint>` SRVs) remains:

```cpp
// Current implementation (separate SRVs for PGRAPH block)
// g_pD3D11PGRegsBuf at t12, uploaded via Map/Unmap per draw (gated by generation counter)
// Shaders use: g_PGRegs[NV_PGRAPH_COMBINECTL >> 2]
```

### 2.5 PGRAPH Sub-State — XFPR, XFCTX, and Derived Fields

Three categories of PGRAPH state require special treatment because they are not directly addressable as single register reads:

**Transform Program RAM (XFPR)** — stores the VSH instruction stream. `NV097_SET_TRANSFORM_PROGRAM` writes 128-bit instructions to the XFPR starting at the slot selected by `NV_PGRAPH_CHEOPS_OFFSET`. In the current implementation, this is stored as `pg->program_data[136][4]` in the `PGRAPHState` struct. For the VS interpreter path, it is uploaded as `StructuredBuffer<uint4>` at t5 (`g_XFPR`). For the JIT path, it is compiled to HLSL at shader-upload time and cached by instruction-stream hash (rapidhash).

**Transform Context RAM (XFCTX)** — stores the 192 VS constant vectors. Written via `NV097_SET_TRANSFORM_CONSTANT_LOAD` + `NV097_SET_TRANSFORM_CONSTANT`. Maintained as `pg->vsh_constants[192][4]` and uploaded as `cbuffer b0` (192 × float4 = 3072 bytes) every draw.

**Software-computed fields** — a small set of values that shaders need but that cannot be derived from a single register read: ColorSign (sign-extension flags for per-vertex color inputs), TexFmtFixup (per-stage texture format quirks), AlphaKill threshold, FrontFaceInfo. These are derived by the CPU at draw time from the relevant PGRAPH registers and placed in a small auxiliary cbuffer (`PSAuxCBLayout`). The goal of eliminating this cbuffer by deriving all fields in-shader from the PGRAPH SRV is achievable but deferred.

### 2.6 The Staging Buffer

A CPU-readable staging texture (not a buffer) is used per RT for readback operations — `CopyResource` requires matching resource types. A staging buffer of the same size as `s_pXboxRAM` also exists for buffer-to-CPU readback paths that do not go through a Texture2D:

```cpp
D3D11_BUFFER_DESC sd = {};
sd.ByteWidth      = XBOX_RAM_SIZE;
sd.Usage          = D3D11_USAGE_STAGING;
sd.BindFlags      = 0;
sd.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
sd.MiscFlags      = 0;
gDevice->CreateBuffer(&sd, nullptr, &s_pStagingBuf);
```

---

## 3. Address Mirrors and Physical Resolution

### 3.1 The Three Windows

```
Guest Virtual Address    Window         Semantics
──────────────────────   ──────────     ──────────────────────────────────────────
0x00000000–0x03FFFFFF    Identity       1:1 physical. CPU + DMA. Direct.
0x80000000–0x83FFFFFF    Contiguous     Same physical, virtually linear.
                                        CPU heap, stack, XDK allocs.
                                        Primary window for the Cxbx write-watch mirror.
0xF0000000–0xF3FFFFFF    Tiled          Same physical. NV2A pre-swizzles addresses
                                        in Morton order before they leave the chip.
                                        CPU accessing this window gets raw linear bytes —
                                        no tiling hardware on the CPU side.
```

All three windows address the same physical DRAM via the MCPX. The tiling transform is applied **inside the NV2A** before the address leaves the chip. The memory controller sees a pre-swizzled linear bus address and performs a straightforward DRAM cell lookup. DRAM has no knowledge of tiling.

### 3.2 What the Tiled Window Actually Means

When the NV2A rasterizer writes a pixel to a render target whose physical address falls within an active PFB_TILE region, the chip's internal addressing unit applies Morton interleave to the pixel's (x, y) screen position, then issues that transformed address on the memory bus. The CPU, lacking this unit, accesses the same physical cells through the Contiguous or Identity window in strictly linear (pitch × row + column) order.

Practical consequence: a CPU read of a tiled RT surface via the Contiguous window sees Morton-encoded data unless the emulator explicitly untiles on readback. The function `CxbxSyncTiledRangeToContiguous()` handles this for explicit cases such as PVIDEO overlay composition, where the CPU-decoded YUV frame is written through the tiled window and must be readable linearly before the display compositor accesses it.

### 3.3 Physical Page Resolution

```cpp
uint32_t GuestVAToPhysical(uint32_t va) {
    if (va < 0x04000000)                     return va;               // Identity
    if (va >= 0x80000000 && va < 0x84000000) return va - 0x80000000;  // Contiguous
    if (va >= 0xF0000000 && va < 0xF4000000) return va - 0xF0000000;  // Tiled (linear from DRAM POV)
    return UINT32_MAX;                                                 // MMIO / unmapped
}
```

From DRAM's perspective all three windows produce linear physical addresses. The difference is only in logical interpretation: NV2A-generated tiled addresses encode spatial (x, y) locality in their physical bits; CPU-generated addresses do not.

---

## 4. Dirty Tracking — Bit-per-Page Arrays

### 4.1 Four Bitmaps (Current Implementation)

Four independent bit-per-page arrays track distinct coherency aspects. All are stored as `uint32_t[512]` (DWORD-granularity, 1 bit per 4KB page × 16384 pages):

```cpp
// Current implementation (Backend_D3D11_PageTracker.cpp)
static uint32_t s_GpuDirtyBitmap[BITMAP_DWORDS];      // GPU wrote → readback on CPU access
static uint32_t s_TextureDirtyBitmap[BITMAP_DWORDS];   // CPU wrote → texture re-upload needed
static uint32_t s_TiledCommittedBitmap[BITMAP_DWORDS]; // 0xF0 pages committed on demand
static uint32_t s_IdentityAllocBitmap[BITMAP_DWORDS];  // physical pages with active contiguous allocs
```

`s_TextureDirtyBitmap` records pages whose content has changed since the last texture upload for a given surface. It is set when CPU-dirty pages are flushed to the GPU mirror, and checked per texture per draw to decide whether a deswizzle/re-upload is needed. It is cleared per-texture after the host texture is updated.

`s_GpuDirtyBitmap` records pages written by the GPU (render targets). Pages marked here are protected `PAGE_NOACCESS`. It gates the VEH readback path.

`s_TiledCommittedBitmap` tracks which 4KB pages in the tiled window (0xF0000000) have been committed on demand by the VEH. At flush time, committed pages are synced back to contiguous memory and decommitted.

`s_IdentityAllocBitmap` tracks which physical pages correspond to active contiguous memory allocations. Only identity-range dirty pages overlapping these allocations are synced to the contiguous range.

**Note:** The document's proposed design (§16) describes a third coherency bitmap (`s_NoaccessArmed`) for tracking which pages are currently `PAGE_NOACCESS`. The current implementation does not maintain this — `PAGE_NOACCESS` state is implicitly derivable from `s_GpuDirtyBitmap` (pages are set NOACCESS exactly when marked GPU-dirty, and restored when the bit is cleared).

### 4.2 Bit Operations

```cpp
// Current implementation uses non-atomic DWORD-granularity operations.
// Single-threaded access to each bitmap makes atomics unnecessary in practice.
static inline void SetBit(uint32_t* bitmap, uint32_t index) {
    bitmap[index >> 5] |= (1u << (index & 31));
}
static inline void ClearBit(uint32_t* bitmap, uint32_t index) {
    bitmap[index >> 5] &= ~(1u << (index & 31));
}
static inline bool TestBit(const uint32_t* bitmap, uint32_t index) {
    return (bitmap[index >> 5] & (1u << (index & 31))) != 0;
}
```

**Proposed improvement:** Cache-line-aligned `uint8_t` bitmaps with `_InterlockedOr8`/`_InterlockedAnd8` for lock-free atomic access from multiple threads. On x86 this emits a single `lock or byte ptr [addr], imm8` — cheaper than a CAS loop for monotone set operations (bits only go 0→1 except during deliberate clear).

### 4.3 Range Scan

The texture-dirty check scans the bitmap for any set bits in a page range. The current implementation uses DWORD-granularity scanning with bitmasks for partial first/last DWORDs:

```cpp
// Current implementation (CxbxPageTrackerIsTextureDirty)
bool CxbxPageTrackerIsTextureDirty(uint32_t offset, uint32_t size) {
    uint32_t firstPage = offset / PAGE_SIZE_;
    uint32_t lastPage  = (offset + size - 1) / PAGE_SIZE_;
    uint32_t firstDW   = firstPage >> 5;
    uint32_t lastDW    = lastPage >> 5;

    if (firstDW == lastDW) {
        uint32_t mask = ((2u << (lastPage & 31)) - 1) & ~((1u << (firstPage & 31)) - 1);
        return (s_TextureDirtyBitmap[firstDW] & mask) != 0;
    }
    // Check partial first/last DWORDs + full DWORDs in between
    ...
}
```

**Proposed improvement:** With AVX2, 32 bytes (256 pages = 1 MiB of Xbox RAM) can be tested per instruction for bulk scanning:

```cpp
bool AnyDirtyInRange(const uint8_t* bitmap, uint32_t firstPage, uint32_t lastPage) {
    uint32_t b0 = firstPage >> 3, b1 = lastPage >> 3;
#ifdef __AVX2__
    for (uint32_t b = b0 & ~31u; b <= b1; b += 32) {
        __m256i v = _mm256_load_si256((const __m256i*)(bitmap + b));
        if (!_mm256_testz_si256(v, v)) return true;
    }
    return false;
#else
    // ... scalar fallback
#endif
}
```

Total coherency metadata for all 64 MiB of Xbox RAM: 4 × 2048 = 8192 bytes (current), fully resident in L1 cache.

---

## 5. CPU-Side Write Detection

### 5.1 GetWriteWatch Bulk Harvest — Once Per Frame

After executing a quantum of guest CPU instructions, once per frame at the first draw after `CxbxPresent()`, the emulator calls `GetWriteWatch` to collect all pages written since the last harvest. The OS maintains this at zero per-access cost in the VMM:

```cpp
// Current implementation (CxbxPageTrackerFlushToGPU)
ULONG_PTR count = PAGE_COUNT;
ULONG granularity;
UINT result = GetWriteWatch(
    WRITE_WATCH_FLAG_RESET,
    (PVOID)CONTIG_BASE, CONTIG_SIZE,
    s_WriteWatchPages, &count, &granularity);

// Mark all flushed pages as texture-dirty (gates texture re-upload)
for (ULONG_PTR i = 0; i < count; i++) {
    uint32_t pageIdx = (uint32_t)((uintptr_t)s_WriteWatchPages[i] - CONTIG_BASE) / PAGE_SIZE_;
    SetBit(s_TextureDirtyBitmap, pageIdx);
}
```

**Once-per-frame is correct.** Because Xbox code runs natively, all guest memory writes are real host stores that the VMM's write-watch tracking captures immediately. The Xbox pushbuffer model additionally guarantees all CPU memory writes complete before the pushbuffer ring is kicked — mid-frame CPU writes to vertex or texture data are structurally impossible in the single-threaded Xbox guest model. The flush is gated by `s_bFirstFlushOfFrame` (reset by `CxbxPageTrackerOnPresent()`), ensuring only one `GetWriteWatch` call per frame rather than one per draw.

**Wine fallback:** When running under Wine (where `GetWriteWatch` may be unreliable), the implementation falls back to a full 64 MiB `memcpy` + marking all pages texture-dirty on every flush.

**Identity range sync:** Before the contiguous-range harvest, dirty pages from the identity-mapped range (0x00010000+) are also synced. Games using `Lock2DSurface` write texture data to physical addresses via identity mapping. Those pages are detected via a second `GetWriteWatch` call on the identity allocation and `memcpy`'d to the corresponding contiguous address, where the normal flush path picks them up.

### 5.2 Flushing Dirty Pages to the GPU Mirror

Dirty pages are flushed from `s_pXboxRAM` to the mirror buffer via `Map`/`Unmap`. Adjacent dirty pages are coalesced into contiguous runs for efficient `memcpy`:

```cpp
// Current implementation (CxbxPageTrackerFlushToGPU)
void CxbxPageTrackerFlushToGPU() {
    if (!s_bFirstFlushOfFrame) return; // Skip redundant mid-frame flushes

    SyncTiledPagesBack();                         // Tiled→contiguous
    CxbxPageTrackerSyncIdentityToContiguous();    // Identity→contiguous

    // Bulk dirty (>25% of pages): DISCARD + full memcpy
    if (count > PAGE_COUNT / 4 && s_bFirstFlushOfFrame) {
        g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
        memcpy(mapped.pData, (void*)CONTIG_BASE, CONTIG_SIZE);
        g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
    } else {
        // Incremental: NO_OVERWRITE, coalesced contiguous runs
        g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_NO_OVERWRITE, 0, &mapped);
        // Walk page list, merge consecutive pages, memcpy per run
        g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
    }
    s_bFirstFlushOfFrame = false;
}
```

Because `s_pMirrorBuf` is `D3D11_USAGE_DYNAMIC`, `WRITE_DISCARD` orphans the old allocation (safe for prior draws that are GPU-pipelined against the old buffer instance) and `WRITE_NO_OVERWRITE` promises the driver that only regions not currently in flight are being written.

---

## 6. GPU-Side Write Detection — NOACCESS VEH

### 6.1 PAGE_NOACCESS, Not PAGE_GUARD

When the GPU renders to an RT, those physical pages are marked GPU-dirty and protected with `PAGE_NOACCESS`. Unlike `PAGE_GUARD` (which fires only on the first access and then self-clears to the underlying protection), `PAGE_NOACCESS` fires on every access — both reads and writes — until explicitly removed:

```cpp
// Current implementation (CxbxPageTrackerMarkGPUDirty)
void CxbxPageTrackerMarkGPUDirty(uint32_t startOffset, uint32_t size) {
    uint32_t firstPage = startOffset / PAGE_SIZE_;
    uint32_t lastPage  = (startOffset + size - 1) / PAGE_SIZE_;

    for (uint32_t p = firstPage; p <= lastPage; p++) {
        SetBit(s_GpuDirtyBitmap, p);
        DWORD oldProtect;
        VirtualProtect((void*)(CONTIG_BASE + p * PAGE_SIZE_),
            PAGE_SIZE_, PAGE_NOACCESS, &oldProtect);
    }

    // Mark overlapping registered RTs as needing readback
    for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
        RegisteredRT& rt = s_RegisteredRTs[i];
        uint32_t rtEnd = rt.offset + rt.pitch * rt.height;
        if (startOffset < rtEnd && rt.offset < startOffset + size)
            rt.needsReadback = true;
    }
}
```

**Note:** Unlike the proposed design, there is no separate `s_NoaccessArmed` bitmap. The GPU-dirty bitmap doubles as the armed-state tracker — pages are set NOACCESS exactly when `SetBit(s_GpuDirtyBitmap, p)` is called, and restored when the bit is cleared.

### 6.2 VEH Handler

The fault handler is integrated into the unified VEH dispatcher (`lleException`), not registered as a standalone VEH. It handles both GPU-dirty pages (contiguous range) and tiled memory commit-on-demand:

```cpp
// Current implementation (CxbxPageTrackerHandleFault)
bool CxbxPageTrackerHandleFault(void* faultAddress, bool isWrite) {
    uintptr_t addr = (uintptr_t)faultAddress;

    // --- Tiled memory redirect (0xF0000000 - 0xF3FFFFFF) ---
    if (addr >= TILED_BASE && addr < (TILED_BASE + TILED_SIZE)) {
        uint32_t pageIdx = (uint32_t)(addr - TILED_BASE) / PAGE_SIZE_;
        if (!TestBit(s_TiledCommittedBitmap, pageIdx)) {
            // Commit page on demand, copy current data from contiguous
            VirtualAlloc((LPVOID)(TILED_BASE + pageIdx * PAGE_SIZE_),
                         PAGE_SIZE_, MEM_COMMIT, PAGE_READWRITE);
            memcpy((void*)(TILED_BASE + pageIdx * PAGE_SIZE_),
                   (void*)(CONTIG_BASE + pageIdx * PAGE_SIZE_), PAGE_SIZE_);
            SetBit(s_TiledCommittedBitmap, pageIdx);
        }
        return true;
    }

    // --- GPU-dirty page handling (0x80000000 - 0x83FFFFFF) ---
    if (addr >= CONTIG_BASE && addr < (CONTIG_BASE + CONTIG_SIZE)) {
        uint32_t pageIdx = (uint32_t)(addr - CONTIG_BASE) / PAGE_SIZE_;
        if (!TestBit(s_GpuDirtyBitmap, pageIdx))
            return false;

        if (!isWrite && g_pD3DDeviceContext != nullptr &&
            TryEnterCriticalSection(&s_D3D11ContextLock)) {
            // Find registered RT covering this page
            const RegisteredRT* pRT = FindRegisteredRTForPage(pageIdx);
            if (pRT && pRT->pTexture) {
                // Clear dirty + restore READWRITE for ENTIRE RT first
                uint32_t rtSize = pRT->pitch * pRT->height;
                uint32_t rtFirstPage = pRT->offset / PAGE_SIZE_;
                uint32_t rtLastPage = (pRT->offset + rtSize - 1) / PAGE_SIZE_;
                for (uint32_t p = rtFirstPage; p <= rtLastPage; p++) {
                    if (TestBit(s_GpuDirtyBitmap, p)) {
                        ClearBit(s_GpuDirtyBitmap, p);
                        DWORD old;
                        VirtualProtect((void*)(CONTIG_BASE + p * PAGE_SIZE_),
                                       PAGE_SIZE_, PAGE_READWRITE, &old);
                    }
                }
                // Staged readback: host RT → staging texture → Map → memcpy
                // Creates staging texture per readback (not cached)
                D3D11_TEXTURE2D_DESC stagingDesc = ...;
                g_pD3DDeviceContext->CopyResource(pStaging, pRT->pTexture);
                g_pD3DDeviceContext->Map(pStaging, 0, D3D11_MAP_READ, 0, &m);
                CopyLinearWithPitch(...); // or swizzle_rect for SWIZZLE surfaces
                g_pD3DDeviceContext->Unmap(pStaging, 0);
                pStaging->Release();
            }
            LeaveCriticalSection(&s_D3D11ContextLock);
            return true;
        }

        // Fallback: write access, lock contended, or no matching RT.
        // Restore page, return stale data (graceful degradation).
        ClearBit(s_GpuDirtyBitmap, pageIdx);
        DWORD old;
        VirtualProtect((void*)(CONTIG_BASE + pageIdx * PAGE_SIZE_),
                       PAGE_SIZE_, PAGE_READWRITE, &old);
        return true;
    }
    return false;
}
```

**Key differences from proposed design:**
- Staging texture is created per-readback rather than cached per-RT (allocation cost per fault).
- The critical section uses `TryEnterCriticalSection` — on contention, the handler degrades gracefully by returning stale data rather than blocking.
- CPU writes to GPU-dirty pages restore access and discard GPU data without readback (CPU wins).

---

## 7. Coherency Synchronization Protocol

### 7.1 Upload: CPU Dirty → Mirror Buffer

For vertex and index data (Draw Arrays path), dirty pages are flushed from contiguous memory to `s_pMirrorBuf` in `CxbxPageTrackerFlushToGPU()` at the first draw after Present. The VS then samples data directly from the mirror buffer via its ByteAddressBuffer SRV (t0).

For cases where a vertex buffer aliases GPU-rendered memory, `CxbxPageTrackerFlushGPUDirtyToMirror()` performs an additional targeted readback: it reads back the D3D11 RT into Xbox RAM, then copies the affected pages to the mirror buffer via `MAP_WRITE_NO_OVERWRITE`.

For texture data, an additional deswizzle/decode pass is required after the flush — see §10.

### 7.2 RT-as-Texture Fast Path (GPU→GPU, Zero CPU Round-Trip)

When a subsequent NV2A draw samples a surface that was previously used as a render target (TEXOFFSET matches a known RT's physBase), the same `ID3D11Texture2D` that served as the RTV is now bound as an SRV. No readback to `s_pXboxRAM` and no re-upload are needed:

```cpp
// Called from CxbxUpdateHostTextures() at draw time
ID3D11Texture2D* CxbxLookupPgraphRTByOffset(uint32_t physBase) {
    auto it = g_PgraphRTCache.find({physBase, currentFmt, currentW, currentH});
    return (it != g_PgraphRTCache.end()) ? it->second.Get() : nullptr;
}
// If found:
//   OMSetRenderTargets(0, nullptr, nullptr)  — unbind RTV (DX11 hazard rule)
//   CreateShaderResourceView(pTex, nullptr, &pSRV)
//   PSSetShaderResources(stage, 1, &pSRV)
```

This is the correct fast path for shadow maps, reflection captures, and post-process intermediates — the common GPU→GPU surface feedback case.

---

## 8. Thread Architecture and PFIFO Synchronization

```
Host Core 0 (pinned)
│
├─ Guest CPU Thread
│    Xbox x86 code runs NATIVELY on the host x86/x64 CPU.
│    Cxbx-Reloaded does NOT JIT-compile guest code. The Pentium III ISA
│    is a strict subset of x86-64; the Xbox executable runs directly once
│    the Xbox address space is mapped into the host process.
│
│    Guest writes reach s_pXboxRAM directly (the Contiguous window mapping
│    is at 0x80000000 in both guest and host VA space). MEM_WRITE_WATCH
│    passively records touched pages without any per-store overhead.
│
│    EMUPATCH hooks (subhook trampolines) intercept specific Xbox D3D API
│    entry points whose behaviour cannot run unmodified on the host:
│    resource creation, swap-chain management, DMA push-buffer submission,
│    and a shrinking set of state-mirrors during the HLE transition period.
│    All other Xbox code — game logic, physics, animation, pushbuffer
│    construction — runs entirely unpatched at native host speed.
│
│    Cxbx also implements the Xbox kernel (NT-derived): MM (virtual memory),
│    IO (device I/O), object manager, thread scheduler, and other system
│    functional areas. Xbox kernel calls from game code are intercepted and
│    handled by this host-side kernel implementation.
│
│    VEH serves two roles:
│      (1) MMIO dispatch — NV2A register accesses (0xFD000000 range) are
│          not backed by real memory. The resulting EXCEPTION_ACCESS_VIOLATION
│          is caught and dispatched to the appropriate NV2A engine handler,
│          which reads/writes s_pNV2AMMIO and fires any side effects (PFIFO
│          kick, interrupt acknowledgement, PCRTC scan-out base change, etc.).
│      (2) GPU-dirty readback — pages protected with PAGE_NOACCESS after
│          GPU rendering trigger a fault on CPU read (or write), handled
│          by the readback path described in §6.
│
│    HarvestCPUWrites() is called once per frame at the pushbuffer kick.
│
├─ PFIFO Pusher Thread
│    Reads USER channel MMIO (DMA PUT pointer kicks from guest code via VEH).
│    Parses the DMA push buffer ring from s_pXboxRAM.
│    Decodes command packet headers: method, subchannel, count, NI/jump type.
│    Enqueues (method, parameter) pairs to CACHE1 (32-entry ring buffer).
│    CRITICAL: DMA_GET must NOT be auto-advanced to PUT when the pusher is
│    enabled, to avoid dropping commands during active DMA.
│
├─ PFIFO Puller / PGRAPH Thread
│    Dequeues from CACHE1, dispatches to pgraph_handle_method().
│    All NV097 methods write decoded parameters to pg->regs[] (PGRAPH struct).
│    Side-effect methods additionally update CPU-side mirrors:
│      pg->program_data[][]   — NV097_SET_TRANSFORM_PROGRAM (XFPR)
│      pg->vsh_constants[][]  — NV097_SET_TRANSFORM_CONSTANT (XFCTX)
│      pg->vertex_attributes[]— NV097_SET_VERTEX_DATA_ARRAY_*
│    On NV097_SET_BEGIN_END(0): triggers pgraph_draw() → D3D11_draw callback.
│    This IS the rendering thread — draws happen directly on the puller.
│    Calls CxbxPageTrackerMarkGPUDirty() after RT writes.
│
└─ PIO / Peripheral Chips Thread
     PAUDIO DSP, PTIMER, interrupt controller ticks.
     Audio DMA writes to s_pXboxRAM → SetTextureDirty() on affected pages.
```

### 8.1 Puller-Driven Draw Architecture (Current)

The current implementation uses the target LLE puller-driven architecture. The puller thread dequeues commands from CACHE1, dispatches them via `pgraph_handle_method()`, and triggers the DX11 draw directly on `NV097_SET_BEGIN_END(0)` via a function-pointer callback (`pgraph_draw` → `D3D11_draw`).

There is no `pfifo_flush_to_pgraph()` synchronization point between an EMUPATCH and the puller — draws are driven entirely by the puller thread processing the pushbuffer ring. EMUPATCH hooks remain for a shrinking set of operations (resource creation, swap-chain management, DMA submission kick) but do not trigger draws directly.

---

## 9. Tiling, Swizzling, and Z Compression

### 9.1 Two Distinct Memory Transforms

**Swizzled textures** (CPU-uploaded static art and read-only surfaces): The NV2A stores these in Morton Z-order at power-of-two dimensions. The XDK's `XGSwizzleRect` / `D3DXLoadSurfaceFromSurface` applies the swizzle when the CPU writes texture data. Dimensions are encoded in TEXFMT as log2(USIZE) / log2(VSIZE) bit fields. Swizzled surfaces are always POT.

**Tiled render targets** (via PFB_TILE registers): Active tile regions cause the NV2A rasterizer to pre-swizzle pixel addresses during rendering. Up to 8 concurrent regions are configurable. This is independent of texture swizzle — it operates at the memory bus level for GPU-generated writes only. Linear-pitch surfaces can fall within tiled regions.

These two mechanisms are independent. A surface can be swizzled (CPU-written texture), tiled via PFB (GPU-rendered RT), linear (pitch-based texture or RT with tiling disabled), or some combination.

### 9.2 PFB Tile Register Set

```
PFB block (offset from NV2A MMIO base, 0xFD000000 on Xbox):
  0x00100240 + N×0x10   PFB_TILE[N]      bits 31:16 = base >> 16; bit 0 = enable
  0x00100244 + N×0x10   PFB_TLIMIT[N]    bits 31:16 = limit >> 16
  0x00100248 + N×0x10   PFB_TSIZE[N]     bits 31:16 = pitch in 256-byte units

PGRAPH block (must be kept in sync with PFB):
  0x00400904 + N×0x10   PGRAPH_TILE[N]
  0x00400908 + N×0x10   PGRAPH_TLIMIT[N]
  0x0040090C + N×0x10   PGRAPH_TSIZE[N]
  0x00400910 + N×0x10   PGRAPH_TSTATUS[N]

Zcomp regions:
  0x00100300 + N×0x08   PFB_ZCOMP[N]
  0x00100304 + N×0x08   PFB_ZCLIMIT[N]
  0x00400980 + N×0x10   PGRAPH_ZCOMP[N]  (mirror)
```

XDK surface setup code writes both PFB and PGRAPH register sets together. Both must be intercepted by the emulator.

### 9.3 Morton Swizzle (Texture Swizzle)

NV2A swizzled textures use a Morton Z-order interleave of bit positions:

```cpp
// Linear (x, y) → swizzled texel index within a POT surface
uint32_t MortonEncode(uint32_t x, uint32_t y) {
    uint32_t result = 0;
    for (int b = 0; b < 16; b++) {
        result |= ((x >> b) & 1u) << (2*b);
        result |= ((y >> b) & 1u) << (2*b + 1);
    }
    return result;
}

// Swizzled index → x component
uint32_t MortonX(uint32_t n) {
    n &= 0x55555555u;
    n = (n | (n >> 1)) & 0x33333333u;
    n = (n | (n >> 2)) & 0x0F0F0F0Fu;
    n = (n | (n >> 4)) & 0x00FF00FFu;
    return (n | (n >> 8)) & 0x0000FFFFu;
}
uint32_t MortonY(uint32_t n) { return MortonX(n >> 1); }
```

For DXT-compressed textures, the swizzle operates at the **block** level — each 4×4 DXT block is treated as a single unit in the Morton ordering, not individual texels.

Deswizzle is the primary pre-draw CPU cost for textures. It is dispatched as a compute shader writing to a `Texture2D UAV` in the texture pool (§10), not applied to the XboxRAM buffer.

### 9.4 PFB Tile Region Swizzle (Render Targets)

The PFB tiling transform operates within 4KB tiles with a Morton interleave of pixel x/y bits:

```cpp
// Linear surface (x, y) → tiled physical byte offset within the tile region
uint32_t LinearToTiled(uint32_t x, uint32_t y, uint32_t pitch, uint32_t bpp) {
    uint32_t tileW       = 256 / bpp;           // pixels per tile row (64 for 32bpp)
    uint32_t tileH       = 64;
    uint32_t tileBytes   = tileW * bpp * tileH;
    uint32_t tilesPerRow = pitch / (tileW * bpp);
    uint32_t mortonIdx   = MortonEncode(x % tileW, y % tileH);
    return (y / tileH * tilesPerRow + x / tileW) * tileBytes + mortonIdx * bpp;
}

uint32_t TiledToLinear(uint32_t tiledOff, uint32_t pitch, uint32_t bpp) {
    uint32_t tileW       = 256 / bpp;
    uint32_t tileH       = 64;
    uint32_t tileBytes   = tileW * bpp * tileH;
    uint32_t tilesPerRow = pitch / (tileW * bpp);
    uint32_t tileIndex   = tiledOff / tileBytes;
    uint32_t mortonIdx   = (tiledOff % tileBytes) / bpp;
    uint32_t inTileX     = MortonX(mortonIdx);
    uint32_t inTileY     = MortonY(mortonIdx);
    uint32_t tileCol     = tileIndex % tilesPerRow;
    uint32_t tileRow     = tileIndex / tilesPerRow;
    return ((tileRow * tileH + inTileY) * pitch) + ((tileCol * tileW + inTileX) * bpp);
}
```

The current RT readback path (§6.2) performs a `CopyLinearWithPitch` memcpy without the tile untransform. This is a known limitation: the correct behavior for tiled RTs read by the CPU via the Contiguous window would require running `TiledToLinear` on each pixel. In practice most CPU reads of RTs (screenshots, save-game thumbnails) tolerate the swizzled data or ignore it, and the GPU→GPU fast path (§7.2) avoids the readback entirely for the common shadow/reflection/post-process case.

### 9.5 Z Compression

The NV2A zcomp is a tag-based lossless scheme: each 4×4 depth block may be stored as a single representative Z value when the block is uniform (most commonly after a depth clear). The compressed store lives in a `PFB_ZCOMP` region; the full-precision Z surface lives in the associated `PFB_TILE` region.

From the host emulator's perspective, the DX11 depth buffer is always maintained in fully-resolved form — DX11 manages its own internal depth compression transparently. A per-ZCOMP-region dirty flag tracks whether a resolve readback is needed before CPU reads:

```cpp
bool gZCompDirty[8];

void ResolveZComp(uint32_t zPhysBase) {
    int idx = FindZCompRegionForPhys(zPhysBase);
    if (idx < 0 || !gZCompDirty[idx]) return;
    ReadbackSurfaceFromAlias(GetDepthSurface(zPhysBase));
    gZCompDirty[idx] = false;
}
```

---

## 10. Host Resource Aliases — Texture Pool and RT Cache

### 10.1 Two Distinct Resource Categories

**Render targets** are created by the GPU and never uploaded from the CPU. They are tracked in `g_PgraphRTCache` keyed by `{physBase, DXGI_FORMAT, width, height}`. Each host `Texture2D` is created with `D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE` so it can serve as both RTV and SRV without re-creation.

**Textures** (static art, dynamic CPU-written surfaces, video frames) are uploaded from Xbox RAM on demand. They live in a separate texture pool keyed by `{physBase, format, dimensions, mipCount}`. Each entry is created on first access and invalidated by the `s_TextureDirtyBitmap` check on every bind.

The two caches are currently separate structures (a unified cache is a planned improvement). The RT-as-texture fast path (§7.2) bridges them for the GPU→GPU surface feedback case.

### 10.2 Texture Upload — Deswizzle Compute Shader

For swizzled textures, a pre-draw CS dispatched per dirty surface reads raw bytes from a **dedicated staging buffer** (uploaded from Xbox RAM for the texture's address range) and writes decoded linear texels into the pool's `Texture2D UAV`. The Morton decode and channel reorder occur in the same dispatch:

```hlsl
// CxbxUnswizzleCS.hlsl — deswizzle a POT Xbox texture into a host Texture2D pool entry
ByteAddressBuffer   g_SrcBuffer : register(t0);  // Staging buffer with raw Xbox texture data
RWTexture2D<uint>   g_DstTexture : register(u0);  // Texture pool entry UAV

cbuffer Params : register(b0) {
    uint g_Width;
    uint g_Height;
    uint g_Depth;
    uint g_BytesPerPixel;
    uint g_RowPitch;
    uint g_SlicePitch;
}

[numthreads(8, 8, 1)]
void main(uint3 tid : SV_DispatchThreadID) {
    if (tid.x >= g_Width || tid.y >= g_Height) return;
    uint swizzledOffset = Swizzle2D(tid.x, tid.y, g_Width, g_Height, g_BytesPerPixel);
    uint rawData = g_SrcBuffer.Load(swizzledOffset);
    uint linearOffset = tid.y * g_RowPitch + tid.x * g_BytesPerPixel;
    g_DstTexture[tid.xy] = rawData;
}
```

**Key difference from proposed design:** The CS reads from a per-texture staging buffer (`g_SrcBuffer`), NOT from the main mirror buffer (`s_pMirrorBuf`). This is because the texture data must first be uploaded to a staging buffer sized to the texture (from Xbox RAM), then the CS deswizzles into the `Texture2D` UAV.

**Fallback:** When the texture format does not support UAV writes (e.g., certain compressed formats), the CPU-side `swizzle_rect()` / `unswizzle_rect()` path is used instead, writing into a staging texture that is then `CopyResource`'d to the final `Texture2D`.

The dirty check that gates dispatch uses `CxbxPageTrackerIsTextureDirty()` — a DWORD-granularity bitmap scan.

For DXT textures, Morton operates at the block level (8 or 16 bytes per block). The CS covers one block per thread group element, reading the full block payload.

For linear textures (NPOT, pitch from TEXCTL1), no Morton decode is needed. Upload is a direct `UpdateSubresource` row by row, handling Xbox pitch vs. host pitch differences.

### 10.3 Palette Resolution from PGRAPH

For P8 (palettized) textures, the palette is resolved from `TEXPALETTE0 + stage*4` via the DMA context A/B base address — no HLE patch involvement. With the struct-based PGRAPH state, the register is read directly from `pg->regs[]`:

```cpp
uint32_t* NV2AResolvePaletteAddress(PGRAPHState* pg, uint32_t stage) {
    // Read PGRAPH TEXPALETTE register from pg->regs[]
    uint32_t palReg  = pg->regs[(NV_PGRAPH_TEXPALETTE0 + stage * 4) >> 2];
    uint32_t ctxSel  = (palReg >> 1) & 1u;  // DMA context A or B
    uint32_t offset  = palReg & ~0xFu;       // byte offset within DMA region
    // RAMHT lookup for the DMA context (lives in PRAMIN, aliases Xbox RAM)
    uint32_t dmaBase = ResolveDMAContextFromPRAMIN(ctxSel);
    return reinterpret_cast<uint32_t*>((uint8_t*)CONTIG_BASE + dmaBase + offset);
}
```

The palette is expanded to RGBA8 during the texture upload pass (CPU-side or via CS).

### 10.4 RT Alias Lifecycle

```
NV097_SET_SURFACE_COLOR_OFFSET(X):
    → GetOrCreateRT(X, pitch, w, h, fmt)
    → if new: CreateTexture2D(BIND_RENDER_TARGET | BIND_SHADER_RESOURCE)
    → CxbxPageTrackerMarkGPUDirty(X, pitch * h)
    → OMSetRenderTargets(pRTV)

[draw, blit, or clear renders into the RT]

Later NV097 texture bind with TEXOFFSET == X:
    → CxbxLookupPgraphRTByOffset(X)  — found in g_PgraphRTCache
    → OMSetRenderTargets(0, nullptr, nullptr)  — unbind RTV first (DX11 hazard rule)
    → CreateShaderResourceView + PSSetShaderResources
    → No CPU round-trip. GPU→GPU zero-copy.

CPU reads from address in [X, X + pitch * h):
    → PAGE_NOACCESS VEH fires
    → TryEnterCriticalSection with puller
    → CopyResource(stagingTex, hostTex)
    → Map + row-by-row memcpy to s_pXboxRAM
    → Restore PAGE_READWRITE for entire RT range
    → Clear GpuDirtyBitmap for RT pages
```

### 10.5 Known Limitations

**Mip tail packing.** The NV2A packs mip levels smaller than one block into a hardware-specific "tail" region with alignment rules that differ from naive `offset += mipSlicePitch` advancement. The current upload path uses naive advancement, which produces incorrect data for the lowest mip levels of high-mipcount textures. Impact is usually invisible (lowest mips are rarely sampled).

**RT cache eviction.** `g_PgraphRTCache` has no eviction policy and grows unbounded. Long-running games with many RT switches leak host `Texture2D` objects until device reset.

---

## 11. NV2A Rendering Operations

### 11.1 Three Vertex Submission Modes

**Draw Arrays** (`NV097_DRAW_ARRAYS`): Vertex data lives in Xbox RAM at addresses in `pg->vertex_attributes[i].offset` with per-attribute stride and format. The VS fetches data directly from `s_pGpuMem` using `SV_VertexID` arithmetic (§12).

**Inline Array** (`NV097_INLINE_ARRAY`): All-attributes-per-vertex packed data is submitted inline through the pushbuffer into `pg->inline_array[]`. Uploaded to a DYNAMIC D3D11 buffer per draw and discarded after.

**Inline Buffer** (via `NV097_SET_VERTEX_DATA2F_M` / `NV097_SET_VERTEX_DATA4F_M` + `NV097_ARRAY_ELEMENT16/32`): Per-vertex float4 values are submitted per-attribute during a Begin/End sequence into `pg->vertex_attributes[i].inline_buffer[]`. Also uploaded per draw to a DYNAMIC buffer.

### 11.2 Draw Call Flow (pgraph_draw / D3D11_draw)

On `NV097_SET_BEGIN_END(0)`, the puller thread calls `pgraph_draw()` which dispatches to the `D3D11_draw` callback:

1. **NOP check** — skip if no color write AND no depth test AND no stencil test.
2. **Surface update** — bind/create host RTs from `pg->surface_color.offset` and `pg->surface_zeta.offset`.
3. **Page flush** — `CxbxPageTrackerFlushToGPU()` for CPU-written vertex pages (once per frame, gated by `s_bFirstFlushOfFrame`).
4. **Texture bind** — per stage 0–3: check `TEXCTL0` enable bit OR `SHADERPROG` non-NONE mode (the second condition handles point sprites using stage 3 without TEXCTL0 set). RT-as-texture fast path if TEXOFFSET matches a known RT.
5. **Shader bind** — VS: JIT-compiled HLSL from `pg->program_data[]`, or interpreter fallback. PS: JIT-compiled combiner HLSL from `pg->regs[]` combiner topology hash, or interpreter fallback.
6. **Upload state** — `pg->vsh_constants[]` → cbuffer b0; `pg->program_data[]` → `g_XFPR` at t5; `pg->regs[]` → `g_PGRegs` at t12 (gated by generation counter).
7. **Pipeline state objects** — blend, depth/stencil, rasterizer decoded from `pg->regs[]` entries for `CONTROL_0/1/2`, `BLEND`, `SETUPRASTER`; looked up from PSO cache keyed on packed state bits. (Most Xbox titles use fewer than 20 distinct PSO combinations per frame.)
8. **Viewport/scissor** — from `NV2ASurfaceState` clip fields and `ZCLIPMIN`/`ZCLIPMAX`.
9. **Draw dispatch** — `Draw(hostVertexCount, 0)` with null IA buffers (§12).
10. **Post-draw dirty** — `CxbxPageTrackerMarkGPUDirty` for color and zeta surface ranges.

### 11.3 Blitting

```cpp
void ExecuteBlit(uint32_t srcPhys, uint32_t dstPhys, ...) {
    CxbxPageTrackerFlushToGPU();
    auto& src = GetOrCreateSurface(srcPhys, ...);
    auto& dst = GetOrCreateSurface(dstPhys, ...);
    EnsureSurfaceUpToDate(src);  // deswizzle CS if needed

    D3D11_BOX srcBox{ srcX, srcY, 0, srcX+w, srcY+h, 1 };
    gD3DCtx->CopySubresourceRegion(dst.texture, 0, dstX, dstY, 0,
                                   src.texture, 0, &srcBox);
    CxbxPageTrackerMarkGPUDirty(dstPhys, dstPitch * (dstY + h));
}
```

For `NV_SCALED_IMAGE_FROM_MEMORY` (scaled blit from linear guest RAM to RT), the source is uploaded directly to a staging texture and scaled via a DX11 compute or PS pass before writing to the destination alias.

### 11.4 Clear Operations

```cpp
void ExecuteClear(uint32_t colorPhys, uint32_t zetaPhys, uint32_t color,
                  float depth, uint8_t stencil, uint32_t mask) {
    if (mask & NV097_CLEAR_SURFACE_COLOR) {
        auto& rt = GetOrCreateRT(colorPhys, ...);
        float c[4] = { ((color>>16)&0xFF)/255.f, ((color>>8)&0xFF)/255.f,
                       ((color    )&0xFF)/255.f, ((color>>24)&0xFF)/255.f };
        gD3DCtx->ClearRenderTargetView(rt.rtv, c);
        CxbxPageTrackerMarkGPUDirty(colorPhys, rt.pitch * rt.height);
    }
    if (mask & (NV097_CLEAR_SURFACE_Z | NV097_CLEAR_SURFACE_STENCIL)) {
        auto& ds = GetOrCreateDepthSurface(zetaPhys, ...);
        UINT flags = ((mask & NV097_CLEAR_SURFACE_Z)       ? D3D11_CLEAR_DEPTH   : 0)
                   | ((mask & NV097_CLEAR_SURFACE_STENCIL) ? D3D11_CLEAR_STENCIL : 0);
        gD3DCtx->ClearDepthStencilView(ds.dsv, flags, depth, stencil);
        int idx = FindZCompRegionForPhys(zetaPhys);
        if (idx >= 0) gZCompDirty[idx] = true;
        CxbxPageTrackerMarkGPUDirty(zetaPhys, ds.pitch * ds.height);
    }
}
```

---

## 12. Vertex Fetch — IA Bypass Architecture

### 12.1 SV_VertexID Approach

The DX11 Input Assembler is bypassed entirely. Every draw is issued as:

```cpp
gD3DCtx->IASetInputLayout(nullptr);
gD3DCtx->IASetVertexBuffers(0, 0, nullptr, nullptr, nullptr);
gD3DCtx->IASetIndexBuffer(nullptr, DXGI_FORMAT_UNKNOWN, 0);
gD3DCtx->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
gD3DCtx->Draw(hostVertexCount, 0);
```

The VS receives only `SV_VertexID` and fetches all attributes from `s_pGpuMem`. No CPU-side attribute conversion and no pre-draw index expansion step are required, except for topology remapping (see §12.3).

### 12.2 Typed SRV Views for Standard Formats

For attributes with DX11 typed SRV equivalents (`F32`, `F16`, `SNORM16`, `UBYTE_OGL`, `UBYTE_D3D`), typed `Buffer<T>` SRVs are created over the XboxRAM region at the PGRAPH-registered base offset, cached by `(baseOffset, stride, dxgiFormat)`:

```hlsl
Buffer<float4> VtxStream_F32   : register(t1);  // F32_4 / F32_3 / F32_2 / F32_1
Buffer<float4> VtxStream_SNORM : register(t2);  // S1 (SNORM16)
Buffer<uint>   VtxStream_UB    : register(t3);  // UB_OGL / UB_D3D
```

For the `CMP` format (11.11.10 packed signed normal — no typed SRV equivalent), five ALU instructions decode inline in the VS:

```hlsl
// NV2A vertex type 6: X in bits 10:0, Y in bits 21:11, Z in bits 31:22
float3 DecodeCMP(uint offset) {
    uint raw = XboxRAM.Load(offset);
    int  x   = (int)(raw << 21u) >> 21;   // sign-extend 11 bits from 10:0
    int  y   = (int)(raw << 10u) >> 21;   // sign-extend 11 bits from 21:11
                                           // NOTE: left shift is 10, not 11
    int  z   = (int)(raw       ) >> 22;   // sign-extend 10 bits from 31:22
    return float3(x / 1023.0f, y / 1023.0f, z / 511.0f);
}
```

The Y-component shift must be **10** (to place bit 21 at the int32 MSB before the arithmetic right shift). Using 11 would discard bit 21 and produce incorrect Y values.

### 12.3 Topology Conversion Inline in VS

Quad lists and triangle fans are remapped to `TRIANGLELIST` entirely within the VS using `SV_VertexID` arithmetic. No pre-draw compute shader and no scratch index buffer:

```hlsl
// Quad list: 4 Xbox verts → 6 host verts (2 triangles per quad)
uint QuadListIndex(uint vertId) {
    static const uint lut[6] = { 0u, 1u, 2u, 0u, 2u, 3u };
    return (vertId / 6u) * 4u + lut[vertId % 6u];
}
// Draw(numQuads * 6, 0)

// Triangle fan: (N-2)*3 host verts; apex always at vertex 0
uint TriFanIndex(uint vertId) {
    uint tri   = vertId / 3u;
    uint local = vertId % 3u;
    return local == 0u ? 0u : (tri + local);
}
// Draw((numVerts - 2) * 3, 0)

// Indexed quads: resolve through Xbox index buffer after topology remap
uint IndexedQuadListIndex(uint vertId, uint xboxIndexBase, bool is16bit) {
    uint logicalIdx = QuadListIndex(vertId);
    if (is16bit) {
        uint word  = XboxRAM.Load((xboxIndexBase + logicalIdx * 2u) & ~3u);
        uint shift = (logicalIdx & 1u) * 16u;
        return (word >> shift) & 0xFFFFu;
    }
    return XboxRAM.Load(xboxIndexBase + logicalIdx * 4u);
}
```

A separate VS permutation is compiled per topology (`PRIM_TYPE` compile-time `#define`) so unused branches compile away entirely. The GPU's post-transform cache handles the repeated vertex 0 fetch in triangle fans.

---

## 13. Vertex and Pixel Shader Emulation

### 13.1 NV2A Vertex Shader ISA

The NV2A VSH is a 136-instruction RISC unit. Each 128-bit instruction encodes up to two parallel operations:

- **ILU**: RCP, RSQ, LOG, EXP, LIT, SGE, SLT
- **MAC**: MOV, MUL, ADD, MAD, DP3, DP4, DPH, DST, MIN, MAX, SFL, SGE, SLT, ARL

Register file: 16 input attributes (v[0]–v[15]), **192** constants (c[0]–c[191], from the XFCTX RAM — `pg->vsh_constants`), 12 temporaries (r[0]–r[11]), 2 address registers. Output: 13 registers (oPos, oD0/oD1 diffuse/specular, oT0–oT3 texcoords, oB0/oB1 back-face colors, oFog, oPts).

#### Interpreter

Evaluates the instruction stream at draw time via `g_XFPR` (StructuredBuffer<uint4> at t5, uploaded from `pg->program_data[]`) and `cbuffer b0` (192 × float4 constants from `pg->vsh_constants`). Used as a reference and as a fallback for JIT-unsupported patterns.

#### JIT Compiler

Translates the NV2A VSH stream to HLSL VS at shader-upload time (`TranslateToHLSL()`). Compiled output is cached in memory by instruction-stream hash (rapidhash) and additionally on disk to avoid recompilation between sessions:

```cpp
void OnSetTransformProgram(const uint32_t* instrs, uint32_t count) {
    uint64_t hash = rapidhash(instrs, count * 16);
    if (auto it = g_VSHCache.find(hash); it != g_VSHCache.end()) {
        gD3DCtx->VSSetShader(it->second.Get(), nullptr, 0);
        return;
    }
    std::string hlsl = TranslateNV2AVSH(instrs, count);
    // D3DCompile → CreateVertexShader → cache (memory + disk)
}
```

### 13.2 NV2A Register Combiner (Pixel Pipeline)

Up to 8 general combiner stages plus a final combiner, parameterized by PGRAPH registers 0x1880–0x1948 (offsets within the PGRAPH block). All register reads in HLSL use `g_PGRegs` (the PGRAPH block `StructuredBuffer<uint>` at t12, uploaded from `pg->regs[]` via `Map/Unmap`), indexed by block-relative offset divided by 4:

```
NV_PGRAPH_COMBINEFACTOR0 + i*4  (0x1880)  C0[i] per stage
NV_PGRAPH_COMBINEFACTOR1 + i*4  (0x18A0)  C1[i] per stage
NV_PGRAPH_COMBINEALPHAI0 + i*4  (0x18C0)  Alpha input selectors
NV_PGRAPH_COMBINEALPHAO0 + i*4  (0x18E0)  Alpha output config
NV_PGRAPH_COMBINECOLORI0 + i*4  (0x1900)  RGB input selectors
NV_PGRAPH_COMBINECOLORO0 + i*4  (0x1920)  RGB output config
NV_PGRAPH_COMBINECTL          (0x1940)    Stage count, MUX-MSB, UniqueC0/C1 flags
NV_PGRAPH_COMBINESPECFOG0     (0x1944)    Final combiner inputs ABCD
NV_PGRAPH_COMBINESPECFOG1     (0x1948)    Final combiner inputs EFG
```

#### Known Correctness Issues in SM3.0 Interpreter (Fixed in SM5 Port)

- **Texture register write**: A `max()` clamp caused all stages to write to T3. Fixed with direct indexed write: `Regs[PS_REGISTER_T0 + stage] = val`.
- **R0 initialization**: R0.rgb started as `float4(1,1,1,x)` (white) instead of zero. Per NV2A hardware: R0.rgb = 0, R0.a = T0.a.
- **PS_CHANNEL_RGB passthrough**: When the alpha-channel-select bit is clear, the full float4 should pass through unchanged. The SM3.0 version replaced `.a` with `.b` (`.rgbb`).
- **DOT_RFLCT_DIFF**: Incorrectly called `reflect()`. The DIFF mode uses the reconstructed normal directly as a cubemap direction; only the SPEC mode reflects.

#### JIT Compiler

Combiner topology (stage count, dot/mux flags, texture modes) is JIT-compiled to a native HLSL PS, cached by combiner state hash (`PixelShaderCache`). Cache lookup avoids hash + mutex + map lookup when combiner state is unchanged between consecutive draws. Dynamic per-frame data (C0/C1 constants, fog parameters, bump matrices) remains in `g_PGRegs` (t12) and is re-uploaded every draw (gated by generation counter) rather than baked into the shader. Eight `NUM_STAGES` variants (1–8) are pre-compiled with the loop bound as a `#define`, eliminating dead stages at compile time.

---

## 14. HLSL SM5.0 Format Conversion Shaders

### 14.1 Swizzled Texture Deswizzle (Compute → Texture2D UAV)

As described in §10.2, the deswizzle CS reads raw bytes from a per-texture staging `ByteAddressBuffer` and writes decoded texels to a `Texture2D` UAV. Before the draw, the UAV binding is released and the texture is rebound as an SRV:

```hlsl
// Current implementation (CxbxUnswizzleCS.hlsl)
ByteAddressBuffer   g_SrcBuffer  : register(t0);  // staging buffer (per-texture upload)
RWTexture2D<uint>   g_DstTexture : register(u0);  // texture pool entry UAV

cbuffer UnswizzleParams : register(b0) {
    uint g_Width;
    uint g_Height;
    uint g_Depth;
    uint g_BytesPerPixel;
    uint g_RowPitch;
    uint g_SlicePitch;
}

[numthreads(8, 8, 1)]
void main(uint3 tid : SV_DispatchThreadID) {
    if (tid.x >= g_Width || tid.y >= g_Height) return;
    uint swizzledOffset = Swizzle2D(tid.x, tid.y, g_Width, g_Height, g_BytesPerPixel);
    uint rawData = g_SrcBuffer.Load(swizzledOffset);
    g_DstTexture[tid.xy] = rawData;
}
```

A BGRA variant handles A8R8G8B8 textures that need channel reordering during deswizzle. Formats that do not support UAV writes fall back to CPU-side `unswizzle_rect()`.

**Proposed improvement:** Read directly from the combined `g_GpuMem` buffer using a `SurfaceBase` offset, eliminating the per-texture staging upload step and reducing CPU→GPU bandwidth.

### 14.2 CMP Vertex Attribute Decode (Inline in VS)

No intermediate buffer needed; five ALU instructions per affected attribute:

```hlsl
// 11.11.10 signed packed normal — NV2A vertex attribute type 6
// Bits 10:0 = X (11-bit signed), bits 21:11 = Y (11-bit signed),
// bits 31:22 = Z (10-bit signed)
float3 DecodeCMP(uint offset) {
    uint raw = XboxRAM.Load(offset);
    int  x   = (int)(raw << 21u) >> 21;  // left 21, right 21 → 11-bit sign-extend
    int  y   = (int)(raw << 10u) >> 21;  // left 10 to put bit21 at MSB, right 21
    int  z   = (int)(raw       ) >> 22;  // right 22 → 10-bit sign-extend
    return float3(x / 1023.0f, y / 1023.0f, z / 511.0f);
}
```

The Y component requires a left shift of **10** (not 11). `<< 10` places bit 21 at the int32 MSB; `<< 11` would discard bit 21 entirely.

### 14.3 YUY2/UYVY Conversion

YUV surfaces are CPU-converted to ARGB8 during the texture upload pass using row converters (`YUY2ToARGBRow_C` / `UYVYToARGBRow_C`). The converted output is uploaded as a standard `BGRA8_UNORM` `Texture2D`. This path is triggered by `EmuXBFormatRequiresConversion()` returning true for YUV format entries.

### 14.4 PFB Retile as Compute (Optional, for CPU-Written Tiled Surfaces)

For surfaces in PFB_TILE regions that are CPU-written and then GPU-read via the tiled aperture, a retile CS converts the linear CPU data to tiled layout in a single dispatch. This avoids a per-page CPU loop:

```hlsl
// cs_5_0 — retile linear CPU data into tiled layout for GPU consumption
ByteAddressBuffer    LinearSrc : register(t0);
RWByteAddressBuffer  TiledDst  : register(u0);  // separate UAV-capable scratch buffer
cbuffer Params : register(b0) { uint Base; uint Pitch; uint Bpp; }

[numthreads(8, 8, 1)]
void Retile(uint3 tid : SV_DispatchThreadID) {
    uint linOff   = tid.y * Pitch + tid.x * Bpp;
    uint tiledOff = LinearToTiledOffset(tid.x, tid.y, Pitch, Bpp);
    TiledDst.Store(Base + tiledOff, LinearSrc.Load(linOff));
}
```

---

## 15. Tile Register Interception

When guest code writes MMIO tile registers (both PFB and PGRAPH must be kept in sync — XDK setup code writes both), the PGRAPH method handler commits the value to `pg->regs[]` in the `NV2AState` struct, and a side-effect handler refreshes the fast `gTileRegions[]` lookup table and invalidates any cached surface aliases whose physical range now falls under changed tiling semantics:

```cpp
void OnTileRegisterWrite(uint32_t mmioOffset, uint32_t value) {
    // Commit to NV2AState struct — pg->regs[] is the ground truth
    pg->regs[mmioOffset >> 2] = value;

    bool isPFB = (mmioOffset >= 0x00100240 && mmioOffset < 0x00100240 + 8*0x10);
    int  n     = isPFB ? (mmioOffset - 0x00100240) / 0x10
                       : (mmioOffset - 0x00400904) / 0x10;

    switch (mmioOffset & 0xCu) {
    case 0x0:  // TILE[N]: base + enable flag
        gTileRegions[n].base    = value & 0xFFFF0000u;
        gTileRegions[n].enabled = (value & 0x1u) != 0;
        break;
    case 0x4:  // TLIMIT[N]
        gTileRegions[n].limit   = (value & 0xFFFF0000u) | 0xFFFFu;
        break;
    case 0x8:  // TSIZE[N]: pitch in 256-byte units stored in bits 31:16
        gTileRegions[n].pitch   = ((value >> 16) & 0xFFFFu) * 256u;
        break;
    }

    // Invalidate host surface aliases covering this tile region range
    for (auto& [physBase, entry] : g_ResourceCache) {
        if (physBase >= gTileRegions[n].base &&
            physBase <= gTileRegions[n].limit) {
            entry.aliasValid = false;
        }
    }
}
```

Because PFB registers are not currently uploaded to the GPU (§2.4), the `gTileRegions[]` fast-lookup table is the primary source for tiling decisions on both the CPU (upload path selection, surface invalidation) and shader compilation (tiling-aware texture fetch). In the proposed combined-buffer design, the PFB block would also be uploaded to the GPU as part of `s_pGpuMem`, enabling shaders to read tile region configuration directly.

---

## 16. Complete Coherency State Machine

Each physical 4KB page has state derived from two independently maintained bits:

```
texDirty  gpuDirty  State        Meaning
────────  ────────  ───────────  ────────────────────────────────────────────
0         0         CLEAN        Both sides coherent. No action.
1         0         CPU_DIRTY    CPU wrote. s_pGpuMem / texture pool stale.
                                 Flush to GPU (Map/deswizzle CS) before read.
0         1         GPU_DIRTY    GPU wrote. s_pXboxRAM stale. NOACCESS armed.
                                 Readback on CPU read (VEH path).
1         1         CONFLICT     Both wrote. CPU wins.
                                 gpuDirty cleared at next flush; GPU data discarded.
```

Transitions at six events:

```
Event                                → texDirty  gpuDirty  Protection   Action
───────────────────────────────────────────────────────────────────────────────
GetWriteWatch fires for page P       → set       clear     unchanged    SetTextureDirty; clear gpuDirty
CxbxPageTrackerFlushToGPU: page done → clear     ─         unchanged    AND-NOT clear; upload done
GPU renders to range [B, B+size)     → ─         set       NOACCESS     MarkGPUDirty; VirtualProtect
VEH: CPU read from gpuDirty page     → ─         clear     READWRITE    CopyResource→Map→memcpy entire RT
VEH: CPU write to gpuDirty page      → ─         clear     READWRITE    GPU data discarded; no readback
Deswizzle CS dispatched for surface  → clear     ─         unchanged    Pool entry updated; texDirty cleared
```

**Key architectural properties (current implementation):**

- `s_TextureDirtyBitmap` and `s_GpuDirtyBitmap` are separate `uint32_t[512]` arrays. In the current single-threaded access model (flush runs on the same thread as draw), non-atomic operations are sufficient.
- The flush is gated by `s_bFirstFlushOfFrame`, ensuring only one `GetWriteWatch` + buffer upload per frame. Mid-frame `CxbxPageTrackerFlushGPUDirtyToMirror()` handles the VB-aliases-RT case with targeted NO_OVERWRITE uploads.
- RT readbacks use `PAGE_NOACCESS` (not `PAGE_GUARD`). Access is restored on the entire RT range before the `CopyResource` + `memcpy` to prevent nested faults. Staging textures are created per-readback (not cached).
- The RT-as-texture fast path (§7.2) eliminates any CPU round-trip for GPU→GPU surface feedback: shadow maps, reflections, post-process intermediates all stay entirely on the host GPU.
- Vertex data flush uses `Map(WRITE_DISCARD)` for bulk (>25% pages dirty, first flush) or `Map(WRITE_NO_OVERWRITE)` for incremental. Adjacent dirty pages are coalesced into contiguous runs before `memcpy`.
- PGRAPH register state is uploaded as a separate `StructuredBuffer<uint>` at t12 via `Map/Unmap`, gated by a generation counter. PFB registers are consumed CPU-side only. VS program instructions (`g_XFPR`) and VS constants (`cbuffer b0`) are uploaded separately.
- Total XboxRAM coherency metadata: 4 × 2048 = 8192 bytes (GpuDirty, TextureDirty, TiledCommitted, IdentityAlloc), fully resident in L1 cache on any modern CPU.

**Proposed improvements (not yet implemented):**

- Combined `s_pGpuMem` buffer (§2.2) with appended MMIO blocks for single-binding shader access.
- Flat 16 MiB `s_pNV2AMMIO` allocation replacing `NV2AState` struct (§2.3) for memcpy-based upload.
- Cache-line-aligned bitmaps with `_InterlockedOr8`/`_InterlockedAnd8` for lock-free multi-threaded access.
- AVX2 bulk bitmap scan for large texture ranges (§4.3).
- Per-RT cached staging textures instead of per-readback allocation (§6.2).
- Deswizzle CS reading directly from the mirror buffer instead of a per-texture staging buffer (§10.2).
