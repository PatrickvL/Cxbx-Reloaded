# NV2A GPU Emulation — DX11 Architecture

## Overview

This document describes a GPU-driven architecture for emulating the Xbox NV2A graphics processor using Direct3D 11. The guiding principle is to push as much of the NV2A emulation workload onto the host GPU as possible, minimising CPU involvement in the per-draw rendering data path.

The NV2A is an Nvidia NV20-family GPU with a fixed 64MB unified memory space shared by textures, vertex buffers, index buffers, push buffers, and the PGRAPH register file. This flat address model maps cleanly onto a single GPU buffer, which becomes the architectural anchor for the entire approach.

---

## Memory Model

### XboxRAM Buffer

The entire 64MB Xbox memory space is represented as a single `ID3D11Buffer` on the host GPU:

```
XboxRAM: ID3D11Buffer  64MB
         BIND_SHADER_RESOURCE
         Usage: DYNAMIC, CPUAccessFlags: WRITE
         MiscFlags: BUFFER_ALLOW_RAW_VIEWS
         Format: DXGI_FORMAT_R32_TYPELESS (allows multiple typed SRV views)
```

> **Implementation note:** The current codebase (`Backend_D3D11_PageTracker.cpp`) uses
> `D3D11_USAGE_DYNAMIC` with `Map`/`Unmap` rather than `UpdateSubresource`. UAV access
> is not bound; deswizzle compute shaders write to `Texture2D` UAVs, not to this buffer.

The CPU maintains a shadow copy of this buffer in system memory. Whenever the Xbox CPU writes to memory, the write goes to the shadow and the affected byte range is marked dirty. Before each draw call, dirty ranges are flushed with `UpdateSubresource`. This is the only CPU→GPU data transfer in the rendering path.

```cpp
void OnXboxMemoryWrite(uint32_t addr, uint32_t value) {
    shadow[addr / 4] = value;
    dirty_lo = min(dirty_lo, addr);
    dirty_hi = max(dirty_hi, addr + 4);
}

void FlushDirtyRanges(ID3D11DeviceContext* ctx) {
    if (dirty_lo >= dirty_hi) return;
    D3D11_BOX box = { dirty_lo, 0, 0, dirty_hi, 1, 1 };
    ctx->UpdateSubresource(xboxRAM_buf, 0, &box,
                           shadow + dirty_lo / 4, 0, 0);
    dirty_lo = UINT_MAX;
    dirty_hi = 0;
}
```

### PGRAPH as a ByteAddressBuffer

The NV2A PGRAPH register file lives at a fixed offset within the 64MB space. All combiner state, vertex array descriptors, texture descriptors, and draw parameters are PGRAPH registers with documented byte offsets. Shaders access them directly via `ByteAddressBuffer.Load` at those offsets — no CPU-side struct construction or field extraction is needed.

```hlsl
ByteAddressBuffer XboxRAM : register(t0);

// Known PGRAPH offsets (from nv2a_regs.h — PGRAPH MMIO space, NOT NV097 method offsets)
static const uint NV_PGRAPH_COMBINEFACTOR0  = 0x1880u;  // C0 per stage
static const uint NV_PGRAPH_COMBINEFACTOR1  = 0x18A0u;  // C1 per stage
static const uint NV_PGRAPH_COMBINEALPHAI0  = 0x18C0u;  // Alpha inputs
static const uint NV_PGRAPH_COMBINEALPHAO0  = 0x18E0u;  // Alpha outputs
static const uint NV_PGRAPH_COMBINECOLORI0  = 0x1900u;  // RGB inputs
static const uint NV_PGRAPH_COMBINECOLORO0  = 0x1920u;  // RGB outputs
static const uint NV_PGRAPH_COMBINECTL      = 0x1940u;  // Combiner control
static const uint NV_PGRAPH_COMBINESPECFOG0 = 0x1944u;  // Final combiner ABCD
static const uint NV_PGRAPH_COMBINESPECFOG1 = 0x1948u;  // Final combiner EFG

uint  RGBInputs(uint stage)  { return XboxRAM.Load(PGRAPH_BASE + NV_PGRAPH_COMBINECOLORI0 + stage * 4u); }
uint  AlphaInputs(uint stage){ return XboxRAM.Load(PGRAPH_BASE + NV_PGRAPH_COMBINEALPHAI0 + stage * 4u); }
uint  RGBOutputs(uint stage) { return XboxRAM.Load(PGRAPH_BASE + NV_PGRAPH_COMBINECOLORO0 + stage * 4u); }
uint  AlphaOutputs(uint stage){ return XboxRAM.Load(PGRAPH_BASE + NV_PGRAPH_COMBINEALPHAO0 + stage * 4u); }
```

Color constants stored as packed ARGB bytes are decoded inline:

```hlsl
float4 LoadColorConstant(uint baseOffset, uint stage) {
    uint packed = XboxRAM.Load(baseOffset + stage * 4u);
    return float4(
        ((packed >> 16u) & 0xFFu) / 255.0f,
        ((packed >>  8u) & 0xFFu) / 255.0f,
        ((packed       ) & 0xFFu) / 255.0f,
        ((packed >> 24u) & 0xFFu) / 255.0f
    );
}
```

### No HostRAM Buffer Required

An earlier design considered a second 64MB `HostRAM` buffer at identity-mapped offsets to hold converted vertex data. Analysis of all NV2A vertex formats showed this is unnecessary:

- **F32, F16, S1 (SNORM16), UB\_OGL, UB\_D3D, S32K** — all have exact DX11 typed SRV equivalents. Multiple typed SRV views can be created over the same `XboxRAM` buffer resource using `DXGI_FORMAT_R32_TYPELESS`, allowing the GPU's vertex fetch unit to decode them natively with zero shader cost.
- **CMP (11.11.10 packed normal)** — no typed SRV equivalent exists, but the format decodes in four ALU instructions inline in the vertex shader. No intermediate buffer is needed.
- **Quad and triangle fan indices** — resolved by arithmetic in the vertex shader from `SV_VertexID` with no index buffer at all (see Topology Conversion).
- **Swizzled textures** — the one remaining case requiring pre-draw work (see Texture Handling).

---

## Vertex Fetch

The DX11 Input Assembler is bypassed entirely. Every draw is issued as `Draw(hostVertexCount, 0)` with null vertex and index buffers. The vertex shader receives only `SV_VertexID` and performs all attribute fetching manually from `XboxRAM`.

### Typed SRV Vertex Fetch

For formats with DX11 equivalents, a typed `Buffer<T>` SRV is created on demand over the XboxRAM region at the PGRAPH-registered base address:

```hlsl
// Declared per active vertex stream, SRV created at draw time
Buffer<float4> VtxStream_F32  : register(t1);
Buffer<float4> VtxStream_SNORM: register(t2);
Buffer<uint>   VtxStream_UB   : register(t3);
```

SRV creation is keyed on `(baseOffset, strideBytes, dxgiFormat)` and cached. The CPU creates a new view only when the format or base address changes.

### CMP Inline Decode

```hlsl
float3 DecodeCMP(uint offset) {
    uint raw = XboxRAM.Load(offset);
    int x = (int)(raw << 21u) >> 21u;  // sign-extend 11 bits
    int y = (int)(raw << 10u) >> 21u;  // sign-extend 11 bits
    int z = (int)(raw       ) >> 22u;  // sign-extend 10 bits
    return float3(x / 1023.0f, y / 1023.0f, z / 511.0f);
}
```

Four bytes in, float3 in registers. The compiler folds the three divisions into multiplies. The cost is approximately five ALU instructions per CMP attribute per vertex — negligible against Xbox-era vertex counts.

### Unified Fetch Function

```hlsl
float4 FetchAttr(uint xboxVtxIdx, uint base, uint stride, uint fmt) {
    uint offset = base + xboxVtxIdx * stride;
    switch (fmt) {
        case NV2A_FMT_F32_4:  return asfloat(XboxRAM.Load4(offset));
        case NV2A_FMT_F32_2:  return float4(asfloat(XboxRAM.Load2(offset)), 0, 1);
        case NV2A_FMT_CMP:    return float4(DecodeCMP(offset), 1);
        // F16, S1, UB paths use typed SRV Buffer<T>.Load(xboxVtxIdx)
        default:              return float4(0, 0, 0, 1);
    }
}
```

---

## Topology Conversion

NV2A supports quad lists and triangle fans, which have no DX11 primitive type equivalent. Rather than converting indices in a pre-draw compute shader and storing the result in a scratch buffer, the topology mapping is computed inline in the vertex shader from `SV_VertexID`. No index buffer, no scratch buffer, no pre-draw dispatch.

### Quad List

An Xbox quad list has quads at vertex indices `[0,1,2,3], [4,5,6,7], ...` Each quad becomes two host triangles using vertices `[0,1,2]` and `[0,2,3]`. Draw `numQuads * 6` null vertices as a triangle list:

```hlsl
uint VSIndexFromQuadList(uint vertId) {
    uint quad  = vertId / 6u;
    uint local = vertId % 6u;
    uint base  = quad * 4u;
    static const uint lut[6] = { 0u, 1u, 2u, 0u, 2u, 3u };
    return base + lut[local];
}
```

### Triangle Fan

A fan with N vertices produces `N-2` triangles, all sharing vertex 0. Draw `(numVerts - 2) * 3` null vertices:

```hlsl
uint VSIndexFromTriFan(uint vertId) {
    uint tri   = vertId / 3u;
    uint local = vertId % 3u;
    // local 0 → always fan apex (vertex 0)
    // local 1 → tri + 1
    // local 2 → tri + 2
    return local == 0u ? 0u : (tri + local);
}
```

The apex vertex (index 0) is fetched repeatedly, but the GPU's post-transform cache handles this — vertex 0 is computed once and reused.

### Indexed Variants

If the Xbox draw uses an index buffer, an extra level of indirection reads the Xbox index from XboxRAM after computing the logical vertex position:

```hlsl
uint VSIndexFromIndexedQuadList(uint vertId, uint xboxIndexBase, bool is16bit) {
    uint quad  = vertId / 6u;
    uint local = vertId % 6u;
    static const uint lut[6] = { 0u, 1u, 2u, 0u, 2u, 3u };
    uint logicalIdx = quad * 4u + lut[local];
    if (is16bit) {
        uint word  = XboxRAM.Load((xboxIndexBase + logicalIdx * 2u) & ~3u);
        uint shift = (logicalIdx & 1u) * 16u;
        return (word >> shift) & 0xFFFFu;
    }
    return XboxRAM.Load(xboxIndexBase + logicalIdx * 4u);
}
```

### Vertex Shader Entry

`PRIM_TYPE` is a compile-time `#define`, so unused branches compile away. A separate VS permutation is compiled per topology:

```hlsl
float4 main(uint vertId : SV_VertexID) : SV_Position {
    uint xboxVtxIdx;
    switch (PRIM_TYPE) {
        case NV2A_QUADS:    xboxVtxIdx = VSIndexFromQuadList(vertId); break;
        case NV2A_TRI_FAN:  xboxVtxIdx = VSIndexFromTriFan(vertId);   break;
        default:            xboxVtxIdx = vertId;                       break;
    }
    return FetchVertex(xboxVtxIdx);
}
```

---

## Register Combiner Interpreter

### Background

The NV2A pixel pipeline is driven by a programmable register combiner — up to eight stages, each performing two products (AB and CD) and a sum or MUX, followed by an optional final combiner. The state is defined entirely by PGRAPH registers: input byte selectors, output destination registers, output mapping (scale/bias), per-stage flags, and texture mode descriptors.

The original Cxbx-Reloaded DX11 branch implements this as an SM3.0 HLSL interpreter (`RegisterCombinerInterpreter.fx`). The port described here upgrades it to DX11/SM5 with three structural optimisations and a full set of bug fixes.

### Bug Fixes

The SM3.0 original contains several correctness bugs:

**Texture register write** — a `max()` clamp on the texture register index caused all four texture stages to always write T3. Fixed by direct indexed write: `RegWrite(Regs, PS_REGISTER_T0 + stage, val)`.

**R0 initialisation** — R0.rgb was initialised to white `float4(1,1,1, T0.a)`. Per NV2A spec, R0.rgb starts at zero and R0.a is initialised from T0.a: `RegWriteA(Regs, PS_REGISTER_R0, Regs[PS_REGISTER_T0].a)`.

**PS\_CHANNEL\_RGB passthrough** — when the channel selector bit is clear (RGB path), the input value should pass through intact. The SM3.0 version used `.rgbb`, replacing alpha with blue. Fixed to use `val.aaaa` only when `useAlphaC` is true; otherwise the full `float4` passes through unchanged.

**CLIPPLANE** — the SM3.0 version called `Sample2D` after the discard test. Fixed to discard only; T\[stage\] stays at `(0,0,0,1)`.

**DOT\_RFLCT\_DIFF** — incorrectly called `reflect()`. The DIFF mode uses the constructed normal directly as the cubemap direction; the SPEC mode uses reflect. Fixed by removing the `reflect()` call from the DIFF path.

### Optimisation Step 2: Register File as float4 Array

The SM3.0 interpreter used `switch` chains with 12+ cases for every register read and write, because SM3.0 could not index registers indirectly. SM5 supports this natively:

```hlsl
float4 Regs[16];  // indexed by PS_REGISTER_* constants

float4 RegRead(float4 Regs[16], uint idx) {
    return Regs[idx & 0xFu];
}

void RegWrite(inout float4 Regs[16], uint idx, float4 val) {
    // Silently drop writes to ZERO, C0, C1, reserved slots
    if (idx == PS_REGISTER_ZERO || idx == PS_REGISTER_C0 ||
        idx == PS_REGISTER_C1   || idx == 6u || idx == 7u)
        return;
    Regs[idx & 0xFu] = val;
}
```

C0 and C1 are resolved directly from the cbuffer in `ResolveInput()` rather than being stored in the array — no slot needed. ZERO is never written, remaining permanently at zero from initialisation.

### Optimisation Step 3: Compile-Time Stage Count

The combiner loop runs 1–8 iterations, but the GPU cannot eliminate dead iterations without knowing the count at compile time. Eight shader variants are compiled with `#define NUM_STAGES N` for N in 1..8:

```hlsl
#ifndef NUM_STAGES
#define NUM_STAGES 8   // safe fallback
#endif

[unroll]
for (uint stage = 0u; stage < (uint)NUM_STAGES; stage++) {
    if (stage < numStages)
        DoCombinerStage(Regs, stage, flagMuxMsb, flagUniqueC0, flagUniqueC1);
}
```

For a variant compiled with `NUM_STAGES=2`, the compiler unrolls to exactly two iterations and eliminates the remaining six stages entirely. The runtime guard `stage < numStages` provides safety in the `NUM_STAGES=8` fallback. The host dispatches the correct variant based on `PSCombinerCount`.

### Optimisation Step 4: Merged RGB and Alpha Passes

The SM3.0 version called `do_color_combiner_stage()` twice per stage — once for RGB and once for alpha — decoding the output control bytes redundantly in each call. The merged version decodes all per-stage state once, batches all eight input reads together (enabling the compiler to schedule them as a unit), computes both RGB and alpha results before any write, then writes RGB outputs followed by alpha outputs:

```hlsl
void DoCombinerStage(inout float4 Regs[16], uint stage,
                     bool flagMuxMsb, bool flagUniqueC0, bool flagUniqueC1) {
    uint rgbIn  = PSRGBInputs[stage];    // read from XboxRAM via PGRAPH offset
    uint aIn    = PSAlphaInputs[stage];
    uint rgbOut = PSRGBOutputs[stage];
    uint aOut   = PSAlphaOutputs[stage];

    // Decode output control bits once for both RGB and alpha
    bool flagCDDot  = GetBit(rgbOut, 12u);
    bool flagABDot  = GetBit(rgbOut, 13u);
    bool flagRGBMux = GetBit(rgbOut, 14u);
    bool flagAMux   = GetBit(aOut,   14u);

    // Fetch all eight inputs in one block
    float4 rgbA = ResolveInput(Regs, GetByte(rgbIn, 3u), false, stage, ...);
    float4 rgbB = ResolveInput(Regs, GetByte(rgbIn, 2u), false, stage, ...);
    float4 rgbC = ResolveInput(Regs, GetByte(rgbIn, 1u), false, stage, ...);
    float4 rgbD = ResolveInput(Regs, GetByte(rgbIn, 0u), false, stage, ...);
    float aA = ResolveInput(Regs, GetByte(aIn, 3u), true, stage, ...).a;
    float aB = ResolveInput(Regs, GetByte(aIn, 2u), true, stage, ...).a;
    float aC = ResolveInput(Regs, GetByte(aIn, 1u), true, stage, ...).a;
    float aD = ResolveInput(Regs, GetByte(aIn, 0u), true, stage, ...).a;

    // Compute AB and CD (dot product for RGB, multiply for alpha)
    float3 rgbAB = flagABDot ? (float3)dot(rgbA.rgb, rgbB.rgb) : rgbA.rgb * rgbB.rgb;
    float3 rgbCD = flagCDDot ? (float3)dot(rgbC.rgb, rgbD.rgb) : rgbC.rgb * rgbD.rgb;
    float  aAB   = aA * aB;
    float  aCD   = aC * aD;

    // SUM or MUX; apply output mapping; write RGB then alpha
    // (RGB writes precede alpha writes so both paths can target the same
    //  destination register without conflict — alpha write patches only .a)
}
```

### Bit Extraction

All `fmod`/`floor`-based bit extraction from the SM3.0 version is replaced with native SM5 bitwise operators:

```hlsl
uint GetByte  (uint v, uint n) { return (v >> (n * 8u)) & 0xFFu; }
uint GetNibble(uint v, uint n) { return (v >> (n * 4u)) & 0x0Fu; }
bool GetBit   (uint v, uint n) { return ((v >> n) & 1u) != 0u; }
```

### Texture Handling in the Combiner

The combiner's texture stage fetch supports all NV2A texture modes: `PROJECT2D`, `PROJECT3D`, `CUBEMAP`, `PASSTHRU`, `CLIPPLANE`, `BUMPENVMAP`, `DPNDNT_AR`, `DPNDNT_GB`, `DOTPRODUCT`, `DOT_ST`, `DOT_ZW`, `DOT_RFLCT_DIFF`, `DOT_RFLCT_SPEC`, `DOT_STR_3D`, `DOT_STR_CUBE`. Each mode reads its input coordinates and prior T register values from the `Regs[]` array and issues the appropriate sample call.

DX11 error X4539 (sampler array elements cannot be used with different intrinsics) is avoided by declaring individual texture objects per stage rather than arrays:

```hlsl
Texture2D   Tex2D_0   : register(t0);  Texture2D   Tex2D_1   : register(t1);
Texture3D   Tex3D_0   : register(t4);  TextureCube TexCube_0 : register(t8);
SamplerState Samp0    : register(s0);  SamplerState Samp1    : register(s1);
```

---

## Texture Handling

Textures are the one component that cannot be handled entirely from XboxRAM without a pre-draw pass. The NV2A stores textures in Morton (Z-order) swizzled layout. DX11's `SamplerState` requires a `Texture2D` resource — it cannot filter arbitrary bytes from a `ByteAddressBuffer`.

### Deswizzle Compute Shader

A pre-draw CS deswizzles from XboxRAM into a pool of `Texture2D` resources. Each texture in the pool corresponds to a surface address registered in PGRAPH. The CS is dispatched only when the corresponding XboxRAM page range is marked dirty:

```hlsl
RWTexture2D<float4> HostTex : register(u0);
ByteAddressBuffer   XboxRAM : register(t0);

[numthreads(8, 8, 1)]
void Deswizzle(uint3 tid : SV_DispatchThreadID) {
    uint2 coord  = tid.xy;
    uint  morton = MortonEncode(coord.x, coord.y);
    uint  xboxOff = SURFACE_BASE + morton * BYTES_PER_TEXEL;
    uint4 texel   = XboxRAM.Load4(xboxOff);
    HostTex[coord] = DecodeTexelFormat(texel, TEXEL_FMT);
}
```

The deswizzle CS writes directly into a `Texture2D UAV`, making the result immediately sample-ready with no intermediate copy.

### Dirty Tracking

Two dirty bits per 4KB page:

| xboxDirty | hostDirty | State |
|-----------|-----------|-------|
| 0 | 0 | Coherent |
| 1 | 0 | CPU modified — needs deswizzle CS before next draw |
| 0 | 1 | GPU modified — needs readback if CPU reads |

The texture pool is keyed on `(surfaceAddress, format, width, height, mipCount)`. Pool entries are evicted LRU when the pool exceeds a size limit.

### Block-Compressed Textures

NV2A DXT textures are Morton-swizzled at the block level — each 4×4 DXT block is treated as a unit for the swizzle. The deswizzle CS encodes block coordinates into Morton order, reads 8-byte (DXT1) or 16-byte (DXT3/5) block payloads, and writes them linearly. The host `Texture2D` is declared with `DXGI_FORMAT_BC1_UNORM` or equivalent so the hardware decompresses during sampling.

---

## Pipeline State Management

DX11 pipeline state (blend equation, blend factors, depth compare, cull mode, alpha test, stencil) must be expressed as `ID3D11BlendState`, `ID3D11DepthStencilState`, and `ID3D11RasterizerState` objects, all created and cached on the CPU.

These are read from the PGRAPH shadow (not from XboxRAM on the GPU) because DX11 has no mechanism for GPU-driven pipeline state. A cache keyed on the packed state fields covers the common case — most Xbox titles use fewer than 20 distinct pipeline state combinations per frame.

This is the one component that cannot be made GPU-driven in DX11. It is addressed in the Vulkan architecture document.

---

## Per-Draw CPU Work Summary

With the complete architecture in place, the CPU's per-draw rendering work is:

| Task | CPU involvement |
|------|-----------------|
| Dirty XboxRAM pages | `UpdateSubresource` for dirty byte ranges |
| Texture deswizzle | Dispatch CS for dirty surface pages only |
| Vertex data | None — VS fetches from XboxRAM directly |
| Index/topology | None — VS computes from `SV_VertexID` |
| Combiner constants | None — PS reads from XboxRAM via PGRAPH offsets |
| Vertex shader constants | None — VS reads from XboxRAM at known offsets |
| Pipeline state | Read PGRAPH shadow, look up or create DX11 state objects |
| Draw call | `Draw(hostVertexCount, 0)` with null IA buffers |

---

## Optimisation Tiers for the Combiner Interpreter

Beyond the structural changes (Steps 2, 3, 4), additional optimisation is possible:

**CPU-side PGRAPH decode (Tier 1):** All constant fields in PGRAPH — texture mode flags, combiner count, final combiner settings — are uniform across all pixels in a draw call. They can be decoded on the CPU and placed in a small `cbuffer`, eliminating repeated per-pixel bit extraction on the hot path.

**Output mapping baked into cbuffer (Tier 1):** The 8-value output mapping is fully determined at draw time. Pre-computed scale and bias floats per stage output eliminate the shader-side decode entirely.

**Partial specialization cache (Tier 3):** A cache of 20–50 pre-compiled shader variants keyed on a hash of the combiner topology (stage count, dot/mux flags per stage, BlueToAlpha flags) covers 95% of real Xbox title combinations. Cache misses fall back to the full interpreter.

---

## Performance Analysis

### Bottleneck

The combiner interpreter, inline vertex decode, and topology arithmetic are all ALU-bound with enormous headroom against Xbox-era geometry and fill rates. The bottleneck is **texture sampling via the deswizzle path at high internal resolutions**.

The texture pool `Texture2D` approach with `SamplerState` uses hardware texture units with 2D-tiled L1 cache. This is the correct approach for DX11 and avoids the bandwidth penalty that would arise from manual `ByteAddressBuffer` texture sampling.

### Expected Performance

| Target resolution | Bottleneck | Headroom vs Xbox GPU |
|------------------|------------|----------------------|
| 640×480 native | Nothing — GPU massively idle | >1000× |
| 2× IR (1280×960) | Nothing | >200× |
| 4× IR (2560×1920) | Deswizzle CS frequency | ~50× |
| 4K (3840×2160) | Texture bandwidth | 5–20× (GPU dependent) |

### ALU Budget

At 4K 60fps with a 2-stage combiner (typical Xbox median):

- Combiner interpreter: ~150 ops/pixel → 74.7 Gops/sec → **0.6% of RTX 3060 ALU**
- Morton + format decode: ~60 ops/pixel → **0.2% of ALU**
- VS topology + CMP decode: ~20 ops/vertex at 3M verts/sec → **<0.01% of ALU**

All ALU-bound components are negligible. The deswizzle texture pool approach provides hardware-filtered sampling at the cost of one CS pass per dirty texture per draw. On modern hardware this is well within budget for any Xbox title.

---

## DX11 Binding Summary

```
t0   ByteAddressBuffer XboxRAM     — all Xbox memory, read by VS and PS
t1   Buffer<float4> VtxStream_*    — typed SRV views over XboxRAM vertex regions
t2+  Texture2D/3D/Cube HostTex[N]  — deswizzled texture pool

u0   RWTexture2D HostTex           — deswizzle CS write target (pre-draw only)

b0   cbuffer DecodedCombinerState  — optional Tier 1 optimisation

s0..s3  SamplerState               — one per active texture stage
```

The VS and PS both read primarily from `t0` (`XboxRAM`). During the pre-draw deswizzle CS pass, `XboxRAM` is bound as an SRV and the target texture as a UAV. Before the draw, the UAV binding is released and the texture is rebound as an SRV.

---

## Upscaling

### Internal Resolution

The NV2A renders to a framebuffer surface whose address and dimensions are registered in PGRAPH. The host render target does not need to match those dimensions. Rendering at a fixed integer multiple (the internal resolution multiplier, IRm) enlarges the render target while leaving all other pipeline components unchanged:

```
Xbox framebuffer:    640 × 480   (from PGRAPH NV_PGRAPH_SURFACE)
Host render target:  2560 × 1920 (IRm = 4)
Texture pool:        native Xbox dimensions  (unchanged)
Deswizzle CS:        native Xbox dimensions  (unchanged)
```

The VS, topology arithmetic, vertex fetch, and combiner interpreter run once per host pixel — more work overall, but identical per-pixel cost to native resolution.

### Surface Type Distinction

The most important architectural consequence of IR rendering is that the texture pool must distinguish between two categories of Xbox surface.

**Static textures** are written by the Xbox CPU and never rendered to by the GPU. These are deswizzled at native Xbox dimensions. Sampling them at any IRm is correct because UV coordinates are normalised — the larger render target samples the same normalised region.

**Render target surfaces** are written by the GPU (framebuffer, shadow maps, reflection captures, post-process intermediates). These must be maintained at IR dimensions so that the geometry rasterised into them matches the host render resolution.

The distinction is made by tracking which Xbox surface addresses appear as NV097 render target registers (`NV097_SET_SURFACE_COLOR_OFFSET`, `NV097_SET_SURFACE_ZETA_OFFSET`). Any address that has been a render target is maintained at IR dimensions; all others are deswizzled at native size.

```cpp
struct SurfaceEntry {
    uint32_t xboxAddr;
    uint32_t nativeW, nativeH;   // Xbox dimensions
    bool     isRenderTarget;     // true when seen as NV097 RT
    ID3D11Texture2D* hostTex;    // native-size if !isRenderTarget
                                 // IR-size    if  isRenderTarget
};
```

### Pixel-Offset Screen-Space Effects at IR

Normalised UV coordinates survive IR scaling correctly. Fixed-pixel-offset effects do not. Xbox games that sample a render target with a hardcoded texel offset (blur kernels, SSAO, edge detection, motion blur accumulation) compute those offsets in terms of native pixels:

```hlsl
// Reconstructed Xbox game shader — native-resolution assumptions
float2 texelSize = float2(1.0 / 640.0, 1.0 / 480.0);
float4 blurred   = tex2D(frameTex, uv + texelSize * 3.0);  // 3-pixel offset
```

At IRm=4, this 3-pixel offset samples 3/2560 of the host surface — one quarter of the intended extent. The effect is weaker than designed.

There is no general solution without game-specific knowledge. Temporal upscalers (FSR 2, XeSS) partially compensate by reconstructing detail across frames, which can mask weakened screen-space effects. In practice most blur-based effects remain visually acceptable at reduced strength.

### LOD Bias at IR

At IRm=4, the GPU's automatic mip level selection (derived from `ddx`/`ddy` of screen-space UV) selects sharper mips than at native resolution — generally desirable, as textures appear more detailed. However it partially counteracts a game-applied positive LOD bias intended to soften a texture.

A global LOD bias offset applied to all samplers controls this behaviour:

```cpp
// 0.0  → IR selects sharper mips naturally (recommended for visual quality)
// +log2(IRm) → restore native mip selection exactly
float globalLodBiasOffset = 0.0f;
```

In DX11, LOD bias is baked into `ID3D11SamplerState`. Sampler objects are created per distinct LOD bias value and cached; the per-stage PGRAPH bias is combined with the global offset at sampler creation time.

### Point Sprite and Line Width Scaling

**Point sprites:** `NV_PGRAPH_POINTSIZE` is expressed in native pixels. The point sprite GS must scale by IRm:

```hlsl
float pointSize = ReadPGRAPHPointSize() * IR_SCALE;  // IR_SCALE = float(IRm) from cbuffer
// ... expand quad by pointSize in screen space
```

**Thick lines:** The thick line GS expands lines to screen-aligned quads by a pixel width. The same `IR_SCALE` multiplier applies to the expansion offset.

### 2D / HUD Pass Detection

Xbox games render the 3D scene first, then composite a 2D HUD or menu using pre-transformed vertices (W-divided, passthrough VS — the NV2A equivalent of `D3DFVF_XYZRHW`). These passes contain UI elements designed for exact native-pixel placement.

**Detection:** The NV2A passthrough VS (`CxbxVertexShaderPassthrough.hlsl`) is active when PGRAPH's programmable VS enable flag (`NV_PGRAPH_CSV0_D`) is clear.

**Option A — Scale in the passthrough VS:** Multiply screen-space XY by `IR_SCALE`. The 2D content renders at IR dimensions. Text and pixel-art UI may appear blurry because the source pixels were designed for native size.

**Option B — Separate native-resolution target:** Render all 2D passes to a native-resolution `ID3D11Texture2D`. After the 3D IR upscale, composite the native-resolution 2D layer onto the display output using nearest-neighbour scaling. This preserves pixel-perfect UI at the cost of one additional composite pass. Correct for games with sharp pixel-art HUDs; optional for games with resolution-independent UI.

### Motion Vector Reconstruction

Temporal upscalers (FSR 2, XeSS) require a screen-space motion vector buffer — a per-pixel 2D velocity indicating where each pixel was in the previous frame. The NV2A produces none; they must be reconstructed in the VS by reprojecting world-space positions through the previous frame's VP matrix.

```hlsl
cbuffer UpscaleData : register(b1) {
    float4x4 prevViewProj;  // VP from previous frame, updated by CPU each frame
    float2   jitter;        // subpixel jitter for this frame
};

struct VS_OUTPUT {
    float4 pos     : SV_Position;
    float4 prevPos : TEXCOORD8;
    // ...
};

VS_OUTPUT main(uint vertId : SV_VertexID) {
    float4 worldPos  = FetchVertex(vertId);
    float4 clip      = mul(currentViewProj, worldPos);
    float4 prevClip  = mul(prevViewProj,    worldPos);
    VS_OUTPUT o;
    o.pos     = clip;
    o.prevPos = prevClip;
    return o;
}
```

A motion vector PS writes the velocity buffer:

```hlsl
float2 main(VS_OUTPUT i) : SV_Target {
    float2 curr = (i.pos.xy     / i.pos.w)     - jitter;  // de-jitter current
    float2 prev = (i.prevPos.xy / i.prevPos.w);
    return (curr - prev) * float2(0.5, -0.5);             // clip → UV-space
}
```

`prevViewProj` is maintained by the CPU: after each frame, the VS constant block is copied from the PGRAPH shadow. Identifying which slot contains the VP matrix uses heuristics — the 4×4 block with the largest column magnitudes and a perspective-characteristic bottom row. This is confirmable per title.

**Limitation:** Motion vectors are incorrect for vertex-shader-animated geometry (skinned characters). FSR 2 and XeSS tolerate this with reduced temporal stability on animated objects rather than hard artefacts.

### Jitter Injection

Temporal upscalers require the projection matrix to carry a sub-pixel jitter offset (Halton sequence) each frame. In the DX11 architecture, the CPU intercepts all VS constant writes and injects jitter into the identified projection matrix before it reaches XboxRAM:

```cpp
void OnVSConstantWrite(uint32_t slot, const float* data, uint32_t floats) {
    if (IsProjectionMatrix(slot, data, floats)) {
        float2 jitter   = HaltonJitter(frameIndex % 8);
        float4x4 jittered = *(float4x4*)data;
        // Add jitter to the clip-space XY translation row
        jittered[2][0] += jitter.x * 2.0f / irWidth;
        jittered[2][1] += jitter.y * 2.0f / irHeight;
        WriteToXboxRAM(VS_CONSTANTS_BASE + slot * 16, &jittered, floats * 4);
    } else {
        WriteToXboxRAM(VS_CONSTANTS_BASE + slot * 16, data, floats * 4);
    }
}
```

The same jitter value is passed to the motion vector PS as `UpscaleData.jitter` so the velocity computation can remove it, ensuring the upscaler sees only true inter-frame motion.

**Projection matrix identification** uses a combination of: the 4×4 float block that contains a value near 1.0 in position [3][2] (perspective divide row), non-zero off-diagonal elements (view rotation), and reasonable near/far plane ratios. A per-title override table handles edge cases.

### Framebuffer Readback at IR

When the Xbox CPU reads from a surface the GPU has rendered to at IR dimensions, a downsample precedes the readback. The IR render target is bilinearly downsampled to native Xbox dimensions before being written into the shadow:

```cpp
void OnXboxCPUReadSurface(uint32_t addr) {
    SurfaceEntry* e = FindSurface(addr);
    if (e && e->isRenderTarget) {
        // Downsample IR → native via bilinear blit PS
        BlitToNativeStaging(e->hostTex, e->nativeStagingTex, e->nativeW, e->nativeH);
        D3D11_MAPPED_SUBRESOURCE mapped;
        ctx->Map(e->nativeStagingTex, 0, D3D11_MAP_READ, 0, &mapped);
        CopyLinearToXboxShadow(shadow + addr, mapped.pData, e->nativeW, e->nativeH);
        ctx->Unmap(e->nativeStagingTex, 0);
    }
}
```

The Xbox CPU always sees native-resolution data.

### Available Upscalers in DX11

| Upscaler | Type | Motion Vectors | DX11 Support | Notes |
|---|---|---|---|---|
| Bilinear blit | Spatial | No | Native | Zero cost; baseline |
| Integer scale | Spatial | No | Custom PS | Nearest-neighbour; pixel-perfect retro output |
| AMD FSR 1 | Spatial | No | HLSL source | Ships as HLSL; PS pass over IR target |
| Nvidia NIS | Spatial | No | HLSL source | Sharpening + upscale |
| AMD FSR 2 | Temporal | **Required** | Official DX11 SDK | Best quality; tolerates imperfect MVs |
| Intel XeSS | Temporal | **Required** | DX11 DP4a path | Works on all SM6-capable hardware |
| AMD FSR 3 (frame gen) | Temporal + frame gen | **Required** | Frame gen DX12 only | Upscale works; frame gen not available in DX11 |
| Nvidia DLSS | Temporal | **Required** | **Not available** | Requires Vulkan or DX12 |

FSR 2 with reconstructed motion vectors is the recommended choice. It produces better results than FSR 1 for the textured 3D scenes typical of Xbox-era games and tolerates imperfect motion vectors on animated geometry better than XeSS.

---

## Texture Replacement

### Concept

A texture replacement pack maps Xbox texture hashes to high-resolution host assets. When a texture is about to be deswizzled and inserted into the pool, its hash is checked against the pack database. On a match, the pack asset is loaded instead. The game renders at IR with high-resolution textures without any modification to the Xbox executable.

### Hashing

The hash is computed over the raw swizzled bytes in XboxRAM at the texture's registered surface address and byte extent. Hashing raw Xbox data (before deswizzle or format conversion) produces a stable identifier independent of host GPU, driver, or IR multiplier.

```cpp
uint64_t HashTexture(uint32_t xboxAddr, uint32_t byteSize) {
    return XXH64(shadow + xboxAddr, byteSize, /* seed */ 0);
}
```

The hash is recomputed only when the texture's dirty page bit is set — once per level load for static art, never again during gameplay. The cost is negligible.

### Pool Integration

```
1. Compute xxHash64 over raw XboxRAM bytes at texture address
2. Look up hash in replacement pack database
3a. Hit  → load replacement asset (async); use native deswizzled texture meanwhile
3b. Miss → deswizzle CS → insert native texture into pool
4. Bind result for draw
```

Replacement assets load on a background thread. The native deswizzled texture is used on first encounter with no frame hitch; the pool entry is atomically swapped when the replacement is ready.

### Render Target Surfaces

Render target surfaces are generated by GPU rendering each frame and cannot be replaced. The pool's `isRenderTarget` flag gates the hash lookup — render target surfaces skip the replacement check entirely.

### Format and Channel Correspondence

Replacement assets must match the channel semantics of the Xbox format, or a conversion pass is applied on load. Most packs ship RGBA PNG or BC7 DDS. The loader verifies channel correspondence with the Xbox PGRAPH format and applies a swizzle where needed (e.g. BGRA Xbox → RGBA replacement).

Normal map replacements in OpenGL convention (Y-up, full XYZ in RGB) require reconstruction from the combiner interpreter's expected encoding. This is handled per-format by detecting an active normal map replacement and routing it through a corrected decode path in the PS.
