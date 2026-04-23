# NV2A GPU Emulation — Vulkan Architecture

## Overview

This document describes a GPU-driven architecture for emulating the Xbox NV2A graphics processor using Vulkan. It extends the DX11 architecture (described in the companion document) by resolving every limitation that DX11's abstraction model imposed, converging on a design where the CPU is entirely absent from the rendering data path.

The three DX11 limitations that Vulkan directly addresses are:

1. **Texture deswizzle pre-pass** — DX11 requires a `Texture2D` resource with hardware tiling; the NV2A stores textures in Morton order. The deswizzle CS was the one remaining pre-draw GPU work item.
2. **PSO creation and caching** — DX11 encodes blend, depth, cull, and stencil state in immutable pipeline state objects. The CPU had to read PGRAPH shadow state and manage a PSO cache per draw.
3. **Push buffer parsing** — DX11's `ExecuteIndirect` is too restricted to replay NV097 command streams from a compute shader.

Vulkan with extensions available on all current desktop drivers eliminates all three.

---

## Memory Model

### Single VkDeviceMemory Allocation

The entire 64MB Xbox address space is represented as one `VkDeviceMemory` allocation. Multiple `VkResource` handles are bound to overlapping ranges of this allocation simultaneously — something Vulkan explicitly supports and DX11 does not.

```cpp
// One allocation backing everything
VkMemoryAllocateInfo allocInfo = {
    .allocationSize  = 64 * 1024 * 1024,
    .memoryTypeIndex = findMemoryType(
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT |
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
        VK_MEMORY_PROPERTY_HOST_COHERENT_BIT) // resizable BAR
};
vkAllocateMemory(device, &allocInfo, nullptr, &xboxMemory);

// Buffer view — vertex fetch, PGRAPH reads, push buffer
VkBufferCreateInfo bufInfo = { .size = 64 * 1024 * 1024,
    .usage = VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
             VK_BUFFER_USAGE_UNIFORM_TEXEL_BUFFER_BIT |
             VK_BUFFER_USAGE_SHADER_DEVICE_ADDRESS_BIT };
vkCreateBuffer(device, &bufInfo, nullptr, &xboxRAM_buf);
vkBindBufferMemory(device, xboxRAM_buf, xboxMemory, 0);
```

The same `xboxMemory` backing is also used for texture `VkImage` resources at identity-mapped offsets (see Texture Handling).

### Persistent Mapped Memory (Resizable BAR)

On resizable BAR hardware — all current discrete GPUs support it via PCIe BAR resize — `VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | DEVICE_LOCAL_BIT` memory can be persistently mapped. The CPU writes directly into GPU memory at cache-line granularity with no driver involvement:

```cpp
void* mapped;
vkMapMemory(device, xboxMemory, 0, 64 * 1024 * 1024, 0, &mapped);
// Keep persistently mapped — no vkUnmapMemory needed

// Xbox CPU emulation writes go here directly
void OnXboxMemoryWrite(uint32_t addr, uint32_t value) {
    ((uint32_t*)mapped)[addr / 4] = value;
    MarkPageDirty(addr >> 12);  // 4KB page granularity
}
```

`UpdateSubresource` (which in DX11 always copies through the driver) becomes a direct `memcpy` into mapped GPU memory. On BAR hardware the cost approaches a CPU L3 cache write. The dirty-range flush that was DX11's per-draw overhead is effectively free.

On non-BAR hardware, a staging buffer is used with `vkCmdCopyBuffer` for dirty ranges — the same cost as DX11's `UpdateSubresource`, but the BAR path covers all modern discrete GPUs.

### Buffer Device Address

`VK_KHR_buffer_device_address` (core in Vulkan 1.2) provides a `uint64_t` GPU virtual address for any `VkBuffer`. Compute shaders can dereference this address directly using pointer-style arithmetic in GLSL:

```glsl
layout(buffer_reference, std430) buffer XboxRAMRef {
    uint data[];
};

// In any shader:
uint64_t xboxRAMAddr = pushConstants.xboxRAMBase;
XboxRAMRef ram = XboxRAMRef(xboxRAMAddr);
uint val = ram.data[offset / 4];
```

This enables a compute shader parsing the push buffer to dereference surface addresses and vertex buffer pointers found in NV097 method arguments directly into XboxRAM without any translation table — the NV2A's flat address model maps one-to-one onto Vulkan's device address space.

---

## Texture Handling

### The DX11 Problem

DX11's `SamplerState` requires a `Texture2D` resource with hardware tiling. The NV2A stores textures in Morton (Z-order) swizzled layout. This forced a pre-draw deswizzle CS writing to a `Texture2D` pool — the one remaining pre-draw pass that prevented a fully push-buffer-driven pipeline.

### Solution: VkImage Over Shared Memory

Vulkan allows binding a `VkImage` and a `VkBuffer` to the same `VkDeviceMemory` region at the same offset. Combined with `VK_EXT_image_drm_format_modifier`, a `VkImage` can be created with a custom memory layout that describes the Morton swizzle pattern. The hardware sampler then traverses the NV2A-format memory directly, with no copy and no deswizzle pass:

```cpp
// Describe the Morton tiling as a DRM format modifier
VkDrmFormatModifierPropertiesEXT modifier = {
    // Morton/Z-order descriptor matching NV2A swizzle for this texture size
};

VkImageDrmFormatModifierExplicitCreateInfoEXT modifierInfo = {
    .sType = VK_STRUCTURE_TYPE_IMAGE_DRM_FORMAT_MODIFIER_EXPLICIT_CREATE_INFO_EXT,
    .drmFormatModifier = NV2A_MORTON_MODIFIER,
    .drmFormatModifierPlaneCount = 1,
    .pPlaneLayouts = &subresourceLayout
};

VkImageCreateInfo imageInfo = {
    .tiling = VK_IMAGE_TILING_DRM_FORMAT_MODIFIER_EXT,
    .pNext  = &modifierInfo,
    // width, height, format from PGRAPH texture descriptor
};
vkCreateImage(device, &imageInfo, nullptr, &texImage);

// Bind to the same VkDeviceMemory as XboxRAM, at the texture's Xbox address
vkBindImageMemory(device, texImage, xboxMemory, textureXboxAddress);
```

The `VkImage` shares physical memory with `XboxRAM`. The CPU's Xbox memory write goes into the persistently-mapped buffer. The `VkImage` view at the same offset reads the same bytes — the hardware sampler decodes the Morton layout described by the modifier. **No deswizzle CS, no texture pool, no copy of any kind.**

### Fallback: Zero-Copy VkImage Write

On drivers that do not support the specific Morton DRM format modifier, the fallback uses a small CS that writes directly into `VkImage` memory (rather than copying to a separate buffer first as in DX11). The result is immediately sample-ready:

```glsl
// Deswizzle CS writes to VkImage UAV directly
layout(set=0, binding=0, rgba8) uniform image2D hostTex;
layout(set=0, binding=1) buffer XboxRAM { uint data[]; };

layout(local_size_x=8, local_size_y=8) in;
void main() {
    uvec2 coord  = gl_GlobalInvocationID.xy;
    uint  morton = MortonEncode(coord.x, coord.y);
    uint  raw    = data[SURFACE_BASE/4 + morton];
    imageStore(hostTex, ivec2(coord), DecodeTexel(raw));
}
```

Even this fallback CS is triggered only on dirty pages. In the common case (texture not modified since last draw), the VkImage view is valid and the sampler reads NV2A memory directly.

### Mip Levels and Cube Maps

Each mip level and cube face has its own base address in PGRAPH. For the modifier path, each level is a separate `VkImage` plane or a separate `VkImage` with its own binding into `xboxMemory` at the level's Xbox base address. The sampler's LOD selection functions normally — no change to the shader.

---

## Pipeline State: Extended Dynamic State

### The DX11 Problem

DX11 encodes blend equation, depth compare, cull mode, stencil op, and fill mode into immutable `ID3D11BlendState`, `ID3D11DepthStencilState`, and `ID3D11RasterizerState` objects. Every unique combination of NV2A pipeline state required the CPU to look up or create a matching set of state objects. While cacheable, this was still a per-draw CPU read of the PGRAPH shadow and a hash table lookup.

### Solution: VK_EXT_extended_dynamic_state3

Available on all current desktop Vulkan drivers (core-promoted in Vulkan 1.3), this extension exposes a `vkCmd*` call for nearly every NV2A pipeline state register:

| NV2A PGRAPH Register | Vulkan Dynamic Command |
|---|---|
| `NV_PGRAPH_BLEND` (equation, factors) | `vkCmdSetColorBlendEquationEXT` |
| `NV_PGRAPH_CONTROL_0` (alpha test) | `vkCmdSetAlphaToOneEnableEXT` |
| `NV_PGRAPH_SETUPRASTER` (cull mode) | `vkCmdSetCullModeEXT` |
| `NV_PGRAPH_ZCOMPRESSOCCLUDE` (depth op) | `vkCmdSetDepthCompareOpEXT` |
| `NV_PGRAPH_CONTROL_2` (stencil op) | `vkCmdSetStencilOpEXT` |
| `NV_PGRAPH_CONTROL_0` (fill mode) | `vkCmdSetPolygonModeEXT` |
| `NV_PGRAPH_CONTROL_0` (depth bias) | `vkCmdSetDepthBiasEnableEXT` + `vkCmdSetDepthBias` |
| `NV_PGRAPH_BLEND` (write mask) | `vkCmdSetColorWriteMaskEXT` |
| `NV_PGRAPH_BLEND` (logic op) | `vkCmdSetLogicOpEXT` |

The CPU reads ~10 PGRAPH shadow fields and issues ~10 `vkCmd*` calls. Zero PSO creation, zero state object cache, zero stalls on first-seen state combinations. A pipeline is created once per render pass configuration (colour format, depth format, sample count) — a handful of variants total — and reused forever.

```cpp
void SetNV2APipelineState(VkCommandBuffer cmd, const uint32_t* pgraphShadow) {
    // Blend
    VkColorBlendEquationEXT blendEq = {
        .colorBlendOp = NV2ABlendOpToVk(pgraphShadow[NV_PGRAPH_BLEND/4]),
        // ...
    };
    vkCmdSetColorBlendEquationEXT(cmd, 0, 1, &blendEq);

    // Depth
    vkCmdSetDepthCompareOpEXT(cmd,
        NV2ADepthFuncToVk(pgraphShadow[NV_PGRAPH_ZCOMPRESSOCCLUDE/4]));

    // Cull
    vkCmdSetCullModeEXT(cmd,
        NV2ACullModeToVk(pgraphShadow[NV_PGRAPH_SETUPRASTER/4]));

    // Stencil, polygon mode, depth bias...
}
```

---

## Push Buffer Replay on GPU

### Background

The NV2A command stream is a push buffer — a sequence of NV097 `(method, value)` pairs processed by PFIFO and dispatched to PGRAPH. In the DX11 architecture, this is walked by the CPU: methods are decoded, PGRAPH state is updated in shadow memory, and draw calls are issued via DX11. This is a fundamental serialisation point — the CPU must finish parsing a push buffer packet before the GPU can act on it.

### VK_EXT_device_generated_commands

This extension (available on Nvidia, AMD, and Intel on current drivers) allows a compute shader to write a sequence of Vulkan command tokens into a `VkIndirectCommandsLayoutEXT`, which the GPU then executes directly. The token types include state changes, descriptor updates, push constant writes, and draw calls.

A CS parses the Xbox push buffer from XboxRAM and emits one token per NV097 method:

```glsl
// Push buffer parser CS — one thread per push buffer DWORD
layout(local_size_x = 64) in;
void main() {
    uint word = xboxRAM.data[PUSH_BUFFER_BASE/4 + gl_GlobalInvocationID.x];

    if (IsDataPacket(word)) {
        uint method = (word >> 2u) & 0x7FFu;
        uint value  = /* next DWORD */ xboxRAM.data[...];

        switch (method) {
        case NV097_SET_BLEND_EQUATION:
            EmitToken(TOKEN_SET_COLOR_BLEND_EQUATION, value);
            break;
        case NV097_DRAW_ARRAYS:
            uint count = (value >> 24u) & 0xFFu;
            uint start = value & 0xFFFFFFu;
            EmitDrawToken(ComputeHostVertexCount(count), start);
            break;
        // ... all NV097 methods
        }
    }
}

// GPU then executes:
// vkCmdExecuteGeneratedCommandsEXT(cmd, indirectCommandsBuffer)
```

The CPU's sequence for a complete frame:

1. Persistently-mapped write of Xbox CPU output into `xboxMemory` (dirty ranges, essentially free on BAR hardware)
2. Dispatch push buffer parser CS
3. `vkCmdExecuteGeneratedCommandsEXT`

The entire frame's state changes, texture binds, and draw calls are determined and executed by the GPU. The CPU issues three Vulkan calls per frame.

### Texture Descriptor Updates Within the Push Buffer

NV097 texture state changes (`NV097_SET_TEXTURE_OFFSET`, `NV097_SET_TEXTURE_FORMAT`) normally require updating descriptor sets — a CPU operation in standard Vulkan. With `VK_EXT_descriptor_buffer` (or `VK_EXT_mutable_descriptor_type`), descriptor data lives in a GPU buffer that the CS can write directly:

```glsl
// Parser CS writes texture descriptor into descriptor buffer
layout(buffer_reference) buffer DescriptorBuffer { uint64_t descriptors[]; };

case NV097_SET_TEXTURE_OFFSET:
    uint64_t texAddr = xboxRAMBase + value; // Xbox texture address → VkImage device address
    descriptorBuffer.descriptors[TEXTURE_SLOT(stage)] = texAddr;
    break;
```

The `VkImage` at that address (bound to `xboxMemory` at the texture's Xbox offset) is referenced by device address — no descriptor set update on the CPU.

---

## Render Pass Architecture: Subpass Dependencies

Xbox games frequently render to a surface and then read it as a texture in the same frame. In DX11, these transitions require explicit `OMSetRenderTargets` / `PSSetShaderResources` swaps and driver-inserted barriers.

Vulkan render passes with explicit subpass dependencies declare these transitions upfront, enabling the driver and hardware to schedule optimally:

```cpp
VkSubpassDependency dep = {
    .srcSubpass      = 0,
    .dstSubpass      = 1,
    .srcStageMask    = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT,
    .dstStageMask    = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT,
    .srcAccessMask   = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT,
    .dstAccessMask   = VK_ACCESS_INPUT_ATTACHMENT_READ_BIT,
    .dependencyFlags = VK_DEPENDENCY_BY_REGION_BIT   // tile-based optimisation
};
```

`VK_DEPENDENCY_BY_REGION_BIT` instructs a tile-based GPU (and Vulkan implementations on tile-based architectures) that the render target data does not need to leave the tile cache between subpasses. For render-to-texture effects, the resolved framebuffer is read back as a texture input without a DRAM round-trip. Xbox games that use this pattern get the optimisation for free on any tile-based Vulkan driver.

---

## Vertex Fetch and Topology Conversion

These components are identical in principle to the DX11 architecture. In GLSL/SPIR-V:

### CMP Inline Decode

```glsl
vec3 DecodeCMP(uint offset) {
    uint raw = xboxRAM.data[offset / 4];
    int x = int(raw << 21u) >> 21;   // sign-extend 11 bits
    int y = int(raw << 10u) >> 21;
    int z = int(raw       ) >> 22;   // sign-extend 10 bits
    return vec3(x / 1023.0, y / 1023.0, z / 511.0);
}
```

### Topology Arithmetic

```glsl
// Quad list
uint VSIndexFromQuadList(uint vertId) {
    uint quad  = vertId / 6u;
    uint local = vertId % 6u;
    const uint lut[6] = uint[6](0u, 1u, 2u, 0u, 2u, 3u);
    return quad * 4u + lut[local];
}

// Triangle fan
uint VSIndexFromTriFan(uint vertId) {
    uint tri   = vertId / 3u;
    uint local = vertId % 3u;
    return local == 0u ? 0u : (tri + local);
}
```

These are the same algorithms as in the DX11 document — the GPU-side implementation is API-agnostic.

---

## Register Combiner Interpreter

The combiner interpreter ports directly from HLSL to GLSL with minor syntax changes. All structural optimisations (Steps 2, 3, 4) and bug fixes carry over:

- **Step 2** — `vec4 Regs[16]` replaces switch chains. GLSL supports indirect indexing in SM5-equivalent compute and fragment stages.
- **Step 3** — `NUM_STAGES` specialisation constant via `layout(constant_id=0) const int NUM_STAGES = 8` — Vulkan's equivalent of `#define NUM_STAGES N`, evaluated at pipeline creation time rather than compile time.
- **Step 4** — Merged RGB and alpha combiner pass, same write ordering.

```glsl
layout(constant_id = 0) const int NUM_STAGES = 8;

void main() {
    // ...
    for (int stage = 0; stage < NUM_STAGES; stage++) {
        if (stage < numStages)
            DoCombinerStage(Regs, stage, flagMuxMsb, flagUniqueC0, flagUniqueC1);
    }
    // ...
}
```

Specialisation constants are set at `vkCreateGraphicsPipelines` time. The push buffer parser CS selects the correct `NUM_STAGES` value from `PSCombinerCount` and writes it into the pipeline's specialisation constant block. The pipeline cache handles deduplication.

---

## GPU-Driven Draws: vkCmdDrawIndirect

In DX11, the CPU computed the host vertex count (adjusting for topology conversion) and called `Draw(hostVertexCount, 0)`. In Vulkan, a CS writes `VkDrawIndirectCommand` structs into a `VkBuffer`, and `vkCmdDrawIndirectCount` executes them:

```cpp
// CS writes this structure; CPU never sees vertex counts
struct VkDrawIndirectCommand {
    uint32_t vertexCount;    // computed from NV2A vertex count + topology
    uint32_t instanceCount;  // always 1
    uint32_t firstVertex;    // from PGRAPH vertex array base
    uint32_t firstInstance;  // always 0
};

// One call covers all draws generated by the push buffer parser
vkCmdDrawIndirectCount(
    cmd,
    drawParamsBuffer, 0,         // buffer of VkDrawIndirectCommand
    drawCountBuffer,  0,         // draw count written by parser CS
    MAX_DRAWS_PER_FRAME,
    sizeof(VkDrawIndirectCommand)
);
```

---

## Full Frame Timeline

```
Frame N:

[CPU]  Xbox CPU emulation runs
       → direct writes into persistently-mapped xboxMemory
       → dirty page bitmap updated

[CPU]  vkBeginCommandBuffer

[CMD]  vkCmdPipelineBarrier  (host write → shader read)
[CMD]  vkCmdDispatch  PushBufferParserCS
         reads NV097 methods from XboxRAM
         writes VkDrawIndirectCommand[] into DrawParamsBuffer
         writes dynamic state tokens into IndirectCommandsBuffer
         writes texture descriptors into DescriptorBuffer (if VK_EXT_descriptor_buffer)

[CMD]  vkCmdBeginRenderPass
[CMD]  vkCmdExecuteGeneratedCommandsEXT  (IndirectCommandsBuffer)
         internally issues:
           vkCmdSetColorBlendEquationEXT  × N
           vkCmdSetDepthCompareOpEXT      × N
           vkCmdSetCullModeEXT            × N
           ... other dynamic state ...
           vkCmdDrawIndirectCount         (all draws from DrawParamsBuffer)
[CMD]  vkCmdEndRenderPass

[CPU]  vkEndCommandBuffer
[CPU]  vkQueueSubmit
```

From `vkBeginCommandBuffer` to `vkQueueSubmit` the CPU issues four Vulkan commands plus one dispatch and one `ExecuteGeneratedCommands`. Every state change, descriptor update, and draw within the frame is driven by the GPU.

---

## Per-Draw CPU Work Summary

| Task | DX11 | Vulkan |
|------|-------|--------|
| Dirty page upload | `UpdateSubresource` (copy) | `memcpy` into mapped BAR memory (near-free) |
| Texture deswizzle | CS dispatch per dirty texture | None — VkImage reads NV2A memory directly (modifier path) |
| Vertex data | None | None |
| Index / topology | None | None |
| Combiner constants | None | None |
| Pipeline state | CPU reads PGRAPH shadow, creates/looks up state objects | GPU emits dynamic state tokens from push buffer |
| Draw calls | CPU computes vertex count, calls `Draw()` | GPU writes and executes `VkDrawIndirectCommand` |

The CPU's rendering contribution in the Vulkan architecture is: dirty page `memcpy`, `vkBeginCommandBuffer`, one CS dispatch, `vkCmdExecuteGeneratedCommandsEXT`, `vkEndCommandBuffer`, `vkQueueSubmit`.

---

## Performance Analysis

### Texture Bottleneck: Eliminated

In the DX11 architecture, the hardware texture sampler required a deswizzle CS to populate a `Texture2D` pool. In Vulkan with the DRM format modifier, the sampler reads Xbox-format memory directly. There is no pre-draw copy, no extra DRAM bandwidth for the copy, and no dirty-tracking delay between a texture write and the next sample.

The performance impact of the modifier path depends on driver implementation. If the driver internally reorders reads to match hardware tile cache geometry, sampling performance is identical to a native `VK_IMAGE_TILING_OPTIMAL` texture. If the driver samples linearly with Morton address computation in the texture unit, the cache efficiency is somewhat lower — but this remains hardware-accelerated filtering, not software sampling.

In either case, the elimination of the deswizzle CS removes the primary bandwidth bottleneck identified in the DX11 analysis.

### Expected Performance

| Target resolution | Bottleneck | Headroom vs Xbox GPU |
|---|---|---|
| 640×480 native | Nothing | >1000× |
| 2× IR | Nothing | >250× |
| 4× IR | Nothing | >100× |
| 4K (3840×2160) | Texture cache efficiency (modifier vs optimal tiling) | 20–80× |
| 4K + anisotropic | Hardware anisotropic in texture unit | 15–60× |

The 4K headroom range is wide because it depends on whether the driver can serve Morton-layout reads from the texture unit's dedicated L1 (upper bound) or falls back to generic cache (lower bound). Either case is substantially better than the DX11 deswizzle path.

### CPU Overhead

The persistently-mapped BAR path reduces the dirty-page flush cost to a `memcpy` with no driver involvement. At 64MB with a typical 10–20% dirty fraction per frame:

```
6–12 MB × 60fps = 360–720 MB/s  ← well within L3 bandwidth on any modern CPU
```

Push buffer parsing on the GPU removes the NV097 method decode loop from the CPU entirely. For complex scenes with thousands of state changes per frame, this was the largest CPU cost in the DX11 path. In Vulkan, the parser CS handles it at GPU compute throughput.

---

## Extension Requirements

| Extension | Status | Purpose |
|---|---|---|
| `VK_KHR_buffer_device_address` | Core 1.2 | Device address for XboxRAM pointer |
| `VK_EXT_extended_dynamic_state3` | Core 1.3 / universal desktop | Replace PSO state objects |
| `VK_EXT_image_drm_format_modifier` | Universal desktop | VkImage over Morton memory |
| `VK_EXT_device_generated_commands` | Nvidia/AMD/Intel current | GPU push buffer replay |
| `VK_EXT_descriptor_buffer` | Universal desktop | GPU texture descriptor updates |
| `VK_KHR_dynamic_rendering` | Core 1.3 | Subpass-free render passes |
| `VK_NV_device_generated_commands` | Nvidia only (older) | Earlier push buffer replay |

All extensions marked "universal desktop" are available on RTX 20+, RX 5000+, and Intel Arc. `VK_EXT_device_generated_commands` became cross-vendor in 2024.

---

## Upscaling

### Internal Resolution

The same surface type distinction from the DX11 architecture applies: static textures are maintained at native Xbox dimensions; render target surfaces are created at IR dimensions. In Vulkan, render target `VkImage` resources are created at IR size and bound to `xboxMemory` only if the modifier path is not active for them — render targets use `VK_IMAGE_TILING_OPTIMAL` since their content is generated by GPU rendering, not read from Morton-ordered Xbox memory.

```
Xbox framebuffer:      640 × 480   (from PGRAPH)
Host VkImage RT:       2560 × 1920 (IRm = 4, VK_IMAGE_TILING_OPTIMAL)
Static texture VkImage: native Xbox dimensions (modifier path, same xboxMemory)
```

Point sprite size, line width, LOD bias, pixel-offset effects, and framebuffer readback at IR all behave identically to the DX11 analysis — these are API-independent problems. See the DX11 document for the detailed treatment of each. The sections below cover where the Vulkan architecture diverges.

### Jitter Injection: The GPU-Driven Problem

This is where the Vulkan path diverges most significantly from DX11. In the DX11 architecture, the CPU intercepts all VS constant writes and injects projection matrix jitter before the data reaches XboxRAM. In the fully GPU-driven Vulkan path, VS constants flow from the Xbox CPU directly into persistently-mapped `xboxMemory` with no CPU inspection. The CPU never sees the constant values — it only manages dirty page tracking.

Injecting jitter therefore requires the **push buffer parser CS** to identify and patch the projection matrix in-flight:

```glsl
// Push buffer parser CS — handling NV097_SET_TRANSFORM_CONSTANT_LOAD
case NV097_SET_TRANSFORM_CONSTANT_LOAD:
    uint constSlot = value;  // which VS constant slot is being loaded
    // Flag: we expect the next 16 DWORDs to be a matrix
    pendingConstLoad = constSlot;
    break;

case NV097_SET_TRANSFORM_CONSTANT:  // data following a LOAD method
    uint offset = method - NV097_SET_TRANSFORM_CONSTANT;
    uint slot   = pendingConstLoad + offset / 4;
    // Write constant to xboxMemory at VS_CONSTANTS_BASE + slot*16
    StoreConstant(slot, value);

    // After every 4th write (completing one float4), check if we just
    // finished writing a full 4×4 matrix and if it looks like a VP matrix
    if (IsMatrixComplete(slot) && IsProjectionMatrix(slot)) {
        float2 jitter = SampleHalton(frameIndex % 8);
        PatchProjectionJitter(slot, jitter);  // write jittered version back
        StoreJitter(jitter);                  // for motion vector PS use
    }
    break;
```

`IsProjectionMatrix` runs in the CS and checks the same heuristics as the DX11 CPU-side version: perspective divide row, rotation column magnitudes, near/far ratio. The CS has full read access to the constant data it just wrote into `xboxMemory` via the device address.

**The fundamental tension** is that the push buffer parser is a compute shader executing thousands of threads simultaneously, processing one or more DWORD per thread. Matrix detection requires correlating 16 sequential DWORDs — a cross-thread communication problem. The practical resolution is a two-pass approach:

- **Pass 1:** All threads write constants to `xboxMemory` unconditionally (fast, parallel)
- **Pass 2:** A single-thread CS scans completed constant writes for projection matrix candidates and applies jitter patches (slow, serial — but only ~8 matrix candidates per frame)

The jitter patch in Pass 2 writes back into `xboxMemory` at the identified slot. The VS subsequently reads the jittered version. The two-pass structure adds one CS dispatch but keeps the main parser fast.

**Alternative: hybrid CPU jitter for Vulkan.** For titles where the projection matrix slot is known (per-title database), the CPU can write the jitter directly into the persistently-mapped `xboxMemory` before the push buffer parser runs. This is a single `memcpy` of 64 bytes — negligible cost — and avoids the two-pass CS entirely. A per-title override table covers 95% of real games; the CS-based detection handles unknowns.

### Motion Vector Reconstruction

Motion vector reconstruction in Vulkan follows the same algorithm as DX11 — reprojecting world-space vertex positions through the previous frame's VP matrix in the VS:

```glsl
// layout(location = 8) out vec4 prevPos;

vec4 worldPos  = FetchVertex(gl_VertexIndex);
vec4 clip      = currentVP * worldPos;
vec4 prevClip  = prevVP    * worldPos;
gl_Position    = clip;
prevPosOut     = prevClip;
```

The `prevVP` matrix is maintained by the push buffer parser: at the end of each frame, the identified projection matrix slot is copied from `xboxMemory` into a `VkBuffer` used as a per-frame uniform. Since `xboxMemory` is persistently mapped, this copy is a direct 64-byte `memcpy` with no GPU readback.

**Skinned geometry limitation** is identical to DX11 — motion vectors are incorrect for VS-animated geometry. All three Vulkan-available temporal upscalers (FSR 2, FSR 3, DLSS 3) handle this gracefully.

### Available Upscalers in Vulkan

| Upscaler | Type | Motion Vectors | Vulkan Support | Notes |
|---|---|---|---|---|
| Bilinear blit | Spatial | No | `vkCmdBlit` | Zero cost; baseline |
| Integer scale | Spatial | No | Custom CS/FS | Nearest-neighbour; retro-correct pixels |
| AMD FSR 1 | Spatial | No | GLSL/SPIR-V port | Ships as HLSL; straightforward GLSL port |
| Nvidia NIS | Spatial | No | GLSL/SPIR-V port | Sharpening + upscale |
| AMD FSR 2 | Temporal | **Required** | Native Vulkan SDK | Best cross-vendor temporal quality |
| AMD FSR 3 | Temporal + frame gen | **Required** | Native Vulkan SDK | Frame generation available (unlike DX11) |
| Intel XeSS | Temporal | **Required** | Native Vulkan SDK | Good quality on all vendors |
| **Nvidia DLSS 3** | Temporal + frame gen | **Required** | Native Vulkan SDK | **Available — not possible in DX11** |

DLSS 3 is the primary advantage of the Vulkan upscaling path over DX11. On Nvidia hardware, DLSS 3 combines temporal upscaling with optical-flow-based frame generation (`VK_NV_optical_flow`), effectively doubling the output frame rate with minimal additional GPU cost. For an emulator rendering at modest IRm on modern hardware, DLSS 3 frame generation turns the large GPU headroom (see Performance section) into a battery/thermal saving rather than wasted cycles.

FSR 3 frame generation is also available in Vulkan — unlike DX11 where the frame generation component requires DX12. This makes cross-vendor frame generation accessible for the first time in this architecture.

### Jitter and Dynamic Rendering

In Vulkan with `VK_KHR_dynamic_rendering` (core 1.3), each render pass begins with `vkCmdBeginRendering` rather than a pre-allocated `VkRenderPass` object. This makes per-frame jitter integration slightly simpler — the render target's `VkImageView` can be a different mip or array layer per frame if needed for accumulation, without creating a new `VkRenderPass` object.

The upscaler's input descriptors (color, depth, motion vectors, reactive mask) are bound via `vkCmdPushDescriptorSetKHR` immediately before the upscale dispatch, using the frame's current VkImage handles. No descriptor pool management per frame.

### Reactive Mask

FSR 2, FSR 3, and DLSS 3 accept a reactive mask — a per-pixel scalar indicating surfaces that should be given less temporal weight (transparent, emissive, or reflective surfaces that change rapidly). Without it, temporal ghosting can occur on reflections and particles.

In the GPU-driven Vulkan path, the reactive mask can be written during the main render pass by having the combiner interpreter PS output a second attachment:

```glsl
// Secondary color attachment: reactive mask
layout(location = 1) out float reactiveMask;

void main() {
    // Primary output: combiner result
    fragColor = DoFinalCombiner(...);

    // Reactive: flag pixels using the EF_PROD or V1R0_SUM registers
    // (these appear in reflections and emissive surfaces in typical NV2A shaders)
    bool usesReflection = (PSFinalCombinerInputsABCD & REFLECTION_REGISTER_MASK) != 0u;
    reactiveMask = usesReflection ? 0.9 : 0.0;
}
```

The reactive mask attachment is declared in `VkPipelineRenderingCreateInfo` alongside the color and depth attachments. No separate render pass is needed.

### 2D / HUD Pass Detection in Vulkan

Detection is identical to DX11 — the programmable VS enable flag in the push buffer stream identifies 2D passes. In the GPU-driven path, the push buffer parser CS detects the `NV097_SET_TRANSFORM_EXECUTION_MODE` method and emits a token that switches the VS permutation to the passthrough variant with `IR_SCALE` applied (Option A) or redirects rendering to a native-resolution attachment (Option B).

The native-resolution HUD attachment in Option B is a separate `VkImage` at 640×480. Its `VkImageView` is swapped in via `vkCmdBeginRendering` parameters when the 2D pass is active. The push buffer parser identifies the transition and emits a `vkCmdEndRendering` / `vkCmdBeginRendering` pair to switch attachments. This is cleaner than the DX11 equivalent because `VK_KHR_dynamic_rendering` allows attachment changes without pre-defined subpass structures.

---

## Texture Replacement

### Hashing in the Vulkan Architecture

The hashing approach is identical in concept to DX11 — xxHash64 over raw XboxRAM bytes at the texture's surface address. In Vulkan with persistently-mapped `xboxMemory`, the CPU can compute the hash directly from the mapped pointer with no staging buffer:

```cpp
uint64_t HashTexture(uint32_t xboxAddr, uint32_t byteSize) {
    return XXH64((uint8_t*)mappedXboxMemory + xboxAddr, byteSize, 0);
}
```

This runs on a CPU thread with no GPU synchronisation — the mapped pointer is always readable by the CPU. The hash is recomputed only when the texture's dirty page bit is set.

### Replacement and the VkImage Modifier Path

The VkImage-over-shared-memory modifier path (where a `VkImage` is bound to `xboxMemory` at the texture's Xbox address) cannot be used for replaced textures — the replacement asset has different content and likely different dimensions than the Xbox original. When a replacement is active:

1. The modifier-path `VkImage` for that address is not created (or is destroyed if it existed)
2. A separate `VkImage` is created at the replacement asset's resolution with `VK_IMAGE_TILING_OPTIMAL`
3. The replacement asset is uploaded via a staging buffer into this `VkImage`
4. The texture pool binds this `VkImage` for the affected surface address

For non-replaced textures, the modifier path remains active — no copy, no separate `VkImage`. The pool maintains both paths simultaneously: per surface address, either a modifier-bound `VkImage` (native) or a standalone `VkImage` (replacement).

### Descriptor Updates for Replacement

When a replacement is swapped in, the descriptor pointing to the old `VkImage` must be updated before the next draw. With `VK_EXT_descriptor_buffer`, the descriptor is a GPU-buffer entry that the replacement loader writes directly:

```cpp
// Replacement loader thread (background)
void OnReplacementReady(uint64_t hash, VkImageView replacementView) {
    // Write new descriptor into descriptor buffer
    VkDescriptorImageInfo imageInfo = { .imageView = replacementView,
                                        .imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL };
    WriteDescriptorToBuffer(textureDescriptorBuffer, slotForHash(hash), imageInfo);
    // Atomic flag for main thread: descriptor for this hash is now live
    descriptorReady[hash] = true;
}
```

The main thread checks `descriptorReady` before each draw that uses a recently-replaced texture. The swap is atomic from the GPU's perspective — the descriptor buffer write is complete before the next `vkQueueSubmit`.

### Render Target Surfaces

Identical to DX11: render target surfaces skip the hash lookup. The `isRenderTarget` flag is set when the surface address appears in a push buffer `NV097_SET_SURFACE_*` method, detected by the push buffer parser CS and stored in a per-surface metadata buffer.

### Async Loading Integration with Dynamic Rendering

Texture replacement assets load on background threads. The native texture (modifier path or deswizzle fallback) is used on first encounter. When the replacement is ready, the descriptor buffer is updated. Because `VK_KHR_dynamic_rendering` requires no pre-compiled render pass objects, the descriptor update does not invalidate any pipeline cache entries — the new `VkImage` is simply referenced by the updated descriptor on the next draw that samples it.

---

## Known Remaining Items

These items exist in both the DX11 and Vulkan architectures and require additional implementation work independent of the API choice:

- **PSCompareMode** — clip-plane comparison mode in the combiner interpreter
- **PSDotMapping** — normal mapping for dot-product texture modes
- **PSInputTexture** — dependent-texture input routing
- **DOT\_RFLCT\_SPEC\_CONST** eye vector — requires wiring `SetEyeVector()` / `D3DRS_PSINPUTTEXTURE`
- **Full BUMPENVMAP** perturbation coordinate routing
- **Full BRDF mode** — requires eye and light sigma vector inputs
- **Framebuffer readback** — GPU→CPU path for Xbox CPU surface reads, triggered on `hostDirty` pages

---

## Architectural Comparison: DX11 vs Vulkan

| Component | DX11 | Vulkan |
|---|---|---|
| XboxRAM upload | `UpdateSubresource` (driver copy) | `memcpy` into BAR-mapped memory |
| Texture access | Deswizzle CS → `Texture2D` pool | `VkImage` over Xbox memory (modifier) |
| Pipeline state | PSO creation + cache | ~10 `vkCmd*` calls, no PSO per state |
| Push buffer | CPU parses NV097, issues draws | CS parses, GPU executes via `GeneratedCommands` |
| Draw call | CPU calls `Draw(hostVertexCount)` | `vkCmdDrawIndirectCount` from GPU buffer |
| Vertex fetch | VS from XboxRAM SRV | VS from XboxRAM device address |
| Topology convert | In-shader `SV_VertexID` arithmetic | In-shader `gl_VertexIndex` arithmetic |
| CMP decode | Inline in VS | Inline in VS |
| Combiner interpreter | SM5 HLSL, Steps 2+3+4 | GLSL, specialisation constants |
| CPU per-frame calls | ~dozens | ~6 |
