# NV2A LLE Migration Plan: From HLE D3D Patches to Pushbuffer-Driven Rendering

## Executive Summary

Migrate Cxbx-Reloaded from the current **HLE D3D patch + D3D11 backend** architecture to a
**pushbuffer-driven NV2A LLE** architecture where:

1. Xbox D3D API calls execute **unpatched** (original Xbox D3D library writes pushbuffer commands)
2. A new NV2A PGRAPH engine **consumes pushbuffer commands directly** (no VEH/MMIO overhead)
3. RC and VS interpreters read from **real NV2A PGRAPH register state**
4. The host GPU backend migrates from **D3D11 → Vulkan**
5. Ultimately, pushbuffer command handling runs as a **Vulkan compute shader** on the host GPU

---

## Current Architecture (Baseline)

```
┌──────────────────────────────────────────────────────────┐
│ Xbox Game (.xbe)                                         │
│   calls Xbox D3D8 API (SetRenderState, DrawVertices...)  │
└──────────────────────┬───────────────────────────────────┘
                       │ subhook intercept
                       ▼
┌──────────────────────────────────────────────────────────┐
│ EMUPATCH Layer (150+ D3D patches)                        │
│   - Calls XB_TRMP() → original Xbox D3D code            │
│   - Mirrors Xbox state into HLE structures               │
│   - Converts resources to host D3D11 format              │
│   - Issues D3D11 draw calls directly                     │
└──────────────────────┬───────────────────────────────────┘
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
┌──────────────┐ ┌──────────┐ ┌───────────────────┐
│ HLE State    │ │ D3D11    │ │ HLSL Shaders      │
│ Mirrors:     │ │ Backend  │ │ - VS Interpreter  │
│ RenderState  │ │ Draw()   │ │ - RC Interpreter  │
│ TextureState │ │ Present()│ │ - Vertex Fetch    │
│ Globals      │ │          │ │ - Fixed-Function  │
└──────────────┘ └──────────┘ └───────────────────┘

Meanwhile, NV2A pushbuffer IS PROCESSED by PFIFO puller thread into PGRAPH regs[],
but HLE draw path reads from HLE state mirrors — NOT from PGRAPH.
MMIO access trapped by VEH → 10-15 layer dispatch (distorm decode per access).

Note: PFIFO puller was previously draining/discarding commands. It now dispatches
methods >= 0x100 to pgraph_handle_method(), populating regs[]. However, this state
is NOT consumed by the D3D11 rendering path due to a race condition: HLE draw calls
happen before the puller has finished processing the pushbuffer for that frame.
```

## Target Architecture

```
┌──────────────────────────────────────────────────────────┐
│ Xbox Game (.xbe)                                         │
│   calls Xbox D3D8 API (unpatched, original code)         │
│   Xbox D3D library writes NV2A pushbuffer commands       │
└──────────────────────┬───────────────────────────────────┘
                       │ DMA to NV2A push buffer memory
                       ▼
┌──────────────────────────────────────────────────────────┐
│ PFIFO Pusher (existing thread, optimized)                │
│   reads push buffer, populates CACHE1                    │
└──────────────────────┬───────────────────────────────────┘
                       ▼
┌──────────────────────────────────────────────────────────┐
│ PFIFO Puller → PGRAPH Method Handler (rewritten)         │
│   - Executes NV097 methods, updates PGRAPH register state│
│   - No OpenGL path — pure state machine                  │
│   - Triggers draw on NV097_SET_BEGIN_END(0)              │
└──────────────────────┬───────────────────────────────────┘
                       │ draw trigger
                       ▼
┌──────────────────────────────────────────────────────────┐
│ Vulkan Rendering Backend                                 │
│   - Reads PGRAPH state (regs[], vertex_attributes, etc.) │
│   - VS Interpreter reads xf.xfpr[] (XFPR RAM)            │
│                       + xf.xfctx[] (XFCTX RAM)           │
│   - RC Interpreter reads combiner regs from PGRAPH       │
│   - Vulkan compute shader for pushbuffer processing      │
│     (ultimate goal)                                      │
└──────────────────────────────────────────────────────────┘
```

---

## Phase 0: Preparation & Infrastructure

### 0.1 — Create Feature Branch & Build Infrastructure
- Branch from current `dx11` or `main`
- Ensure both old (HLE) and new (LLE) paths can coexist via runtime flag
- Add `NV2A_LLE_RENDERING` compile-time and runtime toggle
- Update CMakeLists.txt for Vulkan SDK dependency (Phase 4+)

### 0.2 — Audit & Document Current NV2A PGRAPH State
- **File**: `src/devices/video/nv2a_int.h` — `PGRAPHState` struct
- Catalog every register field and its Xbox hardware meaning
- Cross-reference with NVIDIA NV20/NV2A documentation and xemu's pgraph implementation
- Create a register map document: offset → name → bit fields → semantics
- **Deliverable**: Complete NV2A PGRAPH register reference

### 0.3 — Build NV2A Pushbuffer Command Trace/Replay Tool
- Capture all pushbuffer commands from known-good XDK samples
- Log: `{method, subchannel, parameter}` tuples with timestamps
- Build offline replay tool to feed captured commands to PGRAPH
- Use this for regression testing throughout all phases

---

## Phase 1: Rewrite PGRAPH Command Processing (Remove OpenGL, Remove MMIO Overhead)

This is the foundation. The current `EmuNV2A_PGRAPH.cpp` has dual HLE/LLE paths and OpenGL
rendering interleaved with state management. We need a clean state machine.

### 1.1 — Extract PGRAPH State Machine from OpenGL Rendering
- **Current file**: `src/devices/video/EmuNV2A_PGRAPH.cpp` (~2000+ lines)
- Split into:
  - `NV2A_PGRAPH_State.cpp` — Pure NV2A register state machine (no rendering)
  - `NV2A_PGRAPH_Render.cpp` — Host GPU rendering dispatch (Phase 3+)
- For each NV097 method, extract the **state update** portion:
  - `NV097_SET_SURFACE_*` → update `surface_color`, `surface_zeta`, `surface_shape`
  - `NV097_SET_VERTEX_DATA_ARRAY_*` → update `vertex_attributes[]`
  - `NV097_SET_TRANSFORM_PROGRAM` → update `xf.xfpr[][]` (XFPR RAM mirror)
  - `NV097_SET_TRANSFORM_CONSTANT` → update `xf.xfctx[][]` (XFCTX RAM mirror)
  - `NV097_SET_COMBINER_*` → update combiner registers in `regs[]`
  - `NV097_SET_TEXTURE_*` → update texture state per unit
  - `NV097_SET_BLEND_*`, `NV097_SET_DEPTH_*`, `NV097_SET_STENCIL_*` → update `regs[]`
- The OpenGL draw calls (`glDrawArrays`, `glDrawElements`, shader compilation) are **removed**

### 1.2 — Remove OpenGL LLE Backend Entirely
- **Delete or ifdef-out**:
  - `src/devices/video/nv2a_vsh.cpp` (ARB vertex program translator)
  - `src/devices/video/nv2a_psh.cpp` (GLSL register combiner translator)
  - `src/devices/video/nv2a_shaders.cpp` (OpenGL shader utilities)
  - All `#include <GL/...>` and GLEW references in `nv2a.cpp`, `EmuNV2A_PGRAPH.cpp`
  - `GloContext`, `gl_framebuffer`, `gl_color_buffer`, `gl_zeta_buffer` from `PGRAPHState`
  - `opengl_enabled` flag (replaced by new architecture flag)
- **Keep**: All NV2A register definitions, method constants, push buffer parsing

### 1.3 — Bypass MMIO VEH for NV2A Register Access (Direct Call Path)
- **Current overhead**: Page fault → VEH → distorm x86 decode → PCI bus → block lookup → handler
- **New approach**: Map NV2A MMIO region as **read/write accessible**
  - Back it with actual memory pages containing register values
  - Use `VirtualAlloc` at `0xFD000000` with `PAGE_READWRITE`
  - Xbox code reads/writes NV2A registers directly (fast, no exceptions)
- **Challenge**: Some registers have **side effects** on write (e.g., interrupt acknowledge, FIFO kick)
- **Solution**: Use **write-watch pages** (`PAGE_GUARD` or `MEM_WRITE_WATCH`) only for
  side-effect registers, not for all MMIO
  - Most PGRAPH registers (0xFD400000-0xFD401FFF) are plain state → direct access
  - PFIFO control registers → keep guarded (trigger pusher/puller)
  - PMC interrupt registers → keep guarded (interrupt dispatch)
  - PCRTC registers → keep guarded (vblank timing)
- **Alternative**: For the pushbuffer path, Xbox D3D doesn't read MMIO registers often — it
  writes to the push buffer and kicks PFIFO via USER channel registers. So the hot path
  (USER writes) can be handled with minimal overhead.

### 1.4 — Enable PFIFO Puller to Execute PGRAPH Methods in HLE Mode
- **Current behavior**: Puller dispatches methods >= 0x100 to `pgraph_handle_method()`
  which populates `PGRAPHState.regs[]` — **this is already implemented**
- **Remaining issue**: HLE draw calls execute before the puller finishes processing
  pushbuffer commands (race condition). The RC interpreter was reverted to read from
  PSDef (HLE state) because PGRAPH registers were all zeros at draw time.
- **Solution**: Either (a) flush/synchronize the puller before each HLE draw, or
  (b) move draw triggering to the puller thread itself (Phase 3 target)
- Draw function pointers (`pgraph_draw_arrays` etc.) remain `nullptr` in HLE mode,
  so draw-triggering methods in `pgraph_handle_method` are safe no-ops

### 1.5 — Optimize PFIFO Command Fetch
- Current pusher thread reads 32-bit words one at a time from Xbox memory
- Batch-read: `memcpy` entire push buffer segment (GET→PUT) into local buffer
- Process locally without per-word memory access overhead
- Handle jump/call commands for buffer chaining

---

## Phase 2: Rewire RC/VS Interpreters to Read PGRAPH State

Currently, the HLSL interpreters read from HLE-maintained Xbox state mirrors. They need to
read from the NV2A PGRAPH register state populated by Phase 1.

### 2.1 — RC Interpreter: Switch from Xbox RenderState to PGRAPH Registers
- **Current source**: `CxbxD3D11UploadRCInterpreterState()` in `Backend_D3D11_PixelShader.cpp`
  reads `X_D3DPIXELSHADERDEF` from `XboxRenderStates.GetPixelShaderRenderStatePointer()`
  (HLE state). An earlier attempt to read from PGRAPH `regs[]` was **reverted** due to
  the puller race condition (PGRAPH registers all zeros at draw time).
- **New source**: Read combiner registers directly from `PGRAPHState.regs[]`:
  - `NV_PGRAPH_COMBINECOLORI0..7` → RGB stage inputs
  - `NV_PGRAPH_COMBINECOLORO0..7` → RGB stage outputs
  - `NV_PGRAPH_COMBINEALPHAI0..7` → Alpha stage inputs
  - `NV_PGRAPH_COMBINEALPHAO0..7` → Alpha stage outputs
  - `NV_PGRAPH_COMBINEFACTOR0..1` → C0/C1 constants per stage
  - `NV_PGRAPH_COMBINESPECFOG0/1` → Final combiner inputs
  - `NV_PGRAPH_SHADERCTL` → Texture modes (project, cubemap, dot, bumpenv)
  - `NV_PGRAPH_SHADERPROG` → Shader program config
  - `NV_PGRAPH_TEXCTL0_0..3` → Texture enable/format per unit
  - `NV_PGRAPH_TEXFMT0..3` → Texture format details
- **Constant buffer layout change**: `RCInterpreterCBLayout` struct updated to map
  directly from PGRAPH register offsets
- **The HLSL interpreter code** (`CxbxRegisterCombinerInterpreter.hlsl`) may need minimal
  changes — the register combiner logic is the same, only the source of constants changes

### 2.2 — VS Interpreter: Switch from Xbox VertexShader Slots to PGRAPH Program Data  ✅ DONE
- **Current source**: `CxbxD3D11UploadVSInterpreterState()` in `Backend_D3D11_VertexShader.cpp`
  reads raw NV2A vertex shader microcode from HLE-cached function slots
  (`GetCxbxVertexShaderSlotPtr(g_Xbox_VertexShader_FunctionSlots_StartAddress)`).
  Uploads up to 136 × 4 DWORDs into `VSInterpreterCBLayout` at `b0` (2192 bytes).
- **New source**: Read directly from `PGRAPHState`:
  - `xf.xfpr[NV2A_MAX_TRANSFORM_PROGRAM_LENGTH][4]` — XFPR RAM mirror:
    NV2A Transform Program RAM, on-chip XF SRAM with 136 × 92-bit instructions
    in 128-bit containers. Uploaded via `NV097_SET_TRANSFORM_PROGRAM` with
    `NV_PGRAPH_CHEOPS_OFFSET.PROG_LD_PTR` as auto-incrementing write pointer.
    The RDI (Register Direct Interface) is used for context save/restore.
  - `xf.xfctx[NV2A_VERTEXSHADER_CONSTANTS][4]` — XFCTX RAM mirror:
    NV2A Transform Context RAM, also on-chip XF SRAM behind the RDI interface.
    Holds 192 × float4 constant registers.
  - `NV_PGRAPH_CSV0_D` / `NV_PGRAPH_CSV0_C` — VS program start address, mode bits
- **Constant buffer upload**: `VSInterpreterCBLayout` now populated from PGRAPH state
- **Fixed-function pipeline**: When no vertex shader program is active, the NV2A uses
  a fixed-function transform pipeline. The existing `CxbxFixedFunctionVertexShaderState.hlsli`
  handles this, but state source switches to PGRAPH lighting/transform registers:
  - `NV_PGRAPH_XF_MODELVIEW0..3` → Model-view matrices
  - `NV_PGRAPH_XF_PROJECTION0..3` → Projection matrix
  - `NV_PGRAPH_LIGHT0..7_*` → Light source parameters

### 2.3 — Vertex Attribute Fetch: Switch to PGRAPH Vertex Array State
- **Current**: `Backend_D3D11_VertexFetch.cpp` reads from `g_Xbox_SetStreamSource[]` (HLE state)
- **New**: Read from `PGRAPHState.vertex_attributes[16]`:
  - `dma_select`, `offset`, `format`, `size`, `count`, `stride` per attribute
  - `inline_value[4]` for NV2A sticky/default attribute values
  - `dma_vertex_a` / `dma_vertex_b` for DMA context (physical address base)
- **Cleanup required**: Remove GL-specific fields from `VertexAttribute` struct:
  - `gl_count`, `gl_type`, `gl_normalize`, `gl_converted_buffer`, `gl_inline_buffer`
  - `needs_conversion`, `converted_buffer`, `converted_elements`, `converted_size`
- Vertex data still fetched from Xbox memory via SRVs, but address/format comes from PGRAPH

### 2.4 — Surface/Render Target State from PGRAPH
- **Current**: `g_pXbox_RenderTarget`, `g_pXbox_DepthStencil` (HLE globals)
- **New**: Read from PGRAPH and NV097 methods:
  - `NV097_SET_SURFACE_FORMAT` (0x0208) → color/zeta format
  - `NV097_SET_SURFACE_COLOR_OFFSET` (0x0210) → color buffer address in VRAM
  - `NV097_SET_SURFACE_ZETA_OFFSET` (0x0214) → depth buffer address
  - `NV097_SET_SURFACE_PITCH` (0x020C) → row pitch
  - `NV_PGRAPH_SURFACE` (MMIO 0x0710) → surface type/shape
  - Surface clip from `NV097_SET_SURFACE_CLIP_HORIZONTAL` / `VERTICAL`

### 2.5 — Remaining Render State from PGRAPH
- Blend state: `NV_PGRAPH_BLEND`, `NV_PGRAPH_BLENDCOLOR`
- Depth/Stencil: `NV_PGRAPH_CONTROL_0` (depth test), `NV_PGRAPH_CONTROL_1` (stencil)
- Rasterizer: `NV_PGRAPH_SETUPRASTER` (cull mode, polygon mode, point size)
- Fog: `NV_PGRAPH_FOGCOLOR`, `NV_PGRAPH_FOGPARAM0/1`
- Alpha test: within `NV_PGRAPH_CONTROL_0`
- Viewport: `NV_PGRAPH_WINDOWCLIPX0..7`, `NV_PGRAPH_WINDOWCLIPY0..7`
- Scissor: `NV_PGRAPH_SETUPRASTER` scissor bits

---

## Phase 3: Remove D3D EMUPATCH Layer (Incremental)

With the pushbuffer-driven PGRAPH state machine active and interpreters reading from it,
EMUPATCHes become redundant. Remove them in dependency order.

### 3.1 — Remove State-Setting Patches (No Host GPU Side Effects)
These patches currently mirror Xbox state for the interpreters. Once interpreters read
from PGRAPH, the original Xbox D3D code does this work via pushbuffer:

**Batch 1 — Pure state patches (safest to remove first)**:
- `D3DDevice_SetRenderState_Simple` and all `SetRenderState_*` variants
- `D3DDevice_SetVertexShaderConstant` / `SetVertexShaderConstantNotInline` / `1` / `4`
- `D3DDevice_SetPixelShader` (RC state now comes from PGRAPH)
- `D3DDevice_SetTexture` (texture binding now from PGRAPH)
- `D3DDevice_SetTransform` / `MultiplyTransform`
- `D3DDevice_SetStreamSource` / `SetIndices`
- `D3DDevice_SetVertexShader` / `SelectVertexShader`
- `D3DDevice_SetViewport`
- `D3DDevice_SetScissors`

**Batch 2 — Resource management patches**:
- `D3DDevice_SetPalette`
- `D3DDevice_SetVertexData*` (inline vertex data → NV097_SET_VERTEX* methods)
- `D3DDevice_LoadVertexShader` / `DeleteVertexShader`
- `D3DDevice_SetGammaRamp`

### 3.2 — Remove Draw Call Patches
Once PGRAPH processes `NV097_SET_BEGIN_END` and triggers host draw calls:
- `D3DDevice_DrawVertices` / `DrawVerticesUP`
- `D3DDevice_DrawIndexedVertices` / `DrawIndexedVerticesUP`
- `D3DDevice_Begin` / `End` (inline primitive batching)
- `D3DDevice_DrawRectPatch` / `DrawTriPatch`

### 3.3 — Remove Synchronization Patches
Requires proper GPU event/fence emulation (see Phase 3.5):
- `D3DDevice_BlockOnFence` / `InsertFence` / `IsFencePending`
- `D3DDevice_BlockOnTime`
- `D3DDevice_BeginVisibilityTest` / `EndVisibilityTest` / `GetVisibilityTestResult`
- `D3DResource_BlockUntilNotBusy`

### 3.4 — Remove Presentation Patches
- `D3DDevice_Present` / `Swap` / `PersistDisplay`
- `D3DDevice_EnableOverlay` / `UpdateOverlay`
- These get replaced by PCRTC/PVIDEO handling: read framebuffer from VRAM at vblank

### 3.5 — Implement NV2A GPU Event & Fence Support
The NV2A supports:
- **Semaphore release**: `NV097_SET_SEMAPHORE_OFFSET` + `NV097_BACK_END_WRITE_SEMAPHORE_RELEASE`
  writes a value to memory when all prior commands complete
- **Report queries**: `NV097_GET_REPORT` writes pixel count / Zpass count to memory
- **Notifiers**: DMA notifier objects for command completion
- Implement these in the PGRAPH state machine to unblock removal of sync patches

### 3.6 — Retain Essential Non-GPU Patches
Some patches handle CPU-side concerns that don't go through the pushbuffer:
- `D3DDevice_GetCreationParameters` (info query)
- `D3DDevice_GetDisplayMode` (info query)
- Resource `Lock`/`Unlock` (CPU memory access — may still need interception for
  page-tracking of GPU-readable memory)
- `D3D_CreateDevice` (initialization — always needed for bootstrapping)
- `D3D_KickOffAndWaitForIdle` (CPU-side busy wait — redirect to fence check)

---

## Phase 4: Vulkan Backend

### 4.1 — Vulkan Device Initialization
- Replace D3D11 device creation in `HostDevice.cpp`
- Create `VkInstance`, `VkDevice`, `VkQueue`, `VkSwapchainKHR`
- Implement `VkSurfaceKHR` for the existing Win32 emulator window
- Allocate device memory for VRAM mirror (host-visible, device-local)

### 4.2 — Port HLSL Shaders to SPIR-V
- Use DXC (DirectX Shader Compiler) or glslang to compile HLSL → SPIR-V
  - DXC supports HLSL → SPIR-V directly via `-spirv` flag
  - Most HLSL constructs map cleanly; `ByteAddressBuffer` → SSBO
- Port key shaders:
  - `CxbxVertexShaderInterpreter.hlsl` → SPIR-V compute/vertex shader
  - `CxbxRegisterCombinerInterpreter.hlsl` → SPIR-V fragment shader
  - `CxbxVertexFetch.hlsli` → SPIR-V vertex shader (buffer device address for vertex fetch)
- **Alternative**: Keep HLSL and use `spirv-cross` or DXC for runtime compilation

### 4.3 — Vulkan Render Pass & Pipeline Setup
- Create render passes matching NV2A surface formats:
  - Color: R8G8B8A8, R5G6B5, X1R5G5B5, etc.
  - Depth: D16, D24S8, D24X8, F16, F24S8
- Create pipeline objects for each shader combination
- Implement pipeline cache for shader permutations
- Use dynamic state for viewport, scissor, blend constants, stencil ref

### 4.4 — Vulkan Resource Management
- **VRAM Mapping**: Xbox 64MB VRAM → `VkBuffer` with `VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT`
  - Xbox GPU reads textures/vertices directly from this
  - Use `VkBufferView` for typed access
- **Texture Management**: Create `VkImage`/`VkImageView` on demand from PGRAPH surface state
  - Swizzle conversion (Xbox tile/swizzle → linear) via compute shader
  - Format conversion for Xbox-specific formats (DXT, YUY2, palette)
- **Render Targets**: `VkImage` backed by VRAM regions (color/zeta offsets)
  - Aliased with VRAM buffer for CPU readback

### 4.5 — Vulkan Draw Dispatch
- On `NV097_SET_BEGIN_END(0)` (draw trigger):
  1. Read PGRAPH state → determine pipeline configuration
  2. Bind descriptor sets (vertex SRV, texture samplers, constant buffers)
  3. `vkCmdDraw()` or `vkCmdDrawIndexed()`
  4. Submit command buffer segment
- Frame presentation:
  1. On PCRTC vblank → acquire swapchain image
  2. Blit from VRAM framebuffer → swapchain image
  3. `vkQueuePresentKHR()`

### 4.6 — GPU Synchronization via Vulkan
- `VkFence` for CPU-GPU sync (replaces D3D11_QUERY_EVENT)
- `VkSemaphore` for frame pacing
- `VkEvent` for NV2A semaphore emulation
- Timeline semaphores for fine-grained command ordering

---

## Phase 5: GPU-Side Pushbuffer Processing (Compute Shader)

The ultimate goal: process NV2A pushbuffer commands **on the host GPU** via a Vulkan compute
shader, eliminating the CPU-side puller thread entirely.

### 5.1 — Design GPU Command Processor
- **Compute shader** reads push buffer from shared VRAM buffer
- Decodes command headers: method, subchannel, count, type
- Updates PGRAPH register state stored in a **GPU-side SSBO**
- On draw-trigger methods, emits indirect draw commands

### 5.2 — GPU PGRAPH State Buffer
- Allocate `VkBuffer` (SSBO) for entire PGRAPH register file (~8KB)
- Compute shader writes method parameters to correct register offsets
- Vertex/Pixel shaders read this SSBO for rendering state
- Use `VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT` for maximum bandwidth

### 5.3 — Indirect Draw Dispatch
- Compute shader populates `VkDrawIndirectCommand` / `VkDrawIndexedIndirectCommand`
  structs in a GPU buffer
- After command processing, issue `vkCmdDrawIndirect()` consuming the indirect buffer
- This chains GPU work: command processing → state setup → draw → present

### 5.4 — Synchronization Between Compute and Graphics
- Use `VkEvent` or pipeline barriers between compute (command processing) and
  graphics (rendering) stages
- Memory barriers ensure PGRAPH SSBO writes are visible to vertex/fragment shaders
- `VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT` → `VK_PIPELINE_STAGE_VERTEX_INPUT_BIT`

### 5.5 — Handle Non-Trivial Methods on GPU
Some NV2A methods have complex side effects:
- **DMA object resolution**: Context DMA lookups (RAMHT) may still need CPU assist
- **Surface format changes**: May require pipeline rebinding (GPU can't do this)
- **Interrupt generation**: Must notify CPU (via readback + CPU check)
- **Solution**: Hybrid approach — GPU handles simple register writes (90%+ of commands),
  CPU handles complex side effects via a "deferred command" queue

### 5.6 — Performance Optimization
- Batch pushbuffer segments into single compute dispatch
- Use subgroup operations for parallel command decoding
- Minimize GPU↔CPU synchronization points
- Profile and optimize hot paths (texture setup, vertex array config, combiner state)

---

## Phase 6: Testing & Validation

### 6.1 — Per-Phase Regression Testing
- **Test suite**: XDK samples (BumpLens, BumpEarth, Fur, Lighting, Tiling, etc.)
- **Golden reference**: Capture frame output from current HLE mode
- **Per-phase comparison**: Pixel-diff rendered frames against golden reference
- **Automated**: CI pipeline runs test suite after each merge

### 6.2 — Pushbuffer Trace Comparison
- Record pushbuffer traces from test titles
- Compare PGRAPH register state after trace replay against expected values
- Validate: surface state, combiner config, VS constants, texture setup

### 6.3 — Commercial Title Testing
- Halo, JSRF, Panzer Dragoon Orta, Ninja Gaiden, Dead or Alive 3
- Focus on:
  - Vertex shader programs (complex transforms, skinning)
  - Register combiner setups (multi-texture, bump mapping, specular)
  - GPU sync (fences, visibility tests)
  - Render target switches
  - Video overlay (cutscenes)

### 6.4 — Performance Benchmarking
- Measure frame time per phase
- Target: pushbuffer-driven path ≥ parity with current HLE
- GPU compute command processing should be strictly faster than CPU puller

---

## Dependency Graph

```
Phase 0 (Prep)
    │
    ▼
Phase 1 (PGRAPH State Machine)  ◄── Foundation for everything
    │
    ├──► Phase 2 (Rewire Interpreters to PGRAPH)
    │        │
    │        ▼
    │    Phase 3 (Remove EMUPATCHes)  ◄── Incremental, per-batch
    │
    ▼
Phase 4 (Vulkan Backend)  ◄── Can start in parallel with Phase 2-3 for basic bring-up
    │
    ▼
Phase 5 (GPU Compute Command Processor)  ◄── Requires Phase 4 complete
    │
    ▼
Phase 6 (Testing)  ◄── Continuous throughout, intensifies at end
```

---

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| Xbox D3D internal state diverges from NV2A pushbuffer | High | Trace-compare tool; XDK sample corpus |
| NV2A register semantics undocumented/wrong | High | Cross-reference xemu, MAME NV2A, envytools |
| VEH removal breaks non-GPU MMIO (APU, USB, etc.) | Medium | Only change NV2A MMIO; keep VEH for other devices |
| Vulkan driver differences across vendors | Medium | Test AMD/NVIDIA/Intel; use Vulkan validation layers |
| GPU compute pushbuffer has insufficient indirect draw support | Medium | Hybrid CPU/GPU fallback for complex methods |
| Performance regression during transition | Low | Keep HLE path as fallback; runtime toggle |
| Some games rely on EMUPATCH side effects | Medium | Careful per-patch analysis before removal; whitelist |

---

## Key References

- **xemu** (Xbox emulator with NV2A LLE): https://github.com/xemu-project/xemu
  - Most complete open-source NV2A implementation
  - PGRAPH method handling, register definitions
- **envytools**: https://envytools.readthedocs.io/
  - NVIDIA GPU documentation (NV20 family closely related to NV2A)
- **Mesa Nouveau driver**: NV20/NV2x register definitions and programming model
- **NVIDIA NV_register_combiners spec**: https://www.opengl.org/registry/specs/NV/register_combiners.txt
- **NVIDIA NV_vertex_program spec**: https://www.opengl.org/registry/specs/NV/vertex_program1_1.txt

---

## Estimated Effort Breakdown

| Phase | Scope | Complexity |
|-------|-------|------------|
| Phase 0 | Infrastructure, docs, tooling | Low |
| Phase 1 | PGRAPH rewrite, OpenGL removal, MMIO optimization | **High** |
| Phase 2 | Interpreter rewiring | Medium |
| Phase 3 | EMUPATCH removal (incremental) | Medium-High |
| Phase 4 | Vulkan backend | **High** |
| Phase 5 | GPU compute command processor | **Very High** |
| Phase 6 | Testing & validation | Ongoing |

Phase 1 and Phase 4 are the critical-path items with the highest risk and effort.
Phase 5 is the most ambitious but also the most novel contribution to Xbox emulation.
