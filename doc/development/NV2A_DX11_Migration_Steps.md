# NV2A DX11 Migration: Actionable Step-by-Step Plan

## Goal

Migrate the current HLE D3D patch + D3D11 backend to a pushbuffer-driven PGRAPH
architecture on the existing D3D11 backend. Remove all legacy OpenGL LLE code.
Prepare the codebase for a future Vulkan backend without writing any Vulkan code.

Each step is independently testable. The emulator must remain fully functional after
every step — XDK samples (BumpEarth, BumpLens, Fur, Dolphin, PixelShader, etc.) are
the validation corpus.

---

## Current State Summary

| Component | Status |
|-----------|--------|
| PFIFO puller → PGRAPH `regs[]` | **Active** — methods ≥ 0x100 dispatched, state populated |
| RC interpreter reads from | **HLE PSDef** (PGRAPH attempt reverted — race condition) |
| VS interpreter reads from | **HLE shader slot cache** |
| Vertex fetch reads from | **HLE `g_Xbox_SetStreamSource[]`** |
| Surface/RT state from | **HLE globals** (`g_pXbox_RenderTarget`, etc.) |
| Pipeline state from | **HLE `XboxRenderStates` / `XboxTextureStates`** |
| Draw trigger | **HLE EMUPATCH** draw functions |
| OpenGL LLE path | **Removed** (Step 1 complete) |
| OpenGL types in structs | **Removed** from `PGRAPHState`, `VertexAttribute` |
| GLSL shader translators | **Removed** |
| gloffscreen library | **Removed** |
| GLEW import | **Removed** |
| PFIFO→PGRAPH flush | **Active** — `pfifo_flush_to_pgraph()` called before each HLE draw |

---

## Step 1: Remove OpenGL LLE Rendering Backend  ✅ DONE

**Completed.** All OpenGL/GLEW code removed. 60 files changed, ~38K lines deleted.
Build verified, XDK samples render identically.

### 1.1 — Remove OpenGL draw function implementations

**Files:** `src/devices/video/EmuNV2A_PGRAPH.cpp`

- Delete all `OpenGL_draw_*` functions and `pgraph_init_OpenGL()`:
  - `pgraph_draw_arrays_gl()` / `OpenGL_draw_arrays()`
  - `pgraph_draw_inline_buffer_gl()` / `OpenGL_draw_inline_buffer()`
  - `pgraph_draw_inline_array_gl()` / `OpenGL_draw_inline_array()`
  - `pgraph_draw_inline_elements_gl()` / `OpenGL_draw_inline_elements()`
  - `pgraph_draw_state_update_gl()` / `OpenGL_draw_state_update()`
  - `pgraph_draw_clear_gl()` / `OpenGL_draw_clear()`
  - `pgraph_init_OpenGL()` (where function pointers are set)
- Keep the function pointer fields in `PGRAPHState` (they'll serve as the D3D11
  backend hook point later). Set them to `nullptr` unconditionally in `pgraph_init()`.
- Delete all GL helper functions in the same file:
  - `pgraph_update_surface_part()`, `pgraph_bind_textures()`,
    `pgraph_apply_anti_aliasing_factor()`, `pgraph_get_surface_color_write_mask()`,
    any function containing `gl*` calls
- Remove all `#include <GL/...>` from this file
- Remove `assert(pg->opengl_enabled)` guards (dead asserts)

**Test:** Build succeeds. BumpEarth, PixelShader XDK samples render identically.

### 1.2 — Remove GLSL shader translators

**Delete files:**
- `src/devices/video/nv2a_vsh.cpp` — GLSL vertex shader translator
- `src/devices/video/nv2a_psh.cpp` — GLSL pixel shader translator
- `src/devices/video/nv2a_shaders.cpp` — OpenGL shader compile/link utilities
- `src/devices/video/nv2a_shaders.h` — ShaderBinding struct (GL types)

**Update:** Remove `#include` references to deleted files from:
- `EmuNV2A_PGRAPH.cpp`
- `nv2a.cpp`
- Any CMakeLists.txt `target_sources` lists

**Test:** Build succeeds. No functional change.

### 1.3 — Remove OpenGL types from PGRAPHState

**File:** `src/devices/video/nv2a_int.h`

Remove these fields from `PGRAPHState`:
- `GloContext *gl_context`
- `GLuint gl_framebuffer`, `gl_color_buffer`, `gl_zeta_buffer`
- `GLuint gl_element_buffer`, `gl_memory_buffer`, `gl_vertex_array`
- `GLuint *gl_zpass_pixel_count_queries`
- `ShaderBinding *shader_binding`
- `TextureBinding *texture_binding[4]` (or just the GL fields within)

Remove GL fields from `VertexAttribute`:
- `GLint gl_count`, `GLenum gl_type`, `GLboolean gl_normalize`
- `GLuint gl_converted_buffer`, `gl_inline_buffer`

Remove the `#include <GL/glew.h>` and `#include "gloffscreen.h"` from `nv2a_int.h`.

**Also clean up:**
- `src/devices/video/nv2a_debug.h` — remove GL includes
- `src/devices/video/nv2a.cpp` — remove `#include "glextensions.h"`, remove
  `pgraph_init_OpenGL()` call, remove GL context creation in `pgraph_init()`

**Test:** Build succeeds. All PGRAPHState NV2A-pure fields preserved.

### 1.4 — Remove gloffscreen library and GLEW dependency

**Delete directories/files:**
- `src/common/util/gloffscreen/` (entire directory)
- `import/glew-2.0.0/` (entire directory)

**Update:**
- `CMakeLists.txt`: Remove GLEW find/link, remove gloffscreen sources, remove
  GLEW DLL copy-to-output rules
- `projects/` vcxproj files if they reference GLEW/gloffscreen

**Test:** Clean build succeeds. No GLEW DLL in output. BumpEarth renders.

### 1.5 — Remove `opengl_enabled` / `bLLE_GPU` infrastructure

**Files to update:**
- `src/devices/video/nv2a_int.h` — remove `bool opengl_enabled` from `PGRAPHState`
- `src/devices/video/nv2a.cpp` — remove `bLLE_GPU` usage, remove conditional GL init
- `src/devices/video/EmuNV2A_PFIFO.cpp` — simplify puller: remove
  `opengl_enabled` branching (always take the HLE path that dispatches methods)
- `src/core/hle/Intercept.hpp` / `.cpp` — remove `extern bool bLLE_GPU`
- `src/core/kernel/init/CxbxKrnl.cpp` — remove `bLLE_GPU` assignment
- `src/gui/WndMain.cpp` — remove LLE_GPU menu toggle
- `src/gui/resource/Cxbx.rc` — remove grayed menu item
- `src/core/hle/Patches.cpp` — remove `bLLE_GPU` skip-patches logic
- `src/core/common/imgui/ui.cpp` — remove `bLLE_GPU` overlay check

**Do NOT remove** the `LLE_GPU` flag constant itself yet — other LLE flags (APU, etc.)
use the same enum. Just ensure it's never checked.

**Test:** Build succeeds. No behavioral change. Menu simplified.

---

## Step 2: Solve the PGRAPH Race Condition  ✅ DONE

Implemented `pfifo_flush_to_pgraph()` in `EmuNV2A_PFIFO.cpp`. Called at
the top of `CxbxUpdateNativeD3DResources()` before every HLE draw. Blocks
until the DMA pusher has pushed all commands into CACHE1 and the puller has
dispatched them all to `pgraph_handle_method()`.

This is the critical blocker. The PFIFO puller processes pushbuffer commands
asynchronously. HLE EMUPATCH draw calls execute before the puller catches up,
so PGRAPH `regs[]` are stale/zero at draw time.

### 2.1 — Add PFIFO flush/synchronize primitive

**File:** `src/devices/video/EmuNV2A_PFIFO.cpp`

Implement `pfifo_flush_to_pgraph()`:
- Signal the puller thread to drain all pending CACHE1 entries
- Block until the puller has processed everything up to the current PUT pointer
- Use an event/condition variable: set by puller when CACHE1 is empty and
  GET == PUT, waited on by the caller

**Implementation approach:**
```cpp
// In NV2AState or PFIFOState:
HANDLE hPullerFlushEvent;  // auto-reset event

void pfifo_flush_to_pgraph(NV2AState *d) {
    // If puller is idle (GET == PUT and CACHE1 empty), return immediately
    // Otherwise, signal puller to wake, then WaitForSingleObject(hPullerFlushEvent)
    // Puller sets the event after draining CACHE1 when it sees GET == PUT
}
```

**Test:** Call `pfifo_flush_to_pgraph()` from a test point. Verify PGRAPH `regs[]`
contain expected combiner state after an XDK sample's first draw.

### 2.2 — Insert flush call before each HLE draw

**File:** `src/core/hle/D3D8/Rendering/HostSync.cpp`

At the top of `CxbxUpdateNativeD3DResources()` (the pre-draw sync point), add:
```cpp
extern void pfifo_flush_to_pgraph(NV2AState *d);
if (g_NV2A) pfifo_flush_to_pgraph(g_NV2A);
```

This ensures all pushbuffer commands preceding the HLE draw call have been
processed into PGRAPH state before the interpreters read it.

**Performance note:** This serializes the puller with the HLE thread per draw.
This is acceptable as a transitional measure. The puller processes methods fast
(just register writes). In the target architecture (Step 5+), the draw trigger
moves to the puller itself, eliminating this sync point.

**Test:** BumpEarth — verify PGRAPH combiner registers are non-zero at draw time.
Add a diagnostic dump if needed to confirm.

---

## Step 3: Migrate RC Interpreter to Read PGRAPH State

### 3.1 — Switch core combiner registers to PGRAPH source

**File:** `src/core/hle/D3D8/XbPixelShaderCompiler.cpp`
**Function:** `CxbxD3D11UploadRCInterpreterState()`

Replace PSDef reads with PGRAPH `regs[]` reads for these fields:

| CB Field | Current Source (PSDef) | New Source (PGRAPH) |
|----------|----------------------|---------------------|
| `PSAlphaInputs[i]` | `pPSDef->PSAlphaInputs[i]` | `pg->regs[NV_PGRAPH_COMBINEALPHAI0 + i*4]` |
| `PSAlphaOutputs[i]` | `pPSDef->PSAlphaOutputs[i]` | `pg->regs[NV_PGRAPH_COMBINEALPHAO0 + i*4]` |
| `PSRGBInputs[i]` | `pPSDef->PSRGBInputs[i]` | `pg->regs[NV_PGRAPH_COMBINECOLORI0 + i*4]` |
| `PSRGBOutputs[i]` | `pPSDef->PSRGBOutputs[i]` | `pg->regs[NV_PGRAPH_COMBINECOLORO0 + i*4]` |
| `PSConstant0[i]` | `pPSDef->PSConstant0[i]` | `pg->regs[NV_PGRAPH_COMBINEFACTOR0 + i*4]` |
| `PSConstant1[i]` | `pPSDef->PSConstant1[i]` | `pg->regs[NV_PGRAPH_COMBINEFACTOR1 + i*4]` |
| `PSFinalCombinerInputsABCD` | `pPSDef->PSFinalCombiner...` | `pg->regs[NV_PGRAPH_COMBINESPECFOG0]` |
| `PSFinalCombinerInputsEFG` | `pPSDef->PSFinalCombiner...` | `pg->regs[NV_PGRAPH_COMBINESPECFOG1]` |
| `PSCombinerCount` | `pPSDef->PSCombinerCount` | `pg->regs[NV_PGRAPH_COMBINECTL]` |
| `PSTextureModes` | render state / PSDef | `pg->regs[NV_PGRAPH_SHADERPROG]` |
| `PSInputTexture` | `pPSDef->PSInputTexture` | `pg->regs[NV_PGRAPH_SHADERCTL]` |
| `PSCompareMode` | `pPSDef->PSCompareMode` | `pg->regs[NV_PGRAPH_SHADERCLIPMODE]` |
| `PSDotMapping` | `pPSDef->PSDotMapping` | PGRAPH (verify NV097 method populates regs) |

Keep PSDef as fallback: `if (g_NV2A == nullptr) { /* read from PSDef */ }`.
This ensures compatibility if NV2A init is delayed.

**Test:** BumpEarth, BumpLens, DotProduct3 — compare rendered output before/after.
The combiner state values from PGRAPH and PSDef should be identical (both come from
the same Xbox D3D calls, just via different paths).

### 3.2 — Migrate remaining RC state fields

Fields that don't have direct PGRAPH register equivalents:
- `ColorSign[4]` — derived from texture format; may need to stay HLE-sourced
- `FogColor`, `FogInfo`, `FogEnable` — from render state registers in PGRAPH
- `AlphaTest` — from `NV_PGRAPH_CONTROL_0`
- `BEM[4]`, `LUM[4]` — bump environment map; from PGRAPH bump env matrix
- `ColorKeyOp/Color[4]` — from PGRAPH texture state

Migrate each field individually. Test after each sub-group.

**Test:** PixelShader XDK sample (exercises many combiner configurations).

### 3.3 — Remove PSDef dependency from RC upload path

Once all fields read from PGRAPH:
- Remove the `X_D3DPIXELSHADERDEF` read from `CxbxD3D11UploadRCInterpreterState()`
- Remove `GetPixelShaderRenderStatePointer()` dependency
- Keep `g_NV2A == nullptr` guard as a no-op (return early, skip RC upload)

**Test:** Full XDK sample suite.

---

## Step 4: Migrate VS Interpreter to Read PGRAPH State

### 4.1 — Switch VS microcode source to PGRAPH

**File:** `src/core/hle/D3D8/XbVertexShader.cpp`
**Function:** `CxbxD3D11UploadVSInterpreterState()`

Replace HLE shader slot read with PGRAPH read:
```cpp
// Current: reads from HLE shader slot cache
const uint32_t* pTokens = (const uint32_t*)pXboxMicrocode;

// New: reads from PGRAPHState.program_data[]
PGRAPHState *pg = &g_NV2A->pgraph;
for (int i = 0; i < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH; i++) {
    cb.Instructions[i] = uint4(pg->program_data[i][0],
                               pg->program_data[i][1],
                               pg->program_data[i][2],
                               pg->program_data[i][3]);
}
```

### 4.2 — Switch VS constants source to PGRAPH

**Function:** `CxbxUpdateHostVertexShaderConstants()`

Replace HLE constant buffer read with PGRAPH read:
```cpp
// Current: reads from g_Xbox_VertexShaderConstantMode slots
// New: reads from PGRAPHState.vsh_constants[192][4]
```

The vsh_constants array stores raw uint32_t[4] per constant. Convert to float4
with `reinterpret_cast<float*>`.

### 4.3 — Switch VS mode detection to PGRAPH

**Function:** `CxbxUpdateHostVertexShader()`

Detect fixed-function vs. programmable from PGRAPH:
```cpp
// NV_PGRAPH_CSV0_D (offset 0x0FB4) bit 0: vertex program enable
bool bProgrammable = (pg->regs[NV_PGRAPH_CSV0_D / 4] & 1) != 0;
```

Replace the current check that reads from `g_Xbox_VertexShader_Handle`.

**Test:** Dolphin (programmable VS), BumpEarth (FF VS), PixelShader (passthrough).

---

## Step 5: Migrate Vertex Attribute Fetch to PGRAPH State

### 5.1 — Read vertex array descriptors from PGRAPH

**File:** `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_IABypass.cpp`

Currently reads `g_Xbox_SetStreamSource[]` for per-stream base/stride/offset.
Switch to reading `PGRAPHState.vertex_attributes[16]`:

```cpp
PGRAPHState *pg = &g_NV2A->pgraph;
for (int i = 0; i < 16; i++) {
    const VertexAttribute &attr = pg->vertex_attributes[i];
    layout.attribs[i].offset   = attr.offset;
    layout.attribs[i].stride   = attr.stride;
    layout.attribs[i].format   = MapNV2AFormatToVtxFmt(attr.format, attr.size);
    layout.attribs[i].base     = ResolveDMAAddress(attr.dma_select,
                                                    pg->dma_vertex_a,
                                                    pg->dma_vertex_b);
}
```

The DMA address resolution converts the NV2A DMA context + offset into a physical
Xbox memory address that indexes into the 64MB SRV.

### 5.2 — Read inline/sticky vertex attributes from PGRAPH

`VertexAttribute.inline_value[4]` provides the NV2A "sticky" per-attribute defaults.
Replace `HLE_get_NV2A_vertex_attribute_value_pointer()` with direct PGRAPH read.

**Test:** Fur, SphereMap (multiple vertex streams).

---

## Step 6: Migrate Surface and Pipeline State to PGRAPH

### 6.1 — Surface/render target state

**File:** `src/core/hle/D3D8/Rendering/HostSync.cpp`

Replace `g_pXbox_RenderTarget` / `g_pXbox_DepthStencil` reads with PGRAPH:
- Color offset: `pg->surface_color.offset` (populated by `NV097_SET_SURFACE_COLOR_OFFSET`)
- Zeta offset: `pg->surface_zeta.offset` (populated by `NV097_SET_SURFACE_ZETA_OFFSET`)
- Format: `pg->surface_shape.color_format`, `pg->surface_shape.zeta_format`
- Pitch: from `NV097_SET_SURFACE_PITCH` method handler
- Dimensions: `pg->surface_shape.clip_x` / `clip_y` / `clip_width` / `clip_height`

### 6.2 — Pipeline state (blend, depth, stencil, rasterizer)

**File:** `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_State.cpp`

Replace `XboxRenderStates.Apply()` reads with PGRAPH register reads:
- Blend: `pg->regs[NV_PGRAPH_BLEND]` — equation, factors, enable
- Depth: `pg->regs[NV_PGRAPH_CONTROL_0]` — depth test, write, func
- Stencil: `pg->regs[NV_PGRAPH_CONTROL_1]` — stencil ops, ref, mask
- Rasterizer: `pg->regs[NV_PGRAPH_SETUPRASTER]` — cull mode, fill mode
- Fog: `pg->regs[NV_PGRAPH_FOGCOLOR]`, `NV_PGRAPH_FOGPARAM0/1`
- Alpha test: within `NV_PGRAPH_CONTROL_0`

These map to `D3D11_BLEND_DESC`, `D3D11_DEPTH_STENCIL_DESC`, `D3D11_RASTERIZER_DESC`
objects which are created/cached on the CPU (this is the one part that can't be
GPU-driven in DX11, as noted in the architecture doc).

### 6.3 — Viewport and scissor

Replace `g_Xbox_Viewport` reads with PGRAPH viewport registers.

**Test:** ShadowBuffer (depth/stencil), Glass (blend), ProjectedTexture (viewport).

---

## Step 7: Migrate Texture State to PGRAPH

### 7.1 — Texture binding from PGRAPH

**File:** `src/core/hle/D3D8/Rendering/HostSync.cpp`
**Function:** `CxbxUpdateHostTextures()`

Currently reads `g_pXbox_SetTexture[stage]` (HLE global) and falls back to the
device's internal `m_Textures[]` array.

Switch to reading texture address from PGRAPH:
- Texture offset: `pg->regs[NV_PGRAPH_TEXOFFSET0 + stage * 0x40]`
- Texture format: `pg->regs[NV_PGRAPH_TEXFMT0 + stage * 0x40]`
- Texture control: `pg->regs[NV_PGRAPH_TEXCTL0_0 + stage * 0x40]`
- Image rect: `pg->regs[NV_PGRAPH_TEXIMAGERECT0 + stage * 0x40]`

The texture address resolves via DMA context (`dma_a`/`dma_b`) to a physical
Xbox memory address. The existing deswizzle pipeline and `Texture2D` creation
remain — only the source of the texture descriptor changes.

### 7.2 — Remove texture fallback from device internals

Once PGRAPH provides all texture state, remove:
- The `m_Textures[]` device memory scan fallback in `CxbxUpdateHostTextures()`
- The `D3D_g_pDevice` + `m_Textures_Offset` extraction from `D3DDevice_SetTexture`
  machine code

**Test:** PixelShader (multi-texture), BumpEarth (bump map textures), CubeMap.

---

## Step 8: Move Draw Triggering to PGRAPH (Puller-Driven Draws)

This is the architectural pivot. Instead of HLE EMUPATCH draw functions triggering
D3D11 draws, the PFIFO puller triggers draws when it processes `NV097_SET_BEGIN_END(0)`.

### 8.1 — Register D3D11 draw backend as PGRAPH draw functions

**File:** `src/devices/video/EmuNV2A_PGRAPH.cpp`

Set the draw function pointers to D3D11 backend functions:
```cpp
void pgraph_init_D3D11(PGRAPHState *pg) {
    pg->pgraph_draw_arrays         = D3D11_draw_arrays;
    pg->pgraph_draw_inline_buffer  = D3D11_draw_inline_buffer;
    pg->pgraph_draw_inline_array   = D3D11_draw_inline_array;
    pg->pgraph_draw_inline_elements = D3D11_draw_inline_elements;
    pg->pgraph_draw_state_update   = D3D11_draw_state_update;
    pg->pgraph_draw_clear          = D3D11_draw_clear;
}
```

Each `D3D11_draw_*` function:
1. Reads all state from `PGRAPHState` (combiner, VS, vertex arrays, textures, pipeline)
2. Uploads constant buffers (RC interpreter CB, VS interpreter CB, layout CB)
3. Binds textures, render targets, pipeline state objects
4. Calls `g_pD3DDeviceContext->Draw(hostVertexCount, 0)`

This reuses all existing D3D11 infrastructure — the change is just the trigger
point and state source.

### 8.2 — Handle threading: puller thread → D3D11 device context

The puller thread runs on a separate thread from the D3D11 device creation thread.
D3D11 device contexts are single-threaded. Options:

**Option A (recommended):** Create a deferred context on the puller thread.
Record command lists, execute on the immediate context from the render thread.

**Option B:** Use `ID3D11Multithread` to enable MT access. Simpler but may have
driver-specific performance implications.

**Option C:** Queue draw commands from the puller to the render thread via a
lock-free ring buffer. The render thread drains the queue each frame.

### 8.3 — Remove the `pfifo_flush_to_pgraph()` sync call

Once draws are puller-driven, the flush from Step 2.2 is no longer needed.
Remove it from `CxbxUpdateNativeD3DResources()`.

**Test:** BumpEarth — verify draws come from puller, not from HLE patches.
Frame timing may change (draws happen when puller processes them, not when
HLE patch runs). Watch for tearing / out-of-order issues.

---

## Step 9: Remove HLE EMUPATCH Draw Functions

With puller-driven draws working, EMUPATCH draw functions are redundant.
Remove them in dependency order.

### 9.1 — Remove state-setting patches

These patches mirror Xbox state into HLE structures. Once all state comes from
PGRAPH, the original Xbox D3D code handles this via pushbuffer:

**Batch 1 — Pure state:**
- `D3DDevice_SetRenderState_Simple` and all `SetRenderState_*` variants
- `D3DDevice_SetVertexShaderConstant` / `NotInline` / `1` / `4`
- `D3DDevice_SetPixelShader`
- `D3DDevice_SetTexture`
- `D3DDevice_SetTransform` / `MultiplyTransform`
- `D3DDevice_SetStreamSource` / `SetIndices`
- `D3DDevice_SetVertexShader` / `SelectVertexShader`
- `D3DDevice_SetViewport` / `SetScissors`

**Batch 2 — Resource management:**
- `D3DDevice_SetPalette`
- `D3DDevice_SetVertexData*`
- `D3DDevice_LoadVertexShader` / `DeleteVertexShader`

### 9.2 — Remove draw call patches

- `D3DDevice_DrawVertices` / `DrawVerticesUP`
- `D3DDevice_DrawIndexedVertices` / `DrawIndexedVerticesUP`
- `D3DDevice_Begin` / `End`

### 9.3 — Remove presentation patches

- `D3DDevice_Present` / `Swap`
- Replace with PCRTC/PVIDEO handling: read framebuffer from VRAM at vblank

### 9.4 — Remove sync patches (requires NV2A fence emulation)

- `D3DDevice_BlockOnFence` / `InsertFence` / `IsFencePending`
- `D3DDevice_BeginVisibilityTest` / `EndVisibilityTest` / `GetVisibilityTestResult`
- Implement `NV097_SET_SEMAPHORE_OFFSET` + `NV097_BACK_END_WRITE_SEMAPHORE_RELEASE`
  in pgraph_handle_method

**Test each batch independently.** After each batch, run the full XDK sample suite.

---

## Step 10: Clean Up HLE State Infrastructure

### 10.1 — Remove HLE state globals

Once no code reads from them:
- `XboxRenderStates` / `XboxTextureStates` — render/texture state mirrors
- `g_pXbox_SetTexture[]` — per-stage texture pointers
- `g_Xbox_SetStreamSource[]` — per-stream VB bindings
- `g_pXbox_RenderTarget` / `g_pXbox_DepthStencil`
- `g_Xbox_VertexShader_Handle` / `g_Xbox_VertexShader_FunctionSlots_StartAddress`

### 10.2 — Remove EMUPATCH infrastructure for removed patches

- Clean up `EmuPatches_*.cpp` files
- Remove symbol scan entries for deleted patches
- Remove trampoline slots

### 10.3 — Remove remaining GL fields from VertexAttribute

**File:** `src/devices/video/nv2a_int.h`

Remove `needs_conversion`, `converted_buffer`, `converted_elements`,
`converted_size`, `converted_count` — these were for the OpenGL path's
format conversion. The D3D11 IA bypass handles format decode in shader.

**Test:** Full suite.

---

## Step 11: Prepare for Vulkan Backend (No Vulkan Code)

### 11.1 — Abstract the rendering backend interface

Create a `RenderBackend` interface (pure virtual class or function table):

```cpp
struct RenderBackend {
    void (*DrawArrays)(PGRAPHState *pg);
    void (*DrawInlineBuffer)(PGRAPHState *pg);
    void (*DrawInlineArray)(PGRAPHState *pg);
    void (*DrawInlineElements)(PGRAPHState *pg);
    void (*UpdateState)(PGRAPHState *pg);
    void (*Clear)(PGRAPHState *pg);
    void (*Present)(PGRAPHState *pg);
    void (*Init)(HWND hwnd);
    void (*Shutdown)();
};
```

This maps directly to the existing `pgraph_draw_*` function pointers but is
more formally structured. The D3D11 backend implements this interface.

### 11.2 — Isolate D3D11 code into backend module

Move all D3D11-specific code into `src/core/hle/D3D8/Rendering/Backend/`:
- Ensure no D3D11 types leak outside the Backend directory
- All PGRAPH → backend communication goes through the `RenderBackend` interface
- The HLSL shaders are backend-specific (a Vulkan backend would use SPIR-V or
  cross-compiled HLSL)

### 11.3 — Make PGRAPH state the single source of truth

Verify that all rendering decisions read from `PGRAPHState` only:
- No references to `XboxRenderStates`, `g_pXbox_*`, or `pPSDef` in the render path
- The `RCInterpreterCBLayout` is populated entirely from `pg->regs[]`
- The `VSInterpreterCBLayout` is populated entirely from `pg->program_data[]`
  and `pg->vsh_constants[]`
- The `IABypassLayoutCB` is populated entirely from `pg->vertex_attributes[]`

### 11.4 — Add Vulkan SDK to CMakeLists.txt (optional prep)

```cmake
find_package(Vulkan QUIET)
if(Vulkan_FOUND)
    message(STATUS "Vulkan SDK found: ${Vulkan_LIBRARY}")
    # Don't link yet — just verify availability
endif()
```

**Test:** Full XDK suite. Build with and without Vulkan SDK present.

---

## Dependency Graph

```
Step 1  (Remove OpenGL LLE)         ── no dependencies, safe first
   │
Step 2  (PGRAPH race condition fix) ── independent of Step 1
   │
   ├──► Step 3  (RC interpreter → PGRAPH)
   │       │
   ├──► Step 4  (VS interpreter → PGRAPH)     ── parallel with Step 3
   │       │
   ├──► Step 5  (Vertex fetch → PGRAPH)        ── parallel with 3,4
   │       │
   └──► Step 6  (Surface/pipeline → PGRAPH)    ── parallel with 3,4,5
           │
           ▼
        Step 7  (Texture state → PGRAPH)       ── after 6 (surface state)
           │
           ▼
        Step 8  (Puller-driven draws)          ── after 3,4,5,6,7
           │
           ▼
        Step 9  (Remove EMUPATCH draws)        ── after 8
           │
           ▼
        Step 10 (Clean up HLE state)           ── after 9
           │
           ▼
        Step 11 (Backend abstraction)          ── after 10
```

Steps 3, 4, 5, 6 can be worked in parallel once Step 2 is done.
Step 1 can be done at any time (no functional dependency).

---

## Validation Strategy

After every step:
1. **Build test** — clean compile, no warnings related to changed files
2. **XDK sample smoke test** — run these samples and verify visual output:
   - BumpEarth (FF VS + bump mapping + RC combiners)
   - BumpLens (dot product texture modes)
   - PixelShader (multi-texture, various combiner setups)
   - Fur (programmable VS + multi-pass)
   - Dolphin (programmable VS + water effects)
   - DotProduct3 (dot mapping modes)
   - CubeMap (cube texture + reflection)
   - ShadowBuffer (depth/stencil heavy)
   - Glass (alpha blending)
3. **Regression diff** — capture a reference frame (screenshot) from each sample
   before beginning work. Compare after each step. Pixel differences indicate
   a state-source mismatch.

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| PGRAPH flush stalls hurt frame rate | Profile. Puller processes fast (register writes only). If too slow, batch flushes per frame not per draw. |
| PGRAPH state differs from HLE state | Add diagnostic mode: compare PGRAPH vs PSDef values at draw time. Log mismatches. |
| Removing EMUPATCH breaks LTCG titles | Remove patches incrementally. Keep trampoline stubs that just call Xbox code. |
| D3D11 threading with puller | Use deferred context (Option A). Well-tested D3D11 pattern. |
| Texture state from PGRAPH incomplete | Keep HLE texture fallback as safety net until Step 7.2. |
| Some games rely on EMUPATCH side effects | Per-patch analysis. Whitelist patches that must remain (e.g., `D3D_CreateDevice`). |
