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
| RC interpreter reads from | **PGRAPH `regs[]` raw SRV** (Steps 3.1–3.4 done, CBLayout eliminated) |
| VS interpreter reads from | **PGRAPH `program_data[]` via XFPR SRV** (Steps 4.1–4.3 done, CBLayout eliminated) |
| VS constants reads from | **PGRAPH `vsh_constants[]` (XFCTX RAM mirror)** (already PGRAPH-sourced) |
| Vertex fetch reads from | **PGRAPH `vertex_attributes[]`** (Steps 5.1–5.2 done, inline_value direct) |
| Surface/RT state from | **PGRAPH `surface_color/zeta.offset`** via side-map (Step 6.1 done) |
| Pipeline state from | **PGRAPH `regs[]`** — blend/depth/stencil/rasterizer/viewport (Steps 6.2–6.3 done) |
| Texture state from | **PGRAPH `TEXOFFSET` regs** via side-map (Steps 7.1–7.2 done, m_Textures removed) |
| Draw trigger | **Puller-driven** `pgraph_draw(d)` from `SET_BEGIN_END(END)` (Step 8.1 done) |
| Draw batching | **Cross-bracket squash** — consecutive BEGIN/DRAW_ARRAYS/END merged into single batch |
| OpenGL LLE path | **Removed** (Step 1 complete) |
| OpenGL types in structs | **Removed** from `PGRAPHState`, `VertexAttribute` |
| GLSL shader translators | **Removed** |
| gloffscreen library | **Removed** |
| GLEW import | **Removed** |
| PFIFO→PGRAPH flush | **Active** — skipped when `g_bInPullerContext` (Step 8.3) |

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

### 3.1 — Switch core combiner registers to PGRAPH source ✅ DONE

**File:** `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_PixelShader.cpp` *(moved from XbPixelShaderCompiler.cpp)*
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

### 3.2 — Migrate remaining RC state fields  ✅ DONE

Fields that don't have direct PGRAPH register equivalents:
- `ColorSign[4]` — derived from texture format; may need to stay HLE-sourced
- `FogColor`, `FogInfo`, `FogEnable` — from render state registers in PGRAPH
- `AlphaTest` — from `NV_PGRAPH_CONTROL_0`
- `BEM[4]`, `LUM[4]` — bump environment map; from PGRAPH bump env matrix
- `ColorKeyOp/Color[4]` — from PGRAPH texture state

PSFinalCombinerConstant, FogColor, AlphaTest, BEM, LUM migrated (commit 52b8a152).
FogInfo/FogEnable migrated from PGRAPH CONTROL_3 (commit d9311588).

Migrate each field individually. Test after each sub-group.

**Test:** PixelShader XDK sample (exercises many combiner configurations).

### 3.3 — Remove PSDef dependency from RC upload path  ✅ DONE

PSDef soft dependency only remains for texModeAdjust flag.
`CxbxD3D11UploadRCInterpreterState()` reads all combiner/fog/alpha state from PGRAPH.

**Test:** Full XDK sample suite.

### 3.4 — Replace RCInterpreterCBLayout with raw PGRAPH regs[] buffer

**Architectural change:** Eliminate the `RCInterpreterCBLayout` struct entirely.
Instead of C++ code extracting individual fields from `pg->regs[]` into a custom
struct, upload the raw `pg->regs[]` array (8 KB, 2048 × uint32) as a structured
buffer or constant buffer. The HLSL RC interpreter shader reads registers directly
using their `NV_PGRAPH_*` byte-offset indices.

**Rationale:** The current CBLayout is a hand-maintained C++ ↔ HLSL struct with
84 fields, 1344 bytes, and complex packing (uint registers padded to float4, ABGR
unpacking, etc.). By passing `regs[]` raw, the shader accesses each register at its
hardware-defined offset with its natural type (uint or float). This:
- Eliminates the C++ upload function that extracts and packs fields
- Makes register access self-documenting (shader reads `regs[RI(NV_PGRAPH_COMBINECTL)]`)
- Avoids alignment/padding bugs between C++ and HLSL
- Trivially extends to any new register without touching C++ code

**HLSL side:**
```hlsl
// Replace cbuffer RCInterpreterCBLayout : register(b0) with:
StructuredBuffer<uint> g_PGRegs : register(t4);  // pg->regs[] as uint SRV

// Helper to read a register by its byte offset:
uint  PG_UINT(uint byteOff)  { return g_PGRegs[byteOff >> 2]; }
float PG_FLOAT(uint byteOff) { return asfloat(g_PGRegs[byteOff >> 2]); }

// Example usage:
uint combinerCount = PG_UINT(NV_PGRAPH_COMBINECTL);       // was: PSCombinerCount
uint alphaInput0   = PG_UINT(NV_PGRAPH_COMBINEALPHAI0);   // was: PSAlphaInputs[0]
```

**C++ side:**
- Create a D3D11 structured buffer (8 KB) backed by `pg->regs[]`
- Upload with `Map/memcpy/Unmap` before each draw (or use dirty tracking)
- Bind as SRV to PS slot t4 (avoids conflict with existing t0-t3 texture slots)
- For fields not in `regs[]` (ColorSign, TexFmtFixup, etc.), use a small
  auxiliary cbuffer with only the software-computed fields

**Software-computed fields** that have no PGRAPH register:
- `ColorSign[4]`, `TexFmtFixup`, `ColorKeyOp/Color[4]` — keep in aux cbuffer
- `FogInfo`, `FogEnable` — migrate to PGRAPH regs first (Step 3.2)
- `AlphaKill`, `FrontFaceInfo` — keep in aux cbuffer (runtime-computed)

**For ABGR color registers** (`PSConstant0/1`, `PSFinalCombinerConstant`,
`FogColor`): These are stored as packed ABGR uint32 in `regs[]`. The shader
unpacks them inline: `float4 c = UnpackABGR(PG_UINT(NV_PGRAPH_COMBINEFACTOR0 + i*4))`.

**Migration approach:** Single big-bang replacement.
1. Create the regs[] `StructuredBuffer<uint>` SRV (8 KB) and bind to PS t4
2. Rewrite the HLSL RC interpreter to read all register-sourced fields from
   `g_PGRegs[]` using `PG_UINT()` / `PG_FLOAT()` helpers
3. Replace `RCInterpreterCBLayout` with a small `PSAuxCBLayout` containing
   only the software-computed fields (ColorSign, TexFmtFixup, AlphaKill, etc.)
4. Remove `CxbxD3D11UploadRCInterpreterState()` field-by-field extraction;
   replace with a single `memcpy` of `pg->regs[]` into the SRV
5. Remove `RCInterpreterCBLayout` struct and `CxbxRegisterCombinerInterpreterState.hlsli`

All done in one commit. No intermediate hybrid state.

✅ DONE — Committed as 09a1ef0f. RC interpreter reads raw `regs[]` SRV at t12.
PSAuxCBLayout holds software-computed fields only.

**Test:** BumpEarth, PixelShader, DotProduct3, BumpLens, Dolphin.

---

## Step 4: Migrate VS Interpreter to Read PGRAPH State

### 4.1 — Switch VS microcode source to PGRAPH  ✅ DONE

**File:** `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_VertexShader.cpp` *(moved from XbVertexShader.cpp)*
**Function:** `CxbxD3D11UploadVSInterpreterState()`

Reads `pg->xf.xfpr[startSlot+i][0..3]` using CHEOPS_PROGRAM_START from
`NV_PGRAPH_CSV0_C`. Committed as 55c680c3.

### 4.2 — Switch VS constants source to PGRAPH  ✅ ALREADY DONE

`CxbxUpdateHostVertexShaderConstants()` already reads `pg->xf.xfctx[]`
with dirty tracking. No migration needed.

### 4.3 — Replace VSInterpreterCBLayout with raw PGRAPH buffers  ✅ DONE

**Architectural change:** Eliminate the `VSInterpreterCBLayout` struct entirely.
Instead of C++ code packing `program_data[]` into `Instructions[136]` + `InstructionCount`,
upload `pg->xf.xfpr[]` directly as a structured buffer.

**Hardware background (XFPR RAM):**
On real NV2A, vertex shader microcode lives in the **XFPR (Transform Program RAM)**
— on-chip XF SRAM with 136 slots of 92-bit instructions in 128-bit containers.
This is not VRAM or MMIO-visible; the CPU reaches it through the
`NV097_SET_TRANSFORM_PROGRAM` method range (32 DWORDs per batch), with
`NV_PGRAPH_CHEOPS_OFFSET.PROG_LD_PTR` as the auto-incrementing write pointer.
`NV097_SET_TRANSFORM_PROGRAM_LOAD` resets the write pointer.
The RDI (Register Direct Interface) is used separately for context save/restore.
`pg->xf.xfpr[136][4]` in `PGRAPHState` is the software mirror of XFPR.

**D3D11 implementation:**
```hlsl
// XFPR SRV (Transform Program RAM — pg->xf.xfpr[] mirror):
StructuredBuffer<uint4> g_XFPR : register(t5);

// Program start from regs[]:
// uint startSlot = PG_UINT(NV_PGRAPH_CSV0_C) >> CHEOPS shift;
// Loop until FLD_FINAL — no InstructionCount needed.
```

**C++ side:**
- Create a D3D11 structured buffer for `program_data[136][4]` (2176 bytes)
- Upload with `Map/memcpy/Unmap` (the data is already contiguous in PGRAPH)
- Bind as SRV to VS slot t5
- The program start register is in `regs[]` — accessible via the shared
  `g_PGRegs` SRV at t12 (same buffer bound to PS and VS stages)

**The VS opcode/field constants** (`FLD_ILU`, `FLD_MAC`, etc.) defined in
`CxbxVertexShaderInterpreterState.hlsli` are compile-time constants, not
runtime state — they remain as `#define`s in the HLSL include.

**For `vsh_constants[192][4]` (XFCTX RAM):** Already uploaded separately as
the VS constants cbuffer (b0). On real hardware this is XFCTX (Transform
Context RAM), also on-chip XF SRAM behind the RDI interface. No change needed.

Committed as b8b3b99c. Tested: BumpEarth, Dolphin, Fur.

---

## Step 5: Migrate Vertex Attribute Fetch to PGRAPH State

### 5.1 — Read vertex array descriptors from PGRAPH  ✅ DONE

**File:** `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_VertexFetch.cpp`

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

### 5.2 — Read inline/sticky vertex attributes from PGRAPH  ✅ DONE

`VertexAttribute.inline_value[4]` provides the NV2A "sticky" per-attribute defaults.
Removed `NV2A_get_vertex_attribute_value_pointer()` wrapper; `UploadVertexDefaults()`
now reads `pg->vertex_attributes[i].inline_value` directly.

**Test:** Fur, SphereMap (multiple vertex streams).

---

## Step 6: Migrate Surface and Pipeline State to PGRAPH

### 6.1 — Surface/render target state  ✅ DONE

**File:** `src/core/hle/D3D8/Rendering/HostSync.cpp`

Side-map (VRAM offset → Xbox surface*) populated by EMUPATCH. Committed as 20e012c0.

Replace `g_pXbox_RenderTarget` / `g_pXbox_DepthStencil` reads with PGRAPH:
- Color offset: `pg->surface_color.offset` (populated by `NV097_SET_SURFACE_COLOR_OFFSET`)
- Zeta offset: `pg->surface_zeta.offset` (populated by `NV097_SET_SURFACE_ZETA_OFFSET`)
- Format: `pg->surface_shape.color_format`, `pg->surface_shape.zeta_format`
- Pitch: from `NV097_SET_SURFACE_PITCH` method handler
- Dimensions: `pg->surface_shape.clip_x` / `clip_y` / `clip_width` / `clip_height`

### 6.2 — Pipeline state (blend, depth, stencil, rasterizer)  ✅ DONE

**File:** `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_State.cpp`

Committed as 0a2566a2. `CxbxD3D11UpdatePipelineStateFromPGRAPH()` runs after
`XboxRenderStates.Apply()` to override with authoritative PGRAPH values.

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

### 6.3 — Viewport and scissor  ✅ DONE

Committed as fa49e373. `CxbxD3D11UpdateViewportFromPGRAPH()` derives Xbox viewport
from VPSCL/VPOFF in `vsh_constants[0x3a]/[0x3b]`, depth range from ZCLIPMIN/ZCLIPMAX.

**Test:** ShadowBuffer (depth/stencil), Glass (blend), ProjectedTexture (viewport).

---

## Step 7: Migrate Texture State to PGRAPH

### 7.1 — Texture binding from PGRAPH  ✅ DONE

**File:** `src/core/hle/D3D8/Rendering/HostSync.cpp`
**Function:** `CxbxUpdateHostTextures()`

Committed as c09c2bd7. Side-map (VRAM offset → Xbox texture*) populated by
SetTexture/SwitchTexture/LTCG patches. PGRAPH TEXOFFSET is authoritative source.

Switch to reading texture address from PGRAPH:
- Texture offset: `pg->regs[NV_PGRAPH_TEXOFFSET0 + stage * 0x40]`
- Texture format: `pg->regs[NV_PGRAPH_TEXFMT0 + stage * 0x40]`
- Texture control: `pg->regs[NV_PGRAPH_TEXCTL0_0 + stage * 0x40]`
- Image rect: `pg->regs[NV_PGRAPH_TEXIMAGERECT0 + stage * 0x40]`

The texture address resolves via DMA context (`dma_a`/`dma_b`) to a physical
Xbox memory address. The existing deswizzle pipeline and `Texture2D` creation
remain — only the source of the texture descriptor changes.

### 7.2 — Remove texture fallback from device internals  ✅ DONE

Committed as 16c310e2. Removed m_Textures machine-code scanning and per-frame
device memory read fallback. PGRAPH TEXOFFSET is sole texture source. (-56 lines)

**Test:** PixelShader (multi-texture), BumpEarth (bump map textures), CubeMap.

---

## Step 8: Move Draw Triggering to PGRAPH (Puller-Driven Draws)

This is the architectural pivot. Instead of HLE EMUPATCH draw functions triggering
D3D11 draws, the PFIFO puller triggers draws when it processes `NV097_SET_BEGIN_END(0)`.

### 8.1 — Register D3D11 draw backend as PGRAPH draw functions  ✅ DONE

Committed as a816390d. `D3D11_init_pgraph_plugins()` in `XbPushBuffer.cpp` sets
global function pointers. PGRAPH calls `pgraph_draw(d)` from `SET_BEGIN_END(END)`
for all draw types; the backend internally dispatches by checking which buffer
has data (draw_arrays, inline_buffer, inline_array, inline_elements).

**File:** `src/core/hle/D3D8/XbPushBuffer.cpp`

Global draw function pointers (declared in `nv2a_pgraph_backend.h`):
```cpp
void D3D11_init_pgraph_plugins() {
    g_pgraph_backend.draw              = D3D11_draw;
    g_pgraph_backend.draw_state_update = D3D11_draw_state_update;
    g_pgraph_backend.draw_clear        = D3D11_draw_clear;
    g_pgraph_backend.draw_patch        = D3D11_draw_patch;
    g_pgraph_backend.flip_stall        = D3D11_flip_stall;
    g_pgraph_backend.zpass_begin       = D3D11_zpass_begin;
    g_pgraph_backend.zpass_end         = D3D11_zpass_end;
    g_pgraph_backend.zpass_collect     = D3D11_zpass_collect;
    g_pgraph_backend.launch_transform_program = D3D11_launch_transform_program;
}
```

`D3D11_draw()` dispatches internally:
1. If `draw_arrays_length > 0` → `D3D11_draw_arrays()` (multi-batch Draw calls)
2. If `inline_elements_length > 0` → `D3D11_draw_inline_elements()` (DrawIndexed)
3. If `inline_array_length > 0` → `D3D11_draw_inline_array()` (UP-style draw)
4. If `inline_buffer_length > 0` → `D3D11_draw_inline_buffer()` (Begin/End)

Each path reads state from `PGRAPHState`, uploads constant buffers, binds
textures/render targets/pipeline state, and calls D3D11 Draw/DrawIndexed.

### 8.2 — Handle threading: puller thread → D3D11 device context  ✅ DONE

Resolved via implicit D3D11 runtime MT protection: device created without
`D3D11_CREATE_DEVICE_SINGLETHREADED` flag, so the runtime provides internal
thread safety for the immediate context. Puller thread calls D3D11 directly.

Future optimization: deferred context (Option A) or explicit `ID3D11Multithread`
(Option B) if profiling shows the runtime lock is a bottleneck.

### 8.3 — Puller context flag (partial)  ✅ DONE

Committed as a816390d. `thread_local g_bInPullerContext` prevents `pfifo_flush`
deadlock when `CxbxUpdateNativeD3DResources` is called from puller thread.
Flush is still needed for HLE-triggered draws; will be fully removed when
all draw paths are puller-driven.

**Test:** BumpEarth — verify draws come from puller, not from HLE patches.
Frame timing may change (draws happen when puller processes them, not when
HLE patch runs). Watch for tearing / out-of-order issues.

---

## Step 9: Remove HLE EMUPATCH Draw Functions  ✅ DONE

All D3D EMUPATCH entries disabled in `Patches.cpp`. Xbox native D3D code runs
unpatched; all state/draw/present/sync operations flow through the push buffer
into PFIFO → PGRAPH → D3D11 backend.

### 9.1 — Remove state-setting patches  ✅ DONE

All disabled. Xbox code writes D3D__RenderState[] then calls SetRenderState_Simple
which pushes NV2A methods through the push buffer. PGRAPH regs[] populated natively.

### 9.2 — Remove draw call patches  ✅ DONE

All disabled. Xbox DrawVertices/DrawIndexedVertices/Begin/End push NV097 methods
into the push buffer; PGRAPH `SET_BEGIN_END(END)` triggers `pgraph_draw()`.

### 9.3 — Remove presentation patches  ✅ DONE

Disabled. Native Swap pushes `NV097_FLIP_INCREMENT_WRITE` + `NV097_FLIP_STALL`;
PGRAPH FLIP_STALL handler calls `pgraph_flip_stall` → `D3D11_flip_stall` which
blits the PGRAPH backbuffer to the host swap chain.

### 9.4 — Remove sync patches  ✅ DONE

Disabled. Fences use NV2A reference counter natively. Visibility tests handled
via `pgraph_zpass_begin/end/collect` → `D3D11_zpass_*` (D3D11 occlusion queries).

**Remaining cleanup:** Dead implementation code in `Direct3D9.cpp.unused-patches`
and commented-out PATCH_ENTRY lines can be deleted (Step 10.2).

---

## Step 10: Clean Up HLE State Infrastructure

### 10.1 — Remove HLE state globals (partial)  ✅ DONE

Removed:
- `XboxRenderStates` / `XboxTextureStates` — render/texture state mirrors (deleted classes)
- `g_Xbox_SetStreamSource[]` — already dead (no references)
- `g_pXbox_RenderTarget` / `g_pXbox_DepthStencil` — already dead (no references)
- `g_Xbox_VertexShader_Handle` / `g_Xbox_VertexShader_FunctionSlots_StartAddress` — already dead

Still active (cannot remove yet):
- `g_pXbox_SetTexture[]` — used as texture side-map fallback in HostSync.cpp
  and Backend_D3D11_PixelShader.cpp

### 10.2 — Remove EMUPATCH infrastructure for removed patches  ✅ DONE

Deleted all `EmuPatches_*.cpp`, `EmuPatches_*.h`, `EmuPatches_Unused.h`, and
`Direct3D9.cpp.unused-patches`. Removed entries from CMakeLists.txt and
`#include` directives from RenderGlobals.h. Trampoline infrastructure kept
(still used by active HLE patches in XAPI).

### 10.3 — Remove remaining GL fields from VertexAttribute  ✅ DONE

Committed as bfead6028. All GL fields (`needs_conversion`, `converted_buffer`,
`converted_elements`, `converted_size`, `converted_count`, `gl_count`, `gl_type`,
`gl_normalize`, etc.) removed from `VertexAttribute` and `PGRAPHState`.

---

## Step 11: Prepare for Vulkan Backend (No Vulkan Code)

### 11.1 — Abstract the rendering backend interface  ✅ DONE

Introduced `PgraphBackend` struct in `src/devices/video/nv2a_pgraph_backend.h`
with function pointers for draw, draw_state_update, draw_clear, draw_patch,
flip_stall, zpass_begin/end/collect, and launch_transform_program. Replaced
the scattered `pgraph_draw_*` global function pointers with a single
`g_pgraph_backend` instance. D3D11 backend populates it in
`D3D11_init_pgraph_plugins()`.

### 11.2 — Isolate D3D11 code into backend module  ✅ DONE (partial)

Moved all D3D11-specific implementation files into `Backend/`:
- `HostDevice.cpp`, `HostImGui.cpp`, `HostRender.cpp`, `HostResource.cpp`,
  `HostResourceCreate.cpp`, `HostResourceUpload.cpp`, `HostSync.cpp`,
  `HostWindow.cpp`, `PatchDraw.cpp/.h`

Backend-agnostic files remain in `Rendering/`:
- `RenderGlobals.cpp/.h`, `NV2A_PGRAPH_Helpers.cpp/.h`,
  `IndexBufferConvert.cpp/.h`, `WalkIndexBuffer.cpp/.h`, `EmuD3D8_common.h`

`RenderGlobals.h` split complete: D3D11-typed declarations (device pointers,
query helpers, resource cache types, SetHost*/GetHost* functions) moved to
`Backend/Backend_D3D11.h`.  `RenderGlobals.h` now contains only backend-agnostic
state (Xbox types, resource keys, format arrays, window globals, trampolines).

### 11.3 — Make PGRAPH state the single source of truth  ✅ DONE

Verified:
- `RCInterpreterCBLayout` eliminated — `pg->regs[]` uploaded as raw
  `StructuredBuffer<uint>` SRV; shader indexes with `NV_PGRAPH_*` offsets
- `VSInterpreterCBLayout` eliminated — `pg->xf.xfpr[]` uploaded as
  `StructuredBuffer<uint4>` SRV (`g_XFPR` at t5)
- No references to `XboxRenderStates` in render path (deleted)
- `VertexFetchLayoutCB` remains for per-draw params (PrimType, IndexedDraw,
  etc.) that have no PGRAPH register equivalent
- `g_pXbox_SetTexture[]` remains as texture metadata side-map (PGRAPH doesn't
  carry full Xbox texture header in registers)

### 11.4 — Add Vulkan SDK to CMakeLists.txt (optional prep)  ✅ DONE

Added `find_package(Vulkan QUIET)` with status message. Does not link or
require the SDK — just reports availability. Builds cleanly with and without.

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
