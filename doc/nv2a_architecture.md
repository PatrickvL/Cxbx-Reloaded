# NV2A Architecture & DX11 Renderer Design

## 1. Device Overview

### MMIO Dispatch Table (nv2a.cpp)
- `regions[]` array defines MMIO address space routing
- Each engine has offset, size, read/write handlers
- Lookup via `EmuNV2A_Block(addr)` function

### Emulated Engines
1. **PMC** — Power management, interrupt aggregation
2. **PFIFO** — Command FIFO with pusher/puller threads
3. **PGRAPH** — 3D graphics engine (main emulation target)
4. **PTIMER** — Timer/interrupt handling
5. **PCRTC** — CRT controller
6. **PFB** — Memory interface
7. **PRAMIN** — GPU instance memory access
8. **USER** — PFIFO MMIO/DMA submission area
9. **PVIDEO**, **PTV**, **PVPE**, **PBUS**, **PCOUNTER**, **PRMA**, **PSTRAPS**

---

## 2. Register Block Map

Source files: `src/devices/video/nv2a.cpp` (block table), `src/devices/video/nv2a_regs.h` (defines)

| Block | Name | Base Offset | Size | Handler File | Purpose |
|-------|------|-------------|------|--------------|---------|
| 0 | PMC | 0x000000 | 4K | EmuNV2A_PMC.cpp | Card master control, interrupt routing |
| 1 | PBUS | 0x001000 | 4K | EmuNV2A_PBUS.cpp | Bus control, PCI config |
| 2 | PFIFO | 0x002000 | 8K | EmuNV2A_PFIFO.cpp | MMIO/DMA FIFO, caches |
| 4 | PRMA | 0x007000 | 4K | EmuNV2A_PRMA.cpp | Real mode BAR access |
| 5 | PVIDEO | 0x008000 | 4K | EmuNV2A_PVIDEO.cpp | Video overlay engine |
| 6 | PTIMER | 0x009000 | 4K | EmuNV2A_PTIMER.cpp | Timer, clock counters |
| 7 | PCOUNTER | 0x00a000 | 4K | EmuNV2A_PCOUNTER.cpp | Performance counters |
| 8 | PVPE | 0x00b000 | 4K | EmuNV2A_PVPE.cpp | MPEG2 decode |
| 9 | PTV | 0x00d000 | 4K | EmuNV2A_PTV.cpp | TV encoder |
| 10 | PRMFB | 0x0a0000 | 128K | EmuNV2A_PRMFB.cpp | VGA framebuffer alias |
| 11 | PRMVIO | 0x0c0000 | 32K | EmuNV2A_PRMVIO.cpp | VGA seq/graphics ctrl |
| 12 | PFB | 0x100000 | 4K | EmuNV2A_PFB.cpp | Memory interface |
| 13 | PSTRAPS | 0x101000 | 4K | EmuNV2A_PSTRAPS.cpp | Straps readout |
| 14 | PGRAPH | 0x400000 | 8K | EmuNV2A_PGRAPH.cpp | 2D/3D graphics (helpers: NV2A_PGRAPH_Helpers.h) |
| 15 | PCRTC | 0x600000 | 4K | EmuNV2A_PCRTC.cpp | CRT controller |
| 16 | PRMCIO | 0x601000 | 4K | EmuNV2A_PRMCIO.cpp | VGA CRTC/attr ctrl |
| 17 | PRAMDAC | 0x680000 | 4K | EmuNV2A_PRAMDAC.cpp | RAMDAC, cursor, PLL |
| 18 | PRMDIO | 0x681000 | 4K | EmuNV2A_PRMDIO.cpp | VGA palette |
| 19 | PRAMIN | 0x700000 | 1M | EmuNV2A_PRAMIN.cpp | GPU instance memory |
| 20 | USER | 0x800000 | 4M | EmuNV2A_USER.cpp | PFIFO submission |
| - | UREMAP | 0xC00000 | 4M | EmuNV2A_USER.cpp | Mirror of USER |

---

## 3. PFIFO Architecture

- **Pusher thread**: Reads DMA push buffers, validates commands
- **Puller thread**: Executes commands via `pfifo_run_puller()`
- **RAMHT**: Hash table for object/handle lookups (XOR folding, channel ID mixed in)
- **CACHE1**: Command FIFO cache holding 32 entries

### Key Finding: PFIFO DMA IS Initialized in HLE Mode
- Xbox D3D runtime initializes PFIFO during D3DDevice_Create (PUSH0_ACCESS, DMA_PUSH_ACCESS)
- The DMA pusher/puller threads ARE started and DO process commands
- All unpatched D3D calls write to ring buffer → pusher/puller process to PGRAPH
- **CRITICAL**: DMA_GET drain must NOT advance GET=PUT when pusher is enabled

---

## 4. Registers with Read-Side Effects

### PMC
- **NV_PMC_INTR_0 (0x100)**: Computed live — ORs pending_interrupts & enabled_interrupts from PFIFO, PGRAPH, PCRTC, PVIDEO, PTIMER

### PFIFO
- **NV_PFIFO_CACHE1_DMA_GET (0x1244)**: If GET != PUT and pusher cannot run → auto-advance GET to PUT (deadlock prevention in HLE mode)

### PTIMER
- **NV_PTIMER_TIME_0 (0x08)**: Computed dynamically from hardware clock
- **NV_PTIMER_TIME_1 (0x0C)**: Upper portion of timer

### PGRAPH
- **NV_PGRAPH_INTR (0x100)**: Write-to-clear interrupt register
- **NV_PGRAPH_RDI_DATA (0x754)**: Auto-increments RDI_INDEX address on each access

### USER
- **NV_USER_DMA_GET**: Same fast-path logic as PFIFO CACHE1 DMA_GET

---

## 5. PGRAPH Method Dispatch

- Main handler: `pgraph_handle_method(subchannel, method, parameter)`
- Routes by graphics_class (NV097 = Kelvin 3D, NV062 = Surfaces2D, etc.)
- Data-driven dispatch table (`nv2a_method_table.h`) replaces ~80+ switch cases
- `nv097_init_method_table()` called from `pgraph_init()`
- Dispatch always runs (no early break); switch handles side effects only

### Side-Effect Cases (still in switch)
- SHADOW_ZSLOPE, FOG_PARAMS, TRANSFORM_PROGRAM_LOAD/START, TRANSFORM_CONSTANT_LOAD
- SET_DEPTH_MASK / SET_COLOR_MASK (order-dependent cache read)
- NV097_SET_OBJECT (stores to kelvin.object_instance only, no register equivalent)

---

## 6. Xbox D3D API → NV097 Method Mappings

Determined by push buffer analysis of the Meshes sample.

### Deferred State Emission
Xbox D3D runtime does NOT write NV097 methods for most state-setting APIs. These store to CPU-side data structures and are emitted at draw time.

**APIs that write ZERO push buffer commands:**
- SetTransform, SetLight, SetMaterial, SetBackMaterial, LightEnable, SetStreamSource

**APIs that write commands immediately:**

| API | Methods Written |
|-----|----------------|
| SetVertexShader | SET_VERTEX_DATA4UB, VIEWPORT_SCALE, DEPTH_RANGE, ENGINE, VP_START_FROM_ID, VP_UPLOAD_CONST/INST |
| SetTexture | TX_ENABLE, TX_OFFSET, TX_FORMAT (per stage) |

---

## 7. PGRAPH Register → DX11 Mapping

### Architecture Goal
- `PGRAPHState::regs[]` IS the actual PGRAPH hardware state (ground truth)
- All other PGRAPHState fields are HLE artifacts to be eliminated
- End state: RC/VS interpreter shaders read SOLELY from regs[]
- Method table dispatch writes directly to regs[]

### Raw regs[] Buffer Architecture
Upload `pg->regs[]` as raw `StructuredBuffer<uint>` SRV. HLSL indexes with NV_PGRAPH_* offsets.
- RC+VS: `g_PGRegs : register(t12)` — StructuredBuffer<uint>, 8KB (CxbxPGRAPHRegs.hlsli)
- VS program: `g_XFPR : register(t5)` — StructuredBuffer<uint4>, 136 × uint4 (Transform Program RAM)
- Software-computed fields (ColorSign, TexFmtFixup, etc.) → small aux cbuffer (PSAuxCBLayout)

### Shader Control (0x0FB4–0x0FC8)
| Register | Offset | DX11 Action |
|----------|--------|-------------|
| CSV0_D | 0x0FB4 | VS mode select (FIXED/PROGRAM), FF light config, skin, fog, point |
| CSV0_C | 0x0FB8 | VS program start, specular, normalize, lighting, material sources |
| CSV1_B | 0x0FBC | FF texgen T2/T3 |
| CSV1_A | 0x0FC0 | FF texgen T0/T1 |
| CHEOPS_OFFSET | 0x0FC4 | VS program load offset |

### Blend/Pipeline (0x1800–0x1958)
| Register | Offset | DX11 Action |
|----------|--------|-------------|
| ANTIALIASING | 0x1800 | MSAA sample count |
| BLEND | 0x1804 | D3D11_BLEND_DESC (enable/factors/equation) |
| BLENDCOLOR | 0x1808 | OMSetBlendState factor |
| BORDERCOLOR0-3 | 0x180C-1818 | D3D11_SAMPLER_DESC.BorderColor |
| CONTROL_0 | 0x194C | Alpha test, depth enable/func/write, color write mask |
| CONTROL_1 | 0x1950 | Stencil enable/func/ref/masks |
| CONTROL_2 | 0x1954 | Stencil ops |
| CONTROL_3 | 0x1958 | Fog enable/mode, shade mode, provoking vertex |
| SETUPRASTER | 0x1990 | Cull/fill/frontface, polygon offset enable |

### Register Combiner State (0x1880–0x1948)
| Register | Offset | DX11 Action |
|----------|--------|-------------|
| COMBINEFACTOR0+i*4 | 0x1880 | RC interpreter (PGRegs SRV) |
| COMBINEFACTOR1+i*4 | 0x18A0 | RC interpreter |
| COMBINEALPHAI0+i*4 | 0x18C0 | RC interpreter (alpha input select) |
| COMBINEALPHAO0+i*4 | 0x18E0 | RC interpreter (alpha output config) |
| COMBINECOLORI0+i*4 | 0x1900 | RC interpreter (RGB input select) |
| COMBINECOLORO0+i*4 | 0x1920 | RC interpreter (RGB output config) |
| COMBINECTL | 0x1940 | RC interpreter (stage count) |
| COMBINESPECFOG0 | 0x1944 | RC interpreter (final combiner ABCD) |
| COMBINESPECFOG1 | 0x1948 | RC interpreter (final combiner EFG) |

### PSDef ↔ PGRAPH Register Mapping
| PGRAPH Register | Offset | PSDef Field |
|---|---|---|
| COMBINEALPHAI0+i*4 | 0x18C0 | PSAlphaInputs[i] |
| COMBINEALPHAO0+i*4 | 0x18E0 | PSAlphaOutputs[i] |
| COMBINECOLORI0+i*4 | 0x1900 | PSRGBInputs[i] |
| COMBINECOLORO0+i*4 | 0x1920 | PSRGBOutputs[i] |
| COMBINEFACTOR0+i*4 | 0x1880 | PSConstant0[i] |
| COMBINEFACTOR1+i*4 | 0x18A0 | PSConstant1[i] |
| COMBINESPECFOG0 | 0x1944 | PSFinalCombinerInputsABCD |
| COMBINESPECFOG1 | 0x1948 | PSFinalCombinerInputsEFG |
| COMBINECTL | 0x1940 | PSCombinerCount |
| SHADERPROG | 0x199C | PSTextureModes |
| SHADERCTL | 0x1998 | PSInputTexture |
| SHADERCLIPMODE | 0x1994 | PSCompareMode |
| pg->dot_rgb_mapping | (field) | PSDotMapping |

### Texture State (0x19BC–0x1A40, per stage ×4)
| Register | Offset | DX11 Action |
|----------|--------|-------------|
| TEXADDRESS0-3 | 0x19BC-19C8 | D3D11_SAMPLER_DESC wrap modes |
| TEXCTL0_0-3 | 0x19CC-19D8 | Texture bind + sampler (enable, alpha kill, aniso) |
| TEXCTL1_0-3 | 0x19DC-19E8 | Texture creation (pitch) |
| TEXFILTER0-3 | 0x19F4-1A00 | D3D11_SAMPLER_DESC (min/mag, LOD bias) |
| TEXFMT0-3 | 0x1A04-1A10 | Texture creation (format, dim, mips, cubemap) |
| TEXIMAGERECT0-3 | 0x1A14-1A20 | Texture creation (width/height) |
| TEXOFFSET0-3 | 0x1A24-1A30 | Texture VRAM lookup |
| TEXPALETTE0-3 | 0x1A34-1A40 | Palette texture decode |

### Fog/Depth/Shadow
| Register | Offset | DX11 Action |
|----------|--------|-------------|
| FOGCOLOR | 0x1980 | RC interpreter aux CB |
| FOGPARAM0/1 | 0x1984-1988 | VS fog computation |
| POINTSIZE | 0x198C | VS point size output |
| SHADOWCTL | 0x19A4 | RC interpreter depth compare |
| SPECFOGFACTOR0/1 | 0x19AC-19B0 | RC interpreter final combiner C0/C1 |
| ZCLIPMIN/MAX | 0x1A90/1ABC | D3D11 viewport depth range |
| ZOFFSETBIAS/FACTOR | 0x1AA4-1AA8 | D3D11 rasterizer depth bias |
| EYEVEC0-2 | 0x1AAC-1AB4 | RC interpreter aux CB |

### Non-Register State (PGRAPHState fields)
| Field | DX11 Action |
|-------|-------------|
| surface_color/zeta.offset | RT/DS lookup/creation |
| surface_state.clipX/Y/Width/Height | Viewport/scissor |
| surface_state.zetaFormat | DepthStencil format (Z16/Z24S8) |
| primitive_mode | IASetPrimitiveTopology |
| program_data[][] | VS interpreter SRV (136×uint4) |
| vsh_constants[][4] | cbuffer b0 (192×float4) |
| ltctxa/ltctxb/ltc1 | FF VS light uniforms |
| vertex_attributes[] | IA bypass layout CB |
| draw_arrays_start/count | Draw() calls |
| inline_elements/array/buffer | DrawIndexed() / Draw() |

---

## 8. Rendering Pipeline (draw_begin / draw_end)

Triggered by NV097_SET_BEGIN_END method (same pattern as xemu).

### Data Flow
```
Xbox game → Xbox D3D runtime → NV2A PFIFO pushbuffer → PGRAPH regs[]
                                                           ↓
                                              DX11 Renderer reads regs[]
                                              at draw time (BEGIN_END)
                                                           ↓
                                              VS interpreter / FF VS (HLSL)
                                              RC interpreter / FF PS (HLSL)
                                              IA bypass vertex fetch
                                              DX11 pipeline state objects
                                                           ↓
                                              D3D11 Draw() → host GPU
```

### draw_end (BEGIN_END with parameter == 0)
1. **NOP check**: Skip if no color write AND no depth test AND no stencil test
2. **Surface update**: Bind/create RT+DS from surface_color/zeta.offset
3. **Texture bind**: For each stage 0-3, check TEXCTL0 enable → lookup/create → SRV + sampler
4. **Shader bind**: Select VS (interpreter vs FF) based on CSV0_D MODE; always RC interpreter PS
5. **Upload VS state**: vsh_constants[] → cbuffer; regs → SRV; program_data → SRV
6. **Upload PS state**: regs → SRV; aux CB for software-only fields
7. **Pipeline state**: Blend/DepthStencil/Rasterizer from CONTROL_0/1/2, BLEND, SETUPRASTER
8. **Viewport/scissor**: From NV2ASurfaceState clip + ZCLIPMIN/MAX
9. **Vertex bind**: Fill layout CB from vertex_attributes[]; point mirror SRV
10. **Draw dispatch**: draw_arrays / inline_elements / inline_buffer / inline_array

---

## 9. EMUPATCH Classification

**All D3D EMUPATCHes are disabled** (commented out in `Patches.cpp`). Implementations
moved to the unused/dead code dustbin file. The DX11 renderer operates entirely from
PGRAPH register state populated by the native PFIFO pusher/puller pipeline.

### Previously KEEP (now also disabled)
- Direct3D_CreateDevice — host DX11 device init now triggered by PFIFO
- D3DDevice_Reset — host device reset
- D3DDevice_Present/Swap — replaced by NV097_FLIP_STALL handling
- D3DDevice_GetBackBuffer — backbuffer tracked via RT offset matching
- D3DDevice_EnableOverlay/UpdateOverlay — PVIDEO register-driven
- D3DDevice_BeginPush/EndPush — native PFIFO DMA
- D3DDevice_RunPushBuffer — `pfifo_submit_pushbuffer()` callback
- D3D_DestroyResource — host resource cleanup via cache eviction

### REMOVED (state goes through PFIFO→PGRAPH natively)
- All SetRenderState/SetTextureState patches
- SetTexture, SetTransform, SetLight, SetMaterial, LightEnable
- SetVertexShader, SelectVertexShader, DeleteVertexShader
- SetStreamSource, SetRenderTarget, Clear
- SetPalette, CopyRects
- RunVertexStateShader (now NV097_LAUNCH_TRANSFORM_PROGRAM)
- BlockOnTime (now native PFIFO flush + kernel wait)
- SetVertexShaderInput, InsertCallback, BlockUntilVerticalBlank
- Visibility test patches (now pgraph_zpass callbacks)
- All LTCG variants of above
