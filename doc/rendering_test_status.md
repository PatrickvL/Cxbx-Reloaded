# Rendering Test Status

Consolidated tracking of XDK sample rendering status, known issues, and game compatibility.
This is the single source of truth — other docs should reference this file instead of duplicating status tables.

Last updated: May 2026 (dx11 branch, commit 6cc383454)

> **Note:** ⚠️ and ❌ entries are from automated 20-second captures and need further manual confirmation/investigation.

---

## XDK Sample Status (93 Compiled Samples)

### Legend
- ✅ Renders correctly (may have minor issues)
- ⚠️ Renders with visible bugs
- ❌ Broken (black screen, missing geometry, crash)
- 🔲 Not yet tested on dx11 branch

### Tutorial Samples
| Sample | Status | FPS | Notes |
|--------|--------|-----|-------|
| Vertices | ✅ | 10.32 | RGB triangle |
| Textures | ✅ | 63 | Textured cylinder on blue bg |
| Lights | ✅ | 13.65 | Yellow lit cylinder, specular highlight visible |
| Matrices | ✅ | 9.87 | Rotating colored triangle |
| Meshes | ✅ | 7.96 | Textured tiger on blue bg |
| CreateDevice | ✅ | 6.99 | Magenta clear color (expected behavior) |

### Vertex Shader Samples
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| VertexShaders | ✅ | 35 | VS program execution | Spinning triangle renders correctly |
| Ripple | ✅ | 0.34 | VS program, animated mesh | Yellow ripple pattern on blue bg, circular ripples visible (very slow) |
| CompressedVertices | ✅ | 37 | NORMPACKED3 decode | Yellow teapot renders correctly |
| MatrixPaletteSkinning | ✅ | 58.70 | ARL, blend weights, multi-stream | Green snake/worm with 20 bones, high FPS |
| VertexBlend | ✅ | 31 | Blend weights | "Microsoft" wobble text with vertex blending |
| DisplacementMap | ⚠️ | 63 | D3DCOLOR streams, SGE, ARL | Purple bg + text only, no displaced mesh visible |
| FastVSConstants | ✅ | 7 | Constant upload perf | Triangles + perf counter render correctly |
| MultiShader | ✅ | 0.06 | Multiple VS programs | Color gradient quad rotating, very low FPS |
| ShaderSplicer | ✅ | 29 | VS program switching | Torus with reflections, lighting correct |
| StateShader | ✅ | 7.24 | VS state shaders | RGB gradient colored triangle on gray bg |
| VSLights | ✅ | 12.28 | VS-based lighting | Multi-colored lit teapot (4 lights, Global optimization) |

### Fixed-Function Lighting & Material
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| Dolphin | ✅ | 8.11 | FF lighting, fog, multi-stream morphing | Dolphin+ocean+sky |
| DolphinClassic | ✅ | 4.52 | FF lighting, fog | Dolphin swimming over ocean floor, nice water |
| DolphinHDTV | ✅ | 0.32 | FF lighting, specular, caustics | Same scene, very low FPS (HDTV mode overhead?) |
| TwoSidedLighting | ✅ | 8.50 | Two-sided lighting, back material alpha | Cylinder with front/back colors |
| Strip | ✅ | ~1 | `Power=16`, specular, reflection texgen, aniso | Robot mesh with perf stats |
| Fur | ✅ | 15.87 | `Power=40`, specular, alpha test | Bears with fur shells now rendering |
| Lensflare | ⚠️ | 4.08 | `SPECULARMATERIALSOURCE`, point light | Trees/grass visible, no sun/corona rendered but occlusion detection works |
| Minnaert | ✅ | 47 | Custom lighting model via PS | Female model with Minnaert lighting |
| PerPixelLighting | ✅ | 9.35 | Per-pixel lighting | Earth globe with per-pixel lighting, crosshair cursor |
| PerPixelLightingVS | ⚠️ | 24 | Per-pixel lighting via VS | Globe renders but incorrect bright blue band around it |

### Fog Samples
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| Fog | ✅ | 63 | D3DFOG_LINEAR, EXP, EXP2, table/vertex fog | Columns with visible fog fade |
| AlphaFog | ✅ | 6.54 | Fog + alpha blending | Columns on textured ground with fog effect visible |
| VolumeFog | ✅ | 13.01 | Volume fog via PS | Swamp scene with dead trees, dark fog pool, purple sky |
| HeatShimmer | ⚠️ | 0.06 | D3DFOG_LINEAR + distortion PS | Effect incompletely drawn, lacks backdrop, very low FPS |

### Texgen & Environment Mapping
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| CubeMap | ✅ | 95 | `TCI_CAMERASPACEREFLECTIONVECTOR` | Teapot with full cubemap environment reflections |
| SphereMap | ✅ | 17 | `TCI_CAMERASPACENORMAL` | Sphere-mapped teapot on gradient background |
| FresnelReflect | ✅ | 44 | Reflection texgen | Reflective teapot renders correctly |
| ProjectedTexture | ⚠️ | 12.35 | `TCI_CAMERASPACEPOSITION` (eye-linear texgen) | Blue screen with spotlight thumbnail only, no projected texture on geometry |
| UserClipPlane | ✅ | 30 | Clip planes via texgen | Yellow teapot with clip planes applied correctly |
| MirrorClip | ✅ | 7.21 | Clip planes via texgen | Room with mirror, torus, sphere, ellipsoid visible |

### Pixel Shader & Register Combiners
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| PixelShader | ⚠️ | 28 | RC interpreter, direct + file-based PS | Robot model visible but textures partially black |
| ModifyPixelShader | ⚠️ | — | Runtime PS modification, fog | Blue cylinder on blue bg, no menu bar visible |
| DotProduct3 | ✅ | 0.29 | DOT product in combiners, `D3DFVF_SPECULAR` | 3D face with bump normal map, very low FPS |
| Cartoon | ✅ | 12 | Toon shading via PS | Cel-shading and edge rendering correct |
| QuadLerp | ⚠️ | 5.47 | 4-way lerp blending PS | Blue gradient screen with text only, no visible quad lerp geometry |
| Explosion | ⚠️ | 12 | PROJECT2D, DOT_ZW | Landscape blinks, billboard explosions don't render |
| NoSortAlphaBlend | ✅ | 4.62 | Alpha peel PS, constant PS | Alpha-blended shapes (teapot, spheres, rings) |
| ZSprite | ⚠️ | 13 | Z-sprite PS | Teapots render but depth is disregarded |

### Bump Mapping
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| BumpEarth | ✅ | 11.18 | `TCI_CAMERASPACENORMAL`, BUMPENVMAP, PROJECT2D | Earth globe with bump+cloud |
| bumpearth2 | — | — | Variant of BumpEarth | SKIP: XBE not found |
| BumpDemo | ✅ | 9.65 | BUMPENVMAP, BUMPENVMAP_LUM, DOT_RFLCT_SPEC | Bump-mapped gears with glossmap on blue bg |
| BumpLens | ✅ | 10.53 | `TCI_CAMERASPACEPOSITION`, bump lens distortion | Landscape photo with lens distortion overlay |
| HighQualityBumpMapping | ✅ | 9.91 | 2D dependent texture lookup for specular (1–1000) | Torus with per-pixel bump+specular on blue bg |
| NormalMapGeneration | ⚠️ | 9.01 | Normal map creation | Title+text only, no 3D model visible |
| PolynomialTextureMaps | ✅ | 8.94 | PTM via PS | Rocky moon-like sphere, looks correct |

### Point Sprites & Billboards
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| PointSprites | ✅ | 9.47 | GS-based point sprites | Glowing orange sprites on dark ground |
| Billboard | ✅ | 5.56 | Billboarded trees, alpha test | Forest of billboarded trees |
| PaintEffect | ⚠️ | 27 | Point sprite paint | Only works when pressing A, then fades back to black |
| VolumeSprites | ✅ | 45.68 | Point sprites + floor | Bright particle explosion fountain on blue bg |

### Shadow & Stencil
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| ShadowBuffer | ❌ | — | Shadow compare, depth texture, `.xpu` PS | Cxbx-Reloaded splash screen — sample crashed/didn't load |
| ShadowVolume | ✅ | 7.15 | D3DFOG_LINEAR, stencil volumes | Biplane over mountains with shadow volume info text |
| StencilDepth | ✅ | ~8 | Stencil operations, DPH, SLT, RCC | Helicopter with rotor on blue/purple bg |
| StencilMirror | ✅ | ~9 | D3DFOG_LINEAR, stencil mirror | Biplane with glowing lights |

### Fire, Water & Effects
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| Fire | ✅ | 15.88 | `Fire.xpu` pixel shader | Fire effect with ground plane, purple sky |
| Water | ⚠️ | 35 | Bumpenvmap, reflection, aniso=4, fog, `water.xpu` | First frame has correct textures (no water); once water renders, all whitewashed |
| Glass | ⚠️ | 23 | Alpha test, `Glass.xpu` | Teapot with refraction + reflection, but missing skybox background |
| FocusBlur | ⚠️ | 63 | DPNDNT_AR/GB, 5 pixel shaders | Visible geometry but garbled/blocky, checker patterns |
| MotionBlur | ✅ | 2.94 | Motion blur with alpha test | Moon over purple horizon with motion blur |
| VolumeLight | ⚠️ | 19 | PROJECT2D, PROJECT3D | Stonehenge scene rendered, light beam present but could be brighter |
| XRay | ✅ | 0.06 | `TCI_CAMERASPACENORMAL`, X-ray effect | Translucent blue/purple x-ray robot, very low FPS |
| FuzzyTeapot | ⚠️ | 9.49 | Alpha test, fuzzy material | Teapot visible but looks spiky/exploded, fuzz layers not rendering properly |

### Texture & Format
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| VolumeTexture | ⚠️ | 44 | 3D textures | Volume-textured shape renders |
| Swizzle | ⚠️ | 26 | Texture swizzling | Only draws when model is rotated |
| Tiling | ✅ | 63 | Tiled textures | Text info screen showing tile configuration |
| XPRViewer | ⚠️ | 37 | XPR texture format viewer | No geometry rendered, not even when rotating |
| DynamicGamma | ✅ | 0.06 | Gamma correction | Castle scene with gamma histogram/ramp, very low FPS |

### Rendering Infrastructure
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| PushBuffer | ✅ | 6.31 | PGRAPH attribute path, push buffer draws | Two rainbow gradient triangles on blue bg |
| BeginPush | ✅ | 63.01 | Inline push path | Colorful ribbons/strips rendered correctly |
| BackBufferScale | ✅ | 7.35 | Backbuffer scaling | Blue screen with text (expected — info display only) |
| PersistDisplay | ✅ | 8.44 | Display persistence, fog | Dolphin over ocean with "press button to persist" text |
| SwapCallback | ✅ | ~7 | Swap chain callbacks | VBlank timing info display |
| AntiAlias | ✅ | — | MSAA, `CarOpaque.xpu`, `CarTransparent.xpu` | Renders correctly, may crash on mode switch |
| FieldRender | ✅ | 18 | Interlaced rendering, Dolphin scene, fog, specular | Dolphin scene correct, stutters |
| Patch | ⚠️ | 13 | N-patches / higher-order surfaces | Draws only when rotated, missing triangles (tessellation bug) |
| PerfTest | ✅ | 52.22 | Performance benchmarking | Stonehenge scene rendered, purple sky |
| BenchMark | ✅ | 4.58 | Benchmark suite | Repeating colored parallelogram pattern |
| VisibilityTest | ✅ | 4.69 | Occlusion queries (zpass) | Red textured quad with "Sphere not rendered" text |
| PlayField | ✅ | 7.71 | Aniso=4, gameplay prototype | Green grass field with purple sky |
| HighDynamicRange | ✅ | 6.14 | HDR rendering, hot blur PS | Night village scene with bloom, stars, moon |
| SkyBox | ⚠️ | 49 | Skybox rendering | Sky gradient renders, uncertain if terrain/scene below is complete |

### Misc Samples
| Sample | Status | FPS | Key Features | Notes |
|--------|--------|-----|-------------|-------|
| TrueTypeFont | ⚠️ | 32 | Font rendering | Only labels rendered, no demo font text visible |
| Notifier | ✅ | 3.58 | GPU notification mechanism | Biplane over ground with fence notifier stats |
| Trees | ❌ | — | Billboard trees, fog, alpha test, aniso | Needs more warmup time (killed too soon), black screen |
| Gamepad | ✅ | 9.59 | Input only (no rendering) | Controller diagnostic screen fully rendered |
| Rumble | ✅ | 19.91 | Input only (no rendering) | Motor Test Page with Left/Right 0% indicators |

---

## Known Rendering Issues

### Active Issues
| Issue | Sample(s) | Description | Root Cause | Priority |
|-------|-----------|-------------|-----------|----------|
| ~~Black teapot~~ | ~~CubeMap~~ | ~~Teapot lacks environment reflection~~ | ~~Fixed: PGRAPH RT cache stored faces as ArraySize=1; composed into cubemap~~ | ~~Fixed~~ |
| ~~No sphere~~ | ~~SphereMap~~ | ~~Teapot only visible when rotated~~ | ~~Fixed: FF VS constant upload cache skipped re-upload when shared CB was overwritten by VP draws~~ | ~~Fixed~~ |
| Missing skybox | Glass | Cubemap environment background not rendered | Skybox draw may not go through RT-as-cubemap composition path | Medium |
| Whitewashed scene | Water | First frame correct, then whitewashed once water renders | Unknown — bumpenvmap or render-to-texture issue | Medium |
| Garbled output | FocusBlur | Visible geometry but garbled/blocky checker | DPNDNT_AR/GB dependent texture lookup broken | Medium |
| Blue band | PerPixelLightingVS | Globe renders but has incorrect bright blue band | Unknown VS lighting issue | Medium |
| ZSprite depth | ZSprite | Teapots render but depth is disregarded | Z-sprite depth handling broken | Medium |
| Black textures | PixelShader | Robot visible but textures partially black | Texture binding issue in multi-texture PS sample | Medium |
| No projected texture | ProjectedTexture | Title only, no geometry or projected texture | Eye-linear texgen not projecting correctly | Medium |
| Spiky teapot | FuzzyTeapot | Teapot looks spiky/exploded | Fuzz layers not rendering properly | Low |
| Paint fades | PaintEffect | Only works when pressing A, then fades back to black | Point sprite paint not persisting | Low |
| Missing triangles | Patch | Draws only when rotated, missing triangles | Tessellation bug in N-patch implementation | Low |
| Help screen overlay | All XDK samples | Duke image overlay doesn't render, pause when opening | Unknown HLE/overlay issue | Low |
| Crash | ShadowBuffer | Sample crashes during init | Unknown crash during init | Medium |
| Slow warmup | Trees | Black screen — needs more warmup time | Killed before rendering started | Low |
| Very low FPS | HeatShimmer, XRay, DolphinHDTV | Renders but <1 FPS | Unknown perf issue — possibly shader compilation or fallback path | Low |
| MSAA disabled | All | Aliased edges everywhere | D3D11 MSAA not implemented | Low |

### Fixed Issues
| Issue | Sample(s) | Fix Commit | Description |
|-------|-----------|-----------|-------------|
| Gray dolphin | Dolphin | — | `g_pXbox_PixelShader` was NULL during puller draws → wrong PS path. Fix: COMBINECTL != 0 selects RC interpreter |
| Viewport flash | Dolphin | 2673cf262 | Scene randomly rendered into small top-left box; stale change-detection fast path removed |
| Black screen | Multiple | — | SetHostResource missing D3DUsage → back buffer recreated blank every frame |
| Point sprites broken | PointSprites | cad5d46f8, a44400fd2 | Blend factor, GS, textures, sizing; SETUPRASTER register fix |
| Specular power = 0 | Fur, Strip, all FF lighting | 6cc383454 | Power was never read from NV2A LTC1 → `pow(x,0)=1` for all specular |
| Material alpha = 1 | TwoSidedLighting, Glass | 6cc383454 | Material alpha not read from `ltctxa[CM_COL][3]` |
| Fog depth uninitialized | Fog, AlphaFog | 6cc383454 | `fogDepth` could have garbage if no mode matched |
| Texgen planes ignored | ProjectedTexture, UserClipPlane | 6cc383454 | EYE_LINEAR/OBJECT_LINEAR used raw position without TG*MAT multiply |
| ABS_PLANAR fog | Fog, AlphaFog | 6cc383454 | ABS_PLANAR mapped same as PLANAR without `abs()` |
| DOT_RFLCT_DIFF | — | 70b30da8e | PS texture mode fix |
| CLIPPLANE direction | — | 70b30da8e | PS clip plane direction |
| Shadow compare ordering | ShadowBuffer | 70b30da8e | Compare before PostProcess (note: sample now crashes on init) |
| RT-as-texture rebinding | — | 938ab3e3d | SRV rebinding after RT switch |
| RCC sign preservation | — | 39fd1b808 | VS JIT sign fix |
| LOG(0) handling | — | 39fd1b808 | VS JIT LOG fix |
| PS PASSTHRU PostProcess | PixelShader | 10bc68206 | PostProcess not applied to PASSTHRU mode |
| TCI_OBJECT texgen | — | 3455abb61 | Object-space texgen not implemented |
| RT cubemap composition | CubeMap, Glass | — | PGRAPH RT cache stored cubemap faces as separate ArraySize=1 textures; composed 6 faces into proper TEXTURECUBE |
| FF constant cache stale | SphereMap, Swizzle | — | FF VS constant upload cache (`memcmp`) skipped re-upload when unchanged, but VP draws between frames overwrote the shared D3D11 constant buffer |
| Filter construction | — | c37cc4fa3 | Use D3D11_ENCODE_BASIC_FILTER for correctness |

---

## Key Rendering Commits (dx11 branch)

| Commit | Description | Impact |
|--------|-------------|--------|
| 6cc383454 | FF VS: specular power, material alpha, fog init, texgen planes | Lighting, fog, texgen |
| 10bc68206 | PS JIT: PostProcess for PASSTHRU mode | Pixel shader output |
| 3455abb61 | FF VS: TCI_OBJECT texgen mode | Texgen |
| c37cc4fa3 | Samplers: D3D11_ENCODE_BASIC_FILTER | Texture filtering |
| 70b30da8e | PS: DOT_RFLCT_DIFF, CLIPPLANE, shadow compare | Pixel shader modes |
| 938ab3e3d | D3D11: RT-as-texture SRV rebinding | Render target textures |
| cad5d46f8 | Point sprite rendering (blend, GS, textures, sizing) | Point sprites |
| a44400fd2 | PointSpriteEnable register (SETUPRASTER) | Point sprites |
| 2673cf262 | Remove stale viewport change-detection fast path | Viewport flash |
| 39fd1b808 | VS JIT: RCC sign, LOG(0) | Vertex shader accuracy |
| 5bd6999d4 | GPU→CPU render target readback (VEH fault handler) | Shadow maps, readback |
| 49cd110f5 | DMA context A/B for texture/palette address | Texture addressing |
| d795c73fc | NV097_LAUNCH_TRANSFORM_PROGRAM | Transform programs |
| 09c536aa2 | Shader JIT reorganization with disk cache | Performance |
| 7c8b275c5 | Texture cache dirty check fix | Stale textures |
| a2cfd7374 | YUY2/UYVY texture conversion | FMV/video |
| 160e8b55c | Linear texture pitch < 64 fix | Narrow textures |

---

## Game Rendering Status

### Games with Code Test Cases
These games are referenced in code comments as test cases for specific features:

| Game | Feature(s) Tested | Code Reference |
|------|-------------------|----------------|
| Halo: Combat Evolved | Memory alloc, DES crypto, DS stream, vertex types (NORMSHORT1/2) | VMManager, EmuDes, DirectSound, XbVertexShaderDecoder |
| Halo 2 | Push buffer HLE patch, media loading | EmuPatches_Shader, XFileMediaObject |
| GTA III | Vertex shader constant upload (HLE draws), DS stream sizing | HostRender, DirectSoundStream |
| GTA: San Andreas | Vertex distortion (compact vertex formats via push buffers) | progress.md investigation |
| Panzer Dragoon Orta | File I/O, input polling delay | EmuKrnlNt, EmuKrnlIo, Xapi |
| Jet Set Radio Future | Input enumeration, thread suspend/resume, vertex types (PBYTE4) | Xapi, EmuKrnlKe, XbVertexShaderDecoder |
| Dead or Alive 3 | XBE section loading, compute shader tracking, memory protect | EmuKrnlXe, Backend_D3D11_Compute, EmuKrnlMm |
| Dead or Alive Ultimate | PAL VBlank timing, HAL shutdown | nv2a, EmuKrnlHal |
| Burnout / OutRun 2006 | Surface-as-texture binding, DS buffer conversion | HostSync, HostResource, DirectSoundBuffer |
| Burnout 3 | DS stream pause + SetFormat | DSStream_PacketManager |
| Fable | x86 TLB page fault, memory alloc at XBE_MAX_VA | EmuX86, VMManager |
| Forza Motorsport | NtQueryVirtualMemory 2GB iteration | EmuKrnlNt |
| Splinter Cell 1 & 2 | DES crypto, input enumeration | EmuDes, Xapi |
| Star Wars Battlefront | Memory allocation, SetVertexShaderConstant LTCG | VMManager, EmuPatches_Shader |
| Star Wars: Jedi Academy | XInputSetState | Xapi |
| Lego Star Wars | Input polling delay | Xapi |
| Steel Battalion | Thread APC delivery, custom controller | EmuKrnlNt, EmuKrnlKe, InputDevice |
| Shenmue II | x86 port I/O, screen scale (XYZRHW) | EmuX86, HostDevice |
| Cel Damage | Vertex type NORMSHORT3 | XbVertexShaderDecoder |
| Oddworld: Stranger's Wrath | Input enumeration | Xapi |
| Amped | PsCreateSystemThread, vertex shader program flag | EmuKrnlPs, XbD3D8Types_Resources |
| RalliSport Challenge | RDTSC patching, XeLoadSection | PatchRdtsc, EmuKrnlXe |
| Turok | Inline vertex elements, SHORT3/PBYTE3 types | EmuNV2A_PGRAPH, XbVertexShaderDecoder |
| Baldur's Gate: Dark Alliance 2 | NORMSHORT2 vertex type | XbVertexShaderDecoder |
| Hunter: The Reckoning | Inline vertex elements (ARRAY_ELEMENT16) | EmuNV2A_PGRAPH |
| Otogi | Inline vertex elements | EmuNV2A_PGRAPH |
| Alter Echo | PCRTC raster position read | EmuNV2A_PCRTC |
| Prince of Persia: WW | RunPushBuffer HLE patch | EmuPatches_Shader |
| NASCAR Heat 2002 | SetVertexShader LTCG patch | EmuPatches_Shader |
| Kung Fu Chaos | RunPushBuffer patch | EmuPatches_Shader |
| Aggressive Inline | SetVertexShaderConstant LTCG patch | EmuPatches_Shader |
| Soldier of Fortune II | Vertex shader constant upload (HLE draws) | HostRender |

---

## XDK Sample Feature Coverage

Feature usage from XDK Graphics samples (for test planning):

### Specular Power (Material.Power set explicitly)
Fur (40.0), Strip (16.0)

### Fog Modes (comprehensive)
Fog (LINEAR/EXP/EXP2), AlphaFog, Dolphin, DolphinClassic, DolphinHDTV, HeatShimmer, ModifyPixelShader, ShadowVolume, StencilMirror, Water

### Hardware Texgen (non-passthrough TCI modes)
| Mode | Samples |
|------|---------|
| TCI_CAMERASPACENORMAL | BumpEarth, SphereMap, XRay |
| TCI_CAMERASPACEPOSITION | ProjectedTexture, BumpLens, Tut05_Textures |
| TCI_CAMERASPACEREFLECTIONVECTOR | CubeMap, Strip |
| TCI_OBJECT | No XDK sample (game-specific only) |

### Anisotropic Filtering
PlayField (4), Strip (3), Trees (2), Water (4)

### Alpha Test
Billboard, BumpDemo, DotProduct3, FocusBlur, Fur, FuzzyTeapot, Glass, Lensflare, MotionBlur, PaintEffect, PlayField, PointSprites, Trees, VolumeSprites

---

## Testing Procedure

1. Build: `cmake --build build --config Release --target cxbxr-emu`
2. Deploy HLSL: copy `src/.../Shaders/*.hlsl*` → `build/bin/Release/hlsl/`
3. Clear shader cache: `Remove-Item "build/bin/Release/ShaderCache/*" -Recurse -Force`
4. Capture: `tools/capture_emulator.ps1 -XbePath "path.xbe" -EmulatorPath "build/bin/Release/cxbx.exe" -DelayMs 15000 -StopAfter`
5. Environment: set `CXBX_XBE_SAMPLES` user env var to compiled XBE directory
6. XBE layout: `%CXBX_XBE_SAMPLES%/<SampleName>/<SampleName>.xbe`
