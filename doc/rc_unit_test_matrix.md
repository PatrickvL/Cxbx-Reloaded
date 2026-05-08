# Register Combiner (Pixel Shader) Test Matrix

## Status Key
- ✅ Correct — matches expected output
- ⚠️ Minor issue — renders but with visible artifact
- ❌ Broken — wrong output or crash
- 🔲 Untested
- N/A — sample doesn't exercise this feature

## PS Operations → XDK Sample Mapping

### Texture Modes

| Mode | Samples Using This Mode | JIT Status | Interp Status | Notes |
|------|----------------------------------------------|------------|---------------|-------|
| NONE | Dolphin, DolphinClassic, DolphinHDTV, Explosion, FieldRender, ModifyPixelShader, NoSortAlphaBlend, PerfTest, PersistDisplay, PixelShader, VolumeFog, VolumeLight, ZSprite, BumpMapping, Minnaert | ✅ | ✅ | Basic passthrough |
| PROJECT2D | Dolphin, DolphinClassic, DolphinHDTV, Explosion, FieldRender, ModifyPixelShader, PerfTest, PersistDisplay, PixelShader, VolumeFog, VolumeLight, ZSprite, BumpMapping | ✅ | 🔲 | Projective texturing |
| PROJECT3D | NoSortAlphaBlend, PerfTest, VolumeLight | ✅ | 🔲 | 3D volume projection |
| CUBEMAP | BumpMapping, EnvMapping, Glass, MatrixPaletteSkinning, Minnaert, PerPixelLighting, PerPixelLightingVS, SkyBox, Trees, Water | ✅ | 🔲 | Environment mapping |
| PASSTHRU | BumpMapping (confirmed in code) | ✅ | ⚠️ | JIT skips PostProcess; Interp applies it |
| CLIPPLANE | UserClipPlane, MirrorClip | ✅ | 🔲 | VS clip distance needed for teapot |
| BUMPENVMAP | BumpMapping, Water (D3DTSS_BUMPENVMAT) | ✅ | 🔲 | Bump matrix offset |
| BUMPENVMAP_LUM | BumpMapping (D3DTSS_BUMPENVLSCALE) | ✅ | 🔲 | + luminance scale/bias |
| DOTPRODUCT | BumpMapping, Dolphin, DolphinClassic, DolphinHDTV, Explosion, FieldRender, Minnaert, ModifyPixelShader, NoSortAlphaBlend, PersistDisplay, PixelShader, UserClipPlane, VolumeFog, ZSprite | ✅ | 🔲 | Base for DOT_* modes |
| DOT_ST | Minnaert (PS_TEXTUREMODES_DOT_ST) | ✅ | 🔲 | Dot-product texcoord gen |
| DOT_ZW | Explosion, NoSortAlphaBlend, ZSprite | ✅ | 🔲 | Xemu INCOMPLETE; we have div-by-zero guard |
| DOT_RFLCT_DIFF | (not found in source scan) | ? | 🔲 | Fixed: xemu normal=(prev.x, cur, next) |
| DOT_RFLCT_SPEC | BumpMapping (PS_TEXTUREMODES_DOT_RFLCT_SPEC) | ✅ | 🔲 | Full implementation (xemu UNIMPLEMENTED) |
| DOT_STR_3D | (not found in source scan) | ? | 🔲 | 3-stage dot→3D |
| DOT_STR_CUBE | (not found in source scan) | ? | 🔲 | 3-stage dot→cube |
| DPNDNT_AR | FocusBlur (texREG2AR) | 🔲 | 🔲 | Dependent read via .ar |
| DPNDNT_GB | FocusBlur (texREG2GB) | 🔲 | 🔲 | Dependent read via .gb |
| BRDF | Minnaert (PS_DOTMAPPING_ZERO_TO_ONE used as BRDF) | ✅ | 🔲 | Xemu UNIMPLEMENTED; we have impl |

### Register Combiner Operations

| Operation | Samples Using This Operation | JIT Status | Interp Status | Notes |
|-----------|----------------------------------------------|------------|---------------|-------|
| Basic MUX (MSB select) | Dolphin, DolphinClassic, DolphinHDTV, FieldRender, Minnaert, ModifyPixelShader, NoSortAlphaBlend, PerfTest, PersistDisplay, PixelShader, VolumeFog, VolumeLight | ✅ | 🔲 | |
| DOT product (AB) | PerPixelLighting (via PS_DOTMAPPING) | ✅ | 🔲 | |
| Bias/expand input map | Minnaert (PS_INPUTMAPPING_EXPAND), Fur/Trees (EXPAND in macros) | ✅ | 🔲 | |
| Scale 2x (SHIFTLEFT_1) | VolumeFog (PS_COMBINEROUTPUT_SHIFTLEFT) | ✅ | 🔲 | |
| EF product (final) | Dolphin, DolphinClassic, DolphinHDTV, FieldRender, PersistDisplay, BumpMapping | ✅ | 🔲 | Final combiner E×F |
| V1R0_SUM (specular add) | BumpMapping | ✅ | 🔲 | |
| Final combiner (EFG) | Most PS samples (PS_COMBINERCOUNT implies final combiner) | ✅ | 🔲 | |
| FOG alpha | Dolphin, Fog, HeatShimmer, Water, +6 more (D3DRS_FOGENABLE) | ✅ | ⚠️ | Pre-existing flicker |

### Shadow/Depth Features
| Feature | XDK Sample(s) | JIT Status | Interp Status | Notes |
|---------|---------------|------------|---------------|-------|
| Shadow compare | ShadowBuffer (PS_TEXTUREMODES used with shadow) | ✅ | 🔲 | 🔲 | Fixed: compare before PostProcess |
| Depth texture | ShadowBuffer | ✅ | 🔲 | 🔲 | |
| Stencil buffer | StencilDepth, StencilMirror | ✅ | 🔲 | 🔲 | D3DRS_STENCILENABLE |

### Recommended Minimal PS Test Set

| Sample | XBE Path | PS/RC Features Covered | Status |
|--------|----------|----------------------|--------|
| BumpDemo | `BumpDemo\BumpDemo.xbe` | BUMPENVMAP, BUMPENVMAP_LUM, DOT_RFLCT_SPEC, PASSTHRU, DOTPRODUCT, V1R0_SUM, EF_PROD | 🔲 |
| Dolphin | `Dolphin\Dolphin.xbe` | PROJECT2D, MUX_MSB, EF_PROD, FOG, NONE | ✅ 8s OK |
| Minnaert | `Minnaert\Minnaert.xbe` | DOT_ST, CUBEMAP, DOTPRODUCT, EXPAND mapping, MUX_MSB, BRDF | 🔲 |
| Explosion | `Explosion\Explosion.xbe` | PROJECT2D, DOT_ZW, DOTPRODUCT | 🔲 |
| VolumeFog | `VolumeFog\VolumeFog.xbe` | PROJECT2D, SHIFTLEFT, MUX_MSB | 🔲 |
| VolumeLight | `VolumeLight\VolumeLight.xbe` | PROJECT2D, PROJECT3D | 🔲 |
| UserClipPlane | `UserClipPlane\UserClipPlane.xbe` | CLIPPLANE, DOTPRODUCT | 🔲 |
| PerPixelLighting | `PerPixelLighting\PerPixelLighting.xbe` | CUBEMAP, DOT product in combiners | 🔲 |
| ShadowBuffer | `ShadowBuffer\ShadowBuffer.xbe` | Shadow compare, depth texture | ✅ 8s OK |
| FocusBlur | `FocusBlur\FocusBlur.xbe` | DPNDNT_AR/GB (via texREG2AR/texREG2GB) | 🔲 |

### Known Issues (Pre-existing)
| Issue | Sample | Description | Commit Introduced |
|-------|--------|-------------|------------------|
| Flickering triangles | AlphaFog | Whole triangles disappear per-frame showing background | Pre-existing (before dx11 branch) |
| Missing geometry (flicker) | MatrixPaletteSkinning | Snake body triangles disappear leaving skeleton outline | Pre-existing; same root cause as AlphaFog |
| Missing geometry | VertexBlend | "Microsoft" text mesh has invisible triangles | Pre-existing; same triangle-disappearing pattern |
| Viewport flash | Dolphin | Scene randomly renders into small top-left box then restores | Pre-existing render target/viewport state race |
| White textures | Water (pirate) | Some geometry missing texture data | Pre-existing texture freshness bug |
| Missing teapot | UserClipPlane | VS clip plane computation not implemented | Pre-existing (needs FF texgen) |
| Bright edges | PerPixelLightingVS | Globe edges overbright/blown out | Likely LIT specular or normal issue |
| Black screen | PaintEffect | Point sprite paint not visible | PointSprites fixed (cad5d46f8, a44400fd2) but PaintEffect not yet confirmed |
| Missing floor | VolumeSprites | White sprite fountain visible but no floor texture | Unknown — possibly missing texture or triangle draw |

## Commits Applied
| Hash | Description | Date |
|------|-------------|------|
| 70b30da8e | PS: Fix DOT_RFLCT_DIFF, CLIPPLANE direction, shadow compare ordering | 2026-05-06 |
| 938ab3e3d | D3D11: Fix RT-as-texture SRV rebinding after render target switch | 2026-05-06 |
| cad5d46f8 | D3D11: Fix point sprite rendering (blend, GS, textures, sizing) | 2026-05-07 |

## JIT ↔ Interpreter Remaining Differences
| Area | Difference | Impact |
|------|-----------|--------|
| PASSTHRU PostProcess | JIT skips; Interpreter applies | Unknown — needs test with Water caustics |
| BRDF PostProcess | JIT skips; Interpreter applies | Unknown — no test case yet |
| nv2a_mul zero×inf | Both return 0 (xemu returns NaN) | Correct for NV2A hardware |

## Testing Procedure
1. Clear shader cache: `Remove-Item "$env:APPDATA\Cxbx-Reloaded\ShaderCache" -Recurse -Force`
2. Build: `cmake --build . --config Release --target cxbxr-emu`
3. Launch: `cxbxr-ldr.exe /load "<path>.xbe"`
4. Compare against reference screenshots or expected behavior
5. Toggle PS JIT: No toggle variable exists; comment out the `g_PixelShaderCache.GetShader()` call in XbPixelShaderCompiler.cpp to force interpreter fallback
