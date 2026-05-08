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

### Known Issues & Commits

See [rendering_test_status.md](rendering_test_status.md) for the consolidated known issues, fixed issues, and commit tracking.

## JIT ↔ Interpreter Remaining Differences
| Area | Difference | Impact |
|------|-----------|--------|
| nv2a_mul zero×inf | Both return 0 (xemu returns NaN) | Correct for NV2A hardware |
