# Vertex Shader Test Matrix

## Status Key
- ✅ Correct — matches expected output
- ⚠️ Minor issue — renders but with visible artifact
- ❌ Broken — wrong output or crash
- 🔲 Untested
- N/A — sample doesn't exercise this feature

## VS Operations → XDK Sample Mapping

### Instruction Coverage

| Instruction | Samples Using This Instruction | JIT Status | Interp Status | Notes |
|-------------|----------------------------------------------|------------|---------------|-------|
| MOV | Dolphin, Explosion, Fire, Fur, Glass, MatrixPaletteSkinning, PerPixelLightingVS, ShadowBuffer, StencilBuffer(!), Trees, VertexShader, Water, +28 more | ✅ | ✅ | Ubiquitous |
| MUL | Dolphin, Explosion, Fire, Fur, MatrixPaletteSkinning, ShadowBuffer, StencilBuffer(!), TwoSidedLighting, VertexShader, Water, +20 more | ✅ | 🔲 | nv2a_mul (0×inf=0) |
| ADD | Dolphin, Explosion, Fog, Fur, MatrixPaletteSkinning, MotionBlur, ShadowBuffer, Trees, VertexShader, Water, +33 more | ✅ | 🔲 | |
| MAD | Dolphin, Fog, Fur, MatrixPaletteSkinning, SkyBox, StencilBuffer(!), Trees, VertexShader, Water, +11 more | ✅ | 🔲 | nv2a_mul + add |
| DP3 | Dolphin, Explosion, Fog, Fur, MatrixPaletteSkinning, PerPixelLightingVS, ShadowBuffer, StencilBuffer(!), Trees, TwoSidedLighting, VertexShader, Water, +17 more | ✅ | 🔲 | nv2a_dot3 |
| DP4 | Dolphin, Explosion, Fog, Fur, MatrixPaletteSkinning, ShadowBuffer, StencilBuffer(!), Trees, TwoSidedLighting, VertexShader, +16 more | ✅ | 🔲 | nv2a_dot4 |
| DPH | MatrixPaletteSkinning, StencilBuffer(!) | ✅ | 🔲 | Rare; only 2 samples |
| DST | BumpMapping, Fur, Lensflare, Trees | ✅ | 🔲 | Light attenuation |
| MIN | Fur, MatrixPaletteSkinning, Trees, TwoSidedLighting, +23 more | ✅ | 🔲 | |
| MAX | Fog, Fur, Trees, TwoSidedLighting, VertexShader, +31 more | ✅ | 🔲 | |
| SLT | ShadowBuffer, StencilBuffer(!) | ✅ | 🔲 | Rare; -0 fixup |
| SGE | DisplacementMap, Water | ✅ | 🔲 | Rare |
| ARL | MatrixPaletteSkinning (a0-indexed), DisplacementMap, PolynomialTextureMaps, QuadLerp, ZSprite, BumpMapping, Fur, FocusBlur, ShadowBuffer, PixelShader | ✅ | 🔲 | Bias 0.001 |
| RCP | Dolphin, Explosion, Fur, VertexShader, Water, +5 more | ✅ | 🔲 | |
| RCC | Dolphin, DolphinClassic, DolphinHDTV, AlphaFog, FieldRender, Fog, PersistDisplay, StencilBuffer(!), VertexShader, VolumeFog | ✅ | 🔲 | Clamp [2^-64, 2^64] |
| RSQ | Explosion, Fog, Fur, MatrixPaletteSkinning, ShadowBuffer, StencilBuffer(!), VertexShader, Water, +6 more | ✅ | 🔲 | |
| EXP | Fog, Fur | ✅ | 🔲 | Exact floor (no ARL bias) |
| LOG | Fur | ✅ | 🔲 | IEEE bit extraction |
| LIT | Dolphin, Fire, MatrixPaletteSkinning, Minnaert, PerPixelLighting, PerPixelLightingVS, SkyBox, Tutorials, Water, +3 more | ✅ | 🔲 | kMax=127.9961 |

(!) StencilBuffer in Code\Graphics is named StencilBuffer; compiled variants are StencilDepth and StencilMirror.

### Execution Model
| Feature | XDK Sample(s) | JIT Status | Interp Status | Notes |
|---------|---------------|------------|---------------|-------|
| Paired MAC+ILU | Dolphin (dp4+rcc), Fog (mad+exp), Fur (dp3+rsq, dst+rcp, mul+log), StencilBuffer (dp4+rcc, dp3+rsq) | ✅ | ✅ | Both handle RAW hazards |
| Context writes (c[]) | None known in XDK samples | N/A | ✅ | JIT falls back to interp |
| Relative addressing (a0) | MatrixPaletteSkinning, DisplacementMap, PolynomialTextureMaps, QuadLerp, ZSprite | 🔲 | 🔲 | ARL + c[n+a0] |
| Multi-program slot | StateShader (compiled XBE available) | 🔲 | 🔲 | CHEOPS_PROGRAM_START |

### Recommended Minimal Test Set (covers all instructions)

| Sample | XBE Path | VS Instructions Covered | PS Features | Status |
|--------|----------|------------------------|-------------|--------|
| Fur | `Fur\Fur.xbe` | MOV MUL ADD MAD DP3 DP4 DST MIN MAX RSQ RCP EXP LOG ARL | Fog, expand mapping | ✅ 8s OK |
| MatrixPaletteSkinning | `MatrixPaletteSkinning\MatrixPaletteSkinning.xbe` | MOV MUL MAD DP3 DP4 DPH MIN RSQ LIT ARL | Cubemap, bumpenvmap | ✅ 8s OK |
| ShadowBuffer | `ShadowBuffer\ShadowBuffer.xbe` | MOV MUL ADD DP3 DP4 SLT RCP ARL | Shadow compare, PROJECT2D | ✅ 8s OK |
| Dolphin | `Dolphin\Dolphin.xbe` | MOV MUL ADD MAD DP3 DP4 RCP RCC LIT | Fog, MUX, EF_PROD | ✅ 8s OK |
| Fog | `Fog\Fog.xbe` | MOV ADD MAD DP3 DP4 MAX RSQ RCC EXP | Vertex fog | ✅ 8s OK |
| Water | `Water\Water.xbe` | MOV MUL ADD MAD DP3 DP4 RCP RSQ SGE LIT | Bumpenvmap, dependent tex | ✅ 8s OK |
| StencilDepth | `StencilDepth\StencilDepth.xbe` | MOV MUL MAD DP3 DP4 DPH SLT RSQ RCC | Stencil operations | ✅ 8s OK |
| DisplacementMap | `DisplacementMap\DisplacementMap.xbe` | MOV MUL ADD MAD DP4 SGE ARL | Displacement | ✅ 8s OK |

**Coverage**: This set of 8 samples exercises all 17 MAC ops (MOV MUL ADD MAD DP3 DP4 DPH DST MIN MAX SLT SGE ARL) and all 7 ILU ops (RCP RCC RSQ EXP LOG LIT + ILU MOV via paired).

### Fixed-Function VS Features
| Feature | XDK Sample(s) | JIT/Interp Status | FF VS Status | Notes |
|---------|---------------|-------------------|--------------|-------|
| Vertex transform | Matrices, Vertices, Tutorials | N/A | 🔲 | CxbxFixedFunctionVertexShader.hlsl |
| Lighting (directional) | Lights, VSLights, PerPixelLighting | N/A | 🔲 | 8 lights max |
| Lighting (point) | Lights, Lensflare | N/A | 🔲 | Attenuation (DST+RCP) |
| Lighting (spot) | Lights | N/A | 🔲 | Cone + atten |
| Fog (vertex) | Fog, AlphaFog, Dolphin, FieldRender, Water | N/A | ⚠️ | Pre-existing flicker |
| Two-sided lighting | TwoSidedLighting, HeatShimmer, StencilBuffer | N/A | 🔲 | D3DRS_TWOSIDEDLIGHTING |
| Matrix palette skinning | MatrixPaletteSkinning | N/A | 🔲 | Blend weights + ARL |
| Point sprites | PointSprites, PaintEffect, VolumeSprites | N/A | ✅ | D3DRS_POINTSPRITEENABLE |
| Texgen (object-linear) | UserClipPlane? | N/A | ❌ | NOT IMPLEMENTED |
| Texgen (eye-linear) | UserClipPlane, MirrorClip | N/A | ❌ | NOT IMPLEMENTED |
| Texgen (sphere-map) | SphereMap | N/A | 🔲 | |
| Texgen (reflection) | CubeMap, FresnelReflect, EnvMapping | N/A | 🔲 | |
| Clip planes (NV2A_TX*) | UserClipPlane, MirrorClip | N/A | ❌ | Registers read but not computed in VS |

## JIT ↔ Interpreter Remaining Differences
| Area | Difference | Impact | Priority |
|------|-----------|--------|----------|
| Context writes | JIT blocks; Interpreter supports | Falls back gracefully | N/A |

Note: RCC sign preservation and LOG(0) bugs were fixed in commit 39fd1b808.
Both JIT and interpreter now share identical implementations via CxbxNV2AVshOps.hlsli.

## VS ↔ Xemu Differences
| Area | Cxbx | Xemu | Notes |
|------|------|------|-------|
| nv2a_mul 0×inf | Returns 0 | Returns 0 (via NaNToOne+sign) | Both correct; different implementation |
| nv2a_mul overflow | Clamps to ±FLT_MAX | No clamping | Cxbx more accurate (matches nv2a_vsh_cpu fix_inf_mult) |
| DP3/DP4 overflow | nv2a_mul per-component + sum clamp | Raw GLSL dot() | Cxbx more accurate |
| SLT -0 handling | -0 < +0 = true (bit fixup) | -0 < +0 = false (GLSL) | Cxbx matches nv2a_vsh_cpu nv2a_less_than |
| RCC clamp | [2^-64, 2^64] exact IEEE bits | clampAwayZeroInf (same range) | Both match |
| LOG .y mantissa | IEEE bit extraction | exp2/div decomposition | Both equivalent for normals |
| LIT kMax | 127.9961f | 128 - 1/256 = 127.996094 | Both match (float rounds same) |
| Instruction set | Complete (shared header) | Complete | Both cover all 20 opcodes |
| Fixed-function | Partial (no clip planes) | Full implementation | Cxbx needs FF texgen/clip |

## Known Issues & Commits

See [rendering_test_status.md](rendering_test_status.md) for the consolidated known issues, fixed issues, and commit tracking.

## Testing Procedure
1. Clear shader cache: `Remove-Item "$env:APPDATA\Cxbx-Reloaded\ShaderCache" -Recurse -Force`
2. Build: `cmake --build . --config Release --target cxbxr-emu`
3. Launch: `cxbxr-ldr.exe /load "<path>.xbe"`
4. Toggle JIT: In VertexShaderCache.cpp, set `g_bEnableVSJIT = false` to force interpreter fallback
5. Compare JIT vs interpreter output visually
