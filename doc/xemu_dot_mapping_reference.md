# Xemu NV2A Dot Mapping Reference

Reference implementation from xemu for NV2A register combiner dot product mapping modes.

Source: `hw/xbox/nv2a/pgraph/psh_regs.h` and `hw/xbox/nv2a/pgraph/glsl/psh.c`

## Enum Definition

```c
enum PS_DOTMAPPING {
    PS_DOTMAPPING_ZERO_TO_ONE=         0x00L,
    PS_DOTMAPPING_MINUS1_TO_1_D3D=     0x01L,
    PS_DOTMAPPING_MINUS1_TO_1_GL=      0x02L,
    PS_DOTMAPPING_MINUS1_TO_1=         0x03L,
    PS_DOTMAPPING_HILO_1=              0x04L,
    PS_DOTMAPPING_HILO_HEMISPHERE_D3D= 0x05L,
    PS_DOTMAPPING_HILO_HEMISPHERE_GL=  0x06L,
    PS_DOTMAPPING_HILO_HEMISPHERE=     0x07L,
};
```

## Helper Functions

### sign1 — MINUS1_TO_1_D3D (Mode 1)
`(x*255.0 - 128.0) / 127.0` — D3D convention [0,1] → [-1,1]

### sign2 — MINUS1_TO_1_GL (Mode 2)
```
if (x*255 >= 128): (x*255 - 255.5) / 127.5
else:              (x*255 + 0.5)   / 127.5
```

### sign3 — MINUS1_TO_1 (Mode 3)
```
if (x*255 >= 128): (x*255 - 256.0) / 127.0
else:              (x*255)         / 127.0
```

## Dot Mapping Functions

| Mode | Name | Formula | Notes |
|------|------|---------|-------|
| 0 | ZERO_TO_ONE | `col.rgb` (passthrough) | No conversion |
| 1 | MINUS1_TO_1_D3D | `sign1(col.r/g/b)` | D3D signed |
| 2 | MINUS1_TO_1_GL | `sign2(col.r/g/b)` | GL signed |
| 3 | MINUS1_TO_1 | `sign3(col.r/g/b)` | Generic signed |
| 4 | HILO_1 | hi=(A<<8\|R)/0xFFFF, lo=(G<<8\|B)/0xFFFF, z=1.0 | 16-bit pair |
| 5 | HILO_HEMISPHERE_D3D | **UNIMPLEMENTED** in xemu | Implemented in Cxbx (CxbxPixelShaderFunctions.hlsli) |
| 6 | HILO_HEMISPHERE_GL | **UNIMPLEMENTED** in xemu | Implemented in Cxbx (CxbxPixelShaderFunctions.hlsli) |
| 7 | HILO_HEMISPHERE | **UNIMPLEMENTED** in xemu | Implemented in Cxbx (CxbxPixelShaderFunctions.hlsli) |

Modes 5-7 use signed reconstruction + hemisphere formula: `z = sqrt(max(0, 1 - Hs² - Ls²))`

Used with PS_TEXTUREMODES_DOT_* texture modes (DOT_ST, DOT_ZW, DOT_RFLCT_DIFF, DOT_RFLCT_SPEC, etc.)

Cxbx implementation: `src/core/hle/D3D8/Rendering/Shaders/CxbxPixelShaderFunctions.hlsli` (sign1/sign2/sign3 + dotmap functions)
