// CxbxPGRAPHRegs.hlsli — NV2A PGRAPH register offsets and raw-buffer accessors
//
// HLSL-only header.  Provides register byte offsets matching nv2a_regs.h,
// plus typed accessor functions for reading from the StructuredBuffer<uint>
// bound by the C++ backend (pg->regs[2048]).
//
// Usage in HLSL:
//   uint  val  = PG_UINT(NV_PGRAPH_COMBINECTL);
//   float fval = PG_FLOAT(NV_PGRAPH_BUMPSCALE1);
//
// These are the subset of NV_PGRAPH_* offsets actually used by the RC
// interpreter.  Add more as needed for future shaders.

#ifndef CXBX_PGRAPH_REGS_HLSLI
#define CXBX_PGRAPH_REGS_HLSLI

// ============================================================
// The combined GPU memory buffer — 64 MiB RAM + appended PGRAPH (8 KB).
// Bound as a ByteAddressBuffer at t12 for the PS/VS stages that need
// PGRAPH register access. This is the same physical buffer as t0
// (g_VtxData) but bound separately so the PS stage can access it
// without conflicting with texture slots at t0-t3.
// ============================================================
ByteAddressBuffer g_PGRegs : register(t12);

// PGRAPH block sits at offset 0x04000000 within the combined buffer
#define GPU_PGRAPH_BASE 0x04000000u

// --- Register accessor helpers (byte offset → raw load) ---
uint PG_UINT(uint byteOff)   { return g_PGRegs.Load(GPU_PGRAPH_BASE + byteOff); }
float PG_FLOAT(uint byteOff) { return asfloat(g_PGRegs.Load(GPU_PGRAPH_BASE + byteOff)); }

// ============================================================
// NV_PGRAPH register byte offsets (from nv2a_regs.h)
// Only the subset used by the RC interpreter.
// ============================================================

// Bump environment matrix (stages 1-3, stride = 0x0C per component group)
#define NV_PGRAPH_BUMPMAT00                 0x181C
#define NV_PGRAPH_BUMPMAT01                 0x1828
#define NV_PGRAPH_BUMPMAT10                 0x1834
#define NV_PGRAPH_BUMPMAT11                 0x1840
#define NV_PGRAPH_BUMPOFFSET1               0x184C
#define NV_PGRAPH_BUMPSCALE1                0x1858

// Combiner factor constants (8 stages, stride = 4 bytes per stage)
#define NV_PGRAPH_COMBINEFACTOR0            0x1880
#define NV_PGRAPH_COMBINEFACTOR1            0x18A0

// Combiner alpha inputs/outputs (8 stages, stride = 4)
#define NV_PGRAPH_COMBINEALPHAI0            0x18C0
#define NV_PGRAPH_COMBINEALPHAO0            0x18E0

// Combiner color (RGB) inputs/outputs (8 stages, stride = 4)
#define NV_PGRAPH_COMBINECOLORI0            0x1900
#define NV_PGRAPH_COMBINECOLORO0            0x1920

// Combiner control + specular/fog
#define NV_PGRAPH_COMBINECTL                0x1940
#define NV_PGRAPH_COMBINESPECFOG0           0x1944
#define NV_PGRAPH_COMBINESPECFOG1           0x1948

// Control register 0 (alpha test, depth, etc.)
#define NV_PGRAPH_CONTROL_0                 0x194C
#define NV_PGRAPH_CONTROL_0_ALPHAREF        0x000000FF
#define NV_PGRAPH_CONTROL_0_ALPHAFUNC       0x00000F00
#define NV_PGRAPH_CONTROL_0_ALPHAFUNC_SHIFT 8
#define NV_PGRAPH_CONTROL_0_ALPHATESTENABLE 0x00001000

// Fog color (ABGR packed)
#define NV_PGRAPH_FOGCOLOR                  0x1980

// Shader control registers
#define NV_PGRAPH_SHADERCLIPMODE            0x1994
#define NV_PGRAPH_SHADERCLIPMODE_STAGE_BITS 4       // 4 bits (RSTQ) per stage
#define NV_PGRAPH_SHADERCLIPMODE_STAGE_MASK 0xFu

#define NV_PGRAPH_SHADERCTL                 0x1998
// PSDotMapping: bits [11:0], 3-bit field per stage (stages 1-3), stride 4 bits
#define NV_PGRAPH_SHADERCTL_DOTMAP_STRIDE   4
#define NV_PGRAPH_SHADERCTL_DOTMAP_MASK     0x7u
// PSInputTexture: bits [12:27], source-stage config per texture stage
#define NV_PGRAPH_SHADERCTL_PST2_SHIFT      16      // stage 2: 1-bit field
#define NV_PGRAPH_SHADERCTL_PST2_MASK       0x1u
#define NV_PGRAPH_SHADERCTL_PST3_SHIFT      20      // stage 3: 2-bit field
#define NV_PGRAPH_SHADERCTL_PST3_MASK       0x3u

#define NV_PGRAPH_SHADERPROG                0x199C
// PSTextureModes: 5-bit field per stage, packed sequentially
#define NV_PGRAPH_SHADERPROG_STAGE_BITS     5

// Shadow mapping control
#define NV_PGRAPH_SHADOWCTL                 0x19A4
#define NV_PGRAPH_SHADOWCTL_SHADOW_ZFUNC    0x00000007

// Specular fog factors (final combiner constants, ABGR packed)
#define NV_PGRAPH_SPECFOGFACTOR0            0x19AC
#define NV_PGRAPH_SPECFOGFACTOR1            0x19B0

// Vertex shader control (CHEOPS vertex processor)
#define NV_PGRAPH_CSV0_C                    0x0FB8
#define NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START_SHIFT 8
#define NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START_MASK  0xFF

// ============================================================
// Color unpacking: ABGR uint32 → float4 RGBA [0..1]
//
// Legacy unpacker kept for any PGRAPH registers that happen to
// store colors in ABGR byte order.  Most NV2A color registers
// (combiner factors, specfog factors, fog color) are ARGB —
// use PG_COLOR_ARGB / UnpackARGB for those.
//   bits  0-7  = R
//   bits  8-15 = G
//   bits 16-23 = B
//   bits 24-31 = A
// ============================================================
float4 UnpackABGR(uint c)
{
    return float4(
        float( c        & 0xFFu) / 255.0f,   // R
        float((c >>  8) & 0xFFu) / 255.0f,   // G
        float((c >> 16) & 0xFFu) / 255.0f,   // B
        float((c >> 24) & 0xFFu) / 255.0f    // A
    );
}

// ============================================================
// Color unpacking: ARGB uint32 → float4 RGBA [0..1]
//
// NV_PGRAPH_FOGCOLOR is stored in ARGB byte order (the PGRAPH
// SET_FOG_COLOR handler re-packs the ABGR method parameter):
//   bits  0-7  = B
//   bits  8-15 = G
//   bits 16-23 = R
//   bits 24-31 = A
// ============================================================
float4 UnpackARGB(uint c)
{
    return float4(
        float((c >> 16) & 0xFFu) / 255.0f,   // R
        float((c >>  8) & 0xFFu) / 255.0f,   // G
        float( c        & 0xFFu) / 255.0f,   // B
        float((c >> 24) & 0xFFu) / 255.0f    // A
    );
}

// ============================================================
// Convenience: read a color register and unpack ABGR → float4
// ============================================================
float4 PG_COLOR(uint byteOff) { return UnpackABGR(PG_UINT(byteOff)); }

// Convenience: read FOGCOLOR (stored as ARGB) and unpack → float4
float4 PG_COLOR_ARGB(uint byteOff) { return UnpackARGB(PG_UINT(byteOff)); }

#endif // CXBX_PGRAPH_REGS_HLSLI
