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
// The raw PGRAPH regs[] buffer — 2048 uint32 elements (8 KB).
// Bound as a StructuredBuffer<uint> at t12 (avoids t0-t11 texture slots).
// Shared by both PS (register combiner interpreter) and VS (vertex
// shader interpreter) stages — the same GPU buffer is bound to both.
// ============================================================
StructuredBuffer<uint> g_PGRegs : register(t12);

// --- Register index helper (byte offset → array index) ---
uint PG_UINT(uint byteOff)   { return g_PGRegs[byteOff >> 2]; }
float PG_FLOAT(uint byteOff) { return asfloat(g_PGRegs[byteOff >> 2]); }

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
#define NV_PGRAPH_CONTROL_0_ALPHATESTENABLE 0x00001000

// Fog color (ABGR packed)
#define NV_PGRAPH_FOGCOLOR                  0x1980

// Shader control registers
#define NV_PGRAPH_SHADERCLIPMODE            0x1994
#define NV_PGRAPH_SHADERCTL                 0x1998
#define NV_PGRAPH_SHADERPROG                0x199C

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
// NV2A stores color registers in ABGR byte order:
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
// Convenience: read a color register and unpack ABGR → float4
// ============================================================
float4 PG_COLOR(uint byteOff) { return UnpackABGR(PG_UINT(byteOff)); }

#endif // CXBX_PGRAPH_REGS_HLSLI
