// CxbxPSAuxFromPGRAPH.hlsli — Derive PS auxiliary state directly from PGRAPH SRV
//
// These functions read NV2A registers from the mirror buffer SRV (g_PGRegs at t12)
// and compute values that were previously uploaded in the PSAuxCBLayout cbuffer.
// This eliminates per-draw CPU-side computation for these fields.
//
// Requires: CxbxNV2APixelShaderConstants.hlsli and CxbxPGRAPHRegs.hlsli to be included first.

#ifndef CXBX_PSAUX_FROM_PGRAPH_HLSLI
#define CXBX_PSAUX_FROM_PGRAPH_HLSLI

// ============================================================
// Additional PGRAPH register offsets needed by these helpers
// (supplement those already in CxbxPGRAPHRegs.hlsli)
// ============================================================
#define NV_PGRAPH_TEXCTL0_0             0x19CC
#define NV_PGRAPH_TEXCTL0_1             0x19D0
#define NV_PGRAPH_TEXCTL0_2             0x19D4
#define NV_PGRAPH_TEXCTL0_3             0x19D8
#define NV_PGRAPH_TEXFMT0               0x1A04
#define NV_PGRAPH_TEXFMT1               0x1A08
#define NV_PGRAPH_TEXFMT2               0x1A0C
#define NV_PGRAPH_TEXFMT3               0x1A10
#define NV_PGRAPH_TEXFILTER0            0x19F4
#define NV_PGRAPH_TEXFILTER1            0x19F8
#define NV_PGRAPH_TEXFILTER2            0x19FC
#define NV_PGRAPH_TEXFILTER3            0x1A00
#define NV_PGRAPH_COLORKEYCOLOR0        0x1870
#define NV_PGRAPH_COLORKEYCOLOR1        0x1874
#define NV_PGRAPH_COLORKEYCOLOR2        0x1878
#define NV_PGRAPH_COLORKEYCOLOR3        0x187C
#define NV_PGRAPH_CONTROL_3             0x1958
#define NV_PGRAPH_SETUPRASTER           0x1990
#define NV_PGRAPH_SURFACEFORMAT         0x0714

// Bitmasks
#define TEXCTL0_ENABLE                  0x40000000u
#define TEXCTL0_ALPHAKILLEN             0x00000004u
#define TEXCTL0_COLORKEYMODE            0x00000003u
#define TEXFMT0_CUBEMAPENABLE           0x00000004u
#define TEXFMT0_DIMENSIONALITY          0x000000C0u
#define TEXFMT0_COLOR                   0x00007F00u
#define TEXFMT0_COLOR_SHIFT              8u
#define CONTROL_3_FOGENABLE             0x00000100u
#define CONTROL_3_FOG_MODE              0x00070000u
#define CONTROL_3_FOG_MODE_SHIFT        16u
#define CSV0_C_SPECULAR_ENABLE          0x00010000u
#define CSV0_C_TWO_SIDE_ENABLE          0x20000000u
#define SETUPRASTER_FRONTFACE           0x00800000u
#define SURFACEFORMAT_ZETA              0x000000F0u

// NV2A texture format codes for depth buffers
#define NV2A_TEX_COLOR_D_Z24_S8_FIXED   0x2Au
#define NV2A_TEX_COLOR_D_Z24_S8_FLOAT   0x2Bu
#define NV2A_TEX_COLOR_D_Z16_FIXED      0x2Cu
#define NV2A_TEX_COLOR_D_Z16_FLOAT      0x2Du
#define NV2A_TEX_COLOR_LU_D_Z24_FIXED   0x2Eu
#define NV2A_TEX_COLOR_LU_D_Z24_FLOAT   0x2Fu
#define NV2A_TEX_COLOR_LU_D_Z16_FIXED   0x30u
#define NV2A_TEX_COLOR_LU_D_Z16_FLOAT   0x31u

// PS_TEXTUREMODES constants for DeriveAdjustedPSTextureModes().
// Use _PGAUX suffix to avoid collision with CxbxNV2APixelShaderConstants.hlsli
// which defines the same names as `static const uint`.
#define PGAUX_TEXMODE_PROJECT2D    1u
#define PGAUX_TEXMODE_CUBEMAP      3u
#define PGAUX_TEXMODE_DOT_STR_3D   13u
#define PGAUX_TEXMODE_DOT_STR_CUBE 14u

// PS_REGISTER / PS_CHANNEL constants for final combiner synthesis.
// Use PGAUX_ prefix to avoid collision with CxbxNV2APixelShaderConstants.hlsli.
// Values MUST match the NV2A register encoding.
#define PGAUX_REGISTER_ZERO     0x00u
#define PGAUX_REGISTER_FOG      0x03u
#define PGAUX_REGISTER_R0       0x0Cu
#define PGAUX_REGISTER_V1       0x05u
#define PGAUX_CHANNEL_ALPHA     0x10u

// Surface format zeta values
#define ZETA_Z16                        1u
#define ZETA_Z24S8                      2u

// ============================================================
// Helper: read TEXCTL0 for a given stage (0-3)
// ============================================================
uint PG_TEXCTL0(uint stage)
{
    // Registers are at consecutive 4-byte offsets: 0x19CC, 0x19D0, 0x19D4, 0x19D8
    return PG_UINT(NV_PGRAPH_TEXCTL0_0 + stage * 4u);
}

// Helper: read TEXFMT for a given stage (0-3)
uint PG_TEXFMT(uint stage)
{
    return PG_UINT(NV_PGRAPH_TEXFMT0 + stage * 4u);
}

// Helper: read COLORKEYCOLOR for a given stage (0-3)
uint PG_COLORKEYCOLOR(uint stage)
{
    return PG_UINT(NV_PGRAPH_COLORKEYCOLOR0 + stage * 4u);
}

// ============================================================
// PSTextureModes — adjusted for cubemap/volume
// ============================================================
uint DeriveAdjustedPSTextureModes()
{
    uint modes = PG_UINT(NV_PGRAPH_SHADERPROG);

    [unroll] for (uint i = 0; i < 4; i++) {
        uint mode = (modes >> (i * NV_PGRAPH_SHADERPROG_STAGE_BITS)) & PS_TEXTUREMODES_MASK;
        uint clearMask = ~(PS_TEXTUREMODES_MASK << (i * NV_PGRAPH_SHADERPROG_STAGE_BITS));

        uint texCtl = PG_TEXCTL0(i);
        if (texCtl & TEXCTL0_ENABLE) {
            uint texFmt = PG_TEXFMT(i);
            bool isCubemap = (texFmt & TEXFMT0_CUBEMAPENABLE) != 0;
            if (isCubemap) {
                if (mode == PGAUX_TEXMODE_PROJECT2D)
                    modes = (modes & clearMask) | (PGAUX_TEXMODE_CUBEMAP << (i * NV_PGRAPH_SHADERPROG_STAGE_BITS));
                else if (mode == PGAUX_TEXMODE_DOT_STR_3D)
                    modes = (modes & clearMask) | (PGAUX_TEXMODE_DOT_STR_CUBE << (i * NV_PGRAPH_SHADERPROG_STAGE_BITS));
            }
        }
    }
    return modes;
}

// ============================================================
// PSFinalCombinerInputs — synthesize default if not explicitly set
// ============================================================
void DeriveFinalCombinerInputs(out uint outABCD, out uint outEFG)
{
    uint fcABCD = PG_UINT(NV_PGRAPH_COMBINESPECFOG0);
    uint fcEFG  = PG_UINT(NV_PGRAPH_COMBINESPECFOG1);

    if (fcABCD == 0u && fcEFG == 0u) {
        uint ctrl3 = PG_UINT(NV_PGRAPH_CONTROL_3);
        uint csv0c = PG_UINT(NV_PGRAPH_CSV0_C);
        bool fogEnable = (ctrl3 & CONTROL_3_FOGENABLE) != 0;
        bool specularEnable = (csv0c & CSV0_C_SPECULAR_ENABLE) != 0;

        uint regA = PGAUX_REGISTER_FOG | PGAUX_CHANNEL_ALPHA;
        uint regB = PGAUX_REGISTER_R0;
        uint regC = fogEnable ? PGAUX_REGISTER_FOG : PGAUX_REGISTER_R0;
        uint regD = specularEnable ? PGAUX_REGISTER_V1 : PGAUX_REGISTER_ZERO;
        fcABCD = (regA << 24u) | (regB << 16u) | (regC << 8u) | regD;

        uint regE = PGAUX_REGISTER_ZERO;
        uint regF = PGAUX_REGISTER_ZERO;
        uint regG = PGAUX_REGISTER_R0 | PGAUX_CHANNEL_ALPHA;
        fcEFG = (regE << 24u) | (regF << 16u) | (regG << 8u);
    }

    outABCD = fcABCD;
    outEFG  = fcEFG;
}

// ============================================================
// ColorKeyOp — per-stage color key mode from TEXCTL0
// ============================================================
float4 DeriveColorKeyOp()
{
    return float4(
        float(PG_TEXCTL0(0) & TEXCTL0_COLORKEYMODE),
        float(PG_TEXCTL0(1) & TEXCTL0_COLORKEYMODE),
        float(PG_TEXCTL0(2) & TEXCTL0_COLORKEYMODE),
        float(PG_TEXCTL0(3) & TEXCTL0_COLORKEYMODE)
    );
}

// ============================================================
// ColorKeyColor — per-stage, unpacked from ABGR uint32
// ============================================================
float4 DeriveColorKeyColor(uint stage)
{
    return UnpackABGR(PG_COLORKEYCOLOR(stage));
}

// ============================================================
// AlphaKill — per-stage ALPHAKILLEN bit
// ============================================================
float4 DeriveAlphaKill()
{
    return float4(
        (PG_TEXCTL0(0) & TEXCTL0_ALPHAKILLEN) ? 1.0f : 0.0f,
        (PG_TEXCTL0(1) & TEXCTL0_ALPHAKILLEN) ? 1.0f : 0.0f,
        (PG_TEXCTL0(2) & TEXCTL0_ALPHAKILLEN) ? 1.0f : 0.0f,
        (PG_TEXCTL0(3) & TEXCTL0_ALPHAKILLEN) ? 1.0f : 0.0f
    );
}

// ============================================================
// FogEnable — from CONTROL_3
// ============================================================
uint DeriveFogEnable()
{
    return (PG_UINT(NV_PGRAPH_CONTROL_3) & CONTROL_3_FOGENABLE) ? 1u : 0u;
}

// ============================================================
// FogMode — from CONTROL_3 FOG_MODE field
// ============================================================
uint DeriveFogMode()
{
    return (PG_UINT(NV_PGRAPH_CONTROL_3) & CONTROL_3_FOG_MODE) >> CONTROL_3_FOG_MODE_SHIFT;
}

// ============================================================
// FrontFaceInfo — from CSV0_C (two-sided) + SETUPRASTER (front face)
// ============================================================
float DeriveFrontFaceFactor()
{
    uint csv0c = PG_UINT(NV_PGRAPH_CSV0_C);
    bool twoSided = (csv0c & CSV0_C_TWO_SIDE_ENABLE) != 0;
    if (!twoSided) return 0.0f;

    uint setup = PG_UINT(NV_PGRAPH_SETUPRASTER);
    bool ccwFront = (setup & SETUPRASTER_FRONTFACE) != 0;
    return ccwFront ? -1.0f : 1.0f;
}

// ============================================================
// ShadowCompare — per-stage: is a depth texture bound?
// Check TEXFMT0 COLOR field for NV2A depth format codes.
// ============================================================
bool IsNV2ADepthFormat(uint colorCode)
{
    return (colorCode >= NV2A_TEX_COLOR_D_Z24_S8_FIXED && colorCode <= NV2A_TEX_COLOR_D_Z16_FLOAT)
        || (colorCode >= NV2A_TEX_COLOR_LU_D_Z24_FIXED && colorCode <= NV2A_TEX_COLOR_LU_D_Z16_FLOAT);
}

float4 DeriveShadowCompare()
{
    float4 sc = float4(0, 0, 0, 0);
    [unroll] for (uint i = 0; i < 4; i++) {
        uint texCtl = PG_TEXCTL0(i);
        if (texCtl & TEXCTL0_ENABLE) {
            uint texFmt = PG_TEXFMT(i);
            uint colorCode = (texFmt & TEXFMT0_COLOR) >> TEXFMT0_COLOR_SHIFT;
            if (IsNV2ADepthFormat(colorCode))
                sc[i] = 1.0f;
        }
    }
    return sc;
}

#endif // CXBX_PSAUX_FROM_PGRAPH_HLSLI
