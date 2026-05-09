#ifndef CXBX_VERTEX_SHADER_COMMON_HLSLI
#define CXBX_VERTEX_SHADER_COMMON_HLSLI

// Shared vertex shader input: VS receives only SV_VertexID; all 16 NV2A
// vertex attributes are fetched from a ByteAddressBuffer in the shader
// (see CxbxVertexFetch.hlsli).
struct VS_INPUT
{
	uint vertexId : SV_VertexID;
};

// Output registers — declared identical to pixel shader input (see PS_INPUT)
struct VS_OUTPUT
{
	float4 oPos : SV_Position;  // Homogeneous clip space position (SM4.0+)
	float4 oD0  : COLOR0;    // Primary color (front-facing)
	float4 oD1  : COLOR1;    // Secondary color (front-facing)
	float  oFog : FOG;       // Fog coordinate
	float  oPts : PSIZE;     // Point size
	float4 oB0  : TEXCOORD4; // Back-facing primary color
	float4 oB1  : TEXCOORD5; // Back-facing secondary color
	float4 oT0  : TEXCOORD0; // Texture coordinate set 0
	float4 oT1  : TEXCOORD1; // Texture coordinate set 1
	float4 oT2  : TEXCOORD2; // Texture coordinate set 2
	float4 oT3  : TEXCOORD3; // Texture coordinate set 3
};

// Whether each vertex register is present in the vertex declaration
uniform float4 vRegisterDefaultFlagsPacked[4]  : register(c208);

// Per-stage reciprocal texture coordinate scale factors (1/scale)
// Uploaded from C++ as rcp(scale); multiply is cheaper than divide per vertex.
uniform float4 xboxTextureScaleRcp[4] : register(c214);

// Parameters for the NV2A fog computation.
// CxbxFogInfo: x=fogMode (PGRAPH CONTROL_3 FOG_MODE), y=fogParam0, z=fogParam1, w=unused
// fogParam0/1 are pre-baked coefficients from NV_PGRAPH_FOGPARAM0/1 (set by Xbox D3D
// via NV097_SET_FOG_PARAMS). Their meaning depends on the fog mode:
//   LINEAR: fogParam0 = 1 - end/(end-start), fogParam1 = -1/(end-start)
//   EXP:    fogParam0 = 1.5, fogParam1 = -density/(2*ln(256))
//   EXP2:   fogParam0 = 1.5, fogParam1 = -density/(2*sqrt(ln(256)))
uniform float4 CxbxFogInfo : register(c218); // = CXBX_D3DVS_CONSTREG_FOGINFO

// NV2A-native fog factor computation using FOGPARAM0/1 pre-baked coefficients.
// Matches xemu's implementation exactly.
// fogMode: NV2A PGRAPH FOG_MODE (0=LINEAR, 1=EXP, 3=EXP2, 4=LINEAR_ABS, 5=EXP_ABS, 7=EXP2_ABS)
// fogParam0/1: Pre-baked coefficients from NV_PGRAPH_FOGPARAM0/1
// fogDistance: The fog coordinate (oFog.x from VS, or computed from position for FF)
float CalculateFogFactor(int fogMode, float fogParam0, float fogParam1, float fogDistance)
{
    int baseMode = fogMode & 3;

    // Infinite fog distance: LINEAR = fully fogged, EXP/EXP2 = fully visible.
    if (isinf(fogDistance)) return (baseMode == 0) ? 0.0 : 1.0;

    // Shared bias: LINEAR uses -1.0 (= -1.5 + 0.5), EXP/EXP2 use -1.5.
    float fogFactor    = fogParam0 - 1.5;
    float fogDistParam = fogDistance * fogParam1;

    if (baseMode == 0) {
        // LINEAR / LINEAR_ABS
        fogFactor += 0.5 + fogDistParam;
    } else if (baseMode == 1) {
        // EXP / EXP_ABS
        fogFactor += exp2(fogDistParam * 16.0);
    } else {
        // EXP2 / EXP2_ABS (baseMode == 3)
        fogFactor += exp2(fogDistParam * fogDistParam * -32.0);
    }

    // _ABS variants: bit 2 of fogMode.
    if ((fogMode & 4) != 0)
        fogFactor = abs(fogFactor);

    // exp2 can overflow to Inf then NaN downstream; LINEAR cannot (no transcendentals,
    // and infinite inputs are already handled above). Hardware clamps here too.
    // Bit-pattern check bypasses compiler folding away isinf/isnan on some D3D drivers.
    uint bits = asuint(fogFactor);
    bool nan = ((bits & 0x7F800000u) == 0x7F800000u) && ((bits & 0x007FFFFFu) != 0u);
    return nan ? 1.0 : clamp(fogFactor, -3.4e+38, 3.4e+38);
}

// TEXCOORDINDEX remapping: xyzw = texcoord source index for stages 0-3.
// On NV2A, the texture unit routes interpolated texcoords based on
// D3DTSS_TEXCOORDINDEX. In D3D11 we apply this in the VS footer.
uniform float4 xboxTexCoordIndex : register(c219); // = CXBX_D3DVS_CONSTREG_TEXCOORDINDEX

#endif // CXBX_VERTEX_SHADER_COMMON_HLSLI
