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

// Parameters for mapping the shader's fog output value to a fog factor
uniform float4 CxbxFogInfo : register(c218); // = CXBX_D3DVS_CONSTREG_FOGINFO

// Fog table formula — shared by all VS paths (programmable footer + FF DoFog).
// NV2A evaluates this per-vertex; the rasterizer interpolates the result;
// the PS clamps to [0,1] and blends with the fog color.
// fogMode: 0=NONE (vertex fog passthrough), 1=EXP, 2=EXP2, 3=LINEAR
//          5=EXP_ABS, 6=EXP2_ABS, 7=LINEAR_ABS (apply abs to computed factor)
// Bit 2 (value 4) is the _ABS flag, matching NV2A PGRAPH FOG_MODE encoding.
float CalculateFogFactor(int fogMode, float fogDensity, float fogStart, float fogEnd, float fogDepth)
{
    int baseMode = fogMode & 3;
    float f;
    if (baseMode == 1)       // EXP
        f = exp(-fogDepth * fogDensity);
    else if (baseMode == 2)  // EXP2
    {
        float fd = fogDepth * fogDensity;
        f = exp(-(fd * fd));
    }
    else if (baseMode == 3)  // LINEAR
        f = (fogEnd - fogDepth) / (fogEnd - fogStart);
    else
        return fogDepth;     // 0 = NONE (vertex fog passthrough)

    return (fogMode & 4) ? abs(f) : f;
}

// TEXCOORDINDEX remapping: xyzw = texcoord source index for stages 0-3.
// On NV2A, the texture unit routes interpolated texcoords based on
// D3DTSS_TEXCOORDINDEX. In D3D11 we apply this in the VS footer.
uniform float4 xboxTexCoordIndex : register(c219); // = CXBX_D3DVS_CONSTREG_TEXCOORDINDEX

#endif // CXBX_VERTEX_SHADER_COMMON_HLSLI
