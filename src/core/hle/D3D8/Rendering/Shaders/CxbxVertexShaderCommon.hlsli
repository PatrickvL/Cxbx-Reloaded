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

// Per-stage texture coordinate scale factors
uniform float4 xboxTextureScale[4] : register(c214);

// Parameters for mapping the shader's fog output value to a fog factor
uniform float4 CxbxFogInfo : register(c218); // = CXBX_D3DVS_CONSTREG_FOGINFO

// TEXCOORDINDEX remapping: xyzw = texcoord source index for stages 0-3.
// On NV2A, the texture unit routes interpolated texcoords based on
// D3DTSS_TEXCOORDINDEX. In D3D11 we apply this in the VS footer.
uniform float4 xboxTexCoordIndex : register(c219); // = CXBX_D3DVS_CONSTREG_TEXCOORDINDEX

#endif // CXBX_VERTEX_SHADER_COMMON_HLSLI
