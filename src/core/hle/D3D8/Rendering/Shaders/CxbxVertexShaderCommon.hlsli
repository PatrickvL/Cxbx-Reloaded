#ifndef CXBX_VERTEX_SHADER_COMMON_HLSLI
#define CXBX_VERTEX_SHADER_COMMON_HLSLI

// Shared input layout: flat TEXCOORD array matching the NV2A model of
// 16 generic vertex attribute registers. Used by all vertex shaders
// under D3D11; under D3D9, FixedFunctionVertexShader overrides this
// with a semantic-based layout by defining CXBX_VS_CUSTOM_INPUT before
// including this header.
#ifdef CXBX_IA_BYPASS
// IA bypass mode: VS receives only SV_VertexID, all attributes are
// fetched from a ByteAddressBuffer via CxbxVertexFetch.hlsli
struct VS_INPUT
{
	uint vertexId : SV_VertexID;
};
#elif !defined(CXBX_VS_CUSTOM_INPUT)
struct VS_INPUT
{
	float4 v[16] : TEXCOORD;
};
#endif

// Output registers — declared identical to pixel shader input (see PS_INPUT)
struct VS_OUTPUT
{
#if defined(CXBX_USE_D3D11) || __HLSL_VERSION >= 4
	float4 oPos : SV_Position;  // Homogeneous clip space position (SM4.0+)
#else
	float4 oPos : POSITION;  // Homogeneous clip space position
#endif
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

#ifndef CXBX_USE_D3D11
// Default values for vertex registers, and whether to use them
// D3D9 path: init_v() lerps between vertex data and these defaults
uniform float4 vRegisterDefaultValues[16]  : register(c192);
#endif

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
