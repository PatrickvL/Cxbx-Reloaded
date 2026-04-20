// CxbxPixelShaderInput.hlsli — shared pixel shader input struct
//
// Declares PS_INPUT, which must match the vertex shader output (VS_OUTPUT).
// Included by all pixel shader variants: compiled PS template,
// fixed-function PS, and the RC interpreter ubershader.

#ifndef CXBX_PIXEL_SHADER_INPUT_HLSLI
#define CXBX_PIXEL_SHADER_INPUT_HLSLI

struct PS_INPUT // Declared identical to vertex shader output (see VS_OUTPUT)
{
	float4 iPos : SV_Position; // Screen space position (SM4.0+ requires float4 xyzw)
	float4 iD0  : COLOR0;      // Front-facing primary (diffuse) vertex color (clamped to 0..1)
	float4 iD1  : COLOR1;      // Front-facing secondary (specular) vertex color (clamped to 0..1)
	float  iFog : FOG;
	float  iPts : PSIZE;
	float4 iB0  : TEXCOORD4;   // Back-facing primary (diffuse) vertex color (clamped to 0..1)
	float4 iB1  : TEXCOORD5;   // Back-facing secondary (specular) vertex color (clamped to 0..1)
	float4 iT0  : TEXCOORD0;   // Texture Coord 0
	float4 iT1  : TEXCOORD1;   // Texture Coord 1
	float4 iT2  : TEXCOORD2;   // Texture Coord 2
	float4 iT3  : TEXCOORD3;   // Texture Coord 3
	bool   iFF  : SV_IsFrontFace; // SM4.0+: bool type required
};

#endif // CXBX_PIXEL_SHADER_INPUT_HLSLI
