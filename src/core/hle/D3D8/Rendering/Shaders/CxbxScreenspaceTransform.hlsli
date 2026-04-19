#ifndef CXBX_SCREENSPACE_TRANSFORM_HLSLI
#define CXBX_SCREENSPACE_TRANSFORM_HLSLI

uniform float4 xboxScreenspaceScale : register(c212);
uniform float4 xboxScreenspaceOffset : register(c213);

float4 reverseScreenspaceTransform(float4 oPos)
{
	// On Xbox, oPos should contain the vertex position in screenspace
	// We need to reverse this transformation
	// Conventionally, each Xbox Vertex Shader includes instructions like this
	// mul oPos.xyz, r12, c-38
	// +rcc r1.x, r12.w
	// mad oPos.xyz, r12, r1.x, c-37
	// where c-37 and c-38 are reserved transform values

	// oPos.w and xboxViewportScale.z might be VERY big when a D24 depth buffer is used
	// and multiplying oPos.xyz by oPos.w may cause precision issues.
	// Test case: Burnout 3

	// Reverse screenspace offset
	oPos -= xboxScreenspaceOffset;
	// Reverse screenspace scale
	oPos /= xboxScreenspaceScale;

	// Ensure w is nonzero
	if(oPos.w == 0) oPos.w = 1;
	// Reverse perspective divide
	oPos.xyz *= oPos.w;

	return oPos;
}

#endif // CXBX_SCREENSPACE_TRANSFORM_HLSLI
