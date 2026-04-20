// Xbox HLSL pretransformed vertex shader

#include "CxbxVertexShaderCommon.hlsli"
#include "CxbxVertexFetch.hlsli"
#include "CxbxScreenspaceTransform.hlsli"

VS_OUTPUT main(const VS_INPUT xIn)
{
    // Fetch all 16 vertex input registers from ByteAddressBuffer
    float4 v[16];
    FetchAllAttributes(ResolveVertexIndex(xIn.vertexId), v);

    // For passthrough, map output variables to their corresponding input registers
    float4 oPos = v[0];
    float4 oD0 = v[3];
    float4 oD1 = v[4];
    float4 oFog = v[5];
    float4 oPts = v[6];
    float4 oB0 = v[7];
    float4 oB1 = v[8];
    float4 oT0 = v[9];
    float4 oT1 = v[10];
    float4 oT2 = v[11];
    float4 oT3 = v[12];

    // Copy variables to output struct
    VS_OUTPUT xOut;
#include "CxbxVertexOutputFooter.hlsli"

    return xOut;
}
