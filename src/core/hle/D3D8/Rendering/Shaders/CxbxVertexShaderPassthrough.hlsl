// Xbox HLSL pretransformed vertex shader

#include "CxbxVertexShaderCommon.hlsli"
#include "CxbxVertexFetch.hlsli"
#include "CxbxScreenspaceTransform.hlsli"

VS_OUTPUT main(const VS_INPUT xIn)
{
    // Input registers
    float4 v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

#include "CxbxVertexInputLoad.hlsli"

    // For passthrough, map output variables to their corresponding input registers
    float4 oPos = v0;
    float4 oD0 = v3;
    float4 oD1 = v4;
    float4 oFog = v5;
    float4 oPts = v6;
    float4 oB0 = v7;
    float4 oB1 = v8;
    float4 oT0 = v9;
    float4 oT1 = v10;
    float4 oT2 = v11;
    float4 oT3 = v12;

    // Copy variables to output struct
    VS_OUTPUT xOut;
#include "CxbxVertexOutputFooter.hlsli"

    return xOut;
}
