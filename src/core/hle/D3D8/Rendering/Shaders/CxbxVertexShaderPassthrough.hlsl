// Xbox HLSL pretransformed vertex shader

#include "CxbxVertexShaderCommon.hlsli"
#include "CxbxScreenspaceTransform.hlsli"

// TEXCOORDINDEX remapping: xyzw = texcoord source index for stages 0-3
// On NV2A, the texture unit applies TEXCOORDINDEX after VS output interpolation.
// In D3D11, we must do this in the VS since there's no hardware texcoord routing.
uniform float4 xboxTexCoordIndex : register(c219); // = CXBX_D3DVS_CONSTREG_TEXCOORDINDEX

VS_OUTPUT main(const VS_INPUT xIn)
{
    // Input registers
    float4 v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

#ifdef CXBX_USE_D3D11
    // D3D11: The input assembler delivers correct values for all 16 attributes.
    // Streamed attributes come from their real vertex buffer slots.
    // Non-streamed attributes come from the zero-stride defaults buffer (slot 16)
    // which holds the NV2A's sticky inline_value[] registers.
    #define init_v(i) v##i = xIn.v[i];
#else
    // D3D9: Lerp between vertex data and constant buffer defaults
    float vRegisterDefaultFlags[16] = (float[16])vRegisterDefaultFlagsPacked;
    #define init_v(i) v##i = lerp(xIn.v[i], vRegisterDefaultValues[i], vRegisterDefaultFlags[i]);
#endif
    // Note : unroll manually instead of for-loop, because of the ## concatenation
    init_v( 0); init_v( 1); init_v( 2); init_v( 3);
    init_v( 4); init_v( 5); init_v( 6); init_v( 7);
    init_v( 8); init_v( 9); init_v(10); init_v(11);
    init_v(12); init_v(13); init_v(14); init_v(15);

    // For passthrough, map output variables to their corresponding input registers
    float4 oPos = v0;
    float4 oD0 = v3;
    float4 oD1 = v4;
    float4 oFog = v5;
    float4 oPts = v6;
    float4 oB0 = v7;
    float4 oB1 = v8;

    // Apply TEXCOORDINDEX remapping: on NV2A, the texture unit routes interpolated
    // texcoords to texture stages based on D3DTSS_TEXCOORDINDEX. In D3D11, we must
    // do this in the VS. xboxTexCoordIndex.xyzw holds the source texcoord set index
    // (0-3) for each texture stage.
    float4 texcoordSets[4] = { v9, v10, v11, v12 };
    float4 oT0 = texcoordSets[(int)xboxTexCoordIndex.x];
    float4 oT1 = texcoordSets[(int)xboxTexCoordIndex.y];
    float4 oT2 = texcoordSets[(int)xboxTexCoordIndex.z];
    float4 oT3 = texcoordSets[(int)xboxTexCoordIndex.w];

    // Copy variables to output struct
    VS_OUTPUT xOut;

    xOut.oPos = reverseScreenspaceTransform(oPos);
    xOut.oD0 = saturate(oD0);
    xOut.oD1 = saturate(oD1);
    xOut.oFog = oFog.x; // Note : Xbox clamps fog in pixel shader
    xOut.oPts = oPts.x;
    xOut.oB0 = saturate(oB0);
    xOut.oB1 = saturate(oB1);
    // Scale textures (TODO: or should we apply this to the input register values?)
    xOut.oT0 = oT0 / xboxTextureScale[0];
    xOut.oT1 = oT1 / xboxTextureScale[1];
    xOut.oT2 = oT2 / xboxTextureScale[2];
    xOut.oT3 = oT3 / xboxTextureScale[3];

    return xOut;
}
