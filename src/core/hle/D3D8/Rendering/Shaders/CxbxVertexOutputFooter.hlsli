// CxbxVertexOutputFooter.hlsli — Common VS output processing
//
// Requires in scope: VS_OUTPUT xOut (declared), float4 oPos/oD0/oD1/oB0/oB1/oT0-oT3/oFog/oPts.
// Applies: screen→clip transform, color saturation, TEXCOORDINDEX remapping, texture scale.

    xOut.oPos = reverseScreenspaceTransform(oPos);
    xOut.oD0 = saturate(oD0);
    xOut.oD1 = saturate(oD1);
    xOut.oFog = oFog.x;
    xOut.oPts = oPts.x;
    xOut.oB0 = saturate(oB0);
    xOut.oB1 = saturate(oB1);
    // TEXCOORDINDEX remapping: NV2A texture units route interpolated texcoords
    // to stages based on D3DTSS_TEXCOORDINDEX. In D3D11 we apply this in the VS.
    // For ShaderProgram mode, c219 = {0,1,2,3} (identity, no-op).
    {
        float4 texcoordSets[4] = { oT0, oT1, oT2, oT3 };
        oT0 = texcoordSets[(int)xboxTexCoordIndex.x];
        oT1 = texcoordSets[(int)xboxTexCoordIndex.y];
        oT2 = texcoordSets[(int)xboxTexCoordIndex.z];
        oT3 = texcoordSets[(int)xboxTexCoordIndex.w];
    }
    // Scale textures
    xOut.oT0 = oT0 / xboxTextureScale[0];
    xOut.oT1 = oT1 / xboxTextureScale[1];
    xOut.oT2 = oT2 / xboxTextureScale[2];
    xOut.oT3 = oT3 / xboxTextureScale[3];
