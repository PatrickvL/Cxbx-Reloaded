// CxbxVertexInputLoad.hlsli — Load all 16 vertex input registers (v0..v15)
//
// Requires: VS_INPUT xIn in scope, float4 v0..v15 declared.
// Uses SV_VertexID to fetch all attributes from ByteAddressBuffer.

    // Fetch all attributes from ByteAddressBuffer using SV_VertexID
    {
        uint xboxVtxIdx = ResolveVertexIndex(xIn.vertexId);
        float4 vArr[16];
        FetchAllAttributes(xboxVtxIdx, vArr);
        v0=vArr[0]; v1=vArr[1]; v2=vArr[2]; v3=vArr[3];
        v4=vArr[4]; v5=vArr[5]; v6=vArr[6]; v7=vArr[7];
        v8=vArr[8]; v9=vArr[9]; v10=vArr[10]; v11=vArr[11];
        v12=vArr[12]; v13=vArr[13]; v14=vArr[14]; v15=vArr[15];
    }
