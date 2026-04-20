// CxbxVertexInputLoad.hlsli — Load all 16 vertex input registers (v0..v15)
//
// Requires: VS_INPUT xIn in scope, float4 v0..v15 declared.
// Under IA bypass: uses SV_VertexID to fetch from ByteAddressBuffer.
// Under D3D11 normal: reads from the 16-element TEXCOORD array.
// Under D3D9: lerps between vertex data and constant buffer defaults.

#ifdef CXBX_IA_BYPASS
    // IA bypass: fetch all attributes from ByteAddressBuffer using SV_VertexID
    {
        uint xboxVtxIdx = ResolveVertexIndex(xIn.vertexId);
        float4 vArr[16];
        FetchAllAttributes(xboxVtxIdx, vArr);
        v0=vArr[0]; v1=vArr[1]; v2=vArr[2]; v3=vArr[3];
        v4=vArr[4]; v5=vArr[5]; v6=vArr[6]; v7=vArr[7];
        v8=vArr[8]; v9=vArr[9]; v10=vArr[10]; v11=vArr[11];
        v12=vArr[12]; v13=vArr[13]; v14=vArr[14]; v15=vArr[15];
    }
#elif defined(CXBX_USE_D3D11)
    // D3D11: The input assembler delivers correct values for all 16 attributes.
    #define init_v(i) v##i = xIn.v[i];
    init_v( 0); init_v( 1); init_v( 2); init_v( 3);
    init_v( 4); init_v( 5); init_v( 6); init_v( 7);
    init_v( 8); init_v( 9); init_v(10); init_v(11);
    init_v(12); init_v(13); init_v(14); init_v(15);
    #undef init_v
#else
    // D3D9: Lerp between vertex data and constant buffer defaults
    float vRegisterDefaultFlags[16] = (float[16])vRegisterDefaultFlagsPacked;
    #define init_v(i) v##i = lerp(xIn.v[i], vRegisterDefaultValues[i], vRegisterDefaultFlags[i]);
    init_v( 0); init_v( 1); init_v( 2); init_v( 3);
    init_v( 4); init_v( 5); init_v( 6); init_v( 7);
    init_v( 8); init_v( 9); init_v(10); init_v(11);
    init_v(12); init_v(13); init_v(14); init_v(15);
    #undef init_v
#endif
