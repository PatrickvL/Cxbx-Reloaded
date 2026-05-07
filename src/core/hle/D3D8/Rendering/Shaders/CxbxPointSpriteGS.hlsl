// Point sprite geometry shader — expands point into a screen-aligned quad
#include "CxbxVertexShaderCommon.hlsli"

// GS_INPUT is identical to VS_OUTPUT (signature must match for D3D11 VS→GS linking)
#define GS_INPUT VS_OUTPUT
// GS_OUTPUT must also use VS_OUTPUT so the GS→PS signature matches
#define GS_OUTPUT VS_OUTPUT

cbuffer GSConstants : register(b0) { float4 gsViewportInv; };
[maxvertexcount(4)]
void main(point GS_INPUT input[1], inout TriangleStream<GS_OUTPUT> stream) {
    float4 pos = input[0].oPos;
    GS_OUTPUT o;
    o.oD0  = input[0].oD0;
    o.oD1  = input[0].oD1;
    o.oFog = input[0].oFog;
    o.oPts = input[0].oPts;
    o.oB0  = input[0].oB0;
    o.oB1  = input[0].oB1;
    float ptSize = max(input[0].oPts, 1.0);
    float halfW = ptSize * gsViewportInv.x * pos.w;
    float halfH = ptSize * gsViewportInv.y * pos.w;
    // NV2A generates point sprite UVs for ALL enabled texture stages.
    float4 uv;
    uv = float4(0, 0, 0, 1);
    o.oT0 = uv; o.oT1 = uv; o.oT2 = uv; o.oT3 = uv;
    o.oPos = pos + float4(-halfW, -halfH, 0, 0);
    stream.Append(o);
    uv = float4(1, 0, 0, 1);
    o.oT0 = uv; o.oT1 = uv; o.oT2 = uv; o.oT3 = uv;
    o.oPos = pos + float4(+halfW, -halfH, 0, 0);
    stream.Append(o);
    uv = float4(0, 1, 0, 1);
    o.oT0 = uv; o.oT1 = uv; o.oT2 = uv; o.oT3 = uv;
    o.oPos = pos + float4(-halfW, +halfH, 0, 0);
    stream.Append(o);
    uv = float4(1, 1, 0, 1);
    o.oT0 = uv; o.oT1 = uv; o.oT2 = uv; o.oT3 = uv;
    o.oPos = pos + float4(+halfW, +halfH, 0, 0);
    stream.Append(o);
    stream.RestartStrip();
}
