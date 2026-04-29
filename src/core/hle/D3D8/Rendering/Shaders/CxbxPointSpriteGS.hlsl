// Point sprite geometry shader — expands point into a screen-aligned quad
struct GS_INPUT {
    float4 oPos : SV_Position;
    float4 oD0  : COLOR0;
    float4 oD1  : COLOR1;
    float  oFog : FOG;
    float  oPts : PSIZE;
    float4 oB0  : TEXCOORD4;
    float4 oB1  : TEXCOORD5;
    float4 oT0  : TEXCOORD0;
    float4 oT1  : TEXCOORD1;
    float4 oT2  : TEXCOORD2;
    float4 oT3  : TEXCOORD3;
};
struct GS_OUTPUT {
    float4 oPos : SV_Position;
    float4 oD0  : COLOR0;
    float4 oD1  : COLOR1;
    float  oFog : FOG;
    float  oPts : PSIZE;
    float4 oB0  : TEXCOORD4;
    float4 oB1  : TEXCOORD5;
    float4 oT0  : TEXCOORD0;
    float4 oT1  : TEXCOORD1;
    float4 oT2  : TEXCOORD2;
    float4 oT3  : TEXCOORD3;
};
cbuffer GSConstants : register(b0) { float4 gsViewportInv; };
[maxvertexcount(4)]
void main(point GS_INPUT input[1], inout TriangleStream<GS_OUTPUT> stream) {
    GS_OUTPUT o;
    o.oD0  = input[0].oD0;
    o.oD1  = input[0].oD1;
    o.oFog = input[0].oFog;
    o.oPts = input[0].oPts;
    o.oB0  = input[0].oB0;
    o.oB1  = input[0].oB1;
    o.oT1  = input[0].oT1;
    o.oT2  = input[0].oT2;
    o.oT3  = input[0].oT3;
    float ptSize = max(input[0].oPts, 1.0);
    float4 pos = input[0].oPos;
    float halfW = ptSize * gsViewportInv.x * pos.w;
    float halfH = ptSize * gsViewportInv.y * pos.w;
    o.oPos = pos + float4(-halfW, -halfH, 0, 0);
    o.oT0 = float4(0, 0, 0, 1);
    stream.Append(o);
    o.oPos = pos + float4(+halfW, -halfH, 0, 0);
    o.oT0 = float4(1, 0, 0, 1);
    stream.Append(o);
    o.oPos = pos + float4(-halfW, +halfH, 0, 0);
    o.oT0 = float4(0, 1, 0, 1);
    stream.Append(o);
    o.oPos = pos + float4(+halfW, +halfH, 0, 0);
    o.oT0 = float4(1, 1, 0, 1);
    stream.Append(o);
    stream.RestartStrip();
}
