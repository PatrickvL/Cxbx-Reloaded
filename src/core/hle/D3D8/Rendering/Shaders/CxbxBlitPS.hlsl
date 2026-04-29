// Simple texture sampling pixel shader
Texture2D tex : register(t0);
SamplerState samp : register(s0);
float4 main(float4 pos : SV_Position, float2 uv : TEXCOORD0) : SV_Target {
    return tex.Sample(samp, uv);
}
