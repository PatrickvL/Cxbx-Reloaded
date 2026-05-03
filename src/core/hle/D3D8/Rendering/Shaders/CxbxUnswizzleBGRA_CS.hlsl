// Texture unswizzle compute shader — Morton (Z-order) decode
// Variant for B8G8R8A8_UNORM textures that cannot use R32_UINT UAV reinterpretation.
// Writes float4 to a same-format (B8G8R8A8_UNORM) UAV instead.
ByteAddressBuffer g_SrcBuffer : register(t0);
RWTexture2D<float4> g_DstTexture : register(u0);
cbuffer UnswizzleConstants : register(b0) {
    uint maskX; uint maskY; uint texWidth; uint texHeight; uint bpp;
};
uint MortonIndex(uint x, uint y) {
    uint mx = maskX; uint my = maskY;
    uint result = 0;
    uint xBit = 1, yBit = 1, outBit = 1;
    uint totalMask = mx | my;
    [unroll(20)]
    for (uint i = 0; i < 20; i++) {
        if ((totalMask & outBit) == 0) break;
        if (mx & outBit) { if (x & xBit) result |= outBit; xBit <<= 1; }
        if (my & outBit) { if (y & yBit) result |= outBit; yBit <<= 1; }
        outBit <<= 1;
    }
    return result;
}
[numthreads(8, 8, 1)]
void main(uint3 dtid : SV_DispatchThreadID) {
    uint x = dtid.x; uint y = dtid.y;
    if (x >= texWidth || y >= texHeight) return;
    uint mortonIdx = MortonIndex(x, y);
    uint srcByteOffset = mortonIdx * 4; // Always 4 bpp for BGRA
    uint value = g_SrcBuffer.Load(srcByteOffset);
    // Unpack uint (memory order: B=byte0, G=byte1, R=byte2, A=byte3) to float4.
    // B8G8R8A8_UNORM UAV stores: byte0=.b*255, byte1=.g*255, byte2=.r*255, byte3=.a*255
    float b = float((value >>  0) & 0xFF) / 255.0;
    float g = float((value >>  8) & 0xFF) / 255.0;
    float r = float((value >> 16) & 0xFF) / 255.0;
    float a = float((value >> 24) & 0xFF) / 255.0;
    g_DstTexture[uint2(x, y)] = float4(r, g, b, a);
}
