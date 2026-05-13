// Texture unswizzle compute shader — Morton (Z-order) decode
ByteAddressBuffer g_SrcBuffer : register(t0);
RWTexture2D<uint> g_DstTexture : register(u0);
cbuffer UnswizzleConstants : register(b0) {
    uint maskX; uint maskY; uint texWidth; uint texHeight; uint bpp;
    uint pad0; uint pad1;
    uint srcOffset; // byte offset into g_SrcBuffer (0 = staging, nonzero = mirror)
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
    uint srcByteOffset = srcOffset + mortonIdx * bpp;
    uint value;
    if (bpp == 4) {
        value = g_SrcBuffer.Load(srcByteOffset);
    } else if (bpp == 2) {
        uint dwordAddr = srcByteOffset & ~3u;
        uint shift = (srcByteOffset & 2u) * 8u;
        value = (g_SrcBuffer.Load(dwordAddr) >> shift) & 0xFFFF;
    } else {
        uint dwordAddr = srcByteOffset & ~3u;
        uint shift = (srcByteOffset & 3u) * 8u;
        value = (g_SrcBuffer.Load(dwordAddr) >> shift) & 0xFF;
    }
    g_DstTexture[uint2(x, y)] = value;
}
