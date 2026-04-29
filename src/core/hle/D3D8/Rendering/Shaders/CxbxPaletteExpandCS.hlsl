// Palette texture expansion compute shader
// Combines unswizzle + P8 palette lookup in a single dispatch.
ByteAddressBuffer g_SrcBuffer : register(t0);
ByteAddressBuffer g_Palette : register(t1);
RWTexture2D<uint> g_DstTexture : register(u0);
cbuffer PaletteExpandCB : register(b0) {
    uint maskX; uint maskY; uint texWidth; uint texHeight;
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
    // Read P8 index byte from swizzled source
    uint dwordAddr = mortonIdx & ~3u;
    uint shift = (mortonIdx & 3u) * 8u;
    uint palIdx = (g_SrcBuffer.Load(dwordAddr) >> shift) & 0xFF;
    // Look up ARGB palette entry (stored as D3DCOLOR = 0xAARRGGBB little-endian = [B,G,R,A] bytes)
    uint argb = g_Palette.Load(palIdx * 4);
    // Swap R and B for R8G8B8A8_UNORM output: BGRA bytes -> RGBA bytes
    uint rgba = ((argb & 0x00FF0000u) >> 16) | (argb & 0x0000FF00u) | ((argb & 0x000000FFu) << 16) | (argb & 0xFF000000u);
    g_DstTexture[uint2(x, y)] = rgba;
}
