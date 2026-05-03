// Texture unswizzle compute shader — Morton (Z-order) decode, typed float4 output.
// Used for formats that support same-format typed UAV but NOT R32/R16/R8_UINT
// cross-family reinterpretation (e.g. B8G8R8A8, B4G4R4A4, B5G6R5, B5G5R5A1, R10G10B10A2).
ByteAddressBuffer g_SrcBuffer : register(t0);
RWTexture2D<float4> g_DstTexture : register(u0);
cbuffer UnswizzleConstants : register(b0) {
    uint maskX; uint maskY; uint texWidth; uint texHeight; uint bpp;
    uint fmtDecode; // 0=BGRA8, 1=B4G4R4A4, 2=B5G6R5, 3=B5G5R5A1, 4=R10G10B10A2
    uint pad0; uint pad1;
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

float4 DecodeBGRA8(uint value) {
    float b = float((value >>  0) & 0xFF) / 255.0;
    float g = float((value >>  8) & 0xFF) / 255.0;
    float r = float((value >> 16) & 0xFF) / 255.0;
    float a = float((value >> 24) & 0xFF) / 255.0;
    return float4(r, g, b, a);
}

float4 DecodeB4G4R4A4(uint value) {
    float b = float((value >>  0) & 0xF) / 15.0;
    float g = float((value >>  4) & 0xF) / 15.0;
    float r = float((value >>  8) & 0xF) / 15.0;
    float a = float((value >> 12) & 0xF) / 15.0;
    return float4(r, g, b, a);
}

float4 DecodeB5G6R5(uint value) {
    float b = float((value >>  0) & 0x1F) / 31.0;
    float g = float((value >>  5) & 0x3F) / 63.0;
    float r = float((value >> 11) & 0x1F) / 31.0;
    return float4(r, g, b, 1.0);
}

float4 DecodeB5G5R5A1(uint value) {
    float b = float((value >>  0) & 0x1F) / 31.0;
    float g = float((value >>  5) & 0x1F) / 31.0;
    float r = float((value >> 10) & 0x1F) / 31.0;
    float a = float((value >> 15) & 0x1);
    return float4(r, g, b, a);
}

float4 DecodeR10G10B10A2(uint value) {
    float r = float((value >>  0) & 0x3FF) / 1023.0;
    float g = float((value >> 10) & 0x3FF) / 1023.0;
    float b = float((value >> 20) & 0x3FF) / 1023.0;
    float a = float((value >> 30) & 0x3) / 3.0;
    return float4(r, g, b, a);
}

[numthreads(8, 8, 1)]
void main(uint3 dtid : SV_DispatchThreadID) {
    uint x = dtid.x; uint y = dtid.y;
    if (x >= texWidth || y >= texHeight) return;
    uint mortonIdx = MortonIndex(x, y);
    uint srcByteOffset = mortonIdx * bpp;
    uint value;
    if (bpp == 4) {
        value = g_SrcBuffer.Load(srcByteOffset);
    } else {
        uint dwordAddr = srcByteOffset & ~3u;
        uint shift = (srcByteOffset & 2u) * 8u;
        value = (g_SrcBuffer.Load(dwordAddr) >> shift) & 0xFFFF;
    }
    float4 color;
    switch (fmtDecode) {
    case 1:  color = DecodeB4G4R4A4(value); break;
    case 2:  color = DecodeB5G6R5(value);   break;
    case 3:  color = DecodeB5G5R5A1(value); break;
    case 4:  color = DecodeR10G10B10A2(value); break;
    default: color = DecodeBGRA8(value);    break;
    }
    g_DstTexture[uint2(x, y)] = color;
}
