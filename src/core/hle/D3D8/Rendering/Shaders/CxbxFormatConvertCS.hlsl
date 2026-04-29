// Texture format conversion compute shader
// Combines unswizzle (optional) + format decode -> RGBA output.
ByteAddressBuffer g_SrcBuffer : register(t0);
RWTexture2D<uint> g_DstTexture : register(u0);
cbuffer FormatConvertCB : register(b0) {
    uint maskX; uint maskY; uint texWidth; uint texHeight;
    uint bpp; uint fmtType; uint swizzled; uint srcRowPitch;
};

uint MortonIndex(uint x, uint y) {
    uint mx = maskX; uint my = maskY;
    uint result = 0;
    uint xBit = 1, yBit = 1, outBit = 1;
    uint totalMask = mx | my;
    [loop]
    for (uint i = 0; i < 20; i++) {
        if ((totalMask & outBit) == 0) break;
        if (mx & outBit) { if (x & xBit) result |= outBit; xBit <<= 1; }
        if (my & outBit) { if (y & yBit) result |= outBit; yBit <<= 1; }
        outBit <<= 1;
    }
    return result;
}

uint LoadSrcByte(uint byteOffset) {
    uint dw = g_SrcBuffer.Load(byteOffset & ~3u);
    return (dw >> ((byteOffset & 3u) * 8u)) & 0xFF;
}
// Assumes byteOffset is 2-byte aligned (true for all 16bpp formats with Morton or pitch addressing)
uint LoadSrc16(uint byteOffset) {
    uint dw = g_SrcBuffer.Load(byteOffset & ~3u);
    return (dw >> ((byteOffset & 2u) * 8u)) & 0xFFFF;
}
uint LoadSrc32(uint byteOffset) {
    return g_SrcBuffer.Load(byteOffset);
}

// Unsigned 5-bit expand [0..31] -> 8-bit [0..255]
uint Expand5to8(uint v) { return (v << 3) | (v >> 2); }
// Unsigned 6-bit expand [0..63] -> 8-bit [0..255]
uint Expand6to8(uint v) { return (v << 2) | (v >> 4); }
// Unsigned 4-bit expand [0..15] -> 8-bit [0..255]
uint Expand4to8(uint v) { return (v << 4) | v; }
// Signed 5-bit bias to unsigned 8-bit (signed -> unsigned range)
uint SignedExpand5to8(uint v) {
    int s = int(v); if (s >= 16) s -= 32;
    // Map [-16..+15] to [0..255]: add 16 then scale
    return uint(clamp((s + 16) * 255 / 31, 0, 255));
}

// Pack R,G,B,A (0-255 each) into a uint for R8G8B8A8_UNORM UAV
uint PackRGBA(uint r, uint g, uint b, uint a) {
    return (r & 0xFF) | ((g & 0xFF) << 8) | ((b & 0xFF) << 16) | ((a & 0xFF) << 24);
}

#define CXBX_FMTCONV_L6V5U5   1
#define CXBX_FMTCONV_R6G5B5   2
#define CXBX_FMTCONV_V8U8     3
#define CXBX_FMTCONV_R5G5B5A1 4
#define CXBX_FMTCONV_R4G4B4A4 5
#define CXBX_FMTCONV_A8       6

uint ConvertTexel16(uint raw, uint fmt) {
    if (fmt == CXBX_FMTCONV_L6V5U5) {
        // L6V5U5: bits [4:0]=U(signed), [9:5]=V(signed), [15:10]=L(unsigned)
        uint u8 = SignedExpand5to8(raw & 0x1Fu);
        uint v8 = SignedExpand5to8((raw >> 5u) & 0x1Fu);
        uint l8 = Expand6to8((raw >> 10u) & 0x3Fu);
        return PackRGBA(u8, v8, l8, 255);
    }
    if (fmt == CXBX_FMTCONV_R6G5B5) {
        // R6G5B5: bits [4:0]=B(5), [9:5]=G(5), [15:10]=R(6)
        uint b8 = Expand5to8(raw & 0x1Fu);
        uint g8 = Expand5to8((raw >> 5u) & 0x1Fu);
        uint r8 = Expand6to8((raw >> 10u) & 0x3Fu);
        return PackRGBA(r8, g8, b8, 255);
    }
    if (fmt == CXBX_FMTCONV_V8U8) {
        // V8U8 (G8B8 alias): byte[0]=U(signed), byte[1]=V(signed)
        uint u8 = ((raw & 0xFFu) + 128u) & 0xFFu;
        uint v8 = (((raw >> 8u) & 0xFFu) + 128u) & 0xFFu;
        return PackRGBA(u8, v8, 0, 255);
    }
    if (fmt == CXBX_FMTCONV_R5G5B5A1) {
        // R5G5B5A1: bit[0]=A(1), [5:1]=B(5), [10:6]=G(5), [15:11]=R(5)
        uint a8 = (raw & 1u) ? 255 : 0;
        uint b8 = Expand5to8((raw >> 1u) & 0x1Fu);
        uint g8 = Expand5to8((raw >> 6u) & 0x1Fu);
        uint r8 = Expand5to8((raw >> 11u) & 0x1Fu);
        return PackRGBA(r8, g8, b8, a8);
    }
    if (fmt == CXBX_FMTCONV_R4G4B4A4) {
        // R4G4B4A4: [3:0]=A(4), [7:4]=B(4), [11:8]=G(4), [15:12]=R(4)
        uint a8 = Expand4to8(raw & 0xFu);
        uint b8 = Expand4to8((raw >> 4u) & 0xFu);
        uint g8 = Expand4to8((raw >> 8u) & 0xFu);
        uint r8 = Expand4to8((raw >> 12u) & 0xFu);
        return PackRGBA(r8, g8, b8, a8);
    }
    return PackRGBA(255, 0, 255, 255); // Magenta = unhandled
}

uint ConvertTexel8(uint raw, uint fmt) {
    if (fmt == CXBX_FMTCONV_A8) {
        // A8: Xbox RGB=1, alpha=raw
        return PackRGBA(255, 255, 255, raw & 0xFFu);
    }
    return PackRGBA(255, 0, 255, 255); // Magenta
}

[numthreads(8, 8, 1)]
void main(uint3 dtid : SV_DispatchThreadID) {
    uint x = dtid.x; uint y = dtid.y;
    if (x >= texWidth || y >= texHeight) return;
    uint srcByteOffset;
    if (swizzled != 0) {
        srcByteOffset = MortonIndex(x, y) * bpp;
    } else { // linear: srcRowPitch is the byte stride between rows
        srcByteOffset = y * srcRowPitch + x * bpp;
    }
    uint rgba;
    if (bpp == 2) {
        uint raw16 = LoadSrc16(srcByteOffset);
        rgba = ConvertTexel16(raw16, fmtType);
    } else if (bpp == 1) {
        uint raw8 = LoadSrcByte(srcByteOffset);
        rgba = ConvertTexel8(raw8, fmtType);
    } else if (bpp == 4) {
        // 32-bit formats: raw passthrough; channel swizzle handled by PS TEXFMTFIXUP
        rgba = LoadSrc32(srcByteOffset);
    } else {
        rgba = PackRGBA(255, 0, 255, 255);
    }
    g_DstTexture[uint2(x, y)] = rgba;
}
