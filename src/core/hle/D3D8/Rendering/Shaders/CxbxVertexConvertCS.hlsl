// Vertex format conversion compute shader
// Converts Xbox vertex formats to D3D11-compatible layouts on GPU.
ByteAddressBuffer g_SrcBuffer : register(t0);
RWBuffer<uint> g_DstBuffer : register(u0);
cbuffer VertexConvertCB : register(b0) {
    uint g_VertexCount; uint g_SrcStride; uint g_DstStride; uint g_NumElements;
    uint4 g_Elements[16]; // x=srcOffset, y=dstOffset, z=convType, w=copyDwords
};
#define CONV_COPY        0
#define CONV_NORMSHORT3  1
#define CONV_NORMPACKED3 2
#define CONV_SHORT3      3
#define CONV_PBYTE3      4
#define CONV_FLOAT2H     5
#define CONV_D3DCOLOR    6
#define CONV_NONE        7
uint ReadU32(uint byteOff) {
    uint a = byteOff & ~3u; uint s = (byteOff & 3u) * 8u;
    if (s == 0) return g_SrcBuffer.Load(a);
    return (g_SrcBuffer.Load(a) >> s) | (g_SrcBuffer.Load(a + 4) << (32 - s));
}
uint ReadU16(uint byteOff) {
    uint a = byteOff & ~3u; uint s = (byteOff & 3u) * 8u;
    uint d = g_SrcBuffer.Load(a);
    if (s <= 16) return (d >> s) & 0xFFFF;
    return ((d >> s) | (g_SrcBuffer.Load(a + 4) << (32 - s))) & 0xFFFF;
}
void WriteU32(uint dstByteOff, uint value) { g_DstBuffer[dstByteOff >> 2] = value; }
float PackedIntToFloat(int v, float pm, float nm) {
    return (v >= 0) ? (float(v) / pm) : (float(v) / nm);
}
[numthreads(64, 1, 1)]
void main(uint3 dtid : SV_DispatchThreadID) {
    uint vi = dtid.x;
    if (vi >= g_VertexCount) return;
    uint sb = vi * g_SrcStride;
    uint db = vi * g_DstStride;
    for (uint e = 0; e < g_NumElements; e++) {
        uint so = sb + g_Elements[e].x;
        uint do_ = db + g_Elements[e].y;
        uint ct = g_Elements[e].z;
        uint cw = g_Elements[e].w;
        if (ct == CONV_COPY) {
            for (uint d = 0; d < cw; d++) WriteU32(do_ + d * 4, ReadU32(so + d * 4));
        } else if (ct == CONV_NORMSHORT3) {
            WriteU32(do_, ReadU32(so));
            WriteU32(do_ + 4, ReadU16(so + 4) | 0x7FFF0000u);
        } else if (ct == CONV_NORMPACKED3) {
            int p = asint(ReadU32(so));
            int x = (p << 21) >> 21; int y = (p << 10) >> 21; int z = p >> 22;
            WriteU32(do_, asuint(PackedIntToFloat(x, 1023.0, 1024.0)));
            WriteU32(do_ + 4, asuint(PackedIntToFloat(y, 1023.0, 1024.0)));
            WriteU32(do_ + 8, asuint(PackedIntToFloat(z, 511.0, 512.0)));
        } else if (ct == CONV_SHORT3) {
            WriteU32(do_, ReadU32(so));
            WriteU32(do_ + 4, ReadU16(so + 4) | (1u << 16));
        } else if (ct == CONV_PBYTE3) {
            WriteU32(do_, (ReadU32(so) & 0x00FFFFFFu) | 0xFF000000u);
        } else if (ct == CONV_FLOAT2H) {
            WriteU32(do_, ReadU32(so));
            WriteU32(do_ + 4, ReadU32(so + 4));
            WriteU32(do_ + 8, 0);
            WriteU32(do_ + 12, ReadU32(so + 8));
        } else if (ct == CONV_D3DCOLOR) {
            uint c = ReadU32(so);
            WriteU32(do_, ((c >> 16) & 0xFFu) | (c & 0xFF00FF00u) | ((c & 0xFFu) << 16));
        } // CONV_NONE: skip
    }
}
