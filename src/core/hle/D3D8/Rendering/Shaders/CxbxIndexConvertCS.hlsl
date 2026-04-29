// Index buffer conversion compute shader
// Converts quad-list and triangle-fan indices to triangle-list on GPU.
cbuffer IndexConvertCB : register(b0) {
    uint g_VertexCount;
    uint g_ConversionMode;
    uint g_IsIndexed;
    uint g_Pad;
};
ByteAddressBuffer InputIndices : register(t0);
RWBuffer<uint> OutputIndices : register(u0);
uint ReadIndex16(uint idx) {
    if (g_IsIndexed == 0) return idx;
    uint byteOff = idx * 2;
    uint dwordOff = byteOff & ~3u;
    uint raw = InputIndices.Load(dwordOff);
    return (byteOff & 2u) ? (raw >> 16) : (raw & 0xFFFF);
}
uint Pack16(uint lo, uint hi) { return (lo & 0xFFFF) | (hi << 16); }
[numthreads(64, 1, 1)]
void main(uint3 DTid : SV_DispatchThreadID) {
    uint tid = DTid.x;
    if (g_ConversionMode <= 1) {
        uint numQuads = g_VertexCount / 4;
        if (tid >= numQuads) return;
        uint s = tid * 4;
        uint A = ReadIndex16(s), B = ReadIndex16(s+1), C = ReadIndex16(s+2), D = ReadIndex16(s+3);
        uint d = tid * 3;
        if (g_ConversionMode == 0) {
            OutputIndices[d+0] = Pack16(A, B);
            OutputIndices[d+1] = Pack16(D, B);
            OutputIndices[d+2] = Pack16(C, D);
        } else {
            OutputIndices[d+0] = Pack16(A, D);
            OutputIndices[d+1] = Pack16(B, B);
            OutputIndices[d+2] = Pack16(D, C);
        }
    } else {
        uint numTris = (g_VertexCount >= 3) ? (g_VertexCount - 2) : 0;
        uint triIdx = tid * 2;
        if (triIdx >= numTris) return;
        uint hub = ReadIndex16(0);
        uint i1 = ReadIndex16(triIdx + 1), i2 = ReadIndex16(triIdx + 2);
        uint d = tid * 3;
        if (triIdx + 1 < numTris) {
            uint i3 = ReadIndex16(triIdx + 3);
            OutputIndices[d+0] = Pack16(hub, i1);
            OutputIndices[d+1] = Pack16(i2, hub);
            OutputIndices[d+2] = Pack16(i2, i3);
        } else {
            OutputIndices[d+0] = Pack16(hub, i1);
            OutputIndices[d+1] = Pack16(i2, 0);
        }
    }
}
