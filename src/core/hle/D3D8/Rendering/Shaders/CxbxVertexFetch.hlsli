// CxbxVertexFetch.hlsli — Vertex fetch from ByteAddressBuffer
//
// The vertex shader receives only SV_VertexID. All vertex attribute
// data is fetched from a ByteAddressBuffer containing the raw Xbox
// vertex stream data. Topology conversion (quad→tri, fan→tri) and
// format decode (including NORMPACKED3/CMP 11.11.10) are performed
// inline.

#ifndef CXBX_VERTEX_FETCH_HLSLI
#define CXBX_VERTEX_FETCH_HLSLI

// ---------------------------------------------------------------
// Vertex fetch SRV — raw Xbox vertex data
// ---------------------------------------------------------------
ByteAddressBuffer g_VtxData : register(t0);

// ---------------------------------------------------------------
// Index buffer SRV (optional — only bound for indexed draws)
// ---------------------------------------------------------------
ByteAddressBuffer g_IdxData : register(t1);

// ---------------------------------------------------------------
// Typed SRV views over the same vertex data buffer (hardware format decode)
// These view the same underlying buffer as t0, but with typed formats.
// Element index = byteOffset / 4 (all elements are 4 bytes).
// ---------------------------------------------------------------
Buffer<float2> g_VtxSNorm16x2 : register(t2);  // R16G16_SNORM: 2 signed normalized shorts
Buffer<float4> g_VtxUNorm8x4  : register(t3);  // R8G8B8A8_UNORM: 4 unsigned normalized bytes

// ---------------------------------------------------------------
// Vertex layout constant buffer (b1)
// ---------------------------------------------------------------
// Header: general draw state
// Per-attribute: up to 16 attribute descriptors
cbuffer CxbxVertexLayoutCB : register(b1)
{
    // Header (8 uint = 32 bytes, 2 × uint4 for alignment)
    uint g_PrimType;        // 0=normal, 1=quad, 2=fan, 3=quadstrip, 4=lineloop
    uint g_IndexedDraw;     // 0=non-indexed, 1=indexed 16-bit, 2=indexed 32-bit
    uint g_IndexOffset;     // Byte offset into g_IdxData for the index data start
    uint g_NumAttribs;      // Number of active vertex attributes (1..16)
    uint g_NumVerts;        // Original Xbox vertex count (needed for lineloop wrap)
    uint g_VertexOffset;    // Added to resolved index (StartVertex or BaseVertexIndex)
    uint g_Pad6;
    uint g_Pad7;

    // Per-attribute descriptors (16 × uint4 = 256 bytes)
    // x = byte offset from start of vertex in the stream
    // y = stride (bytes per vertex for this stream)
    // z = format (CXBX_VTXFMT_* constant)
    // w = base byte offset of the stream within g_VtxData
    uint4 g_Attribs[16];
};

// ---------------------------------------------------------------
// Vertex format constants (matches C++ CXBX_VTXFMT_*)
// ---------------------------------------------------------------
#define CXBX_VTXFMT_FLOAT1       0
#define CXBX_VTXFMT_FLOAT2       1
#define CXBX_VTXFMT_FLOAT3       2
#define CXBX_VTXFMT_FLOAT4       3
#define CXBX_VTXFMT_D3DCOLOR     4
#define CXBX_VTXFMT_SHORT2       5
#define CXBX_VTXFMT_SHORT4       6
#define CXBX_VTXFMT_NORMPACKED3  7  // 11.11.10 signed packed (CMP)
#define CXBX_VTXFMT_SHORT2N      8
#define CXBX_VTXFMT_SHORT4N      9
#define CXBX_VTXFMT_PBYTE4       10
#define CXBX_VTXFMT_FLOAT2H      11 // Xbox "half" with W: x,y,0,w stored as 3 floats
#define CXBX_VTXFMT_NONE         12 // Use default value (sticky register)
#define CXBX_VTXFMT_SHORT1N      13 // 1 signed 16-bit normalized (Xbox: 2 bytes)
#define CXBX_VTXFMT_SHORT3N      14 // 3 signed 16-bit normalized (Xbox: 6 bytes, W=1.0)
#define CXBX_VTXFMT_PBYTE1       15 // 1 unsigned byte normalized (Xbox: 1 byte)
#define CXBX_VTXFMT_PBYTE2       16 // 2 unsigned bytes normalized (Xbox: 2 bytes)
#define CXBX_VTXFMT_PBYTE3       17 // 3 unsigned bytes normalized (Xbox: 3 bytes, A=1.0)
#define CXBX_VTXFMT_SHORT1       18 // 1 signed 16-bit unnormalized (Xbox: 2 bytes)
#define CXBX_VTXFMT_SHORT3       19 // 3 signed 16-bit unnormalized (Xbox: 6 bytes, W=1)

// Primitive type constants for topology conversion
#define CXBX_PRIM_NORMAL    0
#define CXBX_PRIM_QUAD      1
#define CXBX_PRIM_FAN       2
#define CXBX_PRIM_QUADSTRIP 3
#define CXBX_PRIM_LINELOOP  4

// ---------------------------------------------------------------
// Read helpers
// ---------------------------------------------------------------

// Read an unaligned uint32 from the vertex data buffer
uint ReadU32(uint byteOff)
{
    uint a = byteOff & ~3u;
    uint s = (byteOff & 3u) * 8u;
    if (s == 0u) return g_VtxData.Load(a);
    return (g_VtxData.Load(a) >> s) | (g_VtxData.Load(a + 4u) << (32u - s));
}

// Read an unaligned uint16 from the vertex data buffer
uint ReadU16(uint byteOff)
{
    uint a = byteOff & ~3u;
    uint s = (byteOff & 3u) * 8u;
    uint d = g_VtxData.Load(a);
    if (s <= 16u) return (d >> s) & 0xFFFFu;
    return ((d >> s) | (g_VtxData.Load(a + 4u) << (32u - s))) & 0xFFFFu;
}

// ---------------------------------------------------------------
// Topology conversion: compute Xbox vertex index from SV_VertexID
// ---------------------------------------------------------------

// Quad list: each quad [0,1,2,3] becomes two triangles [0,1,2] and [0,2,3]
uint QuadVertexIndex(uint vertId)
{
    uint quad  = vertId / 6u;
    uint local = vertId % 6u;
    // LUT: 0,1,2, 0,2,3
    static const uint lut[6] = { 0u, 1u, 2u, 0u, 2u, 3u };
    return quad * 4u + lut[local];
}

// Triangle fan: fan apex is vertex 0, each tri uses (0, i+1, i+2)
uint FanVertexIndex(uint vertId)
{
    uint tri   = vertId / 3u;
    uint local = vertId % 3u;
    // local 0 → apex (0), local 1 → tri+1, local 2 → tri+2
    return local == 0u ? 0u : (tri + local);
}

// Quad strip: each quad uses vertices [2i, 2i+1, 2i+2, 2i+3]
// Triangulated as [2i, 2i+1, 2i+2] and [2i+2, 2i+1, 2i+3]
uint QuadStripVertexIndex(uint vertId)
{
    uint quad  = vertId / 6u;
    uint local = vertId % 6u;
    uint base  = quad * 2u;
    // LUT: +0,+1,+2, +2,+1,+3
    static const uint lut[6] = { 0u, 1u, 2u, 2u, 1u, 3u };
    return base + lut[local];
}

// Line loop: N vertices → N line segments rendered as LINELIST (2N host vertices)
// Segment i uses vertices [i, (i+1) % N]
uint LineLoopVertexIndex(uint vertId, uint numVerts)
{
    uint seg   = vertId / 2u;
    uint local = vertId % 2u;
    // local 0 → seg, local 1 → (seg+1) % numVerts
    return local == 0u ? seg : ((seg + 1u) % numVerts);
}

// Resolve the Xbox vertex index from the host SV_VertexID, applying
// topology conversion and optional index buffer indirection.
uint ResolveVertexIndex(uint hostVertId)
{
    // Step 1: topology remapping
    uint logicalIdx;
    if (g_PrimType == CXBX_PRIM_QUAD) {
        logicalIdx = QuadVertexIndex(hostVertId);
    } else if (g_PrimType == CXBX_PRIM_FAN) {
        logicalIdx = FanVertexIndex(hostVertId);
    } else if (g_PrimType == CXBX_PRIM_QUADSTRIP) {
        logicalIdx = QuadStripVertexIndex(hostVertId);
    } else if (g_PrimType == CXBX_PRIM_LINELOOP) {
        logicalIdx = LineLoopVertexIndex(hostVertId, g_NumVerts);
    } else {
        logicalIdx = hostVertId;
    }

    // Step 2: index buffer indirection (if indexed draw)
    if (g_IndexedDraw == 1u) {
        // 16-bit indices
        uint byteAddr = g_IndexOffset + logicalIdx * 2u;
        uint aligned  = byteAddr & ~3u;
        uint shift    = (byteAddr & 3u) * 8u;
        uint word     = g_IdxData.Load(aligned);
        return ((word >> shift) & 0xFFFFu) + g_VertexOffset;
    } else if (g_IndexedDraw == 2u) {
        // 32-bit indices
        return g_IdxData.Load(g_IndexOffset + logicalIdx * 4u) + g_VertexOffset;
    }

    return logicalIdx + g_VertexOffset;
}

// ---------------------------------------------------------------
// Format decode helpers
// ---------------------------------------------------------------

// Sign-extend a value from 'bits' width to 32-bit int
int SignExtend(uint val, uint bits)
{
    uint signBit = 1u << (bits - 1u);
    return (int)((val ^ signBit) - signBit);
}

// NORMPACKED3 (CMP): 11.11.10 signed packed normal
// Bits [10:0]  = X (11-bit signed)
// Bits [21:11] = Y (11-bit signed)
// Bits [31:22] = Z (10-bit signed)
float3 DecodeNormPacked3(uint raw)
{
    int x = (int)(raw << 21u) >> 21;  // sign-extend 11 bits
    int y = (int)(raw << 10u) >> 21;  // sign-extend 11 bits
    int z = (int)(raw) >> 22;         // sign-extend 10 bits (arithmetic shift)
    return float3(
        (x >= 0) ? ((float)x / 1023.0f) : ((float)x / 1024.0f),
        (y >= 0) ? ((float)y / 1023.0f) : ((float)y / 1024.0f),
        (z >= 0) ? ((float)z / 511.0f)  : ((float)z / 512.0f)
    );
}

// D3DCOLOR: Xbox stores as BGRA (B in low byte), we need RGBA
float4 DecodeD3DColor(uint raw)
{
    float b = (float)((raw      ) & 0xFFu) / 255.0f;
    float g = (float)((raw >> 8u) & 0xFFu) / 255.0f;
    float r = (float)((raw >>16u) & 0xFFu) / 255.0f;
    float a = (float)((raw >>24u) & 0xFFu) / 255.0f;
    return float4(r, g, b, a);
}

// SHORT2N: 2 signed 16-bit normalized values
float4 DecodeShort2N(uint byteOff)
{
    uint lo = ReadU16(byteOff);
    uint hi = ReadU16(byteOff + 2u);
    float x = (float)SignExtend(lo, 16u) / 32767.0f;
    float y = (float)SignExtend(hi, 16u) / 32767.0f;
    return float4(x, y, 0.0f, 1.0f);
}

// SHORT4N: 4 signed 16-bit normalized values
float4 DecodeShort4N(uint byteOff)
{
    float x = (float)SignExtend(ReadU16(byteOff),      16u) / 32767.0f;
    float y = (float)SignExtend(ReadU16(byteOff + 2u), 16u) / 32767.0f;
    float z = (float)SignExtend(ReadU16(byteOff + 4u), 16u) / 32767.0f;
    float w = (float)SignExtend(ReadU16(byteOff + 6u), 16u) / 32767.0f;
    return float4(x, y, z, w);
}

// SHORT2: 2 signed 16-bit unnormalized values (as float)
float4 DecodeShort2(uint byteOff)
{
    float x = (float)SignExtend(ReadU16(byteOff),      16u);
    float y = (float)SignExtend(ReadU16(byteOff + 2u), 16u);
    return float4(x, y, 0.0f, 1.0f);
}

// SHORT4: 4 signed 16-bit unnormalized values (as float)
float4 DecodeShort4(uint byteOff)
{
    float x = (float)SignExtend(ReadU16(byteOff),      16u);
    float y = (float)SignExtend(ReadU16(byteOff + 2u), 16u);
    float z = (float)SignExtend(ReadU16(byteOff + 4u), 16u);
    float w = (float)SignExtend(ReadU16(byteOff + 6u), 16u);
    return float4(x, y, z, w);
}

// PBYTE4: 4 unsigned bytes normalized to [0..1]
float4 DecodePByte4(uint raw)
{
    float x = (float)((raw      ) & 0xFFu) / 255.0f;
    float y = (float)((raw >> 8u) & 0xFFu) / 255.0f;
    float z = (float)((raw >>16u) & 0xFFu) / 255.0f;
    float w = (float)((raw >>24u) & 0xFFu) / 255.0f;
    return float4(x, y, z, w);
}

// FLOAT2H: Xbox stores {x, y, w} as 3 consecutive floats → output {x, y, 0, w}
float4 DecodeFloat2H(uint byteOff)
{
    float x = asfloat(ReadU32(byteOff));
    float y = asfloat(ReadU32(byteOff + 4u));
    float w = asfloat(ReadU32(byteOff + 8u));
    return float4(x, y, 0.0f, w);
}

// ---------------------------------------------------------------
// Main attribute fetch — 2-level dispatch for fast O3 compilation
// Level 1: [branch] selects category (float / typed SRV / raw ALU)
// Level 2: small switch within each category (max 6 cases)
// This prevents 16×20 combinatorial explosion while keeping [unroll].
// ---------------------------------------------------------------
float4 FetchAttribute(uint xboxVtxIdx, uint4 attribDesc, float4 defaultVal)
{
    uint elemOffset = attribDesc.x;  // offset within the vertex
    uint stride     = attribDesc.y;  // bytes per vertex in this stream
    uint fmt        = attribDesc.z;  // CXBX_VTXFMT_* format
    uint streamBase = attribDesc.w;  // base byte offset of stream in g_VtxData

    if (fmt == CXBX_VTXFMT_NONE)
        return defaultVal;

    uint byteOff = streamBase + xboxVtxIdx * stride + elemOffset;

    // Category A: float formats (0-3) — just asfloat loads
    [branch] if (fmt <= CXBX_VTXFMT_FLOAT4) {
        switch (fmt) {
        case CXBX_VTXFMT_FLOAT1:
            return float4(asfloat(ReadU32(byteOff)), 0.0f, 0.0f, 1.0f);
        case CXBX_VTXFMT_FLOAT2:
            return float4(asfloat(ReadU32(byteOff)), asfloat(ReadU32(byteOff + 4u)), 0.0f, 1.0f);
        case CXBX_VTXFMT_FLOAT3:
            return float4(asfloat(ReadU32(byteOff)), asfloat(ReadU32(byteOff + 4u)), asfloat(ReadU32(byteOff + 8u)), 1.0f);
        default: // FLOAT4
            return float4(asfloat(ReadU32(byteOff)), asfloat(ReadU32(byteOff + 4u)), asfloat(ReadU32(byteOff + 8u)), asfloat(ReadU32(byteOff + 12u)));
        }
    }

    // Category B: typed SRV formats — hardware decode via Buffer<T>.Load
    // D3DCOLOR(4), SHORT2N(8), SHORT4N(9), PBYTE4(10), SHORT1N(13), SHORT3N(14), PBYTE1(15), PBYTE2(16), PBYTE3(17)
    [branch] if (fmt == CXBX_VTXFMT_D3DCOLOR || (fmt >= CXBX_VTXFMT_SHORT2N && fmt <= CXBX_VTXFMT_PBYTE4) ||
                 (fmt >= CXBX_VTXFMT_SHORT1N && fmt <= CXBX_VTXFMT_PBYTE3)) {
        uint elem = byteOff / 4u;
        switch (fmt) {
        case CXBX_VTXFMT_D3DCOLOR:
            return g_VtxUNorm8x4.Load(elem).zyxw;
        case CXBX_VTXFMT_SHORT2N:
            return float4(g_VtxSNorm16x2.Load(elem), 0.0f, 1.0f);
        case CXBX_VTXFMT_SHORT4N:
            return float4(g_VtxSNorm16x2.Load(elem), g_VtxSNorm16x2.Load(elem + 1u));
        case CXBX_VTXFMT_PBYTE4:
            return g_VtxUNorm8x4.Load(elem);
        case CXBX_VTXFMT_SHORT1N:
            return float4((float)SignExtend(ReadU16(byteOff), 16u) / 32767.0f, 0.0f, 0.0f, 1.0f);
        case CXBX_VTXFMT_SHORT3N:
            return float4(g_VtxSNorm16x2.Load(elem), g_VtxSNorm16x2.Load(elem + 1u).x, 1.0f);
        case CXBX_VTXFMT_PBYTE1:
            return float4(g_VtxUNorm8x4.Load(elem).x, 0.0f, 0.0f, 1.0f);
        case CXBX_VTXFMT_PBYTE2:
            return float4(g_VtxUNorm8x4.Load(elem).xy, 0.0f, 1.0f);
        default: // PBYTE3
            return float4(g_VtxUNorm8x4.Load(elem).xyz, 1.0f);
        }
    }

    // Category C: raw ALU decode (SHORT2, SHORT4, NORMPACKED3, FLOAT2H, SHORT1, SHORT3)
    switch (fmt) {
    case CXBX_VTXFMT_SHORT2:
        return DecodeShort2(byteOff);
    case CXBX_VTXFMT_SHORT4:
        return DecodeShort4(byteOff);
    case CXBX_VTXFMT_NORMPACKED3:
        return float4(DecodeNormPacked3(ReadU32(byteOff)), 1.0f);
    case CXBX_VTXFMT_FLOAT2H:
        return DecodeFloat2H(byteOff);
    case CXBX_VTXFMT_SHORT1:
        return float4((float)SignExtend(ReadU16(byteOff), 16u), 0.0f, 0.0f, 1.0f);
    case CXBX_VTXFMT_SHORT3:
    {
        float sx = (float)SignExtend(ReadU16(byteOff),      16u);
        float sy = (float)SignExtend(ReadU16(byteOff + 2u), 16u);
        float sz = (float)SignExtend(ReadU16(byteOff + 4u), 16u);
        return float4(sx, sy, sz, 1.0f);
    }
    default:
        return defaultVal;
    }
}

// ---------------------------------------------------------------
// Vertex defaults (NV2A sticky attribute values) — stored in b1
// as part of the layout CB following the attrib descriptors.
// For simplicity, these are uploaded as 16 × float4 at the end
// of the constant buffer.
// ---------------------------------------------------------------
cbuffer CxbxVertexDefaultsCB : register(b2)
{
    float4 g_VtxDefaults[16];
};

// ---------------------------------------------------------------
// Fetch all 16 vertex attributes for a given Xbox vertex index
// ---------------------------------------------------------------
void FetchAllAttributes(uint xboxVtxIdx, out float4 v[16])
{
    [unroll]
    for (uint i = 0u; i < 16u; i++) {
        v[i] = FetchAttribute(xboxVtxIdx, g_Attribs[i], g_VtxDefaults[i]);
    }
}

#endif // CXBX_VERTEX_FETCH_HLSLI
