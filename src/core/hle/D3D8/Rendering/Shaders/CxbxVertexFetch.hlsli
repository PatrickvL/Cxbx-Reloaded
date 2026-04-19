// CxbxVertexFetch.hlsli — IA-bypass vertex fetch from ByteAddressBuffer
//
// When CXBX_IA_BYPASS is defined, the vertex shader receives only
// SV_VertexID. All vertex attribute data is fetched manually from
// a ByteAddressBuffer containing the raw Xbox vertex stream data.
// Topology conversion (quad→tri, fan→tri) and format decode
// (including NORMPACKED3/CMP 11.11.10) are performed inline.

#ifndef CXBX_VERTEX_FETCH_HLSLI
#define CXBX_VERTEX_FETCH_HLSLI

#ifdef CXBX_IA_BYPASS

// ---------------------------------------------------------------
// Vertex fetch SRV — raw Xbox vertex data
// ---------------------------------------------------------------
ByteAddressBuffer g_VtxData : register(t0);

// ---------------------------------------------------------------
// Index buffer SRV (optional — only bound for indexed draws)
// ---------------------------------------------------------------
ByteAddressBuffer g_IdxData : register(t1);

// ---------------------------------------------------------------
// Vertex layout constant buffer (b1)
// ---------------------------------------------------------------
// Header: general draw state
// Per-attribute: up to 16 attribute descriptors
cbuffer CxbxVertexLayoutCB : register(b1)
{
    // Header (4 uint = 16 bytes)
    uint g_PrimType;        // 0=normal, 1=quad, 2=fan
    uint g_IndexedDraw;     // 0=non-indexed, 1=indexed 16-bit, 2=indexed 32-bit
    uint g_IndexOffset;     // Byte offset into g_IdxData for the index data start
    uint g_NumAttribs;      // Number of active vertex attributes (1..16)

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

// Primitive type constants for topology conversion
#define CXBX_PRIM_NORMAL  0
#define CXBX_PRIM_QUAD    1
#define CXBX_PRIM_FAN     2

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
        return (word >> shift) & 0xFFFFu;
    } else if (g_IndexedDraw == 2u) {
        // 32-bit indices
        return g_IdxData.Load(g_IndexOffset + logicalIdx * 4u);
    }

    return logicalIdx;
}

// ---------------------------------------------------------------
// Format decode helpers
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
// Main attribute fetch — read one attribute from the vertex buffer
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

    switch (fmt) {
    case CXBX_VTXFMT_FLOAT1:
        return float4(asfloat(ReadU32(byteOff)), 0.0f, 0.0f, 1.0f);
    case CXBX_VTXFMT_FLOAT2:
        return float4(asfloat(ReadU32(byteOff)), asfloat(ReadU32(byteOff + 4u)), 0.0f, 1.0f);
    case CXBX_VTXFMT_FLOAT3:
        return float4(asfloat(ReadU32(byteOff)), asfloat(ReadU32(byteOff + 4u)), asfloat(ReadU32(byteOff + 8u)), 1.0f);
    case CXBX_VTXFMT_FLOAT4:
        return float4(asfloat(ReadU32(byteOff)), asfloat(ReadU32(byteOff + 4u)), asfloat(ReadU32(byteOff + 8u)), asfloat(ReadU32(byteOff + 12u)));
    case CXBX_VTXFMT_D3DCOLOR:
        return DecodeD3DColor(ReadU32(byteOff));
    case CXBX_VTXFMT_SHORT2:
        return DecodeShort2(byteOff);
    case CXBX_VTXFMT_SHORT4:
        return DecodeShort4(byteOff);
    case CXBX_VTXFMT_NORMPACKED3:
        return float4(DecodeNormPacked3(ReadU32(byteOff)), 1.0f);
    case CXBX_VTXFMT_SHORT2N:
        return DecodeShort2N(byteOff);
    case CXBX_VTXFMT_SHORT4N:
        return DecodeShort4N(byteOff);
    case CXBX_VTXFMT_PBYTE4:
        return DecodePByte4(ReadU32(byteOff));
    case CXBX_VTXFMT_FLOAT2H:
        return DecodeFloat2H(byteOff);
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

#endif // CXBX_IA_BYPASS
#endif // CXBX_VERTEX_FETCH_HLSLI
