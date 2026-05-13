// CxbxVertexFetchLayout.hlsli — shared C++ / HLSL header
//
// Defines the vertex fetch layout constant buffer (b1).
// Must stay in sync between C++ (Backend_D3D11_VertexFetch.cpp) and
// HLSL (CxbxVertexFetch.hlsli).
//
// Layout: 8 packed uint header (32 bytes) + uint4[16] attrib array (256 bytes) = 288 bytes.

#ifdef __cplusplus
#pragma once
#include <cstdint>

#define VFL_BEGIN struct VertexFetchLayoutCB {
#define VFL_END   };
#define VFL_UINT(name)            uint32_t name
#define VFL_UINT4_ARRAY(name, n)  uint32_t name[n][4]

#else
// HLSL cbuffer definition.  8 consecutive scalar uints pack into 2 float4
// registers without padding, followed by a uint4 array.
#define VFL_BEGIN cbuffer CxbxVertexLayoutCB : register(b1) {
#define VFL_END   };
#define VFL_UINT(name)            uint name
#define VFL_UINT4_ARRAY(name, n)  uint4 name[n]

#endif

// ============================================================
// Vertex fetch layout constant buffer.
// Field order MUST match between HLSL and C++ — do not reorder.
// ============================================================
VFL_BEGIN
    VFL_UINT(PrimType);       // 0=normal, 1=quad, 2=fan, 3=quadstrip, 4=lineloop
    VFL_UINT(IndexedDraw);    // 0=non-indexed, 1=indexed 16-bit, 2=indexed 32-bit
    VFL_UINT(IndexOffset);    // Byte offset into index data
    VFL_UINT(NumAttribs);     // Number of active vertex attributes (1..16)
    VFL_UINT(NumVerts);       // Original Xbox vertex count (for lineloop wrap)
    VFL_UINT(VertexOffset);   // Added to resolved index (StartVertex or BaseVertexIndex)
    VFL_UINT(WindingCW);      // 1 = CW front face, 0 = CCW front face
    VFL_UINT(Pad7);
    // Per-attribute descriptors (16 × uint4 = 256 bytes)
    // x = `elemOffset` byte offset from start of vertex in the stream
    // y = `stride` bytes per vertex for this stream
    // z = `format` (CXBX_VTXFMT_* constant)
    // w = `streamBase` base byte offset of the stream within g_VtxData
    VFL_UINT4_ARRAY(Attribs, 16);
VFL_END

// Clean up macros
#undef VFL_BEGIN
#undef VFL_END
#undef VFL_UINT
#undef VFL_UINT4_ARRAY

#ifdef __cplusplus
static_assert(sizeof(VertexFetchLayoutCB) == 288, "VertexFetchLayoutCB size mismatch");
#endif
