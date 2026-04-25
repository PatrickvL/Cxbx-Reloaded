// ******************************************************************
// *
// *  This file is part of the Cxbx project.
// *
// *  Cxbx and Cxbe are free software; you can redistribute them
// *  and/or modify them under the terms of the GNU General Public
// *  License as published by the Free Software Foundation; either
// *  version 2 of the license, or (at your option) any later version.
// *
// *  This program is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// *  GNU General Public License for more details.
// *
// *  You should have recieved a copy of the GNU General Public License
// *  along with this program; see the file COPYING.
// *  If not, write to the Free Software Foundation, Inc.,
// *  59 Temple Place - Suite 330, Bostom, MA 02111-1307, USA.
// *
// *  All rights reserved
// *
// ******************************************************************

// Backend_D3D11_IABypass.cpp — Input Assembler bypass draw path.
//
// This module handles vertex data upload
// and draw calls using SV_VertexID-based vertex fetch in the shader.
// The Input Assembler is not used for vertex/index buffer binding.

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h"
#include "common/AddressRanges.h"
#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\XbPushBuffer.h" // HLE_get_NV2A_vertex_attribute_value_pointer
#include "devices\Xbox.h"              // For extern NV2ADevice* g_NV2A
#include "devices\video\nv2a.h"        // For NV2AState, PGRAPHState, VertexAttribute, nv2a_regs.h

// ******************************************************************
// * Format mapping constants (must match CXBX_VTXFMT_* in CxbxVertexFetch.hlsli)
// ******************************************************************
#define CXBX_VTXFMT_FLOAT1       0
#define CXBX_VTXFMT_FLOAT2       1
#define CXBX_VTXFMT_FLOAT3       2
#define CXBX_VTXFMT_FLOAT4       3
#define CXBX_VTXFMT_D3DCOLOR     4
#define CXBX_VTXFMT_SHORT2       5
#define CXBX_VTXFMT_SHORT4       6
#define CXBX_VTXFMT_NORMPACKED3  7
#define CXBX_VTXFMT_SHORT2N      8
#define CXBX_VTXFMT_SHORT4N      9
#define CXBX_VTXFMT_PBYTE4       10
#define CXBX_VTXFMT_FLOAT2H      11
#define CXBX_VTXFMT_NONE         12
#define CXBX_VTXFMT_SHORT1N      13
#define CXBX_VTXFMT_SHORT3N      14
#define CXBX_VTXFMT_PBYTE1       15
#define CXBX_VTXFMT_PBYTE2       16
#define CXBX_VTXFMT_PBYTE3       17
#define CXBX_VTXFMT_SHORT1       18
#define CXBX_VTXFMT_SHORT3       19

// Prim type constants (must match CXBX_PRIM_* in CxbxVertexFetch.hlsli)
#define CXBX_PRIM_NORMAL    0
#define CXBX_PRIM_QUAD      1
#define CXBX_PRIM_FAN       2
#define CXBX_PRIM_QUADSTRIP 3
#define CXBX_PRIM_LINELOOP  4

// ******************************************************************
// * Persistent GPU resources for IA bypass
// ******************************************************************
// UP draw staging buffer: only used for DrawPrimitiveUP / inline vertex data
// where the source pointer is not in the 64 MiB contiguous mirror.
static ID3D11Buffer*             s_pUPVtxDataBuf = nullptr;
static UINT                      s_UPVtxDataBufSize = 0;
static ID3D11ShaderResourceView* s_pUPVtxDataSRV = nullptr;
static ID3D11ShaderResourceView* s_pUPVtxDataSRV_SNORM16x2 = nullptr;
static ID3D11ShaderResourceView* s_pUPVtxDataSRV_UNORM8x4 = nullptr;

// Index data staging buffer: only for indices not in contiguous memory (pushbuffer inline)
static ID3D11Buffer*             s_pIdxDataBuf = nullptr;
static UINT                      s_IdxDataBufSize = 0;
static ID3D11ShaderResourceView* s_pIdxDataSRV = nullptr;

static ID3D11Buffer*             s_pLayoutCB = nullptr;     // Vertex layout CB (b1)
static ID3D11Buffer*             s_pDefaultsCB = nullptr;   // Vertex defaults CB (b2)

// Optimization: cached last-bound GPU pointers to skip redundant API calls
static ID3D11ShaderResourceView* s_pLastBoundVtxSRV = nullptr;
static ID3D11ShaderResourceView* s_pLastBoundIdxSRV = nullptr;
static ID3D11ShaderResourceView* s_pLastBoundSNormSRV = nullptr;
static ID3D11ShaderResourceView* s_pLastBoundUNormSRV = nullptr;
static ID3D11Buffer*             s_pLastBoundLayoutCB = nullptr;
static ID3D11Buffer*             s_pLastBoundDefaultsCB = nullptr;
static bool                      s_IAAlreadyNull = false;   // IA null-binding elimination

// Layout CB caching: generation counter bumped on state changes
static UINT                      s_LayoutCBGeneration = 0;
static UINT                      s_LastLayoutCBGeneration = UINT_MAX;

// ******************************************************************
// * Layout constant buffer structure (must match CxbxVertexLayoutCB in HLSL)
// ******************************************************************
struct IABypassLayoutCB {
	// Header: 8 uints (32 bytes, matches HLSL CxbxVertexLayoutCB)
	UINT PrimType;        // 0=normal, 1=quad, 2=fan, 3=quadstrip, 4=lineloop
	UINT IndexedDraw;     // 0=non-indexed, 1=indexed 16-bit, 2=indexed 32-bit
	UINT IndexOffset;     // Byte offset into index data
	UINT NumAttribs;      // Number of active attributes
	UINT NumVerts;        // Original Xbox vertex count (for lineloop wrap)
	UINT VertexOffset;    // Added to each resolved vertex index before VB fetch
	                      //   Non-indexed: StartVertex
	                      //   Indexed: BaseVertexIndex
	UINT Pad6;
	UINT Pad7;
	UINT Attribs[16][4];  // Per-attribute: elemOffset, stride, format, streamBase
};

// ******************************************************************
// * Map NV2A hardware format + count to CXBX_VTXFMT_* constant
// * format = NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE (bits 3:0)
// * count  = NV097_SET_VERTEX_DATA_ARRAY_FORMAT_SIZE (bits 7:4)
// ******************************************************************
static UINT NV2AFormatToVtxFmt(unsigned format, unsigned count)
{
	if (count == 0) return CXBX_VTXFMT_NONE;

	switch (format) {
	case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_UB_D3D: // 0 — BGRA unsigned byte normalized
		return CXBX_VTXFMT_D3DCOLOR;
	case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_S1:     // 1 — signed short normalized
		switch (count) {
		case 1: return CXBX_VTXFMT_SHORT1N;
		case 2: return CXBX_VTXFMT_SHORT2N;
		case 3: return CXBX_VTXFMT_SHORT3N;
		case 4: return CXBX_VTXFMT_SHORT4N;
		default: return CXBX_VTXFMT_NONE;
		}
	case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_F:      // 2 — float
		switch (count) {
		case 1: return CXBX_VTXFMT_FLOAT1;
		case 2: return CXBX_VTXFMT_FLOAT2;
		case 3: return CXBX_VTXFMT_FLOAT3;
		case 4: return CXBX_VTXFMT_FLOAT4;
		case 7: return CXBX_VTXFMT_FLOAT2H; // Xbox FLOAT2H: 3 floats (x, y, 1/w)
		default: return CXBX_VTXFMT_NONE;
		}
	case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_UB_OGL: // 4 — RGBA unsigned byte normalized
		switch (count) {
		case 1: return CXBX_VTXFMT_PBYTE1;
		case 2: return CXBX_VTXFMT_PBYTE2;
		case 3: return CXBX_VTXFMT_PBYTE3;
		case 4: return CXBX_VTXFMT_PBYTE4;
		default: return CXBX_VTXFMT_NONE;
		}
	case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_S32K:   // 5 — signed short unnormalized
		switch (count) {
		case 1: return CXBX_VTXFMT_SHORT1;
		case 2: return CXBX_VTXFMT_SHORT2;
		case 3: return CXBX_VTXFMT_SHORT3;
		case 4: return CXBX_VTXFMT_SHORT4;
		default: return CXBX_VTXFMT_NONE;
		}
	case NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_CMP:    // 6 — 11.11.10 packed
		return CXBX_VTXFMT_NORMPACKED3;
	default:
		return CXBX_VTXFMT_NONE;
	}
}

// ******************************************************************
// * Map Xbox vertex format to CXBX_VTXFMT_* constant
// ******************************************************************
static UINT XboxFormatToVtxFmt(UINT xboxType)
{
	switch (xboxType) {
	case 0x12: return CXBX_VTXFMT_FLOAT1;       // X_D3DVSDT_FLOAT1
	case 0x22: return CXBX_VTXFMT_FLOAT2;       // X_D3DVSDT_FLOAT2
	case 0x32: return CXBX_VTXFMT_FLOAT3;       // X_D3DVSDT_FLOAT3
	case 0x42: return CXBX_VTXFMT_FLOAT4;       // X_D3DVSDT_FLOAT4
	case 0x40: return CXBX_VTXFMT_D3DCOLOR;     // X_D3DVSDT_D3DCOLOR
	case 0x25: return CXBX_VTXFMT_SHORT2;       // X_D3DVSDT_SHORT2
	case 0x45: return CXBX_VTXFMT_SHORT4;       // X_D3DVSDT_SHORT4
	case 0x16: return CXBX_VTXFMT_NORMPACKED3;  // X_D3DVSDT_NORMPACKED3
	case 0x11: return CXBX_VTXFMT_SHORT1N;       // X_D3DVSDT_NORMSHORT1 (1 short → float4(x,0,0,1))
	case 0x21: return CXBX_VTXFMT_SHORT2N;       // X_D3DVSDT_NORMSHORT2
	case 0x31: return CXBX_VTXFMT_SHORT3N;       // X_D3DVSDT_NORMSHORT3 (3 shorts → float4(x,y,z,1))
	case 0x41: return CXBX_VTXFMT_SHORT4N;       // X_D3DVSDT_NORMSHORT4
	case 0x14: return CXBX_VTXFMT_PBYTE1;        // X_D3DVSDT_PBYTE1 (1 byte → float4(x,0,0,1))
	case 0x24: return CXBX_VTXFMT_PBYTE2;        // X_D3DVSDT_PBYTE2 (2 bytes → float4(x,y,0,1))
	case 0x34: return CXBX_VTXFMT_PBYTE3;        // X_D3DVSDT_PBYTE3 (3 bytes → float4(x,y,z,1))
	case 0x44: return CXBX_VTXFMT_PBYTE4;        // X_D3DVSDT_PBYTE4
	case 0x15: return CXBX_VTXFMT_SHORT1;        // X_D3DVSDT_SHORT1 (1 short unnormalized)
	case 0x35: return CXBX_VTXFMT_SHORT3;        // X_D3DVSDT_SHORT3 (3 shorts unnormalized)
	case 0x72: return CXBX_VTXFMT_FLOAT2H;      // X_D3DVSDT_FLOAT2H
	case 0x02: return CXBX_VTXFMT_NONE;         // X_D3DVSDT_NONE
	default:   return CXBX_VTXFMT_NONE;
	}
}

// ******************************************************************
// * Initialize IA bypass resources (called once during device init)
// ******************************************************************
void CxbxD3D11IABypassInit()
{
	HRESULT hr;

	// Layout CB (b1) — 32 + 256 = 288 bytes
	hr = CxbxD3D11CreateConstantBuffer(sizeof(IABypassLayoutCB), true, &s_pLayoutCB);
	if (FAILED(hr))
		EmuLog(LOG_LEVEL::WARNING, "IABypassInit: Failed to create layout CB");

	// Defaults CB (b2) — 16 × float4 = 256 bytes
	hr = CxbxD3D11CreateConstantBuffer(16 * 4 * sizeof(float), true, &s_pDefaultsCB);
	if (FAILED(hr))
		EmuLog(LOG_LEVEL::WARNING, "IABypassInit: Failed to create defaults CB");
}

// ******************************************************************
// * Release IA bypass resources
// ******************************************************************
void CxbxD3D11IABypassRelease()
{
	if (s_pUPVtxDataSRV_UNORM8x4) { s_pUPVtxDataSRV_UNORM8x4->Release(); s_pUPVtxDataSRV_UNORM8x4 = nullptr; }
	if (s_pUPVtxDataSRV_SNORM16x2) { s_pUPVtxDataSRV_SNORM16x2->Release(); s_pUPVtxDataSRV_SNORM16x2 = nullptr; }
	if (s_pUPVtxDataSRV) { s_pUPVtxDataSRV->Release(); s_pUPVtxDataSRV = nullptr; }
	if (s_pUPVtxDataBuf) { s_pUPVtxDataBuf->Release(); s_pUPVtxDataBuf = nullptr; }
	s_UPVtxDataBufSize = 0;

	if (s_pIdxDataSRV) { s_pIdxDataSRV->Release(); s_pIdxDataSRV = nullptr; }
	if (s_pIdxDataBuf) { s_pIdxDataBuf->Release(); s_pIdxDataBuf = nullptr; }
	s_IdxDataBufSize = 0;

	if (s_pLayoutCB)   { s_pLayoutCB->Release();   s_pLayoutCB = nullptr; }
	if (s_pDefaultsCB) { s_pDefaultsCB->Release(); s_pDefaultsCB = nullptr; }

	s_pLastBoundVtxSRV = nullptr;
	s_pLastBoundIdxSRV = nullptr;
	s_pLastBoundSNormSRV = nullptr;
	s_pLastBoundUNormSRV = nullptr;
	s_pLastBoundLayoutCB = nullptr;
	s_pLastBoundDefaultsCB = nullptr;
	s_IAAlreadyNull = false;
	s_LayoutCBGeneration = 0;
	s_LastLayoutCBGeneration = UINT_MAX;
}

// Called externally when SetStreamSource or SetVertexShader change
void CxbxD3D11IABypassInvalidateLayout()
{
	s_LayoutCBGeneration++;
}

// ******************************************************************
// * Ensure UP vertex data buffer is large enough
// ******************************************************************
static void EnsureUPVtxDataBuffer(UINT requiredSize)
{
	UINT oldSize = s_UPVtxDataBufSize;
	CxbxD3D11EnsureRawStagingBuffer(requiredSize,
		&s_pUPVtxDataBuf, &s_UPVtxDataBufSize,
		&s_pUPVtxDataSRV, "IABypass_UPVtxData");

	// If the buffer was (re)created, also create typed SRV views for hardware format decode
	if (s_UPVtxDataBufSize != oldSize && s_pUPVtxDataBuf) {
		if (s_pUPVtxDataSRV_SNORM16x2) { s_pUPVtxDataSRV_SNORM16x2->Release(); s_pUPVtxDataSRV_SNORM16x2 = nullptr; }
		if (s_pUPVtxDataSRV_UNORM8x4)  { s_pUPVtxDataSRV_UNORM8x4->Release();  s_pUPVtxDataSRV_UNORM8x4 = nullptr; }

		D3D11_SHADER_RESOURCE_VIEW_DESC typedDesc = {};
		typedDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
		typedDesc.Buffer.FirstElement = 0;
		typedDesc.Buffer.NumElements = s_UPVtxDataBufSize / 4;

		typedDesc.Format = DXGI_FORMAT_R16G16_SNORM;
		g_pD3DDevice->CreateShaderResourceView(s_pUPVtxDataBuf, &typedDesc, &s_pUPVtxDataSRV_SNORM16x2);

		typedDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		g_pD3DDevice->CreateShaderResourceView(s_pUPVtxDataBuf, &typedDesc, &s_pUPVtxDataSRV_UNORM8x4);
	}
}

// ******************************************************************
// * Ensure index data buffer is large enough
// ******************************************************************
static void EnsureIdxDataBuffer(UINT requiredSize)
{
	CxbxD3D11EnsureRawStagingBuffer(requiredSize,
		&s_pIdxDataBuf, &s_IdxDataBufSize,
		&s_pIdxDataSRV, "IABypass_IdxData");
}

// ******************************************************************
// * Upload vertex defaults (NV2A sticky attribute values) to CB b2
// ******************************************************************
// Dirty flag for vertex defaults — set by CxbxSetVertexAttribute, consumed here
bool g_bD3D11IABypassDefaultsDirty = true;

static void UploadVertexDefaults()
{
	if (!s_pDefaultsCB) return;
	if (!g_bD3D11IABypassDefaultsDirty) return;

	g_bD3D11IABypassDefaultsDirty = false;

	D3D11_MAPPED_SUBRESOURCE mapped = {};
	HRESULT hr = g_pD3DDeviceContext->Map(s_pDefaultsCB, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
	if (FAILED(hr)) return;

	float* pDst = (float*)mapped.pData;
	for (int i = 0; i < 16; i++) {
		const float* pSrc = HLE_get_NV2A_vertex_attribute_value_pointer(i);
		pDst[i * 4 + 0] = pSrc[0];
		pDst[i * 4 + 1] = pSrc[1];
		pDst[i * 4 + 2] = pSrc[2];
		pDst[i * 4 + 3] = pSrc[3];
	}

	g_pD3DDeviceContext->Unmap(s_pDefaultsCB, 0);
}

// ******************************************************************
// * Core draw function for IA bypass
// * Returns true if the draw was handled, false to fall back to IA path.
// ******************************************************************
void CxbxD3D11IABypassDraw(CxbxDrawContext& DrawContext)
{
	// When all vertex shaders are compiled with IA bypass, the normal IA
	// fallback path cannot work (shader expects SV_VertexID, not TEXCOORD
	// inputs).
	if (!s_pLayoutCB || !s_pDefaultsCB)
		return;

	CxbxVertexDeclaration* pDecl = CxbxGetVertexDeclaration();
	if (!pDecl || pDecl->NumberOfVertexStreams == 0)
		return;

	// ---------------------------------------------------------------
	// Step 1: Determine topology and host vertex count
	// ---------------------------------------------------------------
	UINT primType = CXBX_PRIM_NORMAL;
	UINT hostVertexCount = DrawContext.dwVertexCount;
	D3D_PRIMITIVE_TOPOLOGY hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;

	switch (DrawContext.XboxPrimitiveType) {
	case xbox::X_D3DPT_QUADLIST:
		primType = CXBX_PRIM_QUAD;
		hostVertexCount = (DrawContext.dwVertexCount / 4) * 6;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_QUADSTRIP:
		primType = CXBX_PRIM_QUADSTRIP;
		// Each pair of vertices adds a quad (2 triangles) after the first 2 vertices
		hostVertexCount = (DrawContext.dwVertexCount >= 4) ? ((DrawContext.dwVertexCount - 2) / 2) * 6 : 0;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_TRIANGLEFAN:
	case xbox::X_D3DPT_POLYGON:
		primType = CXBX_PRIM_FAN;
		hostVertexCount = (DrawContext.dwVertexCount >= 3) ? (DrawContext.dwVertexCount - 2) * 3 : 0;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_TRIANGLELIST:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_TRIANGLESTRIP:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP;
		break;
	case xbox::X_D3DPT_LINELIST:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		break;
	case xbox::X_D3DPT_LINESTRIP:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_LINESTRIP;
		break;
	case xbox::X_D3DPT_POINTLIST:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_POINTLIST;
		break;
	case xbox::X_D3DPT_LINELOOP:
		primType = CXBX_PRIM_LINELOOP;
		// N vertices → N line segments → 2N host vertices as LINELIST
		hostVertexCount = DrawContext.dwVertexCount * 2;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		break;
	default:
		// Unsupported topology — skip draw (can't fall back, shader expects SV_VertexID)
		return;
	}

	if (hostVertexCount == 0)
		return; // Nothing to draw — handled

	// ---------------------------------------------------------------
	// Step 2: Determine data source — mirror (VB draws) or staging (UP draws)
	// ---------------------------------------------------------------
	bool bIsUPDraw = (DrawContext.pXboxVertexStreamZeroData != nullptr);

	// The page-tracked 64 MiB mirror covers all Xbox VBs/IBs in contiguous memory.
	// UP draws use a per-draw staging buffer since data comes from arbitrary pointers.
	ID3D11ShaderResourceView* pMirrorSRV = CxbxPageTrackerGetMirrorSRV();

	// Flush dirty pages so the GPU mirror is current (for VB draws) and
	// s_TextureDirtyBitmap is updated (for texture re-upload detection).
	// Must run for both UP and non-UP draws: UP draws skip the mirror
	// but still need texture dirty tracking via GetWriteWatch().
	CxbxPageTrackerFlushToGPU();

	// ---------------------------------------------------------------
	// Step 2b: Upload UP vertex data to staging buffer
	// ---------------------------------------------------------------
	if (bIsUPDraw) {
		UINT stride = DrawContext.uiXboxVertexStreamZeroStride;
		UINT vtxDataSize = DrawContext.dwVertexCount * stride;
		vtxDataSize = (vtxDataSize + 3) & ~3u; // Align to 4 bytes

		EnsureUPVtxDataBuffer(vtxDataSize);
		if (!s_pUPVtxDataBuf || !s_pUPVtxDataSRV)
			return;

		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pUPVtxDataBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
		if (FAILED(hr)) return;

		memcpy(mapped.pData,
			(const uint8_t*)DrawContext.pXboxVertexStreamZeroData
				+ DrawContext.dwStartVertex * stride,
			vtxDataSize);

		g_pD3DDeviceContext->Unmap(s_pUPVtxDataBuf, 0);
	}

	UINT vertexStart = DrawContext.dwStartVertex;
	UINT numVertices = DrawContext.dwVertexCount;

	// ---------------------------------------------------------------
	// Step 3: Resolve index data (if indexed draw)
	// ---------------------------------------------------------------
	UINT indexedDraw = 0;
	UINT indexOffset = 0;
	bool bIdxFromMirror = false; // true = index data served from the 64 MiB mirror (t0)

	if (DrawContext.pXboxIndexData) {
		indexedDraw = 1; // 16-bit indices

		// Check if index data pointer is in contiguous memory (0x80000000 range).
		// If so, we can read it directly from the mirror buffer — no upload needed.
		uintptr_t idxAddr = (uintptr_t)DrawContext.pXboxIndexData;
		if (idxAddr >= CONTIGUOUS_MEMORY_BASE
			&& idxAddr < (CONTIGUOUS_MEMORY_BASE + XBOX_CONTIGUOUS_MEMORY_SIZE)) {
			// Index data is in the mirror — just pass the byte offset.
			// We bind the mirror SRV as t1 (g_IdxData); the shader reads at g_IndexOffset.
			indexOffset = (UINT)(idxAddr - CONTIGUOUS_MEMORY_BASE);
			bIdxFromMirror = true;
		} else {
			// Index data is NOT in contiguous memory (e.g., pushbuffer inline data).
			// Upload to the per-draw index buffer as before.
			UINT idxDataSize = DrawContext.dwVertexCount * sizeof(INDEX16);
			idxDataSize = (idxDataSize + 3) & ~3u;

			EnsureIdxDataBuffer(idxDataSize);
			if (!s_pIdxDataBuf || !s_pIdxDataSRV)
				return;

			HRESULT hr = CxbxD3D11UpdateDynamicBuffer(s_pIdxDataBuf, DrawContext.pXboxIndexData, DrawContext.dwVertexCount * sizeof(INDEX16));
			if (FAILED(hr)) return;

			indexOffset = 0;
		}
	}

	// ---------------------------------------------------------------
	// Step 4: Fill layout constant buffer (skip if generation unchanged)
	// ---------------------------------------------------------------
	// Bump generation for prim-type or index-mode changes (cheap inline check)
	// The generation is also bumped externally by CxbxD3D11IABypassInvalidateLayout()
	// for SetStreamSource / SetVertexShader changes.
	{
		// Build a local hash of fields that change per-draw but aren't covered by
		// the external invalidation (prim type, indexed mode, vertex range, index offset)
		UINT vertexOffset = indexedDraw ? DrawContext.dwBaseVertexIndex : vertexStart;
		UINT drawLocalKey = primType | (indexedDraw << 2) | (vertexOffset << 4) | (numVertices << 20);
		static UINT s_LastDrawLocalKey = UINT_MAX;
		static UINT s_LastIndexOffset = UINT_MAX;
		bool layoutDirty = (s_LayoutCBGeneration != s_LastLayoutCBGeneration)
		                || (drawLocalKey != s_LastDrawLocalKey)
		                || (indexOffset != s_LastIndexOffset);
		s_LastLayoutCBGeneration = s_LayoutCBGeneration;
		s_LastDrawLocalKey = drawLocalKey;
		s_LastIndexOffset = indexOffset;

		if (!layoutDirty) goto skip_layout_upload;
	}
	{
		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pLayoutCB, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
		if (FAILED(hr)) return;

		IABypassLayoutCB* pCB = (IABypassLayoutCB*)mapped.pData;
		memset(pCB, 0, sizeof(IABypassLayoutCB));

		pCB->PrimType = primType;
		pCB->IndexedDraw = indexedDraw;
		pCB->IndexOffset = indexOffset;
		pCB->NumAttribs = 16; // Always provide all 16 attribute descriptors
		pCB->NumVerts = DrawContext.dwVertexCount; // Original vertex count (for lineloop)

		// VertexOffset: adjusts the resolved vertex index before VB fetch.
		// Non-indexed draws: StartVertex (DrawVertices skips the first N vertices).
		// Indexed draws: BaseVertexIndex (SetIndices offset added to each index).
		if (indexedDraw)
			pCB->VertexOffset = DrawContext.dwBaseVertexIndex;
		else
			pCB->VertexOffset = vertexStart;

		// Fill per-attribute descriptors — default all to NONE (use sticky defaults)
		for (UINT a = 0; a < 16; a++) {
			pCB->Attribs[a][0] = 0;  // elemOffset
			pCB->Attribs[a][1] = 0;  // stride
			pCB->Attribs[a][2] = CXBX_VTXFMT_NONE; // format — default is "use default value"
			pCB->Attribs[a][3] = 0;  // streamBase
		}

		// Two vertex layout paths:
		// 1. PGRAPH path: reads vertex_attributes[] directly from NV2A state.
		//    Used for push buffer draws (HLE_draw_arrays) where PGRAPH is the
		//    authoritative source and g_Xbox_SetStreamSource[] is not populated.
		// 2. HLE path: reads CxbxVertexDeclaration + g_Xbox_SetStreamSource[].
		//    Used for HLE-intercepted draws where SetStreamSource patches populate
		//    the HLE state before DrawPrimitive is called.
		//
		// Strategy: try PGRAPH first for non-UP draws; if no active attributes
		// found, fall back to HLE path.
		bool bUsedPGRAPH = false;
		PGRAPHState* pg = (g_NV2A != nullptr) ? &g_NV2A->GetDeviceState()->pgraph : nullptr;

		if (pg && !bIsUPDraw) {
			// PGRAPH path: slot index = register index, offset = physical address
			for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
				const VertexAttribute& attr = pg->vertex_attributes[i];
				if (attr.count == 0) continue; // inactive attribute → use default

				pCB->Attribs[i][0] = 0;           // elemOffset (baked into offset)
				pCB->Attribs[i][1] = attr.stride;
				pCB->Attribs[i][2] = NV2AFormatToVtxFmt(attr.format, attr.count);
				pCB->Attribs[i][3] = (UINT)attr.offset; // physical addr = SRV byte offset
				bUsedPGRAPH = true;
			}
		}

		if (!bUsedPGRAPH)
		{
			// HLE fallback: walk CxbxVertexDeclaration + g_Xbox_SetStreamSource[]
			for (UINT s = 0; s < pDecl->NumberOfVertexStreams; s++) {
				auto& streamInfo = pDecl->VertexStreams[s];
				UINT streamIdx = streamInfo.XboxStreamIndex;
				auto& streamInput = g_Xbox_SetStreamSource[streamIdx];

				UINT stride;
				if (s == 0 && bIsUPDraw) {
					stride = DrawContext.uiXboxVertexStreamZeroStride;
				} else {
					stride = streamInput.Stride;
					if (stride == 0) stride = streamInfo.HostVertexStride;
				}

				UINT elemOffset = 0;
				UINT hostElemOffset = 0;
				for (UINT e = 0; e < streamInfo.NumberOfVertexElements; e++) {
					auto& elem = streamInfo.VertexElements[e];
					if (elem.XboxType == 0)
						continue;

					UINT regIdx = 0;
					bool found = false;
					if (pDecl->pD3D11InputElements) {
						for (UINT ie = 0; ie < pDecl->D3D11InputElementCount; ie++) {
							auto& inputElem = pDecl->pD3D11InputElements[ie];
							if (inputElem.InputSlot == streamIdx
								&& inputElem.AlignedByteOffset == hostElemOffset) {
								regIdx = inputElem.SemanticIndex;
								found = true;
								break;
							}
						}
					}

					if (found && regIdx < 16) {
						INT streamBase;
						if (bIsUPDraw && s == 0) {
							streamBase = -(INT)vertexStart * (INT)stride;
						} else if (streamInput.VertexBuffer) {
							uintptr_t vbAddr = (uintptr_t)GetDataFromXboxResource(streamInput.VertexBuffer);
							streamBase = (INT)(vbAddr - CONTIGUOUS_MEMORY_BASE) + (INT)streamInput.Offset;
						} else {
							streamBase = 0;
						}

						pCB->Attribs[regIdx][0] = elemOffset;
						pCB->Attribs[regIdx][1] = stride;
						pCB->Attribs[regIdx][2] = XboxFormatToVtxFmt(elem.XboxType);
						pCB->Attribs[regIdx][3] = (UINT)streamBase;
					}

					elemOffset += elem.XboxByteSize;
					hostElemOffset += elem.HostByteSize;
				}
			}
		}

		g_pD3DDeviceContext->Unmap(s_pLayoutCB, 0);
	}
skip_layout_upload:

	// ---------------------------------------------------------------
	// Step 5: Upload vertex defaults (skip if not dirty)
	// ---------------------------------------------------------------
	UploadVertexDefaults();

	// ---------------------------------------------------------------
	// Step 6: Bind resources and issue draw
	// ---------------------------------------------------------------

	// Unbind IA state — null input layout, null vertex/index buffers
	// (skip if already nulled from a prior IA bypass draw)
	if (!s_IAAlreadyNull) {
		g_pD3DDeviceContext->IASetInputLayout(nullptr);
		ID3D11Buffer* nullBufs[17] = {};
		UINT nullStrides[17] = {};
		UINT nullOffsets[17] = {};
		g_pD3DDeviceContext->IASetVertexBuffers(0, 17, nullBufs, nullStrides, nullOffsets);
		g_pD3DDeviceContext->IASetIndexBuffer(nullptr, DXGI_FORMAT_R16_UINT, 0);
		s_IAAlreadyNull = true;
	}
	g_pD3DDeviceContext->IASetPrimitiveTopology(hostTopology);

	// Bind SRVs to VS: t0 = vertex data (raw), t1 = index data,
	// t2 = vertex data (R16G16_SNORM), t3 = vertex data (R8G8B8A8_UNORM)
	// VB draws: SRVs from the 64 MiB page-tracked mirror.
	// UP draws: SRVs from the per-draw staging buffer.
	ID3D11ShaderResourceView* pActiveVtxSRV = bIsUPDraw ? s_pUPVtxDataSRV : pMirrorSRV;
	ID3D11ShaderResourceView* pActiveIdxSRV = bIdxFromMirror ? pMirrorSRV : s_pIdxDataSRV;
	ID3D11ShaderResourceView* pActiveSNormSRV = bIsUPDraw
		? s_pUPVtxDataSRV_SNORM16x2 : CxbxPageTrackerGetMirrorSRV_SNORM16x2();
	ID3D11ShaderResourceView* pActiveUNormSRV = bIsUPDraw
		? s_pUPVtxDataSRV_UNORM8x4 : CxbxPageTrackerGetMirrorSRV_UNORM8x4();

	if (pActiveVtxSRV != s_pLastBoundVtxSRV || pActiveIdxSRV != s_pLastBoundIdxSRV
		|| pActiveSNormSRV != s_pLastBoundSNormSRV || pActiveUNormSRV != s_pLastBoundUNormSRV) {
		ID3D11ShaderResourceView* vsSRVs[4] = { pActiveVtxSRV, pActiveIdxSRV, pActiveSNormSRV, pActiveUNormSRV };
		g_pD3DDeviceContext->VSSetShaderResources(0, 4, vsSRVs);
		s_pLastBoundVtxSRV = pActiveVtxSRV;
		s_pLastBoundIdxSRV = pActiveIdxSRV;
		s_pLastBoundSNormSRV = pActiveSNormSRV;
		s_pLastBoundUNormSRV = pActiveUNormSRV;
	}

	// Bind CBs to VS: b1 = layout, b2 = defaults (skip if unchanged)
	// (b0 is already bound for VS constants by CxbxD3D11FlushVertexShaderConstants)
	if (s_pLayoutCB != s_pLastBoundLayoutCB || s_pDefaultsCB != s_pLastBoundDefaultsCB) {
		ID3D11Buffer* vsCBs[2] = { s_pLayoutCB, s_pDefaultsCB };
		g_pD3DDeviceContext->VSSetConstantBuffers(1, 2, vsCBs);
		s_pLastBoundLayoutCB = s_pLayoutCB;
		s_pLastBoundDefaultsCB = s_pDefaultsCB;
	}

	// Bind thick line GS if needed (only for line primitives that aren't topology-converted)
	if (primType == CXBX_PRIM_NORMAL) {
		CxbxBindThickLineGS(DrawContext.XboxPrimitiveType);
	}

	// Issue the draw
	g_pD3DDeviceContext->Draw(hostVertexCount, 0);

	if (primType == CXBX_PRIM_NORMAL) {
		CxbxUnbindThickLineGS(DrawContext.XboxPrimitiveType);
	}

	// Unbind VS SRVs to avoid conflicts with other passes
	ID3D11ShaderResourceView* nullSRVs[4] = { nullptr, nullptr, nullptr, nullptr };
	g_pD3DDeviceContext->VSSetShaderResources(0, 4, nullSRVs);
	s_pLastBoundVtxSRV = nullptr;
	s_pLastBoundIdxSRV = nullptr;
	s_pLastBoundSNormSRV = nullptr;
	s_pLastBoundUNormSRV = nullptr;

	return;
}

// ******************************************************************
// * Draw inline buffer vertices (Begin/SetVertexData/End path)
// *
// * Inline buffer data is stored as float4 per attribute per vertex
// * in pg->vertex_attributes[i].inline_buffer. This bypasses the
// * normal vertex declaration/stream layout entirely — all 16
// * attributes are packed as FLOAT4 at a fixed 256-byte stride.
// * After drawing, the inline_buffer arrays are freed (same protocol
// * as xemu's pgraph_draw_inline_buffer).
// ******************************************************************
void CxbxD3D11DrawInlineBuffer(PGRAPHState* pg)
{
	if (!s_pLayoutCB || !s_pDefaultsCB)
		return;

	unsigned int vertexCount = pg->inline_buffer_length;
	if (vertexCount == 0) return;

	// ---------------------------------------------------------------
	// Step 1: Pack inline buffer data into contiguous UP vertex buffer
	// ---------------------------------------------------------------
	// All 16 attributes stored as float4 (16 bytes each) = 256 bytes per vertex
	const UINT kAttrSize = 4 * sizeof(float);  // 16 bytes
	const UINT kStride = NV2A_VERTEXSHADER_ATTRIBUTES * kAttrSize; // 256 bytes
	UINT totalSize = vertexCount * kStride;

	EnsureUPVtxDataBuffer(totalSize);
	if (!s_pUPVtxDataBuf || !s_pUPVtxDataSRV)
		return;

	{
		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pUPVtxDataBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
		if (FAILED(hr)) return;

		float* pDst = (float*)mapped.pData;
		for (unsigned int v = 0; v < vertexCount; v++) {
			for (int a = 0; a < NV2A_VERTEXSHADER_ATTRIBUTES; a++) {
				const VertexAttribute& attr = pg->vertex_attributes[a];
				const float* pSrc = attr.inline_buffer
					? &attr.inline_buffer[v * 4]
					: attr.inline_value;
				pDst[0] = pSrc[0];
				pDst[1] = pSrc[1];
				pDst[2] = pSrc[2];
				pDst[3] = pSrc[3];
				pDst += 4;
			}
		}

		g_pD3DDeviceContext->Unmap(s_pUPVtxDataBuf, 0);
	}

	// ---------------------------------------------------------------
	// Step 2: Determine topology and host vertex count
	// ---------------------------------------------------------------
	UINT primType = CXBX_PRIM_NORMAL;
	UINT hostVertexCount = vertexCount;
	D3D_PRIMITIVE_TOPOLOGY hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;

	switch ((xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode) {
	case xbox::X_D3DPT_QUADLIST:
		primType = CXBX_PRIM_QUAD;
		hostVertexCount = (vertexCount / 4) * 6;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_QUADSTRIP:
		primType = CXBX_PRIM_QUADSTRIP;
		hostVertexCount = (vertexCount >= 4) ? ((vertexCount - 2) / 2) * 6 : 0;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_TRIANGLEFAN:
	case xbox::X_D3DPT_POLYGON:
		primType = CXBX_PRIM_FAN;
		hostVertexCount = (vertexCount >= 3) ? (vertexCount - 2) * 3 : 0;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_TRIANGLELIST:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
		break;
	case xbox::X_D3DPT_TRIANGLESTRIP:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP;
		break;
	case xbox::X_D3DPT_LINELIST:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		break;
	case xbox::X_D3DPT_LINESTRIP:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_LINESTRIP;
		break;
	case xbox::X_D3DPT_POINTLIST:
		primType = CXBX_PRIM_NORMAL;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_POINTLIST;
		break;
	case xbox::X_D3DPT_LINELOOP:
		primType = CXBX_PRIM_LINELOOP;
		hostVertexCount = vertexCount * 2;
		hostTopology = D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		break;
	default:
		return;
	}

	if (hostVertexCount == 0)
		return;

	// ---------------------------------------------------------------
	// Step 3: Fill layout CB with float4 layout for all 16 attributes
	// ---------------------------------------------------------------
	{
		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pLayoutCB, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
		if (FAILED(hr)) return;

		IABypassLayoutCB* pCB = (IABypassLayoutCB*)mapped.pData;
		memset(pCB, 0, sizeof(IABypassLayoutCB));

		pCB->PrimType = primType;
		pCB->IndexedDraw = 0;
		pCB->IndexOffset = 0;
		pCB->NumAttribs = 16;
		pCB->NumVerts = vertexCount;
		pCB->VertexOffset = 0;

		for (UINT a = 0; a < 16; a++) {
			pCB->Attribs[a][0] = a * kAttrSize;       // elemOffset
			pCB->Attribs[a][1] = kStride;             // stride
			pCB->Attribs[a][2] = CXBX_VTXFMT_FLOAT4;  // format
			pCB->Attribs[a][3] = 0;                   // streamBase (UP data at offset 0)
		}

		g_pD3DDeviceContext->Unmap(s_pLayoutCB, 0);
	}

	// Invalidate layout cache so the next regular draw refills the CB
	s_LastLayoutCBGeneration = UINT_MAX;

	// ---------------------------------------------------------------
	// Step 4: Upload vertex defaults
	// ---------------------------------------------------------------
	UploadVertexDefaults();

	// ---------------------------------------------------------------
	// Step 5: Bind resources and issue draw
	// ---------------------------------------------------------------
	if (!s_IAAlreadyNull) {
		g_pD3DDeviceContext->IASetInputLayout(nullptr);
		ID3D11Buffer* nullBufs[17] = {};
		UINT nullStrides[17] = {};
		UINT nullOffsets[17] = {};
		g_pD3DDeviceContext->IASetVertexBuffers(0, 17, nullBufs, nullStrides, nullOffsets);
		g_pD3DDeviceContext->IASetIndexBuffer(nullptr, DXGI_FORMAT_R16_UINT, 0);
		s_IAAlreadyNull = true;
	}
	g_pD3DDeviceContext->IASetPrimitiveTopology(hostTopology);

	// Bind UP staging SRVs to VS (t0=raw, t1=null, t2=snorm, t3=unorm)
	{
		ID3D11ShaderResourceView* vsSRVs[4] = {
			s_pUPVtxDataSRV, nullptr,
			s_pUPVtxDataSRV_SNORM16x2, s_pUPVtxDataSRV_UNORM8x4
		};
		g_pD3DDeviceContext->VSSetShaderResources(0, 4, vsSRVs);
		s_pLastBoundVtxSRV = s_pUPVtxDataSRV;
		s_pLastBoundIdxSRV = nullptr;
		s_pLastBoundSNormSRV = s_pUPVtxDataSRV_SNORM16x2;
		s_pLastBoundUNormSRV = s_pUPVtxDataSRV_UNORM8x4;
	}

	// Bind CBs: b1=layout, b2=defaults
	{
		ID3D11Buffer* vsCBs[2] = { s_pLayoutCB, s_pDefaultsCB };
		g_pD3DDeviceContext->VSSetConstantBuffers(1, 2, vsCBs);
		s_pLastBoundLayoutCB = s_pLayoutCB;
		s_pLastBoundDefaultsCB = s_pDefaultsCB;
	}

	if (primType == CXBX_PRIM_NORMAL) {
		CxbxBindThickLineGS((xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode);
	}

	g_pD3DDeviceContext->Draw(hostVertexCount, 0);

	if (primType == CXBX_PRIM_NORMAL) {
		CxbxUnbindThickLineGS((xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode);
	}

	// Unbind VS SRVs
	ID3D11ShaderResourceView* nullSRVs[4] = { nullptr, nullptr, nullptr, nullptr };
	g_pD3DDeviceContext->VSSetShaderResources(0, 4, nullSRVs);
	s_pLastBoundVtxSRV = nullptr;
	s_pLastBoundIdxSRV = nullptr;
	s_pLastBoundSNormSRV = nullptr;
	s_pLastBoundUNormSRV = nullptr;

	// Free per-attribute inline buffers (same protocol as xemu)
	for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
		VertexAttribute& attr = pg->vertex_attributes[i];
		if (attr.inline_buffer) {
			free(attr.inline_buffer);
			attr.inline_buffer = nullptr;
		}
	}
}

