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
// When g_bD3D11IABypass is true, this module handles vertex data upload
// and draw calls using SV_VertexID-based vertex fetch in the shader.
// The Input Assembler is not used for vertex/index buffer binding.

#ifdef CXBX_USE_D3D11

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h"
#include "common/AddressRanges.h"
#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\XbPushBuffer.h" // HLE_get_NV2A_vertex_attribute_value_pointer
#include "../WalkIndexBuffer.h"

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
static ID3D11Buffer*             s_pVtxDataBuf = nullptr;   // ByteAddressBuffer for vertex data
static UINT                      s_VtxDataBufSize = 0;
static ID3D11ShaderResourceView* s_pVtxDataSRV = nullptr;
static ID3D11ShaderResourceView* s_pVtxDataSRV_SNORM16x2 = nullptr; // R16G16_SNORM typed view
static ID3D11ShaderResourceView* s_pVtxDataSRV_UNORM8x4 = nullptr;  // R8G8B8A8_UNORM typed view

static ID3D11Buffer*             s_pIdxDataBuf = nullptr;   // ByteAddressBuffer for index data
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
	UINT Pad5;
	UINT Pad6;
	UINT Pad7;
	UINT Attribs[16][4];  // Per-attribute: elemOffset, stride, format, streamBase
};

// ******************************************************************
// * Map Xbox vertex format to CXBX_VTXFMT_* constant
// ******************************************************************
static UINT XboxFormatToVtxFmt(UINT xboxType)
{
	switch (xboxType) {
	case 0x02: return CXBX_VTXFMT_FLOAT1;       // X_D3DVSDT_FLOAT1
	case 0x12: return CXBX_VTXFMT_FLOAT2;       // X_D3DVSDT_FLOAT2
	case 0x22: return CXBX_VTXFMT_FLOAT3;       // X_D3DVSDT_FLOAT3
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
	case 0:    return CXBX_VTXFMT_NONE;         // X_D3DVSDT_NONE
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
	if (s_pVtxDataSRV_UNORM8x4) { s_pVtxDataSRV_UNORM8x4->Release(); s_pVtxDataSRV_UNORM8x4 = nullptr; }
	if (s_pVtxDataSRV_SNORM16x2) { s_pVtxDataSRV_SNORM16x2->Release(); s_pVtxDataSRV_SNORM16x2 = nullptr; }
	if (s_pVtxDataSRV) { s_pVtxDataSRV->Release(); s_pVtxDataSRV = nullptr; }
	if (s_pVtxDataBuf) { s_pVtxDataBuf->Release(); s_pVtxDataBuf = nullptr; }
	s_VtxDataBufSize = 0;

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
// * Ensure vertex data buffer is large enough
// ******************************************************************
static void EnsureVtxDataBuffer(UINT requiredSize)
{
	UINT oldSize = s_VtxDataBufSize;
	CxbxD3D11EnsureRawStagingBuffer(requiredSize,
		&s_pVtxDataBuf, &s_VtxDataBufSize,
		&s_pVtxDataSRV, "IABypass_VtxData");

	// If the buffer was (re)created, also create typed SRV views for hardware format decode
	if (s_VtxDataBufSize != oldSize && s_pVtxDataBuf) {
		if (s_pVtxDataSRV_SNORM16x2) { s_pVtxDataSRV_SNORM16x2->Release(); s_pVtxDataSRV_SNORM16x2 = nullptr; }
		if (s_pVtxDataSRV_UNORM8x4)  { s_pVtxDataSRV_UNORM8x4->Release();  s_pVtxDataSRV_UNORM8x4 = nullptr; }

		D3D11_SHADER_RESOURCE_VIEW_DESC typedDesc = {};
		typedDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
		typedDesc.Buffer.FirstElement = 0;
		typedDesc.Buffer.NumElements = s_VtxDataBufSize / 4;

		typedDesc.Format = DXGI_FORMAT_R16G16_SNORM;
		g_pD3DDevice->CreateShaderResourceView(s_pVtxDataBuf, &typedDesc, &s_pVtxDataSRV_SNORM16x2);

		typedDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		g_pD3DDevice->CreateShaderResourceView(s_pVtxDataBuf, &typedDesc, &s_pVtxDataSRV_UNORM8x4);
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
bool CxbxD3D11IABypassDraw(CxbxDrawContext& DrawContext)
{
	// When all vertex shaders are compiled with IA bypass, the normal IA
	// fallback path cannot work (shader expects SV_VertexID, not TEXCOORD
	// inputs).  Return true ("handled") to skip the incompatible fallback
	// even if we can't actually draw.
	if (!s_pLayoutCB || !s_pDefaultsCB)
		return true;

	CxbxVertexDeclaration* pDecl = CxbxGetVertexDeclaration();
	if (!pDecl || pDecl->NumberOfVertexStreams == 0)
		return true;

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
		return true;
	}

	if (hostVertexCount == 0)
		return true; // Nothing to draw — handled

	// ---------------------------------------------------------------
	// Step 2: Ensure GPU mirror has latest CPU writes (page tracker path)
	// ---------------------------------------------------------------
	// Flush any pages the game has written since last draw.
	// After this, the 64 MiB mirror ByteAddressBuffer is up-to-date.
	// The mirror path can only be used for VB draws (not UP draws where
	// vertex data comes from an arbitrary user pointer on the stack/heap).
	ID3D11ShaderResourceView* pMirrorSRV = CxbxPageTrackerGetMirrorSRV();
	bool bUsingMirror = (pMirrorSRV != nullptr) && !DrawContext.pXboxVertexStreamZeroData;

	if (bUsingMirror) {
		CxbxPageTrackerFlushToGPU();
	}

	// ---------------------------------------------------------------
	// Step 2b: Upload vertex data (fallback when mirror is not available)
	// ---------------------------------------------------------------
	// Calculate total data size needed across all streams
	UINT streamOffsets[X_VSH_MAX_STREAMS] = {};  // byte offset of each stream within g_VtxData
	UINT totalVtxDataSize = 0;

	// Determine the vertex range to upload
	UINT vertexStart = DrawContext.dwStartVertex;
	UINT numVertices = DrawContext.dwVertexCount;
	if (!bUsingMirror && DrawContext.pXboxIndexData) {
		// For indexed draws in fallback path, we need the range [LowIndex..HighIndex]
		// to know what subset of VB data to upload. (Mirror path doesn't need this —
		// the entire 64 MiB is already available.)
		if (DrawContext.HighIndex == 0) {
			WalkIndexBuffer(DrawContext.LowIndex, DrawContext.HighIndex,
				DrawContext.pXboxIndexData, DrawContext.dwVertexCount);
		}
		vertexStart = DrawContext.LowIndex;
		numVertices = DrawContext.HighIndex - DrawContext.LowIndex + 1;
	}

	bool bIsUPDraw = (DrawContext.pXboxVertexStreamZeroData != nullptr);

	if (!bUsingMirror) {
		for (UINT s = 0; s < pDecl->NumberOfVertexStreams; s++) {
			auto& streamInfo = pDecl->VertexStreams[s];
			UINT streamIdx = streamInfo.XboxStreamIndex;
			auto& streamInput = g_Xbox_SetStreamSource[streamIdx];

			if (!streamInput.VertexBuffer && !bIsUPDraw)
				continue;

			UINT stride;
			if (s == 0 && bIsUPDraw) {
				stride = DrawContext.uiXboxVertexStreamZeroStride;
			} else {
				stride = streamInput.Stride;
				if (stride == 0) stride = streamInfo.HostVertexStride; // Fallback
			}

			UINT streamDataSize = numVertices * stride;
			streamDataSize = (streamDataSize + 3) & ~3u; // Align to 4 bytes

			streamOffsets[s] = totalVtxDataSize;
			totalVtxDataSize += streamDataSize;
		}

		if (totalVtxDataSize == 0)
			return true;

		EnsureVtxDataBuffer(totalVtxDataSize);
		if (!s_pVtxDataBuf || !s_pVtxDataSRV)
			return true;

		// Upload all stream data into the single buffer
		{
			D3D11_MAPPED_SUBRESOURCE mapped = {};
			HRESULT hr = g_pD3DDeviceContext->Map(s_pVtxDataBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
			if (FAILED(hr)) return true;

			uint8_t* pDst = (uint8_t*)mapped.pData;

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

				UINT streamDataSize = numVertices * stride;
				const uint8_t* pSrc = nullptr;

				if (s == 0 && bIsUPDraw) {
					// UP draw: vertex data from user pointer
					pSrc = (const uint8_t*)DrawContext.pXboxVertexStreamZeroData
						+ vertexStart * stride;
				} else if (streamInput.VertexBuffer) {
					// Regular stream: get raw Xbox pointer
					pSrc = (const uint8_t*)GetDataFromXboxResource(streamInput.VertexBuffer)
						+ streamInput.Offset
						+ vertexStart * stride;
				}

				if (pSrc) {
					memcpy(pDst + streamOffsets[s], pSrc, streamDataSize);
				}
			}

			g_pD3DDeviceContext->Unmap(s_pVtxDataBuf, 0);
		}
	} // end if (!bUsingMirror)

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
		if (bUsingMirror && idxAddr >= CONTIGUOUS_MEMORY_BASE
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
				return true;

			HRESULT hr = CxbxD3D11UpdateDynamicBuffer(s_pIdxDataBuf, DrawContext.pXboxIndexData, DrawContext.dwVertexCount * sizeof(INDEX16));
			if (FAILED(hr)) return true;

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
		UINT drawLocalKey = primType | (indexedDraw << 2) | (vertexStart << 8) | (numVertices << 20);
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
		if (FAILED(hr)) return true;

		IABypassLayoutCB* pCB = (IABypassLayoutCB*)mapped.pData;
		memset(pCB, 0, sizeof(IABypassLayoutCB));

		pCB->PrimType = primType;
		pCB->IndexedDraw = indexedDraw;
		pCB->IndexOffset = indexOffset;
		pCB->NumAttribs = 16; // Always provide all 16 attribute descriptors
		pCB->NumVerts = DrawContext.dwVertexCount; // Original vertex count (for lineloop)

		// Fill per-attribute descriptors
		// Walk the vertex declaration's stream info to find each attribute
		for (UINT a = 0; a < 16; a++) {
			pCB->Attribs[a][0] = 0;  // elemOffset
			pCB->Attribs[a][1] = 0;  // stride
			pCB->Attribs[a][2] = CXBX_VTXFMT_NONE; // format — default is "use default value"
			pCB->Attribs[a][3] = 0;  // streamBase
		}

		// Map from the vertex declaration into attribute descriptors
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
			for (UINT e = 0; e < streamInfo.NumberOfVertexElements; e++) {
				auto& elem = streamInfo.VertexElements[e];
				if (elem.XboxType == 0) // X_D3DVSDT_NONE
					continue;

				// Determine which attribute register this element maps to.
				// The register index is encoded in the D3D11 input elements.
				// For now, use the element's position as a sequential register
				// starting from a base that depends on the stream index.
				// Actually — we need the register index from the vertex
				// declaration. Let's use the D3D11InputElements which have SemanticIndex = register.
				// But those might not be exactly parallel... let's find the right approach.

				// The pD3D11InputElements array has SemanticIndex = NV2A register index.
				// We need to find which register this stream element maps to.
				// The simplest way: search pD3D11InputElements for an element with matching
				// InputSlot == streamIdx and AlignedByteOffset == elemOffset.
				UINT regIdx = 0;
				bool found = false;
				if (pDecl->pD3D11InputElements) {
					for (UINT ie = 0; ie < pDecl->D3D11InputElementCount; ie++) {
						auto& inputElem = pDecl->pD3D11InputElements[ie];
						if (inputElem.InputSlot == streamIdx
							&& inputElem.AlignedByteOffset == elemOffset) {
							regIdx = inputElem.SemanticIndex;
							found = true;
							break;
						}
					}
				}

				if (found && regIdx < 16) {
					INT streamBase;
					if (bUsingMirror && streamInput.VertexBuffer) {
						// Mirror path: streamBase is the raw byte offset from CONTIGUOUS_MEMORY_BASE
						// to the start of the VB data. The shader does:
						//   byteOff = streamBase + vtxIdx * stride + elemOffset
						// GetDataFromXboxResource returns a pointer in the 0x80000000 region;
						// subtract CONTIGUOUS_MEMORY_BASE to get the 27-bit offset.
						uintptr_t vbAddr = (uintptr_t)GetDataFromXboxResource(streamInput.VertexBuffer);
						streamBase = (INT)(vbAddr - CONTIGUOUS_MEMORY_BASE) + (INT)streamInput.Offset;
					} else {
						// Fallback path: packed stream data uploaded to s_pVtxDataBuf.
						// streamBase = streamOffsets[s] - vertexStart * stride.
						streamBase = (INT)streamOffsets[s] - (INT)vertexStart * (INT)stride;
					}

					pCB->Attribs[regIdx][0] = elemOffset;
					pCB->Attribs[regIdx][1] = stride;
					pCB->Attribs[regIdx][2] = XboxFormatToVtxFmt(elem.XboxType);
					pCB->Attribs[regIdx][3] = (UINT)streamBase;
				}

				elemOffset += elem.XboxByteSize;
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
	// When using the mirror, all SRVs come from the page tracker's buffer views.
	// Otherwise, use the per-draw fallback buffer's views.
	ID3D11ShaderResourceView* pActiveVtxSRV = bUsingMirror ? pMirrorSRV : s_pVtxDataSRV;
	ID3D11ShaderResourceView* pActiveIdxSRV = bIdxFromMirror ? pMirrorSRV : s_pIdxDataSRV;
	ID3D11ShaderResourceView* pActiveSNormSRV = bUsingMirror
		? CxbxPageTrackerGetMirrorSRV_SNORM16x2() : s_pVtxDataSRV_SNORM16x2;
	ID3D11ShaderResourceView* pActiveUNormSRV = bUsingMirror
		? CxbxPageTrackerGetMirrorSRV_UNORM8x4() : s_pVtxDataSRV_UNORM8x4;

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

	return true;
}

#endif // CXBX_USE_D3D11
