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

// Backend_D3D11_VertexFetch.cpp — Programmable vertex fetching draw path.
//
// This module handles vertex data upload
// and draw calls using SV_VertexID-based vertex fetch in the shader.
// The Input Assembler is not used for vertex/index buffer binding.

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h"
#include "Backend_D3D11_Profiler.h"
#include "common/AddressRanges.h"
#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\XbPushBuffer.h" // NV2A_get_vertex_attribute_value_pointer
#include "core\hle\D3D8\Rendering\IndexBufferConvert.h" // CxbxGetClockWiseWindingOrder
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
// * Persistent GPU resources for vertex fetch
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
static D3D_PRIMITIVE_TOPOLOGY    s_LastTopology = D3D_PRIMITIVE_TOPOLOGY_UNDEFINED; // Topology caching

void CxbxInvalidateTopologyCache()
{
	s_LastTopology = D3D_PRIMITIVE_TOPOLOGY_UNDEFINED;
}

// ******************************************************************
// * Map Xbox primitive type to host topology and compute host vertex count.
// * Returns false if the primitive type is unsupported (caller should skip draw).
// ******************************************************************
static bool ResolveTopology(uint32_t primitiveMode, UINT vertexCount,
	UINT& outPrimType, D3D_PRIMITIVE_TOPOLOGY& outTopology, UINT& outHostVertexCount)
{
	outPrimType = CXBX_PRIM_NORMAL;
	outHostVertexCount = vertexCount;
	outTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;

	switch (primitiveMode) {
	case NV097_SET_BEGIN_END_OP_QUADS:
		outPrimType = CXBX_PRIM_QUAD;
		outHostVertexCount = (vertexCount / 4) * 6;
		break;
	case NV097_SET_BEGIN_END_OP_QUAD_STRIP:
		outPrimType = CXBX_PRIM_QUADSTRIP;
		outHostVertexCount = (vertexCount >= 4) ? ((vertexCount - 2) / 2) * 6 : 0;
		break;
	case NV097_SET_BEGIN_END_OP_TRIANGLE_FAN:
	case NV097_SET_BEGIN_END_OP_POLYGON:
		outPrimType = CXBX_PRIM_FAN;
		outHostVertexCount = (vertexCount >= 3) ? (vertexCount - 2) * 3 : 0;
		break;
	case NV097_SET_BEGIN_END_OP_TRIANGLES:
		break;
	case NV097_SET_BEGIN_END_OP_TRIANGLE_STRIP:
		outTopology = D3D_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP;
		break;
	case NV097_SET_BEGIN_END_OP_LINES:
		outTopology = D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		break;
	case NV097_SET_BEGIN_END_OP_LINE_STRIP:
		outTopology = D3D_PRIMITIVE_TOPOLOGY_LINESTRIP;
		break;
	case NV097_SET_BEGIN_END_OP_POINTS:
		outTopology = D3D_PRIMITIVE_TOPOLOGY_POINTLIST;
		break;
	case NV097_SET_BEGIN_END_OP_LINE_LOOP:
		outPrimType = CXBX_PRIM_LINELOOP;
		outHostVertexCount = vertexCount * 2;
		outTopology = D3D_PRIMITIVE_TOPOLOGY_LINELIST;
		break;
	default:
		return false; // Unsupported
	}
	return outHostVertexCount > 0;
}

// ******************************************************************
// * Shared draw tail: unbind IA, set topology, bind SRVs/CBs, issue draw.
// ******************************************************************
static void BindAndIssueDraw(UINT primType, D3D_PRIMITIVE_TOPOLOGY hostTopology,
	UINT hostVertexCount, uint32_t primitiveMode,
	ID3D11ShaderResourceView* pVtxSRV, ID3D11ShaderResourceView* pIdxSRV,
	ID3D11ShaderResourceView* pSNormSRV, ID3D11ShaderResourceView* pUNormSRV)
{
	// Unbind IA state (skip if already nulled from a prior vertex fetch draw)
	if (!s_IAAlreadyNull) {
		g_pD3DDeviceContext->IASetInputLayout(nullptr);
		ID3D11Buffer* nullBufs[17] = {};
		UINT nullStrides[17] = {};
		UINT nullOffsets[17] = {};
		g_pD3DDeviceContext->IASetVertexBuffers(0, 17, nullBufs, nullStrides, nullOffsets);
		g_pD3DDeviceContext->IASetIndexBuffer(nullptr, DXGI_FORMAT_R16_UINT, 0);
		s_IAAlreadyNull = true;
	}
	if (hostTopology != s_LastTopology) {
		g_pD3DDeviceContext->IASetPrimitiveTopology(hostTopology);
		s_LastTopology = hostTopology;
	}

	// Bind SRVs to VS: t0=vertex data, t1=index data, t2=snorm, t3=unorm
	if (pVtxSRV != s_pLastBoundVtxSRV || pIdxSRV != s_pLastBoundIdxSRV
		|| pSNormSRV != s_pLastBoundSNormSRV || pUNormSRV != s_pLastBoundUNormSRV) {
		ID3D11ShaderResourceView* vsSRVs[4] = { pVtxSRV, pIdxSRV, pSNormSRV, pUNormSRV };
		g_pD3DDeviceContext->VSSetShaderResources(0, 4, vsSRVs);
		s_pLastBoundVtxSRV = pVtxSRV;
		s_pLastBoundIdxSRV = pIdxSRV;
		s_pLastBoundSNormSRV = pSNormSRV;
		s_pLastBoundUNormSRV = pUNormSRV;
	}

	// Bind CBs: b1=layout, b2=defaults (skip if unchanged)
	if (s_pLayoutCB != s_pLastBoundLayoutCB || s_pDefaultsCB != s_pLastBoundDefaultsCB) {
		ID3D11Buffer* vsCBs[2] = { s_pLayoutCB, s_pDefaultsCB };
		g_pD3DDeviceContext->VSSetConstantBuffers(1, 2, vsCBs);
		s_pLastBoundLayoutCB = s_pLayoutCB;
		s_pLastBoundDefaultsCB = s_pDefaultsCB;
	}

	// Bind thick line GS if needed
	if (primType == CXBX_PRIM_NORMAL) {
		CxbxBindThickLineGS(primitiveMode);
	}

	{
		CXBX_PROFILE_SCOPE(PROF_DRAW_CALL);
		g_pD3DDeviceContext->Draw(hostVertexCount, 0);
	}
	g_ProfileDrawCount++;

	if (primType == CXBX_PRIM_NORMAL) {
		CxbxUnbindThickLineGS(primitiveMode);
	}
}

// Layout CB caching: generation counter bumped on state changes
static UINT                      s_LayoutCBGeneration = 0;
static UINT                      s_LastLayoutCBGeneration = UINT_MAX;

// ******************************************************************
// * Layout constant buffer structure (shared with HLSL via CxbxVertexFetchLayout.hlsli)
// ******************************************************************
#include "core\hle\D3D8\Rendering\Shaders\CxbxVertexFetchLayout.hlsli"

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
// * Initialize vertex fetch resources (called once during device init)
// ******************************************************************
void CxbxD3D11VertexFetchInit()
{
	HRESULT hr;

	// Layout CB (b1) — 32 + 256 = 288 bytes
	hr = CxbxD3D11CreateConstantBuffer(sizeof(VertexFetchLayoutCB), false, &s_pLayoutCB);
	if (FAILED(hr))
		EmuLog(LOG_LEVEL::WARNING, "VertexFetchInit: Failed to create layout CB");

	// Defaults CB (b2) — 16 × float4 = 256 bytes
	hr = CxbxD3D11CreateConstantBuffer(16 * 4 * sizeof(float), false, &s_pDefaultsCB);
	if (FAILED(hr))
		EmuLog(LOG_LEVEL::WARNING, "VertexFetchInit: Failed to create defaults CB");
}

// ******************************************************************
// * Release vertex fetch resources
// ******************************************************************
void CxbxD3D11VertexFetchRelease()
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
void CxbxD3D11VertexFetchInvalidateLayout()
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
		&s_pUPVtxDataSRV, "VertexFetch_UPVtxData");

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
		&s_pIdxDataSRV, "VertexFetch_IdxData");
}

// ******************************************************************
// * Upload vertex defaults (NV2A sticky attribute values) to CB b2
// ******************************************************************
// Dirty flag for vertex defaults — set by CxbxSetVertexAttribute, consumed here
bool g_bD3D11VertexFetchDefaultsDirty = true;

static void UploadVertexDefaults()
{
	if (!s_pDefaultsCB) return;
	if (!g_bD3D11VertexFetchDefaultsDirty) return;

	g_bD3D11VertexFetchDefaultsDirty = false;

	// Compare actual attribute values to skip Map/Unmap when nothing changed
	static float s_CachedDefaults[16 * 4] = {};
	static bool s_FirstCall = true;
	bool changed = s_FirstCall;
	for (int i = 0; i < 16 && !changed; i++) {
		const float* pSrc = NV2A_get_vertex_attribute_value_pointer(i);
		if (std::memcmp(pSrc, &s_CachedDefaults[i * 4], sizeof(float) * 4) != 0)
			changed = true;
	}
	if (!changed) return;
	s_FirstCall = false;

	float defaults[16 * 4];
	for (int i = 0; i < 16; i++) {
		// TODO: NV2A hardware updates inline_value[] with the last vertex's data
		// after each streamed draw, so a subsequent draw that doesn't stream an
		// attribute reads the value from the final vertex of the previous draw.
		// Currently, inline_value[] is only written by explicit pushbuffer commands
		// (NV097_SET_VERTEX_DATA4F etc.), not by the streamed vertex path.  To fix
		// this, after each draw we'd need to read back the last vertex's attribute
		// values from the CPU-side vertex data and write them to inline_value[].
		const float* pSrc = NV2A_get_vertex_attribute_value_pointer(i);
		defaults[i * 4 + 0] = pSrc[0];
		defaults[i * 4 + 1] = pSrc[1];
		defaults[i * 4 + 2] = pSrc[2];
		defaults[i * 4 + 3] = pSrc[3];
		std::memcpy(&s_CachedDefaults[i * 4], pSrc, sizeof(float) * 4);
	}

	g_pD3DDeviceContext->UpdateSubresource(s_pDefaultsCB, 0, nullptr, defaults, 0, 0);
}

// ******************************************************************
// * Core draw function for vertex fetch
// * Returns true if the draw was handled, false to fall back to IA path.
// ******************************************************************
void CxbxD3D11VertexFetchDraw(CxbxDrawContext& DrawContext)
{
	// When all vertex shaders are compiled with vertex fetch, the normal IA
	// fallback path cannot work (shader expects SV_VertexID, not TEXCOORD
	// inputs).
	if (!s_pLayoutCB || !s_pDefaultsCB)
		return;

	// ---------------------------------------------------------------
	// Step 1: Determine topology and host vertex count
	// ---------------------------------------------------------------
	UINT primType, hostVertexCount;
	D3D_PRIMITIVE_TOPOLOGY hostTopology;
	uint32_t primitiveMode = DrawContext.XboxPrimitiveType; // NV097_SET_BEGIN_END op value
	if (!ResolveTopology(primitiveMode, DrawContext.dwVertexCount,
		primType, hostTopology, hostVertexCount))
		return; // Unsupported topology

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
	{
		CXBX_PROFILE_SCOPE(PROF_PAGE_FLUSH);
		CxbxPageTrackerFlushToGPU();
	}

	// ---------------------------------------------------------------
	// Step 2a: Flush GPU-dirty RT pages that overlap VB streams
	// ---------------------------------------------------------------
	// When a game renders to a surface then reads it back as a vertex buffer
	// (e.g. DisplacementMap XDK sample using XGSetVertexBufferHeader), the
	// D3D11 RT content must be read back to Xbox RAM and re-uploaded to the
	// GPU mirror before the vertex shader can fetch correct data.
	if (!bIsUPDraw && g_NV2A != nullptr) {
		PGRAPHState* pgVB = &g_NV2A->GetDeviceState()->pgraph;
		for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
			const VertexAttribute& attr = pgVB->vertex_attributes[i];
			if (attr.count == 0) continue;
			// attr.offset is a physical address (27-bit, masked into 64 MiB)
			uint32_t vbOffset = (uint32_t)(attr.offset & 0x07FFFFFF);
			uint32_t vbSize = attr.stride * DrawContext.dwVertexCount;
			if (vbSize == 0) continue;
			CxbxPageTrackerFlushGPUDirtyToMirror(vbOffset, vbSize);
		}
	}

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
	// The generation is also bumped externally by CxbxD3D11VertexFetchInvalidateLayout()
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
		VertexFetchLayoutCB cb = {};

		cb.PrimType = primType;
		cb.IndexedDraw = indexedDraw;
		cb.IndexOffset = indexOffset;
		cb.NumAttribs = 16; // Always provide all 16 attribute descriptors
		cb.NumVerts = DrawContext.dwVertexCount; // Original vertex count (for lineloop)

		// VertexOffset: adjusts the resolved vertex index before VB fetch.
		// Non-indexed draws: StartVertex (DrawVertices skips the first N vertices).
		// Indexed draws: BaseVertexIndex (SetIndices offset added to each index).
		if (indexedDraw)
			cb.VertexOffset = DrawContext.dwBaseVertexIndex;
		else
			cb.VertexOffset = vertexStart;

		// Quad winding: must match NV2A SETUPRASTER front face setting
		cb.WindingCW = CxbxGetClockWiseWindingOrder() ? 1 : 0;

		// Fill per-attribute descriptors — default all to NONE (use sticky defaults)
		for (UINT a = 0; a < 16; a++) {
			cb.Attribs[a][0] = 0;  // elemOffset
			cb.Attribs[a][1] = 0;  // stride
			cb.Attribs[a][2] = CXBX_VTXFMT_NONE; // format — default is "use default value"
			cb.Attribs[a][3] = 0;  // streamBase
		}

		// PGRAPH path: reads vertex_attributes[] directly from NV2A state.
		// All draws come through push buffer → PFIFO → PGRAPH, so PGRAPH is
		// the authoritative source for vertex layout information.
		PGRAPHState* pg = (g_NV2A != nullptr) ? &g_NV2A->GetDeviceState()->pgraph : nullptr;

		if (pg) {
			if (DrawContext.bNV2AInlineData) {
				// NV2A inline_array path: data is packed contiguously per vertex
				// with enabled attributes in register order. Compute element offsets
				// from attribute sizes rather than using physical addresses.
				UINT packedOffset = 0;
				for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
					const VertexAttribute& attr = pg->vertex_attributes[i];
				if (attr.count == 0) continue; // count 0 = disabled (format 0 is valid: UB_D3D/D3DCOLOR)
					cb.Attribs[i][0] = packedOffset; // elemOffset within packed vertex
					cb.Attribs[i][1] = DrawContext.uiXboxVertexStreamZeroStride;
					cb.Attribs[i][2] = NV2AFormatToVtxFmt(attr.format, attr.count);
					cb.Attribs[i][3] = 0; // streamBase = 0 (UP staging buffer)
					packedOffset += attr.count * attr.size;
				}
			} else {
				// VB draw path: slot index = register index, offset = physical address
				for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
					const VertexAttribute& attr = pg->vertex_attributes[i];
					if (attr.count == 0) continue;

					cb.Attribs[i][0] = 0;           // elemOffset (baked into offset)
					cb.Attribs[i][1] = attr.stride;
					cb.Attribs[i][2] = NV2AFormatToVtxFmt(attr.format, attr.count);
					cb.Attribs[i][3] = (UINT)attr.offset; // physical addr = SRV byte offset
				}
			}
		} else {
			return; // No NV2A state — can't draw
		}

		g_pD3DDeviceContext->UpdateSubresource(s_pLayoutCB, 0, nullptr, &cb, 0, 0);
	}
skip_layout_upload:

	// ---------------------------------------------------------------
	// Step 5: Upload vertex defaults (skip if not dirty)
	// ---------------------------------------------------------------
	UploadVertexDefaults();

	// ---------------------------------------------------------------
	// Step 6: Bind resources and issue draw
	// ---------------------------------------------------------------
	// VB draws: SRVs from the 64 MiB page-tracked mirror.
	// UP draws: SRVs from the per-draw staging buffer.
	ID3D11ShaderResourceView* pActiveVtxSRV = bIsUPDraw ? s_pUPVtxDataSRV : pMirrorSRV;
	ID3D11ShaderResourceView* pActiveIdxSRV = bIdxFromMirror ? pMirrorSRV : s_pIdxDataSRV;
	ID3D11ShaderResourceView* pActiveSNormSRV = bIsUPDraw
		? s_pUPVtxDataSRV_SNORM16x2 : CxbxPageTrackerGetMirrorSRV_SNORM16x2();
	ID3D11ShaderResourceView* pActiveUNormSRV = bIsUPDraw
		? s_pUPVtxDataSRV_UNORM8x4 : CxbxPageTrackerGetMirrorSRV_UNORM8x4();

	BindAndIssueDraw(primType, hostTopology, hostVertexCount, primitiveMode,
		pActiveVtxSRV, pActiveIdxSRV, pActiveSNormSRV, pActiveUNormSRV);
}

// ******************************************************************
// * Draw inline buffer vertices (Begin/SetVertexData/End path)
// *
// * Inline buffer data is stored as float4 per attribute per vertex
// * in pg->vertex_attributes[i].inline_buffer. This bypasses the
// * normal vertex declaration/stream layout entirely — all 16
// * attributes are packed as FLOAT4 at a fixed 256-byte stride.
// * After drawing, inline_buffer pointers are reset (pool stays allocated).
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
	UINT primType, hostVertexCount;
	D3D_PRIMITIVE_TOPOLOGY hostTopology;
	uint32_t primitiveMode = pg->primitive_mode; // NV097_SET_BEGIN_END op value
	if (!ResolveTopology(primitiveMode, vertexCount, primType, hostTopology, hostVertexCount))
		return; // Unsupported topology

	// ---------------------------------------------------------------
	// Step 3: Fill layout CB with float4 layout for all 16 attributes
	// ---------------------------------------------------------------
	{
		VertexFetchLayoutCB cb = {};

		cb.PrimType = primType;
		cb.IndexedDraw = 0;
		cb.IndexOffset = 0;
		cb.NumAttribs = 16;
		cb.NumVerts = vertexCount;
		cb.VertexOffset = 0;

		// Quad winding: must match NV2A SETUPRASTER front face setting
		cb.WindingCW = CxbxGetClockWiseWindingOrder() ? 1 : 0;

		for (UINT a = 0; a < 16; a++) {
			cb.Attribs[a][0] = a * kAttrSize;       // elemOffset
			cb.Attribs[a][1] = kStride;             // stride
			cb.Attribs[a][2] = CXBX_VTXFMT_FLOAT4;  // format
			cb.Attribs[a][3] = 0;                   // streamBase (UP data at offset 0)
		}

		g_pD3DDeviceContext->UpdateSubresource(s_pLayoutCB, 0, nullptr, &cb, 0, 0);
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
	BindAndIssueDraw(primType, hostTopology, hostVertexCount, primitiveMode,
		s_pUPVtxDataSRV, nullptr, s_pUPVtxDataSRV_SNORM16x2, s_pUPVtxDataSRV_UNORM8x4);

	// Reset per-attribute inline buffers (pool stays allocated for reuse)
	for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
		VertexAttribute& attr = pg->vertex_attributes[i];
		if (attr.inline_buffer) {
			attr.inline_buffer = nullptr;
		}
	}
}

// ******************************************************************
// * Input Assembler vertex binding / input layout (legacy IA path)
// ******************************************************************

HRESULT CxbxSetStreamSource(UINT HostStreamNumber, ID3D11Buffer* pHostVertexBuffer, UINT VertexStride)
{
	UINT offset = 0;
	g_pD3DDeviceContext->IASetVertexBuffers(HostStreamNumber, 1, &pHostVertexBuffer, &VertexStride, &offset);
	return S_OK;
}

// ******************************************************************
// * Vertex defaults buffer — zero-stride buffer providing NV2A "sticky"
// * attribute values for non-streamed TEXCOORD inputs.
// ******************************************************************
ID3D11Buffer *g_pD3D11VertexDefaultsBuffer = nullptr;

void CxbxD3D11CreateVertexDefaultsBuffer()
{
	// 16 attributes × 4 floats × 4 bytes = 256 bytes
	D3D11_BUFFER_DESC desc = {};
	desc.ByteWidth = X_VSH_MAX_ATTRIBUTES * 4 * sizeof(float);
	desc.Usage = D3D11_USAGE_DYNAMIC;
	desc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

	// Initialize with NV2A default values (0,0,0,1) per attribute
	float initData[X_VSH_MAX_ATTRIBUTES * 4];
	for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
		initData[i * 4 + 0] = 0.0f;
		initData[i * 4 + 1] = 0.0f;
		initData[i * 4 + 2] = 0.0f;
		initData[i * 4 + 3] = 1.0f;
	}

	D3D11_SUBRESOURCE_DATA srd = {};
	srd.pSysMem = initData;

	HRESULT hr = g_pD3DDevice->CreateBuffer(&desc, &srd, &g_pD3D11VertexDefaultsBuffer);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11CreateVertexDefaultsBuffer: CreateBuffer failed (0x%08X)", hr);
		return;
	}

	// Bind to the defaults slot with stride=0 (every vertex reads the same data)
	UINT stride = 0;
	UINT offset = 0;
	g_pD3DDeviceContext->IASetVertexBuffers(CXBX_D3D11_VERTEX_DEFAULTS_SLOT, 1,
		&g_pD3D11VertexDefaultsBuffer, &stride, &offset);
}
