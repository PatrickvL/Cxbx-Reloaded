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

#include "Backend_D3D11_Internal.h"
#include "common/util/hasher.h"
#include <unordered_map>

// ******************************************************************
// * Elastic UAV cache for compute shader dispatch targets
// * Grows on demand up to a hardware-realistic ceiling, then evicts
// * least-recently-used entries.  Sized for the heaviest Xbox titles
// * (e.g. palette-animated textures cycling hundreds of surfaces).
// ******************************************************************

// Xbox has 64 MB unified RAM.  Even the most texture-heavy titles
// (Jet Set Radio Future, Dead or Alive 3) rarely exceed ~200 active
// textures simultaneously.  We cap UAV cache growth at 128 entries
// (each UAV is ~16 bytes of driver state) and evict down to 96.
static constexpr size_t UAV_CACHE_HIGH_WATERMARK = 128;
static constexpr size_t UAV_CACHE_LOW_WATERMARK  = 96;

struct UAVCacheKey {
	ID3D11Resource* pResource;
	DXGI_FORMAT format;
	bool operator==(const UAVCacheKey& other) const {
		return pResource == other.pResource && format == other.format;
	}
};

struct UAVCacheKeyHash {
	size_t operator()(const UAVCacheKey& k) const {
		return static_cast<size_t>(ComputeHash(&k, sizeof(k)));
	}
};

struct UAVCacheValue {
	ID3D11UnorderedAccessView* pUAV;
	uint32_t lastUsed; // monotonic access counter for LRU eviction
};

static std::unordered_map<UAVCacheKey, UAVCacheValue, UAVCacheKeyHash> s_UAVCache;
static uint32_t s_UAVCacheAccessCounter = 0;

static void CxbxPruneUAVCache()
{
	if (s_UAVCache.size() < UAV_CACHE_HIGH_WATERMARK)
		return;

	// Find the eviction threshold (nth_element on lastUsed values).
	// Stack array avoids heap alloc — size is bounded by UAV_CACHE_HIGH_WATERMARK.
	size_t evictCount = s_UAVCache.size() - UAV_CACHE_LOW_WATERMARK;
	uint32_t stamps[UAV_CACHE_HIGH_WATERMARK];
	size_t n = 0;
	for (auto& kv : s_UAVCache)
		stamps[n++] = kv.second.lastUsed;
	std::nth_element(stamps, stamps + evictCount, stamps + n);
	uint32_t threshold = stamps[evictCount];

	for (auto it = s_UAVCache.begin(); it != s_UAVCache.end(); ) {
		if (it->second.lastUsed <= threshold) {
			if (it->second.pUAV) it->second.pUAV->Release();
			it = s_UAVCache.erase(it);
		} else {
			++it;
		}
	}
}

// Cache lookup uses raw resource pointer as key.  This is safe because the
// UAV internally AddRefs the resource, preventing destruction while cached.
// Eviction releases the UAV (and its implicit resource ref) via Release().
static ID3D11UnorderedAccessView* CxbxGetOrCreateTextureUAV(ID3D11Texture2D* pTexture, DXGI_FORMAT format)
{
	UAVCacheKey key = { pTexture, format };
	auto it = s_UAVCache.find(key);
	if (it != s_UAVCache.end()) {
		it->second.lastUsed = ++s_UAVCacheAccessCounter;
		return it->second.pUAV;
	}

	// Create new UAV
	D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
	uavDesc.Format = format;
	uavDesc.ViewDimension = D3D11_UAV_DIMENSION_TEXTURE2D;
	uavDesc.Texture2D.MipSlice = 0;
	ID3D11UnorderedAccessView* pUAV = nullptr;
	HRESULT hr = g_pD3DDevice->CreateUnorderedAccessView(pTexture, &uavDesc, &pUAV);
	if (FAILED(hr))
		return nullptr;

	s_UAVCache[key] = { pUAV, ++s_UAVCacheAccessCounter };

	// Prune if we've grown past the high watermark
	CxbxPruneUAVCache();

	return pUAV;
}

// ******************************************************************
// * GPU unswizzle via compute shader
// ******************************************************************
static void CxbxEnsureUnswizzleStagingBuffer(UINT requiredSize)
{
	CxbxD3D11EnsureRawStagingBuffer(requiredSize,
		&g_pD3D11UnswizzleStagingBuf, &g_UnswizzleStagingBufSize,
		&g_pD3D11UnswizzleSRV, "CxbxEnsureUnswizzleStagingBuffer");
}

// Typed unswizzle format decode constants (must match CxbxUnswizzleBGRA_CS.hlsl)
#define CXBX_UNSW_DECODE_BGRA8      0
#define CXBX_UNSW_DECODE_B4G4R4A4   1
#define CXBX_UNSW_DECODE_B5G6R5     2
#define CXBX_UNSW_DECODE_B5G5R5A1   3
#define CXBX_UNSW_DECODE_R10G10B10A2 4

// Returns the typed-unswizzle decode constant for a format, or -1 if not supported.
static int CxbxGetTypedUnswizzleDecode(DXGI_FORMAT format)
{
	switch (format) {
	case DXGI_FORMAT_B8G8R8A8_UNORM:
	case DXGI_FORMAT_B8G8R8X8_UNORM:     return CXBX_UNSW_DECODE_BGRA8;
	case DXGI_FORMAT_B4G4R4A4_UNORM:     return CXBX_UNSW_DECODE_B4G4R4A4;
	case DXGI_FORMAT_B5G6R5_UNORM:       return CXBX_UNSW_DECODE_B5G6R5;
	case DXGI_FORMAT_B5G5R5A1_UNORM:     return CXBX_UNSW_DECODE_B5G5R5A1;
	case DXGI_FORMAT_R10G10B10A2_UNORM:  return CXBX_UNSW_DECODE_R10G10B10A2;
	default: return -1;
	}
}

bool CxbxD3D11UnswizzleTexture(
	ID3D11Texture2D* pTexture,
	const void* pSwizzledSrc,
	UINT width,
	UINT height,
	UINT bpp,
	DXGI_FORMAT format)
{
	if (!g_pD3D11UnswizzleCS || !g_pD3D11UnswizzleCB)
		return false;

	// Only 2D textures with bpp 1, 2, or 4
	if (bpp != 1 && bpp != 2 && bpp != 4)
		return false;

	UINT dataSize = width * height * bpp;
	// Round up to DWORD alignment for ByteAddressBuffer
	UINT bufferSize = (dataSize + 3) & ~3u;

	// Ensure staging buffer is large enough
	CxbxEnsureUnswizzleStagingBuffer(bufferSize);
	if (!g_pD3D11UnswizzleStagingBuf || !g_pD3D11UnswizzleSRV)
		return false;

	// Upload swizzled source data to the staging buffer
	HRESULT hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11UnswizzleStagingBuf, pSwizzledSrc, dataSize);
	if (FAILED(hr))
		return false;

	// Compute Morton masks (same algorithm as EmuUnswizzleBox)
	UINT maskX = 0, maskY = 0;
	for (UINT i = 1, j = 1; (i < width) || (i < height); i <<= 1) {
		if (i < width)  { maskX |= j; j <<= 1; }
		if (i < height) { maskY |= j; j <<= 1; }
	}

	// Invalidate any cached PS SRV for this texture before binding it as a UAV.
	// Prevents SRV/UAV resource hazards that can trigger GPU TDRs.
	CxbxD3D11InvalidateCachedSRVForTexture(pTexture);

	// Try the fast path: R_UINT UAV reinterpretation (same typeless family)
	DXGI_FORMAT uintFormat;
	switch (bpp) {
	case 4:  uintFormat = DXGI_FORMAT_R32_UINT; break;
	case 2:  uintFormat = DXGI_FORMAT_R16_UINT; break;
	case 1:  uintFormat = DXGI_FORMAT_R8_UINT;  break;
	default: return false;
	}

	ID3D11UnorderedAccessView* pUAV = CxbxGetOrCreateTextureUAV(pTexture, uintFormat);
	if (pUAV) {
		// Fast path: raw uint bit-move via the standard unswizzle CS
		UINT cbData[8] = { maskX, maskY, width, height, bpp, 0, 0, 0 };
		hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11UnswizzleCB, cbData, sizeof(cbData));
		if (FAILED(hr))
			return false;

		UINT groupsX = (width + 7) / 8;
		UINT groupsY = (height + 7) / 8;
		CxbxD3D11DispatchCS(g_pD3D11UnswizzleCS, g_pD3D11UnswizzleCB,
			1, &g_pD3D11UnswizzleSRV, pUAV, groupsX, groupsY, 1);
		return true;
	}

	// R_UINT failed (cross-family format) — use typed float4 CS with same-format UAV
	if (!g_pD3D11UnswizzleBGRA_CS)
		return false;

	int fmtDecode = CxbxGetTypedUnswizzleDecode(format);
	if (fmtDecode < 0) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11UnswizzleTexture: No typed decode for format %u", format);
		return false;
	}

	pUAV = CxbxGetOrCreateTextureUAV(pTexture, format);
	if (!pUAV) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11UnswizzleTexture: Failed to create same-format UAV (format=%u)", format);
		return false;
	}

	UINT cbData[8] = { maskX, maskY, width, height, bpp, (UINT)fmtDecode, 0, 0 };
	hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11UnswizzleCB, cbData, sizeof(cbData));
	if (FAILED(hr))
		return false;

	UINT groupsX = (width + 7) / 8;
	UINT groupsY = (height + 7) / 8;
	CxbxD3D11DispatchCS(g_pD3D11UnswizzleBGRA_CS, g_pD3D11UnswizzleCB,
		1, &g_pD3D11UnswizzleSRV, pUAV, groupsX, groupsY, 1);
	return true;
}

// ******************************************************************
// * GPU palette texture expansion
// ******************************************************************

bool CxbxD3D11ExpandPaletteTexture(
	ID3D11Texture2D* pTexture,
	const void* pSwizzledP8Src,
	UINT width,
	UINT height,
	const void* pPaletteData)
{
	if (!g_pD3D11PaletteExpandCS || !g_pD3D11PaletteExpandCB || !g_pD3D11PaletteBuf || !g_pD3D11PaletteSRV)
		return false;

	UINT dataSize = width * height; // 1 byte per pixel
	UINT bufferSize = (dataSize + 3) & ~3u;

	// Reuse the unswizzle staging buffer for P8 source data
	CxbxEnsureUnswizzleStagingBuffer(bufferSize);
	if (!g_pD3D11UnswizzleStagingBuf || !g_pD3D11UnswizzleSRV)
		return false;

	// Upload swizzled P8 source data
	HRESULT hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11UnswizzleStagingBuf, pSwizzledP8Src, dataSize);
	if (FAILED(hr))
		return false;

	// Upload palette data (256 entries * 4 bytes = 1024 bytes)
	hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11PaletteBuf, pPaletteData, 1024);
	if (FAILED(hr))
		return false;

	// Compute Morton masks
	UINT maskX = 0, maskY = 0;
	for (UINT i = 1, j = 1; (i < width) || (i < height); i <<= 1) {
		if (i < width)  { maskX |= j; j <<= 1; }
		if (i < height) { maskY |= j; j <<= 1; }
	}

	// Update constant buffer
	UINT cbData[4] = { maskX, maskY, width, height };
	hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11PaletteExpandCB, cbData, sizeof(cbData));
	if (FAILED(hr))
		return false;

	// Invalidate any cached PS SRV for this texture before binding it as a UAV.
	// Prevents SRV/UAV resource hazards that can trigger GPU TDRs.
	CxbxD3D11InvalidateCachedSRVForTexture(pTexture);

	// Get or create a cached UAV for destination texture (R8G8B8A8_UNORM → R32_UINT UAV)
	ID3D11UnorderedAccessView* pUAV = CxbxGetOrCreateTextureUAV(pTexture, DXGI_FORMAT_R32_UINT);
	if (!pUAV) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11ExpandPaletteTexture: Failed to create UAV");
		return false;
	}

	// Dispatch 8x8 thread groups: t0=P8 source, t1=palette, u0=destination
	ID3D11ShaderResourceView* srvs[2] = { g_pD3D11UnswizzleSRV, g_pD3D11PaletteSRV };
	UINT groupsX = (width + 7) / 8;
	UINT groupsY = (height + 7) / 8;
	CxbxD3D11DispatchCS(g_pD3D11PaletteExpandCS, g_pD3D11PaletteExpandCB,
		2, srvs, pUAV, groupsX, groupsY, 1);
	return true;
}

// ******************************************************************
// * GPU texture format conversion (unswizzle + decode → RGBA)
// ******************************************************************

bool CxbxD3D11FormatConvertTexture(
	ID3D11Texture2D* pTexture,
	const void* pSrc,
	UINT width,
	UINT height,
	UINT bpp,
	UINT fmtType,
	bool bSwizzled,
	UINT srcRowPitch)
{
	if (!g_pD3D11FormatConvertCS || !g_pD3D11FormatConvertCB)
		return false;

	if (bpp != 1 && bpp != 2 && bpp != 4)
		return false;

	UINT dataSize = bSwizzled ? (width * height * bpp) : (srcRowPitch * height);
	UINT bufferSize = (dataSize + 3) & ~3u;

	// Reuse the unswizzle staging buffer
	CxbxEnsureUnswizzleStagingBuffer(bufferSize);
	if (!g_pD3D11UnswizzleStagingBuf || !g_pD3D11UnswizzleSRV)
		return false;

	HRESULT hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11UnswizzleStagingBuf, pSrc, dataSize);
	if (FAILED(hr))
		return false;

	// Compute Morton masks (only used if swizzled, but always computed)
	UINT maskX = 0, maskY = 0;
	for (UINT i = 1, j = 1; (i < width) || (i < height); i <<= 1) {
		if (i < width)  { maskX |= j; j <<= 1; }
		if (i < height) { maskY |= j; j <<= 1; }
	}

	// Update constant buffer: 8 uints
	UINT cbData[8] = { maskX, maskY, width, height, bpp, fmtType, bSwizzled ? 1u : 0u, srcRowPitch };
	hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11FormatConvertCB, cbData, sizeof(cbData));
	if (FAILED(hr))
		return false;

	// Invalidate any cached PS SRV for this texture before binding it as a UAV.
	// Prevents SRV/UAV resource hazards that can trigger GPU TDRs.
	CxbxD3D11InvalidateCachedSRVForTexture(pTexture);

	// Get or create a cached UAV (R32_UINT view of the R8G8B8A8_UNORM texture)
	ID3D11UnorderedAccessView* pUAV = CxbxGetOrCreateTextureUAV(pTexture, DXGI_FORMAT_R32_UINT);
	if (!pUAV) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11FormatConvertTexture: Failed to create UAV");
		return false;
	}

	UINT groupsX = (width + 7) / 8;
	UINT groupsY = (height + 7) / 8;
	CxbxD3D11DispatchCS(g_pD3D11FormatConvertCS, g_pD3D11FormatConvertCB,
		1, &g_pD3D11UnswizzleSRV, pUAV, groupsX, groupsY, 1);

	return true;
}

// ******************************************************************
// * GPU vertex format conversion
// ******************************************************************

static void CxbxEnsureVertexConvertSrcBuffer(UINT requiredSize)
{
	CxbxD3D11EnsureRawStagingBuffer(requiredSize,
		&g_pD3D11VertexConvertSrcBuf, &g_VertexConvertSrcBufSize,
		&g_pD3D11VertexConvertSrcSRV, "CxbxEnsureVertexConvertSrcBuffer");
}

bool CxbxD3D11ConvertVertexBufferGPU(
	const uint8_t* pSrcVertexData,
	UINT srcDataSize,
	UINT vertexCount,
	UINT srcStride,
	UINT dstStride,
	UINT numElements,
	const UINT* pElementDescriptors, // 4 UINTs per element: srcOffset, dstOffset, convType, copyDwords
	UINT dstBufferSize,
	ID3D11Buffer** ppOutputVB)
{
	if (!g_pD3D11VertexConvertCS || !g_pD3D11VertexConvertCB || vertexCount == 0)
		return false;

	// Require dst stride to be a multiple of 4 (for RWBuffer<uint> writes)
	if ((dstStride & 3) != 0)
		return false;

	UINT bufferSize = (srcDataSize + 3) & ~3u;

	// Upload source vertex data
	CxbxEnsureVertexConvertSrcBuffer(bufferSize);
	if (!g_pD3D11VertexConvertSrcBuf || !g_pD3D11VertexConvertSrcSRV)
		return false;

	HRESULT hr = CxbxD3D11UpdateDynamicBuffer(g_pD3D11VertexConvertSrcBuf, pSrcVertexData, srcDataSize);
	if (FAILED(hr))
		return false;

	// Update constant buffer: header (4 uints) + elements (numElements * 4 uints)
	D3D11_MAPPED_SUBRESOURCE mapped = {};
	hr = g_pD3DDeviceContext->Map(g_pD3D11VertexConvertCB, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
	if (FAILED(hr))
		return false;
	UINT* pCB = (UINT*)mapped.pData;
	pCB[0] = vertexCount;
	pCB[1] = srcStride;
	pCB[2] = dstStride;
	pCB[3] = numElements;
	memcpy(&pCB[4], pElementDescriptors, numElements * 4 * sizeof(UINT));
	// Zero unused element slots
	if (numElements < 16)
		memset(&pCB[4 + numElements * 4], 0, (16 - numElements) * 4 * sizeof(UINT));
	g_pD3DDeviceContext->Unmap(g_pD3D11VertexConvertCB, 0);

	// Create output vertex buffer (DEFAULT + VERTEX_BUFFER + UAV)
	D3D11_BUFFER_DESC outDesc = {};
	outDesc.ByteWidth = dstBufferSize;
	outDesc.Usage = D3D11_USAGE_DEFAULT;
	outDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER | D3D11_BIND_UNORDERED_ACCESS;
	ID3D11Buffer* pOutputBuf = nullptr;
	hr = g_pD3DDevice->CreateBuffer(&outDesc, nullptr, &pOutputBuf);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11ConvertVertexBufferGPU: Failed to create output VB (hr=0x%08X)", hr);
		return false;
	}

	ID3D11UnorderedAccessView* pUAV = nullptr;
	hr = CxbxD3D11CreateTypedBufferUAV(pOutputBuf, dstBufferSize, &pUAV);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11ConvertVertexBufferGPU: Failed to create UAV");
		pOutputBuf->Release();
		return false;
	}

	// Dispatch CS
	UINT numGroups = (vertexCount + 63) / 64;
	CxbxD3D11DispatchCS(g_pD3D11VertexConvertCS, g_pD3D11VertexConvertCB,
		1, &g_pD3D11VertexConvertSrcSRV, pUAV, numGroups, 1, 1);

	pUAV->Release();
	*ppOutputVB = pOutputBuf;
	return true;
}

// ******************************************************************
// * GPU index buffer conversion helpers
// ******************************************************************

static bool CxbxEnsureIndexConvertInputBuffer(UINT requiredSize)
{
	if (g_pD3D11IndexConvertInputBuf && g_IndexConvertInputBufSize >= requiredSize)
		return true;

	// Release old resources
	if (g_pD3D11IndexConvertInputSRV) { g_pD3D11IndexConvertInputSRV->Release(); g_pD3D11IndexConvertInputSRV = nullptr; }
	if (g_pD3D11IndexConvertInputBuf) { g_pD3D11IndexConvertInputBuf->Release(); g_pD3D11IndexConvertInputBuf = nullptr; }

	UINT newSize = (requiredSize + 4095) & ~4095u;
	if (newSize < 4) newSize = 4;

	HRESULT hr = CxbxD3D11CreateRawBuffer(newSize, /*bDynamic=*/false, 0, &g_pD3D11IndexConvertInputBuf);
	if (FAILED(hr)) return false;

	hr = CxbxD3D11CreateRawBufferSRV(g_pD3D11IndexConvertInputBuf, newSize, &g_pD3D11IndexConvertInputSRV);
	if (FAILED(hr)) {
		g_pD3D11IndexConvertInputBuf->Release();
		g_pD3D11IndexConvertInputBuf = nullptr;
		return false;
	}

	g_IndexConvertInputBufSize = newSize;
	return true;
}

static bool CxbxEnsureIndexConvertOutputBuffer(UINT requiredByteSize)
{
	if (g_pD3D11IndexConvertOutputBuf && g_IndexConvertOutputBufSize >= requiredByteSize)
		return true;

	// Release old resources
	if (g_pD3D11IndexConvertOutputUAV) { g_pD3D11IndexConvertOutputUAV->Release(); g_pD3D11IndexConvertOutputUAV = nullptr; }
	if (g_pD3D11IndexConvertOutputBuf) { g_pD3D11IndexConvertOutputBuf->Release(); g_pD3D11IndexConvertOutputBuf = nullptr; }

	UINT newSize = (requiredByteSize + 4095) & ~4095u;
	if (newSize < 4) newSize = 4;

	D3D11_BUFFER_DESC desc = {};
	desc.ByteWidth = newSize;
	desc.Usage = D3D11_USAGE_DEFAULT;
	desc.BindFlags = D3D11_BIND_INDEX_BUFFER | D3D11_BIND_UNORDERED_ACCESS;

	HRESULT hr = g_pD3DDevice->CreateBuffer(&desc, nullptr, &g_pD3D11IndexConvertOutputBuf);
	if (FAILED(hr)) return false;

	hr = CxbxD3D11CreateTypedBufferUAV(g_pD3D11IndexConvertOutputBuf, newSize, &g_pD3D11IndexConvertOutputUAV);
	if (FAILED(hr)) {
		g_pD3D11IndexConvertOutputBuf->Release();
		g_pD3D11IndexConvertOutputBuf = nullptr;
		return false;
	}

	g_IndexConvertOutputBufSize = newSize;
	return true;
}

bool CxbxD3D11ConvertIndexBufferGPU(
	const INDEX16* pSourceIndices,
	UINT sourceVertexCount,
	UINT outputIndexCount,
	int conversionMode)
{
	if (!g_pD3D11IndexConvertCS || !g_pD3D11IndexConvertCB || outputIndexCount == 0)
		return false;

	// 1. Upload source indices to input staging buffer (only if indexed)
	if (pSourceIndices != nullptr) {
		UINT inputByteSize = sourceVertexCount * sizeof(INDEX16);
		if (!CxbxEnsureIndexConvertInputBuffer(inputByteSize))
			return false;

		D3D11_BOX box = { 0, 0, 0, inputByteSize, 1, 1 };
		g_pD3DDeviceContext->UpdateSubresource(g_pD3D11IndexConvertInputBuf, 0, &box, pSourceIndices, inputByteSize, 0);
	}

	// 2. Ensure output buffer is large enough
	// Output: ceil(outputIndexCount / 2) uint32 elements (packed 16-bit pairs)
	UINT outputUintCount = (outputIndexCount + 1) / 2;
	UINT outputByteSize = outputUintCount * sizeof(UINT);
	if (!CxbxEnsureIndexConvertOutputBuffer(outputByteSize))
		return false;

	// 3. Update constant buffer
	struct { UINT VertexCount, ConversionMode, IsIndexed, Pad; } cbData = {
		sourceVertexCount,
		(UINT)conversionMode,
		pSourceIndices != nullptr ? 1u : 0u,
		0
	};
	g_pD3DDeviceContext->UpdateSubresource(g_pD3D11IndexConvertCB, 0, nullptr, &cbData, sizeof(cbData), 0);

	// 4. Bind resources and dispatch
	UINT numSRVs = (pSourceIndices != nullptr) ? 1 : 0;
	UINT numWorkItems;
	if (conversionMode <= CXBX_INDEX_CONVERT_QUAD_CCW) {
		numWorkItems = sourceVertexCount / 4; // one thread per quad
	} else {
		UINT numTris = (sourceVertexCount >= 3) ? (sourceVertexCount - 2) : 0;
		numWorkItems = (numTris + 1) / 2; // one thread per pair of triangles
	}
	UINT numGroups = (numWorkItems + 63) / 64;
	CxbxD3D11DispatchCS(g_pD3D11IndexConvertCS, g_pD3D11IndexConvertCB,
		numSRVs, &g_pD3D11IndexConvertInputSRV, g_pD3D11IndexConvertOutputUAV,
		numGroups, 1, 1);

	// 5. Bind output as index buffer (R16_UINT — the CS packed 16-bit indices into 32-bit writes)
	g_pD3DDeviceContext->IASetIndexBuffer(g_pD3D11IndexConvertOutputBuf, DXGI_FORMAT_R16_UINT, 0);

	return true;
}

