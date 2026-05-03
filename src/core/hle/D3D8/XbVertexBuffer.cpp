// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
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
// *  (c) 2002-2004 Aaron Robinson <caustik@caustik.com>
// *                Kingofc <kingofc@freenet.de>
// *
// *  All rights reserved
// *
// ******************************************************************
#define LOG_PREFIX CXBXR_MODULE::VTXB

#include <unordered_map>
#include "core\kernel\memory-manager\VMManager.h"
#include "common\util\hasher.h"
#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\Rendering\RenderGlobals.h" // For CxbxSetStreamSource, CxbxCreateVertexBuffer, etc.
#include "core\hle\D3D8\Rendering\WalkIndexBuffer.h" // for WalkIndexBuffer
#include "core\hle\D3D8\ResourceTracker.h"
#include "core\hle\D3D8\XbPushBuffer.h" // For CxbxDrawPrimitiveUP
#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11.h"

#include <imgui.h>

#include <ctime>
#include <chrono>
#include <algorithm>

CxbxVertexBufferConverter VertexBufferConverter = {};

// Copy of active Xbox D3D Vertex Streams (and strides), set by [D3DDevice|CxbxImpl]_SetStreamSource*
xbox::X_STREAMINPUT g_Xbox_SetStreamSource[X_VSH_MAX_STREAMS] = { 0 }; // Note : .Offset member is never set (so always 0)

// CxbxSetStreamSource is in Backend_D3D11_Draw.cpp

void CxbxPatchedStream::Activate(CxbxDrawContext *pDrawContext, UINT HostStreamNumber) const
{
	//LOG_INIT // Allows use of DEBUG_D3DRESULT

	// Use the cached stream values on the host
	if (bCacheIsStreamZeroDrawUP) {
		// Set the UserPointer variables in the drawing context
		pDrawContext->pHostVertexStreamZeroData = pCachedHostVertexStreamZeroData;
		pDrawContext->uiHostVertexStreamZeroStride = uiCachedHostVertexStride;
	}
	else {
		HRESULT hRet = CxbxSetStreamSource(
			HostStreamNumber,
			pCachedHostVertexBuffer, 
			uiCachedHostVertexStride);
		if (FAILED(hRet)) {
			CxbxrAbort("Failed to set the type patched buffer as the new stream source!\n");
			// TODO : test-case : XDK Cartoon hits the above case when the vertex cache size is 0.
		}
	}
}

void CxbxPatchedStream::Clear()
{
    if (bCachedHostVertexStreamZeroDataIsAllocated) {
        free(pCachedHostVertexStreamZeroData);
        bCachedHostVertexStreamZeroDataIsAllocated = false;
    }

    pCachedHostVertexStreamZeroData = nullptr;

    if (pCachedHostVertexBuffer != nullptr) {
        pCachedHostVertexBuffer->Release();
        pCachedHostVertexBuffer = nullptr;
    }
}

CxbxPatchedStream::~CxbxPatchedStream()
{
	Clear();
}

// TODO: CountActiveD3DStreams must be removed once we can rely on CxbxGetVertexDeclaration always being set
int CountActiveD3DStreams()
{
	int StreamCount = 0;
	for (int XboxStreamNumber = 0; XboxStreamNumber < X_VSH_MAX_STREAMS; XboxStreamNumber++) {
		if (GetXboxVertexStreamInput(XboxStreamNumber).VertexBuffer != xbox::zeroptr) {
			StreamCount++;
		}
	}

	return StreamCount;
}

UINT CxbxVertexBufferConverter::GetNbrStreams(CxbxDrawContext *pDrawContext) const
{
	// Draw..Up always have one stream
	if (pDrawContext->pXboxVertexStreamZeroData != xbox::zeroptr) {
		return 1;
	}

	CxbxVertexDeclaration *pDecl = CxbxGetVertexDeclaration();
	if (pDecl) {
		return pDecl->NumberOfVertexStreams;
    } 
	
	// TODO: This code and CountActiveD3DStreams must be removed once we can rely on CxbxGetVertexDeclaration always being set
	if (g_Xbox_VertexShader_Handle) {
		return CountActiveD3DStreams();
    }

    return 0;
}

inline FLOAT PackedIntToFloat(const int value, const FLOAT PosFactor, const FLOAT NegFactor)
{
	if (value >= 0) {
		return ((FLOAT)value) / PosFactor;
	}
	else {
		return ((FLOAT)value) / NegFactor;
	}
}

CxbxPatchedStream& CxbxVertexBufferConverter::GetPatchedStream(uint64_t dataKey, uint64_t streamInfoKey)
{
    // First, attempt to fetch an existing patched stream
    const StreamKey key{ dataKey, streamInfoKey };

    auto it = m_PatchedStreams.find(key);
    if (it != m_PatchedStreams.end()) {
        m_TotalLookupSuccesses++;
        m_PatchedStreamUsageList.splice(m_PatchedStreamUsageList.begin(), m_PatchedStreamUsageList, it->second);
        return *it->second;
    }

    // We didn't find an existing patched stream, so we must insert one and get a reference to it
    m_PatchedStreamUsageList.push_front({});
    CxbxPatchedStream& stream = m_PatchedStreamUsageList.front();

    // Insert a reference iterator into the fast lookup map
    m_PatchedStreams[key] = m_PatchedStreamUsageList.begin();

    // If the cache has exceeded it's upper bound, discard the oldest entries in the cache
    if (m_PatchedStreams.size() > (m_MaxCacheSize + m_CacheElasticity)) {
        while (m_PatchedStreams.size() > m_MaxCacheSize) {
            const CxbxPatchedStream& streamToDelete = m_PatchedStreamUsageList.back();

            m_PatchedStreams.erase({ streamToDelete.uiVertexDataHash, streamToDelete.uiVertexStreamInformationHash });
            m_PatchedStreamUsageList.pop_back();
        }
    }
    
    return stream;
}

void CxbxVertexBufferConverter::DrawCacheStats()
{
	const ULONG falsePositives = std::exchange(m_TotalLookupSuccesses, 0) - m_TotalCacheHits;
	const ULONG totalMisses = m_VertexStreamHashMisses + m_DataNotInCacheMisses;

	ImGui::Text("Cache Size: %u", m_PatchedStreams.size());
	ImGui::Text("Hits: %u", std::exchange(m_TotalCacheHits, 0));
	ImGui::Text("Total misses: %u", totalMisses);
	ImGui::Separator();
	ImGui::TextUnformatted("Cache miss details:");
	ImGui::TextWrapped("Vertex stream hash miss: %u", std::exchange(m_VertexStreamHashMisses, 0));
	ImGui::TextWrapped("Data not in cache: %u", std::exchange(m_DataNotInCacheMisses, 0));
}

void CxbxVertexBufferConverter::ConvertStream
(
    CxbxDrawContext *pDrawContext,
    CxbxVertexDeclaration* pCxbxVertexDeclaration,
    UINT             uiStream
)
{
	//X_D3DBaseTexture *pLinearBaseTexture[xbox::X_D3DTS_STAGECOUNT];

	CxbxVertexShaderStreamInfo *pVertexShaderStreamInfo = nullptr;
	UINT XboxStreamNumber = uiStream;
	if (pCxbxVertexDeclaration != nullptr) {
		if (uiStream > pCxbxVertexDeclaration->NumberOfVertexStreams) {
			LOG_TEST_CASE("uiStream > NumberOfVertexStreams");
			return;
		}

		pVertexShaderStreamInfo = &(pCxbxVertexDeclaration->VertexStreams[uiStream]);
		XboxStreamNumber = pVertexShaderStreamInfo->XboxStreamIndex;
	}

	bool bNeedVertexPatching = (pVertexShaderStreamInfo != nullptr && pVertexShaderStreamInfo->NeedPatch);
	bool bNeedStreamCopy = bNeedVertexPatching;

	UINT HostStreamNumber = XboxStreamNumber; // Use Xbox stream index on host
	uint8_t *pXboxVertexData = xbox::zeroptr;
	UINT uiXboxVertexStride = 0;
	UINT uiHostVertexStride = 0;
	uint8_t *pHostVertexData = nullptr;
	ID3D11Buffer *pNewHostVertexBuffer = nullptr;

    if (pDrawContext->pXboxVertexStreamZeroData != xbox::zeroptr) {
		// There should only be one stream (stream zero) in this case
		if (XboxStreamNumber != 0) {
			CxbxrAbort("Trying to patch a Draw..UP with more than stream zero!");
		}

		pXboxVertexData = (uint8_t *)pDrawContext->pXboxVertexStreamZeroData;
		uiXboxVertexStride = pDrawContext->uiXboxVertexStreamZeroStride;
		uiHostVertexStride = (bNeedVertexPatching) ? pVertexShaderStreamInfo->HostVertexStride : uiXboxVertexStride;
	} else {
		xbox::X_STREAMINPUT& XboxStreamInput = GetXboxVertexStreamInput(XboxStreamNumber);
		xbox::X_D3DVertexBuffer *pXboxVertexBuffer = XboxStreamInput.VertexBuffer;
        pXboxVertexData = (uint8_t*)GetDataFromXboxResource(pXboxVertexBuffer);
		if (pXboxVertexData == xbox::zeroptr) {
			HRESULT hRet = CxbxSetStreamSource(
				HostStreamNumber,
				nullptr, 
				0);
			if (FAILED(hRet)) {
				EmuLog(LOG_LEVEL::WARNING, "CxbxSetStreamSource(HostStreamNumber, nullptr, 0)");
			}

			return;
		}

		pXboxVertexData += XboxStreamInput.Offset;
		uiXboxVertexStride = XboxStreamInput.Stride;
		// Dxbx note : Don't overwrite pDrawContext.dwVertexCount with uiVertexCount, because an indexed draw
		// can (and will) use less vertices than the supplied nr of indexes. Thix fixes
		// the missing parts in the CompressedVertices sample (in Vertex shader mode).

		uiHostVertexStride = (bNeedVertexPatching) ? pVertexShaderStreamInfo->HostVertexStride : uiXboxVertexStride;

		// Copy stream for patching and caching.
		bNeedStreamCopy = true;
    }

    // FAST PATH: If this draw is a zerostream based draw, and does not require patching, we can use it directly
    // No need to hash or patch at all in this case!
    if (pDrawContext->pXboxVertexStreamZeroData != xbox::zeroptr && !bNeedStreamCopy) {
        pHostVertexData = pXboxVertexData;

        CxbxPatchedStream stream;
        stream.isValid = true;
        stream.XboxPrimitiveType = pDrawContext->XboxPrimitiveType;
        stream.uiCachedHostVertexStride = uiHostVertexStride;
        stream.bCacheIsStreamZeroDrawUP = true;
        stream.pCachedHostVertexStreamZeroData = pHostVertexData;
        stream.Activate(pDrawContext, HostStreamNumber);
        return;
    }

    // Now we have enough information to hash the existing resource and find it in our cache!
    // To avoid hashing and converting unused vertices, identify the "interesting" region
    // basing on the index/starting vertex data 
    if (pDrawContext->pXboxIndexData != nullptr) {
        pXboxVertexData += (pDrawContext->dwBaseVertexIndex + pDrawContext->LowIndex) * uiXboxVertexStride;
    } else {
        pXboxVertexData += pDrawContext->dwStartVertex * uiXboxVertexStride;
    }

    const UINT uiVertexCount = pDrawContext->NumVerticesToUse;
    const DWORD dwHostVertexDataSize = uiVertexCount * uiHostVertexStride;
    const DWORD xboxVertexDataSize = uiVertexCount * uiXboxVertexStride;
    const uint64_t vertexDataHash = ComputeHash(pXboxVertexData, xboxVertexDataSize);
    const uint64_t pVertexShaderSteamInfoHash = pVertexShaderStreamInfo != nullptr ? ComputeHash(pVertexShaderStreamInfo->VertexElements,
            sizeof(pVertexShaderStreamInfo->VertexElements[0]) * pVertexShaderStreamInfo->NumberOfVertexElements) : 0;

    // Lookup implicity inserts a new entry if not exists, so this always works
    CxbxPatchedStream& patchedStream = GetPatchedStream(vertexDataHash, pVertexShaderSteamInfoHash);

    // We check a few fields of the patched stream to protect against hash collisions (rare)
    // but also to protect against games using the exact same vertex data for different vertex formats (Test Case: Burnout)
    if (patchedStream.isValid && // Check that we found a cached stream
        patchedStream.uiCachedHostVertexStride == patchedStream.uiCachedHostVertexStride && // Make sure the host stride didn't change
        patchedStream.uiCachedXboxVertexStride == uiXboxVertexStride && // Make sure the Xbox Stride didn't change
        patchedStream.uiCachedXboxVertexDataSize == xboxVertexDataSize ) { // Make sure the Xbox Data Size also didn't change
        m_TotalCacheHits++;
        patchedStream.Activate(pDrawContext, HostStreamNumber);
        return;
    }

	// Gather stats
    if (patchedStream.uiVertexStreamInformationHash != pVertexShaderSteamInfoHash)
		m_VertexStreamHashMisses++;
	else
		m_DataNotInCacheMisses++;

    // If execution reaches here, the cached vertex buffer was not valid and we must reconvert the data
    // Free the existing buffers
	patchedStream.Clear();
    assert(pHostVertexData == nullptr);
	assert(pNewHostVertexBuffer == nullptr);

	// If dwHostVertexDataSize is zero, the allocation/creation will fail
	// This can be caused by a stride of 0, and 'other' invalid configurations
	// Test Case :SSX series of games
	if (dwHostVertexDataSize == 0) {
		LOG_TEST_CASE("Attempted to use a 0 sized vertex stream");
		return;
	}

	// Try GPU vertex conversion via compute shader (non-UP, patching case only)
	if (bNeedVertexPatching && pDrawContext->pXboxVertexStreamZeroData == xbox::zeroptr) {
		// Check if all host element sizes are multiples of 4 (required for RWBuffer<uint> writes)
		bool canUseCS = (uiHostVertexStride & 3) == 0;
		UINT elemDescs[16 * 4]; // Up to 16 elements, 4 uints each
		UINT srcOff = 0, dstOff = 0;
		UINT numElems = pVertexShaderStreamInfo->NumberOfVertexElements;
		if (numElems > 16) canUseCS = false;
		if (canUseCS) {
			for (UINT e = 0; e < numElems; e++) {
				UINT hostSize = pVertexShaderStreamInfo->VertexElements[e].HostByteSize;
				UINT xboxSize = pVertexShaderStreamInfo->VertexElements[e].XboxByteSize;
				// All host element sizes and offsets must be 4-byte aligned for CS writes
				if ((hostSize & 3) != 0 || (dstOff & 3) != 0) {
					canUseCS = false;
					break;
				}
				UINT convType, copyDwords = 0;
				switch (pVertexShaderStreamInfo->VertexElements[e].XboxType) {
				case xbox::X_D3DVSDT_NORMSHORT3:  convType = CXBX_VTXCONV_NORMSHORT3; break;
				case xbox::X_D3DVSDT_NORMPACKED3: convType = CXBX_VTXCONV_NORMPACKED3; break;
				case xbox::X_D3DVSDT_SHORT3:      convType = CXBX_VTXCONV_SHORT3; break;
				case xbox::X_D3DVSDT_PBYTE3:      convType = CXBX_VTXCONV_PBYTE3; break;
				case xbox::X_D3DVSDT_FLOAT2H:     convType = CXBX_VTXCONV_FLOAT2H; break;
				case xbox::X_D3DVSDT_D3DCOLOR:    convType = CXBX_VTXCONV_D3DCOLOR; break;
				case xbox::X_D3DVSDT_NONE:        convType = CXBX_VTXCONV_NONE; break;
				default:                          convType = CXBX_VTXCONV_COPY; copyDwords = hostSize / 4; break;
				}
				elemDescs[e * 4 + 0] = srcOff;
				elemDescs[e * 4 + 1] = dstOff;
				elemDescs[e * 4 + 2] = convType;
				elemDescs[e * 4 + 3] = copyDwords;
				srcOff += xboxSize;
				dstOff += hostSize;
			}
		}
		if (canUseCS) {
			ID3D11Buffer* pGPUVB = nullptr;
			if (CxbxD3D11ConvertVertexBufferGPU(
					pXboxVertexData, xboxVertexDataSize, uiVertexCount,
					uiXboxVertexStride, uiHostVertexStride,
					numElems, elemDescs, dwHostVertexDataSize, &pGPUVB)) {
				patchedStream.isValid = true;
				patchedStream.XboxPrimitiveType = pDrawContext->XboxPrimitiveType;
				patchedStream.pCachedXboxVertexData = pXboxVertexData;
				patchedStream.uiCachedXboxVertexDataSize = xboxVertexDataSize;
				patchedStream.uiVertexDataHash = vertexDataHash;
				patchedStream.uiVertexStreamInformationHash = pVertexShaderSteamInfoHash;
				patchedStream.uiCachedXboxVertexStride = uiXboxVertexStride;
				patchedStream.uiCachedHostVertexStride = uiHostVertexStride;
				patchedStream.bCacheIsStreamZeroDrawUP = false;
				patchedStream.pCachedHostVertexBuffer = pGPUVB;
				patchedStream.Activate(pDrawContext, HostStreamNumber);
				return;
			}
		}
	}

    // Allocate new buffers
    if (pDrawContext->pXboxVertexStreamZeroData != xbox::zeroptr) {
        pHostVertexData = (uint8_t*)malloc(dwHostVertexDataSize);

        if (pHostVertexData == nullptr) {
            CxbxrAbort("Couldn't allocate the new stream zero buffer");
        }
    } else {
   	   	HRESULT hRet = CxbxCreateVertexBuffer(dwHostVertexDataSize, &pNewHostVertexBuffer);

        if (FAILED(hRet)) {
            CxbxrAbort("Failed to create vertex buffer");
        }
    }

    // If we need to lock a host vertex buffer, do so now
    if (pHostVertexData == nullptr && pNewHostVertexBuffer != nullptr) {
   	   	pHostVertexData = (uint8_t*)CxbxLockVertexBuffer(pNewHostVertexBuffer);
   	   	if (pHostVertexData == nullptr) {
            CxbxrAbort("Couldn't lock vertex buffer");
        }
    }
	
	if (bNeedVertexPatching) {
	    // assert(bNeedStreamCopy || "bNeedVertexPatching implies bNeedStreamCopy (but copies via conversions");
		for (uint32_t uiVertex = 0; uiVertex < uiVertexCount; uiVertex++) {
			uint8_t *pXboxVertexAsByte = &pXboxVertexData[uiVertex * uiXboxVertexStride];
			uint8_t *pHostVertexAsByte = &pHostVertexData[uiVertex * uiHostVertexStride];
			for (UINT uiElement = 0; uiElement < pVertexShaderStreamInfo->NumberOfVertexElements; uiElement++) {
				FLOAT *pXboxVertexAsFloat = (FLOAT*)pXboxVertexAsByte;
				SHORT *pXboxVertexAsShort = (SHORT*)pXboxVertexAsByte;
				const int XboxElementByteSize = pVertexShaderStreamInfo->VertexElements[uiElement].XboxByteSize;
				FLOAT *pHostVertexAsFloat = (FLOAT*)pHostVertexAsByte;
				SHORT *pHostVertexAsShort = (SHORT*)pHostVertexAsByte;
				// Dxbx note : The following code handles only the D3DVSDT enums that need conversion;
				// All other cases are catched by the memcpy in the default-block.
				switch (pVertexShaderStreamInfo->VertexElements[uiElement].XboxType) {
				case xbox::X_D3DVSDT_NORMSHORT3: { // 0x31:
					// Make it SHORT4N
					pHostVertexAsShort[0] = pXboxVertexAsShort[0];
					pHostVertexAsShort[1] = pXboxVertexAsShort[1];
					pHostVertexAsShort[2] = pXboxVertexAsShort[2];
					pHostVertexAsShort[3] = 32767; // TODO : verify
					break;
				}
				case xbox::X_D3DVSDT_NORMPACKED3: { // 0x16:
					// Make it FLOAT3
					union {
                        int32_t value;
						struct {
							int x : 11;
							int y : 11;
							int z : 10;
						};
					} NormPacked3;

					NormPacked3.value = ((int32_t*)pXboxVertexAsByte)[0];

					pHostVertexAsFloat[0] = PackedIntToFloat(NormPacked3.x, 1023.0f, 1024.f);
					pHostVertexAsFloat[1] = PackedIntToFloat(NormPacked3.y, 1023.0f, 1024.f);
					pHostVertexAsFloat[2] = PackedIntToFloat(NormPacked3.z, 511.0f, 512.f);
					break;
				}
				case xbox::X_D3DVSDT_SHORT3: { // 0x35:
					// Make it a SHORT4 and set the fourth short to 1
					pHostVertexAsShort[0] = pXboxVertexAsShort[0];
					pHostVertexAsShort[1] = pXboxVertexAsShort[1];
					pHostVertexAsShort[2] = pXboxVertexAsShort[2];
					pHostVertexAsShort[3] = 1; // Turok verified (character disappears when this is 32767)
					break;
				}
				case xbox::X_D3DVSDT_PBYTE3: { // 0x34:
					// Make it UBYTE4N
					pHostVertexAsByte[0] = pXboxVertexAsByte[0];
					pHostVertexAsByte[1] = pXboxVertexAsByte[1];
					pHostVertexAsByte[2] = pXboxVertexAsByte[2];
					pHostVertexAsByte[3] = 255; // TODO : Verify
					break;
				}
				case xbox::X_D3DVSDT_FLOAT2H: { // 0x72:
					// Make it FLOAT4 and set the third float to 0.0
					pHostVertexAsFloat[0] = pXboxVertexAsFloat[0];
					pHostVertexAsFloat[1] = pXboxVertexAsFloat[1];
					pHostVertexAsFloat[2] = 0.0f;
					pHostVertexAsFloat[3] = pXboxVertexAsFloat[2];
					break;
				}
				case xbox::X_D3DVSDT_NONE: { // 0x02:
					// Test-case : WWE RAW2
					// Test-case : PetitCopter 
					LOG_TEST_CASE("X_D3DVSDT_NONE");
					// No host element data (but Xbox size can be above zero, when used for X_D3DVSD_MASK_SKIP*
					break;
				}
				case xbox::X_D3DVSDT_D3DCOLOR: { // 0x40: DXGI_FORMAT_R8G8B8A8_UNORM
					// D3DCOLOR is stored as [B,G,R,A] in memory. We use DXGI_FORMAT_R8G8B8A8_UNORM
					// which reads bytes as [R,G,B,A], so swap bytes 0 (B) and 2 (R) to get correct RGBA.
					pHostVertexAsByte[0] = pXboxVertexAsByte[2]; // R
					pHostVertexAsByte[1] = pXboxVertexAsByte[1]; // G
					pHostVertexAsByte[2] = pXboxVertexAsByte[0]; // B
					pHostVertexAsByte[3] = pXboxVertexAsByte[3]; // A
					break;
				}
				case xbox::X_D3DVSDT_FLOAT1: [[fallthrough]]; // 0x12: DXGI_FORMAT_R32_FLOAT
				case xbox::X_D3DVSDT_FLOAT2: [[fallthrough]]; // 0x22: DXGI_FORMAT_R32G32_FLOAT
				case xbox::X_D3DVSDT_FLOAT3: [[fallthrough]]; // 0x32: DXGI_FORMAT_R32G32B32_FLOAT
				case xbox::X_D3DVSDT_FLOAT4: [[fallthrough]]; // 0x42: DXGI_FORMAT_R32G32B32A32_FLOAT
				case xbox::X_D3DVSDT_NORMSHORT1: [[fallthrough]]; // 0x11: DXGI_FORMAT_R16_SNORM
				case xbox::X_D3DVSDT_NORMSHORT2: [[fallthrough]]; // 0x21: DXGI_FORMAT_R16G16_SNORM
				case xbox::X_D3DVSDT_NORMSHORT4: [[fallthrough]]; // 0x41: DXGI_FORMAT_R16G16B16A16_SNORM
				case xbox::X_D3DVSDT_PBYTE1: [[fallthrough]]; // 0x14: DXGI_FORMAT_R8_UNORM
				case xbox::X_D3DVSDT_PBYTE2: [[fallthrough]]; // 0x24: DXGI_FORMAT_R8G8_UNORM
				case xbox::X_D3DVSDT_PBYTE4: [[fallthrough]]; // 0x44: DXGI_FORMAT_R8G8B8A8_UNORM
				case xbox::X_D3DVSDT_SHORT1: [[fallthrough]]; // 0x15: DXGI_FORMAT_R16_SINT
				case xbox::X_D3DVSDT_SHORT2: [[fallthrough]]; // 0x25: DXGI_FORMAT_R16G16_SINT
				case xbox::X_D3DVSDT_SHORT4: [[fallthrough]]; // 0x45: DXGI_FORMAT_R16G16B16A16_SINT
				default: {
					// Generic 'conversion' - just make a copy :
					memcpy(pHostVertexAsByte, pXboxVertexAsByte, XboxElementByteSize);
					break;
				}
				} // switch

				// Increment the Xbox pointer :
				pXboxVertexAsByte += XboxElementByteSize;
				// Increment the host pointer :
				pHostVertexAsByte += pVertexShaderStreamInfo->VertexElements[uiElement].HostByteSize;
			} // for NumberOfVertexElements
		} // for uiVertexCount
    }
    else {
		if (bNeedStreamCopy) {
			memcpy(pHostVertexData, pXboxVertexData, dwHostVertexDataSize);
		}
	}

    patchedStream.isValid = true;
    patchedStream.XboxPrimitiveType = pDrawContext->XboxPrimitiveType;
    patchedStream.pCachedXboxVertexData = pXboxVertexData;
    patchedStream.uiCachedXboxVertexDataSize = xboxVertexDataSize;
    patchedStream.uiVertexDataHash = vertexDataHash;
    patchedStream.uiVertexStreamInformationHash = pVertexShaderSteamInfoHash;
    patchedStream.uiCachedXboxVertexStride = uiXboxVertexStride;
    patchedStream.uiCachedHostVertexStride = uiHostVertexStride;
    patchedStream.bCacheIsStreamZeroDrawUP = (pDrawContext->pXboxVertexStreamZeroData != xbox::zeroptr);
    if (patchedStream.bCacheIsStreamZeroDrawUP) {
        patchedStream.pCachedHostVertexStreamZeroData = pHostVertexData;
        patchedStream.bCachedHostVertexStreamZeroDataIsAllocated = bNeedStreamCopy;
    } else {
        // assert(pNewHostVertexBuffer != nullptr);
   	   	CxbxUnlockVertexBuffer(pNewHostVertexBuffer);
        patchedStream.pCachedHostVertexBuffer = pNewHostVertexBuffer;
    }

	patchedStream.Activate(pDrawContext, HostStreamNumber);
}

void CxbxVertexBufferConverter::Apply(CxbxDrawContext *pDrawContext)
{
	if ((pDrawContext->XboxPrimitiveType < xbox::X_D3DPT_POINTLIST) || (pDrawContext->XboxPrimitiveType > xbox::X_D3DPT_POLYGON))
		CxbxrAbort("Unknown primitive type: 0x%.02X\n", pDrawContext->XboxPrimitiveType);

	CxbxVertexDeclaration* pCxbxVertexDeclaration = CxbxGetVertexDeclaration();

	// When this is an indexed draw, take the index buffer into account
	if (pDrawContext->pXboxIndexData) {
		// Is the highest index in this buffer not set yet?
		if (pDrawContext->HighIndex == 0) {
			// TODO : Instead of calling WalkIndexBuffer here, set LowIndex and HighIndex
			// in all callers that end up here (since they might be able to avoid the call)
			LOG_TEST_CASE("HighIndex == 0"); // TODO : If this is never hit, replace entire block by assert(pDrawContext->HighIndex > 0);
			WalkIndexBuffer(pDrawContext->LowIndex, pDrawContext->HighIndex, pDrawContext->pXboxIndexData, pDrawContext->dwVertexCount);
		}
		// Convert the range of indices into a count
		pDrawContext->NumVerticesToUse = pDrawContext->HighIndex - pDrawContext->LowIndex + 1;
	}
	else {
		// If we are drawing from an offset, we know that the vertex count must have
		// 'offset' vertices before the first drawn vertices
		pDrawContext->NumVerticesToUse = pDrawContext->dwVertexCount;
	}

    // Get the number of streams
    UINT nbrStreams = GetNbrStreams(pDrawContext);
    if (nbrStreams > X_VSH_MAX_STREAMS) {
        LOG_TEST_CASE("nbrStreams count > max number of streams");
        nbrStreams = X_VSH_MAX_STREAMS;
    }

    for(UINT i = 0; i < nbrStreams; i++) {
		ConvertStream(pDrawContext, pCxbxVertexDeclaration, i);
    }

	if (pDrawContext->XboxPrimitiveType == xbox::X_D3DPT_QUADSTRIP) {
		// Quad strip is just like a triangle strip, but requires two vertices per primitive.
		// A quadstrip starts with 4 vertices and adds 2 vertices per additional quad.
		// This is much like a trianglestrip, which starts with 3 vertices and adds
		// 1 vertex per additional triangle, so we use that instead. The planar nature
		// of the quads 'survives' through this change. There's a catch though :
		// In a trianglestrip, every 2nd triangle has an opposing winding order,
		// which would cause backface culling - but this seems to be intelligently
		// handled by d3d :
		// Test-case : XDK Samples (FocusBlur, MotionBlur, Trees, PaintEffect, PlayField)
		// No need to set : pDrawContext->XboxPrimitiveType = X_D3DPT_TRIANGLESTRIP;
		pDrawContext->dwHostPrimitiveCount = ConvertXboxVertexCountToPrimitiveCount(xbox::X_D3DPT_TRIANGLESTRIP, pDrawContext->dwVertexCount);
	} else {
		pDrawContext->dwHostPrimitiveCount = ConvertXboxVertexCountToPrimitiveCount(pDrawContext->XboxPrimitiveType, pDrawContext->dwVertexCount);
	}

	if (pDrawContext->XboxPrimitiveType == xbox::X_D3DPT_POLYGON) {
		// Convex polygon is the same as a triangle fan.
		// No need to set : pDrawContext->XboxPrimitiveType = X_D3DPT_TRIANGLEFAN;
		// Test-case : Panzer Dragoon ORTA (when entering in-game)
		LOG_TEST_CASE("X_D3DPT_POLYGON");
	}
}

void CxbxSetVertexAttribute(int Register, FLOAT a, FLOAT b, FLOAT c, FLOAT d)
{
	if (Register < 0) {
		LOG_TEST_CASE("Register < 0");
		return;
	}
	if (Register >= 16) {
		LOG_TEST_CASE("Register >= 16");
		return;
	}

	// Write these values to the NV2A registers, so that we read them back when needed
	float* attribute_floats = NV2A_get_vertex_attribute_value_pointer(Register);
	attribute_floats[0] = a;
	attribute_floats[1] = b;
	attribute_floats[2] = c;
	attribute_floats[3] = d;

	g_bD3D11VertexFetchDefaultsDirty = true;

	// D3D11: The zero-stride vertex defaults buffer reads inline_value[] directly,
	// so no constant buffer upload is needed for attribute defaults.
}

void CxbxImpl_SetStreamSource(UINT StreamNumber, xbox::X_D3DVertexBuffer* pStreamData, UINT Stride)
{
	if (pStreamData != xbox::zeroptr && Stride == 0) {
		LOG_TEST_CASE("CxbxImpl_SetStreamSource : Stream assigned, and stride set to 0 (might be okay)");
	}

	assert(StreamNumber < X_VSH_MAX_STREAMS);

	g_Xbox_SetStreamSource[StreamNumber].VertexBuffer = pStreamData;
	g_Xbox_SetStreamSource[StreamNumber].Stride = Stride;

	CxbxD3D11VertexFetchInvalidateLayout();
}
