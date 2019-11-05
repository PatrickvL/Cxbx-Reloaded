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
#include "core\hle\D3D8\Direct3D9\Direct3D9.h" // For g_pD3DDevice, g_pXbox_Texture, g_pXbox_RenderTarget, g_pXbox_BackBufferSurface, g_Xbox_VertexShader_Handle
#include "core\hle\D3D8\Direct3D9\WalkIndexBuffer.h" // for WalkIndexBuffer
#include "core\hle\D3D8\ResourceTracker.h"
#include "core\hle\D3D8\XbPushBuffer.h" // for XboxFVF_GetNumberOfTextureCoordinates
#include "core\hle\D3D8\XbVertexBuffer.h"
#include "core\hle\D3D8\XbConvert.h"

#include <ctime>
#include <chrono>
#include <algorithm>

#define MAX_STREAM_NOT_USED_TIME (2 * CLOCKS_PER_SEC) // TODO: Trim the not used time

// Inline vertex buffer emulation
extern XTL::X_D3DPRIMITIVETYPE g_InlineVertexBuffer_PrimitiveType = XTL::X_D3DPT_INVALID;
extern DWORD                   g_InlineVertexBuffer_FVF = 0;
extern struct _D3DIVB         *g_InlineVertexBuffer_Table = nullptr;
extern UINT                    g_InlineVertexBuffer_TableLength = 0;
extern UINT                    g_InlineVertexBuffer_TableOffset = 0;

FLOAT *g_InlineVertexBuffer_pData = nullptr;
UINT   g_InlineVertexBuffer_DataSize = 0;

extern DWORD				g_dwPrimPerFrame = 0;

// Copy of active Xbox D3D Vertex Streams (and strides), set by [D3DDevice|CxbxImpl]_SetStreamSource*
XTL::X_STREAMINPUT g_Xbox_SetStreamSource[X_VSH_MAX_STREAMS] = { 0 }; // Note : .Offset member is never set (so always 0)

// Resource related functions, declared in Direct3D9.cpp :
extern void *GetDataFromXboxResource(XTL::X_D3DResource *pXboxResource);
extern bool GetHostRenderTargetDimensions(DWORD* pHostWidth, DWORD* pHostHeight, IDirect3DSurface* pHostRenderTarget = nullptr);
extern uint32_t GetPixelContainerWidth(XTL::X_D3DPixelContainer* pPixelContainer);
extern uint32_t GetPixelContainerHeight(XTL::X_D3DPixelContainer* pPixelContainer);

// implemented in XbVertexShader.cpp
extern void CxbxUpdateActiveVertexShader(unsigned VerticesInBuffer);

// Reads the active Xbox stream input values (containing VertexBuffer, Offset and Stride) for the given stream number.
// (These values are set through SetStreamSource and can be overridden by SetVertexShaderInput.)
XTL::X_STREAMINPUT& GetXboxVertexStreamInput(unsigned StreamNumber)
{
/* TODO : Enable this once we support SetVertexShaderInput completely (that is, when we convert both vertex buffers AND vertex shaders at drawing time)
	// If SetVertexShaderInput is active, it's arguments overrule those of SetStreamSource
	if (g_Xbox_SetVertexShaderInput_Count > 0) {
		return g_Xbox_SetVertexShaderInput_Data[StreamNumber];
	}
*/
	return g_Xbox_SetStreamSource[StreamNumber];
}

// Returns the Xbox vertex stride from the stream input values that are active on the given stream number.
unsigned GetXboxVertexStreamStride(unsigned StreamNumber)
{
	return GetXboxVertexStreamInput(StreamNumber).Stride;
}

// Returns the Xbox vertex buffer from the stream input values that are active on the given stream number.
XTL::X_D3DVertexBuffer* GetXboxVertexStreamVertexBuffer(unsigned StreamNumber)
{
	return GetXboxVertexStreamInput(StreamNumber).VertexBuffer;
}

// Returns a CPU-readable Xbox data pointer from the given vertex buffer.
uint8_t* GetXboxVertexBufferData(XTL::X_D3DVertexBuffer *pXboxVertexBuffer)
{
	return (uint8_t*)GetDataFromXboxResource(pXboxVertexBuffer);
}

// Returns a CPU-readable Xbox data pointer from the vertex buffer that's active on the given stream number.
uint8_t* GetXboxVertexStreamData(unsigned StreamNumber)
{
	return GetXboxVertexBufferData(GetXboxVertexStreamVertexBuffer(StreamNumber));
}

void CxbxPatchedStream::Activate(CxbxDrawContext *pDrawContext, UINT uiStream) const
{
	//LOG_INIT // Allows use of DEBUG_D3DRESULT

	assert(m_pCachedHostVertexStreamZeroData == nullptr);

	// TODO : Use g_pDummyBuffer when m_pCachedHostVertexBuffer == nullptr?
	// Luke contemplated g_pDummyBuffer might have been intended for cases where Xbox
	// has null streams in between assigned streams, which host might not support.
	// (If not, perhaps we should remove g_pDummyBuffer entirely.)

	// TODO : IF we use g_pDummyBuffer, should we call SetStreamSourceFreq() too,
	// so that it won't run out of it's buffer (although that won't help with
	// indexed draws...)?

	HRESULT hRet = g_pD3DDevice->SetStreamSource(
		/*StreamNumber=*/uiStream,
		/*pStreamData=*/m_pCachedHostVertexBuffer,
		/*OffsetInBytes=*/0,
		/*Stride=*/m_uiCachedHostVertexStride);
	//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetStreamSource");
	if (FAILED(hRet)) {
		CxbxKrnlCleanup("Failed to set the type patched buffer as the new stream source!\n");
		// TODO : test-case : XDK Cartoon hits the above case when the vertex cache size is 0.
	}
}

CxbxPatchedStream::CxbxPatchedStream()
{
    m_bIsCached = false;
}

CxbxPatchedStream::~CxbxPatchedStream()
{
	Clear();
}

void CxbxPatchedStream::Clear()
{
    if (m_bCachedHostVertexStreamZeroDataIsAllocated) {
        free(m_pCachedHostVertexStreamZeroData);
        m_bCachedHostVertexStreamZeroDataIsAllocated = false;
    }

    m_pCachedHostVertexStreamZeroData = nullptr;

    if (m_pCachedHostVertexBuffer != nullptr) {
        // TODO : We could re-use this buffer, if it's large enough...
        m_pCachedHostVertexBuffer->Release();
        m_pCachedHostVertexBuffer = nullptr;
    }
}

CxbxVertexBufferConverter::CxbxVertexBufferConverter()
{
    m_pVertexShaderInfo = nullptr;
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

inline FLOAT NormShortToFloat(const SHORT value)
{
	return PackedIntToFloat((int)value, 32767.0f, 32768.0f);
}

inline FLOAT ByteToFloat(const BYTE value)
{
	return ((FLOAT)value) / 255.0f;
}

CxbxPatchedStream& CxbxVertexBufferConverter::GetPatchedStream(uint64_t key)
{
    // First, attempt to fetch an existing patched stream
    auto it = m_PatchedStreams.find(key);
    if (it != m_PatchedStreams.end()) {
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
            m_PatchedStreams.erase(m_PatchedStreamUsageList.back().m_uiCachedXboxVertexDataHash);
            m_PatchedStreamUsageList.pop_back();
        }
    }
    
    return stream;
}

void CxbxVertexBufferConverter::PrintStats()
{
    printf("Vertex Buffer Cache Status: \n");
    printf("- Cache Size: %d\n", m_PatchedStreams.size());
    printf("- Hits: %d\n", m_TotalCacheHits);
    printf("- Misses: %d\n", m_TotalCacheMisses);
}

void CxbxVertexBufferConverter::PrepareStreamConversion()
{
	m_bVshHandleIsFVF = VshHandleIsFVF(g_Xbox_VertexShader_Handle);
	DWORD XboxFVF = m_bVshHandleIsFVF ? g_Xbox_VertexShader_Handle : 0;

	// Texture normalization can only be set for FVF shaders
	bool bNeedTextureNormalization = false;

	m_uiTextureCoordinatesByteOffsetInVertex = 0;
	m_bMustReset_FVF_XYZRHW = false;
	memset(m_FVFInfo_ActivePixelContainer, 0, sizeof(m_FVFInfo_ActivePixelContainer));

	if (m_bVshHandleIsFVF) {
		DWORD dwTexN = (XboxFVF & X_D3DFVF_TEXCOUNT_MASK) >> X_D3DFVF_TEXCOUNT_SHIFT;
		if (dwTexN > XTL::X_D3DTS_STAGECOUNT) {
			LOG_TEST_CASE("FVF,dwTexN > X_D3DTS_STAGECOUNT");
		}

		// Check for active linear textures.
		//X_D3DBaseTexture *pLinearBaseTexture[XTL::X_D3DTS_STAGECOUNT];
		for (unsigned int i = 0; i < XTL::X_D3DTS_STAGECOUNT; i++) {
			// Only normalize coordinates used by the FVF shader :
			if (i + 1 <= dwTexN) {
				m_FVFInfo_ActivePixelContainer[i].NrTexCoords = XboxFVF_GetNumberOfTextureCoordinates(XboxFVF, i);
				// TODO : Use GetXboxBaseTexture()
				XTL::X_D3DBaseTexture *pXboxBaseTexture = g_pXbox_SetTexture[i];
				if (pXboxBaseTexture != xbnullptr) {
					XTL::X_D3DFORMAT XboxFormat = GetXboxPixelContainerFormat(pXboxBaseTexture);
					if (EmuXBFormatIsLinear(XboxFormat)) {
						// This is often hit by the help screen in XDK samples.
						bNeedTextureNormalization = true;
						// Remember linearity, width and height :
						m_FVFInfo_ActivePixelContainer[i].bTexIsLinear = true;
						// TODO : Use DecodeD3DSize or GetPixelContainerWidth + GetPixelContainerHeight
						m_FVFInfo_ActivePixelContainer[i].Width = (pXboxBaseTexture->Size & X_D3DSIZE_WIDTH_MASK) + 1;
						m_FVFInfo_ActivePixelContainer[i].Height = ((pXboxBaseTexture->Size & X_D3DSIZE_HEIGHT_MASK) >> X_D3DSIZE_HEIGHT_SHIFT) + 1;
						// TODO : Support 3D textures
					}
				}
			}
		}

		if (bNeedTextureNormalization) {
			m_uiTextureCoordinatesByteOffsetInVertex = XboxFVFToVertexSizeInBytes(XboxFVF, /*bIncludeTextures=*/false);
		}

		m_bMustReset_FVF_XYZRHW = ((XboxFVF & X_D3DFVF_POSITION_MASK) == X_D3DFVF_XYZRHW);
	}
}

CxbxPatchedStream& CxbxVertexBufferConverter::ConvertStream
(
	CxbxDrawContext *pDrawContext,
	UINT uiStream
)
{
	extern D3DCAPS g_D3DCaps;

	uint8_t *pXboxVertexData = xbnullptr;

	// Does this draw supply vertex data through a User Pointer?
    if (pDrawContext->pXboxVertexStreamZeroData != xbnullptr) {
		// Only use the supplied user pointer for stream zero
		assert(uiStream == 0);
		pXboxVertexData = (uint8_t*)pDrawContext->pXboxVertexStreamZeroData;
	} else {
        pXboxVertexData = GetXboxVertexStreamData(uiStream);
	}

	// Deactivate this stream on host when there's no Xbox vertex data :
	if (pXboxVertexData == xbnullptr) {
		m_TempPatchedStream.m_pCachedHostVertexBuffer = nullptr; // Passed to SetStreamSource
		m_TempPatchedStream.m_uiCachedHostVertexStride = 0; // Passed to SetStreamSource
		return m_TempPatchedStream;
	}

	CxbxVertexShaderStreamInfo *pVertexShaderStreamInfo = nullptr;
	if (m_pVertexShaderInfo != nullptr) {
		// TODO : Just map streams 1-on-1 instead of checking the count (which might have skipped gaps?)
		if (uiStream > m_pVertexShaderInfo->NumberOfVertexStreams + 1) {
			LOG_TEST_CASE("uiStream > NumberOfVertexStreams");
		} else {
			pVertexShaderStreamInfo = &(m_pVertexShaderInfo->VertexStreams[uiStream]);
		}
	}

	bool bNeedVertexPatching = (pVertexShaderStreamInfo != nullptr && pVertexShaderStreamInfo->NeedPatch);
	bool bNeedTextureNormalization = (m_uiTextureCoordinatesByteOffsetInVertex > 0);
	bool bAllocateHostVertexStreamZero = false;

    // FAST PATH: If this draw is a user-pointer based draw, and does not require patching, we can use it directly
    if (pDrawContext->pXboxVertexStreamZeroData != xbnullptr) {
		if (bNeedVertexPatching || bNeedTextureNormalization || m_bMustReset_FVF_XYZRHW) {
			// Some form of conversion is needed, arrange a buffer where we can convert into
			bAllocateHostVertexStreamZero = true;
		} else {
			// No conversion is needed, we just arrange here that the Xbox stream data will be used on host
			m_TempPatchedStream.m_pCachedHostVertexStreamZeroData = pDrawContext->pXboxVertexStreamZeroData; // Pass Xbox stream data to DrawContext, which in turn gets
			m_TempPatchedStream.m_uiCachedHostVertexStride = pDrawContext->uiXboxVertexStreamZeroStride; // ... passed to host DrawPrimitiveUP()/DrawIndexedPrimitiveUP().
			// No need to hash or patch at all in this case!
			return m_TempPatchedStream;
		}
    }

	// Dxbx note : Don't overwrite pDrawContext.dwVertexCount with uiVertexCount, because an indexed draw
	// can (and will) use less vertices than the supplied nr of indexes. Thix fixes
	// the missing parts in the CompressedVertices sample (in Vertex shader mode).
	UINT uiVertexCount = pDrawContext->VerticesInBuffer; 
	UINT uiXboxVertexStride = (pDrawContext->pXboxVertexStreamZeroData != xbnullptr) ? pDrawContext->uiXboxVertexStreamZeroStride : GetXboxVertexStreamStride(uiStream);
	UINT uiHostVertexStride = (bNeedVertexPatching) ? pVertexShaderStreamInfo->HostVertexStride : uiXboxVertexStride;
	DWORD dwHostVertexDataSize = uiVertexCount * uiHostVertexStride;
	uint8_t *pHostVertexData = nullptr;
	IDirect3DVertexBuffer *pHostVertexBuffer = nullptr;

    // Now we have enough information to hash the existing resource and find it in our cache!
    DWORD XboxVertexDataSize = uiVertexCount * uiXboxVertexStride;
    uint64_t XboxVertexDataHash = ComputeHash(pXboxVertexData, XboxVertexDataSize);
    uint64_t pVertexShaderSteamInfoHash = 0;

    if (pVertexShaderStreamInfo != nullptr) {
        pVertexShaderSteamInfoHash = ComputeHash(pVertexShaderStreamInfo, sizeof(CxbxVertexShaderStreamInfo));
    }

    // Lookup implicity inserts a new entry if not exists, so this always works
    CxbxPatchedStream& patchedStream = GetPatchedStream(XboxVertexDataHash);

    // We check a few fields of the patched stream to protect against hash collisions (rare)
    // but also to protect against games using the exact same vertex data for different vertex formats (Test Case: Burnout)
    if (patchedStream.m_bIsCached) { // Check that we found a cached stream
		if (patchedStream.m_uiCachedVertexStreamInformationHash == pVertexShaderSteamInfoHash && // Check that the vertex conversion is valid
			patchedStream.m_uiCachedHostVertexStride == uiHostVertexStride && // Make sure the host stride didn't change
			patchedStream.m_uiCachedXboxVertexStride == uiXboxVertexStride && // Make sure the Xbox stride didn't change
			patchedStream.m_uiCachedXboxVertexDataSize == XboxVertexDataSize ) { // Make sure the Xbox data size also didn't change
			m_TotalCacheHits++;
			return patchedStream;
	    }

	    // If execution reaches here, the cached vertex buffer was not valid and we must reconvert the data
		patchedStream.Clear();
    }

    m_TotalCacheMisses++;

    // Allocate new buffers
    if (bAllocateHostVertexStreamZero) {
		// The fast path hasn't been taken, so we need to convert the vertex stream zero data user pointer
        pHostVertexData = (uint8_t*)malloc(dwHostVertexDataSize);
        if (pHostVertexData == nullptr) {
            CxbxKrnlCleanup("Couldn't allocate the new stream zero buffer");
        }
    } else {
        HRESULT hRet = g_pD3DDevice->CreateVertexBuffer(
            dwHostVertexDataSize,
            D3DUSAGE_WRITEONLY | D3DUSAGE_DYNAMIC,
            0,
            D3DPOOL_DEFAULT,
            &pHostVertexBuffer,
            nullptr
        );
        if (FAILED(hRet)) {
            CxbxKrnlCleanup("Failed to create vertex buffer");
        }

		if (FAILED(pHostVertexBuffer->Lock(0, 0, (D3DLockData **)&pHostVertexData, D3DLOCK_DISCARD))) {
			CxbxKrnlCleanup("Couldn't lock vertex buffer");
        }
    }
	
	if (bNeedVertexPatching) {
	    // copy via conversions
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
				case XTL::X_D3DVSDT_NORMSHORT1: { // 0x11:
					// Test-cases : Halo - Combat Evolved
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
						// Make it SHORT2N
						pHostVertexAsShort[0] = pXboxVertexAsShort[0];
						pHostVertexAsShort[1] = 0;
					}
					else
					{
						// Make it FLOAT1
						pHostVertexAsFloat[0] = NormShortToFloat(pXboxVertexAsShort[0]);
						//pHostVertexAsFloat[1] = 0.0f; // Would be needed for FLOAT2
					}
					break;
				}
				case XTL::X_D3DVSDT_NORMSHORT2: { // 0x21:
					// Test-cases : Baldur's Gate: Dark Alliance 2, F1 2002, Gun, Halo - Combat Evolved, Scrapland 
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT2N) {
						// No need for patching when D3D9 supports D3DDECLTYPE_SHORT2N
						// TODO : goto default; // ??
						//memcpy(pHostVertexAsByte, pXboxVertexAsByte, XboxElementByteSize);
						// Make it SHORT2N
						pHostVertexAsShort[0] = pXboxVertexAsShort[0];
						pHostVertexAsShort[1] = pXboxVertexAsShort[1];
					}
					else
					{
						// Make it FLOAT2
						pHostVertexAsFloat[0] = NormShortToFloat(pXboxVertexAsShort[0]);
						pHostVertexAsFloat[1] = NormShortToFloat(pXboxVertexAsShort[1]);
					}
					break;
				}
				case XTL::X_D3DVSDT_NORMSHORT3: { // 0x31:
					// Test-cases : Cel Damage, Constantine, Destroy All Humans!
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
						// Make it SHORT4N
						pHostVertexAsShort[0] = pXboxVertexAsShort[0];
						pHostVertexAsShort[1] = pXboxVertexAsShort[1];
						pHostVertexAsShort[2] = pXboxVertexAsShort[2];
						pHostVertexAsShort[3] = 32767; // TODO : verify
					}
					else
					{
						// Make it FLOAT3
						pHostVertexAsFloat[0] = NormShortToFloat(pXboxVertexAsShort[0]);
						pHostVertexAsFloat[1] = NormShortToFloat(pXboxVertexAsShort[1]);
						pHostVertexAsFloat[2] = NormShortToFloat(pXboxVertexAsShort[2]);
					}
					break;
				}
				case XTL::X_D3DVSDT_NORMSHORT4: { // 0x41:
					// Test-cases : Judge Dredd: Dredd vs Death, NHL Hitz 2002, Silent Hill 2, Sneakers, Tony Hawk Pro Skater 4
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_SHORT4N) {
						// No need for patching when D3D9 supports D3DDECLTYPE_SHORT4N
						// TODO : goto default; // ??
						//memcpy(pHostVertexAsByte, pXboxVertexAsByte, XboxElementByteSize);
						// Make it SHORT4N
						pHostVertexAsShort[0] = pXboxVertexAsShort[0];
						pHostVertexAsShort[1] = pXboxVertexAsShort[1];
						pHostVertexAsShort[2] = pXboxVertexAsShort[2];
						pHostVertexAsShort[3] = pXboxVertexAsShort[3];
					}
					else
					{
						// Make it FLOAT4
						pHostVertexAsFloat[0] = NormShortToFloat(pXboxVertexAsShort[0]);
						pHostVertexAsFloat[1] = NormShortToFloat(pXboxVertexAsShort[1]);
						pHostVertexAsFloat[2] = NormShortToFloat(pXboxVertexAsShort[2]);
						pHostVertexAsFloat[3] = NormShortToFloat(pXboxVertexAsShort[3]);
					}
					break;
				}
				case XTL::X_D3DVSDT_NORMPACKED3: { // 0x16:
					// Test-cases : Dashboard
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
				case XTL::X_D3DVSDT_SHORT1: { // 0x15:
					// Make it SHORT2 and set the second short to 0
					pHostVertexAsShort[0] = pXboxVertexAsShort[0];
					pHostVertexAsShort[1] = 0;
					break;
				}
				case XTL::X_D3DVSDT_SHORT3: { // 0x35:
					// Test-cases : Turok
					// Make it a SHORT4 and set the fourth short to 1
					pHostVertexAsShort[0] = pXboxVertexAsShort[0];
					pHostVertexAsShort[1] = pXboxVertexAsShort[1];
					pHostVertexAsShort[2] = pXboxVertexAsShort[2];
					pHostVertexAsShort[3] = 1; // Turok verified (character disappears when this is 32767)
					break;
				}
				case XTL::X_D3DVSDT_PBYTE1: { // 0x14:
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
						// Make it UBYTE4N
						pHostVertexAsByte[0] = pXboxVertexAsByte[0];
						pHostVertexAsByte[1] = 0;
						pHostVertexAsByte[2] = 0;
						pHostVertexAsByte[3] = 255; // TODO : Verify
					}
					else
					{
						// Make it FLOAT1
						pHostVertexAsFloat[0] = ByteToFloat(pXboxVertexAsByte[0]);
					}
					break;
				}
				case XTL::X_D3DVSDT_PBYTE2: { // 0x24:
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
						// Make it UBYTE4N
						pHostVertexAsByte[0] = pXboxVertexAsByte[0];
						pHostVertexAsByte[1] = pXboxVertexAsByte[1];
						pHostVertexAsByte[2] = 0;
						pHostVertexAsByte[3] = 255; // TODO : Verify
					}
					else
					{
						// Make it FLOAT2
						pHostVertexAsFloat[0] = ByteToFloat(pXboxVertexAsByte[0]);
						pHostVertexAsFloat[1] = ByteToFloat(pXboxVertexAsByte[1]);
					}
					break;
				}
				case XTL::X_D3DVSDT_PBYTE3: { // 0x34:
					// Test-cases : Turok
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
						// Make it UBYTE4N
						pHostVertexAsByte[0] = pXboxVertexAsByte[0];
						pHostVertexAsByte[1] = pXboxVertexAsByte[1];
						pHostVertexAsByte[2] = pXboxVertexAsByte[2];
						pHostVertexAsByte[3] = 255; // TODO : Verify
					}
					else
					{
						// Make it FLOAT3
						pHostVertexAsFloat[0] = ByteToFloat(pXboxVertexAsByte[0]);
						pHostVertexAsFloat[1] = ByteToFloat(pXboxVertexAsByte[1]);
						pHostVertexAsFloat[2] = ByteToFloat(pXboxVertexAsByte[2]);
					}
					break;
				}
				case XTL::X_D3DVSDT_PBYTE4: { // 0x44:
					// Test-case : Jet Set Radio Future
					if (g_D3DCaps.DeclTypes & D3DDTCAPS_UBYTE4N) {
						// No need for patching when D3D9 supports D3DDECLTYPE_UBYTE4N
						// TODO : goto default; // ??
						//memcpy(pHostVertexAsByte, pXboxVertexAsByte, XboxElementByteSize);
						// Make it UBYTE4N
						pHostVertexAsByte[0] = pXboxVertexAsByte[0];
						pHostVertexAsByte[1] = pXboxVertexAsByte[1];
						pHostVertexAsByte[2] = pXboxVertexAsByte[2];
						pHostVertexAsByte[3] = pXboxVertexAsByte[3];
					}
					else
					{
						// Make it FLOAT4
						pHostVertexAsFloat[0] = ByteToFloat(pXboxVertexAsByte[0]);
						pHostVertexAsFloat[1] = ByteToFloat(pXboxVertexAsByte[1]);
						pHostVertexAsFloat[2] = ByteToFloat(pXboxVertexAsByte[2]);
						pHostVertexAsFloat[3] = ByteToFloat(pXboxVertexAsByte[3]);
					}
					break;
				}
				case XTL::X_D3DVSDT_FLOAT2H: { // 0x72:
					// Make it FLOAT4 and set the third float to 0.0
					pHostVertexAsFloat[0] = pXboxVertexAsFloat[0];
					pHostVertexAsFloat[1] = pXboxVertexAsFloat[1];
					pHostVertexAsFloat[2] = 0.0f;
					pHostVertexAsFloat[3] = pXboxVertexAsFloat[2];
					break;
				}
				case XTL::X_D3DVSDT_NONE: { // 0x02:
					// Test-case : WWE RAW2
					// Test-case : PetitCopter 
					LOG_TEST_CASE("X_D3DVSDT_NONE");
					// No host element data (but Xbox size can be above zero, when used for X_D3DVSD_MASK_SKIP*
					break;
				}
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
    } else { // !bNeedVertexPatching
		memcpy(pHostVertexData, pXboxVertexData, dwHostVertexDataSize);
	}

	// Xbox FVF shaders are identical to host Direct3D 8.1, however
	// texture coordinates may need normalization if used with linear textures.
	if (bNeedTextureNormalization || m_bMustReset_FVF_XYZRHW) {
		// assert(m_bVshHandleIsFVF);

		// Locate texture coordinate offset in vertex structure.
		if (bNeedTextureNormalization) {
			if (bNeedVertexPatching) {
				LOG_TEST_CASE("Potential xbox vs host texture-offset difference! (bNeedVertexPatching within bNeedTextureNormalization)");
			}
			// As long as vertices aren't resized / patched up until the texture coordinates,
			// the m_uiTextureCoordinatesByteOffsetInVertex on host will match Xbox
			// And since FVF's should not specify unsupported vertex attribute types,
			// patching should not be needed. TODO : Verify this.
		}

        // If for some reason the Xbox Render Target is not set, fallback to the backbuffer
        if (g_pXbox_RenderTarget == xbnullptr) {
            LOG_TEST_CASE("SetRenderTarget fallback to backbuffer");
            g_pXbox_RenderTarget = g_pXbox_BackBufferSurface;
        }

		DWORD HostRenderTarget_Width, HostRenderTarget_Height;
		DWORD XboxRenderTarget_Width = GetPixelContainerWidth(g_pXbox_RenderTarget);
		DWORD XboxRenderTarget_Height = GetPixelContainerHeight(g_pXbox_RenderTarget);
		if (!GetHostRenderTargetDimensions(&HostRenderTarget_Width, &HostRenderTarget_Height)) {
			HostRenderTarget_Width = XboxRenderTarget_Width;
			HostRenderTarget_Height = XboxRenderTarget_Height;
		}

		bool bNeedRHWTransform = XboxRenderTarget_Width < HostRenderTarget_Width && XboxRenderTarget_Height < HostRenderTarget_Height;

		for (uint32_t uiVertex = 0; uiVertex < uiVertexCount; uiVertex++) {
			FLOAT *pVertexDataAsFloat = (FLOAT*)(&pHostVertexData[uiVertex * uiHostVertexStride]);

			// Handle pre-transformed vertices (which bypass the vertex shader pipeline)
			if (m_bMustReset_FVF_XYZRHW) {
                // We need to transform these vertices only if the host render target was upscaled from the Xbox render target
                // Transforming always breaks render to non-upscaled textures: Only surfaces are upscaled, intentionally so
				if (bNeedRHWTransform) {
					pVertexDataAsFloat[0] *= g_RenderScaleFactor;
					pVertexDataAsFloat[1] *= g_RenderScaleFactor;
				}

#if 0
				// Check Z. TODO : Why reset Z from 0.0 to 1.0 ? (Maybe fog-related?)
				if (pVertexDataAsFloat[2] == 0.0f) {
					// LOG_TEST_CASE("D3DFVF_XYZRHW (Z)"); // Test-case : Many XDK Samples (AlphaFog, PointSprites)
					pVertexDataAsFloat[2] = 1.0f;
				}
#endif
#if 1
				// Check RHW. TODO : Why reset from 0.0 to 1.0 ? (Maybe 1.0 indicates that the vertices are not to be transformed)
				if (pVertexDataAsFloat[3] == 0.0f) {
					// LOG_TEST_CASE("D3DFVF_XYZRHW (RHW)"); // Test-case : Many XDK Samples (AlphaFog, PointSprites)
					pVertexDataAsFloat[3] = 1.0f;
				}
#endif
			}

			// Normalize texture coordinates in FVF stream if needed
			if (m_uiTextureCoordinatesByteOffsetInVertex > 0) { // implies bNeedTextureNormalization (using one is more efficient than both)
				FLOAT *pVertexUVData = (FLOAT*)((uintptr_t)pVertexDataAsFloat + m_uiTextureCoordinatesByteOffsetInVertex);
				for (unsigned int i = 0; i < XTL::X_D3DTS_STAGECOUNT; i++) {
					if (m_FVFInfo_ActivePixelContainer[i].bTexIsLinear) {
						switch (m_FVFInfo_ActivePixelContainer[i].NrTexCoords) {
						case 0:
							LOG_TEST_CASE("Normalize 0D?");
							break;
						case 1:
							LOG_TEST_CASE("Normalize 1D");
							pVertexUVData[0] /= m_FVFInfo_ActivePixelContainer[i].Width;
							break;
						case 2:
							pVertexUVData[0] /= m_FVFInfo_ActivePixelContainer[i].Width;
							pVertexUVData[1] /= m_FVFInfo_ActivePixelContainer[i].Height;
							break;
						case 3:
							LOG_TEST_CASE("Normalize 3D");
							// Test case : HeatShimmer
							pVertexUVData[0] /= m_FVFInfo_ActivePixelContainer[i].Width;
							pVertexUVData[1] /= m_FVFInfo_ActivePixelContainer[i].Height;
							pVertexUVData[2] /= m_FVFInfo_ActivePixelContainer[i].Depth;
							break;
						default:
							LOG_TEST_CASE("Normalize ?D");
							break;
						}
					}

					pVertexUVData += m_FVFInfo_ActivePixelContainer[i].NrTexCoords;
				}
			}
		}
	}

    patchedStream.m_bIsCached = true;
    patchedStream.m_uiCachedVertexStreamInformationHash = pVertexShaderSteamInfoHash;
    patchedStream.m_pCachedXboxVertexData = pXboxVertexData;
    patchedStream.m_uiCachedXboxVertexDataSize = XboxVertexDataSize;
    patchedStream.m_uiCachedXboxVertexDataHash = XboxVertexDataHash;
    patchedStream.m_uiCachedXboxVertexStride = uiXboxVertexStride;
    patchedStream.m_uiCachedHostVertexStride = uiHostVertexStride;
	patchedStream.m_bCachedHostVertexStreamZeroDataIsAllocated = bAllocateHostVertexStreamZero;
    if (bAllocateHostVertexStreamZero) {
        patchedStream.m_pCachedHostVertexStreamZeroData = pHostVertexData;
    } else {
        // assert(pHostVertexBuffer != nullptr);
        pHostVertexBuffer->Unlock();
        patchedStream.m_pCachedHostVertexBuffer = pHostVertexBuffer;
    }

	return patchedStream;
}

void CxbxVertexBufferConverter::Apply(CxbxDrawContext *pDrawContext)
{
	if ((pDrawContext->XboxPrimitiveType < XTL::X_D3DPT_POINTLIST) || (pDrawContext->XboxPrimitiveType > XTL::X_D3DPT_POLYGON))
		CxbxKrnlCleanup("Unknown primitive type: 0x%.02X\n", pDrawContext->XboxPrimitiveType);

	// If we are drawing from an offset, we know that the vertex count must have
	// 'offset' vertices before the first drawn vertices
	pDrawContext->VerticesInBuffer = pDrawContext->dwStartVertex + pDrawContext->dwVertexCount;
	// When this is an indexed draw, take the index buffer into account
	if (pDrawContext->pXboxIndexData) {
		// Is the highest index in this buffer not set yet?
		if (pDrawContext->HighIndex == 0) {
			// TODO : Instead of calling WalkIndexBuffer here, set LowIndex and HighIndex
			// in all callers that end up here (since they might be able to avoid the call)
			LOG_TEST_CASE("HighIndex == 0"); // TODO : If this is never hit, replace entire block by assert(pDrawContext->HighIndex > 0);
			WalkIndexBuffer(pDrawContext->LowIndex, pDrawContext->HighIndex, pDrawContext->pXboxIndexData, pDrawContext->dwVertexCount);
		}
		// Convert highest index (including the base offset) into a count
		DWORD dwHighestVertexCount = pDrawContext->dwBaseVertexIndex + pDrawContext->HighIndex + 1;
		// Use the biggest vertex count that can be reached
		if (pDrawContext->VerticesInBuffer < dwHighestVertexCount)
			pDrawContext->VerticesInBuffer = dwHighestVertexCount;
	}

	// TODO : What if we were to check if(pDrawContext->pXboxVertexStreamZeroData != xbnullptr) and in that case
	// call (un)patched D3DDevice_SetStateUP(), otherwise call (un)patched D3DDevice_SetStateVB(g_Xbox_BaseVertexIndex); ?

	CxbxUpdateNativeD3DResources();

	CxbxUpdateActiveVertexShader(pDrawContext->VerticesInBuffer); // Note : Not called from CxbxUpdateNativeD3DResources (lacks access to pDrawContext)

    m_pVertexShaderInfo = nullptr;
    if (VshHandleIsVertexShader(g_Xbox_VertexShader_Handle)) {
        m_pVertexShaderInfo = &(GetCxbxVertexShader(g_Xbox_VertexShader_Handle)->VertexShaderInfo);
    }

	PrepareStreamConversion();

	// Is this a user memory pointer based draws?
	if (pDrawContext->pXboxVertexStreamZeroData != xbnullptr) {
		// For these type of draws, only stream zero needs to be processed
		CxbxPatchedStream& PatchedStream = ConvertStream(pDrawContext, 0);
		// Instead of calling SetStreamSource on host, we pass the patched stream and stride
		// directly to host's DrawPrimitiveUP/DrawIndexedPrimitiveUP via the drawing context :
		pDrawContext->pHostVertexStreamZeroData = PatchedStream.m_pCachedHostVertexStreamZeroData;
		pDrawContext->uiHostVertexStreamZeroStride = PatchedStream.m_uiCachedHostVertexStride;
	} else {
		// Convert and activate all 16 Xbox streams
		for(UINT uiStream = 0; uiStream < X_VSH_MAX_STREAMS; uiStream++) {
			CxbxPatchedStream &PatchedStream = ConvertStream(pDrawContext, uiStream);
			PatchedStream.Activate(pDrawContext, uiStream);
		}
	}

	if (pDrawContext->XboxPrimitiveType == XTL::X_D3DPT_QUADSTRIP) {
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
		pDrawContext->dwHostPrimitiveCount = ConvertXboxVertexCountToPrimitiveCount(XTL::X_D3DPT_TRIANGLESTRIP, pDrawContext->dwVertexCount);
	} else {
		pDrawContext->dwHostPrimitiveCount = ConvertXboxVertexCountToPrimitiveCount(pDrawContext->XboxPrimitiveType, pDrawContext->dwVertexCount);
	}

	if (pDrawContext->XboxPrimitiveType == XTL::X_D3DPT_POLYGON) {
		// Convex polygon is the same as a triangle fan.
		// No need to set : pDrawContext->XboxPrimitiveType = X_D3DPT_TRIANGLEFAN;
		// Test-case : Panzer Dragoon ORTA (when entering in-game)
		LOG_TEST_CASE("X_D3DPT_POLYGON");
	}
}

VOID EmuFlushIVB()
{
    // Parse IVB table with current FVF shader if possible.
    bool bFVF = VshHandleIsFVF(g_Xbox_VertexShader_Handle);
    DWORD dwCurXboxFVF = (bFVF) ? g_Xbox_VertexShader_Handle : g_InlineVertexBuffer_FVF;

    EmuLog(LOG_LEVEL::DEBUG, "g_InlineVertexBuffer_TableOffset := %d", g_InlineVertexBuffer_TableOffset);

	// Check the given FVF
	switch (dwCurXboxFVF & X_D3DFVF_POSITION_MASK) {
	case 0: // No position ?
		if (bFVF) {
			EmuLog(LOG_LEVEL::WARNING, "EmuFlushIVB(): g_Xbox_VertexShader_Handle isn't a valid FVF - using X_D3DFVF_XYZRHW instead!");
			dwCurXboxFVF |= X_D3DFVF_XYZRHW;
		}
		else {
			EmuLog(LOG_LEVEL::WARNING, "EmuFlushIVB(): using g_InlineVertexBuffer_FVF instead of current FVF!");
			dwCurXboxFVF = g_InlineVertexBuffer_FVF;
		}
		break;
	case X_D3DFVF_XYZRHW:
		// X_D3DFVF_NORMAL isn't allowed in combination with X_D3DFVF_XYZRHW 
		if (dwCurXboxFVF & X_D3DFVF_NORMAL) {
			EmuLog(LOG_LEVEL::WARNING, "EmuFlushIVB(): Normal encountered while X_D3DFVF_XYZRHW is given - switching back to X_D3DFVF_XYZ!");
			dwCurXboxFVF &= ~X_D3DFVF_POSITION_MASK;
			dwCurXboxFVF |= X_D3DFVF_XYZ;
		}
		break;
	}

	DWORD dwPos = dwCurXboxFVF & X_D3DFVF_POSITION_MASK;
	DWORD dwTexN = (dwCurXboxFVF & X_D3DFVF_TEXCOUNT_MASK) >> X_D3DFVF_TEXCOUNT_SHIFT;
	size_t TexSize[XTL::X_D3DTS_STAGECOUNT]; // Xbox supports up to 4 textures

	for (unsigned int i = 0; i < dwTexN; i++) {
		TexSize[i] = XboxFVF_GetNumberOfTextureCoordinates(dwCurXboxFVF, i);
	}

	// Use a tooling function to determine the vertex stride :
	UINT uiStride = XboxFVFToVertexSizeInBytes(dwCurXboxFVF, /*bIncludeTextures=*/true);
	// Make sure the output buffer is big enough 
	UINT NeededSize = g_InlineVertexBuffer_TableOffset * uiStride;
	if (g_InlineVertexBuffer_DataSize < NeededSize) {
		g_InlineVertexBuffer_DataSize = NeededSize;
		if (g_InlineVertexBuffer_pData != nullptr) {
			free(g_InlineVertexBuffer_pData);
		}

		g_InlineVertexBuffer_pData = (FLOAT*)malloc(g_InlineVertexBuffer_DataSize);
	}

	FLOAT *pVertexBufferData = g_InlineVertexBuffer_pData;
	for(unsigned int v=0;v<g_InlineVertexBuffer_TableOffset;v++) {
        *pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Position.x;
        *pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Position.y;
        *pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Position.z;
		if (dwPos == X_D3DFVF_XYZRHW) {
            *pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Rhw;
            EmuLog(LOG_LEVEL::DEBUG, "IVB Position := {%f, %f, %f, %f}", g_InlineVertexBuffer_Table[v].Position.x, g_InlineVertexBuffer_Table[v].Position.y, g_InlineVertexBuffer_Table[v].Position.z, g_InlineVertexBuffer_Table[v].Rhw);
		}
		else { // XYZRHW cannot be combined with NORMAL, but the other XYZ formats can :
			switch (dwPos) {
			case X_D3DFVF_XYZ:
				EmuLog(LOG_LEVEL::DEBUG, "IVB Position := {%f, %f, %f}", g_InlineVertexBuffer_Table[v].Position.x, g_InlineVertexBuffer_Table[v].Position.y, g_InlineVertexBuffer_Table[v].Position.z);
				break;
			case X_D3DFVF_XYZB1:
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[0];
				EmuLog(LOG_LEVEL::DEBUG, "IVB Position := {%f, %f, %f, %f}", g_InlineVertexBuffer_Table[v].Position.x, g_InlineVertexBuffer_Table[v].Position.y, g_InlineVertexBuffer_Table[v].Position.z, g_InlineVertexBuffer_Table[v].Blend[0]);
				break;
			case X_D3DFVF_XYZB2:
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[0];
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[1];
				EmuLog(LOG_LEVEL::DEBUG, "IVB Position := {%f, %f, %f, %f, %f}", g_InlineVertexBuffer_Table[v].Position.x, g_InlineVertexBuffer_Table[v].Position.y, g_InlineVertexBuffer_Table[v].Position.z, g_InlineVertexBuffer_Table[v].Blend[0], g_InlineVertexBuffer_Table[v].Blend[1]);
				break;
			case X_D3DFVF_XYZB3:
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[0];
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[1];
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[2];
				EmuLog(LOG_LEVEL::DEBUG, "IVB Position := {%f, %f, %f, %f, %f, %f}", g_InlineVertexBuffer_Table[v].Position.x, g_InlineVertexBuffer_Table[v].Position.y, g_InlineVertexBuffer_Table[v].Position.z, g_InlineVertexBuffer_Table[v].Blend[0], g_InlineVertexBuffer_Table[v].Blend[1], g_InlineVertexBuffer_Table[v].Blend[2]);
				break;
			case X_D3DFVF_XYZB4:
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[0];
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[1];
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[2];
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Blend[3];
				EmuLog(LOG_LEVEL::DEBUG, "IVB Position := {%f, %f, %f, %f, %f, %f, %f}", g_InlineVertexBuffer_Table[v].Position.x, g_InlineVertexBuffer_Table[v].Position.y, g_InlineVertexBuffer_Table[v].Position.z, g_InlineVertexBuffer_Table[v].Blend[0], g_InlineVertexBuffer_Table[v].Blend[1], g_InlineVertexBuffer_Table[v].Blend[2], g_InlineVertexBuffer_Table[v].Blend[3]);
				break;
			default:
				CxbxKrnlCleanup("Unsupported Position Mask (FVF := 0x%.08X dwPos := 0x%.08X)", dwCurXboxFVF, dwPos);
				break;
			}

			if (dwCurXboxFVF & X_D3DFVF_NORMAL) {
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Normal.x;
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Normal.y;
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Normal.z;
				EmuLog(LOG_LEVEL::DEBUG, "IVB Normal := {%f, %f, %f}", g_InlineVertexBuffer_Table[v].Normal.x, g_InlineVertexBuffer_Table[v].Normal.y, g_InlineVertexBuffer_Table[v].Normal.z);
			}
		}

#if 0 // TODO : Was this supported on Xbox from some point in time (pun intended)?
		if (dwCurXboxFVF & X_D3DFVF_PSIZE) {
			*(DWORD*)pVertexBufferData++ = g_InlineVertexBuffer_Table[v].PointSize;
			EmuLog(LOG_LEVEL::DEBUG, "IVB PointSize := 0x%.08X", g_InlineVertexBuffer_Table[v].PointSize);
		}
#endif

        if (dwCurXboxFVF & X_D3DFVF_DIFFUSE) {
            *(DWORD*)pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Diffuse;
            EmuLog(LOG_LEVEL::DEBUG, "IVB Diffuse := 0x%.08X", g_InlineVertexBuffer_Table[v].Diffuse);
        }

		if (dwCurXboxFVF & X_D3DFVF_SPECULAR) {
			*(DWORD*)pVertexBufferData++ = g_InlineVertexBuffer_Table[v].Specular;
			EmuLog(LOG_LEVEL::DEBUG, "IVB Specular := 0x%.08X", g_InlineVertexBuffer_Table[v].Specular);
		}

		for (unsigned int i = 0; i < dwTexN; i++) {
            *pVertexBufferData++ = g_InlineVertexBuffer_Table[v].TexCoord[i].x;
			if (TexSize[i] >= 2) {
				*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].TexCoord[i].y;
				if (TexSize[i] >= 3) {
					*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].TexCoord[i].z;
					if (TexSize[i] >= 4) {
						*pVertexBufferData++ = g_InlineVertexBuffer_Table[v].TexCoord[i].w;
					}
				}
			}

			if (g_bPrintfOn) {
				switch (TexSize[i]) {
				case 1: EmuLog(LOG_LEVEL::DEBUG, "IVB TexCoord%d := {%f}", i + 1, g_InlineVertexBuffer_Table[v].TexCoord[i].x); break;
				case 2: EmuLog(LOG_LEVEL::DEBUG, "IVB TexCoord%d := {%f, %f}", i + 1, g_InlineVertexBuffer_Table[v].TexCoord[i].x, g_InlineVertexBuffer_Table[v].TexCoord[i].y); break;
				case 3: EmuLog(LOG_LEVEL::DEBUG, "IVB TexCoord%d := {%f, %f, %f}", i + 1, g_InlineVertexBuffer_Table[v].TexCoord[i].x, g_InlineVertexBuffer_Table[v].TexCoord[i].y, g_InlineVertexBuffer_Table[v].TexCoord[i].z); break;
				case 4: EmuLog(LOG_LEVEL::DEBUG, "IVB TexCoord%d := {%f, %f, %f, %f}", i + 1, g_InlineVertexBuffer_Table[v].TexCoord[i].x, g_InlineVertexBuffer_Table[v].TexCoord[i].y, g_InlineVertexBuffer_Table[v].TexCoord[i].z, g_InlineVertexBuffer_Table[v].TexCoord[i].w); break;
				}
			}
        }

		if (v == 0) {
			unsigned int VertexBufferUsage = (uintptr_t)pVertexBufferData - (uintptr_t)g_InlineVertexBuffer_pData;
			if (VertexBufferUsage != uiStride) {
				CxbxKrnlCleanup("EmuFlushIVB uses wrong stride!");
			}
		}
	}

	CxbxDrawContext DrawContext = {};

	DrawContext.XboxPrimitiveType = g_InlineVertexBuffer_PrimitiveType;
	DrawContext.dwVertexCount = g_InlineVertexBuffer_TableOffset;
	DrawContext.pXboxVertexStreamZeroData = g_InlineVertexBuffer_pData;
	DrawContext.uiXboxVertexStreamZeroStride = uiStride;

	HRESULT hRet;

	if (bFVF) {
		g_pD3DDevice->SetVertexShader(nullptr);
		hRet = g_pD3DDevice->SetFVF(dwCurXboxFVF); // Note : No comversion needed since host FVF are compatible with Xbox FVF
		//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShader");
	}

	CxbxDrawPrimitiveUP(DrawContext);
	if (bFVF) {
		hRet = g_pD3DDevice->SetFVF(g_Xbox_VertexShader_Handle);
		//DEBUG_D3DRESULT(hRet, "g_pD3DDevice->SetVertexShader");
	}
    g_InlineVertexBuffer_TableOffset = 0; // Might not be needed (also cleared in D3DDevice_Begin)
}

void CxbxImpl_SetStreamSource(UINT StreamNumber, XTL::X_D3DVertexBuffer* pStreamData, UINT Stride)
{
	if (pStreamData != xbnullptr && Stride == 0) {
		LOG_TEST_CASE("CxbxImpl_SetStreamSource : Stream assigned, and stride set to 0 (might be okay)");
	}

	assert(StreamNumber < X_VSH_MAX_STREAMS);

	g_Xbox_SetStreamSource[StreamNumber].VertexBuffer = pStreamData;
	g_Xbox_SetStreamSource[StreamNumber].Stride = Stride;
}
