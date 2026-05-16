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
// *  All rights reserved
// *
// ******************************************************************
#include "../EmuD3D8_common.h"
#include "Backend_D3D11_PageTracker.h"
#include <vector>


xbox::X_D3DRESOURCETYPE GetXboxD3DResourceType(const xbox::X_D3DResource *pXboxResource)
{
	DWORD Type = GetXboxCommonResourceType(pXboxResource);
	switch (Type)
	{
	case X_D3DCOMMON_TYPE_VERTEXBUFFER:
		return xbox::X_D3DRTYPE_VERTEXBUFFER;
	case X_D3DCOMMON_TYPE_INDEXBUFFER:
		return xbox::X_D3DRTYPE_INDEXBUFFER;
	case X_D3DCOMMON_TYPE_PUSHBUFFER:
		return xbox::X_D3DRTYPE_PUSHBUFFER;
	case X_D3DCOMMON_TYPE_PALETTE:
		return xbox::X_D3DRTYPE_PALETTE;
	case X_D3DCOMMON_TYPE_TEXTURE:
	{
		DWORD Format = ((xbox::X_D3DPixelContainer *)pXboxResource)->Format;
		if (Format & X_D3DFORMAT_CUBEMAP)
			return xbox::X_D3DRTYPE_CUBETEXTURE;

		if (GetXboxPixelContainerDimensionCount((xbox::X_D3DPixelContainer *)pXboxResource) > 2)
			return xbox::X_D3DRTYPE_VOLUMETEXTURE;

		return xbox::X_D3DRTYPE_TEXTURE;
	}
	case X_D3DCOMMON_TYPE_SURFACE:
	{
		if (GetXboxPixelContainerDimensionCount((xbox::X_D3DPixelContainer *)pXboxResource) > 2)
			return xbox::X_D3DRTYPE_VOLUME;

		return xbox::X_D3DRTYPE_SURFACE;
	}
	case X_D3DCOMMON_TYPE_FIXUP:
		return xbox::X_D3DRTYPE_FIXUP;
	}

	return xbox::X_D3DRTYPE_NONE;
}

// This can be used to determine if resource Data adddresses
// need the PHYSICAL_MAP_BASE  bit set or cleared
inline bool IsResourceTypeGPUReadable(const DWORD ResourceType)
{
	switch (ResourceType) {
	case X_D3DCOMMON_TYPE_VERTEXBUFFER:
		return true;
	case X_D3DCOMMON_TYPE_INDEXBUFFER:
		/// assert(false); // Index buffers are not allowed to be registered
		break;
	case X_D3DCOMMON_TYPE_PUSHBUFFER:
		return false;
	case X_D3DCOMMON_TYPE_PALETTE:
		return true;
	case X_D3DCOMMON_TYPE_TEXTURE:
		return true;
	case X_D3DCOMMON_TYPE_SURFACE:
		return true;
	case X_D3DCOMMON_TYPE_FIXUP:
		// assert(false); // Fixup's are not allowed to be registered
		break;
	default:
		CxbxrAbort("Unhandled resource type");
	}

	return false;
}

inline bool IsPaletizedTexture(const xbox::dword_xt XboxPixelContainer_Format)
{
	return GetXboxPixelContainerFormat(XboxPixelContainer_Format) == xbox::X_D3DFMT_P8;
}

void *GetDataFromXboxResource(xbox::X_D3DResource *pXboxResource)
{
	// Don't pass in unassigned Xbox resources
	if (pXboxResource == xbox::zeroptr)
		return nullptr;

	xbox::addr_xt pData = pXboxResource->Data;
	if (pData == xbox::zero)
		return nullptr;

	DWORD dwCommonType = GetXboxCommonResourceType(pXboxResource);
	if (IsResourceTypeGPUReadable(dwCommonType))
		pData |= PHYSICAL_MAP_BASE;

	return (uint8_t*)pData;
}

// D3DUSAGE_INVALID is now constexpr in RenderGlobals.h

bool IsResourceAPixelContainer(xbox::dword_xt XboxResource_Common)
{
	DWORD Type = GetXboxCommonResourceType(XboxResource_Common);
	switch (Type)
	{
	case X_D3DCOMMON_TYPE_TEXTURE:
	case X_D3DCOMMON_TYPE_SURFACE:
		return true;
	}

	return false;
}

bool IsResourceAPixelContainer(xbox::X_D3DResource* pXboxResource)
{
	// Don't pass in unassigned Xbox resources
	assert(pXboxResource != xbox::zeroptr);

	return IsResourceAPixelContainer(pXboxResource->Common);
}

// resource_info_t and resource_cache_t are defined in Backend_D3D11.h
resource_cache_t g_Cxbx_Cached_Direct3DResources;
resource_cache_t g_Cxbx_Cached_PaletizedTextures;

// Monotonic counter for LRU eviction (stamped on each cache access)
static uint32_t g_ResourceCacheAccessCounter = 0;

resource_cache_t& GetResourceCache(resource_key_t& key)
{
	return IsResourceAPixelContainer(key.Common) && IsPaletizedTexture(key.Format)
		? g_Cxbx_Cached_PaletizedTextures : g_Cxbx_Cached_Direct3DResources;
}

resource_key_t GetHostResourceKey(xbox::X_D3DResource* pXboxResource, int iTextureStage)
{
	resource_key_t key = {};
	if (pXboxResource != xbox::zeroptr) {
		// Initially, don't base the key on the address of the resource, but on it's uniquely identifying values
		key.Data = pXboxResource->Data;
		key.Common = pXboxResource->Common & CXBX_D3DCOMMON_IDENTIFYING_MASK;
		if (IsResourceAPixelContainer(pXboxResource)) {
			// Pixel containers have more values they must be identified by:
			auto pPixelContainer = (xbox::X_D3DPixelContainer*)pXboxResource;
			key.Format = pPixelContainer->Format;
			key.Size = pPixelContainer->Size;
			// For paletized textures, include the current palette hash as well
			if (IsPaletizedTexture(pPixelContainer->Format)) {
				if (iTextureStage < 0) {
					LOG_TEST_CASE("Unknown texture stage!");
				} else {
					assert(iTextureStage < xbox::X_D3DTS_STAGECOUNT);
					auto paletteState = NV2AGetPaletteState(iTextureStage);
					if (paletteState.data) {
						key.PaletteHash = ComputeHash(paletteState.data, paletteState.size);
					}
				}
			}
		} else {
			// For other resource types, do include their Xbox resource address (TODO : come up with something better)
			key.ResourceAddr = (xbox::addr_xt)pXboxResource;
		}
	}

	return key;
}

void FreeHostResource(resource_key_t key)
{
	// Release the host resource and remove it from the list
	auto& ResourceCache = GetResourceCache(key);
	ResourceCache.erase(key);
}

void ClearResourceCache(resource_cache_t& ResourceCache)
{
	ResourceCache.clear();
}

void ClearAllResourceCaches()
{
	ClearResourceCache(g_Cxbx_Cached_PaletizedTextures);
	ClearResourceCache(g_Cxbx_Cached_Direct3DResources);
}

void PrunePaletizedTexturesCache()
{
	constexpr size_t CACHE_HIGH_WATERMARK = 1500;
	constexpr size_t CACHE_LOW_WATERMARK  = 1000; // Evict down to this level

	if (g_Cxbx_Cached_PaletizedTextures.size() < CACHE_HIGH_WATERMARK)
		return;

	// LRU eviction: find the median access counter and evict entries below it
	// to bring the cache back to the low watermark
	size_t evictCount = g_Cxbx_Cached_PaletizedTextures.size() - CACHE_LOW_WATERMARK;

	// Collect access timestamps and find a threshold to evict the oldest entries
	std::vector<uint32_t> accessTimes;
	accessTimes.reserve(g_Cxbx_Cached_PaletizedTextures.size());
	for (auto& entry : g_Cxbx_Cached_PaletizedTextures) {
		accessTimes.push_back(entry.second.lastAccessFrame);
	}
	std::nth_element(accessTimes.begin(), accessTimes.begin() + evictCount, accessTimes.end());
	uint32_t threshold = accessTimes[evictCount]; // Evict everything at or below this value

	for (auto it = g_Cxbx_Cached_PaletizedTextures.begin(); it != g_Cxbx_Cached_PaletizedTextures.end(); ) {
		if (it->second.lastAccessFrame <= threshold) {
			it = g_Cxbx_Cached_PaletizedTextures.erase(it);
		} else {
			++it;
		}
	}
}

void PruneResourceCache()
{
	constexpr size_t CACHE_HIGH_WATERMARK = 2048;
	constexpr size_t CACHE_LOW_WATERMARK  = 1024;

	if (g_Cxbx_Cached_Direct3DResources.size() < CACHE_HIGH_WATERMARK)
		return;

	size_t evictCount = g_Cxbx_Cached_Direct3DResources.size() - CACHE_LOW_WATERMARK;

	std::vector<uint32_t> accessTimes;
	accessTimes.reserve(g_Cxbx_Cached_Direct3DResources.size());
	for (auto& entry : g_Cxbx_Cached_Direct3DResources) {
		accessTimes.push_back(entry.second.lastAccessFrame);
	}
	std::nth_element(accessTimes.begin(), accessTimes.begin() + evictCount, accessTimes.end());
	uint32_t threshold = accessTimes[evictCount];

	for (auto it = g_Cxbx_Cached_Direct3DResources.begin(); it != g_Cxbx_Cached_Direct3DResources.end(); ) {
		if (it->second.lastAccessFrame <= threshold) {
			it = g_Cxbx_Cached_Direct3DResources.erase(it);
		} else {
			++it;
		}
	}
}

// Forward declaration (defined later in this file)
static void EmuVerifyResourceIsRegistered(xbox::X_D3DResource *pResource, DWORD D3DUsage, int iTextureStage, DWORD dwSize);

ID3D11Resource *GetHostResource(xbox::X_D3DResource *pXboxResource, DWORD D3DUsage = 0, int iTextureStage = -1)
{
	if (pXboxResource == xbox::zeroptr || pXboxResource->Data == xbox::zero)
		return nullptr;

	EmuVerifyResourceIsRegistered(pXboxResource, D3DUsage, iTextureStage, /*dwSize=*/0);

	auto key = GetHostResourceKey(pXboxResource, iTextureStage);
	auto& ResourceCache = GetResourceCache(key);
	auto it = ResourceCache.find(key);
	if (it == ResourceCache.end() || !it->second.pHostResource) {
		EmuLog(LOG_LEVEL::WARNING, "GetHostResource: Resource not registered or does not have a host counterpart!");
		return nullptr;
	}

	it->second.lastAccessFrame = ++g_ResourceCacheAccessCounter;
	return it->second.pHostResource.Get();
}

// Forward declaration of CxbxGetPixelContainerMeasures to prevent
// polluting the diff too much by reshuffling functions around
void CxbxGetPixelContainerMeasures
(
	xbox::X_D3DPixelContainer *pPixelContainer,
	DWORD dwMipMapLevel,
	UINT *pWidth,
	UINT *pHeight,
	UINT *pDepth,
	UINT *pRowPitch,
	UINT *pSlicePitch
);

size_t GetXboxResourceSize(xbox::X_D3DResource* pXboxResource)
{
	// TODO: Smart size calculation based around format of resource
	if (IsResourceAPixelContainer(pXboxResource)) {
		unsigned int Width, Height, Depth, RowPitch, SlicePitch;
		// TODO : Accumulate all mipmap levels!!!
		CxbxGetPixelContainerMeasures(
			(xbox::X_D3DPixelContainer*)pXboxResource,
			0, // dwMipMapLevel
			&Width,
			&Height,
			&Depth,
			&RowPitch,
			&SlicePitch
		);

		return SlicePitch * Depth;
	} else {
		// Fallback to querying the allocation size, if no other calculation was present
		return xbox::MmQueryAllocationSize(GetDataFromXboxResource(pXboxResource));
	}
	
}

bool HostResourceRequiresUpdate(resource_key_t key, xbox::X_D3DResource* pXboxResource, DWORD dwSize)
{
	if (pXboxResource == nullptr) {
		return false;
	}

	// Currently, we only dynamically update Textures and Surfaces, so if our resource
	// isn't of these types, do nothing
	if (!IsResourceAPixelContainer(pXboxResource)) {
		return false;
	}

	auto& ResourceCache = GetResourceCache(key);
	auto it = ResourceCache.find(key);
	if (it == ResourceCache.end()) {
		return false;
	}

	// If the resource size got bigger, we need to re-create it
	if (dwSize > it->second.szXboxDataSize) {
		return true;
	}

	// If the resource type changed, we need to re-create it
	if (it->second.dwXboxResourceType != GetXboxCommonResourceType(pXboxResource)) {
		return true;
	}

	// Render targets have their content managed by the GPU, not by CPU writes
	// to Xbox memory. Don't recreate them based on dirty page tracking, as that
	// would wipe GPU-rendered content (e.g. dynamically rendered cubemap faces).
	if (it->second.HostUsage & D3DUSAGE_RENDERTARGET) {
		return false;
	}

	// Dirty-page-gated texture update: for textures in contiguous memory
	// (0x80000000..0x83FFFFFF), check the page tracker's texture-dirty bitmap
	// instead of hashing. If no pages covering this texture have been written
	// by the CPU since the last host texture upload, the host texture is still valid.
	{
		uintptr_t dataAddr = (uintptr_t)it->second.pXboxData;
		if (dataAddr >= CONTIGUOUS_MEMORY_BASE &&
			dataAddr < (CONTIGUOUS_MEMORY_BASE + XBOX_CONTIGUOUS_MEMORY_SIZE))
		{
			uint32_t offset = (uint32_t)(dataAddr - CONTIGUOUS_MEMORY_BASE);
			uint32_t size = (uint32_t)it->second.szXboxDataSize;
			return CxbxPageTrackerIsTextureDirty(offset, size);
		}
	}

	// Texture not in contiguous memory — can't use page tracking, assume dirty
	return true;
}

void SetHostResource(xbox::X_D3DResource* pXboxResource, ID3D11Resource* pHostResource, int iTextureStage, DWORD D3DUsage, DXGI_FORMAT PCFormat)
{
	auto key = GetHostResourceKey(pXboxResource, iTextureStage);
	auto& ResourceCache = GetResourceCache(key);
	auto& resourceInfo = ResourceCache[key];	// Implicitely inserts a new entry if not already existing

	if (resourceInfo.pHostResource) {
		EmuLog(LOG_LEVEL::WARNING, "SetHostResource: Overwriting an existing host resource");
	}

	// Increments reference count
	resourceInfo.pHostResource = pHostResource;
	resourceInfo.dwXboxResourceType = GetXboxCommonResourceType(pXboxResource);
	resourceInfo.pXboxData = GetDataFromXboxResource(pXboxResource);
	resourceInfo.szXboxDataSize = GetXboxResourceSize(pXboxResource);
	if (PCFormat == EMUFMT_UNKNOWN) {
		D3D11_TEXTURE2D_DESC tex2dDesc = {};
		D3D11_TEXTURE3D_DESC tex3dDesc = {};
		switch (resourceInfo.dwXboxResourceType) {// TODO : Better check pHostResource class type
		case xbox::X_D3DRTYPE_SURFACE:
		case xbox::X_D3DRTYPE_TEXTURE:
		case xbox::X_D3DRTYPE_CUBETEXTURE:
			((ID3D11Texture2D*)pHostResource)->GetDesc(&tex2dDesc);
			PCFormat = tex2dDesc.Format;
			break;
		case xbox::X_D3DRTYPE_VOLUMETEXTURE:
			((ID3D11Texture3D*)pHostResource)->GetDesc(&tex3dDesc);
			PCFormat = tex3dDesc.Format;
			break;
		}
	}

	resourceInfo.HostFormat = PCFormat;
	resourceInfo.HostUsage = D3DUsage;
	resourceInfo.lastAccessFrame = ++g_ResourceCacheAccessCounter;
}

// Inline Set* functions are now in RenderGlobals.h

ID3D11Texture2D *GetHostSurface(xbox::X_D3DResource *pXboxResource, DWORD D3DUsage)
{
	if (pXboxResource == xbox::zeroptr)
		return nullptr;

	if (GetXboxCommonResourceType(pXboxResource) != X_D3DCOMMON_TYPE_SURFACE) // Allows breakpoint below
		assert(GetXboxCommonResourceType(pXboxResource) == X_D3DCOMMON_TYPE_SURFACE);

	return (ID3D11Texture2D*) GetHostResource(pXboxResource, D3DUsage);
}

ID3D11Resource *GetHostBaseTexture(xbox::X_D3DResource *pXboxResource, DWORD D3DUsage, int iTextureStage)
{
	if (pXboxResource == xbox::zeroptr)
		return nullptr;

	if (GetXboxCommonResourceType(pXboxResource) != X_D3DCOMMON_TYPE_TEXTURE) { // Allows breakpoint below
		// test-case : Burnout and Outrun 2006 hit this case (retrieving a surface instead of a texture)
		// TODO : Surfaces can be set in the texture stages, instead of textures - see preparations in CxbxConvertXboxSurfaceToHostTexture
		// We'll need to wrap the surface somehow before using it as a texture
		LOG_TEST_CASE("GetHostBaseTexture called on a non-texture object");
		return nullptr;
		// Note : We'd like to remove the above and do the following instead,
		// but we can't yet since that seems to cause a "CreateCubeTexture Failed!"
		// regression. The root cause for that seems to stem from the X_D3DRTYPE_SURFACE
		// handling in CreateHostResource.
		//assert(GetXboxCommonResourceType(pXboxResource) == X_D3DCOMMON_TYPE_TEXTURE);
	}

	return (ID3D11Resource*)GetHostResource(pXboxResource, D3DUsage, iTextureStage);
}

ID3D11Resource *GetHostBaseTextureWithFormat(xbox::X_D3DResource *pXboxResource, DWORD D3DUsage, int iTextureStage, DXGI_FORMAT *pOutHostFormat)
{
	if (pOutHostFormat) *pOutHostFormat = DXGI_FORMAT_UNKNOWN;

	if (pXboxResource == xbox::zeroptr)
		return nullptr;

	if (GetXboxCommonResourceType(pXboxResource) != X_D3DCOMMON_TYPE_TEXTURE)
		return nullptr;

	if (pXboxResource->Data == xbox::zero)
		return nullptr;

	// Fast path: check cache first without calling EmuVerifyResourceIsRegistered.
	// The verify function does its own lookup, so skipping it when the resource
	// is already cached eliminates one redundant hash map lookup per draw.
	auto key = GetHostResourceKey(pXboxResource, iTextureStage);
	auto& ResourceCache = GetResourceCache(key);
	auto it = ResourceCache.find(key);
	if (it == ResourceCache.end() || !it->second.pHostResource) {
		// Resource not in cache — verify (may create it), then re-lookup
		EmuVerifyResourceIsRegistered(pXboxResource, D3DUsage, iTextureStage, /*dwSize=*/0);
		it = ResourceCache.find(key);
		if (it == ResourceCache.end() || !it->second.pHostResource)
			return nullptr;
	} else if (!(it->second.HostUsage & D3DUSAGE_RENDERTARGET)) {
		// Cache hit — check if CPU has written to this texture since last upload.
		// Render targets are GPU-managed and skip this check.
		uintptr_t dataAddr = (uintptr_t)it->second.pXboxData;
		if (dataAddr >= CONTIGUOUS_MEMORY_BASE &&
			dataAddr < (CONTIGUOUS_MEMORY_BASE + XBOX_CONTIGUOUS_MEMORY_SIZE))
		{
			uint32_t offset = (uint32_t)(dataAddr - CONTIGUOUS_MEMORY_BASE);
			if (CxbxPageTrackerIsTextureDirty(offset, (uint32_t)it->second.szXboxDataSize)) {
				// Texture pages were modified — re-verify (may re-upload)
				EmuVerifyResourceIsRegistered(pXboxResource, D3DUsage, iTextureStage, /*dwSize=*/0);
				it = ResourceCache.find(key);
				if (it == ResourceCache.end() || !it->second.pHostResource)
					return nullptr;
			}
		}
	}

	it->second.lastAccessFrame = ++g_ResourceCacheAccessCounter;
	if (pOutHostFormat)
		*pOutHostFormat = it->second.HostFormat;
	return it->second.pHostResource.Get();
}

ID3D11Texture3D *GetHostVolumeTexture(xbox::X_D3DResource *pXboxResource, int iTextureStage)
{
	return (ID3D11Texture3D *)GetHostBaseTexture(pXboxResource, 0, iTextureStage);

	// TODO : Check for 1 face (and 2 dimensions)?
}

int XboxD3DPaletteSizeToBytes(const xbox::X_D3DPALETTESIZE Size)
{
	return (256 * sizeof(D3DCOLOR)) >> (unsigned)Size;
}

xbox::X_D3DPALETTESIZE GetXboxPaletteSize(const xbox::X_D3DPalette *pPalette)
{
	xbox::X_D3DPALETTESIZE PaletteSize = (xbox::X_D3DPALETTESIZE)
		((pPalette->Common & X_D3DPALETTE_COMMON_PALETTESIZE_MASK) >> X_D3DPALETTE_COMMON_PALETTESIZE_SHIFT);

	return PaletteSize;
}

int GetD3DResourceRefCount(ID3D11Resource *EmuResource)
{
	if (EmuResource != nullptr)
	{
		// Get actual reference count by increasing it using AddRef,
		// and relying on the return value of Release (which is
		// probably more reliable than AddRef)
		EmuResource->AddRef();
		return EmuResource->Release();
	}

	return 0;
}

unsigned int CxbxGetPixelContainerDepth
(
	xbox::X_D3DPixelContainer *pPixelContainer
)
{
	if (pPixelContainer->Size == 0) {
		DWORD l2d = (pPixelContainer->Format & X_D3DFORMAT_PSIZE_MASK) >> X_D3DFORMAT_PSIZE_SHIFT;
		return  1 << l2d;
	}

	return 1;
}

unsigned int CxbxGetPixelContainerMipMapLevels
(
	xbox::X_D3DPixelContainer *pPixelContainer
)
{
	if (pPixelContainer->Size == 0) {
		return (pPixelContainer->Format & X_D3DFORMAT_MIPMAP_MASK) >> X_D3DFORMAT_MIPMAP_SHIFT;
	}

	return 1;
}

uint32_t GetPixelContainerWidth(xbox::X_D3DPixelContainer *pPixelContainer)
{
	DWORD Size = pPixelContainer->Size;
	uint32_t Result;

	if (Size != 0) {
		Result = ((Size & X_D3DSIZE_WIDTH_MASK) /* >> X_D3DSIZE_WIDTH_SHIFT*/) + 1;
	}
	else {
		DWORD l2w = (pPixelContainer->Format & X_D3DFORMAT_USIZE_MASK) >> X_D3DFORMAT_USIZE_SHIFT;

		Result = 1 << l2w;
	}

	return Result;
}

uint32_t GetPixelContainerHeight(xbox::X_D3DPixelContainer *pPixelContainer)
{
	DWORD Size = pPixelContainer->Size;
	uint32_t Result;

	if (Size != 0) {
		Result = ((Size & X_D3DSIZE_HEIGHT_MASK) >> X_D3DSIZE_HEIGHT_SHIFT) + 1;
	}
	else {
		DWORD l2h = (pPixelContainer->Format & X_D3DFORMAT_VSIZE_MASK) >> X_D3DFORMAT_VSIZE_SHIFT;

		Result = 1 << l2h;
	}

	return Result;
}

void GetSurfaceFaceAndLevelWithinTexture(xbox::X_D3DSurface* pSurface, xbox::X_D3DBaseTexture* pTexture, UINT& Level, int& Face)
{
   	auto pSurfaceData = (uintptr_t)GetDataFromXboxResource(pSurface);
   	auto pTextureData = (uintptr_t)GetDataFromXboxResource(pTexture);

   	// Fast path: If the data pointers match, this must be the first surface within the texture
   	if (pSurfaceData == pTextureData) {
   	   	Level = 0;
   	   	Face = 0;
   	   	return;
   	}

   	int numLevels = CxbxGetPixelContainerMipMapLevels(pTexture);
   	int numFaces = pTexture->Format & X_D3DFORMAT_CUBEMAP ? 6 : 1;

   	CxbxGetPixelContainerMipMapLevels(pTexture);

   	// First, we need to fetch the dimensions of both the surface and the texture, for use within our calculations
   	UINT textureWidth, textureHeight, textureDepth, textureRowPitch, textureSlicePitch;
   	CxbxGetPixelContainerMeasures(pTexture, 0, &textureWidth, &textureHeight, &textureDepth, &textureRowPitch, &textureSlicePitch);

   	UINT surfaceWidth, surfaceHeight, surfaceDepth, surfaceRowPitch, surfaceSlicePitch;
   	CxbxGetPixelContainerMeasures(pSurface, 0, &surfaceWidth, &surfaceHeight, &surfaceDepth, &surfaceRowPitch, &surfaceSlicePitch);

   	// Iterate through all faces and levels, until we find a matching pointer
   	bool isCompressed = EmuXBFormatIsCompressed(GetXboxPixelContainerFormat(pTexture));
   	int minSize = (isCompressed) ? 4 : 1;
   	int cubeFaceOffset = 0; int cubeFaceSize = 0;
   	auto pData = pTextureData;

   	for (int face = 0; face < numFaces; face++) {
   	   	int mipWidth = textureWidth;
   	   	int mipHeight = textureHeight;
   	   	int mipDepth = textureDepth;
   	   	int mipRowPitch = textureRowPitch;
   	   	int mipDataOffset = 0;

   	   	for (int level = 0; level < numLevels; level++) {
   	   	   	if (pData + mipDataOffset == pSurfaceData) {
   	   	   	   	Level = level;
   	   	   	   	Face = (int)face;
   	   	   	   	return;
   	   	   	}

   	   	   	// Calculate size of this mipmap level
   	   	   	UINT dwMipSize = mipRowPitch * mipHeight;
   	   	   	if (isCompressed) {
   	   	   	   	dwMipSize /= 4;
   	   	   	}

   	   	   	// If this is the first face, set the cube face size
   	   	   	if (face == 0) {
   	   	   	   	cubeFaceSize = ROUND_UP(textureDepth * dwMipSize, X_D3DTEXTURE_CUBEFACE_ALIGNMENT);
   	   	   	}

   	   	   	// Move to the next mip-map and calculate dimensions for the next iteration
   	   	   	mipDataOffset += dwMipSize;

   	   	   	if (mipWidth > minSize) {
   	   	   	   	mipWidth /= 2;
   	   	   	   	mipRowPitch /= 2;
   	   	   	}

   	   	   	if (mipHeight > minSize) {
   	   	   	   	mipHeight /= 2;
   	   	   	}

   	   	   	if (mipDepth > 1) {
   	   	   	   	mipDepth /= 2;
   	   	   	}
   	   	}

   	   	// Move to the next face
   	   	pData += cubeFaceSize;
   	}

   	LOG_TEST_CASE("Could not find Surface within Texture, falling back to Level = 0, Face = D3DCUBEMAP_FACE_POSITIVE_X");
   	Level = 0;
   	Face = 0;
}

// Wrapper function to allow calling without passing a face
void GetSurfaceFaceAndLevelWithinTexture(xbox::X_D3DSurface* pSurface, xbox::X_D3DBaseTexture* pBaseTexture, UINT& Level)
{
   	int face;
   	GetSurfaceFaceAndLevelWithinTexture(pSurface, pBaseTexture, Level, face);
}

// Direct3D initialization (called before emulation begins)
void CreateHostResource(xbox::X_D3DResource *pResource, DWORD D3DUsage, int iTextureStage, DWORD dwSize); // Forward declartion to prevent restructure of code
static void EmuVerifyResourceIsRegistered(xbox::X_D3DResource *pResource, DWORD D3DUsage, int iTextureStage, DWORD dwSize)
{
	// Skip resources without data
	if (pResource->Data == xbox::zero)
		return;

	auto key = GetHostResourceKey(pResource, iTextureStage);
	auto& ResourceCache = GetResourceCache(key);
	auto it = ResourceCache.find(key);
	if (it != ResourceCache.end()) {
		// TODO : Should this check be (D3DUsage & D3DUSAGE_RENDERTARGET) instead?
		if (D3DUsage == D3DUSAGE_RENDERTARGET && IsResourceAPixelContainer(pResource) && EmuXBFormatIsRenderTarget(GetXboxPixelContainerFormat((xbox::X_D3DPixelContainer*)pResource))) {
   	   	   	// Render targets have special behavior: We can't trash them on guest modification
   	   	   	// this fixes an issue where CubeMaps were broken because the surface Set in GetCubeMapSurface
   	   	   	// would be overwritten by the surface created in SetRenderTarget
   	   	   	// However, if a non-rendertarget surface is used here, we'll need to recreate it as a render target!
   	   	   	auto hostResource = it->second.pHostResource;
			auto xboxSurface = (xbox::X_D3DSurface*)pResource;
			auto xboxTexture = (xbox::X_D3DTexture*)pResource;
			auto xboxResourceType = GetXboxD3DResourceType(pResource);
   	   	   	

   	   	   	// Only continue checking if we were able to get the surface desc, if it failed, we fall-through
   	   	   	// to previous resource management behavior
   	   	   	if (it->second.HostUsage != D3DUSAGE_INVALID) {
   	   	   	   	// If this resource is already created as a render target on the host, simply return
   	   	   	   	if (it->second.HostUsage & D3DUSAGE_RENDERTARGET) {
   	   	   	   	   	return;
   	   	   	   	}

   	   	   	   	// The host resource is not a render target, but the Xbox surface is
   	   	   	   	// We need to re-create it as a render target
   	   	   	   	switch (xboxResourceType) {
   	   	   	   	   	case xbox::X_D3DRTYPE_SURFACE: {
   	   	   	   	   	   	auto xboxSurface = (xbox::X_D3DSurface*)pResource;

   	   	   	   	   	   	// Free the host surface
   	   	   	   	   	   	FreeHostResource(key);

   	   	   	   	   	   	// Free the parent texture, if present
   	   	   	   	   	   	xbox::X_D3DTexture* pParentXboxTexture = (pResource) ? (xbox::X_D3DTexture*)xboxSurface->Parent : xbox::zeroptr;

   	   	   	   	   	   	if (pParentXboxTexture) {
   	   	   	   	   	   	   	// Re-create the texture with D3DUSAGE_RENDERTARGET, this will automatically create any child-surfaces
   	   	   	   	   	   	   	FreeHostResource(GetHostResourceKey(pParentXboxTexture, iTextureStage));
   	   	   	   	   	   	   	CreateHostResource(pParentXboxTexture, D3DUsage, iTextureStage, dwSize);
   	   	   	   	   	   	}

   	   	   	   	   	   	// Re-create the surface with D3DUSAGE_RENDERTARGET
   	   	   	   	   	   	CreateHostResource(pResource, D3DUsage, iTextureStage, dwSize);
   	   	   	   	   	} break;
   	   	   	   	   	case xbox::X_D3DRTYPE_TEXTURE: {
   	   	   	   	   	   	auto xboxTexture = (xbox::X_D3DTexture*)pResource;

   	   	   	   	   	   	// Free the host texture
   	   	   	   	   	   	FreeHostResource(key);

   	   	   	   	   	   	// And re-create the texture with D3DUSAGE_RENDERTARGET
   	   	   	   	   	   	CreateHostResource(pResource, D3DUsage, iTextureStage, dwSize);
   	   	   	   	   	} break;
   	   	   	   	   	case xbox::X_D3DRTYPE_CUBETEXTURE: {
   	   	   	   	   	   	// Free the host cubemap texture
   	   	   	   	   	   	FreeHostResource(key);

   	   	   	   	   	   	// Re-create with D3DUSAGE_RENDERTARGET so faces can be used as render targets
   	   	   	   	   	   	CreateHostResource(pResource, D3DUsage, iTextureStage, dwSize);
   	   	   	   	   	} break;
   	   	   	   	   	default:
   	   	   	   	   	   	LOG_TEST_CASE("Unimplemented rendertarget type");
   	   	   	   	}

   	   	   	   	return;
   	   	   	}
		}

		if (!HostResourceRequiresUpdate(key, pResource, dwSize)) {
			return;
		}

		FreeHostResource(key);
	} else {
		resource_info_t newResource;
		ResourceCache[key] = newResource;
	}

	CreateHostResource(pResource, D3DUsage, iTextureStage, dwSize);

	// After creating/re-creating a texture, clear its texture-dirty bits so that
	// subsequent draws skip re-upload until the CPU modifies those pages again.
	if (IsResourceAPixelContainer(pResource)) {
		uintptr_t dataAddr = (uintptr_t)GetDataFromXboxResource(pResource);
		if (dataAddr >= CONTIGUOUS_MEMORY_BASE &&
			dataAddr < (CONTIGUOUS_MEMORY_BASE + XBOX_CONTIGUOUS_MEMORY_SIZE))
		{
			uint32_t offset = (uint32_t)(dataAddr - CONTIGUOUS_MEMORY_BASE);
			uint32_t texSize = (uint32_t)GetXboxResourceSize(pResource);
			CxbxPageTrackerClearTextureDirty(offset, texSize);
		}
	}
}

