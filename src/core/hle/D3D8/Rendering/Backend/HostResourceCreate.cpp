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
#include "Backend_D3D11.h"


// Forward declaration (defined in HostResourceUpload.cpp)
void UploadPixelContainerMips(
	xbox::X_D3DRESOURCETYPE XboxResourceType,
	const char* ResourceTypeName,
	VAddr VirtualAddr,
	int iTextureStage,
	xbox::X_D3DFORMAT X_Format, DXGI_FORMAT PCFormat, UINT dwBPP,
	UINT xboxWidth, UINT xboxHeight,
	DWORD dwDepth, DWORD dwRowPitch, DWORD dwSlicePitch,
	UINT dwMipMapLevels,
	bool bCubemap, bool bSwizzled, bool bCompressed, bool bConvertTextureFormat,
	ComPtr<ID3D11Resource>& pNewHostResource,
	bool bHostIsDynamic
);

bool IsSupportedFormat(xbox::X_D3DFORMAT X_Format, xbox::X_D3DRESOURCETYPE XboxResourceType, DWORD D3DUsage) {
	// TODO : Nuance the following, because the Direct3D 8 docs states
	// CheckDeviceFormat is needed when D3DUSAGE_RENDERTARGET or
	// D3DUSAGE_DYNAMNIC is specified.
	// Otherwise, lookup resource type and accompanying 'SupportedFormat' array
	bool *pbSupportedFormats = g_bSupportsFormatTexture;

	switch (XboxResourceType) {
		case xbox::X_D3DRTYPE_SURFACE: {
			if (D3DUsage & D3DUSAGE_RENDERTARGET) {
				pbSupportedFormats = g_bSupportsFormatSurfaceRenderTarget;
			} else if (D3DUsage & D3DUSAGE_DEPTHSTENCIL) {
				pbSupportedFormats = g_bSupportsFormatSurfaceDepthStencil;
			} else {
				pbSupportedFormats = g_bSupportsFormatSurface;
			}
			break;
		}
		case xbox::X_D3DRTYPE_VOLUME: {
			pbSupportedFormats = g_bSupportsFormatVolumeTexture;
			break;
		}
		case xbox::X_D3DRTYPE_TEXTURE: {
			if (D3DUsage & D3DUSAGE_RENDERTARGET) {
				pbSupportedFormats = g_bSupportsFormatTextureRenderTarget;
			} else if (D3DUsage & D3DUSAGE_DEPTHSTENCIL) {
				pbSupportedFormats = g_bSupportsFormatTextureDepthStencil;
			} else {
				pbSupportedFormats = g_bSupportsFormatTexture;
			}
			break;
		}
		case xbox::X_D3DRTYPE_VOLUMETEXTURE: {
			pbSupportedFormats = g_bSupportsFormatVolumeTexture;
			break;
		}
		case xbox::X_D3DRTYPE_CUBETEXTURE: {
			pbSupportedFormats = g_bSupportsFormatCubeTexture;
			break;
		}
	} // switch XboxResourceType

	return pbSupportedFormats[X_Format];
}

// ---- Helper: Try to resolve a surface via its parent texture ----
// Returns true if the surface was successfully mapped to a parent host texture
static bool TryResolveParentSurface(
	xbox::X_D3DResource* pResource,
	DWORD D3DUsage,
	int iTextureStage)
{
	LOG_INIT;
	xbox::X_D3DSurface *pXboxSurface = (xbox::X_D3DSurface *)pResource;
	xbox::X_D3DBaseTexture *pParentXboxTexture = (pXboxSurface) ? (xbox::X_D3DBaseTexture*)pXboxSurface->Parent : xbox::zeroptr;

	// Don't init the Parent if the Surface and Surface Parent formats differ
	// Happens in some Outrun 2006 SetRenderTarget calls
	if (!pParentXboxTexture || (pXboxSurface->Format != pParentXboxTexture->Format))
		return false;

	// For surfaces with a parent texture, map these to a host texture first
	// TODO : Investigate how it's possible (and how we could fix) the case when
	// the following call to GetHostBaseTexture would reject non-texture resources,
	// which would seem to trigger a "CreateCubeTexture Failed!" regression.
	// In D3D11, child surfaces are the same texture as the parent - subresource
	// indices (computed from face + mip level) are used when creating views.
	ID3D11Resource *pParentHostBaseTexture = GetHostBaseTexture(pParentXboxTexture, D3DUsage, iTextureStage);
	if (pParentHostBaseTexture) {
		int CubeMapFace = 0;
		UINT SurfaceLevel = 0;
		GetSurfaceFaceAndLevelWithinTexture(pXboxSurface, pParentXboxTexture, SurfaceLevel, CubeMapFace);
		SetHostSurface(pXboxSurface, (ID3D11Texture2D *)pParentHostBaseTexture, iTextureStage);
		EmuLog(LOG_LEVEL::DEBUG, "TryResolveParentSurface : D3D11 mapped child surface to parent texture (Face: %u, Level: %u, pResource: 0x%.08X)",
			CubeMapFace, SurfaceLevel, pResource);
		return true;
	}
	EmuLog(LOG_LEVEL::WARNING, "TryResolveParentSurface : D3D11 failed to get parent host texture - falling through");
	return false;
}

// ---- Helper: Try to resolve a volume via its parent volume texture ----
// Returns true if the volume was successfully mapped to a parent host volume texture
static bool TryResolveParentVolume(
	xbox::X_D3DResource* pResource,
	DWORD D3DUsage,
	int iTextureStage)
{
	LOG_INIT;
	xbox::X_D3DVolume *pXboxVolume = (xbox::X_D3DVolume *)pResource;
	xbox::X_D3DVolumeTexture *pParentXboxVolumeTexture = (pXboxVolume) ? (xbox::X_D3DVolumeTexture *)pXboxVolume->Parent : xbox::zeroptr;
	if (!pParentXboxVolumeTexture)
		return false;

	ID3D11Texture3D *pParentHostVolumeTexture = GetHostVolumeTexture(pParentXboxVolumeTexture, iTextureStage);
	if (pParentHostVolumeTexture) {
		SetHostVolume(pXboxVolume, (ID3D11Texture3D *)pParentHostVolumeTexture, iTextureStage);
		UINT VolumeLevel = 0; // TODO : Derive actual level based on pXboxVolume->Data delta to pParentXboxVolumeTexture->Data
		EmuLog(LOG_LEVEL::DEBUG, "TryResolveParentVolume : D3D11 mapped child volume to parent texture (Level: %u, pResource: 0x%.08X)",
			VolumeLevel, pResource);
		return true;
	}
	EmuLog(LOG_LEVEL::WARNING, "TryResolveParentVolume : D3D11 failed to get parent host volume texture - falling through");
	return false;
}

// ---- Helper: Determine the host format for an Xbox pixel format ----
// Sets bConvertTextureFormat and may clear D3DUSAGE_DEPTHSTENCIL from D3DUsage
static DXGI_FORMAT ResolveHostFormat(
	xbox::X_D3DFORMAT X_Format,
	xbox::X_D3DRESOURCETYPE XboxResourceType,
	DWORD& D3DUsage,
	const char* ResourceTypeName,
	bool& bConvertTextureFormat)
{
	DXGI_FORMAT PCFormat;
	bConvertTextureFormat = false;

	if (EmuXBFormatRequiresConversion(X_Format, /*&*/PCFormat)) {
		bConvertTextureFormat = true;

		// Unset D3DUSAGE_DEPTHSTENCIL: It's not possible for ARGB textures to be depth stencils
		// Fixes CreateTexture error in Virtua Cop 3 (Chihiro)
		D3DUsage &= ~D3DUSAGE_DEPTHSTENCIL;

		// D3D11: For 32-bit formats that differ only in channel order,
		// skip CPU conversion and upload raw bytes as R8G8B8A8_UNORM.
		// The pixel shader applies the correct channel swizzle via TEXFMTFIXUP.
		// Exception: render targets are written by the GPU in correct channel order,
		// so use B8G8R8A8_UNORM (the standard host equivalent) to avoid a
		// spurious TEXFMTFIXUP swizzle when the RT is later sampled as a texture.
		if (X_Format == xbox::X_D3DFMT_B8G8R8A8 || X_Format == xbox::X_D3DFMT_LIN_B8G8R8A8 ||
			X_Format == xbox::X_D3DFMT_R8G8B8A8 || X_Format == xbox::X_D3DFMT_LIN_R8G8B8A8) {
			bConvertTextureFormat = false;
			if (D3DUsage & (D3DUSAGE_RENDERTARGET | D3DUSAGE_DEPTHSTENCIL)) {
				PCFormat = EMUFMT_A8R8G8B8; // = DXGI_FORMAT_B8G8R8A8_UNORM (GPU-rendered data is already correct)
			} else {
				PCFormat = EMUFMT_A8B8G8R8; // = DXGI_FORMAT_R8G8B8A8_UNORM (raw byte upload, PS swizzles via TEXFMTFIXUP)
			}
		}
	}
	else {
		if (IsSupportedFormat(X_Format, XboxResourceType, D3DUsage)) {
			PCFormat = EmuXB2PC_D3DFormat(X_Format);
		}
		else {
			if (D3DUsage & D3DUSAGE_DEPTHSTENCIL) {
				EmuLog(LOG_LEVEL::WARNING, "Xbox %s Format %x will be converted to EMUFMT_D24S8", ResourceTypeName, X_Format);
				PCFormat = EMUFMT_D24S8;
			} else if (EmuXBFormatCanBeConverted(X_Format, /*&*/PCFormat)) {
				EmuLog(LOG_LEVEL::WARNING, "Xbox %s Format %x will be converted to ARGB", ResourceTypeName, X_Format);
				bConvertTextureFormat = true;
			} else {
				/*CxbxrAbort*/EmuLog(LOG_LEVEL::WARNING, "Encountered a completely incompatible %s format!", ResourceTypeName);
				PCFormat = EmuXB2PC_D3DFormat(X_Format);
			}
		}
	}

	return PCFormat;
}

// ---- Helper: Create host Texture2D with appropriate usage/format ----
static HRESULT CreateGpuTexture2D(
	UINT hostWidth, UINT hostHeight, UINT dwMipMapLevels,
	DXGI_FORMAT& PCFormat,
	DWORD D3DUsage,
	bool bSwizzled, bool& bConvertTextureFormat, xbox::X_D3DFORMAT X_Format,
	xbox::X_D3DResource* pResource, const char* ResourceTypeName, int iTextureStage,
	ComPtr<ID3D11Resource>& pNewHostResource, bool& bHostIsDynamic)
{
	LOG_INIT;
	D3D11_TEXTURE2D_DESC desc;
	desc.Width = hostWidth;
	desc.Height = hostHeight;
	desc.MipLevels = dwMipMapLevels;
	desc.ArraySize = 1;
	desc.Format = PCFormat;
	desc.SampleDesc.Count = 1;
	desc.SampleDesc.Quality = 0;
	if ((D3DUsage & D3DUSAGE_DEPTHSTENCIL) || IsDepthFormat(PCFormat)) {
		// Depth textures need typeless format + dual bind flags so they can
		// serve as both depth stencil (DSV) and shader resource (SRV).
		// Games like ShadowBuffer create D16 textures without DEPTHSTENCIL
		// usage, then set them as depth targets via SetRenderTarget.
		desc.Format = GetTypelessDepthFormat(PCFormat);
		desc.BindFlags = D3D11_BIND_DEPTH_STENCIL | D3D11_BIND_SHADER_RESOURCE;
		desc.Usage = D3D11_USAGE_DEFAULT;
		desc.CPUAccessFlags = 0;
	} else if (D3DUsage & D3DUSAGE_RENDERTARGET) {
		desc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
		desc.Usage = D3D11_USAGE_DEFAULT;
		desc.CPUAccessFlags = 0;
	} else if (dwMipMapLevels == 1) {
		// Check if the format supports typed UAV access (needed for GPU unswizzle CS)
		UINT fmtSupport = 0;
		bool bCanUAV = SUCCEEDED(g_pD3DDevice->CheckFormatSupport(PCFormat, &fmtSupport))
			&& (fmtSupport & D3D11_FORMAT_SUPPORT_TYPED_UNORDERED_ACCESS_VIEW);

		if (bSwizzled && !bConvertTextureFormat && bCanUAV) {
			// Swizzled textures use DEFAULT + UAV for GPU compute shader unswizzle
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS;
			desc.CPUAccessFlags = 0;
		} else if (bSwizzled && X_Format == xbox::X_D3DFMT_P8) {
			// P8 textures use DEFAULT + UAV for GPU palette expand CS.
			// Typeless allows R32_UINT UAV reinterpretation; SRV uses R8G8B8A8_UNORM.
			desc.Format = DXGI_FORMAT_R8G8B8A8_TYPELESS;
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS;
			desc.CPUAccessFlags = 0;
		} else if (bConvertTextureFormat && CxbxGetFormatConvertType(X_Format) != 0) {
			// Format-convertible textures use DEFAULT + UAV for GPU CS format conversion.
			// Typeless allows R32_UINT UAV reinterpretation; SRV uses R8G8B8A8_UNORM.
			desc.Format = DXGI_FORMAT_R8G8B8A8_TYPELESS;
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS;
			desc.CPUAccessFlags = 0;
		} else {
			// Block-compressed formats (BC1-BC3 / DXT1-DXT5) do not support
			// D3D11_USAGE_DYNAMIC — use DEFAULT with UpdateSubresource upload.
			bool isBlockCompressed = (PCFormat == EMUFMT_DXT1 || PCFormat == EMUFMT_DXT3 || PCFormat == EMUFMT_DXT5);
			if (isBlockCompressed) {
				desc.Usage = D3D11_USAGE_DEFAULT;
				desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
				desc.CPUAccessFlags = 0;
			} else {
				// D3D11_USAGE_DYNAMIC requires MipLevels == 1
				desc.Usage = D3D11_USAGE_DYNAMIC;
				desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
				desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
				bHostIsDynamic = true;
			}
		}
	} else {
		desc.Usage = D3D11_USAGE_DEFAULT;
		desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
		desc.CPUAccessFlags = 0;
	}
	desc.MiscFlags = 0;

	EmuLog(LOG_LEVEL::DEBUG, "CreateTexture2D: XboxFmt=0x%02X PCFmt=%u (%ux%u mips=%u) "
		"Usage=%u BindFlags=0x%X Swizzled=%d ConvertFmt=%d",
		X_Format, desc.Format, desc.Width, desc.Height, desc.MipLevels,
		desc.Usage, desc.BindFlags, bSwizzled, bConvertTextureFormat);

	HRESULT hRet = g_pD3DDevice->CreateTexture2D(&desc, NULL, reinterpret_cast<ID3D11Texture2D**>(pNewHostResource.ReleaseAndGetAddressOf()));
	DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateTexture2D");

	// If the above failed, we might be able to use an ARGB texture instead
	DXGI_FORMAT TmpPCFormat;
	if ((hRet != S_OK) && (PCFormat != EMUFMT_A8R8G8B8) && EmuXBFormatCanBeConverted(X_Format, TmpPCFormat)) {
		desc.Format = TmpPCFormat;
		hRet = g_pD3DDevice->CreateTexture2D(&desc, NULL, reinterpret_cast<ID3D11Texture2D**>(pNewHostResource.ReleaseAndGetAddressOf()));
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateTexture2D");
		if (hRet == S_OK) {
			// Okay, now this works, make sure the texture gets converted
			bConvertTextureFormat = true;
			PCFormat = TmpPCFormat;
		}
	}

	if (hRet != S_OK) {
		CxbxrAbort("CreateTexture2D Failed!\n\n"
			"Error: 0x%X\nFormat: %d\nDimensions: %dx%d", hRet, PCFormat, hostWidth, hostHeight);
	}
	SetHostResource(pResource, pNewHostResource.Get(), iTextureStage, D3DUsage);
	EmuLog(LOG_LEVEL::DEBUG, "CreateGpuPixelContainerResource : Successfully created %s (0x%.08X, 0x%.08X)",
		ResourceTypeName, pResource, pNewHostResource.Get());
	return hRet;
}

// ---- Helper: Create the host GPU resource for a pixel container ----
// Creates the D3D surface/texture/volume/cube resource based on XboxResourceType
static HRESULT CreateGpuPixelContainerResource(
	xbox::X_D3DRESOURCETYPE XboxResourceType,
	UINT hostWidth, UINT hostHeight, DWORD dwDepth, UINT dwMipMapLevels,
	DXGI_FORMAT& PCFormat,
	DWORD& D3DUsage,
	bool bSwizzled, bool& bConvertTextureFormat, xbox::X_D3DFORMAT X_Format,
	xbox::X_D3DResource* pResource, const char* ResourceTypeName, int iTextureStage,
	ComPtr<ID3D11Resource>& pNewHostResource, bool& bHostIsDynamic
)
{
	LOG_INIT;
	HRESULT hRet = S_OK;

	// Create the surface/volume/(volume/cube/)texture
	switch (XboxResourceType) {
	case xbox::X_D3DRTYPE_SURFACE: {
		D3D11_TEXTURE2D_DESC desc;
		desc.Width = hostWidth;
		desc.Height = hostHeight;
		desc.MipLevels = 1;
		desc.ArraySize = 1;
		desc.Format = PCFormat;
		desc.SampleDesc.Count = 1; // No MSAA for now; enabling requires resolve pass infrastructure
		desc.SampleDesc.Quality = 0;
		desc.Usage = D3D11_USAGE_DYNAMIC;
		if ((D3DUsage & D3DUSAGE_DEPTHSTENCIL) || IsDepthFormat(PCFormat)) {
			desc.Format = GetTypelessDepthFormat(PCFormat);
			desc.BindFlags = D3D11_BIND_DEPTH_STENCIL | D3D11_BIND_SHADER_RESOURCE;
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.CPUAccessFlags = 0;
		} else if (D3DUsage & D3DUSAGE_RENDERTARGET) {
			desc.BindFlags = D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE;
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.CPUAccessFlags = 0;
		} else {
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
			desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
		}
		desc.MiscFlags = 0;

		hRet = g_pD3DDevice->CreateTexture2D(&desc, NULL, reinterpret_cast<ID3D11Texture2D**>(pNewHostResource.ReleaseAndGetAddressOf()));
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateTexture2D");
		// If the fallback failed, show an error and exit execution.
		if (hRet != S_OK) {
			// We cannot safely continue in this state.
			CxbxrAbort("CreateImageSurface Failed!\n\nError: %s\nDesc: %s",
				DXGetErrorString(hRet), DXGetErrorDescription(hRet));
		}

		SetHostResource(pResource, pNewHostResource.Get(), iTextureStage, D3DUsage);
		EmuLog(LOG_LEVEL::DEBUG, "CreateGpuPixelContainerResource : Successfully created %s (0x%.08X, 0x%.08X)",
			ResourceTypeName, pResource, pNewHostResource.Get());
		EmuLog(LOG_LEVEL::DEBUG, "CreateGpuPixelContainerResource : Width : %d, Height : %d, Format : %d",
			hostWidth, hostHeight, PCFormat);
		break;
	}

	case xbox::X_D3DRTYPE_VOLUME: {
		// In D3D11, a standalone volume (no parent VolumeTexture) is backed by
		// a single-depth Texture3D. This is rare but can occur if a game creates
		// a volume resource directly.
		D3D11_TEXTURE3D_DESC desc;
		desc.Width = hostWidth;
		desc.Height = hostHeight;
		desc.Depth = dwDepth > 0 ? dwDepth : 1;
		desc.MipLevels = 1;
		desc.Format = PCFormat;
		desc.Usage = D3D11_USAGE_DYNAMIC;
		desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
		desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
		desc.MiscFlags = 0;

		hRet = g_pD3DDevice->CreateTexture3D(&desc, NULL, reinterpret_cast<ID3D11Texture3D**>(pNewHostResource.ReleaseAndGetAddressOf()));
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateTexture3D (standalone volume)");

		if (hRet != S_OK) {
			CxbxrAbort("CreateTexture3D (standalone volume) Failed!\n\n"
				"Error: 0x%X\nFormat: %d\nDimensions: %dx%dx%d", hRet, PCFormat, hostWidth, hostHeight, desc.Depth);
		}
		SetHostVolume(pResource, (ID3D11Texture3D *)pNewHostResource.Get(), iTextureStage);
		EmuLog(LOG_LEVEL::DEBUG, "CreateGpuPixelContainerResource : Successfully created standalone Volume (0x%.08X, 0x%.08X)",
			pResource, pNewHostResource.Get());
		break;
	}

	case xbox::X_D3DRTYPE_TEXTURE: {
		hRet = CreateGpuTexture2D(
			hostWidth, hostHeight, dwMipMapLevels,
			PCFormat, D3DUsage,
			bSwizzled, bConvertTextureFormat, X_Format,
			pResource, ResourceTypeName, iTextureStage,
			pNewHostResource, bHostIsDynamic);
		break;
	}

	case xbox::X_D3DRTYPE_VOLUMETEXTURE: {
		D3D11_TEXTURE3D_DESC desc;
		desc.Width = hostWidth;
		desc.Height = hostHeight;
		desc.Depth = dwDepth;
		desc.MipLevels = dwMipMapLevels;
		desc.Format = PCFormat;
		if (dwMipMapLevels == 1) {
			desc.Usage = D3D11_USAGE_DYNAMIC;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
			desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
			bHostIsDynamic = true;
		} else {
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
			desc.CPUAccessFlags = 0;
		}
		desc.MiscFlags = 0;

		hRet = g_pD3DDevice->CreateTexture3D(&desc, NULL, reinterpret_cast<ID3D11Texture3D**>(pNewHostResource.ReleaseAndGetAddressOf()));
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateTexture3D");

		if (hRet != S_OK) {
			CxbxrAbort("CreateTexture3D Failed!\n\n"
				"Error: 0x%X\nFormat: %d\nDimensions: %dx%dx%d", hRet, PCFormat, hostWidth, hostHeight, dwDepth);
		}
		SetHostResource(pResource, pNewHostResource.Get(), iTextureStage, D3DUsage);
		EmuLog(LOG_LEVEL::DEBUG, "CreateGpuPixelContainerResource : Successfully created %s (0x%.08X, 0x%.08X)",
			ResourceTypeName, pResource, pNewHostResource.Get());
		break;
	}

	case xbox::X_D3DRTYPE_CUBETEXTURE: {
		EmuLog(LOG_LEVEL::DEBUG, "CreateCubeTexture(%d, %d, 0, %d)", hostWidth,
			dwMipMapLevels, PCFormat);

		D3D11_TEXTURE2D_DESC desc;
		desc.Width = hostWidth;
		desc.Height = hostHeight;
		desc.MipLevels = dwMipMapLevels;
		desc.ArraySize = 6;
		desc.Format = PCFormat;
		desc.SampleDesc.Count = 1;
		desc.SampleDesc.Quality = 0;
		// D3D11_USAGE_DYNAMIC requires ArraySize==1, but cube textures need 6.
		// Always use DEFAULT and upload data via UpdateSubresource.
		desc.Usage = D3D11_USAGE_DEFAULT;
		desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
		if (D3DUsage & D3DUSAGE_RENDERTARGET) {
			desc.BindFlags |= D3D11_BIND_RENDER_TARGET;
		}
		desc.CPUAccessFlags = 0;
		desc.MiscFlags = D3D11_RESOURCE_MISC_TEXTURECUBE;

		hRet = g_pD3DDevice->CreateTexture2D(&desc, NULL, reinterpret_cast<ID3D11Texture2D**>(pNewHostResource.ReleaseAndGetAddressOf()));
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateTexture2D");

		if (hRet != S_OK) {
			CxbxrAbort("CreateCubeTexture Failed!\n\nError: \nDesc: "/*,
				DXGetErrorString(hRet), DXGetErrorDescription(hRet)*/);
		}

		SetHostResource(pResource, pNewHostResource.Get(), iTextureStage, D3DUsage);
		EmuLog(LOG_LEVEL::DEBUG, "CreateGpuPixelContainerResource : Successfully created %s (0x%.08X, 0x%.08X)",
			ResourceTypeName, pResource, pNewHostResource.Get());
		// TODO : Cube face surfaces can be used as a render-target,
		// so we need to associate host surfaces to each surface of this cube texture
   	   	// However, we can't do it here: On Xbox, a new Surface is created on every call to
   	   	// GetCubeMapSurface, so it needs to be done at surface conversion time by looking up
   	   	// the parent CubeTexture
		break;
	}
	} // switch XboxResourceType

	return hRet;
}

// ---- Helper: Create and upload a pixel container (surface/texture/volume) ----
static void CreateHostPixelContainer(
	xbox::X_D3DResource* pResource,
	DWORD D3DUsage,
	int iTextureStage,
	VAddr VirtualAddr,
	const char* ResourceTypeName,
	xbox::X_D3DRESOURCETYPE XboxResourceType)
{
	xbox::X_D3DPixelContainer *pPixelContainer = (xbox::X_D3DPixelContainer*)pResource;
	xbox::X_D3DFORMAT X_Format = GetXboxPixelContainerFormat(pPixelContainer);

	if (EmuXBFormatIsDepthBuffer(X_Format)) {
		D3DUsage |= D3DUSAGE_DEPTHSTENCIL;
	}
	else if (VirtualAddr != 0 && (VirtualAddr & ~PHYSICAL_MAP_BASE) == NV2AGetSurfaceState().colorOffset) {
		// Resource data matches current PGRAPH color render target address
		if (EmuXBFormatIsRenderTarget(X_Format))
			D3DUsage |= D3DUSAGE_RENDERTARGET;
		else
			EmuLog(LOG_LEVEL::WARNING, "Updating RenderTarget %s with an incompatible format!", ResourceTypeName);
	}
	// Determine the host format
	bool bConvertTextureFormat;
	DXGI_FORMAT PCFormat = ResolveHostFormat(X_Format, XboxResourceType, D3DUsage, ResourceTypeName, bConvertTextureFormat);

	// Interpret Width/Height/BPP
	bool bCubemap = pPixelContainer->Format & X_D3DFORMAT_CUBEMAP;
	bool bSwizzled = EmuXBFormatIsSwizzled(X_Format);
	bool bCompressed = EmuXBFormatIsCompressed(X_Format);
	UINT dwBPP = EmuXBFormatBytesPerPixel(X_Format);
	UINT dwMipMapLevels = CxbxGetPixelContainerMipMapLevels(pPixelContainer);
	UINT xboxWidth, xboxHeight, dwDepth, dwRowPitch, dwSlicePitch;

	// Interpret Width/Height/BPP
	CxbxGetPixelContainerMeasures(pPixelContainer, 0, &xboxWidth, &xboxHeight, &dwDepth, &dwRowPitch, &dwSlicePitch);

	// Host width and height dimensions
	UINT hostWidth = xboxWidth;
	UINT hostHeight = xboxHeight;
	// Upscale rendertargets and depth surfaces
	if (D3DUsage & (X_D3DUSAGE_DEPTHSTENCIL | X_D3DUSAGE_RENDERTARGET)) {
		hostWidth *= g_RenderUpscaleFactor;
		hostHeight *= g_RenderUpscaleFactor;
	}

	// Each mip-map level is 1/2 the size of the previous level
	// D3D9 forbids creation of a texture with more mip-map levels than it is divisible
	// EG: A 256x256 texture cannot have more than 8 levels, since that would create a texture smaller than 1x1
	// Because of this, we need to cap dwMipMapLevels when required
	if (dwMipMapLevels > 0) {
		// Calculate how many mip-map levels it takes to get to a texture of 1 pixels in either dimension
		UINT highestMipMapLevel = 1;
		UINT width = xboxWidth; UINT height = xboxHeight;
		while (width > 1 || height > 1) {
			width /= 2;
			height /= 2;
			highestMipMapLevel++;
		}

		// If the desired mip-map level was higher than the maximum possible, cap it
		// Test case: Shin Megami Tensei: Nine
		if (dwMipMapLevels > highestMipMapLevel) {
			LOG_TEST_CASE("Too many mip-map levels");
			dwMipMapLevels = highestMipMapLevel;
		}
	}

	if (dwDepth != 1) {
		// VolumeLight and VolumeTexture samples legitimately use Depth > 1
		// LOG_TEST_CASE("CreateHostPixelContainer : Depth != 1");
	}


	ComPtr<ID3D11Resource> pNewHostResource;
	bool bHostIsDynamic = false;

	CreateGpuPixelContainerResource(
		XboxResourceType,
		hostWidth, hostHeight, dwDepth, dwMipMapLevels,
		PCFormat,
		D3DUsage,
		bSwizzled, bConvertTextureFormat, X_Format,
		pResource, ResourceTypeName, iTextureStage,
		pNewHostResource, bHostIsDynamic
	);

   	// If this resource is a render target or depth stencil, don't attempt to lock/copy it as it won't work anyway
   	// In this case, we simply return
   	if (D3DUsage & D3DUSAGE_RENDERTARGET || D3DUsage & D3DUSAGE_DEPTHSTENCIL) {
   	   	return;
   	}

	UploadPixelContainerMips(
		XboxResourceType, ResourceTypeName, VirtualAddr, iTextureStage,
		X_Format, PCFormat, dwBPP,
		xboxWidth, xboxHeight,
		dwDepth, dwRowPitch, dwSlicePitch,
		dwMipMapLevels,
		bCubemap, bSwizzled, bCompressed, bConvertTextureFormat,
		pNewHostResource, bHostIsDynamic
	);

	// Debug resource dumping
//#define _DEBUG_DUMP_TEXTURE_REGISTER "D:\\"
#ifdef _DEBUG_DUMP_TEXTURE_REGISTER
	bool bDumpConvertedTextures = true;  // TODO : Make this a runtime changeable setting
	if (bDumpConvertedTextures) {
		char szFilePath[MAX_PATH];

		switch (XboxResourceType) {
		case xbox::X_D3DRTYPE_SURFACE: {
			static int dwDumpSurface = 0;
			sprintf(szFilePath, _DEBUG_DUMP_TEXTURE_REGISTER "%.03d-Surface%.03d.dds", X_Format, dwDumpSurface++);
			D3DXSaveSurfaceToFileA(szFilePath, D3DXIFF_DDS, pNewHostSurface, nullptr, nullptr);
			break;
		}
		case xbox::X_D3DRTYPE_VOLUME: {
			// TODO
			break;
		}
		case xbox::X_D3DRTYPE_TEXTURE: {
			static int dwDumpTexure = 0;
			sprintf(szFilePath, _DEBUG_DUMP_TEXTURE_REGISTER "%.03d-Texture%.03d.dds", X_Format, dwDumpTexure++);
			D3DXSaveTextureToFileA(szFilePath, D3DXIFF_DDS, pNewHostTexture, nullptr);
			break;
		}
		case xbox::X_D3DRTYPE_VOLUMETEXTURE: {
			// TODO
			break;
		}
		case xbox::X_D3DRTYPE_CUBETEXTURE: {
			static int dwDumpCubeTexture = 0;
			for (unsigned int face = 0; face <= 5; face++) {
				ID3D11Texture2D *pSurface;
				if (S_OK == pNewHostCubeTexture->GetCubeMapSurface((int)face, 0, &pSurface)) {
					sprintf(szFilePath, _DEBUG_DUMP_TEXTURE_REGISTER "%.03d-CubeTexure%.03d-%d.dds", X_Format, dwDumpCubeTexture, face);
					D3DXSaveSurfaceToFileA(szFilePath, D3DXIFF_DDS, pSurface, nullptr, nullptr);
					pSurface->Release();
				}
			}
			dwDumpCubeTexture++;
			break;
		}
		} // switch XboxResourceType
	}
#endif
}

// Was patch: IDirect3DResource8_Register
void CreateHostResource(xbox::X_D3DResource *pResource, DWORD D3DUsage, int iTextureStage, DWORD dwSize)
{
	if (pResource == xbox::zeroptr)
		return;

	// Determine the resource type name
	const char *ResourceTypeName;
	xbox::X_D3DRESOURCETYPE XboxResourceType = GetXboxD3DResourceType(pResource);

	switch (XboxResourceType) {
	case xbox::X_D3DRTYPE_NONE: ResourceTypeName = "None"; break;
	case xbox::X_D3DRTYPE_SURFACE: ResourceTypeName = "Surface"; break;
	case xbox::X_D3DRTYPE_VOLUME: ResourceTypeName = "Volume"; break;
	case xbox::X_D3DRTYPE_TEXTURE: ResourceTypeName = "Texture"; break;
	case xbox::X_D3DRTYPE_VOLUMETEXTURE: ResourceTypeName = "VolumeTexture"; break;
	case xbox::X_D3DRTYPE_CUBETEXTURE: ResourceTypeName = "CubeTexture"; break;
	case xbox::X_D3DRTYPE_VERTEXBUFFER: ResourceTypeName = "VertexBuffer"; break;
	case xbox::X_D3DRTYPE_INDEXBUFFER: ResourceTypeName = "IndexBuffer"; break;
	case xbox::X_D3DRTYPE_PUSHBUFFER: ResourceTypeName = "PushBuffer"; break;
	case xbox::X_D3DRTYPE_PALETTE: ResourceTypeName = "Palette"; break;
	case xbox::X_D3DRTYPE_FIXUP: ResourceTypeName = "Fixup"; break;
	default:
		EmuLog(LOG_LEVEL::WARNING, "CreateHostResource :-> Unrecognized Xbox Resource Type 0x%.08X", XboxResourceType);
		return;
	}

	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pResource)
		LOG_FUNC_ARG(D3DUsage)
		LOG_FUNC_ARG(iTextureStage)
		LOG_FUNC_ARG(dwSize)
		LOG_FUNC_ARG(ResourceTypeName)
		LOG_FUNC_END;

	// Retrieve and test the xbox resource buffer address
	VAddr VirtualAddr = (VAddr)GetDataFromXboxResource(pResource);
	if ((VirtualAddr & ~PHYSICAL_MAP_BASE) == 0) {
		// TODO: Fix or handle this situation..?
		LOG_TEST_CASE("CreateHostResource : VirtualAddr == 0");
		// This is probably an unallocated resource, mapped into contiguous memory (0x80000000 OR 0xF0000000)
		EmuLog(LOG_LEVEL::WARNING, "CreateHostResource :-> %s carries no data - skipping conversion", ResourceTypeName);
		return;
	}

	switch (XboxResourceType) {
	case xbox::X_D3DRTYPE_NONE: {
		break;
	}

	case xbox::X_D3DRTYPE_SURFACE:
		if (TryResolveParentSurface(pResource, D3DUsage, iTextureStage))
			return;
		[[fallthrough]];
	case xbox::X_D3DRTYPE_VOLUME:
		if (XboxResourceType == xbox::X_D3DRTYPE_VOLUME) {
			if (TryResolveParentVolume(pResource, D3DUsage, iTextureStage))
				return;
		}
		[[fallthrough]];
	case xbox::X_D3DRTYPE_TEXTURE:
	case xbox::X_D3DRTYPE_VOLUMETEXTURE:
	case xbox::X_D3DRTYPE_CUBETEXTURE: {
		CreateHostPixelContainer(pResource, D3DUsage, iTextureStage, VirtualAddr, ResourceTypeName, XboxResourceType);
		break;
	}

	// case X_D3DRTYPE_VERTEXBUFFER: { break; } // TODO
	// case X_D3DRTYPE_INDEXBUFFER: { break; } // TODO
	case xbox::X_D3DRTYPE_PUSHBUFFER: {
		xbox::X_D3DPushBuffer *pPushBuffer = (xbox::X_D3DPushBuffer*)pResource;

		// create push buffer
		dwSize = g_VMManager.QuerySize(VirtualAddr);
		if (dwSize == 0) {
			// TODO: once this is known to be working, remove the warning
			EmuLog(LOG_LEVEL::WARNING, "Push buffer allocation size unknown");
			pPushBuffer->Lock = X_D3DRESOURCE_LOCK_FLAG_NOSIZE;
			break;
		}

		EmuLog(LOG_LEVEL::DEBUG, "CreateHostResource : Successfully created %S (0x%.08X, 0x%.08X, 0x%.08X)", ResourceTypeName, pResource->Data, pPushBuffer->Size, pPushBuffer->AllocationSize);
		break;
	}

	//case X_D3DRTYPE_PALETTE: { break;	}

	case xbox::X_D3DRTYPE_FIXUP: {
		xbox::X_D3DFixup *pFixup = (xbox::X_D3DFixup*)pResource;

		EmuLog(LOG_LEVEL::WARNING, "X_D3DRTYPE_FIXUP is not yet supported\n"
			"0x%.08X (pFixup->Common) \n"
			"0x%.08X (pFixup->Data)   \n"
			"0x%.08X (pFixup->Lock)   \n"
			"0x%.08X (pFixup->Run)    \n"
			"0x%.08X (pFixup->Next)   \n"
			"0x%.08X (pFixup->Size)   \n", pFixup->Common, pFixup->Data, pFixup->Lock, pFixup->Run, pFixup->Next, pFixup->Size);
		break;
	}
   	} // switch XboxResourceType
}

D3DXVECTOR4 toVector(D3DCOLOR color) {
	D3DXVECTOR4 v;
	// ARGB to XYZW
	v.w = (color >> 24 & 0xFF) / 255.f;
	v.x = (color >> 16 & 0xFF) / 255.f;
	v.y = (color >> 8 & 0xFF) / 255.f;
	v.z = (color >> 0 & 0xFF) / 255.f;
	return v;
}

D3DXVECTOR4 toVector(xbox::X_D3DCOLORVALUE val) {
	return D3DXVECTOR4(val.r, val.g, val.b, val.a);
}

