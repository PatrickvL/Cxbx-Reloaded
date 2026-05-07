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
// *  (c) 2002-2003 Aaron Robinson <caustik@caustik.com>
// *
// *  All rights reserved
// *
// ******************************************************************
#include "EmuD3D8_common.h"

XboxRenderStateConverter XboxRenderStates;

FixedFunctionVertexShaderState ffShaderState = {}; // TODO find a home for this and associated code

// Allow use of time duration literals (making 16ms, etc possible)
using namespace std::literals::chrono_literals;

// Global(s)
HWND                         g_hEmuWindow   = NULL; // rendering window
bool                         g_bClipCursor  = false; // indicates that the mouse cursor should be confined inside the rendering window
ID3D11Device                *g_pD3DDevice   = nullptr; // Direct3D Device

// Shared Variable(s)
bool                         g_bSupportsFormatSurface[xbox::X_D3DFMT_LAST + 1]; // Does device support surface format?
bool                         g_bSupportsFormatSurfaceRenderTarget[xbox::X_D3DFMT_LAST + 1]; // Does device support surface format?
bool                         g_bSupportsFormatSurfaceDepthStencil[xbox::X_D3DFMT_LAST + 1]; // Does device support surface format?
bool                         g_bSupportsFormatTexture[xbox::X_D3DFMT_LAST + 1]; // Does device support texture format?
bool                         g_bSupportsFormatTextureRenderTarget[xbox::X_D3DFMT_LAST + 1]; // Does device support texture format?
bool                         g_bSupportsFormatTextureDepthStencil[xbox::X_D3DFMT_LAST + 1]; // Does device support texture format?
bool                         g_bSupportsFormatVolumeTexture[xbox::X_D3DFMT_LAST + 1]; // Does device support surface format?
bool                         g_bSupportsFormatCubeTexture[xbox::X_D3DFMT_LAST + 1]; // Does device support surface format?
bool                         g_bHack_UnlockFramerate = false; // ignore the xbox presentation interval
bool                         g_bHasDepth = false;    // Does device have a Depth Buffer?
bool                         g_bHasStencil = false;  // Does device have a Stencil Buffer?

float                        g_AspectRatioScale = 1.0f;
UINT                         g_AspectRatioScaleWidth = 0;
UINT                         g_AspectRatioScaleHeight = 0;
D3D11_TEXTURE2D_DESC         g_HostBackBufferDesc;
Settings::s_video            g_XBVideo;

bool                         g_bEnableHostQueryVisibilityTest = true;

bool                         g_bHack_DisableHostGPUQueries = false; // TODO : Make configurable
ID3D11Query                 *g_pHostQueryWaitForIdle = nullptr;
int                          g_RenderUpscaleFactor = 1;
xbox::X_D3DMULTISAMPLE_TYPE  g_Xbox_MultiSampleType = xbox::X_D3DMULTISAMPLE_NONE;

// Side-map: VRAM data offset → Xbox texture pointer
static std::unordered_map<xbox::addr_xt, xbox::X_D3DBaseTexture*> g_TexturesByDataAddr;

void CxbxRegisterTextureByDataAddr(xbox::addr_xt dataAddr, xbox::X_D3DBaseTexture *pTexture)
{
	if (dataAddr != xbox::zero)
		g_TexturesByDataAddr[dataAddr] = pTexture;
}

xbox::X_D3DBaseTexture* CxbxLookupTextureByDataAddr(xbox::addr_xt dataAddr)
{
	auto it = g_TexturesByDataAddr.find(dataAddr);
	return (it != g_TexturesByDataAddr.end()) ? it->second : nullptr;
}

DXGI_FORMAT               g_HostTextureFormats[xbox::X_D3DTS_STAGECOUNT]; // Updated by CxbxUpdateHostTextures(), read by CxbxCalcColorSign
xbox::X_D3DBaseTexture       *g_pXbox_SetTexture[xbox::X_D3DTS_STAGECOUNT] = {0,0,0,0}; // Set by our D3DDevice_SetTexture and D3DDevice_SwitchTexture patches

EmuD3D8CreateDeviceProxyData g_EmuCDPD;

// Define trampolines (XB_TRAMPOLINES macro is defined in RenderGlobals.h)
// Non-static so they are accessible from other translation units.
#define XB_trampoline_define(ret, conv, func, arguments) \
   	static const std::string XB_NAME(func) = #func; \
   	XB_TYPE(func) XB_TRMP(func) = nullptr

XB_TRAMPOLINES(XB_trampoline_define);
#undef XB_trampoline_define

void LookupTrampolinesD3D()
{
	XB_TRAMPOLINES(XB_trampoline_lookup);
}

const char *CxbxGetErrorDescription(HRESULT hResult)
{
	// TODO : Use FormatMessage with FORMAT_MESSAGE_FROM_SYSTEM for DirectX errors
	switch (hResult)
	{
	case D3DERR_INVALIDCALL: return "Invalid Call";
	case D3DERR_NOTAVAILABLE: return "Not Available";
	// case D3DERR_OUTOFVIDEOMEMORY: return "Out of Video Memory"; // duplicate of DDERR_OUTOFVIDEOMEMORY

	case S_OK: return "No error occurred.";
	case D3DERR_CONFLICTINGTEXTUREFILTER: return "The current texture filters cannot be used together. ";
	case D3DERR_CONFLICTINGTEXTUREPALETTE: return "The current textures cannot be used simultaneously. This generally occurs when a multitexture device requires that all palettized textures simultaneously enabled also share the same palette. ";
	case D3DERR_CONFLICTINGRENDERSTATE: return "The currently set render states cannot be used together. ";
	case D3DERR_TOOMANYOPERATIONS: return "The application is requesting more texture-filtering operations than the device supports. ";
	case D3DERR_UNSUPPORTEDALPHAARG: return "The device does not support one of the specified texture-blending arguments for the alpha channel. ";
	case D3DERR_UNSUPPORTEDALPHAOPERATION: return "The device does not support one of the specified texture-blending operations for the alpha channel. ";
	case D3DERR_UNSUPPORTEDCOLORARG: return "The device does not support one of the specified texture-blending arguments for color values. ";
	case D3DERR_UNSUPPORTEDCOLOROPERATION: return "The device does not support one of the specified texture-blending operations for color values. ";
	case D3DERR_UNSUPPORTEDFACTORVALUE: return "The specified texture factor value is not supported by the device. ";
	case D3DERR_UNSUPPORTEDTEXTUREFILTER: return "The specified texture filter is not supported by the device. ";
	case D3DERR_WRONGTEXTUREFORMAT: return "The pixel format of the texture surface is not valid. ";
	}

	return nullptr;
}

// A helper function to run any code on a window message thread
// Used for those D3D11 functions which *must* be run on this particular thread
