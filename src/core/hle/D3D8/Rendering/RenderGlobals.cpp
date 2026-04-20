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
XboxTextureStateConverter XboxTextureStates;

D3D8LightState d3d8LightState = D3D8LightState();
D3D8TransformState d3d8TransformState = D3D8TransformState();
FixedFunctionVertexShaderState ffShaderState = {}; // TODO find a home for this and associated code

// Allow use of time duration literals (making 16ms, etc possible)
using namespace std::literals::chrono_literals;

// Global(s)
HWND                                g_hEmuWindow   = NULL; // rendering window
bool                                g_bClipCursor  = false; // indicates that the mouse cursor should be confined inside the rendering window
ID3D11Device                    *g_pD3DDevice   = nullptr; // Direct3D Device

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
float                        g_ZScale = 1.0;
bool                         g_bHasStencil = false;  // Does device have a Stencil Buffer?
DWORD						 g_dwPrimPerFrame = 0;	// Number of primitives within one frame
float                        g_AspectRatioScale = 1.0f;
UINT                         g_AspectRatioScaleWidth = 0;
UINT                         g_AspectRatioScaleHeight = 0;
D3D11_TEXTURE2D_DESC               g_HostBackBufferDesc;
Settings::s_video            g_XBVideo;

bool                         g_bEnableHostQueryVisibilityTest = true;

OverlayProxy g_OverlayProxy;

bool                         g_bHack_DisableHostGPUQueries = false; // TODO : Make configurable
ID3D11Query              *g_pHostQueryWaitForIdle = nullptr;
ID3D11Query              *g_pHostQueryCallbackEvent = nullptr;
int                          g_RenderUpscaleFactor = 1;

xbox::dword_xt                   g_Xbox_PresentationInterval_Default = D3DPRESENT_INTERVAL_IMMEDIATE;
xbox::dword_xt                   g_Xbox_PresentationInterval_Override = 0;
xbox::X_D3DSWAPDATA			g_Xbox_SwapData = {0}; // current swap information
xbox::X_D3DSWAPCALLBACK		g_pXbox_SwapCallback = xbox::zeroptr;	// Swap/Present callback routine
xbox::X_D3DVBLANKDATA			g_Xbox_VBlankData = {0}; // current vertical blank information

xbox::X_D3DSurface           *g_pXbox_BackBufferSurface = xbox::zeroptr;
xbox::X_D3DSurface           *g_pXbox_DefaultDepthStencilSurface = xbox::zeroptr;
xbox::X_D3DSurface           *g_pXbox_RenderTarget = xbox::zeroptr;
xbox::X_D3DSurface           *g_pXbox_DepthStencil = xbox::zeroptr;
xbox::X_D3DMULTISAMPLE_TYPE   g_Xbox_MultiSampleType = xbox::X_D3DMULTISAMPLE_NONE;

xbox::X_VERTEXSHADERCONSTANTMODE g_Xbox_VertexShaderConstantMode = X_D3DSCM_192CONSTANTS; // Set by D3DDevice_SetShaderConstantMode, TODO : Move to XbVertexShader.cpp
xbox::dword_xt                   g_Xbox_BaseVertexIndex = 0; // Set by D3DDevice_SetIndices, read by D3DDevice_DrawIndexedVertices

xbox::X_PixelShader*			g_pXbox_PixelShader = xbox::zeroptr;
xbox::PVOID                   g_pXbox_Palette_Data[xbox::X_D3DTS_STAGECOUNT] = { xbox::zeroptr, xbox::zeroptr, xbox::zeroptr, xbox::zeroptr }; // cached palette pointer
unsigned                     g_Xbox_Palette_Size[xbox::X_D3DTS_STAGECOUNT] = { 0 }; // cached palette size


DXGI_FORMAT               g_HostTextureFormats[xbox::X_D3DTS_STAGECOUNT]; // Updated by CxbxUpdateHostTextures(), read by CxbxCalcColorSign
xbox::X_D3DBaseTexture       *g_pXbox_SetTexture[xbox::X_D3DTS_STAGECOUNT] = {0,0,0,0}; // Set by our D3DDevice_SetTexture and D3DDevice_SwitchTexture patches

xbox::X_D3DVIEWPORT8 g_Xbox_Viewport = { 0 };
float g_Xbox_BackbufferScaleX = 1;
float g_Xbox_BackbufferScaleY = 1;
xbox::X_D3DSWAP g_LastD3DSwap = (xbox::X_D3DSWAP) -1;


static void							UpdateCurrentMSpFAndFPS(); // Used for benchmarking/fps count


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
#if 0
	case D3DERR_BADMAJORVERSION: return "The service that you requested is unavailable in this major version of DirectX. (A major version denotes a primary release, such as DirectX 6.0.) ";
	case D3DERR_BADMINORVERSION: return "The service that you requested is available in this major version of DirectX, but not in this minor version. Get the latest version of the component run time from Microsoft. (A minor version denotes a secondary release, such as DirectX 6.1.) ";
	case D3DERR_COLORKEYATTACHED: return "The application attempted to create a texture with a surface that uses a color key for transparency. ";
#endif
	case D3DERR_CONFLICTINGTEXTUREFILTER: return "The current texture filters cannot be used together. ";
	case D3DERR_CONFLICTINGTEXTUREPALETTE: return "The current textures cannot be used simultaneously. This generally occurs when a multitexture device requires that all palettized textures simultaneously enabled also share the same palette. ";
	case D3DERR_CONFLICTINGRENDERSTATE: return "The currently set render states cannot be used together. ";
#if 0
	case D3DERR_DEVICEAGGREGATED: return "The IDirect3DDevice7::SetRenderTarget method was called on a device that was retrieved from the render target surface. ";
	case D3DERR_EXECUTE_CLIPPED_FAILED: return "The execute buffer could not be clipped during execution. ";
	case D3DERR_EXECUTE_CREATE_FAILED: return "The execute buffer could not be created. This typically occurs when no memory is available to allocate the execute buffer. ";
	case D3DERR_EXECUTE_DESTROY_FAILED: return "The memory for the execute buffer could not be deallocated. ";
	case D3DERR_EXECUTE_FAILED: return "The contents of the execute buffer are invalid and cannot be executed. ";
	case D3DERR_EXECUTE_LOCK_FAILED: return "The execute buffer could not be locked. ";
	case D3DERR_EXECUTE_LOCKED: return "The operation requested by the application could not be completed because the execute buffer is locked. ";
	case D3DERR_EXECUTE_NOT_LOCKED: return "The execute buffer could not be unlocked because it is not currently locked. ";
	case D3DERR_EXECUTE_UNLOCK_FAILED: return "The execute buffer could not be unlocked. ";
	case D3DERR_INBEGIN: return "The requested operation cannot be completed while scene rendering is taking place. Try again after the scene is completed and the IDirect3DDevice7::EndScene method is called. ";
	case D3DERR_INBEGINSTATEBLOCK: return "The operation cannot be completed while recording states for a state block. Complete recording by calling the IDirect3DDevice7::EndStateBlock method, and try again. ";
	case D3DERR_INITFAILED: return "A rendering device could not be created because the new device could not be initialized. ";
	case D3DERR_INVALID_DEVICE: return "The requested device type is not valid. ";
	case D3DERR_INVALIDCURRENTVIEWPORT: return "The currently selected viewport is not valid. ";
	case D3DERR_INVALIDMATRIX: return "The requested operation could not be completed because the combination of the currently set world, view, and projection matrices is invalid (the determinant of the combined matrix is 0). ";
	case D3DERR_INVALIDPALETTE: return "The palette associated with a surface is invalid. ";
	case D3DERR_INVALIDPRIMITIVETYPE: return "The primitive type specified by the application is invalid. ";
	case D3DERR_INVALIDRAMPTEXTURE: return "Ramp mode is being used, and the texture handle in the current material does not match the current texture handle that is set as a render state. ";
	case D3DERR_INVALIDSTATEBLOCK: return "The state block handle is invalid. ";
	case D3DERR_INVALIDVERTEXFORMAT: return "The combination of flexible vertex format flags specified by the application is not valid. ";
	case D3DERR_INVALIDVERTEXTYPE: return "The vertex type specified by the application is invalid. ";
	case D3DERR_LIGHT_SET_FAILED: return "The attempt to set lighting parameters for a light object failed. ";
	case D3DERR_LIGHTHASVIEWPORT: return "The requested operation failed because the light object is associated with another viewport. ";
	case D3DERR_LIGHTNOTINTHISVIEWPORT: return "The requested operation failed because the light object has not been associated with this viewport. ";
	case D3DERR_MATERIAL_CREATE_FAILED: return "The material could not be created. This typically occurs when no memory is available to allocate for the material. ";
	case D3DERR_MATERIAL_DESTROY_FAILED: return "The memory for the material could not be deallocated. ";
	case D3DERR_MATERIAL_GETDATA_FAILED: return "The material parameters could not be retrieved. ";
	case D3DERR_MATERIAL_SETDATA_FAILED: return "The material parameters could not be set. ";
	case D3DERR_MATRIX_CREATE_FAILED: return "The matrix could not be created. This can occur when no memory is available to allocate for the matrix. ";
	case D3DERR_MATRIX_DESTROY_FAILED: return "The memory for the matrix could not be deallocated. ";
	case D3DERR_MATRIX_GETDATA_FAILED: return "The matrix data could not be retrieved. This can occur when the matrix was not created by the current device. ";
	case D3DERR_MATRIX_SETDATA_FAILED: return "The matrix data could not be set. This can occur when the matrix was not created by the current device. ";
	case D3DERR_NOCURRENTVIEWPORT: return "The viewport parameters could not be retrieved because none have been set. ";
	case D3DERR_NOTINBEGIN: return "The requested rendering operation could not be completed because scene rendering has not begun. Call IDirect3DDevice7::BeginScene to begin rendering, and try again. ";
	case D3DERR_NOTINBEGINSTATEBLOCK: return "The requested operation could not be completed because it is only valid while recording a state block. Call the IDirect3DDevice7::BeginStateBlock method, and try again. ";
	case D3DERR_NOVIEWPORTS: return "The requested operation failed because the device currently has no viewports associated with it. ";
	case D3DERR_SCENE_BEGIN_FAILED: return "Scene rendering could not begin. ";
	case D3DERR_SCENE_END_FAILED: return "Scene rendering could not be completed. ";
	case D3DERR_SCENE_IN_SCENE: return "Scene rendering could not begin because a previous scene was not completed by a call to the IDirect3DDevice7::EndScene method. ";
	case D3DERR_SCENE_NOT_IN_SCENE: return "Scene rendering could not be completed because a scene was not started by a previous call to the IDirect3DDevice7::BeginScene method. ";
	case D3DERR_SETVIEWPORTDATA_FAILED: return "The viewport parameters could not be set. ";
	case D3DERR_STENCILBUFFER_NOTPRESENT: return "The requested stencil buffer operation could not be completed because there is no stencil buffer attached to the render target surface. ";
	case D3DERR_SURFACENOTINVIDMEM: return "The device could not be created because the render target surface is not located in video memory. (Hardware-accelerated devices require video-memory render target surfaces.) ";
	case D3DERR_TEXTURE_BADSIZE: return "The dimensions of a current texture are invalid. This can occur when an application attempts to use a texture that has dimensions that are not a power of 2 with a device that requires them. ";
	case D3DERR_TEXTURE_CREATE_FAILED: return "The texture handle for the texture could not be retrieved from the driver. ";
	case D3DERR_TEXTURE_DESTROY_FAILED: return "The device was unable to deallocate the texture memory. ";
	case D3DERR_TEXTURE_GETSURF_FAILED: return "The DirectDraw surface used to create the texture could not be retrieved. ";
	case D3DERR_TEXTURE_LOAD_FAILED: return "The texture could not be loaded. ";
	case D3DERR_TEXTURE_LOCK_FAILED: return "The texture could not be locked. ";
	case D3DERR_TEXTURE_LOCKED: return "The requested operation could not be completed because the texture surface is currently locked. ";
	case D3DERR_TEXTURE_NO_SUPPORT: return "The device does not support texture mapping. ";
	case D3DERR_TEXTURE_NOT_LOCKED: return "The requested operation could not be completed because the texture surface is not locked. ";
	case D3DERR_TEXTURE_SWAP_FAILED: return "The texture handles could not be swapped. ";
	case D3DERR_TEXTURE_UNLOCK_FAILED: return "The texture surface could not be unlocked. ";
#endif
	case D3DERR_TOOMANYOPERATIONS: return "The application is requesting more texture-filtering operations than the device supports. ";
#if 0
	case D3DERR_TOOMANYPRIMITIVES: return "The device is unable to render the provided number of primitives in a single pass. ";
#endif
	case D3DERR_UNSUPPORTEDALPHAARG: return "The device does not support one of the specified texture-blending arguments for the alpha channel. ";
	case D3DERR_UNSUPPORTEDALPHAOPERATION: return "The device does not support one of the specified texture-blending operations for the alpha channel. ";
	case D3DERR_UNSUPPORTEDCOLORARG: return "The device does not support one of the specified texture-blending arguments for color values. ";
	case D3DERR_UNSUPPORTEDCOLOROPERATION: return "The device does not support one of the specified texture-blending operations for color values. ";
	case D3DERR_UNSUPPORTEDFACTORVALUE: return "The specified texture factor value is not supported by the device. ";
	case D3DERR_UNSUPPORTEDTEXTUREFILTER: return "The specified texture filter is not supported by the device. ";
#if 0
	case D3DERR_VBUF_CREATE_FAILED: return "The vertex buffer could not be created. This can happen when there is insufficient memory to allocate a vertex buffer. ";
	case D3DERR_VERTEXBUFFERLOCKED: return "The requested operation could not be completed because the vertex buffer is locked. ";
	case D3DERR_VERTEXBUFFEROPTIMIZED: return "The requested operation could not be completed because the vertex buffer is optimized. (The contents of optimized vertex buffers are driver-specific and considered private.) ";
	case D3DERR_VERTEXBUFFERUNLOCKFAILED: return "The vertex buffer could not be unlocked because the vertex buffer memory was overrun. Be sure that your application does not write beyond the size of the vertex buffer. ";
	case D3DERR_VIEWPORTDATANOTSET: return "The requested operation could not be completed because viewport parameters have not yet been set. Set the viewport parameters by calling the IDirect3DDevice7::SetViewport method, and try again. ";
	case D3DERR_VIEWPORTHASNODEVICE: return "The requested operation could not be completed because the viewport has not yet been associated with a device. Associate the viewport with a rendering device by calling the IDirect3DDevice3::AddViewport method, and try again. ";
#endif
	case D3DERR_WRONGTEXTUREFORMAT: return "The pixel format of the texture surface is not valid. ";
#if 0
	case D3DERR_ZBUFF_NEEDS_SYSTEMMEMORY: return "The requested operation could not be completed because the specified device requires system-memory depth-buffer surfaces. (Software rendering devices require system-memory depth buffers.) ";
	case D3DERR_ZBUFF_NEEDS_VIDEOMEMORY: return "The requested operation could not be completed because the specified device requires video-memory depth-buffer surfaces. (Hardware-accelerated devices require video-memory depth buffers.) ";
	case D3DERR_ZBUFFER_NOTPRESENT: return "The requested operation could not be completed because the render target surface does not have an attached depth buffer. ";
	case DD_OK: return "The request completed successfully.";
#endif
	}

	return nullptr;
}

// A helper function to run any code on a window message thread
// Used for those D3D11 functions which *must* be run on this particular thread
