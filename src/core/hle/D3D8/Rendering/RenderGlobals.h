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
#ifndef RENDERGLOBALS_H
#define RENDERGLOBALS_H

#include "core\hle\XAPI\Xapi.h" // For EMUPATCH

#ifndef VOID
#define VOID void
#endif

#include "core\hle\D3D8\XbD3D8Types.h"
#include "common\Settings.hpp" // For Settings::s_video
#include "common\util\hasher.h" // For ComputeHash

#include <queue>
#include <stack>
#include <map>
#include <unordered_map>
#include <wrl/client.h>

#define CXBX_D3DCOMMON_IDENTIFYING_MASK (X_D3DCOMMON_TYPE_MASK | X_D3DCOMMON_D3DCREATED)

typedef struct resource_key_hash {
	DWORD Common;
	DWORD Data;
	union {
		struct {
			xbox::addr_xt ResourceAddr;
		};
		struct {
			DWORD Format;
			DWORD Size;
			uint64_t PaletteHash;
		};
	};

	bool operator==(const struct resource_key_hash& other) const
	{
		return (Common == other.Common)
			&& (Data == other.Data)
			&& (Format == other.Format)
			&& (Size == other.Size)
			&& (PaletteHash == other.PaletteHash);
	}

	size_t operator()(const struct resource_key_hash& value) const
	{
		return (size_t)ComputeHash((void*)&value, sizeof(value));
	}
} resource_key_t;

// LTCG prologue/epilogue macros for __declspec(naked) patches
#define LTCG_PROLOGUE \
		__asm push ebp \
		__asm mov  ebp, esp \
		__asm sub  esp, __LOCAL_SIZE \
		__asm push esi \
		__asm push edi \
		__asm push ebx

#define LTCG_EPILOGUE \
		__asm pop  ebx \
		__asm pop  edi \
		__asm pop  esi \
		__asm mov  esp, ebp \
		__asm pop  ebp

constexpr UINT WM_CXBXR_RUN_ON_MESSAGE_THREAD = WM_USER+0;

// Shared global variables (defined in RenderGlobals.cpp)
struct FixedFunctionVertexShaderState;
extern FixedFunctionVertexShaderState ffShaderState;
extern HWND                          g_hEmuWindow;
extern bool                          g_bClipCursor;
extern bool                          g_bSupportsFormatSurface[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatSurfaceRenderTarget[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatSurfaceDepthStencil[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatTexture[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatTextureRenderTarget[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatTextureDepthStencil[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatVolumeTexture[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bSupportsFormatCubeTexture[xbox::X_D3DFMT_LAST + 1];
extern bool                          g_bHack_UnlockFramerate;
extern bool                          g_bHasDepth;
extern bool                          g_bHasStencil;
extern float                         g_AspectRatioScale;
extern UINT                          g_AspectRatioScaleWidth;
extern UINT                          g_AspectRatioScaleHeight;
extern Settings::s_video             g_XBVideo;
extern bool                          g_bEnableHostQueryVisibilityTest;
extern bool                          g_bHack_DisableHostGPUQueries;
extern int                           g_RenderUpscaleFactor;
extern xbox::X_D3DMULTISAMPLE_TYPE   g_Xbox_MultiSampleType;

// Side-map: VRAM data offset → Xbox texture pointer.
// Populated by D3DDevice_SetTexture/SwitchTexture patches; queried by
// CxbxUpdateHostTextures to resolve PGRAPH TEXOFFSET values to Xbox textures.
void CxbxRegisterTextureByDataAddr(xbox::addr_xt dataAddr, xbox::X_D3DBaseTexture *pTexture);
xbox::X_D3DBaseTexture* CxbxLookupTextureByDataAddr(xbox::addr_xt dataAddr);

void LookupTrampolinesD3D();

// initialize render window
extern void CxbxInitWindow();

// Save window state (position, size, faux fullscreen) to shared memory for reboot persistence
extern void CxbxSaveWindowStateForReboot();

void CxbxUpdateNativeD3DResources();

// Shader constant helpers
void CxbxSetVertexShaderConstantF(UINT startRegister, const float* pConstantData, UINT Vector4fCount);

// Rendering helpers (implemented in Backend_D3D11*.cpp)
void    CxbxD3DClear(DWORD Count, CONST D3DRECT* pRects, DWORD Flags, D3DCOLOR Color, float Z, DWORD Stencil);
void    CxbxSetScissorRect(CONST RECT *pHostViewportRect);
HRESULT CxbxPresent();
void    CxbxInvalidateActivePixelShader(); // Reset PS state tracking after blit/present


// BeginScene/EndScene: no-ops in D3D11
inline void CxbxBeginScene()
{
}

inline void CxbxEndScene()
{
}

// Xbox resource type helpers (used in HostResource.cpp, HostSync.cpp, etc.)
inline DWORD GetXboxCommonResourceType(const xbox::dword_xt XboxResource_Common)
{
	DWORD dwCommonType = XboxResource_Common & X_D3DCOMMON_TYPE_MASK;
	return dwCommonType;
}

inline DWORD GetXboxCommonResourceType(const xbox::X_D3DResource* pXboxResource)
{
	// Don't pass in unassigned Xbox resources
	assert(pXboxResource != xbox::zeroptr);
	return GetXboxCommonResourceType(pXboxResource->Common);
}

inline int GetXboxPixelContainerDimensionCount(const xbox::X_D3DPixelContainer *pXboxPixelContainer)
{
	// Don't pass in unassigned Xbox pixel container
	assert(pXboxPixelContainer != xbox::zeroptr);
	return (xbox::X_D3DFORMAT)((pXboxPixelContainer->Format & X_D3DFORMAT_DIMENSION_MASK) >> X_D3DFORMAT_DIMENSION_SHIFT);
}

// initialize direct3d
extern void EmuD3DInit();

// cleanup direct3d
extern void EmuD3DCleanup();

extern DXGI_FORMAT g_HostTextureFormats[xbox::X_D3DTS_STAGECOUNT];

extern xbox::X_D3DBaseTexture *g_pXbox_SetTexture[xbox::X_D3DTS_STAGECOUNT];


// Cross-file function declarations (defined in various Host*.cpp files)

// Swap present forward marker
#define CXBX_SWAP_PRESENT_FORWARD (256 + X_D3DSWAP_FINISH + X_D3DSWAP_COPY) // = CxbxPresentForwardMarker + D3DSWAP_FINISH + D3DSWAP_COPY

// HostDevice.cpp
void UpdateHostBackBufferDesc();
void SetAspectRatioScale(const xbox::X_D3DPRESENT_PARAMETERS* pPresentationParameters);
DWORD WINAPI EmuRenderWindow(LPVOID lpParam);
void GetMultiSampleOffset(float& xOffset, float& yOffset);
void GetMultiSampleScaleRaw(float& xScale, float& yScale);
void GetScreenScaleFactors(float& scaleX, float& scaleY);
void GetRenderTargetBaseDimensions(float& x, float& y);
void SetupPresentationParameters(const xbox::X_D3DPRESENT_PARAMETERS *pXboxPresentationParameters);
void DetermineSupportedD3DFormats();

// HostDevice.cpp
void SetXboxMultiSampleType(xbox::X_D3DMULTISAMPLE_TYPE value);

// HostRender.cpp
void CxbxInitHostD3DDevice();
void CreateDefaultDevice(const xbox::X_D3DPRESENT_PARAMETERS* pPresentationParameters);

// HostResource.cpp
xbox::X_D3DRESOURCETYPE GetXboxD3DResourceType(const xbox::X_D3DResource* pXboxResource);
void* GetDataFromXboxResource(xbox::X_D3DResource* pXboxResource);
resource_key_t GetHostResourceKey(xbox::X_D3DResource* pXboxResource, int iTextureStage = -1);
void ClearAllResourceCaches();
uint32_t GetPixelContainerWidth(xbox::X_D3DPixelContainer *pPixelContainer);
uint32_t GetPixelContainerHeight(xbox::X_D3DPixelContainer *pPixelContainer);
unsigned int CxbxGetPixelContainerMipMapLevels(xbox::X_D3DPixelContainer *pPixelContainer);
void GetSurfaceFaceAndLevelWithinTexture(xbox::X_D3DSurface* pSurface, xbox::X_D3DBaseTexture* pTexture, UINT& Level, int& Face);
void GetSurfaceFaceAndLevelWithinTexture(xbox::X_D3DSurface* pSurface, xbox::X_D3DBaseTexture* pBaseTexture, UINT& Level);
void PrunePaletizedTexturesCache();
void PruneResourceCache();
void CxbxGetPixelContainerMeasures(xbox::X_D3DPixelContainer *pPixelContainer, DWORD dwMipMapLevel, UINT *pWidth, UINT *pHeight, UINT *pDepth, UINT *pRowPitch, UINT *pSlicePitch);
void CreateHostResource(xbox::X_D3DResource *pResource, DWORD D3DUsage, int iTextureStage, DWORD dwSize);
int XboxD3DPaletteSizeToBytes(const xbox::X_D3DPALETTESIZE Size);
xbox::X_D3DPALETTESIZE GetXboxPaletteSize(const xbox::X_D3DPalette *pPalette);

// HostRender.cpp
void UpdateFixedFunctionVertexShaderState();
void InvalidateFixedFunctionStateCache(); // Call after VP constants overwrite the shared cbuffer
void CxbxUpdateHostViewPortOffsetAndScaleConstants();

// HostResourceCreate.cpp
D3DXVECTOR4 toVector(D3DCOLOR color);
D3DXVECTOR4 toVector(xbox::X_D3DCOLORVALUE val);

// HostWindow.cpp
void DrawUEM(HWND hWnd);
void CxbxReleaseCursor();
void CxbxClipCursor(HWND hWnd);
void CxbxUpdateCursor(bool forceShow = false);
void RunOnWndMsgThread(const std::function<void()>& func);

const char *CxbxGetErrorDescription(HRESULT hResult);

// CxbxrImpl_GetBackBuffer2 — removed (only used by disabled GetBackBuffer patches).
// xbox::X_D3DSurface* CxbxrImpl_GetBackBuffer2(xbox::int_xt BackBuffer);

// Xbox function trampolines -- defined in RenderGlobals.cpp.
// XB_TRAMPOLINES lists all D3D8 trampolines; we invoke it here with an
// "extern" generator so every TU that includes this header sees the
// declarations.  The actual definitions live in RenderGlobals.cpp.
#define XB_TRAMPOLINES(XB_MACRO)                                                                                                                                                                            \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_CreateVertexShader,                       (CONST xbox::dword_xt*, CONST xbox::dword_xt*, xbox::dword_xt*, xbox::dword_xt)                       );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_DeleteVertexShader,                       (xbox::dword_xt)                                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_DeleteVertexShader_0__LTCG_eax1,          ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_GetBackBuffer,                            (xbox::int_xt, D3DBACKBUFFER_TYPE, xbox::X_D3DSurface**)                                              );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_GetBackBuffer_8__LTCG_eax1,               (D3DBACKBUFFER_TYPE, xbox::X_D3DSurface**)                                                            );  \
   	XB_MACRO(xbox::X_D3DSurface*, WINAPI,     D3DDevice_GetBackBuffer2,                           (xbox::int_xt)                                                                                        );  \
   	XB_MACRO(xbox::X_D3DSurface*, WINAPI,     D3DDevice_GetBackBuffer2_0__LTCG_eax1,              ()                                                                                                    );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_GetDepthStencilSurface,                   (xbox::X_D3DSurface**)                                                                                );  \
   	XB_MACRO(xbox::X_D3DSurface*, WINAPI,     D3DDevice_GetDepthStencilSurface2,                  (xbox::void_xt)                                                                                       );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_GetDisplayMode,                           (xbox::X_D3DDISPLAYMODE*)                                                                             );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_GetRenderTarget,                          (xbox::X_D3DSurface**)                                                                                );  \
   	XB_MACRO(xbox::X_D3DSurface*, WINAPI,     D3DDevice_GetRenderTarget2,                         (xbox::void_xt)                                                                                       );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_LightEnable,                              (xbox::dword_xt, xbox::bool_xt)                                                                       );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_LightEnable_4__LTCG_eax1,                 (xbox::bool_xt)                                                                                       );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_PersistDisplay,                           (xbox::void_xt)                                                                                       );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_Reset,                                    (xbox::X_D3DPRESENT_PARAMETERS*)                                                                      );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_Reset_0__LTCG_edi1,                       ()                                                                                                    );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_Reset_0__LTCG_ebx1,                       ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetBackMaterial,                          (CONST xbox::X_D3DMATERIAL8*)                                                                         );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetIndices,                               (xbox::X_D3DIndexBuffer*, xbox::uint_xt)                                                              );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetIndices_4__LTCG_ebx1,                  (xbox::uint_xt)                                                                                       );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_SetLight,                                 (xbox::dword_xt, CONST xbox::X_D3DLIGHT8*)                                                            );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetMaterial,                              (CONST xbox::X_D3DMATERIAL8*)                                                                         );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetPixelShader,                           (xbox::dword_xt)                                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetPixelShader_0__LTCG_eax1,              ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       __fastcall, D3DDevice_SetRenderState_Simple,                    (xbox::dword_xt, xbox::dword_xt)                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetRenderTarget,                          (xbox::X_D3DSurface*, xbox::X_D3DSurface*)                                                            );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetRenderTarget_0__LTCG_ecx1_eax2,        ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetStreamSource,                          (xbox::uint_xt, xbox::X_D3DVertexBuffer*, xbox::uint_xt)                                              );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetStreamSource_0__LTCG_eax1_edi2_ebx3,   ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetStreamSource_4__LTCG_eax1_ebx2,        (xbox::uint_xt)                                                                                       );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetStreamSource_8__LTCG_eax1,             (xbox::X_D3DVertexBuffer*, xbox::uint_xt)                                                             );  \
   	XB_MACRO(xbox::void_xt,       __fastcall, D3DDevice_SetStreamSource_8__LTCG_edx1,             (void*, xbox::uint_xt, xbox::X_D3DVertexBuffer*, xbox::uint_xt)                                       );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetTexture,                               (xbox::dword_xt, xbox::X_D3DBaseTexture*)                                                             );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetTexture_4__LTCG_eax2,                  (xbox::dword_xt)                                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetTexture_4__LTCG_eax1,                  (xbox::X_D3DBaseTexture*)                                                                             );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetPalette,                               (xbox::dword_xt, xbox::X_D3DPalette*)                                                                 );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetPalette_4__LTCG_eax1,                  (xbox::X_D3DPalette*)                                                                                 );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetVertexShader,                          (xbox::dword_xt)                                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetVertexShader_0__LTCG_ebx1,             ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_LoadVertexShader,                         (xbox::dword_xt, xbox::dword_xt)                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SelectVertexShader,                       (xbox::dword_xt, xbox::dword_xt)                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetVertexShaderInput,                     (xbox::dword_xt, xbox::uint_xt, xbox::X_STREAMINPUT*)                                                 );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetViewport,                              (CONST xbox::X_D3DVIEWPORT8*)                                                                         );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetTransform,                             (xbox::X_D3DTRANSFORMSTATETYPE, CONST xbox::X_D3DMATRIX*)                                             );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_SetTransform_0__LTCG_eax1_edx2,           ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_MultiplyTransform,                        (xbox::X_D3DTRANSFORMSTATETYPE, CONST xbox::X_D3DMATRIX*)                                             );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_MultiplyTransform_0__LTCG_ebx1_eax2,      ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3D_DestroyResource,                                (xbox::X_D3DResource*)                                                                                );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3D_DestroyResource_0__LTCG_edi1,                   (xbox::void_xt)                                                                                       );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     Direct3D_CreateDevice,                              (xbox::uint_xt, D3DDEVTYPE, HWND, xbox::dword_xt, xbox::X_D3DPRESENT_PARAMETERS*, xbox::X_D3DDevice**));  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     Direct3D_CreateDevice_16__LTCG_eax4_ebx6,           (xbox::uint_xt, D3DDEVTYPE, HWND, xbox::X_D3DPRESENT_PARAMETERS*)                                     );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     Direct3D_CreateDevice_16__LTCG_eax4_ecx6,           (xbox::uint_xt, D3DDEVTYPE, HWND, xbox::X_D3DPRESENT_PARAMETERS*)                                     );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     Direct3D_CreateDevice_4__LTCG_eax1_ecx3,            (xbox::X_D3DPRESENT_PARAMETERS*)                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     Lock2DSurface,                                      (xbox::X_D3DPixelContainer*, xbox::X_D3DCUBEMAP_FACES, xbox::uint_xt, xbox::X_D3DLOCKED_RECT*, RECT*, xbox::dword_xt) );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     Lock2DSurface_16__LTCG_esi4_eax5,                   (xbox::X_D3DPixelContainer*, xbox::X_D3DCUBEMAP_FACES, xbox::uint_xt, xbox::dword_xt)                         );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     Lock3DSurface,                                      (xbox::X_D3DPixelContainer*, xbox::uint_xt, xbox::X_D3DLOCKED_BOX*, xbox::X_D3DBOX*, xbox::dword_xt)                  );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     Lock3DSurface_16__LTCG_eax4,                        (xbox::X_D3DPixelContainer*, xbox::uint_xt, xbox::X_D3DLOCKED_BOX*, xbox::dword_xt)                           );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3D_CommonSetRenderTarget,                          (xbox::X_D3DSurface*, xbox::X_D3DSurface*, void*)                                                     );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_BeginPushBuffer,                          (xbox::X_D3DPushBuffer*)                                                                              );  \
   	XB_MACRO(xbox::hresult_xt,    WINAPI,     D3DDevice_EndPushBuffer,                            (xbox::void_xt)                                                                                       );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_DrawVertices,                             (xbox::X_D3DPRIMITIVETYPE, xbox::uint_xt, xbox::uint_xt)                                              );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_DrawIndexedVertices,                      (xbox::X_D3DPRIMITIVETYPE, xbox::uint_xt, CONST PWORD)                                                );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     D3DDevice_DrawVerticesUP,                           (xbox::X_D3DPRIMITIVETYPE, xbox::uint_xt, CONST PVOID, xbox::uint_xt)                                 );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     CDevice_SetStateVB,                                 (xbox::ulong_xt)                                                                                      );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     CDevice_SetStateVB_8,                               (xbox::addr_xt, xbox::ulong_xt)                                                                       );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     CDevice_SetStateUP,                                 ()                                                                                                    );  \
   	XB_MACRO(xbox::void_xt,       WINAPI,     CDevice_SetStateUP_4,                               (xbox::addr_xt)                                                                                       );  \

// Generate extern declarations for trampoline function pointers.
// Each expands to a typedef + extern variable declaration.
#define XB_trampoline_extern(ret, conv, func, arguments) \
   	typedef ret(conv * XB_TRAMPOLINE_##func##_t)arguments; \
   	extern XB_TRAMPOLINE_##func##_t XB_TRAMPOLINE_##func

XB_TRAMPOLINES(XB_trampoline_extern);
#undef XB_trampoline_extern


#endif // RENDERGLOBALS_H
