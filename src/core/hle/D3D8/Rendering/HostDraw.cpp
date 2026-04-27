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
#include "EmuD3D8_common.h"
#include "Backend\Backend_D3D11.h"

void Direct3D_CreateDevice_Start
(
	const xbox::X_D3DPRESENT_PARAMETERS     *pPresentationParameters
)
{
   	CxbxVertexShaderSetFlags();

   	if (!XboxRenderStates.Init()) {
   	   	CxbxrAbort("Failed to init XboxRenderStates");
   	}

   	if (!XboxTextureStates.Init(&XboxRenderStates)) {
   	   	CxbxrAbort("Failed to init XboxTextureStates");
   	}

	SetXboxMultiSampleType(pPresentationParameters->MultiSampleType);

	// create default device *before* calling Xbox Direct3D_CreateDevice trampoline
	// to avoid hitting EMUPATCH'es that need a valid g_pD3DDevice

	if (g_pD3DDevice != nullptr) { // Check to make sure device is null, otherwise no need to create it
		return;
	}

	CreateDefaultDevice(pPresentationParameters);
}

void Direct3D_CreateDevice_End
(
	const xbox::X_D3DPRESENT_PARAMETERS     *pPresentationParameters
)
{

   	UpdateHostBackBufferDesc();
   	SetAspectRatioScale(pPresentationParameters);

	// Reset PGRAPH surface tracking so the first RT bound becomes the tracked backbuffer
	CxbxResetPgraphSurfaceTracking();

   	// Try to determine the Xbox backbuffer and depth stencil surfaces.
	// These are used for side-map registration (helps PGRAPH RT path reuse Xbox resource metadata).
	// With SetRenderTarget patches disabled, the initial SetRenderTarget from CreateDevice won't
	// be intercepted, so we fetch the surfaces via GetRenderTarget/GetDepthStencilSurface trampolines.
	// This is optional — the PGRAPH RT path creates host resources directly from NV2A state if needed.
   	if (g_pXbox_BackBufferSurface == xbox::zeroptr) {
   	   	if (XB_TRMP(D3DDevice_GetRenderTarget)) {
   	   	   	XB_TRMP(D3DDevice_GetRenderTarget)(&g_pXbox_BackBufferSurface);
   	   	}
   	   	else if (XB_TRMP(D3DDevice_GetRenderTarget2)) {
   	   	   	g_pXbox_BackBufferSurface = XB_TRMP(D3DDevice_GetRenderTarget2)();
   	   	}

   	   	if (g_pXbox_BackBufferSurface != xbox::zeroptr) {
   	   	   	CxbxImpl_SetRenderTarget(g_pXbox_BackBufferSurface, xbox::zeroptr);
   	   	} else {
			EmuLog(LOG_LEVEL::WARNING, "Could not determine Xbox backbuffer — PGRAPH path will create host RT directly");
   	   	}
   	}

   	if (g_pXbox_DefaultDepthStencilSurface == xbox::zeroptr) {
   	   	if (XB_TRMP(D3DDevice_GetDepthStencilSurface)) {
   	   	   	XB_TRMP(D3DDevice_GetDepthStencilSurface)(&g_pXbox_DefaultDepthStencilSurface);
   	   	}
   	   	else if (XB_TRMP(D3DDevice_GetDepthStencilSurface2)) {
   	   	   	g_pXbox_DefaultDepthStencilSurface = XB_TRMP(D3DDevice_GetDepthStencilSurface2)();
   	   	}

   	   	if (g_pXbox_DefaultDepthStencilSurface != xbox::zeroptr) {
   	   	   	CxbxImpl_SetRenderTarget(xbox::zeroptr, g_pXbox_DefaultDepthStencilSurface);
   	   	}
   	}
}

// Called by HLE_draw_inline_elements (XbPushBuffer.cpp)
void CxbxDrawIndexed(CxbxDrawContext &DrawContext)
{
	assert(DrawContext.dwStartVertex == 0);
	assert(DrawContext.pXboxIndexData != nullptr);
	assert(DrawContext.dwVertexCount > 0); // TODO : If this fails, make responsible callers do an early-exit

	CxbxD3D11IABypassDraw(DrawContext);
	g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(DrawContext.XboxPrimitiveType, DrawContext.dwVertexCount);
}

// Drawing function specifically for rendering Xbox draw calls supplying a 'User Pointer'.
// Called by HLE_draw_inline_array (XbPushBuffer.cpp)
void CxbxDrawPrimitiveUP(CxbxDrawContext &DrawContext)
{
	assert(DrawContext.dwStartVertex == 0);
	assert(DrawContext.pXboxVertexStreamZeroData != xbox::zeroptr);
	assert(DrawContext.uiXboxVertexStreamZeroStride > 0);
	assert(DrawContext.dwBaseVertexIndex == 0); // No IndexBase under Draw*UP

	CxbxD3D11IABypassDraw(DrawContext);
	g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(DrawContext.XboxPrimitiveType, DrawContext.dwVertexCount);
}

ID3D11Resource* CxbxConvertXboxSurfaceToHostTexture(xbox::X_D3DBaseTexture* pBaseTexture)
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	ID3D11Texture2D* pHostSurface = GetHostSurface(pBaseTexture);
	if (!pHostSurface) {
   	   	LOG_TEST_CASE("Failed to get host surface");
		return nullptr;
	}

	// For D3D11, the surface IS already the texture (ID3D11Texture2D)
	// Just add a reference and return it
	pHostSurface->AddRef();
	return pHostSurface;
}

