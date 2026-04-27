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
#if 0 // Unused :
   	// Set g_Xbox_D3DDevice to point to the Xbox D3D Device
   	auto it = g_SymbolAddresses.find("D3D_g_pDevice");
   	if (it != g_SymbolAddresses.end()) {
   	   	g_Xbox_D3DDevice = (DWORD*)it->second;
   	}
   	else {
   	   	EmuLog(LOG_LEVEL::ERROR2, "D3D_g_pDevice was not found!");
   	}
#endif

   	UpdateHostBackBufferDesc();
   	SetAspectRatioScale(pPresentationParameters);

   	// If the Xbox version of CreateDevice didn't call SetRenderTarget, we must derive the default backbuffer ourselves
   	// This works because CreateDevice always sets the current render target to the Xbox Backbuffer
   	// In later XDKs, it does this inline rather than by calling D3DDevice_SetRenderTarget
   	// meaning our patch doesn't always get called in these cases.
   	// We fix the situation by calling the Xbox GetRenderTarget function, which immediately after CreateDevice
   	// WILL always return the Backbuffer!
   	// Test Case: Shin Megami Tensei: Nine
   	if (g_pXbox_BackBufferSurface == xbox::zeroptr && g_pXbox_DefaultDepthStencilSurface == xbox::zeroptr) {
   	   	// First, log the test case
   	   	LOG_TEST_CASE("Xbox CreateDevice did not call SetRenderTarget");
   	}

   	if (g_pXbox_BackBufferSurface == xbox::zeroptr) {
   	   	if (XB_TRMP(D3DDevice_GetRenderTarget)) {
   	   	   	XB_TRMP(D3DDevice_GetRenderTarget)(&g_pXbox_BackBufferSurface);
   	   	}
   	   	else if (XB_TRMP(D3DDevice_GetRenderTarget2)) {
   	   	   	g_pXbox_BackBufferSurface = XB_TRMP(D3DDevice_GetRenderTarget2)();
   	   	}

   	   	// At this point, g_pXbox_BackBufferSurface should now point to a valid render target
   	   	// if it still doesn't, we cannot continue without crashing at draw time
   	   	if (g_pXbox_BackBufferSurface == xbox::zeroptr) {
   	   	   	CxbxrAbort("Unable to determine default Xbox backbuffer");
   	   	}

   	   	// Set the backbuffer as the initial rendertarget
   	   	CxbxImpl_SetRenderTarget(g_pXbox_BackBufferSurface, xbox::zeroptr);
   	}

   	// Now do the same, but for the default depth stencil surface
   	if (g_pXbox_DefaultDepthStencilSurface == xbox::zeroptr) {
   	   	if (XB_TRMP(D3DDevice_GetDepthStencilSurface)) {
   	   	   	XB_TRMP(D3DDevice_GetDepthStencilSurface)(&g_pXbox_DefaultDepthStencilSurface);
   	   	}
   	   	else if (XB_TRMP(D3DDevice_GetDepthStencilSurface2)) {
   	   	   	g_pXbox_DefaultDepthStencilSurface = XB_TRMP(D3DDevice_GetDepthStencilSurface2)();
   	   	}

   	   	// At this point, g_pXbox_DefaultDepthStencilSurface should now point to a valid depth stencil
   	   	// If it doesn't, just log and carry on: Unlike RenderTarget, this situation is not fatal
   	   	if (g_pXbox_DefaultDepthStencilSurface == xbox::zeroptr) {
   	   	   	LOG_TEST_CASE("Unable to determine default Xbox depth stencil");
   	   	} else {
   	   	   	// Update only the depth stencil
   	   	   	CxbxImpl_SetRenderTarget(xbox::zeroptr, g_pXbox_DefaultDepthStencilSurface);
   	   	}
   	}
}

// Called by D3DDevice_DrawIndexedVertices and EmuExecutePushBufferRaw (twice)
void CxbxDrawIndexed(CxbxDrawContext &DrawContext)
{
	assert(DrawContext.dwStartVertex == 0);
	assert(DrawContext.pXboxIndexData != nullptr);
	assert(DrawContext.dwVertexCount > 0); // TODO : If this fails, make responsible callers do an early-exit

	CxbxD3D11IABypassDraw(DrawContext);
	g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(DrawContext.XboxPrimitiveType, DrawContext.dwVertexCount);
}

// TODO : Move to own file
// Drawing function specifically for rendering Xbox draw calls supplying a 'User Pointer'.
// Called by D3DDevice_DrawVerticesUP and EmuExecutePushBufferRaw
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

