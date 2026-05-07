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

// D3DDevice_GetBackBuffer2, D3DDevice_GetBackBuffer2_0__LTCG_eax1,
// D3DDevice_GetBackBuffer, D3DDevice_GetBackBuffer_8__LTCG_eax1 — disabled.
// Xbox GetBackBuffer2 returns a pointer to the internal backbuffer surface
// structure. Since we no longer intercept SetRenderTarget, the backbuffer
// surface is already known from CreateDevice (g_pXbox_BackBufferSurface).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
// Test-case: NBA 2K2 (LTCG_eax1)

// D3DDevice_Clear — disabled.
// Xbox Clear uses the NV2A 2D engine (SOLID_RECTANGLE class 0x5E or
// GDI_RECTANGLE_TEXT class 0x4A) to fill color/Z/stencil buffers.
// Our NV2A 2D engine implementation now handles these blit classes.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.


// ******************************************************************
// * patch: D3DDevice_CopyRects
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_CopyRects)
(
   	X_D3DSurface*  pSourceSurface,
   	CONST X_RECT*  pSourceRectsArray,
   	uint_xt        cRects,
   	X_D3DSurface*  pDestinationSurface,
   	CONST X_POINT* pDestPointsArray
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(pSourceSurface);
   	   	LOG_FUNC_ARG(pSourceRectsArray);
   	   	LOG_FUNC_ARG(cRects);
   	   	LOG_FUNC_ARG(pDestinationSurface);
   	   	LOG_FUNC_ARG(pDestPointsArray);
   	LOG_FUNC_END;

   	// Copy on the host GPU side for surfaces that have host representations
   	auto pHostSourceSurface = GetHostSurface(pSourceSurface);
   	auto pHostDestSurface = GetHostSurface(pDestinationSurface);

   	if (pHostSourceSurface == nullptr || pHostDestSurface == nullptr) {
   	   	// Test Case: DOA2 attempts to copy from an index buffer resource type
   	   	// TODO: What should we do here?
   	   	LOG_TEST_CASE("D3DDevice-CopyRects: Failed to fetch host surfaces");
   	   	return;
   	}

	D3D11_TEXTURE2D_DESC hostSourceDesc, hostDestDesc;
   	pHostSourceSurface->GetDesc(&hostSourceDesc);
   	pHostDestSurface->GetDesc(&hostDestDesc);

   	// If the source is a render-target and the destination is not, we need force it to be re-created as one
   	// This is because StrechRects cannot copy from a Render-Target to a Non-Render Target
   	// Test Case: Crash Bandicoot: Wrath of Cortex attemps to copy the render-target to a texture
   	// This fixes an issue on the pause screen where the screenshot of the current scene was not displayed correctly
   	if ((hostSourceDesc.Usage & D3DUSAGE_RENDERTARGET) != 0 && (hostDestDesc.Usage & D3DUSAGE_RENDERTARGET) == 0) {
   	   	pHostDestSurface = GetHostSurface(pDestinationSurface, D3DUSAGE_RENDERTARGET);
   	   	pHostDestSurface->GetDesc(&hostDestDesc);
   	}

   	// If no rectangles were given, default to 1 (entire surface)
   	if (cRects == 0) {
   	   	cRects = 1;
   	}

   	// Get Xbox surface dimensions
   	// Host resources may be scaled so we'll account for that later
   	auto xboxSourceWidth = GetPixelContainerWidth(pSourceSurface);
   	auto xboxSourceHeight = GetPixelContainerHeight(pSourceSurface);
   	auto xboxDestWidth = GetPixelContainerWidth(pDestinationSurface);
   	auto xboxDestHeight = GetPixelContainerHeight(pDestinationSurface);

   	for (UINT i = 0; i < cRects; i++) {
   	   	RECT SourceRect, DestRect;

   	   	if (pSourceRectsArray != nullptr) {
   	   	   	SourceRect.left = pSourceRectsArray[i].left;
   	   	   	SourceRect.right = pSourceRectsArray[i].right;
   	   	   	SourceRect.top = pSourceRectsArray[i].top;
   	   	   	SourceRect.bottom = pSourceRectsArray[i].bottom;
   	   	} else {
   	   	   	SourceRect.left = 0;
   	   	   	SourceRect.right = xboxSourceWidth;
   	   	   	SourceRect.top = 0;
   	   	   	SourceRect.bottom = xboxSourceHeight;
   	   	}

   	   	if (pDestPointsArray != nullptr) {
   	   	   	DestRect.left = pDestPointsArray[i].x;
   	   	   	DestRect.right = DestRect.left + (SourceRect.right - SourceRect.left);
   	   	   	DestRect.top = pDestPointsArray[i].y;
   	   	   	DestRect.bottom = DestRect.top + (SourceRect.bottom - SourceRect.top);
   	   	} else if (pSourceRectsArray) {
   	   	   	DestRect = SourceRect;
   	   	} else {
   	   	   	DestRect.left = 0;
   	   	   	DestRect.right = xboxDestWidth;
   	   	   	DestRect.top = 0;
   	   	   	DestRect.bottom = xboxDestHeight;
   	   	}

   	   	// Scale the source and destination rects
   	   	auto sourceScaleX = (uint32_t)hostSourceDesc.Width / xboxSourceWidth;
   	   	auto sourceScaleY = (uint32_t)hostSourceDesc.Height / xboxSourceHeight;

   	   	SourceRect.left *= sourceScaleX;
   	   	SourceRect.right *= sourceScaleX;
   	   	SourceRect.top *= sourceScaleY;
   	   	SourceRect.bottom *= sourceScaleY;

   	   	auto destScaleX = (uint32_t) hostDestDesc.Width / xboxDestWidth;
   	   	auto destScaleY = (uint32_t) hostDestDesc.Height / xboxDestHeight;

   	   	DestRect.left *= destScaleX;
   	   	DestRect.right *= destScaleX;
   	   	DestRect.top *= destScaleY;
   	   	DestRect.bottom *= destScaleY;

   	   	HRESULT hRet;
   	   	hRet = CxbxBltSurface(pHostSourceSurface, &SourceRect, pHostDestSurface, &DestRect, D3DTEXF_LINEAR);
   	   	if (FAILED(hRet)) {
   	   	   	LOG_TEST_CASE("D3DDevice_CopyRects: Failed to copy surface");
   	   	}
   	}
}

// CXBX_SWAP_PRESENT_FORWARD is now defined in RenderGlobals.h

// D3DDevice_Present — disabled.
// Native Swap pushes NV097_FLIP_INCREMENT_WRITE + NV097_FLIP_STALL to push buffer.

// D3DDevice_Swap / D3DDevice_Swap_0__LTCG_eax1 — disabled.
// Native Swap pushes NV097_FLIP_STALL → pgraph_flip_stall → D3D11_flip_stall.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// Lock2DSurface, Lock2DSurface_16__LTCG_esi4_eax5,
// Lock3DSurface, Lock3DSurface_16__LTCG_eax4 — disabled.
// All were trampoline-only (LOG_FUNC + XB_TRMP) with no side effects.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_PersistDisplay — disabled.
// Was an incomplete stub that just called the Xbox trampoline.
// Xbox code saves a framebuffer copy to contiguous memory for dashboard hand-off.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

