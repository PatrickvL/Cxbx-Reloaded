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

// Direct3D_CreateDevice_Start and Direct3D_CreateDevice_End are no longer used.
// Host D3D11 device creation and host-state init moved to CxbxInitHostD3DDevice().
// Xbox state is maintained by native Xbox CreateDevice running unpatched.
// Xbox backbuffer/depth surfaces are resolved by the PGRAPH RT path on first use.

// Called by D3D11_draw_inline_elements (XbPushBuffer.cpp)
void CxbxDrawIndexed(CxbxDrawContext &DrawContext)
{
	assert(DrawContext.dwStartVertex == 0);
	assert(DrawContext.pXboxIndexData != nullptr);
	assert(DrawContext.dwVertexCount > 0); // TODO : If this fails, make responsible callers do an early-exit

	CxbxD3D11IABypassDraw(DrawContext);
	g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(DrawContext.XboxPrimitiveType, DrawContext.dwVertexCount);
}

// Drawing function specifically for rendering Xbox draw calls supplying a 'User Pointer'.
// Called by D3D11_draw_inline_array (XbPushBuffer.cpp)
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

