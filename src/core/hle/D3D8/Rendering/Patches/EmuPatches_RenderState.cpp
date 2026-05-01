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
#include "../IndexBufferConvert.h"

// D3DDevice_SetStreamSource_8__LTCG_edx1 — disabled.
// Xbox native SetStreamSource writes NV097_SET_VERTEX_DATA_ARRAY_OFFSET/FORMAT
// to the push buffer. Host vertex binding reads PGRAPH state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetStreamSource — disabled.
// Same as above. Patch disabled in Patches.cpp.

// D3DDevice_SetVertexShader, D3DDevice_SetVertexShader_0__LTCG_ebx1 — disabled.
// Xbox native SetVertexShader calls LoadVertexShader + SelectVertexShader which
// push NV097_SET_TRANSFORM_PROGRAM and NV097_SET_TRANSFORM_PROGRAM_START.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetRenderTarget, D3DDevice_SetRenderTarget_0__LTCG_ecx1_eax2,
// D3D_CommonSetRenderTarget, D3DDevice_SetRenderTargetFast — disabled.
// Xbox native SetRenderTarget writes NV097_SET_SURFACE_COLOR_OFFSET and related
// PGRAPH surface registers. Host RT is now created from PGRAPH surface state by
// CxbxD3D11UpdateRenderTargetFromPGRAPH. Backbuffer detection handled in CreateDevice.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// CxbxImpl_SetRenderTarget — still used by CreateDevice path (HostDraw.cpp)
// to register backbuffer/depth stencil surfaces for PGRAPH-based RT resolution.
void CxbxImpl_SetRenderTarget
(
	xbox::X_D3DSurface    *pRenderTarget,
	xbox::X_D3DSurface    *pNewZStencil
)
{
	LOG_INIT;

	// In Xbox titles, CreateDevice calls SetRenderTarget for the back buffer
	// We can use this to determine the Xbox backbuffer surface for later use!
	if (g_pXbox_BackBufferSurface == xbox::zeroptr) {
		g_pXbox_BackBufferSurface = pRenderTarget;
		// TODO : Some titles might render to another backbuffer later on,
		// if that happens, we might need to skip the first one or two calls?
	}

	// In Xbox titles, CreateDevice calls SetRenderTarget (our caller) for the depth stencil
	// We can use this to determine the Xbox depth stencil surface for later use!
	if (g_pXbox_DefaultDepthStencilSurface == xbox::zeroptr) {
		g_pXbox_DefaultDepthStencilSurface = pNewZStencil;
		// TODO : Some titles might set another depth stencil later on,
		// if that happens, we might need to skip the first one or two calls?
	}

	// The current render target is only replaced if it's passed in here non-null
	if (pRenderTarget != xbox::zeroptr) {
		g_pXbox_RenderTarget = pRenderTarget;
		CxbxRegisterSurfaceByDataAddr(pRenderTarget->Data, pRenderTarget);
	}
	else if (g_pXbox_RenderTarget == xbox::zeroptr) {
		// No explicit render target and none was set yet.
		// Try the backbuffer (set by the first SetRenderTarget from CreateDevice).
		// When SetRenderTarget patches are disabled and trampolines aren't found,
		// g_pXbox_BackBufferSurface may also be null — that's fine because
		// CxbxD3D11UpdateRenderTargetFromPGRAPH creates host resources directly.
		if (g_pXbox_BackBufferSurface != xbox::zeroptr) {
			g_pXbox_RenderTarget = g_pXbox_BackBufferSurface;
			CxbxRegisterSurfaceByDataAddr(g_pXbox_BackBufferSurface->Data, g_pXbox_BackBufferSurface);
		}
	}

	// The currenct depth stencil is always replaced by whats passed in here (even a null)
	g_pXbox_DepthStencil = pNewZStencil;
	if (pNewZStencil != xbox::zeroptr)
		CxbxRegisterSurfaceByDataAddr(pNewZStencil->Data, pNewZStencil);

	// Host D3D11 render target and depth stencil creation + binding
	// is now fully handled by CxbxD3D11UpdateRenderTargetFromPGRAPH,
	// which creates host resources directly from PGRAPH surface state.
}

// D3DDevice_SetPalette, D3DDevice_SetPalette_4__LTCG_eax1 — disabled.
// Xbox native SetPalette pushes NV097_SET_TEXTURE_PALETTE.
// Host palette lookup uses PGRAPH texture palette registers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_DeleteVertexShader_0__LTCG_eax1 — disabled.
// Vertex shader cache cleanup. Host shader cache uses PGRAPH program store.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetScreenSpaceOffset — disabled.
// g_Xbox_ScreenSpaceOffset was only read by dead CxbxSetVertexShaderPassthroughProgram.
// Patch disabled in Patches.cpp — Xbox code handles it natively.

// D3D_LazySetPointParams — disabled (unimplemented stub, LOG_UNIMPLEMENTED).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetRenderState_Simple — disabled.
// Xbox code already writes D3D__RenderState[] and pushes NV2A methods.
// This patch only mirrored to XboxRenderStates which is redundant.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetTransform / D3DDevice_SetTransform_0__LTCG_eax1_edx2 — disabled.
// Transform state now sourced from PGRAPH XFCTX registers (MMAT0/CMAT/TnMAT).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
// Test case: 25 to Life (MultiplyTransform should call SetTransform internally)
