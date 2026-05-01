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

// D3DDevice_SetBackBufferScale — disabled.
// Host-only concept (upscale factor); Xbox native code doesn't need this.
// Patch disabled in Patches.cpp.

// ******************************************************************
// * patch: D3DDevice_SetGammaRamp
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetGammaRamp)
(
   	dword_xt                   dwFlags,
   	CONST X_D3DGAMMARAMP      *pRamp
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(dwFlags)
		LOG_FUNC_ARG(pRamp)
		LOG_FUNC_END;

	// Use IDXGIOutput::SetGammaControl for D3D11
	if (g_pSwapChain) {
		IDXGIOutput* pOutput = nullptr;
		if (SUCCEEDED(g_pSwapChain->GetContainingOutput(&pOutput))) {
			DXGI_GAMMA_CONTROL gammaControl = {};
			gammaControl.Scale = { 1.0f, 1.0f, 1.0f };
			gammaControl.Offset = { 0.0f, 0.0f, 0.0f };
			for (int v = 0; v < 256; v++) {
				float idx = v / 255.0f * 1024.0f;
				int i = static_cast<int>(idx);
				if (i > 1024) i = 1024;
				gammaControl.GammaCurve[i] = {
					pRamp->red[v] / 255.0f,
					pRamp->green[v] / 255.0f,
					pRamp->blue[v] / 255.0f
				};
			}
			// Interpolate any gaps in the 1025-entry curve
			for (int i = 1; i < 1025; i++) {
				if (gammaControl.GammaCurve[i].Red == 0.0f &&
					gammaControl.GammaCurve[i].Green == 0.0f &&
					gammaControl.GammaCurve[i].Blue == 0.0f && i < 1024) {
					gammaControl.GammaCurve[i] = gammaControl.GammaCurve[i - 1];
				}
			}
			pOutput->SetGammaControl(&gammaControl);
			pOutput->Release();
		}
	}
}

// ******************************************************************
// * patch: D3DDevice_GetGammaRamp
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_GetGammaRamp)
(
   	X_D3DGAMMARAMP     *pRamp
)
{
	LOG_FUNC_ONE_ARG(pRamp);

   	// Use IDXGIOutput::GetGammaControl to retrieve the current gamma ramp
   	bool gotGamma = false;
   	if (g_pSwapChain) {
   	   	IDXGIOutput* pOutput = nullptr;
   	   	if (SUCCEEDED(g_pSwapChain->GetContainingOutput(&pOutput))) {
   	   	   	DXGI_GAMMA_CONTROL gammaControl = {};
   	   	   	if (SUCCEEDED(pOutput->GetGammaControl(&gammaControl))) {
   	   	   	   	for (int v = 0; v < 256; v++) {
   	   	   	   	   	int i = static_cast<int>(v / 255.0f * 1024.0f);
   	   	   	   	   	if (i > 1024) i = 1024;
   	   	   	   	   	pRamp->red[v]   = static_cast<BYTE>(gammaControl.GammaCurve[i].Red * 255.0f);
   	   	   	   	   	pRamp->green[v] = static_cast<BYTE>(gammaControl.GammaCurve[i].Green * 255.0f);
   	   	   	   	   	pRamp->blue[v]  = static_cast<BYTE>(gammaControl.GammaCurve[i].Blue * 255.0f);
   	   	   	   	}
   	   	   	   	gotGamma = true;
   	   	   	}
   	   	   	pOutput->Release();
   	   	}
   	}
   	if (!gotGamma) {
   	   	// Fallback: return a linear ramp
   	   	for (int v = 0; v < 256; v++) {
   	   	   	pRamp->red[v]   = (BYTE)v;
   	   	   	pRamp->green[v] = (BYTE)v;
   	   	   	pRamp->blue[v]  = (BYTE)v;
   	   	}
   	}
}


xbox::X_D3DSurface* CxbxrImpl_GetBackBuffer2
(
   	xbox::int_xt BackBuffer
)
{
	xbox::X_D3DSurface* pXboxBackBuffer = nullptr;

	// Rather than create a new surface, we should forward to the Xbox version of GetBackBuffer,
	// This gives us the correct Xbox surface to update.
	// We get signatures for both backbuffer functions as it changed in later XDKs

	// This also updates the reference count, so we don't need to do this ourselves
	if (XB_TRMP(D3DDevice_GetBackBuffer) != nullptr) {
		XB_TRMP(D3DDevice_GetBackBuffer)(BackBuffer, xbox::X_D3DBACKBUFFER_TYPE_MONO, &pXboxBackBuffer);
	}
	else if (XB_TRMP(D3DDevice_GetBackBuffer_8__LTCG_eax1) != nullptr) {
		__asm {
			lea  eax, pXboxBackBuffer
			push eax
			push D3DBACKBUFFER_TYPE_MONO
			mov  eax, BackBuffer
			call XB_TRMP(D3DDevice_GetBackBuffer_8__LTCG_eax1)
		}
	}
	else if (XB_TRMP(D3DDevice_GetBackBuffer2) != nullptr) {
		pXboxBackBuffer = XB_TRMP(D3DDevice_GetBackBuffer2)(BackBuffer);
	}
	else {
		__asm {
			mov  eax, BackBuffer
			call XB_TRMP(D3DDevice_GetBackBuffer2_0__LTCG_eax1)
			mov  pXboxBackBuffer, eax
		}
	}

	// Now pXboxBackbuffer points to the requested Xbox backbuffer
	if (pXboxBackBuffer == nullptr) {
		CxbxrAbort("D3DDevice_GetBackBuffer2: Could not get Xbox backbuffer");
	}

	return pXboxBackBuffer;
}

// D3DDevice_SetViewport — disabled (trampoline-only after CxbxImpl_SetViewport removal).
// Patch disabled in Patches.cpp — Xbox code runs unpatched.

// CxbxImpl_SetViewport — removed.
// The Xbox trampoline writes viewport state to the NV2A push buffer.
// PGRAPH VPSCL/VPOFF registers are the authority; g_Xbox_Viewport was
// only written here and had no render-thread readers.

// D3DDevice_SetShaderConstantMode_0__LTCG_eax1 — disabled.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetTexture — disabled.
// Xbox native SetTexture pushes NV097_SET_TEXTURE_OFFSET to the push buffer.
// Host texture lookup uses PGRAPH TEXOFFSET registers set by the push buffer.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SwitchTexture — disabled.
// Xbox native SwitchTexture updates texture offset mid-draw via push buffer.
// Host texture lookup uses PGRAPH TEXOFFSET registers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetTransform — disabled.
// Transform state now sourced from PGRAPH XFCTX registers (MMAT0/CMAT/TnMAT).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
// Test case: 25 to Life (MultiplyTransform should call SetTransform internally)

// D3DDevice_MultiplyTransform — disabled.
// Xbox native code calls SetTransform internally which pushes NV2A transform methods.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetStreamSource (all LTCG variants) — disabled.
// Xbox native SetStreamSource writes NV097_SET_VERTEX_DATA_ARRAY_OFFSET/FORMAT
// to the push buffer. Host vertex binding reads PGRAPH state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

