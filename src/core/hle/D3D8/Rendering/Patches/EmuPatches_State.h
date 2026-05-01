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
#ifndef EMUPATCHES_STATE_H
#define EMUPATCHES_STATE_H

#include "core\hle\XAPI\Xapi.h" // For EMUPATCH
#include "core\hle\D3D8\XbD3D8Types.h"

namespace xbox {

// ******************************************************************
// * patch: D3DDevice_SetBackBufferScale — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetBackBufferScale)(float_xt x, float_xt y);

// D3DDevice_SetGammaRamp — disabled (PRMDIO VGA DAC palette emulation handles gamma natively)
// D3DDevice_GetGammaRamp — disabled (PRMDIO VGA DAC palette emulation handles gamma natively)

// ******************************************************************
// * patch: D3DDevice_SetViewport
// ******************************************************************
xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetViewport)
(
   	X_D3DVIEWPORT8 *pViewport
);

xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetShaderConstantMode_0__LTCG_eax1)();

// ******************************************************************
// * patch: D3DDevice_SetTexture — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetTexture)
//(
//	dword_xt           Stage,
//	X_D3DBaseTexture  *pTexture
//);

// ******************************************************************
// * patch: D3DDevice_SwitchTexture — Disabled
// ******************************************************************
//xbox::void_xt __fastcall EMUPATCH(D3DDevice_SwitchTexture)
//(
//	dword_xt           Method,
//	dword_xt           Data,
//	dword_xt           Format
//);

// ******************************************************************
// * patch: IDirect3DResource8_AddRef
// ******************************************************************
xbox::ulong_xt WINAPI EMUPATCH(D3DResource_AddRef)
(
   	X_D3DResource      *pThis
);

// ******************************************************************
// * patch: D3DDevice_SetStreamSource — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetStreamSource)
//(
//	uint_xt                StreamNumber,
//	X_D3DVertexBuffer  *pStreamData,
//	uint_xt                Stride
//);

// ******************************************************************
// * patch: D3DDevice_SetVertexShader — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetVertexShader)
//(
//	dword_xt            Handle
//);

//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetVertexShader_0__LTCG_ebx1)();

// ******************************************************************
// * patch: D3DDevice_SetRenderTarget — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetRenderTarget)
//(
//	X_D3DSurface    *pRenderTarget,
//	X_D3DSurface    *pNewZStencil
//);

// ******************************************************************
// * patch: D3DDevice_SetRenderTarget_0__LTCG_ecx1_eax2 — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetRenderTarget_0__LTCG_ecx1_eax2)();

// ******************************************************************
// * patch: D3D_CommonSetRenderTarget — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3D_CommonSetRenderTarget)
//(
//	X_D3DSurface    *pRenderTarget,
//	X_D3DSurface    *pNewZStencil,
//	void            *unknown
//);

// ******************************************************************
// * patch: D3DDevice_SetPalette — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetPalette)
//(
//	dword_xt         Stage,
//	X_D3DPalette *pPalette
//);

//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetPalette_4__LTCG_eax1)
//(
//	X_D3DPalette *pPalette
//);

//xbox::void_xt WINAPI EMUPATCH(D3DDevice_DeleteVertexShader_0__LTCG_eax1)();

// ******************************************************************
// * patch: D3DDevice_SetScreenSpaceOffset
// ******************************************************************
xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetScreenSpaceOffset)
(
   	float_xt x,
   	float_xt y
);

// ******************************************************************
// * patch: D3DDevice_SetRenderTargetFast — Disabled
// ******************************************************************
//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetRenderTargetFast)
//(
//	X_D3DSurface	*pRenderTarget,
//	X_D3DSurface	*pNewZStencil,
//	dword_xt			Flags
//);

// ******************************************************************
// * patch: D3D::LazySetPointParams
// ******************************************************************
void WINAPI EMUPATCH(D3D_LazySetPointParams)( void* Device );

//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetStreamSource_0__LTCG_eax1_edi2_ebx3)();

//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetStreamSource_4__LTCG_eax1_ebx2)
//(
//	uint_xt                Stride
//);

//xbox::void_xt WINAPI EMUPATCH(D3DDevice_SetStreamSource_8__LTCG_eax1)
//(
//	X_D3DVertexBuffer  *pStreamData,
//	uint_xt             Stride
//);

//xbox::void_xt __fastcall EMUPATCH(D3DDevice_SetStreamSource_8__LTCG_edx1)
//(
//	void*,
//	uint_xt                StreamNumber,
//	X_D3DVertexBuffer  *pStreamData,
//	uint_xt                Stride
//);

} // namespace xbox

#endif // EMUPATCHES_STATE_H
