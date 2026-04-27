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
#include "../Backend\Backend_D3D11.h"

// D3DDevice_Begin, D3DDevice_SetVertexData2f, D3DDevice_SetVertexData2s,
// D3DDevice_SetVertexData4f_16__LTCG_edi1, D3DDevice_SetVertexData4f,
// D3DDevice_SetVertexData4ub, D3DDevice_SetVertexData4s,
// D3DDevice_SetVertexDataColor, D3DDevice_End — disabled.
// Xbox native Begin/End/SetVertexData pushes NV2A methods (NV097_SET_BEGIN_END,
// NV097_ARRAY_ELEMENT16 etc.) through the push buffer → PFIFO → PGRAPH.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_BeginPushBuffer, D3DDevice_BeginPushBuffer_0__LTCG_edi1,
// D3DDevice_EndPushBuffer — disabled.
// These only called the trampoline and toggled g_bRecordingPushBuffer
// which had zero readers. Xbox code runs unpatched.
// Patch disabled in Patches.cpp.

// ******************************************************************
// * patch: D3DDevice_RunPushBuffer
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_RunPushBuffer)
(
   	X_D3DPushBuffer       *pPushBuffer,
   	X_D3DFixup            *pFixup
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pPushBuffer)
		LOG_FUNC_ARG(pFixup)
		LOG_FUNC_END;

	EmuExecutePushBuffer(pPushBuffer, pFixup);    
}

// ******************************************************************
// * patch: D3DDevice_RunPushBuffer_4__LTCG_eax2
// ******************************************************************
// Overload for logging
static void D3DDevice_RunPushBuffer_4__LTCG_eax2
(
   	xbox::X_D3DPushBuffer       *pPushBuffer,
	xbox::X_D3DFixup            *pFixup
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pPushBuffer)
		LOG_FUNC_ARG(pFixup)
		LOG_FUNC_END;
}

// This uses a custom calling convention where parameter is passed in EAX
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_RunPushBuffer_4__LTCG_eax2)
(
   	X_D3DPushBuffer *pPushBuffer
)
{
	X_D3DFixup* pFixup;
	__asm {
		LTCG_PROLOGUE
		mov  pFixup, eax
	}

	// Log
	D3DDevice_RunPushBuffer_4__LTCG_eax2(pPushBuffer, pFixup);

	EmuExecutePushBuffer(pPushBuffer, pFixup);

	__asm {
		LTCG_EPILOGUE
		ret  4
	}
}

// D3DDevice_DrawVertices, DrawVerticesUP, DrawVerticesUP_12__LTCG_ebx3,
// DrawIndexedVertices, DrawIndexedVerticesUP — disabled.
// Xbox native draw calls push NV2A methods (NV097_SET_BEGIN_END, NV097_DRAW_ARRAYS,
// NV097_ARRAY_ELEMENT16, etc.) through the push buffer → PFIFO → PGRAPH.
// Draws are now initiated from NV097 writes in the LLE path.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// CDevice_SetStateVB, CDevice_SetStateVB_8, CDevice_SetStateUP,
// CDevice_SetStateUP_4, CDevice_SetStateUP_0__LTCG_esi1 — disabled.
// These were unimplemented stubs (LOG_UNIMPLEMENTED). Xbox native code pushes
// all necessary NV2A state through the push buffer before draw calls.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// D3DDevice_SetStipple — disabled (unimplemented stub, LOG_IGNORED).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// * patch: D3DDevice_SetSwapCallback
// ******************************************************************
void WINAPI xbox::EMUPATCH(D3DDevice_SetSwapCallback)
(
	X_D3DSWAPCALLBACK		pCallback
)
{
	LOG_FUNC_ONE_ARG(pCallback);

   	g_pXbox_SwapCallback = pCallback;
}

// ******************************************************************
// D3DDevice_PrimeVertexCache — disabled (unimplemented stub, LOG_UNIMPLEMENTED).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// * patch: D3DDevice_DrawRectPatch
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawRectPatch)
(
	uint_xt					Handle,
	CONST float_xt				*pNumSegs,
	CONST X_D3DRECTPATCH_INFO *pRectPatchInfo
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Handle)
		LOG_FUNC_ARG(pNumSegs)
		LOG_FUNC_ARG(pRectPatchInfo)
		LOG_FUNC_END;

	CxbxUpdateNativeD3DResources();

	// D3D11 has no DrawRectPatch - use CPU tessellation
	HRESULT hRet = CxbxDrawRectPatchD3D11(Handle, pNumSegs, pRectPatchInfo);
	DEBUG_D3DRESULT(hRet, "CxbxDrawRectPatchD3D11");

	return hRet;
}

// ******************************************************************
// * patch: D3DDevice_DrawTriPatch
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawTriPatch)
(
	uint_xt					Handle,
	CONST float_xt				*pNumSegs,
	CONST X_D3DTRIPATCH_INFO* pTriPatchInfo
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Handle)
		LOG_FUNC_ARG(pNumSegs)
		LOG_FUNC_ARG(pTriPatchInfo)
		LOG_FUNC_END;

	CxbxUpdateNativeD3DResources();

	// D3D11 has no DrawTriPatch - use CPU tessellation
	HRESULT hRet = CxbxDrawTriPatchD3D11(Handle, pNumSegs, pTriPatchInfo);
	DEBUG_D3DRESULT(hRet, "CxbxDrawTriPatchD3D11");

	return hRet;
}
