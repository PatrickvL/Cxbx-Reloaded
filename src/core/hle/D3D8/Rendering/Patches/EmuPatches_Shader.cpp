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

// Mirror a texture's VRAM offset to the PGRAPH TEXOFFSET register so that
// CxbxUpdateHostTextures() can resolve the Xbox texture even when the
// pushbuffer hasn't been processed yet.
static void CxbxMirrorTexOffsetToPGRAPH(DWORD Stage, xbox::addr_xt dataAddr)
{
	if (Stage < xbox::X_D3DTS_STAGECOUNT) {
		PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;
		pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + Stage * 4)] = dataAddr;
	}
}


// D3DDevice_LoadVertexShader — disabled (trampoline-only after CxbxImpl removal).
// Patch disabled in Patches.cpp — Xbox code runs unpatched.

// Overload for logging
static void D3DDevice_SelectVertexShader_0__LTCG_eax1_ebx2
(
   	xbox::dword_xt                 Handle,
   	xbox::dword_xt                 Address
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(Handle)
   	   	LOG_FUNC_ARG(Address)
   	   	LOG_FUNC_END;
}

// LTCG specific D3DDevice_SelectVertexShader function...
// This uses a custom calling convention where parameter is passed in EAX, EBX
// Test-case: Star Wars - Battlefront
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SelectVertexShader_0__LTCG_eax1_ebx2)
(
)
{
   	dword_xt Handle;
   	dword_xt Address;
   	__asm {
   	   	LTCG_PROLOGUE
   	   	mov  Handle, eax
   	   	mov  Address, ebx
   	}

   	// Log
   	D3DDevice_SelectVertexShader_0__LTCG_eax1_ebx2(Handle, Address);

   	CxbxImpl_SelectVertexShader(Handle, Address);

   	__asm {
   	   	LTCG_EPILOGUE
   	   	ret
   	}
}

// Overload for logging
static void D3DDevice_SelectVertexShader_4__LTCG_eax1
(
   	xbox::dword_xt                 Handle,
   	xbox::dword_xt                 Address
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(Handle)
   	   	LOG_FUNC_ARG(Address)
   	   	LOG_FUNC_END;
}

// LTCG specific D3DDevice_SelectVertexShader function...
// This uses a custom calling convention where parameter is passed in EAX
// Test-case: Aggressive Inline
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SelectVertexShader_4__LTCG_eax1)
(
   	dword_xt                       Address
)
{
   	dword_xt Handle;
   	__asm {
   	   	LTCG_PROLOGUE
   	   	mov  Handle, eax
   	}

   	// Log
   	D3DDevice_SelectVertexShader_4__LTCG_eax1(Handle, Address);

   	CxbxImpl_SelectVertexShader(Handle, Address);

   	__asm {
   	   	LTCG_EPILOGUE
   	   	ret  4
   	}
}

// ******************************************************************
// * patch: D3DDevice_SelectVertexShader
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SelectVertexShader)
(
   	dword_xt                       Handle,
   	dword_xt                       Address
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(Handle)
   	   	LOG_FUNC_ARG(Address)
   	   	LOG_FUNC_END;

	// Call the Xbox trampoline so the NV2A push buffer gets the program start update.
	XB_TRMP(D3DDevice_SelectVertexShader)(Handle, Address);

   	CxbxImpl_SelectVertexShader(Handle, Address);
}

// ******************************************************************
// D3DDevice_SetShaderConstantMode — disabled.
// g_Xbox_VertexShaderConstantMode has no render-thread readers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetVertexShaderConstant variants (8 patches) — disabled.
// Xbox native SetVertexShaderConstant generates NV097_SET_TRANSFORM_CONSTANT
// push buffer commands that flow through PFIFO → PGRAPH. No HLE interception needed.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// Overload for logging
static void D3DDevice_SetTexture_4__LTCG_eax2
(
   	xbox::dword_xt           Stage,
   	xbox::X_D3DBaseTexture  *pTexture
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(Stage)
   	   	LOG_FUNC_ARG(pTexture)
   	   	LOG_FUNC_END;
}

// LTCG specific D3DDevice_SetTexture function...
// This uses a custom calling convention where pTexture is passed in EAX
// Test-case: NASCAR Heat 2002
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetTexture_4__LTCG_eax2)
(
   	dword_xt           Stage
)
{
   	X_D3DBaseTexture *pTexture;
   	__asm {
   	   	LTCG_PROLOGUE
   	   	mov  pTexture, eax
   	}

   	// Log
   	D3DDevice_SetTexture_4__LTCG_eax2(Stage, pTexture);

   	// Call the Xbox implementation of this function, to properly handle reference counting for us
   	__asm {
   	   	mov eax, pTexture
   	   	push Stage
   	   	call XB_TRMP(D3DDevice_SetTexture_4__LTCG_eax2)
   	}

   	g_pXbox_SetTexture[Stage] = pTexture;

   	// Register in VRAM-offset → texture side-map for PGRAPH TEXOFFSET lookup
   	if (pTexture != xbox::zeroptr && pTexture->Data != xbox::zero) {
   	   	CxbxRegisterTextureByDataAddr(pTexture->Data, pTexture);
   	   	CxbxMirrorTexOffsetToPGRAPH(Stage, pTexture->Data);
   	} else {
   	   	CxbxMirrorTexOffsetToPGRAPH(Stage, 0);
   	}

   	__asm {
   	   	LTCG_EPILOGUE
   	   	ret  4
   	}
}

// Overload for logging
static void D3DDevice_SetTexture_4__LTCG_eax1
(
   	xbox::dword_xt           Stage,
   	xbox::X_D3DBaseTexture  *pTexture
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(Stage)
   	   	LOG_FUNC_ARG(pTexture)
   	   	LOG_FUNC_END;
}

// LTCG specific D3DDevice_SetTexture function...
// This uses a custom calling convention where Stage is passed in EAX
// Test-case: Metal Wolf Chaos
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetTexture_4__LTCG_eax1)
(
   	X_D3DBaseTexture  *pTexture
)
{
   	dword_xt Stage;
   	__asm {
   	   	LTCG_PROLOGUE
   	   	mov  Stage, eax
   	}

   	// Log
	D3DDevice_SetTexture_4__LTCG_eax1(Stage, pTexture);

   	// Call the Xbox implementation of this function, to properly handle reference counting for us
   	__asm {
   	   	mov eax, Stage
   	   	push pTexture
   	   	call XB_TRMP(D3DDevice_SetTexture_4__LTCG_eax1)
   	}

   	g_pXbox_SetTexture[Stage] = pTexture;

   	// Register in VRAM-offset → texture side-map for PGRAPH TEXOFFSET lookup
   	if (pTexture != xbox::zeroptr && pTexture->Data != xbox::zero) {
   	   	CxbxRegisterTextureByDataAddr(pTexture->Data, pTexture);
   	   	CxbxMirrorTexOffsetToPGRAPH(Stage, pTexture->Data);
   	} else {
   	   	CxbxMirrorTexOffsetToPGRAPH(Stage, 0);
   	}

   	__asm {
   	   	LTCG_EPILOGUE
   	   	ret  4
   	}
}

// ******************************************************************
// * patch: D3DDevice_SetPixelShader
// D3DDevice_SetPixelShader, D3DDevice_SetPixelShader_0__LTCG_eax1 — disabled.
// These only called CxbxImpl_SetPixelShader (after the trampoline) to write
// g_pXbox_PixelShader, which was only read by the COMBINECTL==0 HLE bridge
// fallback, now removed. Bodies moved to Direct3D9.cpp.unused-patches.

// D3DDevice_DrawVertices_4__LTCG_ecx2_eax3, D3DDevice_DrawVertices_8__LTCG_eax3 — disabled.
// LTCG variants of DrawVertices; patches disabled in Patches.cpp.

// ******************************************************************
// * patch: D3DDevice_DeleteVertexShader
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_DeleteVertexShader)
(
	dword_xt Handle
)
{
	LOG_FUNC_ONE_ARG(Handle);

	CxbxImpl_DeleteVertexShader(Handle);

	// When deleting, call trampoline *after* our implementation,
	// so that we can still access it's fields before it gets deleted!
	XB_TRMP(D3DDevice_DeleteVertexShader)(Handle);
}



// ******************************************************************
// D3DDevice_GetShaderConstantMode — disabled.
// g_Xbox_VertexShaderConstantMode has no render-thread readers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetVertexShader — disabled.
// Getter reads g_Xbox_VertexShader_Handle; Xbox native reads from device struct.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetVertexShaderConstant — disabled.
// Getter reads HLE VS constant shadow; Xbox native reads from device constant table.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// * patch: D3DDevice_SetVertexShaderInput
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexShaderInput)
(
   	dword_xt              Handle,
   	uint_xt               StreamCount,
   	X_STREAMINPUT     *pStreamInputs
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Handle)
		LOG_FUNC_ARG(StreamCount)
		LOG_FUNC_ARG(pStreamInputs)
		LOG_FUNC_END;

	// When this API is in effect, VertexBuffers as set by Xbox SetStreamSource are disregarded,
	// instead, the pStreamInputs[].VertexBuffer streams are used.

	// If Handle is NULL, all VertexShader input state is cleared (after which the VertexBuffers as set by SetStreamSource are used once again).

	// Otherwise, Handle is the address of an Xbox VertexShader struct, or-ed with 1 (X_D3DFVF_RESERVED0)
	// The given pStreamInputs are stored in a global array, and the NV2A is programmed to read
	// each vertex attribute (as defined in the given VertexShader.VertexAttribute.Slots[]) to read
	// the attribute data from the pStreamInputs[slot].VertexBuffer + pStreamInputs[slot].Offset + VertexShader.VertexAttribute.Slots[slot].Offset

	/* LOG_TEST_CASE("SetVertexShaderInput");
	/* Test-cases :
		PushBuffer XDK sample
		Halo 2-3ebe4439.ini:D3DDevice_SetVertexShaderInput = 0x3f7440
		Kung Fu Chaos-d9ab292c.ini:D3DDevice_SetVertexShaderInput = 0x2bc0e0
		NBA LIVE 2005-71d4eeb1.ini:D3DDevice_SetVertexShaderInput = 0x5cf810
		NBA LIVE 2005-71d4eeb1.ini:D3DDevice_SetVertexShaderInputDirect = 0x5ceba0
		Prince of Persia WW-4ccf7369.ini:D3DDevice_SetVertexShaderInput = 0x494830
		Prince of Persia WW-4ccf7369.ini:D3DDevice_SetVertexShaderInputDirect = 0x494280
		Spyro A Hero's Tail-b18e00e5.ini:D3DDevice_SetVertexShaderInput = 0x286cf0
		Spyro A Hero's Tail-b18e00e5.ini:D3DDevice_SetVertexShaderInputDirect = 0x286760
	*/

	CxbxImpl_SetVertexShaderInput(Handle, StreamCount, pStreamInputs);

	// Call trampoline
	if (XB_TRMP(D3DDevice_SetVertexShaderInput))
		XB_TRMP(D3DDevice_SetVertexShaderInput)(Handle, StreamCount, pStreamInputs);
}

// ******************************************************************
// * patch: D3DDevice_RunVertexStateShader
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_RunVertexStateShader)
(
   	dword_xt Address,
   	CONST float_xt *pData
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Address)
		LOG_FUNC_ARG(pData)
		LOG_FUNC_END;

	CxbxrImpl_RunVertexStateShader(Address, pData);
}

// ******************************************************************
// * patch: D3DDevice_RunVertexStateShader_4__LTCG_esi2
// ******************************************************************
// Overload for logging
static void D3DDevice_RunVertexStateShader_4__LTCG_esi2
(
	xbox::dword_xt Address,
	CONST xbox::float_xt* pData
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Address)
		LOG_FUNC_ARG(pData)
		LOG_FUNC_END;
}

// This uses a custom calling convention where parameter is passed in ESI
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_RunVertexStateShader_4__LTCG_esi2)
(
   	dword_xt Address
)
{
	float_xt *pData;
	__asm {
		LTCG_PROLOGUE
		mov  pData, esi
	}

	// Log
	D3DDevice_RunVertexStateShader_4__LTCG_esi2(Address, pData);

	CxbxrImpl_RunVertexStateShader(Address, pData);

	__asm {
		LTCG_EPILOGUE
		ret  4
	}
}

// ******************************************************************
// D3DDevice_SetDepthClipPlanes — disabled (all cases are TODO stubs).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
