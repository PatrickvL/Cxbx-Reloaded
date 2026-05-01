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

// D3DDevice_BeginPush / D3DDevice_EndPush — disabled.
// Xbox native BeginPush returns pointer into GPU ring buffer (CDevice+0x00 = pPut).
// PFIFO processes commands from ring buffer when DMA_PUT is updated by EndPush.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetIndices, D3DDevice_SetIndices_4__LTCG_ebx1 — disabled.
// Trampoline-only after DrawIndexedVertices was disabled (g_Xbox_BaseVertexIndex had no readers).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// * patch: Direct3D_CreateDevice
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(Direct3D_CreateDevice)
(
   	uint_xt                     Adapter,
   	X_D3DDEVTYPE                DeviceType,
	X_HWND                      hFocusWindow,
   	dword_xt                    BehaviorFlags,
   	X_D3DPRESENT_PARAMETERS    *pPresentationParameters,
   	X_D3DDevice               **ppReturnedDeviceInterface
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Adapter)
		LOG_FUNC_ARG(DeviceType)
		LOG_FUNC_ARG(hFocusWindow)
		LOG_FUNC_ARG(BehaviorFlags)
		LOG_FUNC_ARG(pPresentationParameters)
		LOG_FUNC_ARG_OUT(ppReturnedDeviceInterface)
		LOG_FUNC_END;

	Direct3D_CreateDevice_Start(pPresentationParameters);

	// HACK: Disable Software VertexProcessing (Fixes CreateDevice failure in Chihiro titles)
	if (BehaviorFlags & D3DCREATE_SOFTWARE_VERTEXPROCESSING) {
		BehaviorFlags &= ~D3DCREATE_SOFTWARE_VERTEXPROCESSING;
		BehaviorFlags |= D3DCREATE_HARDWARE_VERTEXPROCESSING;
	}

	// Only then call Xbox CreateDevice function
	hresult_xt hRet = XB_TRMP(Direct3D_CreateDevice)(Adapter, DeviceType, hFocusWindow, BehaviorFlags, pPresentationParameters, ppReturnedDeviceInterface);

	Direct3D_CreateDevice_End(pPresentationParameters);

	return hRet;
}

// Overload for logging
static void Direct3D_CreateDevice_4__LTCG_eax1_ecx3
(
   	xbox::dword_xt                BehaviorFlags,
   	xbox::X_D3DPRESENT_PARAMETERS *pPresentationParameters,
   	xbox::X_D3DDevice            **ppReturnedDeviceInterface
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(BehaviorFlags)
		LOG_FUNC_ARG(pPresentationParameters)
		LOG_FUNC_ARG_OUT(ppReturnedDeviceInterface)
		LOG_FUNC_END;
}

// LTCG specific Direct3D_CreateDevice function...
// This uses a custom calling with parameters passed in EAX, ECX and the stack
// Test-case: Ninja Gaiden, Halo 2
__declspec(naked) xbox::hresult_xt WINAPI xbox::EMUPATCH(Direct3D_CreateDevice_4__LTCG_eax1_ecx3)
(
   	X_D3DPRESENT_PARAMETERS     *pPresentationParameters
)
{
   	DWORD BehaviorFlags;
   	xbox::X_D3DDevice **ppReturnedDeviceInterface;
   	__asm {
		LTCG_PROLOGUE
		mov  BehaviorFlags, eax
		mov  ppReturnedDeviceInterface, ecx
	}

	// Log
	Direct3D_CreateDevice_4__LTCG_eax1_ecx3(BehaviorFlags, pPresentationParameters, ppReturnedDeviceInterface);

	Direct3D_CreateDevice_Start(pPresentationParameters);

	// Only then call Xbox CreateDevice function
	hresult_xt hRet;
	__asm {
		push pPresentationParameters
		mov  eax, BehaviorFlags
		mov  ecx, ppReturnedDeviceInterface
		call XB_TRMP(Direct3D_CreateDevice_4__LTCG_eax1_ecx3)
		mov  hRet, eax
	}

	Direct3D_CreateDevice_End(pPresentationParameters);

	__asm {
		mov  eax, hRet
		LTCG_EPILOGUE
		ret  4
	}
}

// Overload for logging
static void Direct3D_CreateDevice_16__LTCG_eax4_ecx6
(
   	xbox::uint_xt                 Adapter,
	xbox::X_D3DDEVTYPE            DeviceType,
	xbox::X_HWND                  hFocusWindow,
   	xbox::dword_xt                BehaviorFlags,
   	xbox::X_D3DPRESENT_PARAMETERS *pPresentationParameters,
   	xbox::X_D3DDevice            **ppReturnedDeviceInterface
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Adapter)
		LOG_FUNC_ARG(DeviceType)
		LOG_FUNC_ARG(hFocusWindow)
		LOG_FUNC_ARG(BehaviorFlags)
		LOG_FUNC_ARG(pPresentationParameters)
		LOG_FUNC_ARG_OUT(ppReturnedDeviceInterface)
		LOG_FUNC_END;
}

// LTCG specific Direct3D_CreateDevice function...
// This uses a custom calling convention passing parameters on stack, in EAX and ECX
// Test-case: Aggressive Inline, Midtown Madness 3
__declspec(naked) xbox::hresult_xt WINAPI xbox::EMUPATCH(Direct3D_CreateDevice_16__LTCG_eax4_ecx6)
(
   	uint_xt                     Adapter,
   	X_D3DDEVTYPE                DeviceType,
	X_HWND                      hFocusWindow,
   	X_D3DPRESENT_PARAMETERS     *pPresentationParameters
)
{
	dword_xt BehaviorFlags;
	xbox::X_D3DDevice **ppReturnedDeviceInterface;
	__asm {
		LTCG_PROLOGUE
		mov  BehaviorFlags, eax
		mov  ppReturnedDeviceInterface, ecx
	}

	// Log
	Direct3D_CreateDevice_16__LTCG_eax4_ecx6(Adapter, DeviceType, hFocusWindow, BehaviorFlags, pPresentationParameters, ppReturnedDeviceInterface);

	Direct3D_CreateDevice_Start(pPresentationParameters);

	// Only then call Xbox CreateDevice function
	hresult_xt hRet;
	__asm {
		push pPresentationParameters
		push hFocusWindow
		push DeviceType
		push Adapter
		mov  eax, BehaviorFlags
		mov  ecx, ppReturnedDeviceInterface
		call XB_TRMP(Direct3D_CreateDevice_16__LTCG_eax4_ecx6)
		mov  hRet, eax
	}

	Direct3D_CreateDevice_End(pPresentationParameters);

	__asm {
		mov  eax, hRet
		LTCG_EPILOGUE
		ret  10h
	}
}

// Overload for logging
static void Direct3D_CreateDevice_16__LTCG_eax4_ebx6
(
   	xbox::uint_xt                 Adapter,
   	xbox::X_D3DDEVTYPE            DeviceType,
	xbox::X_HWND                  hFocusWindow,
   	xbox::dword_xt                BehaviorFlags,
   	xbox::X_D3DPRESENT_PARAMETERS *pPresentationParameters,
   	xbox::X_D3DDevice            **ppReturnedDeviceInterface
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Adapter)
		LOG_FUNC_ARG(DeviceType)
		LOG_FUNC_ARG(hFocusWindow)
		LOG_FUNC_ARG(BehaviorFlags)
		LOG_FUNC_ARG(pPresentationParameters)
		LOG_FUNC_ARG_OUT(ppReturnedDeviceInterface)
		LOG_FUNC_END;
}

// LTCG specific Direct3D_CreateDevice function...
// This uses a custom calling convention passing parameters on stack, in EAX and EBX
// Test-case: NASCAR Heat 2002
__declspec(naked) xbox::hresult_xt WINAPI xbox::EMUPATCH(Direct3D_CreateDevice_16__LTCG_eax4_ebx6)
(
   	uint_xt                     Adapter,
   	X_D3DDEVTYPE                DeviceType,
	X_HWND                      hFocusWindow,
   	X_D3DPRESENT_PARAMETERS     *pPresentationParameters
)
{
	dword_xt BehaviorFlags;
	xbox::X_D3DDevice **ppReturnedDeviceInterface;
	__asm {
		LTCG_PROLOGUE
		mov  BehaviorFlags, eax
		mov  ppReturnedDeviceInterface, ebx
	}

	// Log
	Direct3D_CreateDevice_16__LTCG_eax4_ebx6(Adapter, DeviceType, hFocusWindow, BehaviorFlags, pPresentationParameters, ppReturnedDeviceInterface);

	Direct3D_CreateDevice_Start(pPresentationParameters);

	// Only then call Xbox CreateDevice function
	hresult_xt hRet;
	__asm {
		push pPresentationParameters
		push hFocusWindow
		push DeviceType
		push Adapter
		mov  eax, BehaviorFlags
		mov  ebx, ppReturnedDeviceInterface
		call XB_TRMP(Direct3D_CreateDevice_16__LTCG_eax4_ebx6)
		mov  hRet, eax
	}

	Direct3D_CreateDevice_End(pPresentationParameters);

	__asm {
		mov  eax, hRet
		LTCG_EPILOGUE
		ret  10h
	}
}

// Overload for logging
static void D3DDevice_SetIndices_4__LTCG_ebx1
(
   	xbox::X_D3DIndexBuffer   *pIndexData,
   	xbox::uint_xt             BaseVertexIndex
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(pIndexData)
   	   	LOG_FUNC_ARG(BaseVertexIndex)
   	   	LOG_FUNC_END;
}

// ******************************************************************
// * patch: D3DDevice_SetIndices_4__LTCG_ebx1
// LTCG specific D3DDevice_SetIndices function...
// D3DDevice_SetIndices_4__LTCG_ebx1 — removed (see comment above SetIndices).

// D3DDevice_Reset — disabled.
// Xbox D3DDevice_Reset only does Xbox D3D resource cleanup.
// Host rendering driven by NV2A state needs no reset — host resources
// are managed by PGRAPH RT cache and invalidated by dirty-page detection.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.


// ******************************************************************
// * patch: D3DDevice_GetDisplayFieldStatus
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_GetDisplayFieldStatus)(X_D3DFIELD_STATUS *pFieldStatus)
{
	// NOTE: This can be unpatched only when NV2A does it's own VBlank and HLE _Swap function is unpatched


	LOG_FUNC_ONE_ARG(pFieldStatus);

	// Read AV Flags to determine if Progressive or Interlaced
	// The xbox does this by reading from pDevice->m_Miniport.m_CurrentAvInfo
	// but we don't have an OOVPA for that. Instead, we call the Xbox implementation of 
	// D3DDevice_GetDisplayMode and read the result

	X_D3DDISPLAYMODE displayMode;

	// If we can find the Xbox version of GetDisplayMode, use the real data returned, otherwise
	// use a sensible default
	if (XB_TRMP(D3DDevice_GetDisplayMode) != nullptr) {
		XB_TRMP(D3DDevice_GetDisplayMode)(&displayMode);
	} else {
		// We don't show a warning because doing so pollutes the kernel debug log as this function gets called every frame
		displayMode.Flags = X_D3DPRESENTFLAG_INTERLACED;
	}
	
	// Set the VBlank count
	pFieldStatus->VBlankCount = g_Xbox_VBlankData.VBlank;

	// If we are interlaced, return the current field, otherwise, return progressive scan
	if (displayMode.Flags & X_D3DPRESENTFLAG_INTERLACED) {
		pFieldStatus->Field = (g_Xbox_VBlankData.VBlank % 2 == 0) ? X_D3DFIELD_ODD : X_D3DFIELD_EVEN;
	} else {
		pFieldStatus->Field = X_D3DFIELD_PROGRESSIVE;
	}
}

// ******************************************************************
// * patch: D3DDevice_BeginVisibilityTest
// ******************************************************************
