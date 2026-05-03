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
#include "devices/Xbox.h"              // For extern NV2ADevice* g_NV2A
#include "devices/video/nv2a.h"        // For pfifo_flush_to_pgraph

// D3DDevice_BeginVisibilityTest / EndVisibilityTest / GetVisibilityTestResult — disabled.
// Visibility tests now handled natively via NV2A PGRAPH: Xbox code pushes
// NV097_CLEAR_REPORT_VALUE, NV097_SET_ZPASS_PIXEL_COUNT_ENABLE, and NV097_GET_REPORT
// through PFIFO. PGRAPH wraps D3D11 draws with occlusion queries (D3D11_zpass_begin/end
// in XbPushBuffer.cpp) and accumulates z-passing pixel counts into pg->zpass_pixel_count_result.
// Xbox native D3DDevice_GetVisibilityTestResult reads the report DMA memory directly.

// D3DDevice_EnableOverlay — disabled.
// Native Xbox code programs NV_PVIDEO_STOP/BUFFER; D3D11_flip_stall reads PVIDEO state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_UpdateOverlay — disabled.
// Native Xbox code programs PVIDEO registers (OFFSET, SIZE_IN, FORMAT, POINT_OUT, SIZE_OUT);
// D3D11_flip_stall composites overlay from VRAM using PVIDEO register state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetOverlayUpdateStatus — disabled.
// Hardcoded TRUE stub; Xbox native overlay check is correct.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_InsertFence — disabled.
// Fake 0x8000BEEF stub; Xbox native fence via NV2A reference counter.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_IsFencePending — disabled.
// Hardcoded FALSE stub; Xbox native fence check via NV2A.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_BlockOnFence — disabled.
// Empty stub; Xbox native fence wait via NV2A.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_BlockUntilVerticalBlank — disabled.
// Native Xbox code waits on m_VerticalBlankEvent KEVENT inside the D3D device struct.
// With Direct3D_CreateDevice unpatched, D3D_g_pDevice is valid and VBlank IRQ signals it.

// D3DResource_BlockUntilNotBusy — disabled.
// Empty stub; Xbox native code polls resource state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_InsertCallback — disabled.
// Native InsertCallback pushes NV097_NO_OPERATION(param) to the push buffer.
// PGRAPH raises INTR_ERROR → miniport ISR reads TRAPPED_DATA_LOW → dispatches callback.

// D3DDevice_GetProjectionViewportMatrix — disabled.
// Xbox native code reads projection from D3DDevice struct and builds viewport matrix.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetModelView — disabled.
// SetModelView state now sourced from PGRAPH XFCTX registers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_FlushVertexCache — disabled.
// Unimplemented stub with no side effects.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetModelView — disabled.
// Xbox native code reads WorldView from D3DDevice struct.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3D_SetCommonDebugRegisters — disabled.
// Empty LOG_UNIMPLEMENTED stub; Xbox native code writes harmless debug regs.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_IsBusy — disabled.
// Hardcoded FALSE stub; Xbox native version checks NV_PGRAPH_STATUS.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// * patch: D3D_BlockOnTime
// ******************************************************************
// The Xbox D3D runtime enables the DMA pusher (PUSH0_ACCESS=1,
// DMA_PUSH_ACCESS=1) during CreateDevice.  When the access flags are
// set the PFIFO read fast-path does NOT fake GET=PUT, so the native
// polling loop would spin forever waiting for the pusher to advance
// DMA_GET — which never happens because the Xbox ring buffer is not
// mapped for real FIFO-mode processing in HLE mode.
// We therefore intercept the call and drain any pending PGRAPH
// commands via pfifo_flush_to_pgraph, then return, bypassing the
// native spin loop entirely.
void WINAPI xbox::EMUPATCH(D3D_BlockOnTime)(dword_xt Time, int MakeSpace)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Time)
		LOG_FUNC_ARG(MakeSpace)
		LOG_FUNC_END;

	// Drain pending PFIFO commands so the DMA pusher advances GET,
	// freeing ring buffer space for the caller.
	if (g_NV2A) {
		pfifo_flush_to_pgraph(g_NV2A->GetDeviceState());
	}
}

// ******************************************************************
// * patch: D3D_BlockOnTime_4__LTCG_eax1
// ******************************************************************
__declspec(naked) void WINAPI xbox::EMUPATCH(D3D_BlockOnTime_4__LTCG_eax1)(int MakeSpace)
{
	xbox::dword_xt Time;
	__asm {
		LTCG_PROLOGUE
		mov  Time, eax
	}
	EMUPATCH(D3D_BlockOnTime)(Time, MakeSpace);
	__asm {
		LTCG_EPILOGUE
		ret  4
	}
}

// D3D_DestroyResource — disabled.
// Host resources are NV2A-derived (PGRAPH RT cache, texture cache keyed by VRAM).
// Dirty page tracking invalidates stale host resources.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

