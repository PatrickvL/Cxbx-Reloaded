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
#include "core/hle/D3D8/Rendering/Backend/Backend_D3D11_Profiler.h"

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
		CXBX_PROFILE_SCOPE(PROF_BLOCKONTTIME);
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
