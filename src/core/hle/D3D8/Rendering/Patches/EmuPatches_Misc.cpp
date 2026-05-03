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

// Variables only used in EmuPatches_Misc.cpp
static std::stack<ID3D11Query*> g_HostQueryVisibilityTests;
static std::map<int, ID3D11Query*> g_HostVisibilityTestMap;

xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_BeginVisibilityTest)()
{
	LOG_FUNC();

	if (g_bEnableHostQueryVisibilityTest) {
		D3D11_QUERY_DESC QueryDesc;
		QueryDesc.Query = D3D11_QUERY_OCCLUSION;
		QueryDesc.MiscFlags = 0;
		// Create a D3D occlusion query to handle "visibility test" with
		ID3D11Query* pHostQueryVisibilityTest = nullptr;
		HRESULT hRet = g_pD3DDevice->CreateQuery(&QueryDesc, &pHostQueryVisibilityTest);
		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateQuery (visibility test)");
		if (pHostQueryVisibilityTest != nullptr) {
			CxbxQueryIssueBegin(pHostQueryVisibilityTest);
			{
				g_HostQueryVisibilityTests.push(pHostQueryVisibilityTest);
			}

			pHostQueryVisibilityTest = nullptr;
		}
	}

	return S_OK;
}

// LTCG specific D3DDevice_EndVisibilityTest function...
// This uses a custom calling convention where parameter is passed in EAX
// Test-case: Test Drive: Eve of Destruction
__declspec(naked) xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_EndVisibilityTest_0__LTCG_eax1)
(
)
{
   	dword_xt Index;
	xbox::hresult_xt result;
   	__asm {
   	   	LTCG_PROLOGUE
   	   	mov  Index, eax
   	}

   	result = EMUPATCH(D3DDevice_EndVisibilityTest)(Index);

   	__asm {
   	   	mov  eax, result
   	   	LTCG_EPILOGUE
   	   	ret
   	}
}

// ******************************************************************
// * patch: D3DDevice_EndVisibilityTest
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_EndVisibilityTest)
(
   	dword_xt                       Index
)
{
	LOG_FUNC_ONE_ARG(Index);

	if (g_bEnableHostQueryVisibilityTest) {
		// Check that the dedicated storage for the given Index isn't in use
		if (g_HostVisibilityTestMap[Index] != nullptr) {
			return E_OUTOFMEMORY;
		}

		if (g_HostQueryVisibilityTests.empty()) {
			return 2088; // visibility test incomplete (a prior BeginVisibilityTest call is needed)
		}

		ID3D11Query* pHostQueryVisibilityTest = g_HostQueryVisibilityTests.top();
		g_HostQueryVisibilityTests.pop();
		assert(pHostQueryVisibilityTest != nullptr);

		CxbxQueryIssueEnd(pHostQueryVisibilityTest);
		{
			// Associate the result of this call with the given Index
			g_HostVisibilityTestMap[Index] = pHostQueryVisibilityTest;
		}
	}

   	return S_OK;
}

// ******************************************************************
// * patch: D3DDevice_GetVisibilityTestResult
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_GetVisibilityTestResult)
(
   	dword_xt                       Index,
   	uint_xt                       *pResult,
   	ulonglong_xt                  *pTimeStamp
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Index)
		LOG_FUNC_ARG(pResult)
		LOG_FUNC_ARG(pTimeStamp)
		LOG_FUNC_END;

	if (g_bEnableHostQueryVisibilityTest) {
		ID3D11Query* pHostQueryVisibilityTest = g_HostVisibilityTestMap[Index];
		if (pHostQueryVisibilityTest == nullptr) {
			return E_OUTOFMEMORY;
		}

		// In order to prevent an endless loop if the D3D device becomes lost, we pass
		// the D3DGETDATA_FLUSH flag. This tells GetData to return D3DERR_DEVICELOST if
		// such a situation occurs, and break out of the loop as a result.
		// Note: By Cxbx's design, we cannot do drawing within this while loop in order
		// to further prevent any other endless loop situations.
		UINT64 occlusionData = 0;
		HRESULT hRet;
		while ((hRet = CxbxQueryGetData(pHostQueryVisibilityTest, &occlusionData, sizeof(occlusionData), 0)) == S_FALSE) {
			SwitchToThread(); // Yield CPU while waiting for GPU query result
		}
		if (FAILED(hRet)) {
			EmuLog(LOG_LEVEL::WARNING, "GetVisibilityTestResult: query failed (0x%08X)", hRet);
		}
		if (pResult != xbox::zeroptr)
			*pResult = (uint_xt)occlusionData;

		g_HostVisibilityTestMap[Index] = nullptr;
		pHostQueryVisibilityTest->Release();
	} else {
		// Fallback to old faked result when there's no host occlusion query :
		if (pResult != xbox::zeroptr) {
			*pResult = 640 * 480; // TODO : Use actual backbuffer dimensions
		}
	}

	if (pTimeStamp != xbox::zeroptr) {
		LOG_TEST_CASE("requested value for pTimeStamp");
		*pTimeStamp = sizeof(DWORD); // TODO : This should be an incrementing GPU-memory based DWORD-aligned memory address
	}

   	return S_OK;
}

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

