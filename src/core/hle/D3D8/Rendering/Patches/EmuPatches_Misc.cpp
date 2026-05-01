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
static DWORD g_OverlaySwap = 0; // Set in D3DDevice_UpdateOverlay
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

// ******************************************************************
// * patch: D3DDevice_EnableOverlay
// ******************************************************************
static void CxbxrImpl_EnableOverlay()
{
	// The Xbox D3DDevice_EnableOverlay call merely resets the active
	// NV2A overlay state, it doesn't actually enable or disable anything.
	// Thus, we should just reset our overlay state here too. A title will
	// show new overlay data via D3DDevice_UpdateOverlay (see below).
	g_OverlayProxy = {};
}

xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_EnableOverlay)
(
   	bool_xt Enable
)
{
	LOG_FUNC_ONE_ARG(Enable);

	CxbxrImpl_EnableOverlay();
}

// ******************************************************************
// * patch: D3DDevice_EnableOverlay_0__LTCG
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_EnableOverlay_0__LTCG)()
{
	LOG_FUNC();

	CxbxrImpl_EnableOverlay();
}

static void CxbxrImpl_UpdateOverlay
(
	xbox::X_D3DSurface *pSurface,
	CONST xbox::X_RECT *SrcRect,
	CONST xbox::X_RECT *DstRect,
	xbox::bool_xt       EnableColorKey,
	xbox::X_D3DCOLOR    ColorKey
)
{
	using namespace xbox;

	// Reset and remember the overlay arguments, so our D3DDevice_Swap patch
	// can correctly show this overlay surface data.
	g_OverlayProxy = {};
	if (pSurface) {
		g_OverlayProxy.Surface = *pSurface;
		if (SrcRect)
			g_OverlayProxy.SrcRect = *SrcRect;

		if (DstRect)
			g_OverlayProxy.DstRect = *DstRect;

		g_OverlayProxy.EnableColorKey = EnableColorKey;
		g_OverlayProxy.ColorKey = ColorKey;
		// Update overlay if present was not called since the last call to
		// EmuD3DDevice_UpdateOverlay.
		if (g_OverlaySwap != g_Xbox_SwapData.Swap - 1) {
			EMUPATCH(D3DDevice_Swap)(CXBX_SWAP_PRESENT_FORWARD);
		}

		g_OverlaySwap = g_Xbox_SwapData.Swap;
	}
}

// ******************************************************************
// * patch: D3DDevice_UpdateOverlay
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_UpdateOverlay)
(
	X_D3DSurface *pSurface,
	CONST RECT   *SrcRect,
	CONST RECT   *DstRect,
	bool_xt       EnableColorKey,
	D3DCOLOR      ColorKey
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pSurface)
		LOG_FUNC_ARG(SrcRect)
		LOG_FUNC_ARG(DstRect)
		LOG_FUNC_ARG(EnableColorKey)
		LOG_FUNC_ARG(ColorKey)
		LOG_FUNC_END;

	CxbxrImpl_UpdateOverlay(pSurface, SrcRect, DstRect, EnableColorKey, ColorKey);
}

// ******************************************************************
// * patch: D3DDevice_UpdateOverlay_16__LTCG_eax2
// ******************************************************************
static void D3DDevice_UpdateOverlay_16__LTCG_eax2
(
	xbox::X_D3DSurface *pSurface,
	CONST RECT         *SrcRect,
	CONST RECT         *DstRect,
	xbox::bool_xt       EnableColorKey,
	D3DCOLOR            ColorKey
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pSurface)
		LOG_FUNC_ARG(SrcRect)
		LOG_FUNC_ARG(DstRect)
		LOG_FUNC_ARG(EnableColorKey)
		LOG_FUNC_ARG(ColorKey)
		LOG_FUNC_END;
}

// This uses a custom calling convention where parameter is passed in EAX
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_UpdateOverlay_16__LTCG_eax2)
(
	X_D3DSurface *pSurface,
	CONST RECT   *DstRect,
	bool_xt       EnableColorKey,
	D3DCOLOR      ColorKey
)
{
	RECT* SrcRect;
	__asm {
		LTCG_PROLOGUE
		mov  SrcRect, eax
	}

	// Log
	D3DDevice_UpdateOverlay_16__LTCG_eax2(pSurface, SrcRect, DstRect, EnableColorKey, ColorKey);

	CxbxrImpl_UpdateOverlay(pSurface, SrcRect, DstRect, EnableColorKey, ColorKey);

	__asm {
		LTCG_EPILOGUE
		ret  16
	}
}

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
// Xbox native code waits on VBlank event; NV2A VBlank IRQ signals it correctly.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DResource_BlockUntilNotBusy — disabled.
// Empty stub; Xbox native code polls resource state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// ******************************************************************
// * patch: D3DDevice_InsertCallback
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_InsertCallback)
(
	X_D3DCALLBACKTYPE	Type,
	X_D3DCALLBACK		pCallback,
	dword_xt				Context
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Type)
		LOG_FUNC_ARG(pCallback)
		LOG_FUNC_ARG(Context)
		LOG_FUNC_END;

	CxbxImpl_InsertCallback(Type, pCallback, Context);

	LOG_INCOMPLETE();
}

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

// D3D_BlockOnTime — disabled.
// Xbox native code polls NV_PFIFO_CACHE1_DMA_GET until it equals PUT.
// PFIFO read handler fast-path returns GET=PUT in HLE mode, so the native
// polling loop exits immediately. Implementation moved to Direct3D9.cpp.unused-patches.

// D3D_DestroyResource — disabled.
// Host resources are NV2A-derived (PGRAPH RT cache, texture cache keyed by VRAM).
// Dirty page tracking invalidates stale host resources.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

