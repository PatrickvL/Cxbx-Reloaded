// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
// ******************************************************************
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
// *  (c) 2018 Luke Usher <luke.usher@outlook.com>
// *
// *  All rights reserved
// *
// ******************************************************************

#include "core\kernel\init\CxbxKrnl.h"
#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\Rendering/RenderGlobals.h"
#include "core\hle\JVS\JVS.h"
#include "Patches.hpp"
#include "Intercept.hpp"
#include "devices\Xbox.h"

namespace xbox {
// Minimal HLE: replace DirectSoundDoWork to prevent the native
// implementation from spinning or crashing when the LLE APU
// DSP is not fully functional.
void_xt WINAPI EMUPATCH(DirectSoundDoWork)()
{
	if (g_APU != nullptr) g_APU->SynchronizeAudio();
	if (g_AC97 != nullptr) g_AC97->ServiceAudio();
}
} // namespace xbox

#include <map>
#include <unordered_map>
#include <subhook.h>

typedef struct {
	const void* patchFunc;		// Function pointer of the patch in Cxbx-R codebase
	const uint32_t flags;		// Patch Flags
} xbox_patch_t;

const uint32_t PATCH_ALWAYS = 1 << 0;
const uint32_t PATCH_HLE_D3D = 1 << 1;
const uint32_t PATCH_HLE_OHCI = 1 << 2;
const uint32_t PATCH_IS_FIBER = 1 << 3;
const uint32_t PATCH_HLE_DSOUND = 1 << 4;

#define PATCH_ENTRY(Name, Func, Flags) \
    { Name, xbox_patch_t { (void *)&Func, Flags} }

// Map of Xbox Patch names to Emulator Patches
// A std::string is used as it's possible for a symbol to have multiple names
// This allows for the eventual importing of Dxbx symbol files and even IDA signatures too!
std::map<const std::string, const xbox_patch_t> g_PatchTable = {
	// Direct3D
	// Step 10: Let Xbox native code drive draws through push buffer → PFIFO → PGRAPH → HLE callbacks.
	// CDevice::SetStateVB/UP flush deferred render state to push buffer before draws.
	// Removing these patches lets the native Xbox code write NV2A commands directly.
	//PATCH_ENTRY("CDevice_SetStateUP", xbox::EMUPATCH(CDevice_SetStateUP), PATCH_HLE_D3D),
	//PATCH_ENTRY("CDevice_SetStateUP_4", xbox::EMUPATCH(CDevice_SetStateUP_4), PATCH_HLE_D3D),
	//PATCH_ENTRY("CDevice_SetStateUP_0__LTCG_esi1", xbox::EMUPATCH(CDevice_SetStateUP_0__LTCG_esi1), PATCH_HLE_D3D),
	//PATCH_ENTRY("CDevice_SetStateVB", xbox::EMUPATCH(CDevice_SetStateVB), PATCH_HLE_D3D),
	//PATCH_ENTRY("CDevice_SetStateVB_8", xbox::EMUPATCH(CDevice_SetStateVB_8), PATCH_HLE_D3D),
	// Step 10: Begin/End/SetVertexData flow through push buffer as inline vertex commands.
	//PATCH_ENTRY("D3DDevice_Begin", xbox::EMUPATCH(D3DDevice_Begin), PATCH_HLE_D3D),
	// Disabled: Xbox native BeginPush returns pointer into GPU ring buffer (CDevice+0x00 = pPut).
	// KickOff writes DMA_PUT. PFIFO DMA pusher processes ring commands automatically.
	//PATCH_ENTRY("D3DDevice_BeginPush_4", xbox::EMUPATCH(D3DDevice_BeginPush_4), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_BeginPush_8", xbox::EMUPATCH(D3DDevice_BeginPush_8), PATCH_HLE_D3D),
	// Disabled: trampoline-only, toggles unused g_bRecordingPushBuffer flag
	//PATCH_ENTRY("D3DDevice_BeginPushBuffer", xbox::EMUPATCH(D3DDevice_BeginPushBuffer), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_BeginPushBuffer_0__LTCG_edi1", xbox::EMUPATCH(D3DDevice_BeginPushBuffer_0__LTCG_edi1), PATCH_HLE_D3D),
	// Disabled: Visibility tests now handled natively via NV2A PGRAPH.
	// Xbox code pushes NV097_CLEAR_REPORT_VALUE, SET_ZPASS_PIXEL_COUNT_ENABLE, and
	// NV097_GET_REPORT through PFIFO. PGRAPH wraps D3D11 draws with occlusion queries
	// and accumulates z-passing pixel counts into pg->zpass_pixel_count_result.
	//PATCH_ENTRY("D3DDevice_BeginVisibilityTest", xbox::EMUPATCH(D3DDevice_BeginVisibilityTest), PATCH_HLE_D3D),
	// Disabled: empty LOG_UNIMPLEMENTED stub, Xbox native code writes harmless debug regs
	//PATCH_ENTRY("D3DDevice_BlockOnFence", xbox::EMUPATCH(D3DDevice_BlockOnFence), PATCH_HLE_D3D),
	// Disabled: Native BlockUntilVerticalBlank waits on m_VerticalBlankEvent KEVENT via
	// KeWaitForSingleObject.  The miniport DPC signals it via KeSetEvent.  The kernel
	// missed-wakeup fix (re-checking SignalState in WaitApc) ensures this works correctly.
	//PATCH_ENTRY("D3DDevice_BlockUntilVerticalBlank", xbox::EMUPATCH(D3DDevice_BlockUntilVerticalBlank), PATCH_HLE_D3D),
	// Disabled: Clear now handled by D3D11_draw_clear via NV097_CLEAR_SURFACE → pgraph_handle_method
	//PATCH_ENTRY("D3DDevice_Clear", xbox::EMUPATCH(D3DDevice_Clear), PATCH_HLE_D3D),
	// Disabled: Native CopyRects does CPU memcpy between Xbox surfaces. PVIDEO overlay reads
	// from contiguous memory which is synced from tiled pages by CxbxSyncTiledRangeToContiguous.
	//PATCH_ENTRY("D3DDevice_CopyRects", xbox::EMUPATCH(D3DDevice_CopyRects), PATCH_HLE_D3D),
	// PATCH_ENTRY("D3DDevice_CreateVertexShader", xbox::EMUPATCH(D3DDevice_CreateVertexShader), PATCH_HLE_D3D),
	// ================================================================
	// BATCH DISABLE: All state-caching D3D patches disabled.
	// Xbox native code writes NV2A commands through PFIFO → PGRAPH.
	// Render thread reads authoritative state from PGRAPH registers.
	// ================================================================

	// Vertex shader management
	//PATCH_ENTRY("D3DDevice_DeleteVertexShader", xbox::EMUPATCH(D3DDevice_DeleteVertexShader), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DeleteVertexShader_0__LTCG_eax1", xbox::EMUPATCH(D3DDevice_DeleteVertexShader_0__LTCG_eax1), PATCH_HLE_D3D),
	// Step 10: Draw calls flow through push buffer as NV097_DRAW_ARRAYS, NV097_INLINE_ARRAY,
	// NV097_ARRAY_ELEMENT16/32 commands. PGRAPH SET_BEGIN_END triggers D3D11_draw_* callbacks.
	//PATCH_ENTRY("D3DDevice_DrawIndexedVertices", xbox::EMUPATCH(D3DDevice_DrawIndexedVertices), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawIndexedVerticesUP", xbox::EMUPATCH(D3DDevice_DrawIndexedVerticesUP), PATCH_HLE_D3D),
	// Disabled: Xbox CPU tessellation emits normal NV097 draw commands through push buffer.
	//PATCH_ENTRY("D3DDevice_DrawRectPatch", xbox::EMUPATCH(D3DDevice_DrawRectPatch), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawTriPatch", xbox::EMUPATCH(D3DDevice_DrawTriPatch), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawVertices", xbox::EMUPATCH(D3DDevice_DrawVertices), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawVertices_4__LTCG_ecx2_eax3", xbox::EMUPATCH(D3DDevice_DrawVertices_4__LTCG_ecx2_eax3), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawVertices_8__LTCG_eax3", xbox::EMUPATCH(D3DDevice_DrawVertices_8__LTCG_eax3), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawVerticesUP", xbox::EMUPATCH(D3DDevice_DrawVerticesUP), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_DrawVerticesUP_12__LTCG_ebx3", xbox::EMUPATCH(D3DDevice_DrawVerticesUP_12__LTCG_ebx3), PATCH_HLE_D3D),
	// Disabled: Native EnableOverlay programs NV_PVIDEO_STOP; D3D11_flip_stall reads PVIDEO state.
	//PATCH_ENTRY("D3DDevice_EnableOverlay", xbox::EMUPATCH(D3DDevice_EnableOverlay), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_EnableOverlay_0__LTCG", xbox::EMUPATCH(D3DDevice_EnableOverlay_0__LTCG), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_End", xbox::EMUPATCH(D3DDevice_End), PATCH_HLE_D3D),
	// Disabled: Xbox native EndPush updates pPut in CDevice. PFIFO processes commands from ring buffer.
	//PATCH_ENTRY("D3DDevice_EndPush", xbox::EMUPATCH(D3DDevice_EndPush), PATCH_HLE_D3D),
	// Disabled: trampoline-only, clears unused g_bRecordingPushBuffer flag
	//PATCH_ENTRY("D3DDevice_EndPushBuffer", xbox::EMUPATCH(D3DDevice_EndPushBuffer), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_EndVisibilityTest", xbox::EMUPATCH(D3DDevice_EndVisibilityTest), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_EndVisibilityTest_0__LTCG_eax1", xbox::EMUPATCH(D3DDevice_EndVisibilityTest_0__LTCG_eax1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_FlushVertexCache", xbox::EMUPATCH(D3DDevice_FlushVertexCache), PATCH_HLE_D3D),
	// Disabled: Xbox native GetBackBuffer returns correct Xbox surface pointers.
	// No HLE code depends on intercepting the returned X_D3DSurface*.
	//PATCH_ENTRY("D3DDevice_GetBackBuffer", xbox::EMUPATCH(D3DDevice_GetBackBuffer), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetBackBuffer_8__LTCG_eax1", xbox::EMUPATCH(D3DDevice_GetBackBuffer_8__LTCG_eax1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetBackBuffer2", xbox::EMUPATCH(D3DDevice_GetBackBuffer2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetBackBuffer2_0__LTCG_eax1", xbox::EMUPATCH(D3DDevice_GetBackBuffer2_0__LTCG_eax1), PATCH_HLE_D3D),
	// Disabled: Xbox native reads VBlank count from device struct and field from NV_PCRTC_RASTER bit 20.
	// PCRTC emulation now provides VERT_BLANK and FIELD bits; VBlank ISR increments the kernel counter.
	//PATCH_ENTRY("D3DDevice_GetDisplayFieldStatus", xbox::EMUPATCH(D3DDevice_GetDisplayFieldStatus), PATCH_HLE_D3D),
	// Disabled: PRMDIO VGA DAC palette emulation handles gamma natively
	//PATCH_ENTRY("D3DDevice_GetGammaRamp", xbox::EMUPATCH(D3DDevice_GetGammaRamp), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetMaterial", xbox::EMUPATCH(D3DDevice_GetMaterial), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetModelView", xbox::EMUPATCH(D3DDevice_GetModelView), PATCH_HLE_D3D),
	// Disabled: hardcoded TRUE stub, Xbox native overlay check is correct
	//PATCH_ENTRY("D3DDevice_GetOverlayUpdateStatus", xbox::EMUPATCH(D3DDevice_GetOverlayUpdateStatus), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetProjectionViewportMatrix", xbox::EMUPATCH(D3DDevice_GetProjectionViewportMatrix), PATCH_HLE_D3D),
	// Disabled: g_Xbox_VertexShaderConstantMode has no render-thread readers
	//PATCH_ENTRY("D3DDevice_GetShaderConstantMode", xbox::EMUPATCH(D3DDevice_GetShaderConstantMode), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetTransform", xbox::EMUPATCH(D3DDevice_GetTransform), PATCH_HLE_D3D),
	// Disabled: getter reads g_Xbox_VertexShader_Handle; Xbox native reads from device struct
	//PATCH_ENTRY("D3DDevice_GetVertexShader", xbox::EMUPATCH(D3DDevice_GetVertexShader), PATCH_HLE_D3D),
	// Disabled: getter reads HLE VS constant shadow; Xbox native reads from device constant table
	//PATCH_ENTRY("D3DDevice_GetVertexShaderConstant", xbox::EMUPATCH(D3DDevice_GetVertexShaderConstant), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetVertexShaderDeclaration", xbox::EMUPATCH(D3DDevice_GetVertexShaderDeclaration), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetVertexShaderFunction", xbox::EMUPATCH(D3DDevice_GetVertexShaderFunction), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetVertexShaderInput", xbox::EMUPATCH(D3DDevice_GetVertexShaderInput), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetVertexShaderSize", xbox::EMUPATCH(D3DDevice_GetVertexShaderSize), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetVertexShaderType", xbox::EMUPATCH(D3DDevice_GetVertexShaderType), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetViewportOffsetAndScale", xbox::EMUPATCH(D3DDevice_GetViewportOffsetAndScale), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetViewportOffsetAndScale_0__LTCG_edx1_ecx2", xbox::EMUPATCH(D3DDevice_GetViewportOffsetAndScale_0__LTCG_edx1_ecx2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_GetVisibilityTestResult", xbox::EMUPATCH(D3DDevice_GetVisibilityTestResult), PATCH_HLE_D3D),
	// Disabled: Native InsertCallback pushes NV097_NO_OPERATION(param) to the push buffer.
	// PGRAPH raises INTR_ERROR → miniport ISR reads TRAPPED_DATA_LOW → dispatches callback.
	// The HLE version stored callbacks in g_Xbox_CallbackQueue but CxbxHandleXboxCallbacks()
	// was never called, so callbacks were never dispatched anyway.
	//PATCH_ENTRY("D3DDevice_InsertCallback", xbox::EMUPATCH(D3DDevice_InsertCallback), PATCH_HLE_D3D),
	// Disabled: fake 0x8000BEEF stub, Xbox native fence via NV2A reference counter
	//PATCH_ENTRY("D3DDevice_InsertFence", xbox::EMUPATCH(D3DDevice_InsertFence), PATCH_HLE_D3D),
	// Disabled: hardcoded FALSE stub, Xbox native fence check via NV2A
	//PATCH_ENTRY("D3DDevice_IsBusy", xbox::EMUPATCH(D3DDevice_IsBusy), PATCH_HLE_D3D),
	// Disabled: hardcoded FALSE stub, Xbox native fence check via NV2A
	//PATCH_ENTRY("D3DDevice_IsFencePending", xbox::EMUPATCH(D3DDevice_IsFencePending), PATCH_HLE_D3D),
	// Disabled: trampoline-only, Xbox native pushes NV2A light methods directly
	//PATCH_ENTRY("D3DDevice_LightEnable", xbox::EMUPATCH(D3DDevice_LightEnable), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_LightEnable_4__LTCG_eax1", xbox::EMUPATCH(D3DDevice_LightEnable_4__LTCG_eax1), PATCH_HLE_D3D),
	// Disabled: trampoline-only after CxbxImpl_LoadVertexShader removal.
	// Xbox code runs unpatched and writes to push buffer directly.
	//PATCH_ENTRY("D3DDevice_LoadVertexShader", xbox::EMUPATCH(D3DDevice_LoadVertexShader), PATCH_HLE_D3D),
	// Step 11: LoadVertexShaderProgram flows through push buffer → PFIFO → pg->xf.xfpr[].
	// LTCG LoadVertexShader variants kept for logging only (no trampoline available).
	//PATCH_ENTRY("D3DDevice_LoadVertexShaderProgram", xbox::EMUPATCH(D3DDevice_LoadVertexShaderProgram), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_LoadVertexShader_0__LTCG_ecx1_eax2", xbox::EMUPATCH(D3DDevice_LoadVertexShader_0__LTCG_ecx1_eax2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_LoadVertexShader_0__LTCG_edx1_eax2", xbox::EMUPATCH(D3DDevice_LoadVertexShader_0__LTCG_edx1_eax2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_LoadVertexShader_4__LTCG_eax1", xbox::EMUPATCH(D3DDevice_LoadVertexShader_4__LTCG_eax1), PATCH_HLE_D3D),
	// Disabled: Xbox native MultiplyTransform calls SetTransform which pushes NV2A transform methods
	//PATCH_ENTRY("D3DDevice_MultiplyTransform", xbox::EMUPATCH(D3DDevice_MultiplyTransform), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_MultiplyTransform_0__LTCG_ebx1_eax2", xbox::EMUPATCH(D3DDevice_MultiplyTransform_0__LTCG_ebx1_eax2), PATCH_HLE_D3D),
	// Disabled: pure trampoline with LOG_INCOMPLETE. Xbox native runs unpatched.
	//PATCH_ENTRY("D3DDevice_PersistDisplay", xbox::EMUPATCH(D3DDevice_PersistDisplay), PATCH_HLE_D3D),
	// Disabled: Native Swap pushes NV097_FLIP_INCREMENT_WRITE + NV097_FLIP_STALL to push buffer.
	// PGRAPH FLIP_STALL handler calls pgraph_flip_stall → D3D11_flip_stall which blits
	// the PGRAPH backbuffer to the host swap chain and presents.
	//PATCH_ENTRY("D3DDevice_Present", xbox::EMUPATCH(D3DDevice_Present), PATCH_HLE_D3D),
	// Disabled: unimplemented stub (LOG_UNIMPLEMENTED), intercepting does nothing useful
	//PATCH_ENTRY("D3DDevice_PrimeVertexCache", xbox::EMUPATCH(D3DDevice_PrimeVertexCache), PATCH_HLE_D3D),
	// Disabled: Xbox D3DDevice_Reset only does Xbox D3D resource cleanup.
	// Host rendering driven by NV2A state needs no reset — host resources
	// are managed by PGRAPH RT cache and invalidated by dirty-page detection.
	//PATCH_ENTRY("D3DDevice_Reset", xbox::EMUPATCH(D3DDevice_Reset), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_Reset_0__LTCG_edi1", xbox::EMUPATCH(D3DDevice_Reset_0__LTCG_edi1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_Reset_0__LTCG_ebx1", xbox::EMUPATCH(D3DDevice_Reset_0__LTCG_ebx1), PATCH_HLE_D3D),
	// Disabled: native RunPushBuffer applies fixups and pushes commands through
	// the real GPU FIFO.  D3D_BlockOnTime waits for ring space via semaphores.
	//PATCH_ENTRY("D3DDevice_RunPushBuffer", xbox::EMUPATCH(D3DDevice_RunPushBuffer), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_RunPushBuffer_4__LTCG_eax2", xbox::EMUPATCH(D3DDevice_RunPushBuffer_4__LTCG_eax2), PATCH_HLE_D3D),
	// Disabled: PGRAPH now handles NV097_SET_TRANSFORM_DATA + NV097_LAUNCH_TRANSFORM_PROGRAM natively
	//PATCH_ENTRY("D3DDevice_RunVertexStateShader", xbox::EMUPATCH(D3DDevice_RunVertexStateShader), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_RunVertexStateShader_4__LTCG_esi2", xbox::EMUPATCH(D3DDevice_RunVertexStateShader_4__LTCG_esi2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SelectVertexShader", xbox::EMUPATCH(D3DDevice_SelectVertexShader), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SelectVertexShaderDirect", xbox::EMUPATCH(D3DDevice_SelectVertexShaderDirect), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SelectVertexShaderDirect_0__LTCG_eax1_ebx2", xbox::EMUPATCH(D3DDevice_SelectVertexShaderDirect_0__LTCG_eax1_ebx2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SelectVertexShader_0__LTCG_eax1_ebx2", xbox::EMUPATCH(D3DDevice_SelectVertexShader_0__LTCG_eax1_ebx2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SelectVertexShader_4__LTCG_eax1", xbox::EMUPATCH(D3DDevice_SelectVertexShader_4__LTCG_eax1), PATCH_HLE_D3D),
	// Disabled: Xbox native stores scale in device struct. Default 1.0 is correct for most titles.
	// TODO: Read from Xbox device struct in Swap path if non-1.0 scale needed.
	//PATCH_ENTRY("D3DDevice_SetBackBufferScale", xbox::EMUPATCH(D3DDevice_SetBackBufferScale), PATCH_HLE_D3D),
	// Disabled: trampoline-only, Xbox native pushes NV2A back material methods directly
	//PATCH_ENTRY("D3DDevice_SetBackMaterial", xbox::EMUPATCH(D3DDevice_SetBackMaterial), PATCH_HLE_D3D),
	// Disabled: all cases are TODO stubs, intercepting does nothing useful
	//PATCH_ENTRY("D3DDevice_SetDepthClipPlanes", xbox::EMUPATCH(D3DDevice_SetDepthClipPlanes), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetFlickerFilter", xbox::EMUPATCH(D3DDevice_SetFlickerFilter), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetFlickerFilter_0__LTCG_esi1", xbox::EMUPATCH(D3DDevice_SetFlickerFilter_0__LTCG_esi1), PATCH_HLE_D3D),
	// Disabled: PRMDIO VGA DAC palette emulation handles gamma natively
	//PATCH_ENTRY("D3DDevice_SetGammaRamp", xbox::EMUPATCH(D3DDevice_SetGammaRamp), PATCH_HLE_D3D),
	// Disabled: trampoline-only, g_Xbox_BaseVertexIndex was only reader (DrawIndexedVertices disabled)
	//PATCH_ENTRY("D3DDevice_SetIndices", xbox::EMUPATCH(D3DDevice_SetIndices), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetIndices_4__LTCG_ebx1", xbox::EMUPATCH(D3DDevice_SetIndices_4__LTCG_ebx1), PATCH_HLE_D3D),
	// Disabled: trampoline-only, Xbox native pushes NV2A light/material methods directly
	//PATCH_ENTRY("D3DDevice_SetLight", xbox::EMUPATCH(D3DDevice_SetLight), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetMaterial", xbox::EMUPATCH(D3DDevice_SetMaterial), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetModelView", xbox::EMUPATCH(D3DDevice_SetModelView), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetPalette", xbox::EMUPATCH(D3DDevice_SetPalette), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetPalette_4__LTCG_eax1", xbox::EMUPATCH(D3DDevice_SetPalette_4__LTCG_eax1), PATCH_HLE_D3D),
	// SetPixelShader patches disabled: CxbxImpl_SetPixelShader only wrote g_pXbox_PixelShader
	// which was only read by the COMBINECTL==0 HLE bridge fallback, now removed.
	// Xbox native SetPixelShader pushes all combiner state through PFIFO → PGRAPH.
	//PATCH_ENTRY("D3DDevice_SetPixelShader", xbox::EMUPATCH(D3DDevice_SetPixelShader), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetPixelShaderConstant", xbox::EMUPATCH(D3DDevice_SetPixelShaderConstant), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetPixelShaderConstant_4__LTCG_ecx1_eax3", xbox::EMUPATCH(D3DDevice_SetPixelShaderConstant_4__LTCG_ecx1_eax3), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetPixelShader_0__LTCG_eax1", xbox::EMUPATCH(D3DDevice_SetPixelShader_0__LTCG_eax1), PATCH_HLE_D3D),
	// Step 9.1: PGRAPH regs[] is authoritative for render state (Step 6.2). Xbox native code
	// already writes D3D__RenderState[] before calling SetRenderState_Simple, so this patch
	// was only doing a redundant write. The pushbuffer method goes to PGRAPH via the puller.
	//PATCH_ENTRY("D3DDevice_SetRenderState_Simple", xbox::EMUPATCH(D3DDevice_SetRenderState_Simple), PATCH_HLE_D3D),
	// Disabled: CxbxD3D11UpdateRenderTargetFromPGRAPH creates host RTs/DSs directly from PGRAPH state.
	// Present uses g_pHostPgraphBackBuffer instead of g_pXbox_BackBufferSurface.
	//PATCH_ENTRY("D3DDevice_SetRenderTarget", xbox::EMUPATCH(D3DDevice_SetRenderTarget), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetRenderTargetFast", xbox::EMUPATCH(D3DDevice_SetRenderTargetFast), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetRenderTarget_0__LTCG_ecx1_eax2", xbox::EMUPATCH(D3DDevice_SetRenderTarget_0__LTCG_ecx1_eax2), PATCH_HLE_D3D),
	// Step 12: SetScreenSpaceOffset only wrote g_Xbox_ScreenSpaceOffset which was only
	// read by the dead CxbxSetVertexShaderPassthroughProgram. Xbox code handles it natively.
	//PATCH_ENTRY("D3DDevice_SetScreenSpaceOffset", xbox::EMUPATCH(D3DDevice_SetScreenSpaceOffset), PATCH_HLE_D3D),
	// Disabled: g_Xbox_VertexShaderConstantMode has no render-thread readers
	//PATCH_ENTRY("D3DDevice_SetShaderConstantMode", xbox::EMUPATCH(D3DDevice_SetShaderConstantMode), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetShaderConstantMode_0__LTCG_eax1", xbox::EMUPATCH(D3DDevice_SetShaderConstantMode_0__LTCG_eax1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetSoftDisplayFilter", xbox::EMUPATCH(D3DDevice_SetSoftDisplayFilter), PATCH_HLE_D3D),
	// Disabled: unimplemented stub (LOG_IGNORED), intercepting does nothing useful
	//PATCH_ENTRY("D3DDevice_SetStipple", xbox::EMUPATCH(D3DDevice_SetStipple), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetStreamSource", xbox::EMUPATCH(D3DDevice_SetStreamSource), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetStreamSource_0__LTCG_eax1_edi2_ebx3", xbox::EMUPATCH(D3DDevice_SetStreamSource_0__LTCG_eax1_edi2_ebx3), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetStreamSource_4__LTCG_eax1_ebx2", xbox::EMUPATCH(D3DDevice_SetStreamSource_4__LTCG_eax1_ebx2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetStreamSource_8__LTCG_eax1", xbox::EMUPATCH(D3DDevice_SetStreamSource_8__LTCG_eax1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetStreamSource_8__LTCG_edx1", xbox::EMUPATCH(D3DDevice_SetStreamSource_8__LTCG_edx1), PATCH_HLE_D3D),
	// Disabled: Xbox native stores callback in device struct. Swap callback not invoked without this.
	// TODO: Read callback pointer from Xbox device struct in Swap path if needed.
	//PATCH_ENTRY("D3DDevice_SetSwapCallback", xbox::EMUPATCH(D3DDevice_SetSwapCallback), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetTexture", xbox::EMUPATCH(D3DDevice_SetTexture), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetTexture_4__LTCG_eax2", xbox::EMUPATCH(D3DDevice_SetTexture_4__LTCG_eax2), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetTexture_4__LTCG_eax1", xbox::EMUPATCH(D3DDevice_SetTexture_4__LTCG_eax1), PATCH_HLE_D3D),
	// Disabled: transform state now sourced from PGRAPH XFCTX registers (MMAT0/CMAT/TnMAT)
	//PATCH_ENTRY("D3DDevice_SetTransform", xbox::EMUPATCH(D3DDevice_SetTransform), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetTransform_0__LTCG_eax1_edx2", xbox::EMUPATCH(D3DDevice_SetTransform_0__LTCG_eax1_edx2), PATCH_HLE_D3D),
	// Step 10: SetVertexData writes NV097_SET_VERTEX_DATA* push buffer commands → inline_buffer path
	//PATCH_ENTRY("D3DDevice_SetVertexData2f", xbox::EMUPATCH(D3DDevice_SetVertexData2f), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexData2s", xbox::EMUPATCH(D3DDevice_SetVertexData2s), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexData4f", xbox::EMUPATCH(D3DDevice_SetVertexData4f), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexData4f_16__LTCG_edi1", xbox::EMUPATCH(D3DDevice_SetVertexData4f_16__LTCG_edi1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexData4s", xbox::EMUPATCH(D3DDevice_SetVertexData4s), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexData4ub", xbox::EMUPATCH(D3DDevice_SetVertexData4ub), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexDataColor", xbox::EMUPATCH(D3DDevice_SetVertexDataColor), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShader", xbox::EMUPATCH(D3DDevice_SetVertexShader), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShader_0__LTCG_ebx1", xbox::EMUPATCH(D3DDevice_SetVertexShader_0__LTCG_ebx1), PATCH_HLE_D3D),
	// Step 9.2: VS constants flow through pushbuffer → puller → pg->xf.xfctx[]
	// with dirty tracking. pfifo_flush before draws ensures constants are up to date.
	// Xbox native SetVertexShaderConstant generates NV097_SET_TRANSFORM_CONSTANT push buffer
	// commands, which pgraph_handle_method writes to pg->xf.xfctx[] with dirty flags.
	// Removing these patches lets the push buffer path be the sole source of truth.
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstant", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstant), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstant1", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstant1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstant1Fast", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstant1Fast), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstant4", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstant4), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstantNotInline", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstantNotInline), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstantNotInline_0__LTCG_ebx1_edx2_eax3", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstantNotInline_0__LTCG_ebx1_edx2_eax3), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstantNotInlineFast", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstantNotInlineFast), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderConstant_8__LTCG_edx3", xbox::EMUPATCH(D3DDevice_SetVertexShaderConstant_8__LTCG_edx3), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderInput", xbox::EMUPATCH(D3DDevice_SetVertexShaderInput), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVertexShaderInputDirect", xbox::EMUPATCH(D3DDevice_SetVertexShaderInputDirect), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SetVerticalBlankCallback", xbox::EMUPATCH(D3DDevice_SetVerticalBlankCallback), PATCH_HLE_D3D),
	// Disabled: trampoline-only after CxbxImpl_SetViewport removal.
	// Xbox code runs unpatched and writes viewport to push buffer directly.
	//PATCH_ENTRY("D3DDevice_SetViewport", xbox::EMUPATCH(D3DDevice_SetViewport), PATCH_HLE_D3D),
	// Disabled: Native Swap pushes NV097_FLIP_STALL → pgraph_flip_stall → D3D11_flip_stall.
	//PATCH_ENTRY("D3DDevice_Swap", xbox::EMUPATCH(D3DDevice_Swap), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_Swap_0__LTCG_eax1", xbox::EMUPATCH(D3DDevice_Swap_0__LTCG_eax1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_SwitchTexture", xbox::EMUPATCH(D3DDevice_SwitchTexture), PATCH_HLE_D3D),
	// Disabled: Native UpdateOverlay programs PVIDEO registers; D3D11_flip_stall composites overlay.
	//PATCH_ENTRY("D3DDevice_UpdateOverlay", xbox::EMUPATCH(D3DDevice_UpdateOverlay), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3DDevice_UpdateOverlay_16__LTCG_eax2", xbox::EMUPATCH(D3DDevice_UpdateOverlay_16__LTCG_eax2), PATCH_HLE_D3D),
	// Disabled: empty LOG_UNIMPLEMENTED stub, Xbox native code polls resource state
	//PATCH_ENTRY("D3DResource_BlockUntilNotBusy", xbox::EMUPATCH(D3DResource_BlockUntilNotBusy), PATCH_HLE_D3D),
	// D3D_BlockOnTime: unpatched — native D3D_BlockOnTime uses semaphore,
	// wait-for-idle, and nop methods pushed into the command buffer.
	// The DMA_PUT write handler processes commands inline (pfifo_run_pusher),
	// so completion signals fire and the native wait returns naturally.
	//PATCH_ENTRY("D3D_BlockOnTime", xbox::EMUPATCH(D3D_BlockOnTime), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3D_BlockOnTime_4__LTCG_eax1", xbox::EMUPATCH(D3D_BlockOnTime_4__LTCG_eax1), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3D_CommonSetRenderTarget", xbox::EMUPATCH(D3D_CommonSetRenderTarget), PATCH_HLE_D3D),
	// Disabled: host resources are NV2A-derived (PGRAPH RT cache, texture cache keyed by VRAM).
	// Dirty page tracking invalidates stale host resources. FreeHostResource on Xbox D3D keys is vestigial.
    //PATCH_ENTRY("D3D_DestroyResource", xbox::EMUPATCH(D3D_DestroyResource), PATCH_HLE_D3D),
	//PATCH_ENTRY("D3D_DestroyResource_0__LTCG_edi1", xbox::EMUPATCH(D3D_DestroyResource_0__LTCG_edi1), PATCH_HLE_D3D),
	// Disabled: unimplemented stub (LOG_UNIMPLEMENTED), intercepting does nothing useful
	//PATCH_ENTRY("D3D_LazySetPointParams", xbox::EMUPATCH(D3D_LazySetPointParams), PATCH_HLE_D3D),
	// Disabled: empty LOG_UNIMPLEMENTED stub, Xbox native code writes harmless debug regs
	//PATCH_ENTRY("D3D_SetCommonDebugRegisters", xbox::EMUPATCH(D3D_SetCommonDebugRegisters), PATCH_HLE_D3D),
	// Disabled: Host D3D11 device now created early by CxbxInitHostD3DDevice().
	// Native Xbox CreateDevice runs unpatched, populating D3D_g_pDevice and all internal state.
	//PATCH_ENTRY("Direct3D_CreateDevice", xbox::EMUPATCH(Direct3D_CreateDevice), PATCH_HLE_D3D),
	//PATCH_ENTRY("Direct3D_CreateDevice_16__LTCG_eax4_ebx6", xbox::EMUPATCH(Direct3D_CreateDevice_16__LTCG_eax4_ebx6), PATCH_HLE_D3D),
	//PATCH_ENTRY("Direct3D_CreateDevice_16__LTCG_eax4_ecx6", xbox::EMUPATCH(Direct3D_CreateDevice_16__LTCG_eax4_ecx6), PATCH_HLE_D3D),
	//PATCH_ENTRY("Direct3D_CreateDevice_4__LTCG_eax1_ecx3", xbox::EMUPATCH(Direct3D_CreateDevice_4__LTCG_eax1_ecx3), PATCH_HLE_D3D),
	// Disabled: trampoline-only (LOG_FUNC + XB_TRMP), no side effects
	//PATCH_ENTRY("Lock2DSurface", xbox::EMUPATCH(Lock2DSurface), PATCH_HLE_D3D),
	//PATCH_ENTRY("Lock2DSurface_16__LTCG_esi4_eax5", xbox::EMUPATCH(Lock2DSurface_16__LTCG_esi4_eax5), PATCH_HLE_D3D),
	//PATCH_ENTRY("Lock3DSurface", xbox::EMUPATCH(Lock3DSurface), PATCH_HLE_D3D),
	//PATCH_ENTRY("Lock3DSurface_16__LTCG_eax4", xbox::EMUPATCH(Lock3DSurface_16__LTCG_eax4), PATCH_HLE_D3D),

	// OHCI
	PATCH_ENTRY("XGetDeviceChanges", xbox::EMUPATCH(XGetDeviceChanges), PATCH_HLE_OHCI),
	PATCH_ENTRY("XGetDeviceEnumerationStatus", xbox::EMUPATCH(XGetDeviceEnumerationStatus), PATCH_HLE_OHCI),
	PATCH_ENTRY("XGetDevices", xbox::EMUPATCH(XGetDevices), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInitDevices", xbox::EMUPATCH(XInitDevices), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputClose", xbox::EMUPATCH(XInputClose), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputGetCapabilities", xbox::EMUPATCH(XInputGetCapabilities), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputGetDeviceDescription", xbox::EMUPATCH(XInputGetDeviceDescription), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputGetState", xbox::EMUPATCH(XInputGetState), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputOpen", xbox::EMUPATCH(XInputOpen), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputPoll", xbox::EMUPATCH(XInputPoll), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputSetLightgunCalibration", xbox::EMUPATCH(XInputSetLightgunCalibration), PATCH_HLE_OHCI),
	PATCH_ENTRY("XInputSetState", xbox::EMUPATCH(XInputSetState), PATCH_HLE_OHCI),

	// XAPI
	PATCH_ENTRY("ConvertThreadToFiber", xbox::EMUPATCH(ConvertThreadToFiber), PATCH_IS_FIBER),
	PATCH_ENTRY("CreateFiber", xbox::EMUPATCH(CreateFiber), PATCH_IS_FIBER),
	PATCH_ENTRY("DeleteFiber", xbox::EMUPATCH(DeleteFiber), PATCH_IS_FIBER),
	//PATCH_ENTRY("GetExitCodeThread", xbox::EMUPATCH(GetExitCodeThread), PATCH_ALWAYS),
	//PATCH_ENTRY("GetThreadPriority", xbox::EMUPATCH(GetThreadPriority), PATCH_ALWAYS),
	//PATCH_ENTRY("OutputDebugStringA", xbox::EMUPATCH(OutputDebugStringA), PATCH_ALWAYS), // Game's original code calls DbgPrint which is fully implemented
	//PATCH_ENTRY("RaiseException", xbox::EMUPATCH(RaiseException), PATCH_ALWAYS),
	//PATCH_ENTRY("SetThreadPriority", xbox::EMUPATCH(SetThreadPriority), PATCH_ALWAYS),
	//PATCH_ENTRY("SetThreadPriorityBoost", xbox::EMUPATCH(SetThreadPriorityBoost), PATCH_ALWAYS),
	//PATCH_ENTRY("SignalObjectAndWait", xbox::EMUPATCH(SignalObjectAndWait), PATCH_ALWAYS), // Game's original code calls NtSignalAndWaitForSingleObjectEx which is fully implemented
	PATCH_ENTRY("SwitchToFiber", xbox::EMUPATCH(SwitchToFiber), PATCH_IS_FIBER),
	PATCH_ENTRY("XMountMUA", xbox::EMUPATCH(XMountMUA), PATCH_ALWAYS),
	PATCH_ENTRY("XMountMURootA", xbox::EMUPATCH(XMountMURootA), PATCH_ALWAYS),
	//PATCH_ENTRY("XSetProcessQuantumLength", xbox::EMUPATCH(XSetProcessQuantumLength), PATCH_ALWAYS),
	//PATCH_ENTRY("timeKillEvent", xbox::EMUPATCH(timeKillEvent), PATCH_ALWAYS),
	//PATCH_ENTRY("timeSetEvent", xbox::EMUPATCH(timeSetEvent), PATCH_ALWAYS),
	PATCH_ENTRY("XReadMUMetaData", xbox::EMUPATCH(XReadMUMetaData), PATCH_ALWAYS),
	PATCH_ENTRY("XUnmountMU", xbox::EMUPATCH(XUnmountMU), PATCH_ALWAYS),

	// JVS Functions
	PATCH_ENTRY("JVS_SendCommand", xbox::EMUPATCH(JVS_SendCommand), PATCH_ALWAYS),
	PATCH_ENTRY("JVS_SendCommand2", xbox::EMUPATCH(JVS_SendCommand), PATCH_ALWAYS),
	PATCH_ENTRY("JVS_SendCommand3", xbox::EMUPATCH(JVS_SendCommand), PATCH_ALWAYS),
	PATCH_ENTRY("JvsBACKUP_Read", xbox::EMUPATCH(JvsBACKUP_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsBACKUP_Read2", xbox::EMUPATCH(JvsBACKUP_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsBACKUP_Read3", xbox::EMUPATCH(JvsBACKUP_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsBACKUP_Write", xbox::EMUPATCH(JvsBACKUP_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsBACKUP_Write2", xbox::EMUPATCH(JvsBACKUP_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsEEPROM_Read", xbox::EMUPATCH(JvsEEPROM_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsEEPROM_Read2", xbox::EMUPATCH(JvsEEPROM_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsEEPROM_Read3", xbox::EMUPATCH(JvsEEPROM_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsEEPROM_Write", xbox::EMUPATCH(JvsEEPROM_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsEEPROM_Write2", xbox::EMUPATCH(JvsEEPROM_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsEEPROM_Write3", xbox::EMUPATCH(JvsEEPROM_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareDownload", xbox::EMUPATCH(JvsFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareDownload2", xbox::EMUPATCH(JvsFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareDownload3", xbox::EMUPATCH(JvsFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareDownload4", xbox::EMUPATCH(JvsFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareUpload", xbox::EMUPATCH(JvsFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareUpload2", xbox::EMUPATCH(JvsFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareUpload3", xbox::EMUPATCH(JvsFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsFirmwareUpload4", xbox::EMUPATCH(JvsFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsNodeReceivePacket", xbox::EMUPATCH(JvsNodeReceivePacket), PATCH_ALWAYS),
	PATCH_ENTRY("JvsNodeReceivePacket2", xbox::EMUPATCH(JvsNodeReceivePacket), PATCH_ALWAYS),
	PATCH_ENTRY("JvsNodeSendPacket", xbox::EMUPATCH(JvsNodeSendPacket), PATCH_ALWAYS),
	PATCH_ENTRY("JvsNodeSendPacket2", xbox::EMUPATCH(JvsNodeSendPacket), PATCH_ALWAYS),
	PATCH_ENTRY("JvsRTC_Read", xbox::EMUPATCH(JvsRTC_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsRTC_Read2", xbox::EMUPATCH(JvsRTC_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsRTC_Read3", xbox::EMUPATCH(JvsRTC_Read), PATCH_ALWAYS),
	PATCH_ENTRY("JvsRTC_Write", xbox::EMUPATCH(JvsRTC_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsRTC_Write2", xbox::EMUPATCH(JvsRTC_Write), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareDownload", xbox::EMUPATCH(JvsScFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareDownload2", xbox::EMUPATCH(JvsScFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareDownload3", xbox::EMUPATCH(JvsScFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareDownload4", xbox::EMUPATCH(JvsScFirmwareDownload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareUpload", xbox::EMUPATCH(JvsScFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareUpload2", xbox::EMUPATCH(JvsScFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScFirmwareUpload3", xbox::EMUPATCH(JvsScFirmwareUpload), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScReceiveMidi", xbox::EMUPATCH(JvsScReceiveMidi), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScReceiveMidi2", xbox::EMUPATCH(JvsScReceiveMidi), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScReceiveRs323c", xbox::EMUPATCH(JvsScReceiveRs323c), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScReceiveRs323c2", xbox::EMUPATCH(JvsScReceiveRs323c), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScSendMidi", xbox::EMUPATCH(JvsScSendMidi), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScSendMidi2", xbox::EMUPATCH(JvsScSendMidi), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScSendRs323c", xbox::EMUPATCH(JvsScSendRs323c), PATCH_ALWAYS),
	PATCH_ENTRY("JvsScSendRs323c2", xbox::EMUPATCH(JvsScSendRs323c), PATCH_ALWAYS),

	// DirectSound HLE — replace DoWork to prevent crash/spin
	PATCH_ENTRY("DirectSoundDoWork", xbox::EMUPATCH(DirectSoundDoWork), PATCH_ALWAYS),
};

std::unordered_map<std::string, subhook::Hook> g_FunctionHooks;

inline bool TitleRequiresUnpatchedFibers()
{
    static bool detected = false;
    static bool result = false;

    // Prevent running the check every time this function is called
    if (detected) {
        return result;
    }

    // Array of known games that require the fiber unpatch hack
    DWORD titleIds[] = {
        0x46490002, // Futurama PAL
        0x56550008, // Futurama NTSC
        0
    };

    DWORD* pTitleId = &titleIds[0];
    while (*pTitleId != 0) {
        if (g_pCertificate->dwTitleId == *pTitleId) {
            result = true;
            break;
        }

        pTitleId++;
    }

    detected = true;
    return result;
}


// NOTE: EmuInstallPatch do not get to be in XbSymbolDatabase, do the patches in Cxbx project only.
inline void EmuInstallPatch(const std::string FunctionName, const xbox::addr_xt FunctionAddr)
{
	auto it = g_PatchTable.find(FunctionName);
	if (it == g_PatchTable.end()) {
		return;
	}

	auto patch = it->second;

	if ((patch.flags & PATCH_HLE_OHCI) && bLLE_USB) {
		printf("HLE: %s: Skipped (LLE OHCI Enabled)\n", FunctionName.c_str());
		return;
	}

    // HACK: Some titles require unpatched Fibers, otherwise they enter an infinite loop
    // while others require patched Fibers, otherwise they outright crash
    // This is caused by limitations of Direct Code Execution and Cxbx-R's threading model
    if ((patch.flags & PATCH_IS_FIBER) && TitleRequiresUnpatchedFibers()) {
        printf("HLE: %s: Skipped (Game requires unpatched Fibers)\n", FunctionName.c_str());
        return;
    }

	auto success = g_FunctionHooks[FunctionName].Install((void*)(FunctionAddr), (void*)patch.patchFunc);
	if (success) {
		printf("HLE: %s Patched\n", FunctionName.c_str());
	}
	else {
		printf("HLE: %s Patch Failed\n", FunctionName.c_str());
	}
	
}

void EmuInstallPatches()
{
	for (const auto& it : g_SymbolAddresses) {
		EmuInstallPatch(it.first, it.second);
	}

	LookupTrampolinesD3D();
	LookupTrampolinesXAPI();
}

void* GetPatchedFunctionTrampoline(const std::string functionName)
{
	auto it = g_FunctionHooks.find(functionName);
	if (it != g_FunctionHooks.end()) {
		auto trampoline = it->second.GetTrampoline();
		if (trampoline == nullptr) {
			EmuLogEx(CXBXR_MODULE::HLE, LOG_LEVEL::WARNING, "Failed to get XB_Trampoline for %s", functionName.c_str());
		}

		return trampoline;
	}

	return nullptr;
}
