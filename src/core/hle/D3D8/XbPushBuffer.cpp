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
// *  (c) 2002-2003 Aaron Robinson <caustik@caustik.com>
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::PSHB

#include <assert.h> // For assert()

#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\XbD3D8Types.h" // For X_D3DFORMAT
#include "core\hle\D3D8\ResourceTracker.h"
#include "core\hle\D3D8\Rendering\RenderGlobals.h" // For g_Xbox_VertexShader_Handle
#include "core\hle\D3D8\XbPushBuffer.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11.h" // For CxbxD3D11VertexFetchDraw
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11_Profiler.h"
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11_PageTracker.h"
#include "core\hle\D3D8\Rendering\PatchDraw.h" // For D3D11_draw_patch
#include "core\hle\D3D8\XbVertexShader.h" // For D3D11_launch_transform_program
#include "common/AddressRanges.h" // For CONTIGUOUS_MEMORY_BASE
#include "core/common/video/RenderBase.hpp" // For g_renderbase
#include "devices/video/nv2a.h" // For g_NV2A, PGRAPHState
#include "devices/video/nv2a_int.h" // For NV** defines
#include "core/hle/D3D8/Rendering/NV2A_PGRAPH_Helpers.h"
#include "Logging.h"

// TODO: Find somewhere to put this that doesn't conflict with xbox::
extern void CxbxUpdateHostTextures();

// Persistent PVIDEO overlay texture — DYNAMIC so YUY2→ARGB writes directly into mapped GPU memory
static UINT g_OverlayTexWidth = 0, g_OverlayTexHeight = 0;
ID3D11Texture2D *g_pOverlayTex = nullptr;

void CxbxReleaseOverlayResources()
{
	g_OverlayTexWidth = 0;
	g_OverlayTexHeight = 0;
	if (g_pOverlayTex) { g_pOverlayTex->Release(); g_pOverlayTex = nullptr; }
}

const char *NV2AMethodToString(DWORD dwMethod); // forward

static void D3D11_draw_arrays(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	for (unsigned int i = 0; i < pg->draw_arrays_length; i++) {
		CxbxDrawContext DrawContext = {};

		DrawContext.XboxPrimitiveType = (xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode;
		DrawContext.dwStartVertex = pg->gl_draw_arrays_start[i];
		DrawContext.dwVertexCount = pg->gl_draw_arrays_count[i];

		CxbxD3D11VertexFetchDraw(DrawContext);
	}
}

static void D3D11_draw_inline_buffer(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	if (pg->inline_buffer_length == 0)
		return;

	CxbxD3D11DrawInlineBuffer(pg);
}

static void D3D11_draw_inline_array(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	// Compute per-vertex stride from NV2A vertex attribute format registers.
	// Inline array data packs all enabled attributes contiguously per vertex,
	// unlike array-based draws which use the stride field from the format register.
	unsigned int nv2a_stride = 0;
	for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
		if (pg->vertex_attributes[i].count != 0) { // count 0 = disabled (format 0 is valid: UB_D3D/D3DCOLOR)
			nv2a_stride += pg->vertex_attributes[i].count * pg->vertex_attributes[i].size;
		}
	}

	UINT VertexCount = (nv2a_stride > 0 && pg->inline_array_length > 0)
		? (pg->inline_array_length * sizeof(DWORD)) / nv2a_stride : 0;

	if (nv2a_stride == 0 || pg->inline_array_length == 0) {
		return;
	}

	if (VertexCount == 0) {
		return;
	}

	CxbxDrawContext DrawContext = {};
	DrawContext.XboxPrimitiveType = (xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode;
	DrawContext.dwVertexCount = VertexCount;
	DrawContext.pXboxVertexStreamZeroData = pg->inline_array;
	DrawContext.uiXboxVertexStreamZeroStride = nv2a_stride;
	DrawContext.bNV2AInlineData = true;

	CxbxD3D11VertexFetchDraw(DrawContext);
}

static void D3D11_draw_inline_elements(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	unsigned int uiIndexCount = pg->inline_elements_length;
	CxbxDrawContext DrawContext = {};

	DrawContext.XboxPrimitiveType = (xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode;
	DrawContext.dwVertexCount = uiIndexCount;
	DrawContext.pXboxIndexData = d->pgraph.inline_elements;

	CxbxD3D11VertexFetchDraw(DrawContext);
}

// Unified draw callback — dispatches to the appropriate draw method based on
// which vertex data buffer was filled between BEGIN and END.
void D3D11_draw(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	if (pg->draw_arrays_length) {
		D3D11_draw_arrays(d);
	} else if (pg->inline_buffer_length) {
		D3D11_draw_inline_buffer(d);
	} else if (pg->inline_array_length) {
		D3D11_draw_inline_array(d);
	} else if (pg->inline_elements_length) {
		D3D11_draw_inline_elements(d);
	}
}

void D3D11_draw_state_update(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	// Vertex attribute inline_value may have changed via NV2A push buffer
	// (SET_VERTEX_DATA4F/4UB/2S) since the last draw. Mark defaults dirty
	// so the vertex fetch re-uploads them before the next draw.
	g_bD3D11VertexFetchDefaultsDirty = true;

	// Only invalidate layout CB when vertex_attributes actually changed
	// (FORMAT or OFFSET methods bumped the generation counter in PGRAPH).
	// This avoids expensive per-draw layout CB re-upload for consecutive
	// draws that share the same vertex format.
	{
		static uint32_t s_LastVertexAttribGeneration = ~0u;
		if (pg->vertex_attributes_generation != s_LastVertexAttribGeneration) {
			s_LastVertexAttribGeneration = pg->vertex_attributes_generation;
			CxbxD3D11VertexFetchInvalidateLayout();
		}
	}

	CxbxUpdateNativeD3DResources();

	LOG_INCOMPLETE(); // TODO : Read state from pgraph, convert to D3D
}

// ---- NV2A Zpass pixel count (visibility test) via D3D11 occlusion queries ----

// Persistent occlusion query reused across draw calls.
// Created on first use; Begin/End bracket each draw when zpass counting is enabled.
static ID3D11Query* g_pZpassQuery = nullptr;
static bool g_bZpassQueryPending = false; // true = query ended but result not yet collected

// Collect any pending zpass query result (non-blocking first, blocking if forced)
static void CollectPendingZpassResult(PGRAPHState *pg, bool bBlock)
{
	if (!g_bZpassQueryPending)
		return;

	UINT64 pixelCount = 0;
	HRESULT hr = g_pD3DDeviceContext->GetData(g_pZpassQuery, &pixelCount, sizeof(pixelCount),
		D3D11_ASYNC_GETDATA_DONOTFLUSH);
	if (hr == S_OK) {
		pg->zpass_pixel_count_result += (unsigned int)pixelCount;
		g_bZpassQueryPending = false;
		return;
	}

	if (!bBlock)
		return;

	// Result not ready — flush command buffer and poll until available
	g_pD3DDeviceContext->Flush();
	while (g_pD3DDeviceContext->GetData(g_pZpassQuery, &pixelCount, sizeof(pixelCount), 0) == S_FALSE) {
		// Yield briefly — the flush above ensures GPU is processing
		SwitchToThread();
	}
	pg->zpass_pixel_count_result += (unsigned int)pixelCount;
	g_bZpassQueryPending = false;
}

void D3D11_zpass_begin(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	if (pg->zpass_pixel_count_active)
		return; // already between Begin/End

	if (!g_pD3DDevice || !g_pD3DDeviceContext)
		return;

	// Create the occlusion query on first use
	if (g_pZpassQuery == nullptr) {
		D3D11_QUERY_DESC desc = {};
		desc.Query = D3D11_QUERY_OCCLUSION;
		HRESULT hr = g_pD3DDevice->CreateQuery(&desc, &g_pZpassQuery);
		if (FAILED(hr) || g_pZpassQuery == nullptr)
			return;
	}

	// Collect any still-pending result before reusing the query object
	CollectPendingZpassResult(pg, true);

	g_pD3DDeviceContext->Begin(g_pZpassQuery);
	pg->zpass_pixel_count_active = true;
}

void D3D11_zpass_end(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	if (!pg->zpass_pixel_count_active)
		return; // no query in flight

	if (!g_pD3DDeviceContext || !g_pZpassQuery)
		return;

	g_pD3DDeviceContext->End(g_pZpassQuery);
	pg->zpass_pixel_count_active = false;

	// Mark as pending — result will be collected lazily when needed
	// (at next zpass_begin or GET_REPORT). This avoids stalling the CPU
	// immediately after the draw, giving the GPU time to finish.
	g_bZpassQueryPending = true;
}

static void D3D11_zpass_collect(NV2AState *d)
{
	CollectPendingZpassResult(&d->pgraph, true);
}

// ---- End zpass ----

void D3D11_draw_clear(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	CxbxUpdateNativeD3DResources();

	// Read clear parameters from PGRAPH registers (set by method table dispatch
	// of NV097_SET_CLEAR_RECT_HORIZONTAL/VERTICAL, NV097_SET_COLOR_CLEAR_VALUE,
	// NV097_SET_ZSTENCIL_CLEAR_VALUE before NV097_CLEAR_SURFACE triggers this).
	uint32_t flags = pg->clear_surface_flags;

	// Map NV097 clear flags to host D3D flags.
	// NV097 flags match X_D3DCLEAR values exactly, but host D3DCLEAR_TARGET
	// is a single bit while NV097 has per-channel RGBA bits.
	DWORD hostFlags = 0;
	if (flags & NV097_CLEAR_SURFACE_COLOR)
		hostFlags |= D3DCLEAR_TARGET;
	if (flags & NV097_CLEAR_SURFACE_Z)
		hostFlags |= D3DCLEAR_ZBUFFER;
	if (flags & NV097_CLEAR_SURFACE_STENCIL)
		hostFlags |= D3DCLEAR_STENCIL;

	if (hostFlags == 0)
		return;

	D3DCOLOR color = pg->regs[RI(NV_PGRAPH_COLORCLEARVALUE)];
	uint32_t zstencil = pg->regs[RI(NV_PGRAPH_ZSTENCILCLEARVALUE)];

	// Decode Z and stencil based on the surface zeta format:
	// Z16 (format 1): 16-bit depth in bits [15:0], no stencil
	// Z24S8 (format 2): 24-bit depth in bits [31:8], 8-bit stencil in bits [7:0]
	float z;
	DWORD stencil;
	unsigned int zeta_format = NV2AGetSurfaceState(pg).zetaFormat;
	if (zeta_format == NV097_SET_SURFACE_FORMAT_ZETA_Z16) {
		z = (float)(zstencil & 0xFFFF) / (float)0xFFFF;
		stencil = 0;
		// Z16 has no stencil — strip stencil clear flag
		hostFlags &= ~D3DCLEAR_STENCIL;
	} else {
		// Z24S8 (default)
		z = (float)(zstencil >> 8) / (float)0xFFFFFF;
		stencil = zstencil & 0xFF;
	}

	// Read clear rect from PGRAPH.  Use a single rect covering the clear area.
	uint32_t rectx = pg->regs[RI(NV_PGRAPH_CLEARRECTX)];
	uint32_t recty = pg->regs[RI(NV_PGRAPH_CLEARRECTY)];

	D3DRECT rect;
	rect.left   = rectx & NV_PGRAPH_CLEARRECTX_XMIN;
	rect.right  = (rectx & NV_PGRAPH_CLEARRECTX_XMAX) >> 16;
	rect.top    = recty & NV_PGRAPH_CLEARRECTY_YMIN;
	rect.bottom = (recty & NV_PGRAPH_CLEARRECTY_YMAX) >> 16;

	// NV2A clear rect right/bottom are inclusive; D3D expects exclusive
	rect.right  += 1;
	rect.bottom += 1;

	// Scale for upscale factor and MSAA
	float aaX, aaY;
	GetMultiSampleScaleRaw(aaX, aaY);
	float Xscale = aaX * g_RenderUpscaleFactor;
	float Yscale = aaY * g_RenderUpscaleFactor;
	rect.left   = static_cast<LONG>(rect.left   * Xscale);
	rect.right  = static_cast<LONG>(rect.right  * Xscale);
	rect.top    = static_cast<LONG>(rect.top    * Yscale);
	rect.bottom = static_cast<LONG>(rect.bottom * Yscale);

	CxbxD3DClear(1, &rect, hostFlags, color, z, stencil);
}

// Import pgraph_draw_* variables, declared in EmuNV2A_PGRAPH.cpp :
extern void(*pgraph_draw)(NV2AState *d);
extern void(*pgraph_draw_state_update)(NV2AState *d);
extern void(*pgraph_draw_clear)(NV2AState *d);
extern void(*pgraph_draw_patch)(NV2AState *d);
extern void(*pgraph_flip_stall)(NV2AState *d);
extern void(*pgraph_zpass_begin)(NV2AState *d);
extern void(*pgraph_zpass_end)(NV2AState *d);
extern void(*pgraph_zpass_collect)(NV2AState *d);
extern void(*pgraph_launch_transform_program)(NV2AState *d, unsigned int program_start);

extern void CxbxImGui_RenderD3D(ImGuiUI* m_imgui, ID3D11Texture2D* renderTarget);

// D3D11_flip_stall: Triggered by NV097_FLIP_STALL in the push buffer.
// Blits the PGRAPH-tracked backbuffer to the host swap chain and presents.
static void D3D11_flip_stall(NV2AState *d)
{
	// Get host swap chain backbuffer
	ID3D11Texture2D *pHostBackBuffer = nullptr;
	HRESULT hRet = CxbxGetBackBuffer(&pHostBackBuffer);
	if (hRet != S_OK || !pHostBackBuffer)
		return;

	// Save and restore the game's render target around the present blit.
	// CxbxD3D11Blt manages its own RT state internally, so the PGRAPH RT
	// is effectively unbound during the copy (no resource hazard).
	ID3D11Texture2D* pExistingRT = CxbxGetCurrentRenderTarget();

	// Clear host backbuffer to black (prevents artifacts on aspect ratio change)
	(void)CxbxSetRenderTarget(pHostBackBuffer);
	CxbxD3DClear(0, nullptr, D3DCLEAR_TARGET, 0xFF000000, 1.0f, 0);
	if (pExistingRT) {
		(void)CxbxSetRenderTarget(pExistingRT);
	}

	// Calculate destination rect (centered, aspect-ratio aware)
	float width, height;
	if (g_XBVideo.bMaintainAspect && g_AspectRatioScaleWidth > 0 && g_AspectRatioScaleHeight > 0) {
		width = g_AspectRatioScaleWidth * g_AspectRatioScale;
		height = g_AspectRatioScaleHeight * g_AspectRatioScale;
	} else {
		width = (float)g_HostBackBufferDesc.Width;
		height = (float)g_HostBackBufferDesc.Height;
	}

	// Blit PGRAPH backbuffer to host backbuffer
	auto pXboxBackBufferHostSurface = g_pHostPgraphBackBuffer;
	if (pXboxBackBufferHostSurface) {
		RECT dest{};
		dest.top = (LONG)((g_HostBackBufferDesc.Height - height) / 2);
		dest.left = (LONG)((g_HostBackBufferDesc.Width - width) / 2);
		dest.right = (LONG)(dest.left + width);
		dest.bottom = (LONG)(dest.top + height);

		{
			CXBX_PROFILE_SCOPE(PROF_PRESENT_BLIT);
			CxbxBltSurface(pXboxBackBufferHostSurface, nullptr, pHostBackBuffer, &dest, D3DTEXF_LINEAR);
		}
	}

	// Composite PVIDEO overlay (if enabled by Xbox D3DDevice_UpdateOverlay → PVIDEO registers)
	if (d->enable_overlay) {
		CXBX_PROFILE_SCOPE(PROF_PRESENT_OVERLAY);
		// Determine which buffer is active (bit 0 = buffer 0, bit 4 = buffer 1)
		uint32_t pvideo_buffer = d->pvideo.regs[RI(NV_PVIDEO_BUFFER)];
		int buf = (pvideo_buffer & NV_PVIDEO_BUFFER_0_USE) ? 0 : 1;

		uint32_t pvideo_base = d->pvideo.regs[RI(NV_PVIDEO_BASE(buf))];
		uint32_t pvideo_offset = d->pvideo.regs[RI(NV_PVIDEO_OFFSET(buf))];
		uint32_t pvideo_size_in = d->pvideo.regs[RI(NV_PVIDEO_SIZE_IN(buf))];
		uint32_t pvideo_format = d->pvideo.regs[RI(NV_PVIDEO_FORMAT(buf))];
		uint32_t pvideo_point_out = d->pvideo.regs[RI(NV_PVIDEO_POINT_OUT(buf))];
		uint32_t pvideo_size_out = d->pvideo.regs[RI(NV_PVIDEO_SIZE_OUT(buf))];

		UINT overlayWidth = GET_MASK(pvideo_size_in, NV_PVIDEO_SIZE_IN_WIDTH);
		UINT overlayHeight = GET_MASK(pvideo_size_in, NV_PVIDEO_SIZE_IN_HEIGHT);
		UINT overlayPitch = GET_MASK(pvideo_format, NV_PVIDEO_FORMAT_PITCH);

		if (overlayWidth > 0 && overlayHeight > 0 && overlayPitch > 0) {
			// Sync any tiled pages (0xF0000000 WC mapping) back to contiguous memory
			// so the overlay data is visible at the contiguous address we read below.
			// Games write decoded video frames via the tiled/WC mapping for performance,
			// but the page tracker only syncs during FlushToGPU (which requires draw calls).
			// During FMV-only playback (no 3D rendering), the sync never happens otherwise.
			uint32_t overlayPhysAddr = pvideo_base + pvideo_offset;
			uint32_t overlayBytes = overlayPitch * overlayHeight;
			CxbxSyncTiledRangeToContiguous(overlayPhysAddr, overlayBytes);

			uint8_t *pOverlayData = (uint8_t *)(CONTIGUOUS_MEMORY_BASE + pvideo_base + pvideo_offset);

			// Calculate output rectangle (PVIDEO coordinates → host backbuffer)
			int out_x = GET_MASK(pvideo_point_out, NV_PVIDEO_POINT_OUT_X);
			int out_y = GET_MASK(pvideo_point_out, NV_PVIDEO_POINT_OUT_Y);
			int out_w = GET_MASK(pvideo_size_out, NV_PVIDEO_SIZE_OUT_WIDTH);
			int out_h = GET_MASK(pvideo_size_out, NV_PVIDEO_SIZE_OUT_HEIGHT);

			// Scale overlay output rect from Xbox framebuffer coords to host backbuffer coords
			DWORD XboxBackBufferWidth = g_PgraphBackBufferWidth;
			DWORD XboxBackBufferHeight = g_PgraphBackBufferHeight;
			if (XboxBackBufferWidth == 0) XboxBackBufferWidth = 640;
			if (XboxBackBufferHeight == 0) XboxBackBufferHeight = 480;

			float xScale = width / (float)XboxBackBufferWidth;
			float yScale = height / (float)XboxBackBufferHeight;
			float offsetX = (g_HostBackBufferDesc.Width - width) / 2.0f;
			float offsetY = (g_HostBackBufferDesc.Height - height) / 2.0f;

			RECT destRect;
			destRect.left = (LONG)(out_x * xScale + offsetX);
			destRect.top = (LONG)(out_y * yScale + offsetY);
			destRect.right = (LONG)((out_x + out_w) * xScale + offsetX);
			destRect.bottom = (LONG)((out_y + out_h) * yScale + offsetY);

			// Clamp to host backbuffer
			if (destRect.right > (LONG)g_HostBackBufferDesc.Width)
				destRect.right = (LONG)g_HostBackBufferDesc.Width;
			if (destRect.bottom > (LONG)g_HostBackBufferDesc.Height)
				destRect.bottom = (LONG)g_HostBackBufferDesc.Height;

			// Reallocate overlay texture only when the current one is too small
			if (overlayWidth > g_OverlayTexWidth || overlayHeight > g_OverlayTexHeight) {
				if (g_pOverlayTex) { g_pOverlayTex->Release(); g_pOverlayTex = nullptr; }

				D3D11_TEXTURE2D_DESC texDesc = {};
				texDesc.Width = overlayWidth;
				texDesc.Height = overlayHeight;
				texDesc.MipLevels = 1;
				texDesc.ArraySize = 1;
				texDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
				texDesc.SampleDesc.Count = 1;
				texDesc.Usage = D3D11_USAGE_DYNAMIC;
				texDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
				texDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

				g_pD3DDevice->CreateTexture2D(&texDesc, nullptr, &g_pOverlayTex);
				g_OverlayTexWidth = overlayWidth;
				g_OverlayTexHeight = overlayHeight;
			}

			// Check NV_PVIDEO_FORMAT_DISPLAY bit for destination color keying.
			// When enabled, the overlay should only replace framebuffer pixels
			// whose RGB matches NV_PVIDEO_COLOR_KEY (destination color key).
			// This requires reading back destination pixels — not yet implemented.
			bool colorKeyEnabled = (pvideo_format & NV_PVIDEO_FORMAT_DISPLAY) != 0;
			if (colorKeyEnabled) {
				LOG_TEST_CASE("PVIDEO destination color key enabled");
			}

			// Map texture, convert YUY2→ARGB directly into GPU memory, unmap
			if (g_pOverlayTex) {
				D3D11_MAPPED_SUBRESOURCE mapped;
				HRESULT hr = g_pD3DDeviceContext->Map(g_pOverlayTex, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
				if (SUCCEEDED(hr)) {
					const uint8_t *pSrcRow = pOverlayData;
					uint8_t *pDstRow = (uint8_t *)mapped.pData;
					for (UINT row = 0; row < overlayHeight; row++) {
						____YUY2ToARGBRow_C(pSrcRow, pDstRow, overlayWidth);
						pSrcRow += overlayPitch;
						pDstRow += mapped.RowPitch;
					}
					g_pD3DDeviceContext->Unmap(g_pOverlayTex, 0);
					RECT srcRect = { 0, 0, (LONG)overlayWidth, (LONG)overlayHeight };
					CxbxBltSurface(g_pOverlayTex, &srcRect, pHostBackBuffer, &destRect, D3DTEXF_LINEAR);
				}
			}
		}
	}

	// Render ImGui overlay
	if (g_renderbase) {
		static std::function<void(ImGuiUI*, ID3D11Texture2D*)> internal_render = &CxbxImGui_RenderD3D;
		g_renderbase->Render(internal_render, pHostBackBuffer);
	}

	pHostBackBuffer->Release();

	// Present to display.
	{
		CXBX_PROFILE_SCOPE(PROF_PRESENT_SWAP);
		CxbxPresent();
	}

	// Update FPS counter
	g_renderbase->UpdateFPSCounter();

	// Profiler: tick frame and dump timing breakdown once per second
	CxbxProfilerFrameTick();
}

void D3D11_init_pgraph_plugins()
{
	/* attach HLE Direct3D render plugins */
	pgraph_draw = D3D11_draw;
	pgraph_draw_state_update = D3D11_draw_state_update;
	pgraph_draw_clear = D3D11_draw_clear;
	pgraph_draw_patch = D3D11_draw_patch;
	pgraph_flip_stall = D3D11_flip_stall;
	pgraph_zpass_begin = D3D11_zpass_begin;
	pgraph_zpass_end = D3D11_zpass_end;
	pgraph_zpass_collect = D3D11_zpass_collect;
	pgraph_launch_transform_program = D3D11_launch_transform_program;
}

extern void pgraph_handle_method(
	NV2AState *d,
	unsigned int subchannel,
	unsigned int method,
	uint32_t parameter);

uint32_t NV2A_read_pgraph_register(const int reg)
{
	NV2AState* dev = g_NV2A->GetDeviceState();
	PGRAPHState *pg = &(dev->pgraph);
	return pg->regs[RI(reg)];
}

float *NV2A_get_vertex_attribute_value_pointer(unsigned slot)
{
	NV2AState* dev = g_NV2A->GetDeviceState();
	PGRAPHState *pg = &(dev->pgraph);

	// See CASE_16(NV097_SET_VERTEX_DATA4UB, 4) in LLE pgraph_handle_method()
	VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
	return vertex_attribute->inline_value;
}

const char *NV2AMethodToString(DWORD dwMethod)
{
	switch (dwMethod) {

#define ENUM_RANGED_ToString_N(Name, Method, Pitch, N) \
	case Name(N): return #Name "((" #N ")*" #Pitch ")";

#define ENUM_RANGED_ToString_1(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 0)

#define ENUM_RANGED_ToString_2(Name, Method, Pitch) \
	ENUM_RANGED_ToString_1(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 1)

#define ENUM_RANGED_ToString_3(Name, Method, Pitch) \
	ENUM_RANGED_ToString_2(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 2)

#define ENUM_RANGED_ToString_4(Name, Method, Pitch) \
	ENUM_RANGED_ToString_3(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 3) 

#define ENUM_RANGED_ToString_6(Name, Method, Pitch) \
	ENUM_RANGED_ToString_4(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 4) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 5)

#define ENUM_RANGED_ToString_8(Name, Method, Pitch) \
	ENUM_RANGED_ToString_6(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 6) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 7)

#define ENUM_RANGED_ToString_10(Name, Method, Pitch) \
	ENUM_RANGED_ToString_8(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 8) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 9) \

#define ENUM_RANGED_ToString_16(Name, Method, Pitch) \
	ENUM_RANGED_ToString_10(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 10) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 11) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 12) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 13) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 14) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 15)

#define ENUM_RANGED_ToString_32(Name, Method, Pitch) \
	ENUM_RANGED_ToString_16(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 16) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 17) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 18) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 19) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 20) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 21) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 22) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 23) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 24) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 25) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 26) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 27) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 28) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 29) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 30) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 31)

#define ENUM_METHOD_ToString(Name, Method) case Method: return #Name;
#define ENUM_RANGED_ToString(Name, Method, Pitch, Repeat) ENUM_RANGED_ToString_##Repeat(Name, Method, Pitch)
#define ENUM_BITFLD_Ignore(Name, Value)
#define ENUM_VALUE_Ignore(Name, Value)

	ENUM_NV2A(ENUM_METHOD_ToString, ENUM_RANGED_ToString, ENUM_BITFLD_Ignore, ENUM_VALUE_Ignore)

	default:
		return "UNLABLED";
	}
}
