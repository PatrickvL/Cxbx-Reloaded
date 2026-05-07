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
#include "EmuD3D8_common.h"
#include "Backend\Backend_D3D11.h"
#include "Backend\Backend_D3D11_Profiler.h"
#include <algorithm> // std::min

// NV2A-native linear format check.  On NV2A, linear (pitch-based) textures
// use format color codes with the LU_IMAGE or LC_IMAGE prefix.
// These occupy specific ranges in the 8-bit color code field.
static inline bool IsNV2AColorFormatLinear(uint32_t colorFmt) {
	// LU_IMAGE range 1: 0x10..0x20
	if (colorFmt >= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A1R5G5B5
		&& colorFmt <= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8Y8)
		return true;
	// LC_IMAGE (YUV): 0x24..0x26
	if (colorFmt >= NV097_SET_TEXTURE_FORMAT_COLOR_LC_IMAGE_CR8YB8CB8YA8
		&& colorFmt <= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_A8CR8CB8Y8)
		return true;
	// LU_IMAGE depth: 0x2E..0x31
	if (colorFmt >= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_X8_Y24_FIXED
		&& colorFmt <= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_Y16_FLOAT)
		return true;
	// LU_IMAGE range 2: 0x35..0x40
	if (colorFmt >= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_Y16
		&& colorFmt <= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_B8G8R8A8)
		return true;
	return false;
}

// Thread-local flag: true when executing on the PFIFO puller thread.
// When set, CxbxUpdateNativeD3DResources skips pfifo_flush_to_pgraph
// because PGRAPH registers are already current (we ARE the puller).
thread_local bool g_bInPullerContext = false;

void CxbxSetPullerContext(bool active) { g_bInPullerContext = active; }

// Synthetic Xbox texture objects constructed from PGRAPH registers.
// Used when SetTexture patches are disabled: the Xbox D3D runtime writes
// texture format/offset/size to the NV2A pushbuffer, so PGRAPH has all
// the information needed to reconstruct the Xbox texture descriptor.
// One per texture stage; updated each draw by CxbxUpdateHostTextures.
static xbox::X_D3DBaseTexture s_SyntheticTextures[NV2A_MAX_TEXTURES] = {};

// Per-stage SRV cache: avoids recreating SRVs every frame for the same resource.
// Promoted to file scope so CxbxD3D11InvalidateCachedSRVForTexture can access them.
static ID3D11Resource*           s_CachedResource[NV2A_MAX_TEXTURES] = {};
static ID3D11ShaderResourceView* s_CachedSRV[NV2A_MAX_TEXTURES] = {};
static D3D11_SRV_DIMENSION       s_CachedDim[NV2A_MAX_TEXTURES] = {};

// Shared texture state generation counter — incremented by CxbxUpdateHostTextures
// when TEXOFFSET/TEXCTL0/TEXFMT change.  CxbxUpdateHostTextureScaling uses this to
// avoid redundantly re-reading the same 12 registers for its own fast path.
static uint32_t s_TextureStateGeneration = 0;

// Cached per-stage register values from the fast-path check in CxbxUpdateHostTextures.
// CxbxUpdateHostTextureScaling reuses these instead of re-reading PGRAPH.
static uint32_t s_CachedTexOff[4] = {};
static uint32_t s_CachedTexCtl[4] = {};
static uint32_t s_CachedTexFmt[4] = {};

// Force CxbxUpdateHostTextures to rebind SRVs on the next draw call.
// Must be called whenever PS SRV bindings are disturbed externally
// (e.g., blit/present unbinding slot 0, CS dispatch unbinding all slots).
static bool s_TextureSRVsDirty = false;

// Lightweight: only forces SRV rebinding on next draw without trashing
// the register cache (avoids expensive full texture re-upload).
// Clears cached resource pointers so the full lookup path runs and
// correctly resolves RT-as-texture scenarios.
void CxbxMarkTextureSRVsDirty()
{
    s_TextureSRVsDirty = true;
}

void CxbxInvalidateTextureStateCache()
{
    // Setting cached values to ~0 guarantees the fast-path check will detect a "change"
    for (int i = 0; i < 4; i++) {
        s_CachedTexOff[i] = ~0u;
        s_CachedTexCtl[i] = ~0u;
        s_CachedTexFmt[i] = ~0u;
    }
    s_TextureSRVsDirty = true;
}

// Invalidate any cached SRV that wraps pTexture and unbind it from all PS slots.
// Must be called before binding pTexture as a UAV for a compute shader dispatch
// to eliminate SRV/UAV resource hazards that can trigger GPU TDRs.
void CxbxD3D11InvalidateCachedSRVForTexture(ID3D11Resource* pTexture)
{
	for (int stage = 0; stage < NV2A_MAX_TEXTURES; stage++) {
		if (s_CachedResource[stage] == pTexture) {
			if (s_CachedSRV[stage]) {
				ID3D11ShaderResourceView* pNullSRV = nullptr;
				g_pD3DDeviceContext->PSSetShaderResources(stage, 1, &pNullSRV);
				g_pD3DDeviceContext->PSSetShaderResources(4 + stage, 1, &pNullSRV);
				g_pD3DDeviceContext->PSSetShaderResources(8 + stage, 1, &pNullSRV);
				s_CachedSRV[stage]->Release();
				s_CachedSRV[stage] = nullptr;
			}
			s_CachedResource[stage] = nullptr;
		}
	}
}

void CxbxUpdateHostTextures()
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	auto d = g_NV2A->GetDeviceState();
	auto pg = &d->pgraph;

	// Fast path: skip entire function if texture-related registers unchanged.
	// This avoids hash map lookups, format decoding, and SRV creation.
	{
		bool anyChanged = false;
		for (int i = 0; i < 4; i++) {
			uint32_t off = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + i * 4)];
			uint32_t ctl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
			uint32_t fmt = pg->regs[RI(NV_PGRAPH_TEXFMT0 + i * 4)];
			if (off != s_CachedTexOff[i] || ctl != s_CachedTexCtl[i] || fmt != s_CachedTexFmt[i]) {
				s_CachedTexOff[i] = off;
				s_CachedTexCtl[i] = ctl;
				s_CachedTexFmt[i] = fmt;
				anyChanged = true;
			}
		}
		if (!anyChanged && !s_TextureSRVsDirty)
			return; // All texture state unchanged and SRVs still bound — skip expensive work
		if (anyChanged)
			s_TextureStateGeneration++; // Signal CxbxUpdateHostTextureScaling
	}

	// Set the host texture for each stage
	for (int stage = 0; stage < NV2A_MAX_TEXTURES; stage++) {
		auto pXboxBaseTexture = g_pXbox_SetTexture[stage];

		// Check PGRAPH TEXCTL0 enable bit.  When the Xbox D3D runtime disables
		// a texture stage, it writes CONTROL0 with the enable bit (bit 30) cleared.
		// We must respect this: disabled stages should not have textures bound,
		// otherwise we may create D3D11 resource hazards (e.g., the same texture
		// bound as both RTV and SRV) or sample stale data from a previous draw.
		//
		// Exception: NV2A SHADERPROG mode overrides TEXCTL0.  When SHADERPROG
		// specifies a non-NONE mode (e.g. PROJECT2D) for a stage, the texture
		// unit IS active regardless of TEXCTL0.  This matters for point sprites
		// where the game may not explicitly enable TEXCTL0 for stage 3.
		uint32_t texCtl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];
		bool bTextureEnabled = (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) != 0;
		if (!bTextureEnabled) {
			uint32_t shaderProg = pg->regs[RI(NV_PGRAPH_SHADERPROG)];
			uint32_t stageMode = (shaderProg >> (stage * 5)) & 0x1Fu;
			if (stageMode != 0) // PS_TEXTUREMODES_NONE = 0
				bTextureEnabled = true;
		}

		if (!bTextureEnabled) {
			// Texture stage is disabled in PGRAPH — unbind and skip
			if (s_CachedSRV[stage]) {
				s_CachedSRV[stage]->Release();
				s_CachedSRV[stage] = nullptr;
			}
			s_CachedResource[stage] = nullptr;
			ID3D11ShaderResourceView* pNullSRV = nullptr;
			g_pD3DDeviceContext->PSSetShaderResources(stage, 1, &pNullSRV);
			g_pD3DDeviceContext->PSSetShaderResources(4 + stage, 1, &pNullSRV);
			g_pD3DDeviceContext->PSSetShaderResources(8 + stage, 1, &pNullSRV);
			continue;
		}

		// Read texture VRAM offset from PGRAPH — authoritative source.
		// The Xbox D3D runtime always writes SET_TEXTURE_OFFSET to the
		// pushbuffer, so PGRAPH has the physical address of the texture.
		// The offset is relative to DMA context A or B (selected by TEXFMT0 bit 1).
		ID3D11Resource* pHostBaseTexture = nullptr;
		bool bNeedRelease = false;
		bool bIsRenderTargetTexture = false;
		uint32_t texOffsetRaw = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage * 4)];

		// Resolve DMA base: CONTEXT_DMA bit in TEXFMT selects dma_a (0) or dma_b (1)
		uint32_t texFmtReg = pg->regs[RI(NV_PGRAPH_TEXFMT0 + stage * 4)];
		bool texDmaSelect = (texFmtReg & NV_PGRAPH_TEXFMT0_CONTEXT_DMA) != 0;
		uint32_t texDmaBase = NV2ADevice::ResolveDmaBaseAddress(
			d, texDmaSelect ? pg->dma_b : pg->dma_a);
		uint32_t texOffset = texDmaBase + texOffsetRaw;

		if (texOffset != 0) {
			// Check if this texture offset corresponds to a render target:
			// the game may render caustics/shadows to an offscreen RT, then
			// sample that RT as a texture in a later draw.
			// Only treat as RT-texture if the offset is NOT the currently bound
			// color surface (sampling the active RT is undefined) and is not the
			// current depth surface (can't sample while bound as DSV).
			if (texOffset != pg->regs[RI(NV_PGRAPH_BOFFSET4)]
				&& texOffset != pg->regs[RI(NV_PGRAPH_BOFFSET3)]) {
				auto pPgraphRT = CxbxLookupPgraphRTByOffset(texOffset);
				if (pPgraphRT) {
					pHostBaseTexture = pPgraphRT;
					bIsRenderTargetTexture = true;
					// D3D11 will unbind the RTV when this resource is bound as SRV.
					// Invalidate RT tracking so the next draw rebinds the RTV.
					CxbxInvalidatePgraphRTBinding();
				}
			}

			// For non-RT textures, try the texture side-map first (populated
			// by SetTexture patches if they're enabled), then fall back to
			// the HLE-tracked texture, then to constructing a synthetic
			// Xbox texture from PGRAPH registers.
			if (!bIsRenderTargetTexture) {
				auto pgTex = CxbxLookupTextureByDataAddr(texOffset);
				if (pgTex != nullptr) {
					pXboxBaseTexture = pgTex;
				} else if (pXboxBaseTexture != xbox::zeroptr
				           && pXboxBaseTexture != &s_SyntheticTextures[stage]) {
					// No side-map entry, but HLE has a genuine texture for
					// this stage (set via SetTexture/SwitchTexture patches).
					// Prefer it and register in the side-map for future lookups.
					// Exclude stale synthetic textures: their Data/Format fields
					// may belong to the previous draw's texture, not the current
					// one identified by texOffset (PGRAPH TEXOFFSET).
					CxbxRegisterTextureByDataAddr(texOffset, pXboxBaseTexture);
				} else {
					// No side-map entry: build/update synthetic X_D3DBaseTexture
					// from current PGRAPH registers.  Must re-derive every draw
					// because the same stage may bind different textures across
					// draws (e.g., ocean floor then font overlay in Dolphin).
					// The Xbox D3D runtime writes pTexture->Format directly as
					// the NV097_SET_TEXTURE_FORMAT argument, so PGRAPH TEXFMT
					// contains the exact Xbox Format field value.
					auto& synth = s_SyntheticTextures[stage];
					synth.Common = X_D3DCOMMON_TYPE_TEXTURE | X_D3DCOMMON_D3DCREATED | 1; // type + d3d-created + refcount
					synth.Data = texOffset;
					synth.Lock = 0;
					synth.Format = pg->regs[RI(NV_PGRAPH_TEXFMT0 + stage * 4)];

					// Reconstruct the Size field for linear textures.
					// Swizzled textures use Size=0 (dimensions from Format log2 bits).
					uint32_t fmtColor = GET_MASK(synth.Format, NV097_SET_TEXTURE_FORMAT_COLOR);
					if (IsNV2AColorFormatLinear(fmtColor)) {
						uint32_t texImageRect = pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + stage * 4)];
						uint32_t texCtl1 = pg->regs[RI(NV_PGRAPH_TEXCTL1_0 + stage * 4)];
						uint32_t width = (texImageRect >> 16) & 0x1FFF;
						uint32_t height = texImageRect & 0x1FFF;
						uint32_t pitch = (texCtl1 >> 16) & 0xFFFF;
						// NV2A minimum pitch alignment is 64 bytes; clamp to
						// avoid synth.Size=0 which would force swizzled decode.
						if (pitch < 64)
							pitch = 64;
						if (width > 0 && height > 0)
							synth.Size = ((width - 1) & 0xFFF)
								| (((height - 1) & 0xFFF) << X_D3DSIZE_HEIGHT_SHIFT)
								| ((((pitch / 64) - 1) & 0xFF) << X_D3DSIZE_PITCH_SHIFT);
						else
							synth.Size = 0;
					} else {
						synth.Size = 0;
					}

					pXboxBaseTexture = &synth;
					// Publish so downstream code (CxbxGetTexFmtFixup,
					// CxbxUpdateHostTextureScaling) can resolve this stage.
					g_pXbox_SetTexture[stage] = &synth;
				}
			}
		}

		if (!bIsRenderTargetTexture && pXboxBaseTexture != xbox::zeroptr) {
			DWORD XboxResourceType = GetXboxCommonResourceType(pXboxBaseTexture);
			switch (XboxResourceType) {
			case X_D3DCOMMON_TYPE_TEXTURE: {
				DXGI_FORMAT hostFormat = DXGI_FORMAT_UNKNOWN;
				pHostBaseTexture = GetHostBaseTextureWithFormat(pXboxBaseTexture, /*D3DUsage=*/0, stage, &hostFormat);
				if (hostFormat != DXGI_FORMAT_UNKNOWN)
					g_HostTextureFormats[stage] = hostFormat;
				break;
			}
			case X_D3DCOMMON_TYPE_SURFACE:
				// Surfaces can be set in the texture stages, instead of textures
				LOG_TEST_CASE("ActiveTexture set to a surface (non-texture) resource"); // Test cases : Burnout, Outrun 2006
				// For D3D11, the surface IS already the texture (ID3D11Texture2D)
				{
					ID3D11Texture2D* pHostSurface = GetHostSurface(pXboxBaseTexture);
					if (pHostSurface) {
						pHostSurface->AddRef();
						pHostBaseTexture = pHostSurface;
					} else {
						LOG_TEST_CASE("Failed to get host surface");
					}
				}
				// Release this texture (after SetTexture) when we succeeded in creating it :
				bNeedRelease = pHostBaseTexture != nullptr;
				break;
			default:
				LOG_TEST_CASE("ActiveTexture set to an unhandled resource type!");
				break;
			}
		}

		if (pHostBaseTexture != nullptr) {
			// Reuse cached SRV if the underlying resource hasn't changed
			if (s_CachedResource[stage] == pHostBaseTexture && s_CachedSRV[stage] != nullptr) {
				// SRV already cached and resource unchanged — rebind only if externally dirtied
				if (s_TextureSRVsDirty) {
					g_pD3DDeviceContext->PSSetShaderResources(stage, 1, &s_CachedSRV[stage]);
					if (s_CachedDim[stage] == D3D11_SRV_DIMENSION_TEXTURE3D)
						g_pD3DDeviceContext->PSSetShaderResources(4 + stage, 1, &s_CachedSRV[stage]);
					else if (s_CachedDim[stage] == D3D11_SRV_DIMENSION_TEXTURECUBE)
						g_pD3DDeviceContext->PSSetShaderResources(8 + stage, 1, &s_CachedSRV[stage]);
				}
			} else {
				// Release old cached SRV
				if (s_CachedSRV[stage]) {
					s_CachedSRV[stage]->Release();
					s_CachedSRV[stage] = nullptr;
				}
				s_CachedResource[stage] = nullptr;

				// Create a shader resource view for the texture
				D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
				D3D11_RESOURCE_DIMENSION dim;
				pHostBaseTexture->GetType(&dim);

				switch (dim) {
				case D3D11_RESOURCE_DIMENSION_TEXTURE2D: {
					D3D11_TEXTURE2D_DESC texDesc = {};
					((ID3D11Texture2D*)pHostBaseTexture)->GetDesc(&texDesc);
					// Depth textures use typeless format — map to SRV-compatible format for sampling
					srvDesc.Format = IsDepthFormat(texDesc.Format) ? GetDepthSRVFormat(texDesc.Format) : texDesc.Format;
					if (texDesc.ArraySize == 6) {
						srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURECUBE;
						srvDesc.TextureCube.MipLevels = texDesc.MipLevels;
						srvDesc.TextureCube.MostDetailedMip = 0;
					} else {
						srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
						srvDesc.Texture2D.MipLevels = texDesc.MipLevels;
						srvDesc.Texture2D.MostDetailedMip = 0;
					}
					break;
				}
				case D3D11_RESOURCE_DIMENSION_TEXTURE3D: {
					D3D11_TEXTURE3D_DESC texDesc = {};
					((ID3D11Texture3D*)pHostBaseTexture)->GetDesc(&texDesc);
					srvDesc.Format = texDesc.Format;
					srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE3D;
					srvDesc.Texture3D.MipLevels = texDesc.MipLevels;
					srvDesc.Texture3D.MostDetailedMip = 0;
					break;
				}
				default:
					// Unsupported resource type
					if (bNeedRelease) pHostBaseTexture->Release();
					continue;
				}

				ID3D11ShaderResourceView* pSRV = nullptr;
				HRESULT hRet = g_pD3DDevice->CreateShaderResourceView(pHostBaseTexture, &srvDesc, &pSRV);
				DEBUG_D3DRESULT(hRet, "g_pD3DDevice->CreateShaderResourceView");
				if (FAILED(hRet)) {
					if (hRet == DXGI_ERROR_DEVICE_REMOVED) {
						HRESULT reason = g_pD3DDevice->GetDeviceRemovedReason();
						CxbxrAbort("D3D11 device removed (DXGI_ERROR_DEVICE_REMOVED).\n"
							"Reason: 0x%08X\n\n"
							"This is usually caused by a GPU driver crash (TDR) triggered by\n"
							"an invalid compute shader dispatch or resource hazard.\n"
							"Enable D3D11 debug layer for more details.", reason);
					}
					EmuLog(LOG_LEVEL::WARNING, "CxbxUpdateHostTextures : g_pD3DDevice->CreateShaderResourceView "
						"D3D error (0x%08X: format=%u)", hRet, srvDesc.Format);
				}

				if (SUCCEEDED(hRet) && pSRV != nullptr) {
					s_CachedResource[stage] = pHostBaseTexture;
					s_CachedSRV[stage] = pSRV; // Keep ref for cache
					s_CachedDim[stage] = srvDesc.ViewDimension;
					// Always bind to the base slot (for compiled PS path)
					g_pD3DDeviceContext->PSSetShaderResources(stage, 1, &pSRV);
					// All pixel shaders use separate Texture2D/3D/Cube declarations
					// at t0-3/t4-7/t8-11; bind to the type-appropriate slot too
					if (srvDesc.ViewDimension == D3D11_SRV_DIMENSION_TEXTURE3D)
						g_pD3DDeviceContext->PSSetShaderResources(4 + stage, 1, &pSRV);
					else if (srvDesc.ViewDimension == D3D11_SRV_DIMENSION_TEXTURECUBE)
						g_pD3DDeviceContext->PSSetShaderResources(8 + stage, 1, &pSRV);
				}
			}
		} else {
			// Clear cache and unbind
			if (s_CachedSRV[stage]) {
				s_CachedSRV[stage]->Release();
				s_CachedSRV[stage] = nullptr;
			}
			s_CachedResource[stage] = nullptr;
			ID3D11ShaderResourceView* pNullSRV = nullptr;
			g_pD3DDeviceContext->PSSetShaderResources(stage, 1, &pNullSRV);
			g_pD3DDeviceContext->PSSetShaderResources(4 + stage, 1, &pNullSRV);
			g_pD3DDeviceContext->PSSetShaderResources(8 + stage, 1, &pNullSRV);
		}
		if (bNeedRelease) {
			pHostBaseTexture->Release();
		}
	}
	s_TextureSRVsDirty = false;
}

void CxbxUpdateHostTextureScaling()
{
	auto d = g_NV2A->GetDeviceState();
	auto pg = &d->pgraph;

	// Fast path: skip if texture state hasn't changed since last call.
	// CxbxUpdateHostTextures (called immediately before us) already checks
	// TEXOFFSET/TEXCTL0/TEXFMT and bumps s_TextureStateGeneration on change.
	// We only need to additionally check TEXIMAGERECT and surface_color.offset.
	{
		static uint32_t s_LastTexGen = ~0u;
		static uint32_t s_LastRect[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastSurfColor = ~0u;
		bool anyChanged = false;
		if (s_TextureStateGeneration != s_LastTexGen) {
			s_LastTexGen = s_TextureStateGeneration;
			anyChanged = true;
		}
		uint32_t surfColor = pg->regs[RI(NV_PGRAPH_BOFFSET3)];
		if (surfColor != s_LastSurfColor) { s_LastSurfColor = surfColor; anyChanged = true; }
		for (int i = 0; i < 4; i++) {
			uint32_t rect = pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + i * 4)];
			if (rect != s_LastRect[i]) { s_LastRect[i] = rect; anyChanged = true; }
		}
		if (!anyChanged) return;
	}

	// Xbox works with "Linear" and "Swizzled" texture formats
	// Linear formats are not addressed with normalized coordinates (similar to https://www.khronos.org/opengl/wiki/Rectangle_Texture?)
	// We want to use normalized coordinates in our shaders, so need to be able to scale the coordinates back
	// Note texcoords aren't only used for texture lookups
	// TODO store scaling per texture instead of per stage, and scale during lookup in the pixel shader

	// Each texture stage has one texture coordinate set associated with it
	// We'll store scale factors for each texture coordinate set
	std::array<std::array<float, 4>, NV2A_MAX_TEXTURES> texcoordScales;
	texcoordScales.fill({ 1, 1, 1, 1 });

	for (int stage = 0; stage < NV2A_MAX_TEXTURES; stage++) {
		// Reuse cached register values from CxbxUpdateHostTextures (avoids re-reading PGRAPH)
		uint32_t texFmt = s_CachedTexFmt[stage];
		uint32_t texOffset = s_CachedTexOff[stage];
		uint32_t texCtl0 = s_CachedTexCtl[stage];

		// No texture bound or disabled — skip
		bool texEnabled = (texCtl0 & (1 << 30)) != 0;
		if (!texEnabled || texOffset == 0) {
			continue;
		}

		uint32_t colorFmt = GET_MASK(texFmt, NV097_SET_TEXTURE_FORMAT_COLOR);

		// Texcoord index. Just the texture stage unless fixed function mode
		int texCoordIndex = stage;
		if (NV2AIsFixedFunctionMode(pg)) {
			// Read texgen mode from PGRAPH CSV1_A/CSV1_B to determine if
			// coordinates are generated (no HLE dependency).
			unsigned int csvReg = (stage < 2) ? NV_PGRAPH_CSV1_A : NV_PGRAPH_CSV1_B;
			unsigned int sMask  = (stage % 2) ? NV_PGRAPH_CSV1_A_T1_S : NV_PGRAPH_CSV1_A_T0_S;
			uint32_t texgenS = GET_MASK(pg->regs[RI(csvReg)], sMask);

			// If coordinates are generated, we don't have to worry about the coordinates coming from the title
			bool isGenerated = (texgenS != NV_PGRAPH_CSV1_A_T0_S_DISABLE);
			if (isGenerated) {
				continue;
			}

			// On NV2A, texcoord routing is identity for FF (stage i uses TEXCOORD i)
			texCoordIndex = stage;
		}

		auto texCoordScale = &texcoordScales[texCoordIndex];

		// Check for active linear textures.
		// NV2A linear formats use LU_IMAGE/LC_IMAGE color codes (ranges 0x10-0x20, 0x24-0x26, 0x2E-0x31, 0x35-0x40).
		if (IsNV2AColorFormatLinear(colorFmt)) {
			// Test-case : This is often hit by the help screen in XDK samples.
			// Set scaling factor for this texture, which will be applied to
			// all texture-coordinates in the vertex shader
			// Note : Linear textures are two-dimensional at most (right?)
			// Read dimensions from PGRAPH TEXIMAGERECT (authoritative, replaces HLE reads)
			uint32_t texImageRect = pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + stage * 4)];
			float width  = (float)((texImageRect >> 16) & 0x1FFF);
			float height = (float)(texImageRect & 0x1FFF);

			// Account for MSAA when texture is the current render target (backbuffer)
			if (texOffset == pg->regs[RI(NV_PGRAPH_BOFFSET3)]) {
				// Test case: Max Payne 2 (bullet time)
				// Use PGRAPH anti_aliasing directly instead of g_Xbox_MultiSampleType HLE global
				auto surf = NV2AGetSurfaceState(pg);
				if (surf.antiAliasing != NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_1) {
					float aaX = 1.0f, aaY = 1.0f;
					switch (surf.antiAliasing) {
					case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_CENTER_CORNER_2:
						aaX = 2.0f; break;
					case NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_SQUARE_OFFSET_4:
						aaX = 2.0f; aaY = 2.0f; break;
					}
					width /= aaX;
					height /= aaY;
				}
			}

			*texCoordScale = {
				width,
				height,
				1.0f, // TODO should this be mip levels for volume textures?
				1.0f
			};
		}

		// When a depth buffer is used as a texture
		// We do 'Native Shadow Mapping'
		// https://aras-p.info/texts/D3D9GPUHacks.html
		// The z texture coordinate component holds a depth value, which needs to be normalized
		// TODO implement handling for
		// - X_D3DRS_SHADOWFUNC
		// - X_D3DRS_POLYGONOFFSETZSLOPESCALE
		// - X_D3DRS_POLYGONOFFSETZOFFSET
		// NV2A depth texture format codes are 0x2A-0x31 (swizzled and linear depth)
		if (colorFmt >= NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_X8_Y24_FIXED
			&& colorFmt <= NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_Y16_FLOAT) {
			// Derive Z scale from the NV2A texture format color code
			float zScale = 1.0f;
			switch (colorFmt) {
				case NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_Y16_FIXED:
				case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_Y16_FIXED:
					zScale = 65535.0f;    break;
				case NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_X8_Y24_FIXED:
				case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_X8_Y24_FIXED:
					zScale = 16777215.0f; break;
				case NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_Y16_FLOAT:
				case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_Y16_FLOAT:
					zScale = 511.9375f;   break;
				case NV097_SET_TEXTURE_FORMAT_COLOR_SZ_DEPTH_X8_Y24_FLOAT:
				case NV097_SET_TEXTURE_FORMAT_COLOR_LU_IMAGE_DEPTH_X8_Y24_FLOAT:
					zScale = 1.0e30f;     break;
				default: break;
			}
			(*texCoordScale)[2] = zScale;
		}
	}
	// Convert texture scales to reciprocals for GPU-side multiply (cheaper than divide).
	// Upload as xboxTextureScaleRcp[4] at c214.
	std::array<std::array<float, 4>, NV2A_MAX_TEXTURES> texcoordScaleRcp;
	for (int i = 0; i < NV2A_MAX_TEXTURES; i++) {
		for (int j = 0; j < 4; j++) {
			texcoordScaleRcp[i][j] = 1.0f / texcoordScales[i][j];
		}
	}
	CxbxSetVertexShaderConstantF(CXBX_D3DVS_TEXTURES_SCALE_BASE, (float*)texcoordScaleRcp.data(), CXBX_D3DVS_TEXTURES_SCALE_SIZE);

	// Upload TEXCOORDINDEX remapping for the vertex shader output footer.
	// On NV2A, the texture unit applies D3DTSS_TEXCOORDINDEX after VS output
	// interpolation. In D3D11, we must do this remapping in the VS.
	// On NV2A, texcoord routing is always identity (stage i uses TEXCOORD i),
	// so we always upload the identity mapping.
	{
		float defaultIndices[4] = { 0.0f, 1.0f, 2.0f, 3.0f };
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_TEXCOORDINDEX, defaultIndices, 1);
	}
}

void CxbxUpdateDirtyVertexShaderConstants(const float* constants, bool* dirty) {
	// Reduce the number of calls by updating contiguous "batches" of dirty states
	int batchStartIndex = -1; // -1 means we aren't in a batch

	for (int i = 0; i < X_D3DVS_CONSTREG_COUNT; i++) {
		if (batchStartIndex == -1 && dirty[i]) {
			batchStartIndex = i; // Start a batch
		}
		else if (batchStartIndex != -1 && !dirty[i]) {
			// Finish the batch
			int count = i - batchStartIndex;
			CxbxSetVertexShaderConstantF(batchStartIndex, &constants[batchStartIndex * 4], count);
			batchStartIndex = -1;
		}

		// Constant is no longer dirty
		dirty[i] = false;
	}

	// Send the final batch
	if (batchStartIndex != -1) {
		int count = X_D3DVS_CONSTREG_COUNT - batchStartIndex;
		CxbxSetVertexShaderConstantF(batchStartIndex, &constants[batchStartIndex * 4], count);
	}
}

// D3DDevice_SetVertexShaderConstant patches have been removed;
// Xbox native code pushes NV097_SET_TRANSFORM_CONSTANT through PFIFO → PGRAPH.
void CxbxUpdateHostVertexShaderConstants()
{
	// Track which constants are currently written
	// So we can skip updates
	static bool isXboxConstants = false;
	auto pg = &(g_NV2A->GetDeviceState()->pgraph);

	if (NV2AIsFixedFunctionMode(pg)) {
		UpdateFixedFunctionVertexShaderState();
		isXboxConstants = false;
	}
	else {
		auto constant_floats = (float*)pg->vsh_constants;

		if (isXboxConstants) {
			CxbxUpdateDirtyVertexShaderConstants(constant_floats, pg->vsh_constants_dirty);
		}
		else {
			CxbxSetVertexShaderConstantF(0, constant_floats, X_D3DVS_CONSTREG_COUNT);
		}

		isXboxConstants = true;
		CxbxUpdateHostViewPortOffsetAndScaleConstants();
	}

	// Upload NV2A fog parameters from PGRAPH registers.
	// Only update if fog-related registers changed.
	{
		auto *pg = &g_NV2A->GetDeviceState()->pgraph;
		uint32_t ctl3 = pg->regs[RI(NV_PGRAPH_CONTROL_3)];
		uint32_t fogP0 = pg->regs[RI(NV_PGRAPH_FOGPARAM0)];
		uint32_t fogP1 = pg->regs[RI(NV_PGRAPH_FOGPARAM1)];

		static uint32_t s_LastFogCtl3 = ~0u;
		static uint32_t s_LastFogP0 = ~0u;
		static uint32_t s_LastFogP1 = ~0u;

		if (ctl3 != s_LastFogCtl3 || fogP0 != s_LastFogP0 || fogP1 != s_LastFogP1) {
			s_LastFogCtl3 = ctl3;
			s_LastFogP0 = fogP0;
			s_LastFogP1 = fogP1;
			float fogMode = (float)GET_MASK(ctl3, NV_PGRAPH_CONTROL_3_FOG_MODE);
			float fogParam0; std::memcpy(&fogParam0, &fogP0, sizeof(float));
			float fogParam1; std::memcpy(&fogParam1, &fogP1, sizeof(float));
			float fogStuff[4] = { fogMode, fogParam0, fogParam1, 0.0f };
			CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_FOGINFO, fogStuff, 1);
		}
	}
}

extern void CxbxUpdateHostVertexDeclaration(); // TMP glue
extern void CxbxUpdateHostVertexShader(); // TMP glue

void CxbxUpdateNativeD3DResources()
{
	// Drain all pending pushbuffer commands so PGRAPH regs[] are current.
	// This closes the race between the async PFIFO puller and the HLE
	// interpreters that read register state at draw time.
	// Skip when called from the puller thread itself (registers are
	// already current, and calling flush would deadlock on pfifo_lock).
	{
		CXBX_PROFILE_SCOPE(PROF_PFIFO_FLUSH);
		if (!g_bInPullerContext) {
			pfifo_flush_to_pgraph(g_NV2A->GetDeviceState());
		}
	}

	// Hold pgraph_lock while reading PGRAPH registers for this draw.
	// After pfifo_flush_to_pgraph returns, the puller thread is free to
	// process NEW commands pushed by the game thread.  Without this lock,
	// the puller can overwrite PGRAPH registers (VS constants, combiner
	// state, viewport, etc.) mid-draw-setup, causing intermittent flicker
	// (e.g., dolphin drawn at wrong position, seafloor going black).
	// The lock is released after all PGRAPH reads and before the D3D11 draw.
	bool pgraph_locked = false;
	if (!g_bInPullerContext) {
		CXBX_PROFILE_SCOPE(PROF_PGRAPH_LOCK_WAIT);
		qemu_mutex_lock(&g_NV2A->GetDeviceState()->pgraph.pgraph_lock);
		pgraph_locked = true;
	}

	// Single pg pointer for the entire per-draw state update sequence.
	PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;

	// Before we start, make sure our resource cache stays limited in size
	PrunePaletizedTexturesCache(); // TODO : Could we move this to Swap instead?

	// NOTE: Vertex shader must be updated before vertex declaration,
	// because D3D11 input layout creation depends on compiled VS bytecode
	CxbxUpdateHostVertexShader();

	CxbxUpdateHostVertexDeclaration();

	{
		CXBX_PROFILE_SCOPE(PROF_VS_CONSTANTS);
		CxbxUpdateHostVertexShaderConstants();
	}

	// Bind render target from PGRAPH surface offsets BEFORE viewport setup.
	// The viewport dimensions are clamped to the render target size, so the
	// correct RT must be bound first. Otherwise, if the RT switches from a
	// small offscreen target (e.g. 256x256 caustic texture) to the backbuffer
	// (640x480), GetHostRenderTargetDimensions returns the old (small) size,
	// causing the scissor rect to clip the viewport incorrectly.
	{
		CXBX_PROFILE_SCOPE(PROF_RENDER_TARGET);
		CxbxD3D11UpdateRenderTargetFromPGRAPH(pg);
	}

	// Set viewport from PGRAPH registers (authoritative).
	{
		CXBX_PROFILE_SCOPE(PROF_VIEWPORT);
		CxbxD3D11UpdateViewportFromPGRAPH(pg);
	}

	{
		CXBX_PROFILE_SCOPE(PROF_TEXTURES);
		CxbxUpdateHostTextures();
		CxbxUpdateHostTextureScaling();
	}

	// Pipeline state and sampler states from PGRAPH registers.
	// This replaces the former XboxRenderStates.Apply() (blend/depth/stencil/rasterizer)
	// and XboxTextureStates.Apply() (sampler configuration) which read from Xbox D3D
	// runtime memory. All state is now sourced from NV2A PGRAPH registers directly.
	{
		CXBX_PROFILE_SCOPE(PROF_PIPELINE_STATE);
		CxbxD3D11UpdatePipelineStateFromPGRAPH(pg);
	}
	{
		CXBX_PROFILE_SCOPE(PROF_SAMPLERS);
		CxbxD3D11UpdateSamplersFromPGRAPH(pg);
	}
	{
		extern float g_fLineWidth;
		g_fLineWidth = pg->line_width;
	}

	// Point sprite texture swap: NV2A uses stage 3 for point sprite textures.
	// Copy the SRV from slot 3 to slot 0 so the GS-generated UVs on TEXCOORD0
	// sample the correct texture.
	extern bool g_bPointSpriteEnabled;
	if (g_bPointSpriteEnabled) {
		ID3D11ShaderResourceView* pSRV = nullptr;
		g_pD3DDeviceContext->PSGetShaderResources(3, 1, &pSRV);
		g_pD3DDeviceContext->PSSetShaderResources(0, 1, &pSRV);
		if (pSRV) pSRV->Release();
	}

	// If Pixel Shaders are not disabled, process them
	{
		CXBX_PROFILE_SCOPE(PROF_PIXEL_SHADER);
		if (!g_DisablePixelShaders) {
			CxbxUpdateActivePixelShader();
		}
	}

	// Note: Vertex defaults upload (NV2A sticky attribute values) is handled
	// internally by CxbxD3D11VertexFetchDraw's UploadVertexDefaults() at draw
	// time. The IA-path CxbxD3D11UpdateVertexDefaultsBuffer() is no longer needed
	// since all draws use the vertex fetch CS path (SV_VertexID-based fetching).

	// Release pgraph_lock — all PGRAPH register reads for this draw are done.
	// The puller thread is now free to process new commands for the next draw.
	if (pgraph_locked) {
		qemu_mutex_unlock(&pg->pgraph_lock);
		pgraph_locked = false;
	}

	// Apply any pending D3D11 state object changes before drawing
	CxbxD3D11ApplyDirtyStates();
}

// This function should be called in tight idle-wait loops.
// It's purpose is to lower CPU cost in such a way that the
// caller will still repond quickly, without actually waiting
// or giving up it's time-slice.
// See https://docs.microsoft.com/en-us/windows/win32/api/winnt/nf-winnt-yieldprocessor
// and https://software.intel.com/en-us/cpp-compiler-developer-guide-and-reference-pause-intrinsic
inline void CxbxCPUIdleWait() // TODO : Apply wherever applicable
{
	YieldProcessor();
}

// This function indicates whether Cxbx can flush host GPU commands.
bool CxbxCanFlushHostGPU()
{
	return (g_pHostQueryWaitForIdle != nullptr);
}

// Wait until host GPU finished processing it's command queue
bool CxbxFlushHostGPU()
{
	// The following can only work when host GPU queries are available
	if (!CxbxCanFlushHostGPU()) {
		// If we can't query host GPU, return failure
		return false;
	}

	// Add an end marker to the command buffer queue.
	// This, so that the next GetData will always have at least one
	// final query event to flush out, after which GPU will be done.
	CxbxQueryIssueEnd(g_pHostQueryWaitForIdle);

	// Empty the command buffer and wait until host GPU is idle.
	BOOL queryData = FALSE;
	while (CxbxQueryGetData(g_pHostQueryWaitForIdle, &queryData, sizeof(queryData), 0) == S_FALSE)
		CxbxCPUIdleWait();

	// Signal caller that host GPU has been flushed
	return true;
}

// CxbxHandleXboxCallbacks and CxbxImpl_InsertCallback — removed.
// Native InsertCallback pushes NV097_NO_OPERATION(param) to the push buffer.
// PGRAPH raises INTR_ERROR → miniport ISR reads TRAPPED_DATA_LOW → dispatches callback.

// ******************************************************************
// * patch: D3DDevice_SetPixelShader
// ******************************************************************
