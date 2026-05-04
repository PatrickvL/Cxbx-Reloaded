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
#include <algorithm> // std::min

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
static xbox::X_D3DBaseTexture s_SyntheticTextures[xbox::X_D3DTS_STAGECOUNT] = {};

// Per-stage SRV cache: avoids recreating SRVs every frame for the same resource.
// Promoted to file scope so CxbxD3D11InvalidateCachedSRVForTexture can access them.
static ID3D11Resource*           s_CachedResource[xbox::X_D3DTS_STAGECOUNT] = {};
static ID3D11ShaderResourceView* s_CachedSRV[xbox::X_D3DTS_STAGECOUNT] = {};
static D3D11_SRV_DIMENSION       s_CachedDim[xbox::X_D3DTS_STAGECOUNT] = {};

// Invalidate any cached SRV that wraps pTexture and unbind it from all PS slots.
// Must be called before binding pTexture as a UAV for a compute shader dispatch
// to eliminate SRV/UAV resource hazards that can trigger GPU TDRs.
void CxbxD3D11InvalidateCachedSRVForTexture(ID3D11Resource* pTexture)
{
	for (int stage = 0; stage < xbox::X_D3DTS_STAGECOUNT; stage++) {
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

	auto pg = &(g_NV2A->GetDeviceState()->pgraph);

	// Fast path: skip entire function if texture-related registers unchanged.
	// This avoids hash map lookups, format decoding, and SRV creation.
	{
		static uint32_t s_LastTexOff[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastTexCtl[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastTexFmt[4] = { ~0u, ~0u, ~0u, ~0u };
		bool anyChanged = false;
		for (int i = 0; i < 4; i++) {
			uint32_t off = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + i * 4)];
			uint32_t ctl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
			uint32_t fmt = pg->regs[RI(NV_PGRAPH_TEXFMT0 + i * 4)];
			if (off != s_LastTexOff[i] || ctl != s_LastTexCtl[i] || fmt != s_LastTexFmt[i]) {
				s_LastTexOff[i] = off;
				s_LastTexCtl[i] = ctl;
				s_LastTexFmt[i] = fmt;
				anyChanged = true;
			}
		}
		if (!anyChanged)
			return; // All texture state unchanged — skip expensive work
	}

	// Set the host texture for each stage
	for (int stage = 0; stage < xbox::X_D3DTS_STAGECOUNT; stage++) {
		auto pXboxBaseTexture = g_pXbox_SetTexture[stage];

		// Check PGRAPH TEXCTL0 enable bit.  When the Xbox D3D runtime disables
		// a texture stage, it writes CONTROL0 with the enable bit (bit 30) cleared.
		// We must respect this: disabled stages should not have textures bound,
		// otherwise we may create D3D11 resource hazards (e.g., the same texture
		// bound as both RTV and SRV) or sample stale data from a previous draw.
		uint32_t texCtl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];
		bool bTextureEnabled = (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) != 0;

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
		ID3D11Resource* pHostBaseTexture = nullptr;
		bool bNeedRelease = false;
		bool bIsRenderTargetTexture = false;
		uint32_t texOffset = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage * 4)];

		if (texOffset != 0) {
			// Check if this offset corresponds to a render target or depth stencil surface.
			// Exclude the current depth/stencil surface (identified by PGRAPH
			// surface_zeta.offset) — it cannot be sampled while bound as depth.
			auto pXboxSurface = CxbxLookupSurfaceByDataAddr(texOffset);
			if (pXboxSurface && texOffset != pg->surface_zeta.offset) {
				auto pHostRT = GetHostSurface(pXboxSurface, D3DUSAGE_RENDERTARGET);
				if (!pHostRT) {
					// Also try depth stencil — shadow mapping binds a depth surface as texture
					pHostRT = GetHostSurface(pXboxSurface, D3DUSAGE_DEPTHSTENCIL);
				}
				if (pHostRT) {
					pHostBaseTexture = pHostRT;
					bIsRenderTargetTexture = true;
				}
			}

			// If no Xbox surface registered (SetRenderTarget patches disabled),
			// check the PGRAPH RT cache for render targets created directly
			// from PGRAPH surface state.  This enables render-to-texture:
			// the game renders caustics/shadows to an offscreen RT, then
			// samples that RT as a texture in a later draw.
			if (!bIsRenderTargetTexture && texOffset != pg->surface_zeta.offset) {
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
					synth.Common = X_D3DCOMMON_TYPE_TEXTURE | 1; // type + refcount
					synth.Data = texOffset;
					synth.Lock = 0;
					synth.Format = pg->regs[RI(NV_PGRAPH_TEXFMT0 + stage * 4)];

					// Reconstruct the Size field for linear textures.
					// Swizzled textures use Size=0 (dimensions from Format log2 bits).
					xbox::X_D3DFORMAT xboxFmt = GetXboxPixelContainerFormat(synth.Format);
					if (EmuXBFormatIsLinear(xboxFmt)) {
						uint32_t texImageRect = pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + stage * 4)];
						uint32_t texCtl1 = pg->regs[RI(NV_PGRAPH_TEXCTL1_0 + stage * 4)];
						uint32_t width = (texImageRect >> 16) & 0x1FFF;
						uint32_t height = texImageRect & 0x1FFF;
						uint32_t pitch = (texCtl1 >> 16) & 0xFFFF;
						if (width > 0 && height > 0 && pitch >= 64)
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
			case X_D3DCOMMON_TYPE_TEXTURE:
				pHostBaseTexture = GetHostBaseTexture(pXboxBaseTexture, /*D3DUsage=*/0, stage);
				break;
			case X_D3DCOMMON_TYPE_SURFACE:
				// Surfaces can be set in the texture stages, instead of textures
				LOG_TEST_CASE("ActiveTexture set to a surface (non-texture) resource"); // Test cases : Burnout, Outrun 2006
				// We must wrap the surface before using it as a texture
				pHostBaseTexture = CxbxConvertXboxSurfaceToHostTexture(pXboxBaseTexture);
				// Release this texture (after SetTexture) when we succeeded in creating it :
				bNeedRelease = pHostBaseTexture != nullptr;
				break;
			default:
				LOG_TEST_CASE("ActiveTexture set to an unhandled resource type!");
				break;
			}

			// Read HostFormat from GetResourceCache :
			// TODO : Optimize this, as we're doing the lookup twice (once in GetHostBaseTexture, once here)
			auto key = GetHostResourceKey(pXboxBaseTexture, stage);
			auto& ResourceCache = GetResourceCache(key);
			auto it = ResourceCache.find(key);
			if (it != ResourceCache.end()) {
				g_HostTextureFormats[stage] = it->second.HostFormat;
			}
		}

		if (pHostBaseTexture != nullptr) {
			// Reuse cached SRV if the underlying resource hasn't changed
			if (s_CachedResource[stage] == pHostBaseTexture && s_CachedSRV[stage] != nullptr) {
				// SRV already cached and resource unchanged — skip rebind
				// (the PSSetShaderResources call from the last time is still in effect)
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
}

void CxbxUpdateHostTextureScaling()
{
	auto pg = &(g_NV2A->GetDeviceState()->pgraph);

	// Fast path: skip if texture format/offset/ctl/imagerect registers unchanged.
	// This avoids format decoding, division, and VS constant upload every draw.
	{
		static uint32_t s_LastFmt[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastOff[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastCtl[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastRect[4] = { ~0u, ~0u, ~0u, ~0u };
		static uint32_t s_LastSurfColor = ~0u;
		bool anyChanged = false;
		uint32_t surfColor = pg->surface_color.offset;
		if (surfColor != s_LastSurfColor) { s_LastSurfColor = surfColor; anyChanged = true; }
		for (int i = 0; i < 4; i++) {
			uint32_t fmt = pg->regs[RI(NV_PGRAPH_TEXFMT0 + i * 4)];
			uint32_t off = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + i * 4)];
			uint32_t ctl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + i * 4)];
			uint32_t rect = pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + i * 4)];
			if (fmt != s_LastFmt[i] || off != s_LastOff[i] || ctl != s_LastCtl[i] || rect != s_LastRect[i]) {
				s_LastFmt[i] = fmt; s_LastOff[i] = off; s_LastCtl[i] = ctl; s_LastRect[i] = rect;
				anyChanged = true;
			}
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
	std::array<std::array<float, 4>, xbox::X_D3DTS_STAGECOUNT> texcoordScales;
	texcoordScales.fill({ 1, 1, 1, 1 });

	for (int stage = 0; stage < xbox::X_D3DTS_STAGECOUNT; stage++) {
		// Read texture format directly from PGRAPH (authoritative, no HLE dependency)
		uint32_t texFmt = pg->regs[RI(NV_PGRAPH_TEXFMT0 + stage * 4)];
		uint32_t texOffset = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage * 4)];
		uint32_t texCtl0 = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];

		// No texture bound or disabled — skip
		bool texEnabled = (texCtl0 & (1 << 30)) != 0;
		if (!texEnabled || texOffset == 0) {
			continue;
		}

		xbox::X_D3DFORMAT XboxFormat = GetXboxPixelContainerFormat(texFmt);

		// Texcoord index. Just the texture stage unless fixed function or passthrough mode
		int texCoordIndex = stage;
		if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction
			|| g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough) {
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
		if (EmuXBFormatIsLinear(XboxFormat)) {
			// Test-case : This is often hit by the help screen in XDK samples.
			// Set scaling factor for this texture, which will be applied to
			// all texture-coordinates in the vertex shader
			// Note : Linear textures are two-dimensional at most (right?)
			// Read dimensions from PGRAPH TEXIMAGERECT (authoritative, replaces HLE reads)
			uint32_t texImageRect = pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + stage * 4)];
			float width  = (float)((texImageRect >> 16) & 0x1FFF);
			float height = (float)(texImageRect & 0x1FFF);

			// Account for MSAA when texture is the current render target (backbuffer)
			if (texOffset == pg->surface_color.offset) {
				// Test case: Max Payne 2 (bullet time)
				if (g_Xbox_MultiSampleType & xbox::X_D3DMULTISAMPLE_SAMPLING_MULTI) {
					float aaX, aaY;
					GetMultiSampleScaleRaw(aaX, aaY);
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
		if (EmuXBFormatIsDepthBuffer(XboxFormat)) {
			// Derive Z scale from PGRAPH-sourced format (no Xbox object needed)
			float zScale = 1.0f;
			switch (XboxFormat) {
				case xbox::X_D3DFMT_D16:
				case xbox::X_D3DFMT_LIN_D16:     zScale = 65535.0f;    break;
				case xbox::X_D3DFMT_D24S8:
				case xbox::X_D3DFMT_LIN_D24S8:   zScale = 16777215.0f; break;
				case xbox::X_D3DFMT_F16:
				case xbox::X_D3DFMT_LIN_F16:     zScale = 511.9375f;   break;
				case xbox::X_D3DFMT_F24S8:
				case xbox::X_D3DFMT_LIN_F24S8:   zScale = 1.0e30f;     break;
				default: break;
			}
			(*texCoordScale)[2] = zScale;
		}
	}
	// Convert texture scales to reciprocals for GPU-side multiply (cheaper than divide).
	// Upload as xboxTextureScaleRcp[4] at c214.
	std::array<std::array<float, 4>, xbox::X_D3DTS_STAGECOUNT> texcoordScaleRcp;
	for (int i = 0; i < xbox::X_D3DTS_STAGECOUNT; i++) {
		for (int j = 0; j < 4; j++) {
			texcoordScaleRcp[i][j] = 1.0f / texcoordScales[i][j];
		}
	}
	CxbxSetVertexShaderConstantF(CXBX_D3DVS_TEXTURES_SCALE_BASE, (float*)texcoordScaleRcp.data(), CXBX_D3DVS_TEXTURES_SCALE_SIZE);

	// Upload TEXCOORDINDEX remapping for the passthrough vertex shader.
	// On NV2A, the texture unit applies D3DTSS_TEXCOORDINDEX after VS output
	// interpolation. In D3D11 passthrough mode, we must do this remapping in the VS.
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

	if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction) {
		UpdateFixedFunctionVertexShaderState();
		isXboxConstants = false;
	}
	else {
		auto pg = &(g_NV2A->GetDeviceState()->pgraph);
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
	if (!g_bInPullerContext) {
		pfifo_flush_to_pgraph(g_NV2A->GetDeviceState());
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
		qemu_mutex_lock(&g_NV2A->GetDeviceState()->pgraph.pgraph_lock);
		pgraph_locked = true;
	}

	// Derive the vertex shader mode entirely from PGRAPH state.
	// g_Xbox_VertexShaderMode is ONLY written here on the render thread;
	// the game-thread patches no longer touch it, eliminating the race.
	//
	// Detection strategy:
	// - FIXED mode + CMAT ≈ identity → Passthrough (XYZRHW in fixed pipeline)
	// - FIXED mode + CMAT ≠ identity → FixedFunction (normal W*V*P transform)
	// - PROGRAM mode + VPSCL ≈ (1, ±1, ...) → Passthrough (XYZRHW via VS program)
	// - PROGRAM mode + VPSCL has large values → ShaderProgram (real VS program)
	//
	// Rationale: On Xbox, SetVertexShader(D3DFVF_XYZRHW|...) sets MODE=PROGRAM
	// and loads a trivial passthrough VS.  CMAT is NOT updated (stale from
	// previous draws), so we cannot use CMAT for PROGRAM mode.  Instead, the
	// runtime sets VPSCL to identity (1,1,1,0) since the VS outputs screen-space
	// positions directly.  Normal VS programs have VPSCL = (W/2, -H/2, zScale, 0).
	{
		PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;
		uint32_t pgraph_mode = GET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_MODE);

		if (pgraph_mode == NV097_SET_TRANSFORM_EXECUTION_MODE_MODE_PROGRAM) {
			// For PROGRAM mode, use VPSCL to distinguish passthrough from real VS.
			// XYZRHW passthrough: VPSCL.x ≈ 1, VPSCL.y ≈ ±1
			// Real VS programs:   VPSCL.x = Width/2 (≥ 4), VPSCL.y = -Height/2
			float vpscl[4];
			std::memcpy(vpscl, pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL], 16);

			if (fabsf(vpscl[0]) <= 1.5f && fabsf(vpscl[1]) <= 1.5f) {
				g_Xbox_VertexShaderMode = VertexShaderMode::Passthrough;
			} else {
				g_Xbox_VertexShaderMode = VertexShaderMode::ShaderProgram;
			}
		} else {
			// FIXED mode: distinguish passthrough from normal FF via CMAT.
			// For passthrough (XYZRHW), CMAT ≈ identity.
			// For normal FF, CMAT = World*View*Proj*Viewport with large values.
			float cmat[4][4];
			for (int row = 0; row < 4; row++)
				std::memcpy(&cmat[row][0], &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_CMAT0 + row][0], 16);

			bool isIdentity = true;
			for (int r = 0; r < 4 && isIdentity; r++) {
				for (int c = 0; c < 4 && isIdentity; c++) {
					float expected = (r == c) ? 1.0f : 0.0f;
					if (fabsf(cmat[r][c] - expected) > 0.01f)
						isIdentity = false;
				}
			}

			if (isIdentity) {
				g_Xbox_VertexShaderMode = VertexShaderMode::Passthrough;
			} else {
				g_Xbox_VertexShaderMode = VertexShaderMode::FixedFunction;
			}
		}
	}

	// Before we start, make sure our resource cache stays limited in size
	PrunePaletizedTexturesCache(); // TODO : Could we move this to Swap instead?

	// NOTE: Vertex shader must be updated before vertex declaration,
	// because D3D11 input layout creation depends on compiled VS bytecode
	CxbxUpdateHostVertexShader();

	CxbxUpdateHostVertexDeclaration();

	CxbxUpdateHostVertexShaderConstants();

	// Bind render target from PGRAPH surface offsets BEFORE viewport setup.
	// The viewport dimensions are clamped to the render target size, so the
	// correct RT must be bound first. Otherwise, if the RT switches from a
	// small offscreen target (e.g. 256x256 caustic texture) to the backbuffer
	// (640x480), GetHostRenderTargetDimensions returns the old (small) size,
	// causing the scissor rect to clip the viewport incorrectly.
	CxbxD3D11UpdateRenderTargetFromPGRAPH(&g_NV2A->GetDeviceState()->pgraph);

	// Set viewport from PGRAPH registers (authoritative).
	CxbxD3D11UpdateViewportFromPGRAPH(&g_NV2A->GetDeviceState()->pgraph);

	CxbxUpdateHostTextures();
	CxbxUpdateHostTextureScaling();

	// Pipeline state and sampler states from PGRAPH registers.
	// This replaces the former XboxRenderStates.Apply() (blend/depth/stencil/rasterizer)
	// and XboxTextureStates.Apply() (sampler configuration) which read from Xbox D3D
	// runtime memory. All state is now sourced from NV2A PGRAPH registers directly.
	{
		auto pg = &g_NV2A->GetDeviceState()->pgraph;
		CxbxD3D11UpdatePipelineStateFromPGRAPH(pg);
		CxbxD3D11UpdateSamplersFromPGRAPH(pg);
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
   	if (!g_DisablePixelShaders) {
   	   	CxbxUpdateActivePixelShader();
   	}

	// Refresh the zero-stride vertex defaults buffer with current NV2A sticky
	// attribute values before every draw, not just on vertex declaration changes.
	// This ensures non-streamed attributes (e.g. texcoords not in the vertex
	// declaration) always read the latest inline_value[] data.
	CxbxD3D11UpdateVertexDefaultsBuffer();

	// Release pgraph_lock — all PGRAPH register reads for this draw are done.
	// The puller thread is now free to process new commands for the next draw.
	if (pgraph_locked) {
		qemu_mutex_unlock(&g_NV2A->GetDeviceState()->pgraph.pgraph_lock);
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
