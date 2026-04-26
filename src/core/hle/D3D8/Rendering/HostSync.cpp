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
#include <algorithm> // std::min

// Thread-local flag: true when executing on the PFIFO puller thread.
// When set, CxbxUpdateNativeD3DResources skips pfifo_flush_to_pgraph
// because PGRAPH registers are already current (we ARE the puller).
thread_local bool g_bInPullerContext = false;

void CxbxSetPullerContext(bool active) { g_bInPullerContext = active; }

static std::queue<s_Xbox_Callback> g_Xbox_CallbackQueue;

void CxbxUpdateHostTextures()
{
	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	// Per-stage SRV cache: avoids recreating SRVs every frame for the same resource
	static ID3D11Resource*           s_CachedResource[xbox::X_D3DTS_STAGECOUNT] = {};
	static ID3D11ShaderResourceView* s_CachedSRV[xbox::X_D3DTS_STAGECOUNT] = {};
	static D3D11_SRV_DIMENSION       s_CachedDim[xbox::X_D3DTS_STAGECOUNT] = {};

	// Set the host texture for each stage
	for (int stage = 0; stage < xbox::X_D3DTS_STAGECOUNT; stage++) {
		auto pXboxBaseTexture = g_pXbox_SetTexture[stage];

		// Check PGRAPH TEXCTL0 enable bit.  When the Xbox D3D runtime disables
		// a texture stage, it writes CONTROL0 with the enable bit (bit 30) cleared.
		// We must respect this: disabled stages should not have textures bound,
		// otherwise we may create D3D11 resource hazards (e.g., the same texture
		// bound as both RTV and SRV) or sample stale data from a previous draw.
		bool bTextureEnabled = true; // default enabled for non-PGRAPH path
		if (g_NV2A) {
			auto pg = &(g_NV2A->GetDeviceState()->pgraph);
			uint32_t texCtl = pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];
			bTextureEnabled = (texCtl & NV_PGRAPH_TEXCTL0_0_ENABLE) != 0;
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

		// Authoritative: read the texture VRAM offset from PGRAPH
		// registers and resolve to an Xbox texture via the side-map
		// populated by SetTexture/SwitchTexture patches.  This covers
		// all cases including direct pushbuffer writes and inlined LTCG
		// code (e.g. XDK CXBFont/CXBHelp) — the Xbox D3D runtime always
		// writes SET_TEXTURE_OFFSET to the pushbuffer, so PGRAPH has it.
		// Fallback: when PGRAPH texOffset is 0, keep pXboxBaseTexture
		// from g_pXbox_SetTexture[stage] (set by HLE patches).  The IVB
		// (Begin/End) path draws synchronously before the puller thread
		// processes the push buffer, so PGRAPH may not yet have the
		// texture offset even though SetTexture was called.
		ID3D11Resource* pHostBaseTexture = nullptr;
		bool bNeedRelease = false;
		bool bIsRenderTargetTexture = false;
		if (g_NV2A) {
			auto pg = &(g_NV2A->GetDeviceState()->pgraph);
			uint32_t texOffset = pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage * 4)];
			if (texOffset != 0) {
				// When a render target is used as a texture, the TEXOFFSET
				// contains the surface's Data address.  The texture side-map
				// (populated by SetTexture) may contain a DIFFERENT Xbox
				// object than the surface side-map (populated by SetRenderTarget).
				// GetHostBaseTexture(pXboxTexture) creates a separate D3D11
				// texture from Xbox memory, which has stale/uninitialized data
				// for a render target.  Instead, look up the surface side-map
				// and use GetHostSurface to get the D3D11 texture that holds
				// the actual rendered content (created with both RT and SRV flags).
				auto pXboxSurface = CxbxLookupSurfaceByDataAddr(texOffset);
				if (pXboxSurface && pXboxSurface != g_pXbox_DepthStencil) {
					// This surface is a known render target (not depth stencil).
					// Use the host RT texture which has the rendered content.
					auto pHostRT = GetHostSurface(pXboxSurface, D3DUSAGE_RENDERTARGET);
					if (pHostRT) {
						pHostBaseTexture = pHostRT;
						bIsRenderTargetTexture = true;
					}
				}

				// Fallback: use the texture side-map for non-RT textures
				if (!bIsRenderTargetTexture) {
					auto pgTex = CxbxLookupTextureByDataAddr(texOffset);
					if (pgTex != nullptr)
						pXboxBaseTexture = pgTex;
				}
			}
			// When texOffset == 0 and g_pXbox_SetTexture[stage] is also null,
			// pXboxBaseTexture stays zeroptr — the texture will be unbound.
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
				g_pD3DDeviceContext->PSSetShaderResources(stage, 1, &s_CachedSRV[stage]);
				// All pixel shaders use separate Texture2D/3D/Cube declarations
				// at t0-3/t4-7/t8-11; bind to the type-appropriate slot too
				if (s_CachedDim[stage] == D3D11_SRV_DIMENSION_TEXTURE3D)
					g_pD3DDeviceContext->PSSetShaderResources(4 + stage, 1, &s_CachedSRV[stage]);
				else if (s_CachedDim[stage] == D3D11_SRV_DIMENSION_TEXTURECUBE)
					g_pD3DDeviceContext->PSSetShaderResources(8 + stage, 1, &s_CachedSRV[stage]);
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
					srvDesc.Format = texDesc.Format;
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
		auto pXboxBaseTexture = g_pXbox_SetTexture[stage];

		// No texture, no scaling to do
		if (pXboxBaseTexture == xbox::zeroptr) {
			continue;
		}

		// Texcoord index. Just the texture stage unless fixed function or passthrough mode
		int texCoordIndex = stage;
		if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction
			|| g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough) {
			// Get TEXCOORDINDEX for the current texture stage's state
			// Stores both the texture stage index and information for generating coordinates
			// See D3DTSS_TEXCOORDINDEX
			auto texCoordIndexState = XboxTextureStates.Get(stage, xbox::X_D3DTSS_TEXCOORDINDEX);

			// If coordinates are generated, we don't have to worry about the coordinates coming from the title
			bool isGenerated = texCoordIndexState >= X_D3DTSS_TCI_CAMERASPACENORMAL;
			if (isGenerated) {
				continue;
			}

			// Determine the texture coordinate addressing this texture stage
			texCoordIndex = (texCoordIndexState & 0x3); // 0 - 3
		}

		auto texCoordScale = &texcoordScales[texCoordIndex];

		// Check for active linear textures.
		xbox::X_D3DFORMAT XboxFormat = GetXboxPixelContainerFormat(pXboxBaseTexture);
		if (EmuXBFormatIsLinear(XboxFormat)) {
			// Test-case : This is often hit by the help screen in XDK samples.
			// Set scaling factor for this texture, which will be applied to
			// all texture-coordinates in CxbxVertexShaderTemplate.hlsl
			// Note : Linear textures are two-dimensional at most (right?)
			float width, height;
			if ((xbox::X_D3DSurface*)pXboxBaseTexture == g_pXbox_BackBufferSurface) {
				// Account for MSAA
				// Test case: Max Payne 2 (bullet time)
				GetBackBufferPixelDimensions(width, height);
			}
			else {
				width = (float)GetPixelContainerWidth(pXboxBaseTexture);
				height = (float)GetPixelContainerHeight(pXboxBaseTexture);
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
			(*texCoordScale)[2] = (float)GetZScaleForPixelContainer(pXboxBaseTexture);
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
	if (g_Xbox_VertexShaderMode == VertexShaderMode::Passthrough) {
		float texCoordIndices[4];
		for (int stage = 0; stage < xbox::X_D3DTS_STAGECOUNT; stage++) {
			auto texCoordIndexState = XboxTextureStates.Get(stage, xbox::X_D3DTSS_TEXCOORDINDEX);
			texCoordIndices[stage] = (float)(texCoordIndexState & 0x3); // 0 - 3
		}
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_TEXCOORDINDEX, texCoordIndices, 1);
	} else {
		// Default: each stage uses its own texcoord set (identity mapping)
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

// TODO : Once we're able to flush the NV2A push buffer
// remove our patches on D3DDevice_SetVertexShaderConstant (and CxbxImpl_SetVertexShaderConstant)
void CxbxUpdateHostVertexShaderConstants()
{
	// For Xbox vertex shader programs, the Xbox vertex shader constants
	// are mirrored on the host.
	// Otherwise, the same set of constants is used for the fixed function vertex shader
	// implementation instead

	// Track which constants are currently written
	// So we can skip updates
	static bool isXboxConstants = false;

	if (g_Xbox_VertexShaderMode == VertexShaderMode::FixedFunction) {
		// Write host FF shader state
		// TODO dirty tracking like for Xbox constants?
		UpdateFixedFunctionVertexShaderState();
		isXboxConstants = false;
	}
	else {
		// Write Xbox constants
		auto pg = &(g_NV2A->GetDeviceState()->pgraph);
		auto constant_floats = (float*)pg->vsh_constants;

		if (isXboxConstants) {
			// Only need to overwrite what's changed
			CxbxUpdateDirtyVertexShaderConstants(constant_floats, pg->vsh_constants_dirty);
		}
		else {
			// We need to update everything
			CxbxSetVertexShaderConstantF(0, constant_floats, X_D3DVS_CONSTREG_COUNT);
		}

		// We've written the Xbox constants
		isXboxConstants = true;

		// FIXME our viewport constants don't match Xbox values
		// If we write them to pgraph constants, like we do with constants set by the title,
		// the Xbox could overwrite them (at any time?) and we get flickering geometry.
		// For now, set our viewport constants directly in the call below,
		// overwriting whatever was in pgraph
		// Test case:
		// Xbox dashboard (during initial fade from black)
		// Need for Speed: Hot Pursuit 2 (car select)
		CxbxUpdateHostViewPortOffsetAndScaleConstants();
	}

	// Placed this here until we find a better place
	float fogTableMode = static_cast<float>(XboxRenderStates.GetXboxRenderState(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGTABLEMODE));
	// When table fog is active, check PGRAPH for _ABS fog mode variants.
	// NV2A PGRAPH fog modes 4/5/7 apply abs() to the computed fog factor;
	// bit 2 of the PGRAPH fog mode field is the _ABS flag.
	if (fogTableMode > 0.0f && g_NV2A) {
		auto *pg = &g_NV2A->GetDeviceState()->pgraph;
		if (pg->regs[RI(NV_PGRAPH_CONTROL_3)] & 0x00040000u) { // bit 2 of FOG_MODE field
			fogTableMode += 4.0f; // Promote to _ABS variant (5=EXP_ABS, 6=EXP2_ABS, 7=LINEAR_ABS)
		}
	}
	const float fogDensity = XboxRenderStates.GetXboxRenderStateAsFloat(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGDENSITY);
	const float fogStart = XboxRenderStates.GetXboxRenderStateAsFloat(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGSTART);
	const float fogEnd = XboxRenderStates.GetXboxRenderStateAsFloat(xbox::_X_D3DRENDERSTATETYPE::X_D3DRS_FOGEND);
	float fogStuff[4] = {fogTableMode, fogDensity, fogStart, fogEnd};
	CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_FOGINFO, fogStuff, 1);
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
	if (g_NV2A && !g_bInPullerContext) {
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
	if (g_NV2A && !g_bInPullerContext) {
		qemu_mutex_lock(&g_NV2A->GetDeviceState()->pgraph.pgraph_lock);
		pgraph_locked = true;
	}

	// Derive the vertex shader mode entirely from PGRAPH state.
	// g_Xbox_VertexShaderMode was previously set by the HLE SetVertexShader
	// patch on the game thread, which runs AHEAD of the puller — a race.
	// Now we read CSV0_D MODE for Program vs Fixed, and detect Passthrough
	// (XYZRHW pre-transformed vertices) via the VPSCL/VPOFF sign: the Xbox
	// D3D runtime maps screen coords to clip space such that the derived
	// X,Y origin is negative (e.g. -320,-240 for 640x480).  Normal fixed-
	// function viewports always produce X,Y >= 0.
	if (g_NV2A) {
		PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;
		uint32_t pgraph_mode = GET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_MODE);
		if (pgraph_mode == NV097_SET_TRANSFORM_EXECUTION_MODE_MODE_PROGRAM) {
			g_Xbox_VertexShaderMode = VertexShaderMode::ShaderProgram;
		} else {
			// MODE_FIXED: distinguish true fixed-function from passthrough
			// by checking the PGRAPH viewport constants.
			float vpoff0, vpoff1, vpscl0, vpscl1;
			std::memcpy(&vpoff0, &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPOFF][0], sizeof(float));
			std::memcpy(&vpoff1, &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPOFF][1], sizeof(float));
			std::memcpy(&vpscl0, &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL][0], sizeof(float));
			std::memcpy(&vpscl1, &pg->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL][1], sizeof(float));
			float xboxX = vpoff0 - vpscl0;
			float xboxY = vpoff1 + vpscl1;
			if (xboxX < 0.0f || xboxY < 0.0f) {
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
	if (g_NV2A) {
		CxbxD3D11UpdateRenderTargetFromPGRAPH(&g_NV2A->GetDeviceState()->pgraph);
	}

	// Set viewport from PGRAPH registers (authoritative).
	if (g_NV2A) {
		CxbxD3D11UpdateViewportFromPGRAPH(&g_NV2A->GetDeviceState()->pgraph);
	}

	// NOTE: Order is important here
   	// Some Texture States depend on RenderState values (Point Sprites)
   	// And some Pixel Shaders depend on Texture State values (BumpEnvMat, etc)
	CxbxUpdateHostTextures();
	CxbxUpdateHostTextureScaling();
   	XboxRenderStates.Apply();
   	XboxTextureStates.Apply();

	// Override blend/depth-stencil/rasterizer state from PGRAPH registers.
	if (g_NV2A) {
		auto pg = &g_NV2A->GetDeviceState()->pgraph;
		if (pg->surface_color.offset != 0) {
			CxbxD3D11UpdatePipelineStateFromPGRAPH(pg);
		}
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

// This function mimicks NV2A software callback events.
// Normally, these would be handled by actual push-buffer
// command handling at the point where they where inserted.
// Since our HLE mostly circumvents the NV2A pushbuffer,
// this function has to be called after 'pushing' functions.
void CxbxHandleXboxCallbacks()
{
	// The following can only work when host GPU queries are available
	if (g_pHostQueryCallbackEvent != nullptr) {
		// Query whether host GPU encountered a callback event already
		BOOL queryData = FALSE;
		if (S_FALSE == CxbxQueryGetData(g_pHostQueryCallbackEvent, &queryData, sizeof(queryData), 0)) {
			// If not, don't handle callbacks
			return;
		}
	}

	// Process inserted callbacks
	while (!g_Xbox_CallbackQueue.empty()) {
		// Fetch a callback from the FIFO callback queue
		s_Xbox_Callback XboxCallback = g_Xbox_CallbackQueue.front();
		g_Xbox_CallbackQueue.pop();

		// Differentiate between write and read callbacks
		if (XboxCallback.Type == xbox::X_D3DCALLBACK_WRITE) {
			// Write callbacks should wait until GPU is idle
			if (!CxbxFlushHostGPU()) {
				// Host GPU can't be flushed. In the old behaviour, we made the callback anyway
				// TODO : Should we keep doing that?
			}
		} else {
			assert(XboxCallback.Type == xbox::X_D3DCALLBACK_READ);
			// Should we mimick Read callback old behaviour?
			if (g_bHack_DisableHostGPUQueries) {
				// Note : Previously, we only processed Write, and ignored Read callbacks
				continue;
			} else {
				// New behaviour does place Read callbacks too
			}
		}

		// Make the callback
		XboxCallback.pCallback(XboxCallback.Context);
	}
}

// On Xbox, this function inserts push-buffer commands that
// will trigger the software handler to perform the callback
// when the GPU processes these commands.
// The type X_D3DCALLBACK_WRITE callbacks are prefixed with an
// wait-for-idle command, but otherwise they're identical.
// (Software handlers are triggered on NV2A via NV097_NO_OPERATION) 
void CxbxImpl_InsertCallback
(
	xbox::X_D3DCALLBACKTYPE	Type,
	xbox::X_D3DCALLBACK		pCallback,
	xbox::dword_xt				Context
)
{
	if (Type > xbox::X_D3DCALLBACK_WRITE) {
		LOG_TEST_CASE("Illegal callback type!");
		return;
	}

	if (pCallback == xbox::zeroptr) {
		LOG_TEST_CASE("pCallback == xbox::zeroptr!");
		return;
	}

	// Should we mimick old behaviour?
	if (g_bHack_DisableHostGPUQueries) {
		// Mimick old behaviour, in which only the final callback event
		// was remembered, by emptying the callback queue entirely :
		while (!g_Xbox_CallbackQueue.empty()) {
			g_Xbox_CallbackQueue.pop();
		}
	}

	// Push this callback's arguments into the callback queue :
	s_Xbox_Callback XboxCallback = { pCallback, Type, Context };
	g_Xbox_CallbackQueue.push(XboxCallback); // g_Xbox_CallbackQueue.emplace(pCallback, Type, Context); doesn't compile?

	// Does host supports GPU queries?
	if (g_pHostQueryCallbackEvent != nullptr) {
		// Insert a callback event on host GPU,
		// which will be handled by CxbxHandleXboxCallback
		CxbxQueryIssueEnd(g_pHostQueryCallbackEvent);
	}
}

xbox::void_xt CxbxImpl_SetPixelShader(xbox::dword_xt Handle)
{
   	// Cache the active shader handle
   	g_pXbox_PixelShader = (xbox::X_PixelShader*)Handle;

   	// Copy the Pixel Shader data to our RenderState handler (this includes values for pixel shader constants)
   	// This mirrors the fact that unpatched SetPixelShader does the same thing!
   	// This shouldn't be necessary anymore, but shaders still break if we don't do this
	// This breakage might be caused by our push-buffer processing could be "trailing behind" what our patches do;
	// By writing to render state during this patch, we avoid missing out on updates that push buffer commands would perform.
	// However, any updates that occur mid-way can overwrite what we store here, and still cause problems!
	// The only viable solution for that would be to draw entirely based on push-buffer handling (which might require removing possibly all D3D patches!)
   	if (g_pXbox_PixelShader != nullptr) {
   	   	// TODO : If D3DDevice_SetPixelShader() in XDKs don't overwrite the X_D3DRS_PS_RESERVED slot with PSDef.PSTextureModes,
   	   	// store it here and restore after memcpy, or alternatively, perform two separate memcpy's (the halves before, and after the reserved slot).
   	   	memcpy(XboxRenderStates.GetPixelShaderRenderStatePointer(), g_pXbox_PixelShader->pPSDef, sizeof(xbox::X_D3DPIXELSHADERDEF) - 3 * sizeof(DWORD));
   	   	// Copy the PSDef.PSTextureModes field to its dedicated slot, which lies outside the range of PixelShader render state slots.
   	   	// Always write from the new PSDef — a subsequent SetRenderState can override this if needed.
   	   	XboxRenderStates.SetXboxRenderState(xbox::X_D3DRS_PSTEXTUREMODES, g_pXbox_PixelShader->pPSDef->PSTextureModes);
		// NOTE: PGRAPH combiner registers are bridged at draw time in
		// CxbxD3D11UploadRCInterpreterState() to also catch subsequent
		// SetPixelShaderConstant / SetRenderState changes.
   	} else {
		// When clearing the pixel shader (handle=0), sync PSTextureModes from
		// the PSDef area that the native trampoline just wrote (default combiner
		// program). Without this, X_D3DRS_PSTEXTUREMODES retains the previous
		// shader's value — e.g. CUBEMAP mode leaks into subsequent fixed-function
		// draws, causing them to sample from the wrong SRV slot and produce black.
		const xbox::X_D3DPIXELSHADERDEF *pRSPSDef = (const xbox::X_D3DPIXELSHADERDEF*)(XboxRenderStates.GetPixelShaderRenderStatePointer());
		if (pRSPSDef) {
			XboxRenderStates.SetXboxRenderState(xbox::X_D3DRS_PSTEXTUREMODES, pRSPSDef->PSTextureModes);
		}
	}
}

// ******************************************************************
// * patch: D3DDevice_SetPixelShader
// ******************************************************************
