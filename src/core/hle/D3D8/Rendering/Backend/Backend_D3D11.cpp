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

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h"
#include "Backend_D3D11_Profiler.h"

// ******************************************************************
// * D3D11 device globals — definitions
// ******************************************************************
IDXGISwapChain                     *g_pSwapChain   = nullptr;
ID3D11DeviceContext                *g_pD3DDeviceContext = nullptr;
bool                                g_bTearingSupported = false;

// D3D11 render target and depth stencil views for the back buffer
ID3D11RenderTargetView             *g_pD3DBackBufferView = nullptr;
ID3D11DepthStencilView             *g_pD3DDepthStencilView = nullptr;
// D3D11 currently bound render target view (backbuffer or offscreen)
ID3D11RenderTargetView             *g_pD3DCurrentRTV = nullptr;
// D3D11 depth/stencil buffer texture
ID3D11Texture2D                    *g_pD3DDepthStencilBuffer = nullptr;
// D3D11 back buffer texture (used as fallback for dimension queries)
ID3D11Texture2D                    *g_pD3DBackBufferSurface = nullptr;
// D3D11 current host render target surface (used for dimension queries)
ID3D11Texture2D                    *g_pD3DCurrentHostRenderTarget = nullptr;

// ******************************************************************
// * D3D11 state descriptors — definitions
// * Initialize with valid D3D11 defaults so that the first
// * CreateXxxState() call (triggered by the dirty flags below)
// * succeeds even if the game hasn't set every render state yet.
// ******************************************************************
D3D11_RASTERIZER_DESC    g_D3D11RasterizerDesc = {
	/* FillMode              */ D3D11_FILL_SOLID,
	/* CullMode              */ D3D11_CULL_BACK,
	/* FrontCounterClockwise */ FALSE,
	/* DepthBias             */ 0,
	/* DepthBiasClamp        */ 0.0f,
	/* SlopeScaledDepthBias  */ 0.0f,
	/* DepthClipEnable       */ TRUE,
	/* ScissorEnable         */ FALSE,
	/* MultisampleEnable     */ FALSE,
	/* AntialiasedLineEnable */ FALSE,
};
D3D11_DEPTH_STENCIL_DESC g_D3D11DepthStencilDesc = {
	/* DepthEnable           */ TRUE,
	/* DepthWriteMask        */ D3D11_DEPTH_WRITE_MASK_ALL,
	/* DepthFunc             */ D3D11_COMPARISON_LESS,
	/* StencilEnable         */ FALSE,
	/* StencilReadMask       */ D3D11_DEFAULT_STENCIL_READ_MASK,
	/* StencilWriteMask      */ D3D11_DEFAULT_STENCIL_WRITE_MASK,
	/* FrontFace             */ { D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_COMPARISON_ALWAYS },
	/* BackFace              */ { D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_STENCIL_OP_KEEP, D3D11_COMPARISON_ALWAYS },
};
D3D11_BLEND_DESC         g_D3D11BlendDesc = {
	/* AlphaToCoverageEnable  */ FALSE,
	/* IndependentBlendEnable */ FALSE,
	/* RenderTarget[0]        */ {
		/* BlendEnable           */ FALSE,
		/* SrcBlend              */ D3D11_BLEND_ONE,
		/* DestBlend             */ D3D11_BLEND_ZERO,
		/* BlendOp               */ D3D11_BLEND_OP_ADD,
		/* SrcBlendAlpha         */ D3D11_BLEND_ONE,
		/* DestBlendAlpha        */ D3D11_BLEND_ZERO,
		/* BlendOpAlpha          */ D3D11_BLEND_OP_ADD,
		/* RenderTargetWriteMask */ D3D11_COLOR_WRITE_ENABLE_ALL,
	},
};

// Dirty flags
bool  g_bD3D11RasterizerStateDirty = true;
bool  g_bD3D11DepthStencilStateDirty = true;
bool  g_bD3D11BlendStateDirty = true;

// Additional OMSet* parameters
UINT  g_D3D11StencilRef = 0;
FLOAT g_D3D11BlendFactor[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
UINT  g_D3D11SampleMask = 0xFFFFFFFF;

// ******************************************************************
// * D3D11 state objects (internal — only used by ApplyDirtyStates)
// ******************************************************************
ComPtr<ID3D11RasterizerState>   g_pD3DRasterizerState;
ComPtr<ID3D11DepthStencilState> g_pD3DDepthStencilState;
ComPtr<ID3D11BlendState>        g_pD3DBlendState;

// ******************************************************************
// * Vertex shader constant buffer
// ******************************************************************
   	   ID3D11Buffer *g_pD3D11VSConstantBuffer = nullptr;
float         g_D3D11VSConstants[CXBX_D3D11_VS_CB_COUNT][4] = {};
UINT          g_D3D11VSConstantsDirtyMin = 0;
UINT          g_D3D11VSConstantsDirtyMax = CXBX_D3D11_VS_CB_COUNT; // full upload on first frame

// ******************************************************************
// * Blit shader resources (StretchRect replacement)
// ******************************************************************
ID3D11VertexShader  *g_pD3D11BlitVS = nullptr;
ID3D11PixelShader   *g_pD3D11BlitPS = nullptr;
ID3D11SamplerState  *g_pD3D11BlitSamplerLinear = nullptr;
ID3D11SamplerState  *g_pD3D11BlitSamplerPoint = nullptr;

// ******************************************************************
// * Point sprite geometry shader resources
// ******************************************************************
ID3D11GeometryShader *g_pD3D11PointSpriteGS = nullptr;
ID3D11Buffer         *g_pD3D11GSConstantBuffer = nullptr;
bool                  g_bPointSpriteEnabled = false;

// ******************************************************************
// * Thick line geometry shader resources
// ******************************************************************
ID3D11GeometryShader *g_pD3D11ThickLineGS = nullptr;
float                 g_fLineWidth = 1.0f;

// ******************************************************************
// * Compute shader unswizzle resources
// ******************************************************************
ID3D11ComputeShader  *g_pD3D11UnswizzleCS = nullptr;
ID3D11ComputeShader  *g_pD3D11UnswizzleBGRA_CS = nullptr; // float4 variant for B8G8R8A8_UNORM UAV
ID3D11Buffer         *g_pD3D11UnswizzleCB = nullptr; // constant buffer: maskX, maskY, width, height, bpp
ID3D11Buffer         *g_pD3D11UnswizzleStagingBuf = nullptr; // reusable ByteAddressBuffer for upload
UINT                  g_UnswizzleStagingBufSize = 0;
ID3D11ShaderResourceView *g_pD3D11UnswizzleSRV = nullptr; // SRV for staging buffer

// ******************************************************************
// * Compute shader index buffer conversion resources
// ******************************************************************
ID3D11ComputeShader       *g_pD3D11IndexConvertCS = nullptr;
ID3D11Buffer              *g_pD3D11IndexConvertCB = nullptr; // constant buffer: vertexCount, mode, isIndexed, pad
ID3D11Buffer              *g_pD3D11IndexConvertInputBuf = nullptr; // ByteAddressBuffer for source indices
UINT                       g_IndexConvertInputBufSize = 0;
ID3D11ShaderResourceView  *g_pD3D11IndexConvertInputSRV = nullptr;
ID3D11Buffer              *g_pD3D11IndexConvertOutputBuf = nullptr; // output buffer (INDEX_BUFFER + UAV)
UINT                       g_IndexConvertOutputBufSize = 0;
ID3D11UnorderedAccessView *g_pD3D11IndexConvertOutputUAV = nullptr;

// ******************************************************************
// * Compute shader palette texture expansion resources
// ******************************************************************
ID3D11ComputeShader       *g_pD3D11PaletteExpandCS = nullptr;
ID3D11Buffer              *g_pD3D11PaletteExpandCB = nullptr; // constant buffer: maskX, maskY, width, pad
ID3D11Buffer              *g_pD3D11PaletteBuf = nullptr; // 256-entry palette upload buffer
ID3D11ShaderResourceView  *g_pD3D11PaletteSRV = nullptr;

// ******************************************************************
// * Compute shader format conversion resources
// ******************************************************************
ID3D11ComputeShader       *g_pD3D11FormatConvertCS = nullptr;
ID3D11Buffer              *g_pD3D11FormatConvertCB = nullptr; // constant buffer: maskX, maskY, width, height, bpp, fmtType, swizzled, pad

// ******************************************************************
// * Register combiner interpreter (PS ubershader) resources
// ******************************************************************
ID3D11PixelShader         *g_pD3D11RCInterpreterPS = nullptr;
ID3D11Buffer              *g_pD3D11RCInterpreterAuxCB = nullptr; // PSAuxCBLayout (software-computed fields)

// ******************************************************************
// * Vertex shader interpreter (VS ubershader) resources
// ******************************************************************
bool                       g_bUseVSInterpreter = true; // default on — ubershader path
ID3D11VertexShader        *g_pD3D11VSInterpreterVS = nullptr;
ID3DBlob                  *g_pD3D11VSInterpreterBytecode = nullptr; // kept for input layout creation
ID3D11Buffer              *g_pD3D11XFPRBuf = nullptr;           // XFPR (Transform Program RAM) structured buffer — pg->xf.xfpr[]
ID3D11ShaderResourceView  *g_pD3D11XFPRSRV = nullptr;           // SRV for g_XFPR : register(t5)

// ******************************************************************
// * Compute shader vertex format conversion resources
// ******************************************************************
ID3D11ComputeShader       *g_pD3D11VertexConvertCS = nullptr;
ID3D11Buffer              *g_pD3D11VertexConvertCB = nullptr; // constant buffer: header + 16 element descriptors
ID3D11Buffer              *g_pD3D11VertexConvertSrcBuf = nullptr; // staging ByteAddressBuffer for source vertices
UINT                       g_VertexConvertSrcBufSize = 0;
ID3D11ShaderResourceView  *g_pD3D11VertexConvertSrcSRV = nullptr;

// ******************************************************************
// * D3D11 buffer/view creation helpers
// ******************************************************************

// Create a constant buffer. If bDynamic, uses DYNAMIC + CPU_ACCESS_WRITE
// (for Map/Unmap); otherwise DEFAULT (for UpdateSubresource).
HRESULT CxbxD3D11CreateConstantBuffer(UINT byteWidth, bool bDynamic, ID3D11Buffer** ppBuffer)
{
	D3D11_BUFFER_DESC desc = {};
	desc.ByteWidth = byteWidth;
	desc.Usage = bDynamic ? D3D11_USAGE_DYNAMIC : D3D11_USAGE_DEFAULT;
	desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	desc.CPUAccessFlags = bDynamic ? D3D11_CPU_ACCESS_WRITE : 0;
	return g_pD3DDevice->CreateBuffer(&desc, nullptr, ppBuffer);
}

// Create a ByteAddressBuffer (BUFFER_ALLOW_RAW_VIEWS) with optional CPU write access.
// If bDynamic, uses DYNAMIC + CPU_ACCESS_WRITE; otherwise DEFAULT.
HRESULT CxbxD3D11CreateRawBuffer(UINT byteWidth, bool bDynamic, UINT extraBindFlags, ID3D11Buffer** ppBuffer)
{
	D3D11_BUFFER_DESC desc = {};
	desc.ByteWidth = byteWidth;
	desc.Usage = bDynamic ? D3D11_USAGE_DYNAMIC : D3D11_USAGE_DEFAULT;
	desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | extraBindFlags;
	desc.CPUAccessFlags = bDynamic ? D3D11_CPU_ACCESS_WRITE : 0;
	desc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS;
	return g_pD3DDevice->CreateBuffer(&desc, nullptr, ppBuffer);
}

// Create a R32_TYPELESS raw SRV (for ByteAddressBuffer access in shaders).
HRESULT CxbxD3D11CreateRawBufferSRV(ID3D11Buffer* pBuffer, UINT byteWidth, ID3D11ShaderResourceView** ppSRV)
{
	D3D11_SHADER_RESOURCE_VIEW_DESC desc = {};
	desc.Format = DXGI_FORMAT_R32_TYPELESS;
	desc.ViewDimension = D3D11_SRV_DIMENSION_BUFFEREX;
	desc.BufferEx.Flags = D3D11_BUFFEREX_SRV_FLAG_RAW;
	desc.BufferEx.NumElements = byteWidth / 4;
	return g_pD3DDevice->CreateShaderResourceView(pBuffer, &desc, ppSRV);
}

// Create a R32_UINT typed UAV for a buffer.
HRESULT CxbxD3D11CreateTypedBufferUAV(ID3D11Buffer* pBuffer, UINT byteWidth, ID3D11UnorderedAccessView** ppUAV)
{
	D3D11_UNORDERED_ACCESS_VIEW_DESC desc = {};
	desc.Format = DXGI_FORMAT_R32_UINT;
	desc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
	desc.Buffer.FirstElement = 0;
	desc.Buffer.NumElements = byteWidth / 4;
	return g_pD3DDevice->CreateUnorderedAccessView(pBuffer, &desc, ppUAV);
}

// Ensure a dynamic raw staging buffer + SRV pair are at least requiredSize bytes.
// Rounds up to 4KB. Releases and re-creates if too small.
void CxbxD3D11EnsureRawStagingBuffer(
	UINT requiredSize,
	ID3D11Buffer** ppBuffer, UINT* pCurrentSize,
	ID3D11ShaderResourceView** ppSRV,
	const char* debugName)
{
	requiredSize = (requiredSize + 4095) & ~4095u;
	if (*ppBuffer && *pCurrentSize >= requiredSize)
		return;

	if (*ppSRV) { (*ppSRV)->Release(); *ppSRV = nullptr; }
	if (*ppBuffer) { (*ppBuffer)->Release(); *ppBuffer = nullptr; }

	HRESULT hr = CxbxD3D11CreateRawBuffer(requiredSize, /*bDynamic=*/true, 0, ppBuffer);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "%s: Failed to create buffer (%u bytes)", debugName, requiredSize);
		return;
	}
	*pCurrentSize = requiredSize;

	hr = CxbxD3D11CreateRawBufferSRV(*ppBuffer, requiredSize, ppSRV);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "%s: Failed to create SRV", debugName);
	}
}

// Bind CS pipeline state, dispatch, and unbind resources to avoid hazards.
void CxbxD3D11DispatchCS(
	ID3D11ComputeShader* pShader,
	ID3D11Buffer* pCB,
	UINT numSRVs,
	ID3D11ShaderResourceView* const* ppSRVs,
	ID3D11UnorderedAccessView* pUAV,
	UINT groupsX, UINT groupsY, UINT groupsZ)
{
	// Unbind all PS SRV slots to prevent SRV/UAV hazard on any resource
	// that may be simultaneously bound as a UAV for the CS dispatch.
	// 12 = 3 slot ranges × 4 stages: base (0-3), 3D (4-7), cube (8-11).
	static ID3D11ShaderResourceView* const nullPSSRVs[12] = {};
	g_pD3DDeviceContext->PSSetShaderResources(0, 12, nullPSSRVs);

	// Invalidate texture state cache so the next draw rebinds PS SRVs
	extern void CxbxInvalidateTextureStateCache();
	CxbxInvalidateTextureStateCache();

	g_pD3DDeviceContext->CSSetShader(pShader, nullptr, 0);
	g_pD3DDeviceContext->CSSetConstantBuffers(0, 1, &pCB);
	if (numSRVs > 0)
		g_pD3DDeviceContext->CSSetShaderResources(0, numSRVs, ppSRVs);
	g_pD3DDeviceContext->CSSetUnorderedAccessViews(0, 1, &pUAV, nullptr);
	g_pD3DDeviceContext->Dispatch(groupsX, groupsY, groupsZ);
	InterlockedIncrement(&g_ProfileCSDispatchCount);

	// Unbind CS resources to avoid hazards
	ID3D11ShaderResourceView* nullSRVs[2] = { nullptr, nullptr };
	ID3D11UnorderedAccessView* pNullUAV = nullptr;
	UINT unbindCount = numSRVs > 0 ? numSRVs : 1;
	g_pD3DDeviceContext->CSSetShaderResources(0, unbindCount, nullSRVs);
	g_pD3DDeviceContext->CSSetUnorderedAccessViews(0, 1, &pNullUAV, nullptr);
	g_pD3DDeviceContext->CSSetShader(nullptr, nullptr, 0);

#ifdef _DEBUG
	// Flush and check for device removal to catch TDRs at the exact dispatch
	g_pD3DDeviceContext->Flush();
	HRESULT hrRemoved = g_pD3DDevice->GetDeviceRemovedReason();
	if (FAILED(hrRemoved)) {
		EmuLog(LOG_LEVEL::ERROR2, "CxbxD3D11DispatchCS: Device removed after Dispatch! Reason: 0x%08X", hrRemoved);
	}
#endif
}

// Map a DYNAMIC buffer with WRITE_DISCARD, memcpy data, and Unmap.
HRESULT CxbxD3D11UpdateDynamicBuffer(ID3D11Buffer* pBuffer, const void* pData, size_t dataSize)
{
	D3D11_MAPPED_SUBRESOURCE mapped = {};
	HRESULT hr = g_pD3DDeviceContext->Map(pBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
	if (FAILED(hr))
		return hr;
	memcpy(mapped.pData, pData, dataSize);
	g_pD3DDeviceContext->Unmap(pBuffer, 0);
	return S_OK;
}

// ******************************************************************
// * Shader constant functions
// ******************************************************************
void CxbxSetVertexShaderConstantF(UINT startRegister, const float* pConstantData, UINT Vector4fCount)
{
	if (!g_pD3D11VSConstantBuffer || !pConstantData || Vector4fCount == 0)
		return;

	UINT endRegister = startRegister + Vector4fCount;
	if (endRegister > CXBX_D3D11_VS_CB_COUNT)
		endRegister = CXBX_D3D11_VS_CB_COUNT;

	UINT count = endRegister - startRegister;
	memcpy(g_D3D11VSConstants[startRegister], pConstantData, count * sizeof(float) * 4);

	// Expand dirty range to cover the written registers
	if (startRegister < g_D3D11VSConstantsDirtyMin)
		g_D3D11VSConstantsDirtyMin = startRegister;
	if (endRegister > g_D3D11VSConstantsDirtyMax)
		g_D3D11VSConstantsDirtyMax = endRegister;
}

void CxbxD3D11FlushVertexShaderConstants()
{
	if (!g_pD3D11VSConstantBuffer || g_D3D11VSConstantsDirtyMin >= g_D3D11VSConstantsDirtyMax)
		return;

	g_pD3DDeviceContext->UpdateSubresource(g_pD3D11VSConstantBuffer, 0, nullptr, g_D3D11VSConstants, 0, 0);

	g_D3D11VSConstantsDirtyMin = CXBX_D3D11_VS_CB_COUNT;
	g_D3D11VSConstantsDirtyMax = 0;
}

// ******************************************************************
// * Register combiner interpreter — init (loads precompiled CSO)
// ******************************************************************

bool CxbxD3D11InitRCInterpreter()
{
	if (g_pD3D11RCInterpreterPS)
		return true; // Already initialized

	// Load precompiled shader from embedded CSO blob
	ID3DBlob* pBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxRCInterpreterPS", &pBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "RC Interpreter: failed to load precompiled CSO");
		return false;
	}

	HRESULT hr = g_pD3DDevice->CreatePixelShader(pBlob->GetBufferPointer(), pBlob->GetBufferSize(),
		nullptr, &g_pD3D11RCInterpreterPS);
	pBlob->Release();
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "RC Interpreter CreatePixelShader failed: 0x%08X", hr);
		return false;
	}

	// Create the auxiliary constant buffer (software-computed fields only).
	// Kept as DEFAULT + UpdateSubresource: benchmarking showed Map/WRITE_DISCARD
	// was ~12% slower for this 112-byte buffer (544-616 fps vs 656-708 fps).
	// Small isolated CB updates are faster via UpdateSubresource (driver can
	// DMA-copy from the command buffer without allocation overhead).
	hr = CxbxD3D11CreateConstantBuffer(sizeof(PSAuxCBLayout), false, &g_pD3D11RCInterpreterAuxCB);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "RC Interpreter CreateConstantBuffer (aux) failed: 0x%08X", hr);
		g_pD3D11RCInterpreterPS->Release();
		g_pD3D11RCInterpreterPS = nullptr;
		return false;
	}

	EmuLog(LOG_LEVEL::INFO, "RC Interpreter ubershader loaded successfully (%u byte aux cbuffer)",
		(unsigned)sizeof(PSAuxCBLayout));
	return true;
}


// ******************************************************************
// * Vertex shader interpreter — init (loads precompiled CSO)
// ******************************************************************

bool CxbxD3D11InitVSInterpreter()
{
	if (g_pD3D11VSInterpreterVS)
		return true; // Already initialized

	// Load precompiled shader from embedded CSO blob
	ID3DBlob* pBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxVSInterpreterVS", &pBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "VS Interpreter: failed to load precompiled CSO");
		return false;
	}

	HRESULT hr = g_pD3DDevice->CreateVertexShader(pBlob->GetBufferPointer(), pBlob->GetBufferSize(),
		nullptr, &g_pD3D11VSInterpreterVS);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "VS Interpreter CreateVertexShader failed: 0x%08X", hr);
		pBlob->Release();
		return false;
	}

	// Keep the bytecode alive for input layout creation
	g_pD3D11VSInterpreterBytecode = pBlob;

	// Create the XFPR StructuredBuffer<uint4> (136 slots × 16 bytes = 2176 bytes)
	// XFPR = NV2A Transform Program RAM (on-chip XF SRAM, 136 × 92-bit instructions in 128-bit containers)
	{
		D3D11_BUFFER_DESC desc = {};
		desc.ByteWidth = XFPR_LENGTH * 4 * sizeof(uint32_t); // 136 × 16 = 2176
		desc.Usage = D3D11_USAGE_DYNAMIC;
		desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
		desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
		desc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
		desc.StructureByteStride = 4 * sizeof(uint32_t); // 16 bytes per uint4
		hr = g_pD3DDevice->CreateBuffer(&desc, nullptr, &g_pD3D11XFPRBuf);
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "VS Interpreter CreateBuffer (XFPR) failed: 0x%08X", hr);
			g_pD3D11VSInterpreterVS->Release();
			g_pD3D11VSInterpreterVS = nullptr;
			g_pD3D11VSInterpreterBytecode->Release();
			g_pD3D11VSInterpreterBytecode = nullptr;
			return false;
		}

		D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
		srvDesc.Format = DXGI_FORMAT_UNKNOWN; // structured buffer
		srvDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
		srvDesc.Buffer.FirstElement = 0;
		srvDesc.Buffer.NumElements = XFPR_LENGTH; // 136
		hr = g_pD3DDevice->CreateShaderResourceView(g_pD3D11XFPRBuf, &srvDesc, &g_pD3D11XFPRSRV);
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "VS Interpreter CreateSRV (XFPR) failed: 0x%08X", hr);
			g_pD3D11XFPRBuf->Release();
			g_pD3D11XFPRBuf = nullptr;
			g_pD3D11VSInterpreterVS->Release();
			g_pD3D11VSInterpreterVS = nullptr;
			g_pD3D11VSInterpreterBytecode->Release();
			g_pD3D11VSInterpreterBytecode = nullptr;
			return false;
		}
	}

	EmuLog(LOG_LEVEL::INFO, "VS Interpreter ubershader loaded successfully (%u byte XFPR SRV, shared mirror SRV at t%u)",
		(unsigned)(XFPR_LENGTH * 4 * sizeof(uint32_t)),
		(unsigned)CXBX_D3D11_VS_PGREGS_SRV_SLOT);
	return true;
}

// ******************************************************************
// * Blit shader init
// ******************************************************************
void CxbxD3D11InitBlit()
{
	LOG_INIT;

	// Load precompiled blit vertex shader
	ID3DBlob* pVSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxBlitVS", &pVSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load blit VS CSO");
		return;
	}
	HRESULT hr = g_pD3DDevice->CreateVertexShader(pVSBlob->GetBufferPointer(), pVSBlob->GetBufferSize(), nullptr, &g_pD3D11BlitVS);
	pVSBlob->Release();
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create blit VS");
		return;
	}

	// Load precompiled blit pixel shader
	ID3DBlob* pPSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxBlitPS", &pPSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load blit PS CSO");
		return;
	}
	hr = g_pD3DDevice->CreatePixelShader(pPSBlob->GetBufferPointer(), pPSBlob->GetBufferSize(), nullptr, &g_pD3D11BlitPS);
	pPSBlob->Release();
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create blit PS");
		return;
	}

	// Create linear sampler
	D3D11_SAMPLER_DESC sd = {};
	sd.Filter = D3D11_FILTER_MIN_MAG_MIP_LINEAR;
	sd.AddressU = sd.AddressV = sd.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
	sd.MaxLOD = D3D11_FLOAT32_MAX;
	hr = g_pD3DDevice->CreateSamplerState(&sd, &g_pD3D11BlitSamplerLinear);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create linear sampler");
	}

	// Create point sampler
	sd.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT;
	hr = g_pD3DDevice->CreateSamplerState(&sd, &g_pD3D11BlitSamplerPoint);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create point sampler");
	}

	// ************************************************************
	// Point sprite geometry shader
	// ************************************************************
	ID3DBlob* pGSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxPointSpriteGS", &pGSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load point sprite GS CSO");
	} else {
		hr = g_pD3DDevice->CreateGeometryShader(pGSBlob->GetBufferPointer(), pGSBlob->GetBufferSize(), nullptr, &g_pD3D11PointSpriteGS);
		pGSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create point sprite GS");
		}
	}

	// ************************************************************
	// Thick line geometry shader
	// ************************************************************
	pGSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxThickLineGS", &pGSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load thick line GS CSO");
	} else {
		hr = g_pD3DDevice->CreateGeometryShader(pGSBlob->GetBufferPointer(), pGSBlob->GetBufferSize(), nullptr, &g_pD3D11ThickLineGS);
		pGSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create thick line GS");
		}
	}

	// Create GS constant buffer (1 float4: inverse viewport dimensions + line width)
	hr = CxbxD3D11CreateConstantBuffer(16, false, &g_pD3D11GSConstantBuffer);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create GS constant buffer");
	}

	// ************************************************************
	// Compute shader for texture unswizzle
	// ************************************************************
	ID3DBlob* pCSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxUnswizzleCS", &pCSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load unswizzle CS CSO");
	} else {
		hr = g_pD3DDevice->CreateComputeShader(pCSBlob->GetBufferPointer(), pCSBlob->GetBufferSize(), nullptr, &g_pD3D11UnswizzleCS);
		pCSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create unswizzle CS");
		}
	}

	// BGRA variant of the unswizzle CS (writes float4 to B8G8R8A8_UNORM UAV)
	pCSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxUnswizzleBGRA_CS", &pCSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load unswizzle BGRA CS CSO");
	} else {
		hr = g_pD3DDevice->CreateComputeShader(pCSBlob->GetBufferPointer(), pCSBlob->GetBufferSize(), nullptr, &g_pD3D11UnswizzleBGRA_CS);
		pCSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create unswizzle BGRA CS");
		}
	}

	// Create unswizzle constant buffer (5 uints: maskX, maskY, width, height, bpp)
	hr = CxbxD3D11CreateConstantBuffer(32, true, &g_pD3D11UnswizzleCB);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create unswizzle CB");
	}

	// ---------------------------------------------------------------
	// Index buffer conversion compute shader (cs_5_0)
	// ---------------------------------------------------------------
	ID3DBlob* pIdxCSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxIndexConvertCS", &pIdxCSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load index convert CS CSO");
	} else {
		hr = g_pD3DDevice->CreateComputeShader(pIdxCSBlob->GetBufferPointer(), pIdxCSBlob->GetBufferSize(), nullptr, &g_pD3D11IndexConvertCS);
		pIdxCSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create index convert CS");
		}
	}

	// Create index convert constant buffer (4 uints: vertexCount, mode, isIndexed, pad)
	hr = CxbxD3D11CreateConstantBuffer(16, false, &g_pD3D11IndexConvertCB);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create index convert CB");
	}

	// ---------------------------------------------------------------
	// Palette texture expansion compute shader (cs_5_0)
	// ---------------------------------------------------------------
	ID3DBlob* pPalCSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxPaletteExpandCS", &pPalCSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load palette expand CS CSO");
	} else {
		hr = g_pD3DDevice->CreateComputeShader(pPalCSBlob->GetBufferPointer(), pPalCSBlob->GetBufferSize(), nullptr, &g_pD3D11PaletteExpandCS);
		pPalCSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create palette expand CS");
		}
	}

	// Create palette expand constant buffer (4 uints: maskX, maskY, width, pad)
	hr = CxbxD3D11CreateConstantBuffer(16, true, &g_pD3D11PaletteExpandCB);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create palette expand CB");
	}

	// Create palette data buffer (256 entries * 4 bytes = 1024 bytes)
	hr = CxbxD3D11CreateRawBuffer(1024, true, 0, &g_pD3D11PaletteBuf);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create palette buffer");
	} else {
		hr = CxbxD3D11CreateRawBufferSRV(g_pD3D11PaletteBuf, 1024, &g_pD3D11PaletteSRV);
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create palette SRV");
		}
	}

	// ---------------------------------------------------------------
	// Texture format conversion compute shader (cs_5_0)
	// ---------------------------------------------------------------
	ID3DBlob* pFmtCSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxFormatConvertCS", &pFmtCSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load format convert CS CSO");
	} else {
		hr = g_pD3DDevice->CreateComputeShader(pFmtCSBlob->GetBufferPointer(), pFmtCSBlob->GetBufferSize(), nullptr, &g_pD3D11FormatConvertCS);
		pFmtCSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create format convert CS");
		}
	}

	// Create format convert constant buffer (8 uints)
	hr = CxbxD3D11CreateConstantBuffer(32, true, &g_pD3D11FormatConvertCB);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create format convert CB");
	}

	// ---------------------------------------------------------------
	// Vertex format conversion compute shader (cs_5_0)
	// ---------------------------------------------------------------
	ID3DBlob* pVtxCSBlob = nullptr;
	if (!LoadPrecompiledCSO("CxbxVertexConvertCS", &pVtxCSBlob)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to load vertex convert CS CSO");
	} else {
		hr = g_pD3DDevice->CreateComputeShader(pVtxCSBlob->GetBufferPointer(), pVtxCSBlob->GetBufferSize(), nullptr, &g_pD3D11VertexConvertCS);
		pVtxCSBlob->Release();
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create vertex convert CS");
		}
	}

	// Create vertex convert constant buffer (header 16 bytes + 16 elements * 16 bytes = 272 bytes)
	hr = CxbxD3D11CreateConstantBuffer(272, true, &g_pD3D11VertexConvertCB);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "CxbxD3D11InitBlit: Failed to create vertex convert CB");
	}

	// Initialize page-tracked 64 MiB mirror (must precede vertex fetch init)
	CxbxPageTrackerInit();

	// Initialize vertex fetch resources
	CxbxD3D11VertexFetchInit();
}

// ******************************************************************
// * Blit implementation
// ******************************************************************
HRESULT CxbxD3D11Blt(
	ID3D11Texture2D* pSrc, const RECT* pSrcRect,
	ID3D11Texture2D* pDst, const RECT* pDstRect,
	D3DTEXTUREFILTERTYPE Filter)
{
	D3D11_TEXTURE2D_DESC srcDesc, dstDesc;
	pSrc->GetDesc(&srcDesc);
	pDst->GetDesc(&dstDesc);

	// Determine source region
	UINT srcX = pSrcRect ? pSrcRect->left : 0;
	UINT srcY = pSrcRect ? pSrcRect->top : 0;
	UINT srcW = pSrcRect ? (pSrcRect->right - pSrcRect->left) : srcDesc.Width;
	UINT srcH = pSrcRect ? (pSrcRect->bottom - pSrcRect->top) : srcDesc.Height;

	// Determine dest region
	UINT dstX = pDstRect ? pDstRect->left : 0;
	UINT dstY = pDstRect ? pDstRect->top : 0;
	UINT dstW = pDstRect ? (pDstRect->right - pDstRect->left) : dstDesc.Width;
	UINT dstH = pDstRect ? (pDstRect->bottom - pDstRect->top) : dstDesc.Height;

	// Fast path: same-size, same-format copy
	if (srcW == dstW && srcH == dstH && srcDesc.Format == dstDesc.Format) {
		D3D11_BOX srcBox = { srcX, srcY, 0, srcX + srcW, srcY + srcH, 1 };
		g_pD3DDeviceContext->CopySubresourceRegion(pDst, 0, dstX, dstY, 0, pSrc, 0, &srcBox);
		return S_OK;
	}

	// Scaled blit path
	if (!g_pD3D11BlitVS || !g_pD3D11BlitPS) {
		g_pD3DDeviceContext->CopyResource(pDst, pSrc);
		return S_OK;
	}

	// Create temporary SRV for source
	ID3D11ShaderResourceView* pSRV = nullptr;
	D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
	srvDesc.Format = (srcDesc.Format == DXGI_FORMAT_R8G8B8A8_TYPELESS) ? DXGI_FORMAT_R8G8B8A8_UNORM : srcDesc.Format;
	srvDesc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
	srvDesc.Texture2D.MipLevels = 1;
	HRESULT hr = g_pD3DDevice->CreateShaderResourceView(pSrc, &srvDesc, &pSRV);
	if (FAILED(hr)) return hr;

	// Create temporary RTV for destination
	ID3D11RenderTargetView* pRTV = nullptr;
	D3D11_RENDER_TARGET_VIEW_DESC rtvDesc = {};
	rtvDesc.Format = dstDesc.Format;
	rtvDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
	hr = g_pD3DDevice->CreateRenderTargetView(pDst, &rtvDesc, &pRTV);
	if (FAILED(hr)) {
		pSRV->Release();
		return hr;
	}

	// Save current pipeline state
	ID3D11RenderTargetView* pOldRTV = nullptr;
	ID3D11DepthStencilView* pOldDSV = nullptr;
	g_pD3DDeviceContext->OMGetRenderTargets(1, &pOldRTV, &pOldDSV);
	D3D11_VIEWPORT oldVP;
	UINT numVP = 1;
	g_pD3DDeviceContext->RSGetViewports(&numVP, &oldVP);
	// Save blend state so the game's blend settings don't leak into the blit
	ComPtr<ID3D11BlendState> pOldBlendState;
	FLOAT oldBlendFactor[4];
	UINT oldSampleMask;
	g_pD3DDeviceContext->OMGetBlendState(&pOldBlendState, oldBlendFactor, &oldSampleMask);
	// Save rasterizer state so the game's scissor rect doesn't clip the blit
	ComPtr<ID3D11RasterizerState> pOldRasterizerState;
	g_pD3DDeviceContext->RSGetState(&pOldRasterizerState);

	// Set blit pipeline state
	D3D11_VIEWPORT vp = { (FLOAT)dstX, (FLOAT)dstY, (FLOAT)dstW, (FLOAT)dstH, 0.0f, 1.0f };
	g_pD3DDeviceContext->RSSetViewports(1, &vp);
	g_pD3DDeviceContext->RSSetState(nullptr); // Default rasterizer: no scissor, no culling
	g_pD3DDeviceContext->OMSetRenderTargets(1, &pRTV, nullptr);
	g_pD3DDeviceContext->OMSetBlendState(nullptr, nullptr, 0xFFFFFFFF); // Default blend: no blending, all channels written
	g_pD3DDeviceContext->VSSetShader(g_pD3D11BlitVS, nullptr, 0);
	g_pD3DDeviceContext->PSSetShader(g_pD3D11BlitPS, nullptr, 0);
	g_pD3DDeviceContext->GSSetShader(nullptr, nullptr, 0); // Ensure point sprite GS doesn't interfere
	g_pD3DDeviceContext->PSSetShaderResources(0, 1, &pSRV);
	ID3D11SamplerState* pSampler = (Filter == D3DTEXF_LINEAR) ? g_pD3D11BlitSamplerLinear : g_pD3D11BlitSamplerPoint;
	g_pD3DDeviceContext->PSSetSamplers(0, 1, &pSampler);
	g_pD3DDeviceContext->IASetInputLayout(nullptr);
	g_pD3DDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	// Draw full-screen triangle
	g_pD3DDeviceContext->Draw(3, 0);

	// Restore previous state
	g_pD3DDeviceContext->OMSetRenderTargets(1, &pOldRTV, pOldDSV);
	g_pD3DDeviceContext->RSSetViewports(1, &oldVP);
	g_pD3DDeviceContext->RSSetState(pOldRasterizerState.Get());
	g_pD3DDeviceContext->OMSetBlendState(pOldBlendState.Get(), oldBlendFactor, oldSampleMask);
	if (pOldRTV) pOldRTV->Release();
	if (pOldDSV) pOldDSV->Release();

	// Unbind source SRV to avoid hazard
	ID3D11ShaderResourceView* pNullSRV = nullptr;
	g_pD3DDeviceContext->PSSetShaderResources(0, 1, &pNullSRV);

	// Invalidate texture state cache so the next draw rebinds PS SRVs
	extern void CxbxInvalidateTextureStateCache();
	CxbxInvalidateTextureStateCache();

	// Invalidate the PS state tracking since we bypassed CxbxSetPixelShader
	// to bind the blit PS directly. Without this, the next CxbxSetPixelShader
	// call would skip rebinding the game's pixel shader.
	CxbxInvalidateActivePixelShader();

	// Invalidate VS/GS/topology caches since we bound blit shaders directly
	extern void CxbxInvalidateVertexShaderCache();
	CxbxInvalidateVertexShaderCache();

	pRTV->Release();
	pSRV->Release();

	return S_OK;
}


// ******************************************************************
// * Release all backend resources (called from device release lambda)
// ******************************************************************
void CxbxD3D11ReleaseBackendResources()
{
	g_pD3DBlendState.Reset();
	g_pD3DDepthStencilState.Reset();
	g_pD3DRasterizerState.Reset();
	if (g_pD3D11VSConstantBuffer) { g_pD3D11VSConstantBuffer->Release(); g_pD3D11VSConstantBuffer = nullptr; }
	if (g_pD3D11BlitVS) { g_pD3D11BlitVS->Release(); g_pD3D11BlitVS = nullptr; }
	if (g_pD3D11BlitPS) { g_pD3D11BlitPS->Release(); g_pD3D11BlitPS = nullptr; }
	if (g_pD3D11BlitSamplerLinear) { g_pD3D11BlitSamplerLinear->Release(); g_pD3D11BlitSamplerLinear = nullptr; }
	if (g_pD3D11BlitSamplerPoint) { g_pD3D11BlitSamplerPoint->Release(); g_pD3D11BlitSamplerPoint = nullptr; }
	if (g_pD3D11PointSpriteGS) { g_pD3D11PointSpriteGS->Release(); g_pD3D11PointSpriteGS = nullptr; }
	if (g_pD3D11ThickLineGS) { g_pD3D11ThickLineGS->Release(); g_pD3D11ThickLineGS = nullptr; }
	if (g_pD3D11GSConstantBuffer) { g_pD3D11GSConstantBuffer->Release(); g_pD3D11GSConstantBuffer = nullptr; }
	if (g_pD3D11UnswizzleCS) { g_pD3D11UnswizzleCS->Release(); g_pD3D11UnswizzleCS = nullptr; }
	if (g_pD3D11UnswizzleBGRA_CS) { g_pD3D11UnswizzleBGRA_CS->Release(); g_pD3D11UnswizzleBGRA_CS = nullptr; }
	if (g_pD3D11UnswizzleCB) { g_pD3D11UnswizzleCB->Release(); g_pD3D11UnswizzleCB = nullptr; }
	if (g_pD3D11UnswizzleSRV) { g_pD3D11UnswizzleSRV->Release(); g_pD3D11UnswizzleSRV = nullptr; }
	if (g_pD3D11UnswizzleStagingBuf) { g_pD3D11UnswizzleStagingBuf->Release(); g_pD3D11UnswizzleStagingBuf = nullptr; }
	g_UnswizzleStagingBufSize = 0;
	if (g_pD3D11IndexConvertCS) { g_pD3D11IndexConvertCS->Release(); g_pD3D11IndexConvertCS = nullptr; }
	if (g_pD3D11IndexConvertCB) { g_pD3D11IndexConvertCB->Release(); g_pD3D11IndexConvertCB = nullptr; }
	if (g_pD3D11IndexConvertInputSRV) { g_pD3D11IndexConvertInputSRV->Release(); g_pD3D11IndexConvertInputSRV = nullptr; }
	if (g_pD3D11IndexConvertInputBuf) { g_pD3D11IndexConvertInputBuf->Release(); g_pD3D11IndexConvertInputBuf = nullptr; }
	g_IndexConvertInputBufSize = 0;
	if (g_pD3D11IndexConvertOutputUAV) { g_pD3D11IndexConvertOutputUAV->Release(); g_pD3D11IndexConvertOutputUAV = nullptr; }
	if (g_pD3D11IndexConvertOutputBuf) { g_pD3D11IndexConvertOutputBuf->Release(); g_pD3D11IndexConvertOutputBuf = nullptr; }
	g_IndexConvertOutputBufSize = 0;
	if (g_pD3D11FormatConvertCS) { g_pD3D11FormatConvertCS->Release(); g_pD3D11FormatConvertCS = nullptr; }
	if (g_pD3D11FormatConvertCB) { g_pD3D11FormatConvertCB->Release(); g_pD3D11FormatConvertCB = nullptr; }
	if (g_pD3D11RCInterpreterPS) { g_pD3D11RCInterpreterPS->Release(); g_pD3D11RCInterpreterPS = nullptr; }
	if (g_pD3D11RCInterpreterAuxCB) { g_pD3D11RCInterpreterAuxCB->Release(); g_pD3D11RCInterpreterAuxCB = nullptr; }
	if (g_pD3D11VSInterpreterVS) { g_pD3D11VSInterpreterVS->Release(); g_pD3D11VSInterpreterVS = nullptr; }
	if (g_pD3D11VSInterpreterBytecode) { g_pD3D11VSInterpreterBytecode->Release(); g_pD3D11VSInterpreterBytecode = nullptr; }
	if (g_pD3D11XFPRSRV) { g_pD3D11XFPRSRV->Release(); g_pD3D11XFPRSRV = nullptr; }
	if (g_pD3D11XFPRBuf) { g_pD3D11XFPRBuf->Release(); g_pD3D11XFPRBuf = nullptr; }
	ClearRTVCache();
	CxbxReleaseOverlayResources();
}

