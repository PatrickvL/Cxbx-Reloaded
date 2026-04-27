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


/* Unused :
static xbox::dword_xt                  *g_Xbox_D3DDevice; // TODO: This should be a D3DDevice structure
*/

static void DrawInitialBlackScreen
(
)
{
   	// initially, show a black screen
   	// Only clear depth buffer and stencil if present
   	//
   	// Avoids following DirectX Debug Runtime error report
   	//    [424] Direct3D8: (ERROR) :Invalid flag D3DCLEAR_ZBUFFER: no zbuffer is associated with device. Clear failed. 
   	//

	CxbxD3DClear(
		/*Count=*/0,
		/*pRects=*/nullptr,
		D3DCLEAR_TARGET | (g_bHasDepth ? D3DCLEAR_ZBUFFER : 0) | (g_bHasStencil ? D3DCLEAR_STENCIL : 0),
		/*Color=*/0xFF000000, // TODO : Use constant for this
		/*Z=*/g_bHasDepth ? 1.0f : 0.0f,
		/*Stencil=*/0);

	CxbxBeginScene();

	CxbxPresent();
}

void CreateDefaultDevice
(
   	const xbox::X_D3DPRESENT_PARAMETERS     *pPresentationParameters
)
{
   	LOG_INIT;

   	// only one device should be created at once
   	if (g_pD3DDevice != nullptr) {
   	   	EmuLog(LOG_LEVEL::DEBUG, "CreateDefaultDevice releasing old Device.");

		CxbxEndScene();

   	   	ClearAllResourceCaches();

   	   	// TODO: ensure all other resources are cleaned up too

   	   	// Final release of IDirect3DDevice9 must be called from the window message thread
   	   	// See https://docs.microsoft.com/en-us/windows/win32/direct3d9/multithreading-issues
   	   	RunOnWndMsgThread([] {
   	   	   	// We only need to call bundled device release once here.
   	   	   	g_renderbase->DeviceRelease();
   	   	});
   	}

   	// Apply render scale factor for high-resolution rendering
   	g_RenderUpscaleFactor = g_XBVideo.renderScaleFactor;

   	// Setup the HostPresentationParameters
   	SetupPresentationParameters(pPresentationParameters);

	// This flag adds support for surfaces with a different color channel 
	// ordering than the API default. It is required for compatibility with
	// Direct2D.
	UINT creationFlags = D3D11_CREATE_DEVICE_BGRA_SUPPORT; // See enum D3D11_CREATE_DEVICE_FLAG
#if defined(_DEBUG)
	// If the project is in a debug build, enable debugging via SDK Layers.
	creationFlags |= D3D11_CREATE_DEVICE_DEBUG;
#endif
	// only use feature level 10.0
	D3D_FEATURE_LEVEL featureLevels[] = {
		D3D_FEATURE_LEVEL_11_0, // Required for cs_5_0, typed UAV access, ByteAddressBuffer
		D3D_FEATURE_LEVEL_10_0,
	};

	// Create the Direct3D 11 API device object and a corresponding context.
	ComPtr<ID3D11Device> device;
	ComPtr<ID3D11DeviceContext> context;
	HRESULT hr = D3D11CreateDevice(
		g_EmuCDPD.Adapter,
		g_EmuCDPD.DeviceType,
		nullptr,
		creationFlags,
		featureLevels,
		ARRAYSIZE(featureLevels),
		D3D11_SDK_VERSION, // UWP apps must set this to D3D11_SDK_VERSION.
		&device, // Returns the Direct3D device created.
		nullptr, // pFeatureLevel
		&context // Returns the device immediate context.
	);
#if defined(_DEBUG)
	// If debug layer failed (SDK not installed), retry without it
	if (FAILED(hr) && (creationFlags & D3D11_CREATE_DEVICE_DEBUG)) {
		EmuLog(LOG_LEVEL::WARNING, "D3D11CreateDevice failed with debug layer (hr=0x%08X), retrying without", hr);
		creationFlags &= ~D3D11_CREATE_DEVICE_DEBUG;
		hr = D3D11CreateDevice(
			g_EmuCDPD.Adapter,
			g_EmuCDPD.DeviceType,
			nullptr,
			creationFlags,
			featureLevels,
			ARRAYSIZE(featureLevels),
			D3D11_SDK_VERSION,
			&device,
			nullptr,
			&context
		);
	}
#endif
	// If device creation failed with a specific adapter or non-hardware driver type,
	// fall back to default adapter with D3D_DRIVER_TYPE_HARDWARE (most compatible, works with DXVK/Proton)
	if (FAILED(hr) && (g_EmuCDPD.Adapter != nullptr || g_EmuCDPD.DeviceType != D3D_DRIVER_TYPE_HARDWARE)) {
		EmuLog(LOG_LEVEL::WARNING, "D3D11CreateDevice failed (hr=0x%08X), retrying with default adapter and HARDWARE driver type", hr);
		hr = D3D11CreateDevice(
			nullptr,
			D3D_DRIVER_TYPE_HARDWARE,
			nullptr,
			creationFlags,
			featureLevels,
			ARRAYSIZE(featureLevels),
			D3D11_SDK_VERSION,
			&device,
			nullptr,
			&context
		);
	}
	// Last resort: fall back to WARP software rasterizer (Windows only, not available under Wine/Proton)
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "D3D11CreateDevice failed (hr=0x%08X), falling back to WARP", hr);
		hr = D3D11CreateDevice(
			nullptr,
			D3D_DRIVER_TYPE_WARP,
			nullptr,
			creationFlags,
			featureLevels,
			ARRAYSIZE(featureLevels),
			D3D11_SDK_VERSION,
			&device,
			nullptr,
			&context
		);
	}
   	DEBUG_D3DRESULT(hr, "D3D11CreateDevice");
	if (FAILED(hr))
		CxbxrAbort("D3D11CreateDevice failed (hr=0x%08X)", hr);

	// Store pointers to the Direct3D 11 API device and immediate context.
	device->QueryInterface(__uuidof(ID3D11Device), reinterpret_cast<void**>(&g_pD3DDevice));
	context->QueryInterface(__uuidof(ID3D11DeviceContext), reinterpret_cast<void**>(&g_pD3DDeviceContext));

	// Create a swap chain using the HWND (Win32 window)
	// Get DXGI objects from device
	ComPtr<IDXGIDevice1> dxgiDevice;
	g_pD3DDevice->QueryInterface(__uuidof(IDXGIDevice1), reinterpret_cast<void**>(dxgiDevice.GetAddressOf()));

	ComPtr<IDXGIAdapter> dxgiAdapter;
	dxgiDevice->GetAdapter(dxgiAdapter.GetAddressOf());

	ComPtr<IDXGIFactory2> dxgiFactory;
	dxgiAdapter->GetParent(__uuidof(IDXGIFactory2), reinterpret_cast<void**>(dxgiFactory.GetAddressOf()));

	// Configure swap chain description for Win32 HWND
	DXGI_SWAP_CHAIN_DESC1 SwapChainDesc = {};
	SwapChainDesc.Width = g_EmuCDPD.HostPresentationParameters.BackBufferWidth;
	SwapChainDesc.Height = g_EmuCDPD.HostPresentationParameters.BackBufferHeight;
	SwapChainDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM; // Common back buffer format
	SwapChainDesc.Stereo = FALSE;
	SwapChainDesc.SampleDesc.Count = 1;
	SwapChainDesc.SampleDesc.Quality = 0;
	SwapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	SwapChainDesc.BufferCount = 2;
	SwapChainDesc.Scaling = DXGI_SCALING_STRETCH;
	SwapChainDesc.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
	SwapChainDesc.AlphaMode = DXGI_ALPHA_MODE_UNSPECIFIED;
	SwapChainDesc.Flags = 0;

	DXGI_SWAP_CHAIN_FULLSCREEN_DESC fullscreenDesc = {};
	fullscreenDesc.RefreshRate.Numerator = g_EmuCDPD.HostPresentationParameters.FullScreen_RefreshRateInHz;
	fullscreenDesc.RefreshRate.Denominator = 1;
	fullscreenDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED;
	fullscreenDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;
	fullscreenDesc.Windowed = g_EmuCDPD.HostPresentationParameters.Windowed;

	ComPtr<IDXGISwapChain1> swapChain1;
	hr = dxgiFactory->CreateSwapChainForHwnd(
		g_pD3DDevice,
		g_hEmuWindow,
		&SwapChainDesc,
		fullscreenDesc.Windowed ? nullptr : &fullscreenDesc,
		nullptr, // pRestrictToOutput
		&swapChain1
	);
	DEBUG_D3DRESULT(hr, "IDXGIFactory2::CreateSwapChainForHwnd");
	if (FAILED(hr))
		CxbxrAbort("IDXGIFactory2::CreateSwapChainForHwnd failed");

	swapChain1->QueryInterface(__uuidof(IDXGISwapChain), reinterpret_cast<void**>(&g_pSwapChain));

	// Prevent DXGI from interfering with ALT+ENTER fullscreen toggle
	dxgiFactory->MakeWindowAssociation(g_hEmuWindow, DXGI_MWA_NO_ALT_ENTER);

	dxgiDevice->SetMaximumFrameLatency(1);

	// Configure the back buffer as a render target
	ComPtr<ID3D11Texture2D> backBuffer;
	hr = g_pSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), reinterpret_cast<void**>(backBuffer.GetAddressOf()));
	DEBUG_D3DRESULT(hr, "IDXGISwapChain::GetBuffer");

	// Create a render target view on the back buffer.
	hr = g_pD3DDevice->CreateRenderTargetView(backBuffer.Get(), nullptr, &g_pD3DBackBufferView);
	DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateRenderTargetView");

	// Keep a reference to the back buffer texture for dimension queries
	g_pD3DBackBufferSurface = backBuffer.Get();
	g_pD3DBackBufferSurface->AddRef();

	D3D11_TEXTURE2D_DESC backBufferDesc = {};
	backBuffer->GetDesc(&backBufferDesc);

	// Create a depth/stencil buffer and view to match the back buffer dimensions
	D3D11_TEXTURE2D_DESC depthDesc = {};
	depthDesc.Width = backBufferDesc.Width;
	depthDesc.Height = backBufferDesc.Height;
	depthDesc.MipLevels = 1;
	depthDesc.ArraySize = 1;
	depthDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
	depthDesc.SampleDesc.Count = 1;
	depthDesc.SampleDesc.Quality = 0;
	depthDesc.Usage = D3D11_USAGE_DEFAULT;
	depthDesc.BindFlags = D3D11_BIND_DEPTH_STENCIL;
	depthDesc.CPUAccessFlags = 0;
	depthDesc.MiscFlags = 0;

	hr = g_pD3DDevice->CreateTexture2D(&depthDesc, nullptr, &g_pD3DDepthStencilBuffer);
	DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateTexture2D (depth stencil)");

	if (SUCCEEDED(hr)) {
		hr = g_pD3DDevice->CreateDepthStencilView(g_pD3DDepthStencilBuffer, nullptr, &g_pD3DDepthStencilView);
		DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateDepthStencilView");
	}

	// Bind render target and depth stencil views to the output merger stage
	g_pD3DCurrentRTV = g_pD3DBackBufferView;
	g_pD3DCurrentHostRenderTarget = g_pD3DBackBufferSurface;
	g_pD3DDeviceContext->OMSetRenderTargets(1, &g_pD3DBackBufferView, g_pD3DDepthStencilView);

	// Store back buffer description for later use
	g_HostBackBufferDesc = backBufferDesc;

	// Set up default viewport to match back buffer size
	D3D11_VIEWPORT viewport = {};
	viewport.TopLeftX = 0.0f;
	viewport.TopLeftY = 0.0f;
	viewport.Width = static_cast<float>(backBufferDesc.Width);
	viewport.Height = static_cast<float>(backBufferDesc.Height);
	viewport.MinDepth = 0.0f;
	viewport.MaxDepth = 1.0f;
	g_pD3DDeviceContext->RSSetViewports(1, &viewport);

	// Initialize default D3D11 rasterizer state desc
	g_D3D11RasterizerDesc.FillMode = D3D11_FILL_SOLID;
	g_D3D11RasterizerDesc.CullMode = D3D11_CULL_BACK;
	g_D3D11RasterizerDesc.FrontCounterClockwise = FALSE;
	g_D3D11RasterizerDesc.DepthBias = 0;
	g_D3D11RasterizerDesc.SlopeScaledDepthBias = 0.0f;
	g_D3D11RasterizerDesc.DepthBiasClamp = 0.0f;
	g_D3D11RasterizerDesc.DepthClipEnable = TRUE;
	g_D3D11RasterizerDesc.ScissorEnable = FALSE;
	g_D3D11RasterizerDesc.MultisampleEnable = FALSE;
	g_D3D11RasterizerDesc.AntialiasedLineEnable = FALSE;

	// Initialize default D3D11 depth stencil state desc (Z test enabled, write enabled)
	g_D3D11DepthStencilDesc.DepthEnable = TRUE;
	g_D3D11DepthStencilDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
	g_D3D11DepthStencilDesc.DepthFunc = D3D11_COMPARISON_LESS_EQUAL;
	g_D3D11DepthStencilDesc.StencilEnable = FALSE;
	g_D3D11DepthStencilDesc.StencilReadMask = D3D11_DEFAULT_STENCIL_READ_MASK;
	g_D3D11DepthStencilDesc.StencilWriteMask = D3D11_DEFAULT_STENCIL_WRITE_MASK;
	g_D3D11DepthStencilDesc.FrontFace.StencilFunc = D3D11_COMPARISON_ALWAYS;
	g_D3D11DepthStencilDesc.FrontFace.StencilPassOp = D3D11_STENCIL_OP_KEEP;
	g_D3D11DepthStencilDesc.FrontFace.StencilFailOp = D3D11_STENCIL_OP_KEEP;
	g_D3D11DepthStencilDesc.FrontFace.StencilDepthFailOp = D3D11_STENCIL_OP_KEEP;
	g_D3D11DepthStencilDesc.BackFace = g_D3D11DepthStencilDesc.FrontFace;

	// Initialize default blend state (no blending)
	g_D3D11BlendDesc.AlphaToCoverageEnable = FALSE;
	g_D3D11BlendDesc.IndependentBlendEnable = FALSE;
	g_D3D11BlendDesc.RenderTarget[0].BlendEnable = FALSE;
	g_D3D11BlendDesc.RenderTarget[0].SrcBlend = D3D11_BLEND_ONE;
	g_D3D11BlendDesc.RenderTarget[0].DestBlend = D3D11_BLEND_ZERO;
	g_D3D11BlendDesc.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
	g_D3D11BlendDesc.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE;
	g_D3D11BlendDesc.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ZERO;
	g_D3D11BlendDesc.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
	g_D3D11BlendDesc.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;

	// Create the vertex shader constant buffer for D3D11
	{
		HRESULT cbHr = CxbxD3D11CreateConstantBuffer(CXBX_D3D11_VS_CB_COUNT * sizeof(float) * 4, true, &g_pD3D11VSConstantBuffer);
		DEBUG_D3DRESULT(cbHr, "g_pD3DDevice->CreateBuffer (VS constant buffer)");
		if (SUCCEEDED(cbHr)) {
			g_pD3DDeviceContext->VSSetConstantBuffers(CXBX_D3D11_VS_CB_SLOT, 1, &g_pD3D11VSConstantBuffer);
		}
	}

	// Create the zero-stride vertex defaults buffer for NV2A sticky attribute emulation
	CxbxD3D11CreateVertexDefaultsBuffer();

	// Initialize TEXCOORDINDEX remapping to identity (each stage reads its own texcoord set).
	// This ensures c219 is valid before the first passthrough draw, even if
	// CxbxUpdateHostTextureScaling() hasn't run yet.
	{
		float defaultTexCoordIndices[4] = { 0.0f, 1.0f, 2.0f, 3.0f };
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_TEXCOORDINDEX, defaultTexCoordIndices, 1);
	}

   	// Which texture formats does this device support?
   	DetermineSupportedD3DFormats();

	D3D11_QUERY_DESC QueryDesc;
	QueryDesc.Query = D3D11_QUERY_EVENT;
	QueryDesc.MiscFlags = 0;
   	// Can host driver create event queries?
   	if (SUCCEEDED(g_pD3DDevice->CreateQuery(&QueryDesc, nullptr))) {
   	   	// Is host GPU query creation enabled?
   	   	if (!g_bHack_DisableHostGPUQueries) {
   	   	   	// Create a D3D event query to handle "wait-for-idle" with
   	   	   	hr = g_pD3DDevice->CreateQuery(&QueryDesc, &g_pHostQueryWaitForIdle);
   	   	   	DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateQuery (wait for idle)");

   	   	   	// Create a D3D event query to handle "callback events" with
   	   	   	hr = g_pD3DDevice->CreateQuery(&QueryDesc, &g_pHostQueryCallbackEvent);
   	   	   	DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateQuery (callback event)");
   	   	}
   	} else {
   	   	LOG_TEST_CASE("Can't CreateQuery(D3DQUERYTYPE_EVENT) on host!");
   	}

   	// Can host driver create occlusion queries?
	QueryDesc.Query = D3D11_QUERY_OCCLUSION;
   	g_bEnableHostQueryVisibilityTest = false;
   	if (SUCCEEDED(g_pD3DDevice->CreateQuery(&QueryDesc, nullptr))) {
   	   	// Is host GPU query creation enabled?
   	   	if (!g_bHack_DisableHostGPUQueries) {
   	   	   	g_bEnableHostQueryVisibilityTest = true;
   	   	} else {
   	   	   	LOG_TEST_CASE("Disabled D3DQUERYTYPE_OCCLUSION on host!");
   	   	}
   	} else {
   	   	LOG_TEST_CASE("Can't CreateQuery(D3DQUERYTYPE_OCCLUSION) on host!");
   	}

   	DrawInitialBlackScreen();

   	// Set up ImGui's render backend
	ImGui_ImplDX11_Init(g_pD3DDevice, g_pD3DDeviceContext);
	CxbxD3D11InitBlit();
   	g_renderbase->SetDeviceRelease([] {
   	   	ImGui_ImplDX11_Shutdown();
   	   	g_VertexShaderCache.Clear();
   	   	CxbxD3D11ReleaseBackendResources(); // Also resets g_pD3DCurrentRTV for cached entries
   	   	if (g_pD3DDepthStencilView) { g_pD3DDepthStencilView->Release(); g_pD3DDepthStencilView = nullptr; }
   	   	if (g_pD3DDepthStencilBuffer) { g_pD3DDepthStencilBuffer->Release(); g_pD3DDepthStencilBuffer = nullptr; }
   	   	// Reset g_pD3DCurrentRTV before releasing back buffer view (it may equal g_pD3DBackBufferView)
   	   	g_pD3DCurrentRTV = nullptr;
   	   	if (g_pD3DBackBufferView) { g_pD3DBackBufferView->Release(); g_pD3DBackBufferView = nullptr; }
   	   	if (g_pD3DBackBufferSurface) { g_pD3DBackBufferSurface->Release(); g_pD3DBackBufferSurface = nullptr; }
   	   	if (g_pSwapChain) { g_pSwapChain->Release(); g_pSwapChain = nullptr; }
   	   	if (g_pD3DDeviceContext) { g_pD3DDeviceContext->Release(); g_pD3DDeviceContext = nullptr; }
   	   	g_pD3DDevice->Release();
   	});
}


// check if a resource has been registered yet (if not, register it)
bool GetHostRenderTargetDimensions(DWORD *pHostWidth, DWORD *pHostHeight, ID3D11Texture2D* pHostRenderTarget)
{
	if (pHostRenderTarget == nullptr) {
		pHostRenderTarget = CxbxGetCurrentRenderTarget();
	}

	// The following can only work if we could retrieve a host render target
	if (!pHostRenderTarget) {
		return false;
	}

	// Get current host render target dimensions
	D3D11_TEXTURE2D_DESC HostRenderTarget_Desc;
	pHostRenderTarget->GetDesc(&HostRenderTarget_Desc);

	*pHostWidth = HostRenderTarget_Desc.Width;
	*pHostHeight = HostRenderTarget_Desc.Height;

	return true;
}

DWORD ScaleDWORD(DWORD Value, DWORD FromMax, DWORD ToMax)
{
	uint64_t tmp = Value;
	tmp *= ToMax;
	tmp /= FromMax;
	return (DWORD)tmp;
}

void ValidateRenderTargetDimensions(DWORD HostRenderTarget_Width, DWORD HostRenderTarget_Height, DWORD XboxRenderTarget_Width, DWORD XboxRenderTarget_Height)
{
   	// This operation is often used to change the display resolution without calling SetRenderTarget!
   	// This works by updating the underlying Width & Height of the Xbox surface, without reallocating the data
   	// Because of this, we need to validate that the associated host resource still matches the dimensions of the Xbox Render Target
   	// If not, we must force them to be re-created
   	// TEST CASE: Chihiro Factory Test Program
   	DWORD XboxRenderTarget_Width_Scaled = XboxRenderTarget_Width * g_RenderUpscaleFactor;
   	DWORD XboxRenderTarget_Height_Scaled = XboxRenderTarget_Height * g_RenderUpscaleFactor;
   	if (HostRenderTarget_Width != XboxRenderTarget_Width_Scaled || HostRenderTarget_Height != XboxRenderTarget_Height_Scaled) {
   	   	LOG_TEST_CASE("Existing RenderTarget width/height changed");

   	   	FreeHostResource(GetHostResourceKey(g_pXbox_RenderTarget)); CxbxSetRenderTarget(GetHostSurface(g_pXbox_RenderTarget, D3DUSAGE_RENDERTARGET));
		FreeHostResource(GetHostResourceKey(g_pXbox_DepthStencil));
		CxbxSetDepthStencilSurface(GetHostSurface(g_pXbox_DepthStencil, D3DUSAGE_DEPTHSTENCIL));
   	}
}

float GetZScaleForPixelContainer(xbox::X_D3DPixelContainer* pSurface)
{
   	// If no surface was present, fallback to 1
   	if (pSurface == xbox::zeroptr) {
   	   	return 1.0f;
   	}

   	auto format = GetXboxPixelContainerFormat(pSurface);
   	switch (format) {
   	   	case xbox::X_D3DFMT_D16:
   	   	case xbox::X_D3DFMT_LIN_D16:
   	   	   	return 65535.0f;

   	   	case xbox::X_D3DFMT_D24S8:
   	   	case xbox::X_D3DFMT_LIN_D24S8:
   	   	   	return 16777215.0f;

   	   	case xbox::X_D3DFMT_F16:
   	   	case xbox::X_D3DFMT_LIN_F16:
   	   	   	return 511.9375f;

   	   	case xbox::X_D3DFMT_F24S8:
   	   	case xbox::X_D3DFMT_LIN_F24S8:
   	   	   	// 24bit floating point is close to precision maximum, so a lower value is used
   	   	   	// We can't use a double here since the vertex shader is only at float precision
   	   	   	return 1.0e30f; 
   	}

   	// Default to 1 if unknown depth format
   	LOG_TEST_CASE("GetZScaleForSurface: Unknown Xbox Depth Format");
   	return 1.0f;
}

void CxbxUpdateHostViewPortOffsetAndScaleConstants()
{
	// Xbox outputs vertex positions in rendertarget pixel coordinate space, with non-normalized Z
	// e.g. 0 < x < 640 and 0 < y < 480
	// We want to scale it back to normalized device coordinates i.e. XY are (-1, +1) and Z is (0, 1)

	// The screenspace is a combination of the rendertarget
	// and various scale factors
	// Get the rendertarget width and height
	float xboxRenderTargetWidth;
	float xboxRenderTargetHeight;
	GetRenderTargetBaseDimensions(xboxRenderTargetWidth, xboxRenderTargetHeight);

	float screenScaleX, screenScaleY;
	float aaOffsetX, aaOffsetY;
	GetScreenScaleFactors(screenScaleX, screenScaleY);
	GetMultiSampleOffset(aaOffsetX, aaOffsetY);

	// No half-pixel offset needed (D3D11 has pixel-center-at-0.5 natively)
	// https://aras-p.info/blog/2016/04/08/solving-dx9-half-pixel-offset/

	float xboxScreenspaceWidth = xboxRenderTargetWidth * screenScaleX;
	float xboxScreenspaceHeight = xboxRenderTargetHeight * screenScaleY;

	// Passthrough should range 0 to 1, instead of 0 to zbuffer depth
	// Test case: DoA3 character select
	// Detect passthrough from PGRAPH VPSCL/VPOFF sign to avoid racing
	// g_Xbox_VertexShaderMode which the game thread writes.
	bool isPassthrough = false;
	{
		auto pg_z = &(g_NV2A->GetDeviceState()->pgraph);
		float vpoff0, vpoff1, vpscl0, vpscl1;
		std::memcpy(&vpoff0, &pg_z->vsh_constants[NV_IGRAPH_XF_XFCTX_VPOFF][0], sizeof(float));
		std::memcpy(&vpoff1, &pg_z->vsh_constants[NV_IGRAPH_XF_XFCTX_VPOFF][1], sizeof(float));
		std::memcpy(&vpscl0, &pg_z->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL][0], sizeof(float));
		std::memcpy(&vpscl1, &pg_z->vsh_constants[NV_IGRAPH_XF_XFCTX_VPSCL][1], sizeof(float));
		float xboxX = vpoff0 - vpscl0;
		float xboxY = vpoff1 + vpscl1;
		isPassthrough = (xboxX < 0.0f || xboxY < 0.0f);
	}
	float zOutputScale = isPassthrough ? 1 : g_ZScale;

	float screenspaceScale[4] = { xboxScreenspaceWidth / 2,  -xboxScreenspaceHeight / 2, zOutputScale, 1 };
	float screenspaceOffset[4] = { xboxScreenspaceWidth / 2 + aaOffsetX, xboxScreenspaceHeight / 2 + aaOffsetY, 0, 0 };
	CxbxSetVertexShaderConstantF(CXBX_D3DVS_SCREENSPACE_SCALE_BASE, screenspaceScale, CXBX_D3DVS_NORMALIZE_SCALE_SIZE);
	CxbxSetVertexShaderConstantF(CXBX_D3DVS_SCREENSPACE_OFFSET_BASE, screenspaceOffset, CXBX_D3DVS_NORMALIZE_OFFSET_SIZE);

	// Reserved constants c[-38] (slot 58) and c[-37] (slot 59) hold viewport
	// scale/offset for screen-space transformation in programmable VS programs.
	// In the NV2A-driven path, PGRAPH already has the correct values written
	// by NV097_SET_VIEWPORT_SCALE/OFFSET through the push buffer.  Overwriting
	// them with HLE-computed values from g_Xbox_Viewport causes mismatches:
	// the Xbox thread updates g_Xbox_Viewport ahead of pfifo processing, so
	// mid-frame viewport changes produce stale overwrite values for earlier
	// draws.  Additionally, the HLE computation may not exactly match the
	// Xbox D3D runtime's internal NV2A register values.
	//
	// With draws going through pfifo, the PGRAPH values ARE the source of
	// truth — they are written before the draw in the push buffer and
	// processed sequentially by the puller.  CxbxUpdateHostVertexShaderConstants
	// already uploads them from pg->vsh_constants[58/59].
	//
	// TODO: Re-enable this overwrite if HLE draw patches are restored.
	// Test Case: GTA III, Soldier of Fortune II (needed when HLE draws are active)

}

// ******************************************************************
// * patch: D3DDevice_SetViewport
// ******************************************************************
void UpdateFixedFunctionShaderLight(int d3dLightIndex, Light* pShaderLight, D3DXVECTOR4* pLightAmbient) {
	if (d3dLightIndex == -1) {
		pShaderLight->Type = 0; // Disable the light
		return;
	}

	auto d3dLight = &d3d8LightState.Lights[d3dLightIndex];
	auto viewTransform = (D3DXMATRIX)d3d8TransformState.Transforms[xbox::X_D3DTS_VIEW];

	// TODO remove D3DX usage
	// Pre-transform light position to viewspace
	D3DXVECTOR4 positionV;
	D3DXVec3Transform(&positionV, (D3DXVECTOR3*)&d3dLight->Position, &viewTransform);
	pShaderLight->PositionV = (D3DXVECTOR3)positionV;

	// Pre-transform light direction to viewspace and normalize
	D3DXVECTOR4 directionV;
	D3DXMATRIX viewTransform3x3;
	D3DXMatrixIdentity(&viewTransform3x3);
	for (int y = 0; y < 3; y++) {
		for (int x = 0; x < 3; x++) {
			viewTransform3x3.m[x][y] = viewTransform.m[x][y];
		}
	}

	D3DXVec3Transform(&directionV, (D3DXVECTOR3*)&d3dLight->Direction, &viewTransform3x3);
	D3DXVec3Normalize((D3DXVECTOR3*)&pShaderLight->DirectionVN, (D3DXVECTOR3*)&directionV);

	bool SpecularEnable = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_SPECULARENABLE) != FALSE;

	// Map D3D light to state struct
	pShaderLight->Type = (int)d3dLight->Type;
	pShaderLight->Diffuse = toVector(d3dLight->Diffuse);
	pShaderLight->Specular = SpecularEnable ? toVector(d3dLight->Specular) : toVector(0);
	pShaderLight->Range = d3dLight->Range;
	pShaderLight->Falloff = d3dLight->Falloff;
	pShaderLight->Attenuation.x = d3dLight->Attenuation0;
	pShaderLight->Attenuation.y = d3dLight->Attenuation1;
	pShaderLight->Attenuation.z = d3dLight->Attenuation2;

	pLightAmbient->x += d3dLight->Ambient.r;
	pLightAmbient->y += d3dLight->Ambient.g;
	pLightAmbient->z += d3dLight->Ambient.b;

	auto cosHalfPhi = cos(d3dLight->Phi / 2);
	pShaderLight->CosHalfPhi = cosHalfPhi;
	pShaderLight->SpotIntensityDivisor = cos(d3dLight->Theta / 2) - cos(d3dLight->Phi / 2);
}

void UpdateFixedFunctionVertexShaderState()
{
	extern xbox::X_VERTEXATTRIBUTEFORMAT* GetXboxVertexAttributeFormat(); // TMP glue
	using namespace xbox;

	// Vertex blending
	// Prepare vertex blending mode variables used in transforms, below
	auto VertexBlend = XboxRenderStates.GetXboxRenderState(X_D3DRS_VERTEXBLEND);
	// Xbox and host D3DVERTEXBLENDFLAGS :
	//     D3DVBF_DISABLE           = 0 : 1 matrix,   0 weights => final weight 1
	//     D3DVBF_1WEIGHTS          = 1 : 2 matrices, 1 weights => final weight calculated
	//     D3DVBF_2WEIGHTS          = 3 : 3 matrices, 2 weights => final weight calculated
	//     D3DVBF_3WEIGHTS          = 5 : 4 matrices, 3 weights => final weight calculated
	// Xbox X_D3DVERTEXBLENDFLAGS :
	//   X_D3DVBF_2WEIGHTS2MATRICES = 2 : 2 matrices, 2 weights
	//   X_D3DVBF_3WEIGHTS3MATRICES = 4 : 3 matrices, 3 weights
	//   X_D3DVBF_4WEIGHTS4MATRICES = 6 : 4 matrices, 4 weights
	//
	if (VertexBlend > xbox::X_D3DVBF_4WEIGHTS4MATRICES) LOG_TEST_CASE("X_D3DRS_VERTEXBLEND out of range");
	// Calculate the number of matrices, by adding the LSB to turn (0,1,3,5) and (0,2,4,6) into (0,2,4,6); Then divide by 2 to get (0,1,2,3), and add 1 to get 1, 2, 3 or 4 matrices :
	auto NrBlendMatrices = ((VertexBlend + (VertexBlend & 1)) / 2) + 1;
	// Looking at the above values, 0 or the LSB of VertexBlend signals that the final weight needs to be calculated from all previous weigths (deducting them all from an initial 1) :
	auto CalcLastBlendWeight = (VertexBlend == xbox::X_D3DVBF_DISABLE) || (VertexBlend & 1);
	// Copy the resulting values over to shader state :
	ffShaderState.Modes.VertexBlend_NrOfMatrices = NrBlendMatrices;
	ffShaderState.Modes.VertexBlend_CalcLastWeight = CalcLastBlendWeight;

	// Transforms
	// Read transform matrices from PGRAPH XFCTX constants.
	// The Xbox D3D runtime writes World*View to MMAT, Inverse(World*View) to IMMAT,
	// and the full composite (World*View*Projection*Viewport) to CMAT.
	// PMAT is NOT written by the Xbox D3D runtime.
	//
	// CMAT includes the NV2A viewport transform (baked in by the Xbox D3D runtime for FF mode).
	// We must strip the viewport to get a pure clip-space projection matrix for D3D11.
	//
	// NV2A uses column-vector convention (clipPos = M * v), while our HLSL FF shader
	// uses row-vector convention (result = mul(v, M) = v * M). For NV2A matrices stored
	// register-per-row, a direct copy (no C++ transpose) makes the HLSL column-major
	// interpretation give transpose(M_pgraph), and mul(v, transpose(M)) = M * v.
	{
		PGRAPHState* pg = &g_NV2A->GetDeviceState()->pgraph;

		// Helper: direct-copy a 4x4 matrix from vsh_constants[base..base+3] — NO transpose
		auto ReadXFCTXMatrix = [&](D3DXMATRIX* pDst, int base) {
			for (int row = 0; row < 4; row++) {
				std::memcpy(&pDst->m[row][0], &pg->vsh_constants[base + row][0], 16);
			}
		};

		// Read MMAT0 (ModelView) and CMAT (Composite with viewport)
		D3DXMATRIX mmat, cmat;
		ReadXFCTXMatrix(&mmat, NV_IGRAPH_XF_XFCTX_MMAT0);
		ReadXFCTXMatrix(&cmat, NV_IGRAPH_XF_XFCTX_CMAT0);

		// Derive Projection-with-viewport = CMAT * inverse(MMAT)
		D3DXMATRIX mmatInv, projVP;
		bool validMMAT = (D3DXMatrixInverse(&mmatInv, nullptr, &mmat) != nullptr);
		if (validMMAT) {
			D3DXMatrixMultiply(&projVP, &cmat, &mmatInv);
		} else {
			// Singular MMAT (e.g., first few frames before state is populated)
			D3DXMatrixIdentity(&projVP);
			D3DXMatrixIdentity(&mmat);
		}

		// Strip the NV2A viewport from the projection.
		// The NV2A viewport matrix (column-vector):
		//   VP = [[sx,  0,  0, ox],  with sx = ox = W/2
		//         [ 0, sy,  0, oy],       sy = -H/2, oy = H/2
		//         [ 0,  0, sz, oz],       sz = z-buffer max, oz = 0
		//         [ 0,  0,  0,  1]]
		// PureProj = VP^-1 * projVP, computed element-by-element:
		//   row 0: (projVP[0][j] - ox * projVP[3][j]) / sx
		//   row 1: (projVP[1][j] - oy * projVP[3][j]) / sy
		//   row 2: projVP[2][j] / sz  (since oz = 0)
		//   row 3: projVP[3][j]
		D3DXMATRIX pureProj;
		float vpWidth = 0, vpHeight = 0;

		if (validMMAT && projVP.m[3][2] != 0.0f) {
			// Extract viewport offsets from the projection.
			// For standard perspective: projVP[3] = [0, 0, 1, 0], so
			// projVP[0][2] = ox (viewport X offset) and projVP[1][2] = oy (viewport Y offset).
			float ox = projVP.m[0][2] / projVP.m[3][2]; // typically W/2
			float oy = projVP.m[1][2] / projVP.m[3][2]; // typically H/2
			float sx = ox;      // centered viewport: sx = ox
			float sy = -oy;     // Y-flip: sy = -oy

			vpWidth  = 2.0f * ox;
			vpHeight = 2.0f * oy;

			// Z-buffer depth scale from surface format
			float sz = 1.0f;
			switch (pg->surface_shape.zeta_format) {
				case NV097_SET_SURFACE_FORMAT_ZETA_Z16:   sz = 65535.0f;    break;
				case NV097_SET_SURFACE_FORMAT_ZETA_Z24S8: sz = 16777215.0f; break;
				default:                                  sz = 65535.0f;    break;
			}

			for (int j = 0; j < 4; j++) {
				pureProj.m[0][j] = (projVP.m[0][j] - ox * projVP.m[3][j]) / sx;
				pureProj.m[1][j] = (projVP.m[1][j] - oy * projVP.m[3][j]) / sy;
				pureProj.m[2][j] =  projVP.m[2][j] / sz;
				pureProj.m[3][j] =  projVP.m[3][j];
			}
		} else {
			D3DXMatrixIdentity(&pureProj);
		}

		// Upload Projection (direct copy, no C++ transpose)
		std::memcpy(&ffShaderState.Transforms.Projection, &pureProj, sizeof(pureProj));

		// Set D3D11 viewport for FF mode (since VPSCL/VPOFF are zero, the
		// normal viewport update skips FF mode — we must set it here).
		if (vpWidth > 0 && vpHeight > 0) {
			D3D11_VIEWPORT hostViewport;
			hostViewport.TopLeftX = 0;
			hostViewport.TopLeftY = 0;
			hostViewport.Width    = vpWidth  * g_RenderUpscaleFactor;
			hostViewport.Height   = vpHeight * g_RenderUpscaleFactor;
			hostViewport.MinDepth = 0.0f;
			hostViewport.MaxDepth = 1.0f;
			CxbxSetViewport(&hostViewport);
		}

		// View matrix: PGRAPH XFCTX doesn't store View separately (only combined ModelView).
		// Set View to identity — the WorldView matrices already include it.
		D3DXMatrixIdentity((D3DXMATRIX*)&ffShaderState.Transforms.View);

		// Texture transforms (T0MAT..T3MAT) — direct copy
		static const int TnMAT[] = {
			NV_IGRAPH_XF_XFCTX_T0MAT, NV_IGRAPH_XF_XFCTX_T1MAT,
			NV_IGRAPH_XF_XFCTX_T2MAT, NV_IGRAPH_XF_XFCTX_T3MAT
		};
		for (unsigned i = 0; i < 4; i++) {
			ReadXFCTXMatrix((D3DXMATRIX*)&ffShaderState.Transforms.Texture[i], TnMAT[i]);
		}

		// WorldView matrices (MMAT0..MMAT3) — already pre-combined World*View, direct copy
		static const int MMATn[] = {
			NV_IGRAPH_XF_XFCTX_MMAT0, NV_IGRAPH_XF_XFCTX_MMAT1,
			NV_IGRAPH_XF_XFCTX_MMAT2, NV_IGRAPH_XF_XFCTX_MMAT3
		};
		for (unsigned i = 0; i < (unsigned)ffShaderState.Modes.VertexBlend_NrOfMatrices; i++) {
			ReadXFCTXMatrix((D3DXMATRIX*)&ffShaderState.Transforms.WorldView[i], MMATn[i]);
		}

		// WorldView inverse transpose — for normal transformation in lighting.
		// Compute from WorldView: upload = (MMAT^-1)^T so HLSL sees MMAT^-1,
		// and mul(normal, MMAT^-1) gives the correct (M^-1)^T * n transform.
		for (unsigned i = 0; i < (unsigned)ffShaderState.Modes.VertexBlend_NrOfMatrices; i++) {
			D3DXMATRIX wv, wvInv, wvInvT;
			ReadXFCTXMatrix(&wv, MMATn[i]);
			if (D3DXMatrixInverse(&wvInv, nullptr, &wv)) {
				D3DXMatrixTranspose(&wvInvT, &wvInv);
			} else {
				D3DXMatrixIdentity(&wvInvT);
			}
			std::memcpy(&ffShaderState.Transforms.WorldViewInverseTranspose[i], &wvInvT, sizeof(wvInvT));
		}
	}

	// Lighting
	// Point sprites aren't lit - 'each point is always rendered with constant colors.'
	// https://docs.microsoft.com/en-us/windows/win32/direct3d9/point-sprites
	bool PointSpriteEnable = XboxRenderStates.GetXboxRenderState(X_D3DRS_POINTSPRITEENABLE);
	bool LightingEnable = XboxRenderStates.GetXboxRenderState(X_D3DRS_LIGHTING);
	ffShaderState.Modes.Lighting = LightingEnable && !PointSpriteEnable;
	ffShaderState.Modes.TwoSidedLighting = XboxRenderStates.GetXboxRenderState(X_D3DRS_TWOSIDEDLIGHTING) ? 1 : 0;
	ffShaderState.Modes.LocalViewer = XboxRenderStates.GetXboxRenderState(X_D3DRS_LOCALVIEWER) ? 1 : 0;

	// Material sources
	bool ColorVertex = XboxRenderStates.GetXboxRenderState(X_D3DRS_COLORVERTEX) != FALSE;
	ffShaderState.Modes.AmbientMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_AMBIENTMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.DiffuseMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_DIFFUSEMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.SpecularMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_SPECULARMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.EmissiveMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_EMISSIVEMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.BackAmbientMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_BACKAMBIENTMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.BackDiffuseMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_BACKDIFFUSEMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.BackSpecularMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_BACKSPECULARMATERIALSOURCE) : D3DMCS_MATERIAL;
	ffShaderState.Modes.BackEmissiveMaterialSource = ColorVertex ? XboxRenderStates.GetXboxRenderState(X_D3DRS_BACKEMISSIVEMATERIALSOURCE) : D3DMCS_MATERIAL;

	// Point sprites; Fetch required variables
	float pointSize = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_POINTSIZE);
	float pointSize_Min = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_POINTSIZE_MIN);
	float pointSize_Max = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_POINTSIZE_MAX);
	bool PointScaleEnable = XboxRenderStates.GetXboxRenderState(X_D3DRS_POINTSCALEENABLE);
	float pointScale_A = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_POINTSCALE_A);
	float pointScale_B = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_POINTSCALE_B);
	float pointScale_C = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_POINTSCALE_C);
	float renderTargetHeight = (float)GetPixelContainerHeight(g_pXbox_RenderTarget);
	// Make sure to disable point scaling when point sprites are not enabled
	PointScaleEnable &= PointSpriteEnable;
	// Set variables in shader state
	ffShaderState.PointSprite.PointSize = PointSpriteEnable ? pointSize : 1.0f;
	ffShaderState.PointSprite.PointSize_Min = PointSpriteEnable ? pointSize_Min : 1.0f;
	ffShaderState.PointSprite.PointSize_Max = PointSpriteEnable ? pointSize_Max : 1.0f;
	ffShaderState.PointSprite.PointScaleABC.x = PointScaleEnable ? pointScale_A : 1.0f;
	ffShaderState.PointSprite.PointScaleABC.y = PointScaleEnable ? pointScale_B : 0.0f;
	ffShaderState.PointSprite.PointScaleABC.z = PointScaleEnable ? pointScale_C : 0.0f;
	ffShaderState.PointSprite.XboxRenderTargetHeight = PointScaleEnable ? renderTargetHeight : 1.0f;
	ffShaderState.PointSprite.RenderUpscaleFactor = (float)g_RenderUpscaleFactor;

	// Fog
	// Determine how the fog depth is transformed into the fog factor
	auto fogEnable = XboxRenderStates.GetXboxRenderState(X_D3DRS_FOGENABLE);
	auto fogTableMode = XboxRenderStates.GetXboxRenderState(X_D3DRS_FOGTABLEMODE);
	ffShaderState.Fog.Enable = fogEnable ? 1 : 0;
	ffShaderState.Fog.TableMode = fogTableMode;

	// Determine how fog depth is calculated
	if (fogEnable && fogTableMode != D3DFOG_NONE) {
		D3DXMATRIX projMtx = ffShaderState.Transforms.Projection;

		if (XboxRenderStates.GetXboxRenderState(X_D3DRS_RANGEFOGENABLE)) {
			LOG_TEST_CASE("Using RANGE fog");
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_RANGE;
		}
		else if (projMtx._14 == 0 &&
			projMtx._24 == 0 &&
			projMtx._34 == 0 &&
			projMtx._44 == 1) {
			LOG_TEST_CASE("Using Z fog");
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_Z;
		}
		else {
			// Test case:
			// Fog sample
			// JSRF (non-compliant projection matrix)
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_W;
		}

		ffShaderState.Fog.Density = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_FOGDENSITY);
		ffShaderState.Fog.Start = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_FOGSTART);
		ffShaderState.Fog.End = XboxRenderStates.GetXboxRenderStateAsFloat(X_D3DRS_FOGEND);
	}
	else {
		ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_NONE;
	}

	// Texture state
	for (int i = 0; i < xbox::X_D3DTS_STAGECOUNT; i++) {
		auto transformFlags = XboxTextureStates.Get(i, X_D3DTSS_TEXTURETRANSFORMFLAGS);
		ffShaderState.TextureStates[i].TextureTransformFlagsCount = transformFlags & ~D3DTTFF_PROJECTED;
		ffShaderState.TextureStates[i].TextureTransformFlagsProjected = transformFlags & D3DTTFF_PROJECTED;

		auto texCoordIndex = XboxTextureStates.Get(i, X_D3DTSS_TEXCOORDINDEX);
		ffShaderState.TextureStates[i].TexCoordIndex = texCoordIndex & 0x7; // 8 coords
		ffShaderState.TextureStates[i].TexCoordIndexGen = texCoordIndex >> 16; // D3DTSS_TCI flags
	}

	// Read current TexCoord component counts
	xbox::X_VERTEXATTRIBUTEFORMAT* pXboxVertexAttributeFormat = GetXboxVertexAttributeFormat();
	// Note : There seem to be other ways to access this, but we can use only this one;
	// This, because CxbxGetVertexDeclaration() can't be used, since it doesn't track VertexAttributes
	// (plus, it contains the overhead of shader lookup).
	// Another, GetXboxVertexShader(), can't be used, because it doesn't honor vertex attribute overrides
	// like those that apply for active SetVertexShaderInput.
	// Also, the xbox::X_D3DVertexShader.Dimensionality[] field contains somewhat strange values.
	for (int i = 0; i < xbox::X_D3DTS_STAGECOUNT; i++) {
		auto vertexDataFormat = pXboxVertexAttributeFormat->Slots[xbox::X_D3DVSDE_TEXCOORD0 + i].Format;
		reinterpret_cast<float*>(&ffShaderState.TexCoordComponentCount)[i] = (float)GetXboxVertexDataComponentCount(vertexDataFormat);
	}

	// Update lights from PGRAPH registers.
	// The NV2A light enable mask is in CSV0_D (2 bits per light: 0=off, 1=infinite/directional, 2=local/point, 3=spot).
	// Light colors in ltctxb[] are pre-multiplied by material by the Xbox D3D runtime, so we set
	// material to white to let the shader's (material × light) give the correct pre-multiplied result.
	{
		PGRAPHState* pg = &g_NV2A->GetDeviceState()->pgraph;
		uint32_t lightMask = pg->regs[NV_PGRAPH_CSV0_D / 4] & NV_PGRAPH_CSV0_D_LIGHTS;

		auto LightAmbient = D3DXVECTOR4(0.f, 0.f, 0.f, 0.f);

		// Helper to reinterpret uint32_t bit pattern as float
		auto AsFloat = [](uint32_t u) -> float { float f; std::memcpy(&f, &u, 4); return f; };

		for (size_t i = 0; i < ffShaderState.Lights.size(); i++) {
			Light* pShaderLight = &ffShaderState.Lights[i];
			unsigned nv2aType = (lightMask >> (i * 2)) & 0x3;

			if (nv2aType == 0) {
				pShaderLight->Type = 0; // Disabled
				continue;
			}

			// Map NV2A light type to shader type:
			//   NV2A 1 (INFINITE) → shader 3 (DIRECTIONAL)
			//   NV2A 2 (LOCAL)    → shader 1 (POINT)
			//   NV2A 3 (SPOT)     → shader 2 (SPOT)
			static const int typeMap[] = { 0, 3, 1, 2 };
			pShaderLight->Type = typeMap[nv2aType];

			// Diffuse color from ltctxb (3 floats stored as uint32_t bit patterns)
			int base = NV_IGRAPH_XF_LTCTXB_L0_DIF + (int)i * 6;
			pShaderLight->Diffuse = D3DXVECTOR4(
				AsFloat(pg->ltctxb[base][0]),
				AsFloat(pg->ltctxb[base][1]),
				AsFloat(pg->ltctxb[base][2]),
				1.0f);

			// Specular color
			bool SpecularEnable = XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_SPECULARENABLE) != FALSE;
			base = NV_IGRAPH_XF_LTCTXB_L0_SPC + (int)i * 6;
			if (SpecularEnable) {
				pShaderLight->Specular = D3DXVECTOR4(
					AsFloat(pg->ltctxb[base][0]),
					AsFloat(pg->ltctxb[base][1]),
					AsFloat(pg->ltctxb[base][2]),
					1.0f);
			} else {
				pShaderLight->Specular = D3DXVECTOR4(0, 0, 0, 0);
			}

			// Accumulate per-light ambient
			base = NV_IGRAPH_XF_LTCTXB_L0_AMB + (int)i * 6;
			LightAmbient.x += AsFloat(pg->ltctxb[base][0]);
			LightAmbient.y += AsFloat(pg->ltctxb[base][1]);
			LightAmbient.z += AsFloat(pg->ltctxb[base][2]);

			// Direction (for directional lights — already in view-space, normalized)
			pShaderLight->DirectionVN = D3DXVECTOR3(
				pg->light_infinite_direction[i][0],
				pg->light_infinite_direction[i][1],
				pg->light_infinite_direction[i][2]);

			// Position (for point/spot lights — already in view-space)
			pShaderLight->PositionV = D3DXVECTOR3(
				pg->light_local_position[i][0],
				pg->light_local_position[i][1],
				pg->light_local_position[i][2]);

			// Attenuation
			pShaderLight->Attenuation = D3DXVECTOR3(
				pg->light_local_attenuation[i][0],
				pg->light_local_attenuation[i][1],
				pg->light_local_attenuation[i][2]);

			// Range (stored in ltc1)
			pShaderLight->Range = AsFloat(pg->ltc1[NV_IGRAPH_XF_LTC1_r0 + i][0]);

			// Spot parameters from ltctxa
			int spotBase = NV_IGRAPH_XF_LTCTXA_L0_K + (int)i * 2;
			pShaderLight->Falloff = AsFloat(pg->ltctxa[spotBase][2]); // falloff stored in K[2]
			pShaderLight->CosHalfPhi = AsFloat(pg->ltctxa[spotBase][0]);
			pShaderLight->SpotIntensityDivisor = AsFloat(pg->ltctxa[spotBase][1]);
		}

		// Scene ambient from PGRAPH ltctxa[FR_AMB] (3 floats)
		D3DXVECTOR4 SceneAmbient(
			AsFloat(pg->ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][0]),
			AsFloat(pg->ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][1]),
			AsFloat(pg->ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][2]),
			0.f);
		D3DXVECTOR4 BackSceneAmbient(
			AsFloat(pg->ltctxa[NV_IGRAPH_XF_LTCTXA_BR_AMB][0]),
			AsFloat(pg->ltctxa[NV_IGRAPH_XF_LTCTXA_BR_AMB][1]),
			AsFloat(pg->ltctxa[NV_IGRAPH_XF_LTCTXA_BR_AMB][2]),
			0.f);

		ffShaderState.TotalLightsAmbient.Front = (D3DXVECTOR3)(LightAmbient + SceneAmbient);
		ffShaderState.TotalLightsAmbient.Back = (D3DXVECTOR3)(LightAmbient + BackSceneAmbient);

		// Material: set to white since NV2A ltctxb values are pre-multiplied by material.
		// The shader computes (material * light), so white material preserves the pre-multiplied values.
		ffShaderState.Materials[0].Diffuse  = D3DXVECTOR4(1, 1, 1, 1);
		ffShaderState.Materials[0].Ambient  = D3DXVECTOR4(1, 1, 1, 1);
		ffShaderState.Materials[0].Specular = D3DXVECTOR4(1, 1, 1, 1);
		ffShaderState.Materials[0].Emissive = D3DXVECTOR4(0, 0, 0, 0);
		ffShaderState.Materials[0].Power    = 0.0f;
		ffShaderState.Materials[1] = ffShaderState.Materials[0]; // back material
	}

	// Misc flags
	ffShaderState.Modes.NormalizeNormals = XboxRenderStates.GetXboxRenderState(X_D3DRS_NORMALIZENORMALS) ? 1 : 0;

	// Write fixed function state to shader constants
	const int slotSize = 16;
	const int fixedFunctionStateSize = (sizeof(FixedFunctionVertexShaderState) + slotSize - 1) / slotSize;
	CxbxSetVertexShaderConstantF(0, (float*)&ffShaderState, fixedFunctionStateSize);
}

