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
#include <dxgi1_5.h> // IDXGIFactory5, DXGI_FEATURE_PRESENT_ALLOW_TEARING
#include "Backend/Backend_D3D11_PageTracker.h"
#include "devices\video\nv2a.h" // NV2AState


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

void CxbxInitHostD3DDevice()
{
	// Create the host D3D11 device before emulation starts.
	// This allows native Xbox Direct3D_CreateDevice to run unpatched,
	// which populates D3D_g_pDevice and all internal D3D device fields
	// (including the VBlank KEVENT that BlockUntilVerticalBlank waits on).
	if (g_pD3DDevice != nullptr) {
		return; // Already created
	}

	CreateDefaultDevice(nullptr);

	// Host-side init that was formerly in Direct3D_CreateDevice_End.
	CxbxResetPgraphSurfaceTracking();
}

// ---- Helper: Create D3D11 device with fallback chain ----
static void CreateD3D11DeviceWithFallbacks(UINT creationFlags)
{
   	LOG_INIT;

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
}

// ---- Helper: Create a DXGI swap chain using the HWND (Win32 window) for the emulator window ----
static void CreateSwapChainForWindow()
{
	LOG_INIT;
	ComPtr<IDXGIDevice1> dxgiDevice;
	g_pD3DDevice->QueryInterface(__uuidof(IDXGIDevice1), reinterpret_cast<void**>(dxgiDevice.GetAddressOf()));

	ComPtr<IDXGIAdapter> dxgiAdapter;
	dxgiDevice->GetAdapter(dxgiAdapter.GetAddressOf());

	ComPtr<IDXGIFactory2> dxgiFactory;
	dxgiAdapter->GetParent(__uuidof(IDXGIFactory2), reinterpret_cast<void**>(dxgiFactory.GetAddressOf()));

	// Get DXGI objects from device
	// Check if the system supports tearing (variable refresh rate / no-vsync fast path)
	bool bTearingSupported = false;
	{
		IDXGIFactory5* factory5 = nullptr;
		if (SUCCEEDED(dxgiFactory->QueryInterface(__uuidof(IDXGIFactory5), reinterpret_cast<void**>(&factory5)))) {
			BOOL allowTearing = FALSE;
			if (SUCCEEDED(factory5->CheckFeatureSupport(
					DXGI_FEATURE_PRESENT_ALLOW_TEARING, &allowTearing, sizeof(allowTearing)))) {
				bTearingSupported = (allowTearing == TRUE);
			}
			factory5->Release();
		}
	}

	// Configure swap chain description for Win32 HWND
	DXGI_SWAP_CHAIN_DESC1 SwapChainDesc = {};
	SwapChainDesc.Width = g_EmuCDPD.HostPresentationParameters.BackBufferWidth;
	SwapChainDesc.Height = g_EmuCDPD.HostPresentationParameters.BackBufferHeight;
	SwapChainDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM; // Common back buffer format
	SwapChainDesc.Stereo = FALSE;
	SwapChainDesc.SampleDesc.Count = 1;
	SwapChainDesc.SampleDesc.Quality = 0;
	SwapChainDesc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	SwapChainDesc.BufferCount = 3;  // Triple buffer to avoid Present blocking on buffer availability
	SwapChainDesc.Scaling = DXGI_SCALING_STRETCH;
	SwapChainDesc.SwapEffect = DXGI_SWAP_EFFECT_FLIP_DISCARD;
	SwapChainDesc.AlphaMode = DXGI_ALPHA_MODE_UNSPECIFIED;
	SwapChainDesc.Flags = bTearingSupported ? DXGI_SWAP_CHAIN_FLAG_ALLOW_TEARING : 0;

	DXGI_SWAP_CHAIN_FULLSCREEN_DESC fullscreenDesc = {};
	fullscreenDesc.RefreshRate.Numerator = g_EmuCDPD.HostPresentationParameters.FullScreen_RefreshRateInHz;
	fullscreenDesc.RefreshRate.Denominator = 1;
	fullscreenDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED;
	fullscreenDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;
	fullscreenDesc.Windowed = g_EmuCDPD.HostPresentationParameters.Windowed;

	ComPtr<IDXGISwapChain1> swapChain1;
	HRESULT hr = dxgiFactory->CreateSwapChainForHwnd(
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
	g_bTearingSupported = bTearingSupported;
	if (bTearingSupported) {
		printf("[CXBX] DXGI: Tearing (ALLOW_TEARING) is supported and enabled\n");
	}

	// Prevent DXGI from interfering with ALT+ENTER fullscreen toggle
	dxgiFactory->MakeWindowAssociation(g_hEmuWindow, DXGI_MWA_NO_ALT_ENTER);

	// Allow up to 2 frames in the present queue so the CPU can prepare
	// the next frame while the GPU composites the previous one.
	// (Default is 3, but 2 keeps input latency reasonable while avoiding
	// the Present-blocks-every-frame bottleneck of latency=1.)
	dxgiDevice->SetMaximumFrameLatency(2);
}

// ---- Helper: Set up back buffer, depth stencil, and default render target ----
static void SetupBackBufferAndDepthStencil()
{
	LOG_INIT;
	ComPtr<ID3D11Texture2D> backBuffer;
	HRESULT hr = g_pSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), reinterpret_cast<void**>(backBuffer.GetAddressOf()));
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
}

// ---- Helper: Initialize default D3D11 pipeline state descriptors ----
static void InitializeDefaultPipelineState()
{
	g_D3D11RasterizerDesc.FillMode = D3D11_FILL_SOLID;
	g_D3D11RasterizerDesc.CullMode = D3D11_CULL_BACK;
	g_D3D11RasterizerDesc.FrontCounterClockwise = FALSE;
	g_D3D11RasterizerDesc.DepthBias = 0;
	g_D3D11RasterizerDesc.SlopeScaledDepthBias = 0.0f;
	g_D3D11RasterizerDesc.DepthBiasClamp = 0.0f;
	g_D3D11RasterizerDesc.DepthClipEnable = FALSE; // NV2A has no depth clipping, only depth testing
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

	CreateD3D11DeviceWithFallbacks(creationFlags);
	CreateSwapChainForWindow();
	SetupBackBufferAndDepthStencil();
	InitializeDefaultPipelineState();

	// Create the vertex shader constant buffer for D3D11
	{
		HRESULT cbHr = CxbxD3D11CreateConstantBuffer(CXBX_D3D11_VS_CB_COUNT * sizeof(float) * 4, false, &g_pD3D11VSConstantBuffer);
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
   	   	   	HRESULT hr = g_pD3DDevice->CreateQuery(&QueryDesc, &g_pHostQueryWaitForIdle);
   	   	   	DEBUG_D3DRESULT(hr, "g_pD3DDevice->CreateQuery (wait for idle)");
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

	// Cache RT dimensions — GetDesc is a COM virtual call; skip when RT unchanged
	static ID3D11Texture2D* s_CachedRT = nullptr;
	static DWORD s_CachedWidth = 0, s_CachedHeight = 0;
	if (pHostRenderTarget != s_CachedRT) {
		D3D11_TEXTURE2D_DESC HostRenderTarget_Desc;
		pHostRenderTarget->GetDesc(&HostRenderTarget_Desc);
		s_CachedRT = pHostRenderTarget;
		s_CachedWidth = HostRenderTarget_Desc.Width;
		s_CachedHeight = HostRenderTarget_Desc.Height;
	}

	*pHostWidth = s_CachedWidth;
	*pHostHeight = s_CachedHeight;

	return true;
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

	// Z output scale: read directly from NV2A viewport scale Z component (VPSCL.z).
	// The game's VS multiplies clip-space Z by this value (via reserved constant c-38.z).
	// We divide by it to reverse back to normalized [0,1] for D3D11.
	// Clamp to host depth format range: since we use D24_UNORM (max representable = 1.0),
	// extremely large scales (e.g. 1e30 for F24S8) would lose precision in the dp4 math.
	// The host stores depth as UNORM regardless of Xbox float format, so the effective
	// max is limited by the host format's integer depth range.
	float zOutputScale = 1.0f;
	{
		auto pg_z = &(g_NV2A->GetDeviceState()->pgraph);
		float vpscl_z;
		std::memcpy(&vpscl_z, &pg_z->xf.xfctx[NV_IGRAPH_XF_XFCTX_VPSCL][2], sizeof(float));
		if (vpscl_z != 0.0f) {
			zOutputScale = vpscl_z;
		} else {
			// Fallback: use format-derived default if VPSCL not yet programmed
			auto surf = NV2AGetSurfaceState(pg_z);
			switch (surf.zetaFormat) {
				case NV097_SET_SURFACE_FORMAT_ZETA_Z16:   zOutputScale = 65535.0f;    break;
				case NV097_SET_SURFACE_FORMAT_ZETA_Z24S8: zOutputScale = 16777215.0f; break;
				default:                                  zOutputScale = 65535.0f;    break;
			}
		}
	}

	float screenspaceScale[4] = { xboxScreenspaceWidth / 2,  -xboxScreenspaceHeight / 2, zOutputScale, 1 };
	float screenspaceOffset[4] = { xboxScreenspaceWidth / 2 + aaOffsetX, xboxScreenspaceHeight / 2 + aaOffsetY, 0, 0 };

	// Skip VS constant upload if the computed values haven't changed
	static float s_LastScale[4] = { 0, 0, 0, 0 };
	static float s_LastOffset[4] = { 0, 0, 0, 0 };
	if (std::memcmp(screenspaceScale, s_LastScale, sizeof(s_LastScale)) != 0 ||
	    std::memcmp(screenspaceOffset, s_LastOffset, sizeof(s_LastOffset)) != 0) {
		std::memcpy(s_LastScale, screenspaceScale, sizeof(s_LastScale));
		std::memcpy(s_LastOffset, screenspaceOffset, sizeof(s_LastOffset));
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_SCREENSPACE_SCALE_BASE, screenspaceScale, CXBX_D3DVS_NORMALIZE_SCALE_SIZE);
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_SCREENSPACE_OFFSET_BASE, screenspaceOffset, CXBX_D3DVS_NORMALIZE_OFFSET_SIZE);
	}

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
	// already uploads them from pg->xf.xfctx[58/59].
	//
	// TODO: Re-enable this overwrite if HLE draw patches are restored.
	// Test Case: GTA III, Soldier of Fortune II (needed when HLE draws are active)

}

// ******************************************************************
// * Reconstruct specular power from NV2A specular params (6 floats).
// * Ported from xemu (pgraph.c: reconstruct_specular_power).
// * The Xbox D3D runtime converts a single material power (shininess)
// * into a piecewise polynomial stored in LTC1 specular params.
// * This function reverses that conversion to recover the original power.
// ******************************************************************
static float ReconstructSpecularPower(const float* params)
{
	// Values < 1.0 produce a positive c1 and can be recovered directly
	if (params[1] > 0.0f && params[2] < 1.0f) {
		return params[2] - (params[0] * 2.0f);
	}

	float c0 = params[0];
	float c3 = params[3];

	// Positive first coefficients shouldn't occur from DirectX-generated values
	if (c0 > 0.0f || c3 > 0.0f) {
		return 0.0001f;
	}

	// Lookup table: breakpoints in c0 space, each corresponding to a power range
	struct CurveCoeffs { float a, b, c; };
	static const float kStepPoints[] = {
		-0.022553957999f,  // power =    1.25
		-0.421539008617f,  // power =    4.00
		-0.678715527058f,  // power =    9.00
		-0.838916420937f,  // power =   20.00
		-0.961754500866f,  // power =   90.00
		-0.990773200989f,  // power =  375.00
		-0.994858562946f,  // power =  650.00
		-0.996561050415f,  // power = 1000.00
		-0.999547004700f,  // power = 1250.00
	};
	static const CurveCoeffs kCoeffs[] = {
		{ 1.000108475163f,     -9.838607076280f,      54.829089549713f  },
		{ 1.199164441703f,     -3.292603784852f,       7.799987995214f  },
		{ 8.653441252033f,     29.189473787191f,      43.586027561823f  },
		{-531.307758450301f,  117.398468683934f,     113.155490738338f  },
		{ -4.662713151292f,    1.221108944572f,        1.217360986939f  },
		{-124.435242105211f,  35.401219563514f,       35.408114377045f  },
		{ 10672560.259502f,   21565843.555824f,      10894794.336297f   },
		{-51973801.463934f, -104199997.554353f,     -52225454.356278f   },
		{ 972270.324080f,     2025882.096547f,        1054898.052467f   },
	};
	// Reconstruction methods: quadratic for most ranges, saturation growth rate for mid-ranges
	// Quadratic: a + b*x + c*x^2
	// Saturation: (a*x) / (b + c*x)
	enum { QUAD, SGR }; // reconstruction method
	static const int kMethod[] = { QUAD, QUAD, QUAD, SGR, SGR, SGR, QUAD, QUAD, QUAD };
	constexpr int kCount = sizeof(kStepPoints) / sizeof(kStepPoints[0]);

	auto Reconstruct = [](int method, float x, const CurveCoeffs& c) -> float {
		if (method == QUAD) return c.a + c.b * x + c.c * x * x;
		else                return (c.a * x) / (c.b + c.c * x);
	};

	float power = 0.f;
	for (int i = 0; i < kCount; ++i) {
		if (c0 > kStepPoints[i]) {
			power = Reconstruct(kMethod[i], c0, kCoeffs[i]);
			break;
		}
	}

	float halfPower = 0.f;
	for (int i = 0; i < kCount; ++i) {
		if (c3 > kStepPoints[i]) {
			halfPower = Reconstruct(kMethod[i], c3, kCoeffs[i]);
			break;
		}
	}

	// Extend range beyond 1250 using the half-power params
	if (power == 0.f || (halfPower > power && c0 < -0.1f)) {
		return halfPower * 2.f;
	}

	return power;
}

// ---- Helper: Update vertex blending state from PGRAPH ----
static uint32_t UpdateFFState_VertexBlending(PGRAPHState* pg, uint32_t csv0d)
{
	// Vertex blending — read from PGRAPH CSV0_D SKIN field
	// SKIN values 0..6 map directly to D3D VertexBlend (DISABLE, 1WEIGHTS, 2W2M, 2WEIGHTS, 3W3M, 3WEIGHTS, 4W4M)
	uint32_t skinMode = GET_MASK(csv0d, NV_PGRAPH_CSV0_D_SKIN);
	// NV2A SKIN field values:
	//   OFF  = 0 : 1 matrix,   0 weights => final weight 1
	//   2G   = 1 : 2 matrices, 1 weight,  generate last
	//   2    = 2 : 2 matrices, 2 weights
	//   3G   = 3 : 3 matrices, 2 weights, generate last
	//   3    = 4 : 3 matrices, 3 weights
	//   4G   = 5 : 4 matrices, 3 weights, generate last
	//   4    = 6 : 4 matrices, 4 weights
	//
	if (skinMode > NV_PGRAPH_CSV0_D_SKIN_4) LOG_TEST_CASE("PGRAPH CSV0_D SKIN out of range");
	// Calculate the number of matrices, by adding the LSB to turn (0,1,3,5) and (0,2,4,6) into (0,2,4,6); Then divide by 2 to get (0,1,2,3), and add 1 to get 1, 2, 3 or 4 matrices :
	auto NrBlendMatrices = ((skinMode + (skinMode & 1)) / 2) + 1;
	// Looking at the above values, 0 or the LSB of skinMode signals that the final weight needs to be calculated from all previous weigths (deducting them all from an initial 1) :
	auto CalcLastBlendWeight = (skinMode == NV_PGRAPH_CSV0_D_SKIN_OFF) || (skinMode & 1);
	// Copy the resulting values over to shader state :
	ffShaderState.Modes.VertexBlend_NrOfMatrices = NrBlendMatrices;
	ffShaderState.Modes.VertexBlend_CalcLastWeight = CalcLastBlendWeight;
	return skinMode;
}

// ---- Helper: Update transform matrices from PGRAPH XFCTX ----
static void UpdateFFState_Transforms(PGRAPHState* pg, uint32_t skinMode)
{
	using namespace xbox;

	// Transforms
	// Read transform matrices from PGRAPH XFCTX constants.
	// The Xbox D3D runtime writes World*View to MMAT, Inverse(World*View) to IMMAT,
	// and the full composite to CMAT.  PMAT is NOT written by the Xbox D3D runtime.
	//
	// CMAT contents depend on whether vertex blending (skinning) is active:
	//   Skinning OFF:  CMAT = VP * Proj * View * World  (everything baked in)
	//   Skinning ON:   CMAT = VP * Proj                 (MV is per-bone in MMAT0..3)
	// (This matches xemu's vsh-ff.c: "If skinning is off the composite matrix
	// already includes the MV matrix".)
	//
	// We need projVP = VP * Proj to strip the viewport and get pure projection.
	//   Skinning OFF:  projVP = CMAT * inv(MMAT)
	//   Skinning ON:   projVP = CMAT  (already VP * Proj)
	//
	// NV2A uses column-vector convention (clipPos = M * v), while our HLSL FF shader
	// uses row-vector convention (result = mul(v, M) = v * M). For NV2A matrices stored
	// register-per-row, a direct copy (no C++ transpose) makes the HLSL column-major
	// interpretation give transpose(M_pgraph), and mul(v, transpose(M)) = M * v.
	{
		// Helper: direct-copy a 4x4 matrix from vsh_constants[base..base+3] — NO transpose
		auto ReadXFCTXMatrix = [&](D3DXMATRIX* pDst, int base) {
			for (int row = 0; row < 4; row++) {
				std::memcpy(&pDst->m[row][0], &pg->xf.xfctx[base + row][0], 16);
			}
		};

		// Read MMAT0 (ModelView) and CMAT (Composite with viewport)
		D3DXMATRIX mmat, cmat;
		ReadXFCTXMatrix(&mmat, NV_IGRAPH_XF_XFCTX_MMAT0);
		ReadXFCTXMatrix(&cmat, NV_IGRAPH_XF_XFCTX_CMAT0);

		// Read viewport offset from PGRAPH XFCTX (half-pixel bias applied by Xbox runtime)
		float vpoff[4];
		for (int i = 0; i < 4; i++) {
			std::memcpy(&vpoff[i], &pg->xf.xfctx[NV_IGRAPH_XF_XFCTX_VPOFF][i], sizeof(float));
		}

		// NV2A FF pipeline: CMAT × position → screen-space coordinates (viewport baked in).
		// We upload CMAT directly and do the screen→NDC conversion in the HLSL shader,
		// matching xemu's approach:
		//   1. screenPos = CMAT × position (or CMAT × blended_eye_pos when skinning)
		//   2. screenPos.xy /= screenPos.w  (perspective divide for xy)
		//   3. screenPos.xy += VPOFF.xy     (half-pixel offset)
		//   4. clipPos.xy = (2 * screenPos.xy - surfaceSize) / surfaceSize * w  (screen→NDC)
		//   5. clipPos.z = screenPos.z / clipRange  (depth normalization)
		// D3D11 viewport is set to full surface size by CxbxD3D11UpdateViewportFromPGRAPH().
		std::memcpy(&ffShaderState.Transforms.Projection, &cmat, sizeof(cmat));
		ffShaderState.Modes.UseDirectComposite = (skinMode == NV_PGRAPH_CSV0_D_SKIN_OFF) ? 1 : 0;

		// Pass surface size and viewport offset to shader for screen→NDC conversion.
		// Guarantee >= 1 so the shader can unconditionally divide by these.
		auto surf = NV2AGetSurfaceState(pg);
		ffShaderState.Modes.SurfaceWidth = std::max(1.0f, static_cast<float>(surf.clipWidth));
		ffShaderState.Modes.SurfaceHeight = std::max(1.0f, static_cast<float>(surf.clipHeight));
		ffShaderState.Modes.ViewportOffsetX = vpoff[0];
		ffShaderState.Modes.ViewportOffsetY = vpoff[1];

		// Depth max (zmax) for Z normalization — read from NV2A viewport scale Z (VPSCL.z).
		// The viewport transform bakes VPSCL.z into CMAT's Z column, so we must reverse it.
		{
			float vpscl_z;
			std::memcpy(&vpscl_z, &pg->xf.xfctx[NV_IGRAPH_XF_XFCTX_VPSCL][2], sizeof(float));
			if (vpscl_z > 0.0f) {
				ffShaderState.Modes.DepthMax = vpscl_z;
			} else {
				// Fallback: format-based default if VPSCL not yet programmed
				switch (surf.zetaFormat) {
					case NV097_SET_SURFACE_FORMAT_ZETA_Z16:   ffShaderState.Modes.DepthMax = 65535.0f;    break;
					default:                                  ffShaderState.Modes.DepthMax = 16777215.0f; break;
				}
			}
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

		// Texgen plane matrices (TG0MAT..TG3MAT) — used for EYE_LINEAR/OBJECT_LINEAR texgen
		static const int TGnMAT[] = {
			NV_IGRAPH_XF_XFCTX_TG0MAT, NV_IGRAPH_XF_XFCTX_TG1MAT,
			NV_IGRAPH_XF_XFCTX_TG2MAT, NV_IGRAPH_XF_XFCTX_TG3MAT
		};
		for (unsigned i = 0; i < 4; i++) {
			ReadXFCTXMatrix((D3DXMATRIX*)&ffShaderState.Transforms.TexgenMatrix[i], TGnMAT[i]);
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
}

// ---- Helper: Update lighting and material source modes from PGRAPH ----
static void UpdateFFState_LightingModes(PGRAPHState* pg, uint32_t csv0c, bool PointSpriteEnable)
{
	bool LightingEnable = (csv0c & NV_PGRAPH_CSV0_C_LIGHTING) != 0;
	ffShaderState.Modes.Lighting = LightingEnable && !PointSpriteEnable;
	ffShaderState.Modes.TwoSidedLighting = (csv0c & NV_PGRAPH_CSV0_C_TWO_SIDE_LIGHTING) ? 1 : 0;
	ffShaderState.Modes.LocalViewer = (csv0c & NV_PGRAPH_CSV0_C_LOCALEYE) ? 1 : 0;

	// Material sources — read from PGRAPH CSV0_C (the Xbox D3D runtime bakes ColorVertex logic
	// into SET_COLOR_MATERIAL: when ColorVertex=FALSE, all sources are set to FROM_MATERIAL=0)
	// NV2A source values: 0=MATERIAL, 1=COLOR1, 2=COLOR2 — matches D3DMCS_* directly
	ffShaderState.Modes.AmbientMaterialSource  = GET_MASK(csv0c, NV_PGRAPH_CSV0_C_AMBIENT);
	ffShaderState.Modes.DiffuseMaterialSource  = GET_MASK(csv0c, NV_PGRAPH_CSV0_C_DIFFUSE);
	ffShaderState.Modes.SpecularMaterialSource = GET_MASK(csv0c, NV_PGRAPH_CSV0_C_SPECULAR);
	ffShaderState.Modes.EmissiveMaterialSource = GET_MASK(csv0c, NV_PGRAPH_CSV0_C_EMISSION);
	// NV2A has no separate back material source bits — back face uses same sources as front
	ffShaderState.Modes.BackAmbientMaterialSource  = ffShaderState.Modes.AmbientMaterialSource;
	ffShaderState.Modes.BackDiffuseMaterialSource  = ffShaderState.Modes.DiffuseMaterialSource;
	ffShaderState.Modes.BackSpecularMaterialSource = ffShaderState.Modes.SpecularMaterialSource;
	ffShaderState.Modes.BackEmissiveMaterialSource = ffShaderState.Modes.EmissiveMaterialSource;
}

// ---- Helper: Update point sprite state from PGRAPH ----
static void UpdateFFState_PointSprites(PGRAPHState* pg, uint32_t ctl3, bool PointSpriteEnable)
{
	// Point sprites — read from PGRAPH registers using NV2A's native formula.
	// NV_PGRAPH_POINTSIZE is a fixed-point integer (value / 8.0 = size in pixels).
	// NV2A point_params[0..7] are used directly (the Xbox D3D runtime pre-bakes
	// PointSize * RTHeight into params[3], attenuation into [0,1,2], etc.).
	// Formula (from NV2A hardware / xemu):
	//   Scaled:    oPts = clamp(rsqrt(A + B*d + C*d²) * p3 + p7, min, max)
	//   Non-scaled: oPts = max(1, POINTSIZE_REG / 8)
	// We unify by choosing constants so the same formula works for both paths.
	bool PointScaleEnable = (ctl3 & NV_PGRAPH_CONTROL_3_POINTPARAMSENABLE) != 0;
	PointScaleEnable &= PointSpriteEnable;

	if (PointScaleEnable) {
		// Scaled path: use NV2A params directly
		ffShaderState.PointSprite.PointScaleABC.x = pg->point_params[0]; // A
		ffShaderState.PointSprite.PointScaleABC.y = pg->point_params[1]; // B
		ffShaderState.PointSprite.PointScaleABC.z = pg->point_params[2]; // C
		float p3 = pg->point_params[3]; // scale (PointSize * RTHeight, set by Xbox D3D runtime)
		float p7 = pg->point_params[7]; // bias / min size
		float ptMin = std::min(p7, 63.875f);
		float ptMax = std::min(p3 + ptMin, 63.875f);
		ffShaderState.PointSprite.PointSize = p7;               // bias added after multiply
		ffShaderState.PointSprite.PointSize_Min = ptMin;        // clamp min
		ffShaderState.PointSprite.PointSize_Max = ptMax;        // clamp max
		ffShaderState.PointSprite.XboxRenderTargetHeight = p3;  // scale multiplier
	} else if (PointSpriteEnable) {
		// Non-scaled path: fixed size from POINTSIZE register (fixed-point / 8)
		float ptSize = (float)pg->regs[RI(NV_PGRAPH_POINTSIZE)] / 8.0f;
		ffShaderState.PointSprite.PointScaleABC.x = 1.0f; // rsqrt(1) = 1
		ffShaderState.PointSprite.PointScaleABC.y = 0.0f;
		ffShaderState.PointSprite.PointScaleABC.z = 0.0f;
		ffShaderState.PointSprite.PointSize = 0.0f;             // bias = 0
		ffShaderState.PointSprite.PointSize_Min = 1.0f;         // min 1 pixel
		ffShaderState.PointSprite.PointSize_Max = 63.875f;      // hardware max
		ffShaderState.PointSprite.XboxRenderTargetHeight = ptSize; // scale = the size
	} else {
		// Point sprites disabled: output 1.0 (doesn't matter, GS not bound)
		ffShaderState.PointSprite.PointScaleABC.x = 1.0f;
		ffShaderState.PointSprite.PointScaleABC.y = 0.0f;
		ffShaderState.PointSprite.PointScaleABC.z = 0.0f;
		ffShaderState.PointSprite.PointSize = 0.0f;
		ffShaderState.PointSprite.PointSize_Min = 1.0f;
		ffShaderState.PointSprite.PointSize_Max = 1.0f;
		ffShaderState.PointSprite.XboxRenderTargetHeight = 1.0f;
	}
	ffShaderState.PointSprite.RenderUpscaleFactor = (float)g_RenderUpscaleFactor;
}

// ---- Helper: Update fog state from PGRAPH ----
static void UpdateFFState_Fog(PGRAPHState* pg, uint32_t csv0d, uint32_t ctl3)
{
	bool fogEnable = (ctl3 & NV_PGRAPH_CONTROL_3_FOGENABLE) != 0;
	uint32_t fogMode = GET_MASK(ctl3, NV_PGRAPH_CONTROL_3_FOG_MODE);
	ffShaderState.Fog.Enable = fogEnable ? 1 : 0;
	ffShaderState.Fog.FogMode = fogMode;

	// Determine how fog depth is calculated from PGRAPH FOGGENMODE
	if (fogEnable) {
		uint32_t fogGenMode = GET_MASK(csv0d, NV_PGRAPH_CSV0_D_FOGGENMODE);

		// Map NV2A foggen to our depth mode enum
		switch (fogGenMode) {
		case NV_PGRAPH_CSV0_D_FOGGENMODE_SPEC_ALPHA:
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_NONE; // specular.a
			break;
		case NV_PGRAPH_CSV0_D_FOGGENMODE_RADIAL:
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_RANGE;
			break;
		case NV_PGRAPH_CSV0_D_FOGGENMODE_PLANAR:
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_W;
			break;
		case NV_PGRAPH_CSV0_D_FOGGENMODE_ABS_PLANAR:
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_W_ABS;
			break;
		case NV_PGRAPH_CSV0_D_FOGGENMODE_FOG_X:
		default:
			ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_NONE; // fog coord passthrough
			break;
		}

		float fogParam0; std::memcpy(&fogParam0, &pg->regs[RI(NV_PGRAPH_FOGPARAM0)], sizeof(float));
		float fogParam1; std::memcpy(&fogParam1, &pg->regs[RI(NV_PGRAPH_FOGPARAM1)], sizeof(float));
		ffShaderState.Fog.FogParam0 = fogParam0;
		ffShaderState.Fog.FogParam1 = fogParam1;
	}
	else {
		ffShaderState.Fog.DepthMode = FixedFunctionVertexShader::FOG_DEPTH_NONE;
		ffShaderState.Fog.FogParam0 = 0.0f;
		ffShaderState.Fog.FogParam1 = 0.0f;
	}
}

// ---- Helper: Update texture coordinate state from PGRAPH ----
static void UpdateFFState_TextureStates(PGRAPHState* pg)
{
	// Texture state — read from PGRAPH (authoritative, no HLE dependency)
	for (int i = 0; i < NV2A_MAX_TEXTURES; i++) {
		// TextureTransformFlags: derived from PGRAPH texture_matrix_enable[]
		// and the shader stage program mode (for PROJECTED).
		// When texture_matrix_enable is true, the NV2A always does a full 4x4
		// matrix multiply (equivalent to D3DTTFF_COUNT4). When false, count=0
		// (D3DTTFF_DISABLE — no texture matrix transform).
		bool matrixEnabled = pg->texture_matrix_enable[i];
		ffShaderState.TextureStates[i].TextureTransformFlagsCount = matrixEnabled ? 4 : 0;

		// PROJECTED comes from the shader stage program mode (2D/3D_PROJECTIVE).
		// NV_PGRAPH_SHADERPROG stores all 4 stages; each stage has 5 bits.
		static const uint32_t stageMasks[4] = {
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE0,
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE1,
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE2,
			NV097_SET_SHADER_STAGE_PROGRAM_STAGE3
		};
		uint32_t shaderProg = pg->regs[RI(NV_PGRAPH_SHADERPROG)];
		uint32_t stageMode = GET_MASK(shaderProg, stageMasks[i]);
		// 2D_PROJECTIVE=1 and 3D_PROJECTIVE=2 are the projective modes
		bool projected = (stageMode == 1 || stageMode == 2);
		ffShaderState.TextureStates[i].TextureTransformFlagsProjected = projected ? D3DTTFF_PROJECTED : 0;

		// TexCoordIndex: low bits = which texcoord set, high bits = texgen mode.
		// On NV2A, texcoord routing is identity for FF (stage i uses TEXCOORD i).
		// Texgen mode is stored in CSV1_A (stages 0,1) / CSV1_B (stages 2,3).
		unsigned int csvReg = (i < 2) ? NV_PGRAPH_CSV1_A : NV_PGRAPH_CSV1_B;
		unsigned int sMask  = (i % 2) ? NV_PGRAPH_CSV1_A_T1_S : NV_PGRAPH_CSV1_A_T0_S;
		uint32_t texgenS = GET_MASK(pg->regs[RI(csvReg)], sMask);
		// Map NV2A texgen values to D3DTSS_TCI values (high 16 bits of TEXCOORDINDEX)
		unsigned int tci = 0; // TCI_PASSTHRU
		switch (texgenS) {
		case NV_PGRAPH_CSV1_A_T0_S_DISABLE:        tci = 0; break; // TCI_PASSTHRU
		case NV_PGRAPH_CSV1_A_T0_S_EYE_LINEAR:     tci = 2; break; // TCI_CAMERASPACEPOSITION
		case NV_PGRAPH_CSV1_A_T0_S_OBJECT_LINEAR:  tci = 4; break; // TCI_OBJECT (Xbox ext)
		case NV_PGRAPH_CSV1_A_T0_S_SPHERE_MAP:     tci = 5; break; // TCI_SPHERE (Xbox ext)
		case NV_PGRAPH_CSV1_A_T0_S_NORMAL_MAP:     tci = 1; break; // TCI_CAMERASPACENORMAL
		case NV_PGRAPH_CSV1_A_T0_S_REFLECTION_MAP: tci = 3; break; // TCI_CAMERASPACEREFLECTIONVECTOR
		}
		ffShaderState.TextureStates[i].TexCoordIndex = i; // identity routing
		ffShaderState.TextureStates[i].TexCoordIndexGen = tci;
	}

	// Read current TexCoord component counts from PGRAPH vertex attributes.
	// PGRAPH is the authoritative source since the puller processes attribute
	// format commands from the push buffer before each draw.
	for (int i = 0; i < NV2A_MAX_TEXTURES; i++) {
		int attrIdx = NV2A_VERTEX_ATTR_TEXTURE0 + i;
		const VertexAttribute& attr = pg->vertex_attributes[attrIdx];
		float componentCount;
		if (attr.count > 0) {
			componentCount = (float)attr.count;
		} else {
			componentCount = 2.0f; // Safe default (2D texcoords)
		}
		reinterpret_cast<float*>(&ffShaderState.TexCoordComponentCount)[i] = componentCount;
	}
}

// ---- Helper: Update lights and materials from PGRAPH ----
// The NV2A light enable mask is in CSV0_D (2 bits per light: 0=off, 1=infinite/directional, 2=local/point, 3=spot).
// Light colors in ltctxb[] are pre-multiplied by material by the Xbox D3D runtime, so we set
// material to white to let the shader's (material × light) give the correct pre-multiplied result.
static void UpdateFFState_Lighting(PGRAPHState* pg, uint32_t csv0c)
{
	// Skip if no lighting-related state changed since last call.
	// NV2A_DIRTY_LIGHTING covers ltctxa/ltctxb/ltc1/light writes;
	// NV2A_DIRTY_PGRAPH covers CSV0_D (light enable mask) and CSV0_C (specular enable).
	static uint32_t s_lastLightingGen = 0;
	static uint32_t s_lastPgraphGen = 0;
	if (pg->dirty[NV2A_DIRTY_LIGHTING] == s_lastLightingGen && pg->dirty[NV2A_DIRTY_PGRAPH] == s_lastPgraphGen)
		return;
	s_lastLightingGen = pg->dirty[NV2A_DIRTY_LIGHTING];
	s_lastPgraphGen = pg->dirty[NV2A_DIRTY_PGRAPH];

	// Helper to reinterpret uint32_t bit pattern as float
	auto AsFloat = [](uint32_t u) -> float { float f; std::memcpy(&f, &u, 4); return f; };

	uint32_t lightMask = pg->regs[RI(NV_PGRAPH_CSV0_D)] & NV_PGRAPH_CSV0_D_LIGHTS;
	auto LightAmbient = D3DXVECTOR4(0.f, 0.f, 0.f, 0.f);

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
			AsFloat(pg->xf.ltctxb[base][0]),
			AsFloat(pg->xf.ltctxb[base][1]),
			AsFloat(pg->xf.ltctxb[base][2]),
			1.0f);

		// Specular color
		bool SpecularEnable = (csv0c & NV_PGRAPH_CSV0_C_SPECULAR_ENABLE) != 0;
		base = NV_IGRAPH_XF_LTCTXB_L0_SPC + (int)i * 6;
		if (SpecularEnable) {
			pShaderLight->Specular = D3DXVECTOR4(
				AsFloat(pg->xf.ltctxb[base][0]),
				AsFloat(pg->xf.ltctxb[base][1]),
				AsFloat(pg->xf.ltctxb[base][2]),
				1.0f);
		} else {
			pShaderLight->Specular = D3DXVECTOR4(0, 0, 0, 0);
		}

		// Accumulate per-light ambient
		base = NV_IGRAPH_XF_LTCTXB_L0_AMB + (int)i * 6;
		LightAmbient.x += AsFloat(pg->xf.ltctxb[base][0]);
		LightAmbient.y += AsFloat(pg->xf.ltctxb[base][1]);
		LightAmbient.z += AsFloat(pg->xf.ltctxb[base][2]);

		// Direction (for directional lights — already in view-space, normalized)
		pShaderLight->DirectionVN = D3DXVECTOR3(
			pg->light[i].infinite_direction[0],
			pg->light[i].infinite_direction[1],
			pg->light[i].infinite_direction[2]);

		// Position (for point/spot lights — already in view-space)
		pShaderLight->PositionV = D3DXVECTOR3(
			pg->light[i].local_position[0],
			pg->light[i].local_position[1],
			pg->light[i].local_position[2]);

		// Attenuation
		pShaderLight->Attenuation = D3DXVECTOR3(
			pg->light[i].local_attenuation[0],
			pg->light[i].local_attenuation[1],
			pg->light[i].local_attenuation[2]);

		// Range (stored in ltc1)
		pShaderLight->Range = AsFloat(pg->xf.ltc1[NV_IGRAPH_XF_LTC1_r0 + i][0]);

		// Spot parameters from ltctxa
		int spotBase = NV_IGRAPH_XF_LTCTXA_L0_K + (int)i * 2;
		pShaderLight->Falloff = AsFloat(pg->xf.ltctxa[spotBase][2]); // falloff stored in K[2]
		pShaderLight->CosHalfPhi = AsFloat(pg->xf.ltctxa[spotBase][0]);
		pShaderLight->SpotIntensityDivisor = AsFloat(pg->xf.ltctxa[spotBase][1]);
	}

	// Scene ambient from PGRAPH ltctxa[FR_AMB] (3 floats)
	D3DXVECTOR4 SceneAmbient(
		AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][0]),
		AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][1]),
		AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_FR_AMB][2]),
		0.f);
	D3DXVECTOR4 BackSceneAmbient(
		AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_BR_AMB][0]),
		AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_BR_AMB][1]),
		AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_BR_AMB][2]),
		0.f);

	ffShaderState.TotalLightsAmbient.Front = (D3DXVECTOR3)(LightAmbient + SceneAmbient);
	ffShaderState.TotalLightsAmbient.Back = (D3DXVECTOR3)(LightAmbient + BackSceneAmbient);

	// Material: set to white since NV2A ltctxb values are pre-multiplied by material.
	// The shader computes (material * light), so white material preserves the pre-multiplied values.
	// Emission is already baked into the scene ambient register (FR_AMB/BR_AMB) by the Xbox D3D runtime.
	// Material alpha comes from NV097_SET_MATERIAL_ALPHA → ltctxa[CM_COL][3].
	float materialAlpha     = AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_CM_COL][3]);
	float backMaterialAlpha = AsFloat(pg->xf.ltctxa[NV_IGRAPH_XF_LTCTXA_BCM_COL][3]);

	ffShaderState.Materials[0].Diffuse  = D3DXVECTOR4(1, 1, 1, materialAlpha);
	ffShaderState.Materials[0].Ambient  = D3DXVECTOR4(1, 1, 1, 1);
	ffShaderState.Materials[0].Specular = D3DXVECTOR4(1, 1, 1, 1);
	ffShaderState.Materials[0].Emissive = D3DXVECTOR4(0, 0, 0, 0);

	// Reconstruct specular power from NV2A's 6 polynomial coefficients (LTC1).
	// Front specular params: ltc1[l0][0..3] + ltc1[l0+1][0..1]
	float frontParams[6];
	for (int j = 0; j < 4; j++) frontParams[j]     = AsFloat(pg->xf.ltc1[NV_IGRAPH_XF_LTC1_l0][j]);
	for (int j = 0; j < 2; j++) frontParams[4 + j]  = AsFloat(pg->xf.ltc1[NV_IGRAPH_XF_LTC1_l0 + 1][j]);
	ffShaderState.Materials[0].Power = ReconstructSpecularPower(frontParams);

	ffShaderState.Materials[1] = ffShaderState.Materials[0]; // back material (start from front)
	ffShaderState.Materials[1].Diffuse.w = backMaterialAlpha;

	// Back specular params: ltc1[Bl0][0..3] + ltc1[Bl0+1][0..1]
	float backParams[6];
	for (int j = 0; j < 4; j++) backParams[j]     = AsFloat(pg->xf.ltc1[NV_IGRAPH_XF_LTC1_Bl0][j]);
	for (int j = 0; j < 2; j++) backParams[4 + j]  = AsFloat(pg->xf.ltc1[NV_IGRAPH_XF_LTC1_Bl0 + 1][j]);
	ffShaderState.Materials[1].Power = ReconstructSpecularPower(backParams);
}

// ******************************************************************
// * UpdateFixedFunctionVertexShaderState
// ******************************************************************
void UpdateFixedFunctionVertexShaderState()
{
	PGRAPHState* pg = &g_NV2A->GetDeviceState()->pgraph;
	uint32_t csv0c = pg->regs[RI(NV_PGRAPH_CSV0_C)];
	uint32_t csv0d = pg->regs[RI(NV_PGRAPH_CSV0_D)];
	uint32_t ctl3  = pg->regs[RI(NV_PGRAPH_CONTROL_3)];

	uint32_t skinMode = UpdateFFState_VertexBlending(pg, csv0d);
	UpdateFFState_Transforms(pg, skinMode);

	// Point sprite enable comes from NV_PGRAPH_SETUPRASTER (D3DRS_POINTSPRITEENABLE →
	// NV097_SET_POINT_SMOOTH_ENABLE → NV_PGRAPH_SETUPRASTER_POINTSMOOTHENABLE).
	// Point scale/params enable comes from NV_PGRAPH_CONTROL_3_POINTPARAMSENABLE
	// (D3DRS_POINTSCALEENABLE → NV097_SET_POINT_PARAMETERS_ENABLE).
	// These are distinct Xbox render states and must be derived from separate PGRAPH bits.
	uint32_t setupRaster = pg->regs[RI(NV_PGRAPH_SETUPRASTER)];
	bool PointSpriteEnable = (setupRaster & NV_PGRAPH_SETUPRASTER_POINTSMOOTHENABLE) != 0;

	UpdateFFState_LightingModes(pg, csv0c, PointSpriteEnable);
	UpdateFFState_PointSprites(pg, ctl3, PointSpriteEnable);
	UpdateFFState_Fog(pg, csv0d, ctl3);
	UpdateFFState_TextureStates(pg);
	UpdateFFState_Lighting(pg, csv0c);

	// Misc flags
	ffShaderState.Modes.NormalizeNormals = (csv0c & NV_PGRAPH_CSV0_C_NORMALIZATION_ENABLE) ? 1 : 0;

	// Write fixed function state to shader constants.
	// Must always upload because VP draws between FF draws overwrite the shared constant buffer.
	const int slotSize = 16;
	const int fixedFunctionStateSize = (sizeof(FixedFunctionVertexShaderState) + slotSize - 1) / slotSize;
	CxbxSetVertexShaderConstantF(0, (float*)&ffShaderState, fixedFunctionStateSize);
}

// ******************************************************************
// * Present / display helpers
// ******************************************************************

// Build a DXGI_GAMMA_CONTROL from the NV2A VGA DAC palette (256-entry CLUT),
// linearly interpolating into the 1025-entry DXGI curve.
static void CxbxApplyNV2AGamma(NV2AState* d)
{
	IDXGIOutput* pOutput = nullptr;
	if (FAILED(g_pSwapChain->GetContainingOutput(&pOutput))) return;

	auto clut = d->puserdac.palette;

	DXGI_GAMMA_CONTROL gammaControl = {};
	gammaControl.Scale  = { 1.0f, 1.0f, 1.0f };
	gammaControl.Offset = { 0.0f, 0.0f, 0.0f };

	for (int j = 0; j <= 1024; ++j) {
		float x    = (j / 1024.0f) * 255.0f;
		int   lo   = (int)x;
		int   hi   = lo < 255 ? lo + 1 : 255;
		float frac = x - (float)lo;

		auto lerp = [&](uint8_t a, uint8_t b) -> float {
			return (a + frac * (float)(b - a)) / 255.0f;
		};

		gammaControl.GammaCurve[j] = {
			lerp(clut[lo * 3 + 0], clut[hi * 3 + 0]),
			lerp(clut[lo * 3 + 1], clut[hi * 3 + 1]),
			lerp(clut[lo * 3 + 2], clut[hi * 3 + 2])
		};
	}

	pOutput->SetGammaControl(&gammaControl);
	pOutput->Release();
}

HRESULT CxbxPresent()
{
	LOG_INIT;
	CxbxEndScene();

	// Apply NV2A gamma LUT (PRMDIO VGA DAC palette) to DXGI output
	NV2AState* d = g_NV2A->GetDeviceState();
	if (d->puserdac.dirty) {
		d->puserdac.dirty = false;
		CxbxApplyNV2AGamma(d);
	}

	HRESULT hRet = g_pSwapChain->Present(0, g_bTearingSupported ? DXGI_PRESENT_ALLOW_TEARING : 0);
	DEBUG_D3DRESULT(hRet, "g_pSwapChain->Present");
	// Allow the next page tracker flush to use DISCARD (safe at frame boundary
	// since no draw calls from this frame are still referencing the buffer)
	CxbxPageTrackerOnPresent();
	CxbxBeginScene();
	return hRet;
}

HRESULT CxbxGetBackBuffer(ID3D11Texture2D** ppBackBuffer)
{
	return g_pSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), reinterpret_cast<void**>(ppBackBuffer));
}

