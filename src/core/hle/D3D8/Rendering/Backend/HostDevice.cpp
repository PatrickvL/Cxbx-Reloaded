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

// Variables only used in HostDevice.cpp
static HBRUSH g_hBgBrush = NULL; // Background Brush
static bool g_bIsFauxFullscreen = false;
static int g_iWireframe = 0; // wireframe toggle
static HWINEVENTHOOK g_parentEventHook = NULL; // WinEvent hook for tracking parent window moves

void CxbxSaveWindowStateForReboot()
{
	if (g_hEmuWindow == NULL)
		return;

	RECT rect;
	if (GetWindowRect(g_hEmuWindow, &rect)) {
		g_EmuShared->SetSavedWindowState(&rect, g_bIsFauxFullscreen);
	}
}

// Forward declarations (defined later in this file)
LRESULT WINAPI EmuMsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);
void ToggleFauxFullscreen(HWND hWnd);

// Reposition the owned popup to cover the GUI parent's client area exactly.
static void RepositionToParentClientArea()
{
	if (!CxbxKrnl_hEmuParent || !g_hEmuWindow || g_bIsFauxFullscreen)
		return;

	RECT clientRect;
	if (GetClientRect(CxbxKrnl_hEmuParent, &clientRect)) {
		MapWindowPoints(CxbxKrnl_hEmuParent, NULL, (LPPOINT)&clientRect, 2);
		SetWindowPos(g_hEmuWindow, NULL,
			clientRect.left, clientRect.top,
			clientRect.right - clientRect.left,
			clientRect.bottom - clientRect.top,
			SWP_NOZORDER | SWP_NOACTIVATE);
	}
}

// WinEvent hook callback: fires when the GUI parent window moves or resizes.
// WINEVENT_OUTOFCONTEXT means the callback runs on the emu's message pump thread.
static void CALLBACK ParentLocationChanged(
	HWINEVENTHOOK hWinEventHook, DWORD event, HWND hwnd,
	LONG idObject, LONG idChild, DWORD dwEventThread, DWORD dwmsEventTime)
{
	if (hwnd == CxbxKrnl_hEmuParent && idObject == OBJID_WINDOW) {
		RepositionToParentClientArea();
	}
}

// Forward declarations (defined in XbPushBuffer.cpp)
extern void D3D11_init_pgraph_plugins();

void EmuD3DInit()
{
	D3D11_init_pgraph_plugins(); // TODO : Hook more nv_dma_map() result uses in EmuNV2A_PGRAPH.cpp

	// Initialise CreateDevice Proxy Data struct
	{
		g_EmuCDPD.Adapter = nullptr; // Will be set below if user selected a non-default adapter
		g_EmuCDPD.DeviceType = (g_XBVideo.direct3DDevice == 0) ? D3D_DRIVER_TYPE_HARDWARE : D3D_DRIVER_TYPE_REFERENCE;
		// Enumerate DXGI adapters to select the one matching g_XBVideo.adapter index
		if (g_XBVideo.adapter != 0) {
			IDXGIFactory1 *pFactory = nullptr;
			if (SUCCEEDED(CreateDXGIFactory1(__uuidof(IDXGIFactory1), reinterpret_cast<void**>(&pFactory)))) {
				IDXGIAdapter *pAdapter = nullptr;
				if (SUCCEEDED(pFactory->EnumAdapters(g_XBVideo.adapter, &pAdapter))) {
					g_EmuCDPD.Adapter = pAdapter; // Takes ownership; released when device is destroyed
					// When specifying an adapter, driver type must be D3D_DRIVER_TYPE_UNKNOWN
					g_EmuCDPD.DeviceType = D3D_DRIVER_TYPE_UNKNOWN;
				}
				pFactory->Release();
			}
		}
	}
}

// cleanup Direct3D
void EmuD3DCleanup() {}

// window message processing thread
DWORD WINAPI EmuRenderWindow(LPVOID lpParam)
{
	CxbxSetThreadName("Cxbx Render Window");

   	// register window class
   	{
   	   	LOGBRUSH logBrush = {BS_SOLID, RGB(0,0,0)};

   	   	g_hBgBrush = CreateBrushIndirect(&logBrush);

   	   	WNDCLASSEX wc =
   	   	{
   	   	   	sizeof(WNDCLASSEX),
   	   	   	CS_CLASSDC,
   	   	   	EmuMsgProc,
   	   	   	0, 0, hActiveModule, // Was GetModuleHandle(nullptr),
			0, // TODO : LoadIcon(hmodule, ?)
   	   	   	LoadCursor(NULL, IDC_ARROW),
   	   	   	(HBRUSH)(g_hBgBrush), NULL,
   	   	   	"CxbxRender",
			nullptr
   	   	};

   	   	RegisterClassEx(&wc);
   	}

   	// create the window
   	{
		// Peform selection if running in GUI or kernel mode first.
		HWND hwndParent = (!CxbxKrnl_hEmuParent ? GetDesktopWindow() : CxbxKrnl_hEmuParent);
		DWORD dwStyle = WS_POPUP;
		RECT windowRect = { 0 };

		// Obtain the selected resolution from GUI or full desktop screen in kernel mode.
		if (!GetWindowRect(hwndParent, &windowRect)) {
			// Fall back resolution if failed
			windowRect = { 0, 0, 640, 480 };
		}

		// Then perform additional checks if not running in full screen.
		if (!g_XBVideo.bFullScreen) {

			// If running as kernel mode, force use the xbox's default resolution.
			if (!CxbxKrnl_hEmuParent) {
				// Xbox default resolution (standalone window is resizable by the way)
				windowRect.right = 640;
				windowRect.bottom = 480;
				dwStyle = WS_OVERLAPPEDWINDOW;
			}
			else {
				// Use WS_POPUP (owned by parent) instead of WS_CHILD.
				// This gives the emu process its own top-level window so it can be
				// the foreground process (getting OS scheduling priority, proper timer
				// resolution, and unthrottled DXGI Present without cross-process
				// DWM composition synchronization).
				dwStyle = WS_POPUP;
				// Position the popup over the GUI's client area (screen coordinates).
				RECT clientRect;
				if (GetClientRect(hwndParent, &clientRect)) {
					MapWindowPoints(hwndParent, NULL, (LPPOINT)&clientRect, 2);
					windowRect = clientRect;
				}
			}
		}

   	   	g_hEmuWindow = CreateWindow
   	   	(
   	   	   	"CxbxRender", "Cxbx-Reloaded",
   	   	   	dwStyle, 
			windowRect.left,
			windowRect.top,
			windowRect.right - windowRect.left,
			windowRect.bottom - windowRect.top,
   	   	   	hwndParent, nullptr, hActiveModule, // Was GetModuleHandle(nullptr),
   	   	   	nullptr
   	   	);
   	}

   	ShowWindow(g_hEmuWindow, ((CxbxKrnl_hEmuParent == 0) || g_XBVideo.bFullScreen) ? SW_SHOWDEFAULT : SW_SHOW);
   	UpdateWindow(g_hEmuWindow);

	// Restore window state from a previous reboot (secondary XBE load)
	{
		RECT savedRect;
		bool bSavedFaux, bValid;
		g_EmuShared->GetSavedWindowState(&savedRect, &bSavedFaux, &bValid);
		if (bValid) {
			g_EmuShared->ClearSavedWindowState();
			if (bSavedFaux) {
				// Restore faux fullscreen mode
				ToggleFauxFullscreen(g_hEmuWindow);
			} else if (!g_XBVideo.bFullScreen && !CxbxKrnl_hEmuParent) {
				// Restore standalone window position and size
				SetWindowPos(g_hEmuWindow, NULL,
					savedRect.left, savedRect.top,
					savedRect.right - savedRect.left,
					savedRect.bottom - savedRect.top,
					SWP_NOZORDER | SWP_NOACTIVATE);
			}
		}
	}

   	if(!g_XBVideo.bFullScreen && (CxbxKrnl_hEmuParent != NULL))
   	{
   	   	// Notify the GUI of our window handle (owned popup doesn't trigger
   	   	// automatic WM_PARENTNOTIFY/WM_CREATE like WS_CHILD did).
   	   	ipc_send_gui_update(IPC_UPDATE_GUI::WINDOW_HANDLE,
   	   	   	static_cast<unsigned int>(reinterpret_cast<uintptr_t>(g_hEmuWindow)));

   	   	// Install a WinEvent hook to track when the GUI parent window moves or
   	   	// resizes, so we can reposition the popup to stay aligned.
   	   	DWORD parentThreadId = GetWindowThreadProcessId(CxbxKrnl_hEmuParent, nullptr);
   	   	g_parentEventHook = SetWinEventHook(
   	   	   	EVENT_OBJECT_LOCATIONCHANGE, EVENT_OBJECT_LOCATIONCHANGE,
   	   	   	NULL, ParentLocationChanged,
   	   	   	0, parentThreadId,
   	   	   	WINEVENT_OUTOFCONTEXT);

   	   	// Claim foreground status for the emulation process.
   	   	// The GUI must have called AllowSetForegroundWindow for this to succeed.
   	   	SetForegroundWindow(g_hEmuWindow);
   	}

   	EmuLog(LOG_LEVEL::DEBUG, "Message-Pump thread is running.");

   	SetFocus(g_hEmuWindow);

   	SetEvent(*reinterpret_cast<PHANDLE>(lpParam));

   	// message processing loop
   	{
   	   	MSG msg;
   	   	BOOL bRet;
   	   	while((bRet = GetMessage(&msg, NULL, 0U, 0U)) != FALSE)
   	   	{
   	   	   	if(bRet == -1)
   	   	   	{
   	   	   	   	CxbxrAbort("GetMessage failed!");
   	   	   	}

   	   	   	TranslateMessage(&msg);
   	   	   	DispatchMessage(&msg);
   	   	}

   	   	CxbxrAbort(nullptr);
   	}

   	return 0;
}

// simple helper function
void ToggleFauxFullscreen(HWND hWnd)
{
   	if(g_XBVideo.bFullScreen) {
   	   	return;
   	}

   	// Remember last known position and style before go into faux full screen mode.
   	static RECT lRect = {};
   	static LONG gwl_style = {};

   	// Require to toggle before start process due to WM_SETFOCUS will get trigger earlier.
   	g_bIsFauxFullscreen = !g_bIsFauxFullscreen;
   	if (g_bIsFauxFullscreen) {
   	   	GetWindowRect(hWnd, &lRect);
   	   	gwl_style = GetWindowLong(hWnd, GWL_STYLE);
   	   	SetWindowLong(hWnd, GWL_STYLE, WS_POPUP);
   	   	if (CxbxKrnl_hEmuParent) {
   	   	   	// Window is already WS_POPUP (owned), just go topmost and maximize.
   	   	   	// No SetParent(NULL) needed since we're not a child window.
   	   	}
   	   	SetWindowPos(hWnd, HWND_TOPMOST, lRect.left, lRect.top, 0, 0, SWP_NOSIZE);
   	   	ShowWindow(hWnd, SW_MAXIMIZE);
   	}
   	else {
   	   	SetWindowLong(hWnd, GWL_STYLE, gwl_style);
   	   	if(CxbxKrnl_hEmuParent) {
   	   	   	// Restore popup position over the GUI's client area.
   	   	   	SetWindowPos(hWnd, HWND_NOTOPMOST, 0, 0, 0, 0, SWP_NOSIZE | SWP_NOMOVE);
   	   	   	RepositionToParentClientArea();
   	   	   	ShowWindow(hWnd, SW_SHOW);
   	   	}
   	   	else {
   	   	   	ShowWindow(hWnd, SW_RESTORE);
   	   	   	SetWindowPos(hWnd, HWND_NOTOPMOST, lRect.left, lRect.top, lRect.right - lRect.left, lRect.bottom - lRect.top, 0);
   	   	   	SetFocus(hWnd);
   	   	}
   	}
}

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

// rendering window message procedure
LRESULT WINAPI EmuMsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	const LRESULT imguiResult = ImGui_ImplWin32_WndProcHandler(hWnd, msg, wParam, lParam);
	if (imguiResult != 0) return imguiResult;

   	switch(msg)
   	{
   	   	case WM_DESTROY:
   	   	{
   	   	   	// Notify GUI that our window is gone (owned popup doesn't trigger
   	   	   	// automatic WM_PARENTNOTIFY/WM_DESTROY like WS_CHILD did).
   	   	   	if (CxbxKrnl_hEmuParent) {
   	   	   	   	ipc_send_gui_update(IPC_UPDATE_GUI::WINDOW_DESTROYED, 0);
   	   	   	}

   	   	   	// Unhook parent position tracking
   	   	   	if (g_parentEventHook) {
   	   	   	   	UnhookWinEvent(g_parentEventHook);
   	   	   	   	g_parentEventHook = NULL;
   	   	   	}

   	   	   	CxbxReleaseCursor();
   	   	   	DeleteObject(g_hBgBrush);
   	   	   	PostQuitMessage(0);
   	   	   	return S_OK; // = 0
   	   	}
   	   	break;

		case WM_PAINT:
		{
			if (g_CxbxPrintUEM)
			{
				DrawUEM(hWnd);
			}
			else
			{
				DefWindowProc(hWnd, msg, wParam, lParam);
			}
		}
		break;

		case WM_COMMAND:
		{
			switch (LOWORD(wParam))
			{
				case ID_SYNC_CONFIG_LOGGING:
				{
					// Sync run-time config log settings from GUI process.
					log_sync_config();
					log_generate_active_filter_output(CXBXR_MODULE::CXBXR);
				}
				break;

				case ID_SYNC_CONFIG_INPUT:
				{
					SDL_Event UpdateInputEvent;
					SDL_memset(&UpdateInputEvent, 0, sizeof(SDL_Event));
					UpdateInputEvent.type = Sdl::UpdateInputEvent_t;
					UpdateInputEvent.user.data1 = new int(lParam);
					SDL_PushEvent(&UpdateInputEvent);
				}
				break;

				case ID_SYNC_TIME_CHANGE:
				{
					// Sent by the GUI when it detects WM_TIMECHANGE
					xbox::KeSystemTimeChanged.test_and_set();
				}
				break;

				default:
					break;
			}
		}
		break;

		case WM_INPUT_DEVICE_CHANGE:
		{
			// sent by rawinput when it detects changes for the registered device types

			if (wParam == GIDC_ARRIVAL) {
				g_InputDeviceManager.HotplugHandler(false);
			}
			return DefWindowProc(hWnd, msg, wParam, lParam);
		}
		break;

		case WM_TIMECHANGE:
		{
			// NOTE: this is only received if the loader was launched from the command line without the GUI
			xbox::KeSystemTimeChanged.test_and_set();
			return DefWindowProc(hWnd, msg, wParam, lParam);
		}
		break;

   	   	case WM_SYSKEYDOWN:
   	   	{
   	   	   	if(wParam == VK_RETURN)
   	   	   	{
   	   	   	   	ToggleFauxFullscreen(hWnd);
   	   	   	}
   	   	   	else if(wParam == VK_F4)
   	   	   	{
   	   	   	   	PostMessage(hWnd, WM_CLOSE, 0, 0);
   	   	   	}
   	   	   	// NOTE: Windows does not send F10 key message to WM_KEYDOWN.
   	   	   	// Source: https://docs.microsoft.com/en-us/windows/desktop/inputdev/wm-syskeydown
   	   	   	else if(wParam == VK_F10)
   	   	   	{
   	   	   	   	ToggleFauxFullscreen(hWnd);
   	   	   	}
   	   	   	else
   	   	   	{
   	   	   	   	return DefWindowProc(hWnd, msg, wParam, lParam);
   	   	   	}
   	   	}
   	   	break;

   	   	case WM_KEYDOWN:
   	   	{
   	   	   	/*! disable fullscreen if we are set to faux mode, and faux fullscreen is active */
   	   	   	if(wParam == VK_ESCAPE)
   	   	   	{
   	   	   	   	if(g_XBVideo.bFullScreen)
   	   	   	   	{
   	   	   	   	   	SendMessage(hWnd, WM_CLOSE, 0, 0);
   	   	   	   	}
   	   	   	   	else if(g_bIsFauxFullscreen)
   	   	   	   	{
   	   	   	   	   	ToggleFauxFullscreen(hWnd);
   	   	   	   	}
   	   	   	}
   	   	   	else if (wParam == VK_F1)
   	   	   	{
				if (g_renderbase)
					g_renderbase->ToggleImGui();

   	   	   	   	CxbxUpdateCursor();
   	   	   	}

   	   	   	else if (wParam == VK_F3)
   	   	   	{
   	   	   	   	g_bClipCursor = !g_bClipCursor;
   	   	   	   	g_EmuShared->SetClipCursorFlag(g_bClipCursor);

   	   	   	   	if (g_bClipCursor) {
   	   	   	   	   	CxbxClipCursor(hWnd);
   	   	   	   	}
   	   	   	   	else {
   	   	   	   	   	CxbxReleaseCursor();
   	   	   	   	}
   	   	   	}
   	   	   	else if (wParam == VK_F6)
   	   	   	{
   	   	   	   	// For some unknown reason, F6 isn't handled in WndMain::WndProc
   	   	   	   	// sometimes, so detect it and stop emulation from here too :
   	   	   	   	SendMessage(hWnd, WM_CLOSE, 0, 0); // See StopEmulation();
   	   	   	}
			else if (wParam == VK_F7)
			{
				// F7 key - previously toggled passthrough HLSL (removed)
			}
   	   	   	else if(wParam == VK_F8)
   	   	   	{
   	   	   	   	g_bPrintfOn = !g_bPrintfOn;
   	   	   	   	g_EmuShared->SetIsKrnlLogEnabled(g_bPrintfOn);
   	   	   	   	LOG_THREAD_INIT;
   	   	   	   	std::cout << _logThreadPrefix << g_EnumModules2String[to_underlying(CXBXR_MODULE::CXBXR)] << "Enable log is " << g_bPrintfOn << std::endl;
   	   	   	   	ipc_send_gui_update(IPC_UPDATE_GUI::LOG_ENABLED, static_cast<UINT>(g_bPrintfOn));
   	   	   	}
			else if (wParam == VK_F9)
			{
				// Toggle frame-limiting
				g_bHack_UnlockFramerate = !g_bHack_UnlockFramerate;
			}
   	   	   	else if(wParam == VK_F11)
   	   	   	{
   	   	   	   	if (g_iWireframe++ == 2) {
   	   	   	   	   	g_iWireframe = 0;
   	   	   	   	}

   	   	   	   	// TODO: implement debug wireframe via PGRAPH register override
   	   	   	}
   	   	   	else
   	   	   	{
   	   	   	   	return DefWindowProc(hWnd, msg, wParam, lParam);
   	   	   	}
   	   	}
   	   	break;

   	   	case WM_SIZE:
   	   	{
   	   	   	switch(wParam)
   	   	   	{
   	   	   	   	case SIZE_MINIMIZED:
   	   	   	   	{
   	   	   	   	   	if(g_XBVideo.bFullScreen)
   	   	   	   	   	   	CxbxrAbort(nullptr);
   	   	   	   	}
   	   	   	   	break;
   	   	   	}

   	   	   	if (g_bClipCursor) {
   	   	   	   	CxbxClipCursor(hWnd);
   	   	   	}
   	   	}
   	   	break;

   	   	case WM_MOVE:
   	   	{
   	   	   	if (g_bClipCursor) {
   	   	   	   	CxbxClipCursor(hWnd);
   	   	   	}
   	   	}
   	   	break;

   	   	case WM_MOUSELEAVE:
   	   	{
   	   	   	DInput::mo_leave_wnd = true;
   	   	   	g_bIsTrackingMoLeave = false;
   	   	   	g_bIsTrackingMoMove = true;
   	   	   	CxbxUpdateCursor(true);
   	   	}
   	   	break;

   	   	case WM_MOUSEMOVE:
   	   	{
   	   	   	if (g_bClipCursor) {
   	   	   	   	CxbxClipCursor(hWnd);
   	   	   	}

   	   	   	if (!g_bIsTrackingMoLeave) {
   	   	   	   	TRACKMOUSEEVENT tme;
   	   	   	   	tme.cbSize = sizeof(TRACKMOUSEEVENT);
   	   	   	   	tme.hwndTrack = hWnd;
   	   	   	   	tme.dwFlags = TME_LEAVE;
   	   	   	   	TrackMouseEvent(&tme);
   	   	   	   	g_bIsTrackingMoLeave = true;
   	   	   	   	CxbxUpdateCursor();

   	   	   	   	if (g_bIsTrackingMoMove) {
   	   	   	   	   	DInput::mo_leave_wnd = false;
   	   	   	   	   	g_bIsTrackingMoMove = false;
   	   	   	   	}
   	   	   	}
   	   	}
   	   	break;

   	   	case WM_ACTIVATEAPP:
   	   	{
   	   	   	if (wParam == FALSE) {
   	   	   	   	DInput::mo_leave_wnd = true;
   	   	   	} else {
   	   	   	   	DInput::mo_leave_wnd = false;
   	   	   	}
   	   	}
   	   	break;

   	   	case WM_CLOSE:
   	   	   	CxbxReleaseCursor();
   	   	   	DestroyWindow(hWnd);
   	   	   	CxbxrShutDown();
   	   	   	break;

   	   	case WM_SETFOCUS:
   	   	{
   	   	   	// With WS_POPUP, let the emu window keep focus so the emulation
   	   	   	// process stays foreground (better scheduling + timer resolution).
   	   	}
   	   	break;

   	   	case WM_SETCURSOR:
   	   	{
   	   	   	if(g_XBVideo.bFullScreen || g_bIsFauxFullscreen)
   	   	   	{
   	   	   	   	SetCursor(NULL);
   	   	   	   	return S_OK; // = Is not part of D3D8 handling.
   	   	   	}
   	   	   	return DefWindowProc(hWnd, msg, wParam, lParam);
   	   	}
   	   	break;

		case WM_CXBXR_RUN_ON_MESSAGE_THREAD:
		{
			auto func = reinterpret_cast<const std::function<void()>*>(wParam);
			(*func)();
			return 0;
		}

   	   	default:
   	   	   	return DefWindowProc(hWnd, msg, wParam, lParam);
   	}

   	return S_OK; // = Is not part of D3D8 handling.
}

void UpdateDepthStencilFlags(ID3D11Texture2D *pDepthStencilSurface)
{
	g_bHasDepth = false;
	g_bHasStencil = false;
	if (pDepthStencilSurface != nullptr) {
		D3D11_TEXTURE2D_DESC Desc;
		pDepthStencilSurface->GetDesc(&Desc);

		switch (Desc.Format) {
		// Under D3D11, multiple D3D9 depth formats map to the same DXGI format,
		// so we switch on the unique DXGI values only :
		case DXGI_FORMAT_D16_UNORM: // Covers D16, D15S1 (stencil bit lost), D16_LOCKABLE
			g_bHasDepth = true;
			break;
		case DXGI_FORMAT_D24_UNORM_S8_UINT: // Covers D24S8, D24X8, D24X4S4 (all have stencil in D3D11)
			g_bHasDepth = true;
			g_bHasStencil = true;
			break;
		case DXGI_FORMAT_D32_FLOAT: // Covers D32
			g_bHasDepth = true;
			break;
		}
	}
}

void SetupPresentationParameters
(
   	const xbox::X_D3DPRESENT_PARAMETERS     *pXboxPresentationParameters
)
{
   	auto& params = g_EmuCDPD.HostPresentationParameters;

   	params.Windowed = !g_XBVideo.bFullScreen;

    (void)pXboxPresentationParameters; // Presentation interval handled by PGRAPH flip_stall

   	// retrieve resolution from configuration
   	char szBackBufferFormat[16] = {};
   	const char* resolution = g_XBVideo.szVideoResolution;
   	if (4 != sscanf(resolution, "%u x %u %*dbit %s (%u hz)",
   	   	&params.BackBufferWidth,
   	   	&params.BackBufferHeight,
   	   	szBackBufferFormat,
   	   	&params.FullScreen_RefreshRateInHz)) {
   		if (pXboxPresentationParameters) {
   	   		EmuLog(LOG_LEVEL::DEBUG, "EmuCreateDeviceProxy: Couldn't parse resolution : %s. Using Xbox Default (%d, %d @ %uhz)", resolution,
   	   	   		pXboxPresentationParameters->BackBufferWidth, pXboxPresentationParameters->BackBufferHeight,
   	   	   		pXboxPresentationParameters->FullScreen_RefreshRateInHz);
   	   		params.BackBufferWidth = pXboxPresentationParameters->BackBufferWidth;
   	   		params.BackBufferHeight = pXboxPresentationParameters->BackBufferHeight;
   	   		// Field rendering uses half-height backbuffers (one per field); the host
   	   		// display should be full frame height so both fields can be presented.
   	   		if (pXboxPresentationParameters->Flags & X_D3DPRESENTFLAG_FIELD) {
   	   			params.BackBufferHeight *= 2;
   	   		}
   	   		params.FullScreen_RefreshRateInHz = pXboxPresentationParameters->FullScreen_RefreshRateInHz;
   		} else {
   			// Early init without Xbox params — use sensible defaults
   			EmuLog(LOG_LEVEL::DEBUG, "SetupPresentationParameters: No Xbox params, using 640x480 @ 60hz defaults");
   			params.BackBufferWidth = 640;
   			params.BackBufferHeight = 480;
   			params.FullScreen_RefreshRateInHz = 60;
   		}
   	}

   	if (params.Windowed) {
   	   	// Refresh rate must be 0 in windowed mode, as is documented
   	   	params.FullScreen_RefreshRateInHz = 0;
   	}
}

void DetermineSupportedD3DFormats
(
)
{
   	memset(g_bSupportsFormatSurface, false, sizeof(g_bSupportsFormatSurface));
   	memset(g_bSupportsFormatSurfaceRenderTarget, false, sizeof(g_bSupportsFormatSurfaceRenderTarget));
   	memset(g_bSupportsFormatSurfaceDepthStencil, false, sizeof(g_bSupportsFormatSurfaceDepthStencil));
   	memset(g_bSupportsFormatTexture, false, sizeof(g_bSupportsFormatTexture));
   	memset(g_bSupportsFormatTextureRenderTarget, false, sizeof(g_bSupportsFormatTextureRenderTarget));
   	memset(g_bSupportsFormatTextureDepthStencil, false, sizeof(g_bSupportsFormatTextureDepthStencil));
   	memset(g_bSupportsFormatVolumeTexture, false, sizeof(g_bSupportsFormatVolumeTexture));
   	memset(g_bSupportsFormatCubeTexture, false, sizeof(g_bSupportsFormatCubeTexture));
   	for (int X_Format = xbox::X_D3DFMT_FIRST; X_Format <= xbox::X_D3DFMT_LAST; X_Format++) {
   	   	// Only process Xbox formats that are directly mappable to host
		DXGI_FORMAT PCFormat;
   	   	if (!EmuXBFormatRequiresConversion((xbox::X_D3DFORMAT)X_Format, /*&*/PCFormat)) {
   	   	   	// Convert the Xbox format into host format (without warning, thanks to the above restriction)
   	   	   	PCFormat = EmuXB2PC_D3DFormat((xbox::X_D3DFORMAT)X_Format);
   	   	   	if (PCFormat != EMUFMT_UNKNOWN) {
   	   	   	   	// Index with Xbox D3DFormat, because host FourCC codes are too big to be used as indices
   	   	   	   	UINT FormatSupport = 0;
				g_pD3DDevice->CheckFormatSupport(PCFormat, &FormatSupport);

   	   	   	   	g_bSupportsFormatSurface[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_TEXTURE2D;
   	   	   	   	g_bSupportsFormatSurfaceRenderTarget[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_RENDER_TARGET;
   	   	   	   	g_bSupportsFormatSurfaceDepthStencil[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_DEPTH_STENCIL;
   	   	   	   	// In D3D11, surfaces are textures, so the same format support applies
   	   	   	   	// Regular textures must support both creation AND shader sampling;
   	   	   	   	// formats like YUY2 may support TEXTURE2D but not SHADER_SAMPLE on some hardware,
   	   	   	   	// which causes CreateShaderResourceView to fail later.
   	   	   	   	g_bSupportsFormatTexture[X_Format] = (FormatSupport & D3D11_FORMAT_SUPPORT_TEXTURE2D)
   	   	   	   	   	&& (FormatSupport & D3D11_FORMAT_SUPPORT_SHADER_SAMPLE);
   	   	   	   	g_bSupportsFormatTextureRenderTarget[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_RENDER_TARGET;
   	   	   	   	g_bSupportsFormatTextureDepthStencil[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_DEPTH_STENCIL;

   	   	   	   	g_bSupportsFormatVolumeTexture[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_TEXTURE3D;
   	   	   	   	g_bSupportsFormatCubeTexture[X_Format] = FormatSupport & D3D11_FORMAT_SUPPORT_TEXTURECUBE;
   	   	   	}
   	   	}
   	}
}

void UpdateHostBackBufferDesc()
{
   	ID3D11Texture2D *pCurrentHostBackBuffer = nullptr;

   	HRESULT hRet = CxbxGetBackBuffer(&pCurrentHostBackBuffer);
   	if (hRet != S_OK) {
   	   	CxbxrAbort("Unable to get host backbuffer surface");
   	}

   	pCurrentHostBackBuffer->GetDesc(&g_HostBackBufferDesc);

   	pCurrentHostBackBuffer->Release();
}

void SetAspectRatioScale(const xbox::X_D3DPRESENT_PARAMETERS* pPresentationParameters)
{
   	// NOTE: Some games use anamorphic widesceen (expecting a 4:3 surface to be displayed at 16:9)
   	// For those, we *lie* about the default width, for the scaler
   	// 720p / 1080i are *always* widescreen, and will have the correct backbuffer size, so we only
   	// apply this 'hack' for non-hd resolutions
   	g_AspectRatioScaleWidth = pPresentationParameters->BackBufferWidth;
   	g_AspectRatioScaleHeight = pPresentationParameters->BackBufferHeight;

   	if (pPresentationParameters->Flags & X_D3DPRESENTFLAG_WIDESCREEN &&
   	   	pPresentationParameters->BackBufferHeight < 720) {
   	   	// Lie and pretend we are 1280x720, this works because this ratio is only used in calculations
   	   	// and not used as actual raw input values
   	   	g_AspectRatioScaleWidth = 1280;
   	   	g_AspectRatioScaleHeight = 720;
   	}
   	
   	const auto imageAspect = (float)g_AspectRatioScaleWidth / (float)g_AspectRatioScaleHeight;
   	const auto screenAspect = (float)g_HostBackBufferDesc.Width / (float)g_HostBackBufferDesc.Height;
   	g_AspectRatioScale = screenAspect > imageAspect ? (float)g_HostBackBufferDesc.Height / (float)g_AspectRatioScaleHeight : (float)g_HostBackBufferDesc.Width / (float)g_AspectRatioScaleWidth;
}

#define CXBX_D3DMULTISAMPLE_XSCALE(type) (((type) & xbox::X_D3DMULTISAMPLE_XSCALE_MASK) >> xbox::X_D3DMULTISAMPLE_XSCALE_SHIFT)
#define CXBX_D3DMULTISAMPLE_YSCALE(type) (((type) & xbox::X_D3DMULTISAMPLE_YSCALE_MASK) >> xbox::X_D3DMULTISAMPLE_YSCALE_SHIFT)

void SetXboxMultiSampleType(xbox::X_D3DMULTISAMPLE_TYPE value)
{
	// Validate & correct input, to detect test cases and avoid trouble when using g_Xbox_MultiSampleType :
	if (value == 0) {
		// Correcting zero to X_D3DMULTISAMPLE_NONE
		value = xbox::X_D3DMULTISAMPLE_NONE;
	}
	if (value & ~xbox::X_D3DMULTISAMPLE_KNOWN_MASK) {
		LOG_TEST_CASE("Removing unknown X_D3DMULTISAMPLE bits");
		value &= xbox::X_D3DMULTISAMPLE_KNOWN_MASK;
	}
	if (CXBX_D3DMULTISAMPLE_XSCALE(value) == 0) {
		LOG_TEST_CASE("Correcting XSCALE 0 to 1");
		value |= 1 << xbox::X_D3DMULTISAMPLE_XSCALE_SHIFT;
	}
	if (CXBX_D3DMULTISAMPLE_YSCALE(value) == 0) {
		LOG_TEST_CASE("Correcting YSCALE 0 to 1");
		value |= 1 << xbox::X_D3DMULTISAMPLE_YSCALE_SHIFT;
	}
	if ((CXBX_D3DMULTISAMPLE_XSCALE(value) == 1) && (CXBX_D3DMULTISAMPLE_YSCALE(value) == 1)) {
		if (value & xbox::X_D3DMULTISAMPLE_ALGO_MASK) {
			LOG_TEST_CASE("Removing algorithm for 1 by 1 scale");
			value &= ~xbox::X_D3DMULTISAMPLE_ALGO_MASK;
		}
		if (value & xbox::X_D3DMULTISAMPLE_SAMPLING_MASK) {
			LOG_TEST_CASE("Removing sampling for 1 by 1 scale");
			value &= ~xbox::X_D3DMULTISAMPLE_SAMPLING_MASK;
		}
	}
	if (value == xbox::X_D3DMULTISAMPLE_9_SAMPLES_MULTISAMPLE_GAUSSIAN) {
		LOG_TEST_CASE("X_D3DMULTISAMPLE_9_SAMPLES_MULTISAMPLE_GAUSSIAN"); // This used to use scale 1.5f
	}
	// Set the now-validated g_Xbox_MultiSampleType value :
	g_Xbox_MultiSampleType = value;
}

float GetMultiSampleOffsetDelta()
{
	// TODO : What impact does X_D3DMULTISAMPLE_SAMPLING_SUPER have on offset?
	return (g_Xbox_MultiSampleType & xbox::X_D3DMULTISAMPLE_SAMPLING_MULTI) ? 0.0f : 0.5f;
}

void GetMultiSampleOffset(float& xOffset, float& yOffset)
{
	xOffset = yOffset = GetMultiSampleOffsetDelta();
}

// Get the raw X and Y scale of the multisample mode
// Both MSAA and SSAA result in increased rendertarget size
void GetMultiSampleScaleRaw(float& xScale, float& yScale) {
	xScale = static_cast<float>(CXBX_D3DMULTISAMPLE_XSCALE(g_Xbox_MultiSampleType));
	yScale = static_cast<float>(CXBX_D3DMULTISAMPLE_YSCALE(g_Xbox_MultiSampleType));
}

// Get the "screen" scale factors
// Multisampling modes and backbufferscale affect the screen coordinate space
// Titles can render pre-transformed vertices to screen space (using XYZRHW vertex position data or otherwise)
// so we need to know the space they are in to interpret it correctly
void GetScreenScaleFactors(float& scaleX, float& scaleY) {
	scaleX = 1;
	scaleY = 1;

	// With fixed-function mode, titles don't have to account for these scale factors,
	// so we don't have reverse them.
	// Test cases:
	// Fixed-func XYZRHW, title does not apply SSAA scale:
	// - Shenmue II (Menu)
	// Fixed-func XYZRHW, title does not apply backbuffer scale:
	// - Antialias sample(background gradient)
	// Only apply scale factors for PROGRAM mode (shader programs account for them).
	if (NV2AIsFixedFunctionMode()) {
		return;
	}

	// Example:
	// NFS HP2 renders in-game at 640*480
	// The title uses MSAA, which increases the rendertarget size, but leaves the screen scale unaffected
	// The backbuffer scale is (0.75, 1), which contributes to the screen scale
	// So the Xbox expects vertices in 480*480 coordinate space

	// SSAA increases the screen scale (but MSAA does not)
	// Read from PGRAPH ANTIALIASING register (authoritative, no HLE dependency)
	bool isMultiSampleEnabled = false;
	if (g_NV2A) {
		PGRAPHState* pg = &g_NV2A->GetDeviceState()->pgraph;
		isMultiSampleEnabled = (pg->regs[RI(NV_PGRAPH_ANTIALIASING)] & NV_PGRAPH_ANTIALIASING_ENABLE) != 0;
	}
	bool isSuperSampleMode = g_Xbox_MultiSampleType & xbox::X_D3DMULTISAMPLE_SAMPLING_SUPER;

	// Apply multisample scale if supersampling is enabled
	// Test cases:
	// Antialias sample (text overlay is drawn with antialiasing disabled)
	if (isMultiSampleEnabled && isSuperSampleMode) {
		GetMultiSampleScaleRaw(scaleX, scaleY);
	}

}

// GetRenderTargetRawDimensions, GetBackBufferPixelDimensions removed:
// dead code (never called). PGRAPH NV2ASurfaceState.clipWidth/clipHeight is used instead.

// Get the base rendertarget dimensions excluding multisample scaling
// e.g. a raw 1280*960 rendertarget with 2x MSAA would be have a base 640*480
void GetRenderTargetBaseDimensions(float& x, float& y) {
	// Read render target dimensions from PGRAPH surface clip (replaces HLE g_pXbox_RenderTarget lookup)
	auto surf = NV2AGetSurfaceState();
	x = (float)surf.clipWidth;
	y = (float)surf.clipHeight;

	// NV2A clip registers (SURFACECLIPX/Y) contain LOGICAL dimensions.
	// The AA factor (CENTER_CORNER_2, SQUARE_OFFSET_4) is a separate
	// hardware register that scales the surface physically — it is NOT
	// baked into the clip rect. So no AA division needed here.
	// For SSAA, the D3D runtime adjusts the viewport transform via
	// GetScreenScaleFactors, which multiplies by the AA factor.
}
