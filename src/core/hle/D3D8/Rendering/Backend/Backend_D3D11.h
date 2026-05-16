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
#ifndef BACKEND_D3D11_H
#define BACKEND_D3D11_H

#include "core\hle\D3D8\XbD3D8Types.h"
#include "../RenderGlobals.h" // resource_key_t, resource_key_hash
#include <vector>
#include <unordered_map>
#include <wrl/client.h>

struct _CxbxVertexDeclaration;
typedef struct _CxbxVertexDeclaration CxbxVertexDeclaration;

// ******************************************************************
// * D3D11 device globals
// ******************************************************************
extern ID3D11Device                *g_pD3DDevice;
extern IDXGISwapChain              *g_pSwapChain;
extern bool                         g_bTearingSupported;
extern ID3D11DeviceContext         *g_pD3DDeviceContext;
extern ID3D11RenderTargetView      *g_pD3DBackBufferView;
extern ID3D11DepthStencilView      *g_pD3DDepthStencilView;
extern ID3D11RenderTargetView      *g_pD3DCurrentRTV;
extern ID3D11Texture2D             *g_pD3DDepthStencilBuffer;
extern ID3D11Texture2D             *g_pD3DBackBufferSurface;
extern ID3D11Texture2D             *g_pD3DCurrentHostRenderTarget;

// PGRAPH-tracked backbuffer: set by CxbxD3D11UpdateRenderTargetFromPGRAPH
// when the first color RT is bound (which is always the backbuffer from CreateDevice).
extern ID3D11Texture2D             *g_pHostPgraphBackBuffer;
extern UINT                         g_PgraphBackBufferWidth;
extern UINT                         g_PgraphBackBufferHeight;
extern D3D11_TEXTURE2D_DESC         g_HostBackBufferDesc;
extern ID3D11Query                 *g_pHostQueryWaitForIdle;
void CxbxResetPgraphSurfaceTracking();
ID3D11Texture2D* CxbxLookupPgraphRTByOffset(xbox::addr_xt offset);
void CxbxPgraphRTCacheEvict();
void CxbxInvalidatePgraphRTBinding();

// ******************************************************************
// * Constant buffer sizing
// ******************************************************************
static const UINT CXBX_D3D11_VS_CB_SLOT = 0;
static const UINT CXBX_D3D11_VS_CB_COUNT = 256;
static const UINT CXBX_D3D11_PS_CB_SLOT = 0;
static const UINT CXBX_D3D11_PS_PGREGS_SRV_SLOT = 12; // ByteAddressBuffer (mirror SRV) : register(t12) — PGRAPH at offset 0x04000000
static const UINT CXBX_D3D11_VS_PGREGS_SRV_SLOT = 12; // Same mirror SRV shared with VS : register(t12)
static const UINT CXBX_D3D11_VS_XFPR_SRV_SLOT = 5;     // StructuredBuffer<uint4> g_XFPR : register(t5) — NV2A XFPR (Transform Program RAM)

// ******************************************************************
// * Vertex defaults buffer — provides all 16 TEXCOORD attributes via
// * a zero-stride vertex buffer bound to a dedicated input slot.
// * Non-streamed attributes read NV2A "sticky" values from this buffer.
// * This satisfies DXVK's requirement that every ISGN entry has a
// * matching input layout element.
// ******************************************************************
static const UINT CXBX_D3D11_VERTEX_DEFAULTS_SLOT = 16; // Input slot for zero-stride defaults buffer
extern ID3D11Buffer *g_pD3D11VertexDefaultsBuffer;
// Constant buffers (created in RenderGlobals.cpp device init, used by Backend_D3D11.cpp)
extern ID3D11Buffer *g_pD3D11VSConstantBuffer;

// ******************************************************************
// * D3D11 state descriptors (modified by RenderStates.cpp, applied by CxbxD3D11ApplyDirtyStates)
// ******************************************************************
extern D3D11_RASTERIZER_DESC    g_D3D11RasterizerDesc;
extern D3D11_DEPTH_STENCIL_DESC g_D3D11DepthStencilDesc;
extern D3D11_BLEND_DESC         g_D3D11BlendDesc;

// Dirty flags to trigger D3D11 state object recreation
extern bool  g_bD3D11RasterizerStateDirty;
extern bool  g_bD3D11DepthStencilStateDirty;
extern bool  g_bD3D11BlendStateDirty;

// Additional state parameters passed to OMSet* calls
extern UINT  g_D3D11StencilRef;
extern FLOAT g_D3D11BlendFactor[4];
extern UINT  g_D3D11SampleMask;

// ******************************************************************
// * Depth format helpers for typeless creation + SRV/DSV views
// ******************************************************************

// Map a depth format to its typeless equivalent for dual-bind (DSV+SRV) texture creation.
inline DXGI_FORMAT GetTypelessDepthFormat(DXGI_FORMAT depthFormat) {
	switch (depthFormat) {
		case DXGI_FORMAT_D16_UNORM:           return DXGI_FORMAT_R16_TYPELESS;
		case DXGI_FORMAT_D24_UNORM_S8_UINT:   return DXGI_FORMAT_R24G8_TYPELESS;
		default: return depthFormat;
	}
}

// Map a typeless or depth format to the SRV-compatible format for sampling depth textures.
inline DXGI_FORMAT GetDepthSRVFormat(DXGI_FORMAT format) {
	switch (format) {
		case DXGI_FORMAT_R16_TYPELESS:
		case DXGI_FORMAT_D16_UNORM:           return DXGI_FORMAT_R16_UNORM;
		case DXGI_FORMAT_R24G8_TYPELESS:
		case DXGI_FORMAT_D24_UNORM_S8_UINT:   return DXGI_FORMAT_R24_UNORM_X8_TYPELESS;
		default: return format;
	}
}

// Map a typeless format back to the typed depth format for DSV creation.
inline DXGI_FORMAT GetDepthDSVFormat(DXGI_FORMAT format) {
	switch (format) {
		case DXGI_FORMAT_R16_TYPELESS:        return DXGI_FORMAT_D16_UNORM;
		case DXGI_FORMAT_R24G8_TYPELESS:      return DXGI_FORMAT_D24_UNORM_S8_UINT;
		default: return format;
	}
}

// Check if a format is a depth/stencil or typeless-depth format.
inline bool IsDepthFormat(DXGI_FORMAT format) {
	switch (format) {
		case DXGI_FORMAT_D16_UNORM:
		case DXGI_FORMAT_D24_UNORM_S8_UINT:
		case DXGI_FORMAT_R16_TYPELESS:
		case DXGI_FORMAT_R24G8_TYPELESS:
			return true;
		default:
			return false;
	}
}

// ******************************************************************
// * D3D11 backend functions
// ******************************************************************

// Compile blit shaders and create sampler states (called once at device init)
void CxbxD3D11InitBlit();

// Create a D3D11 constant buffer (DYNAMIC or DEFAULT).
HRESULT CxbxD3D11CreateConstantBuffer(UINT byteWidth, bool bDynamic, ID3D11Buffer** ppBuffer);

// Release all D3D11 backend resources (blit shaders, samplers, constant buffers, state objects)
void CxbxD3D11ReleaseBackendResources();
void CxbxReleaseOverlayResources();

// Read NV2A PGRAPH registers and update D3D11 blend/depth-stencil/rasterizer descriptors.
struct PGRAPHState; // forward decl
void CxbxD3D11UpdatePipelineStateFromPGRAPH(PGRAPHState *pg);

// Read PGRAPH texture registers (TEXADDRESS, TEXFILTER, TEXCTL0, BORDERCOLOR)
// and create/bind D3D11 sampler states. Replaces XboxTextureStates.Apply() for samplers.
void CxbxD3D11UpdateSamplersFromPGRAPH(PGRAPHState *pg);

// Read viewport offset/scale and window clip from PGRAPH and set D3D11 viewport/scissor.
// Called from CxbxUpdateNativeD3DResources() to set viewport from PGRAPH register data.
void CxbxD3D11UpdateViewportFromPGRAPH(PGRAPHState *pg);

// Read PGRAPH surface_color/zeta offsets and rebind D3D11 render target / depth-stencil
// if they differ from what's currently bound. Creates host surfaces directly from PGRAPH state.
void CxbxD3D11UpdateRenderTargetFromPGRAPH(PGRAPHState *pg);

// Recreate D3D11 state objects that have been marked dirty, and flush constant buffers
void CxbxD3D11ApplyDirtyStates();

// Flush vertex shader constant buffer to GPU if dirty
void CxbxD3D11FlushVertexShaderConstants();

// ******************************************************************
// * Register combiner interpreter (PS ubershader)
// ******************************************************************

// Compile the RC interpreter ubershader and create its constant buffer.
// Called lazily on first use.  Returns true if compilation succeeded.
bool CxbxD3D11InitRCInterpreter();

// Upload Xbox register combiner state to the RC interpreter constant buffer
// and bind it to the pixel shader stage.  Called each frame when the
// interpreter is active, before the draw call.
void CxbxD3D11UploadRCInterpreterState();

// ******************************************************************
// * Vertex shader interpreter (VS ubershader)
// ******************************************************************

// When true, use the vertex shader interpreter ubershader instead of
// per-program recompiled vertex shaders. Set via user option or debug toggle.
extern bool g_bUseVSInterpreter;

// Compile the VS interpreter ubershader and create its constant buffer.
// Called lazily on first use.  Returns true if compilation succeeded.
bool CxbxD3D11InitVSInterpreter();

// D3D11 blit: copy source texture region to dest texture region with optional filtering
// Fast path for same-size copies, shader-based path for scaled copies
HRESULT CxbxD3D11Blt(
	ID3D11Texture2D* pSrc, const RECT* pSrcRect,
	ID3D11Texture2D* pDst, const RECT* pDstRect,
	D3DTEXTUREFILTERTYPE Filter);

// Bind/unbind thick line geometry shader around line draw calls.
// Call CxbxBindThickLineGS before and CxbxUnbindThickLineGS after direct
// Draw/DrawIndexed calls that may use line primitives (UP draw paths).
void CxbxBindThickLineGS(uint32_t primitiveMode);
void CxbxUnbindThickLineGS(uint32_t primitiveMode);

// GPU-accelerated texture unswizzle via compute shader.
// Returns true if the CS path was used, false if caller should fall back to CPU.
// Texture must be D3D11_USAGE_DEFAULT with D3D11_BIND_UNORDERED_ACCESS.
bool CxbxD3D11UnswizzleTexture(
	ID3D11Texture2D* pTexture,
	const void* pSwizzledSrc,
	UINT width,
	UINT height,
	UINT bpp,
	DXGI_FORMAT format);

// Index conversion mode constants
#define CXBX_INDEX_CONVERT_QUAD_CW  0
#define CXBX_INDEX_CONVERT_QUAD_CCW 1
#define CXBX_INDEX_CONVERT_FAN      2

// GPU-accelerated index buffer topology conversion via compute shader.
// Converts quad-list or triangle-fan indices to triangle-list on the GPU.
// Returns true if successful and sets the output buffer as the active index buffer (R16_UINT).
// pSourceIndices: Xbox 16-bit index data, or nullptr for non-indexed (sequential) draws.
// sourceVertexCount: number of source vertices/indices.
// outputIndexCount: number of output triangle-list indices.
// conversionMode: CXBX_INDEX_CONVERT_QUAD_CW, CXBX_INDEX_CONVERT_QUAD_CCW, or CXBX_INDEX_CONVERT_FAN.
bool CxbxD3D11ConvertIndexBufferGPU(
	const INDEX16* pSourceIndices,
	UINT sourceVertexCount,
	UINT outputIndexCount,
	int conversionMode);

// GPU-accelerated palette texture expansion via compute shader.
// Combines unswizzle + palette lookup in a single CS dispatch.
// Texture must be DXGI_FORMAT_R8G8B8A8_UNORM with DEFAULT+UAV.
// pSwizzledP8Src: raw swizzled P8 pixel data (1 byte per pixel).
// pPaletteData: 256-entry ARGB palette (1024 bytes).
bool CxbxD3D11ExpandPaletteTexture(
	ID3D11Texture2D* pTexture,
	const void* pSwizzledP8Src,
	UINT width,
	UINT height,
	const void* pPaletteData);

// Format conversion type constants (must match CS shader)
#define CXBX_FMTCONV_L6V5U5      1
#define CXBX_FMTCONV_R6G5B5      2  // Same bit layout as L6V5U5, unsigned interpretation
#define CXBX_FMTCONV_V8U8        3  // V8U8/G8B8 are aliases (Xbox fmt 0x28); signedness via colorsign
#define CXBX_FMTCONV_R5G5B5A1    4
#define CXBX_FMTCONV_R4G4B4A4    5
#define CXBX_FMTCONV_A8          6

// Map Xbox texture format to CS format conversion type.
// Returns 0 if the format is not handled by the CS.
inline UINT CxbxGetFormatConvertType(xbox::X_D3DFORMAT X_Format) {
	switch (X_Format) {
	case xbox::X_D3DFMT_L6V5U5:     case xbox::X_D3DFMT_LIN_L6V5U5:   return CXBX_FMTCONV_L6V5U5;
	case xbox::X_D3DFMT_V8U8: /* = X_D3DFMT_G8B8 */                   return CXBX_FMTCONV_V8U8;
	case xbox::X_D3DFMT_R5G5B5A1:   case xbox::X_D3DFMT_LIN_R5G5B5A1: return CXBX_FMTCONV_R5G5B5A1;
	case xbox::X_D3DFMT_R4G4B4A4:   case xbox::X_D3DFMT_LIN_R4G4B4A4: return CXBX_FMTCONV_R4G4B4A4;
	case xbox::X_D3DFMT_A8:         case xbox::X_D3DFMT_LIN_A8:        return CXBX_FMTCONV_A8;
	// YUY2/UYVY are 4:2:2 packed — chroma shared across texel pairs, can't be decoded per-texel.
	// They need a separate dispatch path or fall back to CPU conversion.
	default: return 0;
	}
}

// GPU-accelerated texture format conversion via compute shader.
// Combines optional unswizzle + format decode → R8G8B8A8_UNORM output.
// Texture must be DXGI_FORMAT_R8G8B8A8_UNORM with DEFAULT+UAV.
// fmtType: one of CXBX_FMTCONV_* constants.
// srcRowPitch: row pitch of source data (only used when bSwizzled=false).
bool CxbxD3D11FormatConvertTexture(
	ID3D11Texture2D* pTexture,
	const void* pSrc,
	UINT width,
	UINT height,
	UINT bpp,
	UINT fmtType,
	bool bSwizzled,
	UINT srcRowPitch);

// Vertex format conversion type constants
#define CXBX_VTXCONV_COPY       0
#define CXBX_VTXCONV_NORMSHORT3 1
#define CXBX_VTXCONV_NORMPACKED3 2
#define CXBX_VTXCONV_SHORT3     3
#define CXBX_VTXCONV_PBYTE3     4
#define CXBX_VTXCONV_FLOAT2H    5
#define CXBX_VTXCONV_D3DCOLOR   6
#define CXBX_VTXCONV_NONE       7

// GPU-accelerated vertex format conversion via compute shader.
// Converts Xbox vertex data to D3D11-compatible layouts on the GPU.
// pElementDescriptors: array of numElements * 4 UINTs (srcOffset, dstOffset, convType, copyDwords).
// Returns true if successful, with *ppOutputVB set to the converted vertex buffer.
bool CxbxD3D11ConvertVertexBufferGPU(
	const uint8_t* pSrcVertexData,
	UINT srcDataSize,
	UINT vertexCount,
	UINT srcStride,
	UINT dstStride,
	UINT numElements,
	const UINT* pElementDescriptors,
	UINT dstBufferSize,
	ID3D11Buffer** ppOutputVB);

// Filter D3D11 input layout elements to only include semantics present in
// the shader's input signature (parsed from the DXBC ISGN chunk).

// ******************************************************************
// * vertex fetch — manual vertex fetch from ByteAddressBuffer via SV_VertexID
// ******************************************************************
typedef struct _CxbxDrawContext CxbxDrawContext; // forward decl
struct PGRAPHState; // forward decl for inline buffer draw
void CxbxD3D11VertexFetchInit();
void CxbxD3D11VertexFetchRelease();
void CxbxD3D11VertexFetchDraw(CxbxDrawContext& DrawContext);
void CxbxD3D11DrawInlineBuffer(PGRAPHState* pg);
void CxbxD3D11VertexFetchInvalidateLayout();  // Bump layout CB generation counter
extern bool g_bD3D11VertexFetchDefaultsDirty; // Set true when vertex defaults change

// Create the vertex defaults buffer and bind it to the defaults slot.
// Called once during device initialization.
void CxbxD3D11CreateVertexDefaultsBuffer();

// Lazily create the D3D11 input layout for a vertex declaration (if not yet
// created) and bind it via IASetInputLayout.  Encapsulates all device access
// so callers outside the Rendering folder never touch g_pD3DDevice directly.
void CxbxD3D11SetVertexDeclaration(CxbxVertexDeclaration* pCxbxVertexDeclaration);

// ******************************************************************
// * Device creation data (uses DXGI/D3D11 types)
// ******************************************************************

// information passed to the create device proxy thread
struct EmuD3D8CreateDeviceProxyData
{
	IDXGIAdapter*  Adapter;
	D3D_DRIVER_TYPE   DeviceType;
	struct {
		UINT BackBufferWidth;
		UINT BackBufferHeight;
		UINT FullScreen_RefreshRateInHz;
		BOOL Windowed;
	} HostPresentationParameters;
};

extern EmuD3D8CreateDeviceProxyData  g_EmuCDPD;

// ******************************************************************
// * Rendering helpers (implemented in Backend_D3D11*.cpp)
// ******************************************************************
HRESULT CxbxSetRenderTarget(ID3D11Texture2D* pHostRenderTarget, UINT mipSlice = 0, UINT arraySlice = 0);
void    CxbxSetViewport(D3D11_VIEWPORT *pHostViewport);
HRESULT CxbxBltSurface(ID3D11Texture2D* pSrc, const RECT* pSrcRect, ID3D11Texture2D* pDst, const RECT* pDstRect, D3DTEXTUREFILTERTYPE Filter);
void    CxbxSetDepthStencilSurface(ID3D11Texture2D* pHostDepthStencil);
ID3D11Texture2D* CxbxGetCurrentRenderTarget(); // Returns current RT (non-owning pointer)
HRESULT CxbxGetBackBuffer(ID3D11Texture2D** ppBackBuffer); // Returns back buffer (caller owns ref)
HRESULT CxbxSetStreamSource(UINT HostStreamNumber, ID3D11Buffer* pHostVertexBuffer, UINT VertexStride);
HRESULT CxbxSetVertexShader(ID3D11VertexShader* pHostVertexShader);

// ******************************************************************
// * Query helpers
// ******************************************************************
inline void CxbxQueryIssueBegin(ID3D11Query* pQuery)
{
	g_pD3DDeviceContext->Begin(pQuery);
}

inline void CxbxQueryIssueEnd(ID3D11Query* pQuery)
{
	g_pD3DDeviceContext->End(pQuery);
}

inline HRESULT CxbxQueryGetData(ID3D11Query* pQuery, void* pData, DWORD dwSize, DWORD dwGetDataFlags)
{
	return g_pD3DDeviceContext->GetData(pQuery, pData, dwSize, dwGetDataFlags);
}

// ******************************************************************
// * Resource cache types (defined in HostResource.cpp)
// ******************************************************************
constexpr DWORD D3DUSAGE_INVALID = 0xFFFFFFFF;

typedef struct _resource_info_t {
	Microsoft::WRL::ComPtr<ID3D11Resource> pHostResource;
	DXGI_FORMAT HostFormat = EMUFMT_UNKNOWN;
	DWORD HostUsage = D3DUSAGE_INVALID;
	DWORD dwXboxResourceType = 0;
	void* pXboxData = xbox::zeroptr;
	size_t szXboxDataSize = 0;
	uint32_t lastAccessFrame = 0; // Frame counter for LRU eviction
} resource_info_t;

typedef std::unordered_map<resource_key_t, resource_info_t, resource_key_hash> resource_cache_t;

// HostDevice.cpp
void UpdateDepthStencilFlags(ID3D11Texture2D *pDepthStencilSurface);

// HostResource.cpp
void SetHostResource(xbox::X_D3DResource* pXboxResource, ID3D11Resource* pHostResource, int iTextureStage = -1, DWORD D3DUsage = 0xFFFFFFFFu, DXGI_FORMAT PCFormat = EMUFMT_UNKNOWN);
void FreeHostResource(resource_key_t key);
ID3D11Texture2D* GetHostSurface(xbox::X_D3DResource* pXboxResource, DWORD D3DUsage = 0);
ID3D11Resource* GetHostBaseTexture(xbox::X_D3DResource* pXboxResource, DWORD D3DUsage = 0, int iTextureStage = 0);
ID3D11Resource* GetHostBaseTextureWithFormat(xbox::X_D3DResource* pXboxResource, DWORD D3DUsage, int iTextureStage, DXGI_FORMAT* pOutHostFormat);
ID3D11Texture3D* GetHostVolumeTexture(xbox::X_D3DResource* pXboxResource, int iTextureStage = 0);
resource_cache_t& GetResourceCache(resource_key_t& key);

// Inline resource setters
inline void SetHostSurface(xbox::X_D3DResource* pXboxResource, ID3D11Texture2D* pHostSurface, int iTextureStage = -1)
{
	SetHostResource(pXboxResource, (ID3D11Resource*)pHostSurface, iTextureStage);
}

inline void SetHostVolume(xbox::X_D3DResource* pXboxResource, ID3D11Texture3D* pHostVolume, int iTextureStage = -1)
{
	SetHostResource(pXboxResource, (ID3D11Resource*)pHostVolume, iTextureStage);
}

inline void SetHostTexture(xbox::X_D3DResource* pXboxResource, ID3D11Texture2D* pHostTexture, int iTextureStage = -1)
{
	SetHostResource(pXboxResource, (ID3D11Resource*)pHostTexture, iTextureStage);
}

inline void SetHostVolumeTexture(xbox::X_D3DResource* pXboxResource, ID3D11Texture3D* pHostVolumeTexture, int iTextureStage = -1)
{
	SetHostResource(pXboxResource, (ID3D11Resource*)pHostVolumeTexture, iTextureStage);
}

inline void SetHostCubeTexture(xbox::X_D3DResource* pXboxResource, ID3D11Texture2D* pHostCubeTexture, int iTextureStage = -1)
{
	SetHostResource(pXboxResource, (ID3D11Resource*)pHostCubeTexture, iTextureStage);
}

// HostRender.cpp
bool GetHostRenderTargetDimensions(DWORD *pHostWidth, DWORD *pHostHeight, ID3D11Texture2D* pHostRenderTarget = nullptr);

// HostSync.cpp
void CxbxD3D11InvalidateCachedSRVForTexture(ID3D11Resource* pTexture);

#endif // BACKEND_D3D11_H
