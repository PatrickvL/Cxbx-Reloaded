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

// Internal shared state for Backend_D3D11 split files.
// Not for use outside the Backend/ folder.
#ifndef BACKEND_D3D11_INTERNAL_H
#define BACKEND_D3D11_INTERNAL_H

#undef LOG_PREFIX
#define LOG_PREFIX CXBXR_MODULE::D3D8

#include "Backend_D3D11.h"
#include "../RenderGlobals.h"
#include "core\kernel\init\CxbxKrnl.h"
#include "core\hle\D3D8\XbD3D8Logging.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\Rendering\Backend\Shading\Shader.h"

#include <cstring>
#include <unordered_map>
#include <utility>
#include <vector>
#include <wrl/client.h>
using namespace Microsoft::WRL;

// ******************************************************************
// * D3D11 state objects (shared across split files)
// ******************************************************************
extern ComPtr<ID3D11RasterizerState>   g_pD3DRasterizerState;
extern ComPtr<ID3D11DepthStencilState> g_pD3DDepthStencilState;
extern ComPtr<ID3D11BlendState>        g_pD3DBlendState;

// ******************************************************************
// * Constant buffer shadow arrays
// ******************************************************************
extern float g_D3D11VSConstants[CXBX_D3D11_VS_CB_COUNT][4];
extern UINT  g_D3D11VSConstantsDirtyMin;
extern UINT  g_D3D11VSConstantsDirtyMax;

// ******************************************************************
// * Blit shader resources
// ******************************************************************
extern ID3D11VertexShader  *g_pD3D11BlitVS;
extern ID3D11PixelShader   *g_pD3D11BlitPS;
extern ID3D11SamplerState  *g_pD3D11BlitSamplerLinear;
extern ID3D11SamplerState  *g_pD3D11BlitSamplerPoint;

// ******************************************************************
// * Geometry shader resources
// ******************************************************************
extern ID3D11GeometryShader *g_pD3D11PointSpriteGS;
extern ID3D11Buffer         *g_pD3D11GSConstantBuffer;
extern bool                  g_bPointSpriteEnabled;
extern ID3D11GeometryShader *g_pD3D11ThickLineGS;
extern float                 g_fLineWidth;

// ******************************************************************
// * Compute shader resources — unswizzle
// ******************************************************************
extern ID3D11ComputeShader       *g_pD3D11UnswizzleCS;
extern ID3D11ComputeShader       *g_pD3D11UnswizzleBGRA_CS; // float4 variant for B8G8R8A8_UNORM UAV
extern ID3D11Buffer              *g_pD3D11UnswizzleCB;
extern ID3D11Buffer              *g_pD3D11UnswizzleStagingBuf;
extern UINT                       g_UnswizzleStagingBufSize;
extern ID3D11ShaderResourceView  *g_pD3D11UnswizzleSRV;

// ******************************************************************
// * Compute shader resources — index convert
// ******************************************************************
extern ID3D11ComputeShader       *g_pD3D11IndexConvertCS;
extern ID3D11Buffer              *g_pD3D11IndexConvertCB;
extern ID3D11Buffer              *g_pD3D11IndexConvertInputBuf;
extern UINT                       g_IndexConvertInputBufSize;
extern ID3D11ShaderResourceView  *g_pD3D11IndexConvertInputSRV;
extern ID3D11Buffer              *g_pD3D11IndexConvertOutputBuf;
extern UINT                       g_IndexConvertOutputBufSize;
extern ID3D11UnorderedAccessView *g_pD3D11IndexConvertOutputUAV;

// ******************************************************************
// * Compute shader resources — palette expand
// ******************************************************************
extern ID3D11ComputeShader       *g_pD3D11PaletteExpandCS;
extern ID3D11Buffer              *g_pD3D11PaletteExpandCB;
extern ID3D11Buffer              *g_pD3D11PaletteBuf;
extern ID3D11ShaderResourceView  *g_pD3D11PaletteSRV;

// ******************************************************************
// * Compute shader resources — format conversion
// ******************************************************************
extern ID3D11ComputeShader       *g_pD3D11FormatConvertCS;
extern ID3D11Buffer              *g_pD3D11FormatConvertCB;

// ******************************************************************
// * Register combiner interpreter (PS ubershader)
// ******************************************************************
extern ID3D11PixelShader         *g_pD3D11RCInterpreterPS;       // RC interpreter ubershader
extern ID3D11Buffer              *g_pD3D11RCInterpreterAuxCB;  // PSAuxCBLayout (software-computed fields)

// RC interpreter constant buffer layout — shared with the HLSL cbuffer
// definition in CxbxRegisterCombinerInterpreterState.hlsli.
#include "../Shaders/CxbxNV2APixelShaderConstants.hlsli"
#include "../Shaders/CxbxRegisterCombinerInterpreterState.hlsli"

// ******************************************************************
// * Vertex shader interpreter (VS ubershader)
// ******************************************************************
extern bool                       g_bUseVSInterpreter;
extern ID3D11VertexShader        *g_pD3D11VSInterpreterVS;
extern ID3DBlob                  *g_pD3D11VSInterpreterBytecode;
extern ID3D11Buffer              *g_pD3D11XFPRBuf;              // XFPR (Transform Program RAM) StructuredBuffer — pg->xf.xfpr[]
extern ID3D11ShaderResourceView  *g_pD3D11XFPRSRV;              // SRV for g_XFPR : register(t5)

// VS interpreter instruction field constants — shared with HLSL.
#include "../Shaders/CxbxVertexShaderInterpreterState.hlsli"

// ******************************************************************
// * Compute shader resources — vertex convert
// ******************************************************************
extern ID3D11ComputeShader       *g_pD3D11VertexConvertCS;
extern ID3D11Buffer              *g_pD3D11VertexConvertCB;
extern ID3D11Buffer              *g_pD3D11VertexConvertSrcBuf;
extern UINT                       g_VertexConvertSrcBufSize;
extern ID3D11ShaderResourceView  *g_pD3D11VertexConvertSrcSRV;

// ******************************************************************
// * Internal helper functions (shared across split files)
// ******************************************************************

// CXBX_INDEX_CONVERT_* and CXBX_VTXCONV_* are defined in Backend_D3D11.h

HRESULT CxbxD3D11CreateConstantBuffer(UINT byteWidth, bool bDynamic, ID3D11Buffer** ppBuffer);
HRESULT CxbxD3D11CreateRawBuffer(UINT byteWidth, bool bDynamic, UINT extraBindFlags, ID3D11Buffer** ppBuffer);
HRESULT CxbxD3D11CreateRawBufferSRV(ID3D11Buffer* pBuffer, UINT byteWidth, ID3D11ShaderResourceView** ppSRV);
HRESULT CxbxD3D11CreateTypedBufferUAV(ID3D11Buffer* pBuffer, UINT byteWidth, ID3D11UnorderedAccessView** ppUAV);
void CxbxD3D11EnsureRawStagingBuffer(UINT requiredSize, ID3D11Buffer** ppBuffer, UINT* pCurrentSize, ID3D11ShaderResourceView** ppSRV, const char* debugName);
void CxbxD3D11DispatchCS(ID3D11ComputeShader* pShader, ID3D11Buffer* pCB, UINT numSRVs, ID3D11ShaderResourceView* const* ppSRVs, ID3D11UnorderedAccessView* pUAV, UINT groupsX, UINT groupsY, UINT groupsZ);
HRESULT CxbxD3D11UpdateDynamicBuffer(ID3D11Buffer* pBuffer, const void* pData, size_t dataSize);

// RTV cache helpers
// Key: (texture*, mipSlice, arraySlice) — arraySlice selects cubemap face (0-5) or 0 for 2D
using RTVCacheKey = std::tuple<ID3D11Texture2D*, UINT, UINT>;
struct RTVCacheKeyHash {
	size_t operator()(const RTVCacheKey& k) const {
		struct { const void* p; UINT a; UINT b; } packed = { std::get<0>(k), std::get<1>(k), std::get<2>(k) };
		return static_cast<size_t>(ComputeHash(&packed, sizeof(packed)));
	}
};
void ClearRTVCache();
extern std::unordered_map<RTVCacheKey, ID3D11RenderTargetView*, RTVCacheKeyHash> g_RTVCache;

// Scene begin/end (defined in HostDevice.cpp or HostRender.cpp)
extern void CxbxBeginScene();
extern void CxbxEndScene();

#endif // BACKEND_D3D11_INTERNAL_H
