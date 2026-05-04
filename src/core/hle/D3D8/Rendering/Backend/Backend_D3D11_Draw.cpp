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

// ******************************************************************
// * Rendering helpers (D3D11 implementations)
// ******************************************************************

void CxbxD3DClear(DWORD Count, CONST D3DRECT* pRects, DWORD Flags, D3DCOLOR Color, float Z, DWORD Stencil)
{
	LOG_INIT;
	FLOAT clearColor[4];
	clearColor[0] = ((Color >> 16) & 0xFF) / 255.0f;
	clearColor[1] = ((Color >>  8) & 0xFF) / 255.0f;
	clearColor[2] = ((Color >>  0) & 0xFF) / 255.0f;
	clearColor[3] = ((Color >> 24) & 0xFF) / 255.0f;

	// Diagnostic: log clear color (first 3 occurrences only)
	{
		static int s_clearDiag = 0;
		if (s_clearDiag < 3) {
			s_clearDiag++;
			EmuLog(LOG_LEVEL::INFO, "Clear diag [%d]: Color=0x%08X -> RGBA(%f,%f,%f,%f) Flags=0x%X",
				s_clearDiag, Color, clearColor[0], clearColor[1], clearColor[2], clearColor[3], Flags);
		}
	}

	if ((Flags & D3DCLEAR_TARGET) && g_pD3DCurrentRTV != nullptr) {
		if (Count > 0 && pRects != nullptr) {
			ComPtr<ID3D11DeviceContext1> context1;
			if (SUCCEEDED(g_pD3DDeviceContext->QueryInterface(IID_PPV_ARGS(context1.GetAddressOf())))) {
				context1->ClearView(g_pD3DCurrentRTV, clearColor, (const D3D11_RECT*)pRects, Count);
			}
		} else {
			g_pD3DDeviceContext->ClearRenderTargetView(g_pD3DCurrentRTV, clearColor);
		}
	}

	if (g_pD3DDepthStencilView != nullptr) {
		UINT clearFlags = 0;
		if (Flags & D3DCLEAR_ZBUFFER) clearFlags |= D3D11_CLEAR_DEPTH;
		if (Flags & D3DCLEAR_STENCIL) clearFlags |= D3D11_CLEAR_STENCIL;
		if (clearFlags != 0) {
			g_pD3DDeviceContext->ClearDepthStencilView(g_pD3DDepthStencilView, clearFlags, Z, (UINT8)Stencil);
		}
	}
}

HRESULT CxbxBltSurface(ID3D11Texture2D* pSrc, const RECT* pSrcRect, ID3D11Texture2D* pDst, const RECT* pDstRect, D3DTEXTUREFILTERTYPE Filter)
{
	return CxbxD3D11Blt(pSrc, pSrcRect, pDst, pDstRect, Filter);
}

ID3D11Buffer *CxbxDynBuffer::Update(const void *pData, UINT size)
{
	if (size == 0)
		return nullptr;

	// Grow the buffer if it's too small (or doesn't exist yet)
	if (pBuffer == nullptr || capacity < size) {
		Release();
		// Round up to next power-of-two-ish size to avoid frequent re-allocs
		UINT newCap = (size < 4096) ? 4096 : size;
		// Round up to next multiple of 4096
		newCap = (newCap + 4095) & ~4095u;

		D3D11_BUFFER_DESC desc = {};
		desc.ByteWidth = newCap;
		desc.Usage = D3D11_USAGE_DYNAMIC;
		desc.BindFlags = bindFlags;
		desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

		HRESULT hr = g_pD3DDevice->CreateBuffer(&desc, nullptr, &pBuffer);
		if (FAILED(hr) || pBuffer == nullptr)
			return nullptr;
		capacity = newCap;
	}

	// Map-discard and upload
	if (FAILED(CxbxD3D11UpdateDynamicBuffer(pBuffer, pData, size)))
		return nullptr;
	return pBuffer;
}

void CxbxDynBuffer::Release()
{
	if (pBuffer != nullptr) {
		pBuffer->Release();
		pBuffer = nullptr;
	}
	capacity = 0;
}

void CxbxRawSetPixelShader(ID3D11PixelShader* pPixelShader)
{
	g_pD3DDeviceContext->PSSetShader(pPixelShader, nullptr, 0);
}

// ******************************************************************
// * Dual-backend wrappers — D3D11 implementations
// ******************************************************************

static ID3D11VertexShader* s_LastBoundVS = nullptr;

HRESULT CxbxSetVertexShader(ID3D11VertexShader* pHostVertexShader)
{
	if (pHostVertexShader != s_LastBoundVS) {
		g_pD3DDeviceContext->VSSetShader(pHostVertexShader, nullptr, 0);
		s_LastBoundVS = pHostVertexShader;
	}
	return S_OK;
}

void CxbxInvalidateVertexShaderCache()
{
	s_LastBoundVS = nullptr;
	// Also invalidate topology — blit sets TRIANGLELIST directly
	extern void CxbxInvalidateTopologyCache();
	CxbxInvalidateTopologyCache();
	// Also invalidate GS — blit sets GS=nullptr directly
	extern void CxbxInvalidateGSCache();
	CxbxInvalidateGSCache();
}
