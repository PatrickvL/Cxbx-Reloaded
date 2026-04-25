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
#include "../IndexBufferConvert.h"
#include "../Backend\Backend_D3D11.h"

// ******************************************************************
// * patch: D3DDevice_Begin
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_Begin)
(
   	X_D3DPRIMITIVETYPE     PrimitiveType
)
{
	LOG_FUNC_ONE_ARG(PrimitiveType);

	CxbxImpl_Begin(PrimitiveType);
}

// ******************************************************************
// * patch: D3DDevice_SetVertexData2f
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexData2f)
(
   	int_xt     Register,
   	float_xt   a,
   	float_xt   b
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(a)
		LOG_FUNC_ARG(b)
		LOG_FUNC_END;

	CxbxImpl_SetVertexData4f(Register, a, b, 0.0f, 1.0f);
}

static inline DWORD FtoDW(FLOAT f) { return *((DWORD*)&f); }
static inline FLOAT DWtoF(DWORD f) { return *((FLOAT*)&f); }

// ******************************************************************
// * patch: D3DDevice_SetVertexData2s
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexData2s)
(
   	int_xt   Register,
   	short_xt a,
   	short_xt b
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(a)
		LOG_FUNC_ARG(b)
		LOG_FUNC_END;

	// Test case: Halo
	// Note : XQEMU verified that the int16_t arguments
	// must be mapped to floats in the range [-32768.0, 32767.0]
	// (See https://github.com/xqemu/xqemu/pull/176)
	const float fa = static_cast<float>(a);
	const float fb = static_cast<float>(b);


	CxbxImpl_SetVertexData4f(Register, a, b, 0.0f, 1.0f);
}

// Overload for logging
static void D3DDevice_SetVertexData4f_16__LTCG_edi1
(
   	xbox::int_xt     Register,
   	xbox::float_xt   a,
   	xbox::float_xt   b,
   	xbox::float_xt   c,
   	xbox::float_xt   d
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(a)
		LOG_FUNC_ARG(b)
		LOG_FUNC_ARG(c)
		LOG_FUNC_ARG(d)
		LOG_FUNC_END;
}

// ******************************************************************
// * patch: D3DDevice_SetVertexData4f_16__LTCG_edi1
// ******************************************************************
// This is an LTCG specific version of SetVertexData4f where the first param is passed in EDI
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexData4f_16__LTCG_edi1)
(
	float_xt   a,
	float_xt   b,
	float_xt   c,
	float_xt   d
)
{
	int_xt Register;

	__asm {
		LTCG_PROLOGUE
		mov  Register, edi
	}

	// Log
	D3DDevice_SetVertexData4f_16__LTCG_edi1(Register, a, b, c, d);

	CxbxImpl_SetVertexData4f(Register, a, b, c, d);

	_asm {
		LTCG_EPILOGUE
		ret  10h
	}
}

// ******************************************************************
// * patch: D3DDevice_SetVertexData4f
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexData4f)
(
   	int_xt     Register,
   	float_xt   a,
   	float_xt   b,
   	float_xt   c,
   	float_xt   d
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(a)
		LOG_FUNC_ARG(b)
		LOG_FUNC_ARG(c)
		LOG_FUNC_ARG(d)
		LOG_FUNC_END;

	CxbxImpl_SetVertexData4f(Register, a, b, c, d);
}

// ******************************************************************
// * patch: D3DDevice_SetVertexData4ub
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexData4ub)
(
	int_xt	Register,
	byte_xt	a,
	byte_xt	b,
	byte_xt	c,
	byte_xt	d
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(a)
		LOG_FUNC_ARG(b)
		LOG_FUNC_ARG(c)
		LOG_FUNC_ARG(d)
		LOG_FUNC_END;

	const float fa = a / 255.0f;
	const float fb = b / 255.0f;
	const float fc = c / 255.0f;
	const float fd = d / 255.0f;

   	CxbxImpl_SetVertexData4f(Register, fa, fb, fc, fd);
}

// ******************************************************************
// * patch: D3DDevice_SetVertexData4s
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexData4s)
(
	int_xt	 Register,
	short_xt a,
	short_xt b,
	short_xt c,
	short_xt d
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(a)
		LOG_FUNC_ARG(b)
		LOG_FUNC_ARG(c)
		LOG_FUNC_ARG(d)
		LOG_FUNC_END;

	// Test case: Halo
	// Note : XQEMU verified that the int16_t arguments
	// must be mapped to floats in the range [-32768.0, 32767.0]
	// (See https://github.com/xqemu/xqemu/pull/176)
	const float fa = static_cast<float>(a);
	const float fb = static_cast<float>(b);
	const float fc = static_cast<float>(c);
	const float fd = static_cast<float>(d);

   	CxbxImpl_SetVertexData4f(Register, fa, fb, fc, fd);
}

// ******************************************************************
// * patch: D3DDevice_SetVertexDataColor
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_SetVertexDataColor)
(
   	int_xt      Register,
   	X_D3DCOLOR  Color
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Register)
		LOG_FUNC_ARG(Color)
		LOG_FUNC_END;

   	const D3DXCOLOR XColor = Color;

   	CxbxImpl_SetVertexData4f(Register, XColor.r, XColor.g, XColor.b, XColor.a);
}

// ******************************************************************
// * patch: D3DDevice_End
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_End)()
{
	LOG_FUNC();

	CxbxImpl_End();
}

// ******************************************************************
// * patch: D3DDevice_BeginPushBuffer
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_BeginPushBuffer)
(
	X_D3DPushBuffer *pPushBuffer
)
{
	LOG_FUNC_ONE_ARG(pPushBuffer);

	// Call through to original Xbox code to redirect the push buffer pointer.
	// Only set recording flag if the trampoline was available and the original
	// code actually ran.
	if (XB_TRMP(D3DDevice_BeginPushBuffer) != nullptr) {
		XB_TRMP(D3DDevice_BeginPushBuffer)(pPushBuffer);
		g_bRecordingPushBuffer = true;
	} else {
		LOG_TEST_CASE("BeginPushBuffer trampoline not available");
	}
}

// ******************************************************************
// * patch: D3DDevice_BeginPushBuffer_0__LTCG_edi1
// ******************************************************************
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_BeginPushBuffer_0__LTCG_edi1)()
{
	X_D3DPushBuffer* pPushBuffer;
	__asm {
		LTCG_PROLOGUE
		mov  pPushBuffer, edi
	}

	EMUPATCH(D3DDevice_BeginPushBuffer)(pPushBuffer);

	__asm {
		LTCG_EPILOGUE
		ret  0
	}
}

// ******************************************************************
// * patch: D3DDevice_EndPushBuffer
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_EndPushBuffer)()
{
	LOG_FUNC();

	// Always clear recording flag
	g_bRecordingPushBuffer = false;

	// Call through to original Xbox code to finalize the push buffer
	if (XB_TRMP(D3DDevice_EndPushBuffer) != nullptr) {
		return XB_TRMP(D3DDevice_EndPushBuffer)();
	} else {
		LOG_TEST_CASE("EndPushBuffer trampoline not available");
		return (HRESULT)0;
	}
}

// ******************************************************************
// * patch: D3DDevice_RunPushBuffer
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_RunPushBuffer)
(
   	X_D3DPushBuffer       *pPushBuffer,
   	X_D3DFixup            *pFixup
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pPushBuffer)
		LOG_FUNC_ARG(pFixup)
		LOG_FUNC_END;

	EmuExecutePushBuffer(pPushBuffer, pFixup);    
}

// ******************************************************************
// * patch: D3DDevice_RunPushBuffer_4__LTCG_eax2
// ******************************************************************
// Overload for logging
static void D3DDevice_RunPushBuffer_4__LTCG_eax2
(
   	xbox::X_D3DPushBuffer       *pPushBuffer,
	xbox::X_D3DFixup            *pFixup
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(pPushBuffer)
		LOG_FUNC_ARG(pFixup)
		LOG_FUNC_END;
}

// This uses a custom calling convention where parameter is passed in EAX
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_RunPushBuffer_4__LTCG_eax2)
(
   	X_D3DPushBuffer *pPushBuffer
)
{
	X_D3DFixup* pFixup;
	__asm {
		LTCG_PROLOGUE
		mov  pFixup, eax
	}

	// Log
	D3DDevice_RunPushBuffer_4__LTCG_eax2(pPushBuffer, pFixup);

	EmuExecutePushBuffer(pPushBuffer, pFixup);

	__asm {
		LTCG_EPILOGUE
		ret  4
	}
}

// ******************************************************************
// * patch: D3DDevice_DrawVertices
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawVertices)
(
   	X_D3DPRIMITIVETYPE PrimitiveType,
   	uint_xt            StartVertex,
   	uint_xt            VertexCount
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(PrimitiveType)
		LOG_FUNC_ARG(StartVertex)
		LOG_FUNC_ARG(VertexCount)
		LOG_FUNC_END;

	// Dxbx Note : In DrawVertices and DrawIndexedVertices, PrimitiveType may not be D3DPT_POLYGON

	if (!IsValidXboxVertexCount(PrimitiveType, VertexCount)) {
		LOG_TEST_CASE("Invalid VertexCount");
		return;
	}

	// During push buffer recording, call through to the original Xbox code
	// so it writes NV2A commands to the user's push buffer. Skip HLE drawing
	// to avoid rendering during the recording phase.
	if (g_bRecordingPushBuffer && XB_TRMP(D3DDevice_DrawVertices) != nullptr) {
		XB_TRMP(D3DDevice_DrawVertices)(PrimitiveType, StartVertex, VertexCount);
		return;
	}

	// TODO : Call unpatched CDevice_SetStateVB[_8](0);

	CxbxUpdateNativeD3DResources();

	CxbxDrawContext DrawContext = {};

	DrawContext.XboxPrimitiveType = PrimitiveType;
	DrawContext.dwVertexCount = VertexCount;
	DrawContext.dwStartVertex = StartVertex;

	CxbxD3D11IABypassDraw(DrawContext);
	g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(PrimitiveType, VertexCount);

	CxbxHandleXboxCallbacks();
}

// ******************************************************************
// * patch: D3DDevice_DrawVerticesUP
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawVerticesUP)
(
   	X_D3DPRIMITIVETYPE  PrimitiveType,
   	uint_xt             VertexCount,
   	CONST PVOID         pVertexStreamZeroData,
   	uint_xt             VertexStreamZeroStride
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(PrimitiveType)
		LOG_FUNC_ARG(VertexCount)
		LOG_FUNC_ARG(pVertexStreamZeroData)
		LOG_FUNC_ARG(VertexStreamZeroStride)
		LOG_FUNC_END;

	if (!IsValidXboxVertexCount(PrimitiveType, VertexCount)) {
		LOG_TEST_CASE("Invalid VertexCount");
		return;
	}

	// During push buffer recording, call through to the original Xbox code
	if (g_bRecordingPushBuffer && XB_TRMP(D3DDevice_DrawVerticesUP) != nullptr) {
		XB_TRMP(D3DDevice_DrawVerticesUP)(PrimitiveType, VertexCount, pVertexStreamZeroData, VertexStreamZeroStride);
		return;
	}

	// TODO : Call unpatched CDevice_SetStateUP();

	CxbxUpdateNativeD3DResources();

	CxbxDrawContext DrawContext = {};

	DrawContext.XboxPrimitiveType = PrimitiveType;
	DrawContext.dwVertexCount = VertexCount;
	DrawContext.pXboxVertexStreamZeroData = pVertexStreamZeroData;
	DrawContext.uiXboxVertexStreamZeroStride = VertexStreamZeroStride;

	CxbxDrawPrimitiveUP(DrawContext);

	CxbxHandleXboxCallbacks();
}

// LTCG specific D3DDevice_DrawVerticesUP function...
// This uses a custom calling convention where pVertexStreamZeroData is passed in EBX
// Test-case: NASCAR Heat 20002
__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawVerticesUP_12__LTCG_ebx3)
(
   	X_D3DPRIMITIVETYPE  PrimitiveType,
   	uint_xt             VertexCount,
   	uint_xt             VertexStreamZeroStride
)
{
   	PVOID pVertexStreamZeroData;
   	__asm {
   	   	LTCG_PROLOGUE
   	   	mov  pVertexStreamZeroData, ebx
   	}

   	EMUPATCH(D3DDevice_DrawVerticesUP)(PrimitiveType, VertexCount, pVertexStreamZeroData, VertexStreamZeroStride);

   	__asm {
   	   	LTCG_EPILOGUE
   	   	ret  0Ch
   	}
}

// ******************************************************************
// * patch: D3DDevice_DrawIndexedVertices
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawIndexedVertices)
(
   	X_D3DPRIMITIVETYPE  PrimitiveType,
   	uint_xt             VertexCount,
   	CONST PWORD         pIndexData
)
{
	// Test-cases : XDK samples (Cartoon, Gamepad)
	// Note : In gamepad.xbe, the gamepad is drawn by D3DDevice_DrawIndexedVertices
	// Dxbx Note : In DrawVertices and DrawIndexedVertices, PrimitiveType may not be D3DPT_POLYGON

	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(PrimitiveType)
		LOG_FUNC_ARG(VertexCount)
		LOG_FUNC_ARG(pIndexData)
		LOG_FUNC_END;

	if (!IsValidXboxVertexCount(PrimitiveType, VertexCount)) {
		LOG_TEST_CASE("Invalid VertexCount");
		return;
	}

	// During push buffer recording, call through to the original Xbox code
	if (g_bRecordingPushBuffer && XB_TRMP(D3DDevice_DrawIndexedVertices) != nullptr) {
		XB_TRMP(D3DDevice_DrawIndexedVertices)(PrimitiveType, VertexCount, pIndexData);
		return;
	}

	// TODO : Call unpatched CDevice_SetStateVB[_8](g_Xbox_BaseVertexIndex);

	CxbxUpdateNativeD3DResources();

	CxbxDrawContext DrawContext = {};

	DrawContext.XboxPrimitiveType = PrimitiveType;
	DrawContext.dwVertexCount = VertexCount;
	DrawContext.dwBaseVertexIndex = g_Xbox_BaseVertexIndex; // Multiplied by vertex stride and added to the vertex buffer start
	DrawContext.pXboxIndexData = pIndexData; // Used to derive VerticesInBuffer

	// Test case JSRF draws all geometry through this function (only sparks are drawn via another method)
	// using X_D3DPT_TRIANGLELIST and X_D3DPT_TRIANGLESTRIP PrimitiveType
	CxbxDrawIndexed(DrawContext);

	CxbxHandleXboxCallbacks();
}

// ******************************************************************
// * patch: D3DDevice_DrawIndexedVerticesUP
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawIndexedVerticesUP)
(
   	X_D3DPRIMITIVETYPE  PrimitiveType,
   	uint_xt                VertexCount,
   	CONST PVOID         pIndexData,
   	CONST PVOID         pVertexStreamZeroData,
   	uint_xt                VertexStreamZeroStride
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(PrimitiveType)
		LOG_FUNC_ARG(VertexCount)
		LOG_FUNC_ARG(pIndexData)
		LOG_FUNC_ARG(pVertexStreamZeroData)
		LOG_FUNC_ARG(VertexStreamZeroStride)
		LOG_FUNC_END;

	if (!IsValidXboxVertexCount(PrimitiveType, VertexCount)) {
		LOG_TEST_CASE("Invalid VertexCount");
		return;
	}

	// TODO : Call unpatched CDevice_SetStateUP();

	CxbxUpdateNativeD3DResources();

		CxbxDrawContext DrawContext = {};
		INDEX16* pXboxIndexData = (INDEX16*)pIndexData;

		DrawContext.XboxPrimitiveType = PrimitiveType;
		DrawContext.dwVertexCount = VertexCount;
		DrawContext.pXboxIndexData = pXboxIndexData; // Used to derive VerticesInBuffer
		// Note : D3DDevice_DrawIndexedVerticesUP does NOT use g_Xbox_BaseVertexIndex, so keep DrawContext.dwBaseVertexIndex at 0!
		DrawContext.pXboxVertexStreamZeroData = pVertexStreamZeroData;
		DrawContext.uiXboxVertexStreamZeroStride = VertexStreamZeroStride;

		// Determine LowIndex and HighIndex *before* VerticesInBuffer gets derived
		WalkIndexBuffer(DrawContext.LowIndex, DrawContext.HighIndex, pXboxIndexData, VertexCount);

		VertexBufferConverter.Apply(&DrawContext);

		INDEX16* pHostIndexData;
		UINT PrimitiveCount = DrawContext.dwHostPrimitiveCount;

		bool bConvertQuadListToTriangleList = (DrawContext.XboxPrimitiveType == X_D3DPT_QUADLIST);
		bool bConvertTriFanToTriangleList = (DrawContext.XboxPrimitiveType == X_D3DPT_TRIANGLEFAN
			|| DrawContext.XboxPrimitiveType == xbox::X_D3DPT_POLYGON);
		bool bConvertedPrimitive = bConvertQuadListToTriangleList || bConvertTriFanToTriangleList;
		if (bConvertQuadListToTriangleList) {
			LOG_TEST_CASE("X_D3DPT_QUADLIST");
			// Test-case : Buffy: The Vampire Slayer
			// Test-case : XDK samples : FastLoad, BackBufferScale, DisplacementMap, Donuts3D, VolumeLight, PersistDisplay, PolynomialTextureMaps, SwapCallback, Tiling, VolumeFog, DebugKeyboard, Gamepad
			// Convert draw arguments from quads to triangles :
			pHostIndexData = CxbxCreateQuadListToTriangleListIndexData(pXboxIndexData, VertexCount);
			PrimitiveCount *= TRIANGLES_PER_QUAD;
			// Note, that LowIndex and HighIndex won't change due to this quad-to-triangle conversion,
			// so it's less work to WalkIndexBuffer over the input instead of the converted index buffer.
		} else if (bConvertTriFanToTriangleList) {
			pHostIndexData = CxbxCreateTriFanToTriangleListIndexData(pXboxIndexData, VertexCount);
			PrimitiveCount = (VertexCount >= 3) ? VertexCount - 2 : 0;
		} else {
			// LOG_TEST_CASE("DrawIndexedPrimitiveUP"); // Test-case : Burnout, Namco Museum 50th Anniversary
			pHostIndexData = pXboxIndexData;
		}

		HRESULT hRet;
		// D3D11 has no DrawIndexedPrimitiveUP - use reusable dynamic buffers
		UINT vertexDataSize = DrawContext.dwVertexCount * DrawContext.uiHostVertexStreamZeroStride;

		static CxbxDynBuffer s_IdxUpVB = { nullptr, 0, D3D11_BIND_VERTEX_BUFFER };
		ID3D11Buffer* pVB = s_IdxUpVB.Update(DrawContext.pHostVertexStreamZeroData, vertexDataSize);

		if (pVB != nullptr) {
			UINT stride = DrawContext.uiHostVertexStreamZeroStride;
			UINT offset = 0;
			g_pD3DDeviceContext->IASetVertexBuffers(0, 1, &pVB, &stride, &offset);

			UINT indexCount = bConvertedPrimitive ? PrimitiveCount * 3 : DrawContext.dwVertexCount;
			bool bCsHandled = false;

			// Try GPU index conversion for converted primitives
			if (bConvertedPrimitive) {
				int mode = bConvertQuadListToTriangleList ?
					(CxbxGetClockWiseWindingOrder() ? CXBX_INDEX_CONVERT_QUAD_CW : CXBX_INDEX_CONVERT_QUAD_CCW) :
					CXBX_INDEX_CONVERT_FAN;
				bCsHandled = CxbxD3D11ConvertIndexBufferGPU(pXboxIndexData, VertexCount, indexCount, mode);
			}

			if (!bCsHandled) {
				// CPU path: upload converted (or original) index data
				UINT indexDataSize = indexCount * sizeof(INDEX16);
				static CxbxDynBuffer s_IdxUpIB = { nullptr, 0, D3D11_BIND_INDEX_BUFFER };
				ID3D11Buffer* pIB = s_IdxUpIB.Update(pHostIndexData, indexDataSize);
				if (pIB != nullptr) {
					g_pD3DDeviceContext->IASetIndexBuffer(pIB, DXGI_FORMAT_R16_UINT, 0);
				}
			}

			D3D_PRIMITIVE_TOPOLOGY topology = bConvertedPrimitive ?
				D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST :
				EmuXB2PC_D3D11PrimitiveTopology(DrawContext.XboxPrimitiveType);
			g_pD3DDeviceContext->IASetPrimitiveTopology(topology);
			CxbxBindThickLineGS(DrawContext.XboxPrimitiveType);
			g_pD3DDeviceContext->DrawIndexed(indexCount, 0, 0);
			CxbxUnbindThickLineGS(DrawContext.XboxPrimitiveType);
			hRet = S_OK;
		}

		DEBUG_D3DRESULT(hRet, "g_pD3DDevice->DrawIndexedPrimitiveUP");

		if (bConvertQuadListToTriangleList) {
			CxbxReleaseQuadListToTriangleListIndexData(pHostIndexData);
		}
		else if (bConvertTriFanToTriangleList) {
			free(pHostIndexData);
		}

		g_dwPrimPerFrame += PrimitiveCount;
		if (DrawContext.XboxPrimitiveType == X_D3DPT_LINELOOP) {
			// Close line-loops using a final single line, drawn from the end to the start vertex
			LOG_TEST_CASE("X_D3DPT_LINELOOP"); // TODO : Which titles reach this test-case?
			// Read the end and start index from the supplied index data
			INDEX16 LowIndex = pXboxIndexData[0];
			INDEX16 HighIndex = pXboxIndexData[DrawContext.dwHostPrimitiveCount];
			// If needed, swap so highest index is higher than lowest (duh)
			if (HighIndex < LowIndex) {
				std::swap(HighIndex, LowIndex);
			}

			// Close line-loops using a final single line, drawn from the end to the start vertex :
			CxbxDrawIndexedClosingLineUP(
				LowIndex,
				HighIndex,
				DrawContext.pHostVertexStreamZeroData,
				DrawContext.uiHostVertexStreamZeroStride
			);
		}

	CxbxHandleXboxCallbacks();
}

// ******************************************************************
// * patch: CDevice_SetStateVB
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(CDevice_SetStateVB)(ulong_xt Unknown1)
{
	addr_xt _this;
	__asm mov _this, ecx;

	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(_this)
		LOG_FUNC_ARG(Unknown1)
		LOG_FUNC_END;

	// During push buffer recording, call through to the original Xbox code
	// so that NV2A state commands are written to the user's push buffer.
	// This is critical: DrawVertices calls SetStateVB internally, and the
	// original DrawVertices code depends on SetStateVB to manage the push
	// buffer write pointer and space.
	if (g_bRecordingPushBuffer && XB_TRMP(CDevice_SetStateVB) != nullptr) {
		__asm mov ecx, _this;
		XB_TRMP(CDevice_SetStateVB)(Unknown1);
		return;
	}

	LOG_UNIMPLEMENTED();
}

xbox::void_xt WINAPI xbox::EMUPATCH(CDevice_SetStateVB_8)(addr_xt _this, ulong_xt Unknown1)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(_this)
		LOG_FUNC_ARG(Unknown1)
		LOG_FUNC_END;

	if (g_bRecordingPushBuffer && XB_TRMP(CDevice_SetStateVB_8) != nullptr) {
		XB_TRMP(CDevice_SetStateVB_8)(_this, Unknown1);
		return;
	}

	LOG_UNIMPLEMENTED();
}

// ******************************************************************
// * patch: CDevice_SetStateUP (D3D::CDevice::SetStateUP)
// ******************************************************************
xbox::void_xt CxbxrImpl_CDevice_SetStateUP(xbox::addr_xt _this)
{
	// During push buffer recording, call through to original Xbox code
	if (g_bRecordingPushBuffer && XB_TRMP(CDevice_SetStateUP) != nullptr) {
		__asm mov ecx, _this;
		XB_TRMP(CDevice_SetStateUP)();
		return;
	}

	LOG_UNIMPLEMENTED();

	// TODO: Anything?
	//__asm int 3;
}

xbox::void_xt WINAPI xbox::EMUPATCH(CDevice_SetStateUP)()
{
	addr_xt _this;
	__asm mov _this, ecx;

	LOG_FUNC_ONE_ARG(_this);

	CxbxrImpl_CDevice_SetStateUP(_this);
}

xbox::void_xt WINAPI xbox::EMUPATCH(CDevice_SetStateUP_4)(xbox::addr_xt _this)
{
	LOG_FUNC_ONE_ARG(_this);

	if (g_bRecordingPushBuffer && XB_TRMP(CDevice_SetStateUP_4) != nullptr) {
		XB_TRMP(CDevice_SetStateUP_4)(_this);
		return;
	}

	CxbxrImpl_CDevice_SetStateUP(_this);
}
static void CDevice_SetStateUP_0__LTCG_esi1(xbox::addr_xt _this)
{
	LOG_FUNC_ONE_ARG(_this);
}

__declspec(naked) xbox::void_xt WINAPI xbox::EMUPATCH(CDevice_SetStateUP_0__LTCG_esi1)()
{
	addr_xt _this;
	__asm {
		LTCG_PROLOGUE
		mov  _this, esi
	}

	// Log
	CDevice_SetStateUP_0__LTCG_esi1(_this);

	CxbxrImpl_CDevice_SetStateUP(_this);

	__asm {
		LTCG_EPILOGUE
		ret
	}
}

// ******************************************************************
// * patch: D3DDevice_SetStipple
// ******************************************************************
void WINAPI xbox::EMUPATCH(D3DDevice_SetStipple)( dword_xt* pPattern )
{
	LOG_FUNC_ONE_ARG(pPattern);

	// We need an OpenGL port... badly

	LOG_IGNORED();
}

// ******************************************************************
// * patch: D3DDevice_SetSwapCallback
// ******************************************************************
void WINAPI xbox::EMUPATCH(D3DDevice_SetSwapCallback)
(
	X_D3DSWAPCALLBACK		pCallback
)
{
	LOG_FUNC_ONE_ARG(pCallback);

   	g_pXbox_SwapCallback = pCallback;
}

// ******************************************************************
// * patch: D3DDevice_PrimeVertexCache
// ******************************************************************
xbox::void_xt WINAPI xbox::EMUPATCH(D3DDevice_PrimeVertexCache)
(
	uint_xt  VertexCount,
	WORD *pIndexData
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(VertexCount)
		LOG_FUNC_ARG(pIndexData)
		LOG_FUNC_END;

	// TODO: Implement
	LOG_UNIMPLEMENTED();
}

// ******************************************************************
// * patch: D3DDevice_DrawRectPatch
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawRectPatch)
(
	uint_xt					Handle,
	CONST float_xt				*pNumSegs,
	CONST X_D3DRECTPATCH_INFO *pRectPatchInfo
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Handle)
		LOG_FUNC_ARG(pNumSegs)
		LOG_FUNC_ARG(pRectPatchInfo)
		LOG_FUNC_END;

	CxbxUpdateNativeD3DResources();

	// D3D11 has no DrawRectPatch - use CPU tessellation
	HRESULT hRet = CxbxDrawRectPatchD3D11(Handle, pNumSegs, pRectPatchInfo);
	DEBUG_D3DRESULT(hRet, "CxbxDrawRectPatchD3D11");

	return hRet;
}

// ******************************************************************
// * patch: D3DDevice_DrawTriPatch
// ******************************************************************
xbox::hresult_xt WINAPI xbox::EMUPATCH(D3DDevice_DrawTriPatch)
(
	uint_xt					Handle,
	CONST float_xt				*pNumSegs,
	CONST X_D3DTRIPATCH_INFO* pTriPatchInfo
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(Handle)
		LOG_FUNC_ARG(pNumSegs)
		LOG_FUNC_ARG(pTriPatchInfo)
		LOG_FUNC_END;

	CxbxUpdateNativeD3DResources();

	// D3D11 has no DrawTriPatch - use CPU tessellation
	HRESULT hRet = CxbxDrawTriPatchD3D11(Handle, pNumSegs, pTriPatchInfo);
	DEBUG_D3DRESULT(hRet, "CxbxDrawTriPatchD3D11");

	return hRet;
}
