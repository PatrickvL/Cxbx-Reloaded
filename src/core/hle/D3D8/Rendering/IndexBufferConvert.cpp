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
#include "IndexBufferConvert.h"
#include "devices/video/nv2a.h" // For NV2ADevice, g_NV2A
#include "devices/video/nv2a_int.h" // For PGRAPHState
#include "devices/video/nv2a_regs.h" // For NV_PGRAPH_SETUPRASTER

// ******************************************************************
// * Triangle fan to triangle list conversion
// ******************************************************************

// Convert a triangle fan (N vertices) to triangle list indices.
// Fan vertex 0 is the hub; for each i in [1..N-2], triangle = (0, i, i+1).
// When pFanIndexData != null, reads original fan indices; otherwise generates sequential.
// Output: (N-2)*3 triangle indices written to pTriangleIndexData.
void CxbxConvertTriFanToTriangleListIndices(
	INDEX16* pFanIndexData,
	unsigned uNrOfFanVertices,
	INDEX16* pTriangleIndexData)
{
	assert(uNrOfFanVertices >= 3);
	assert(pTriangleIndexData);

	unsigned out = 0;
	for (unsigned i = 1; i + 1 < uNrOfFanVertices; i++) {
		pTriangleIndexData[out++] = pFanIndexData ? pFanIndexData[0]     : 0;
		pTriangleIndexData[out++] = pFanIndexData ? pFanIndexData[i]     : (INDEX16)i;
		pTriangleIndexData[out++] = pFanIndexData ? pFanIndexData[i + 1] : (INDEX16)(i + 1);
	}
}

UINT FanToTriangleVertexCount(UINT NrOfFanVertices)
{
	return (NrOfFanVertices >= 3) ? (NrOfFanVertices - 2) * VERTICES_PER_TRIANGLE : 0;
}

INDEX16* CxbxCreateTriFanToTriangleListIndexData(INDEX16* pFanIndexData, unsigned FanVertexCount)
{
	UINT NrOfTriangleIndices = FanToTriangleVertexCount(FanVertexCount);
	INDEX16* pBuffer = (INDEX16*)malloc(NrOfTriangleIndices * sizeof(INDEX16));
	CxbxConvertTriFanToTriangleListIndices(pFanIndexData, FanVertexCount, pBuffer);
	return pBuffer;
}

// ******************************************************************
// * Quad list to triangle list conversion
// ******************************************************************

// Determine winding order from NV2A PGRAPH SETUPRASTER register.
// FRONTFACE bit 23: 0 = CW, 1 = CCW
bool CxbxGetClockWiseWindingOrder()
{
	extern NV2ADevice* g_NV2A;
	if (g_NV2A) {
		PGRAPHState* pg = &g_NV2A->GetDeviceState()->pgraph;
		uint32_t setupraster = pg->regs[NV_PGRAPH_SETUPRASTER / 4];
		return (setupraster & NV_PGRAPH_SETUPRASTER_FRONTFACE) == 0; // 0 = CW
	}
	return XboxRenderStates.GetXboxRenderState(xbox::X_D3DRS_FRONTFACE) == 0x900;
}

UINT QuadToTriangleVertexCount(UINT NrOfQuadVertices)
{
	return (NrOfQuadVertices * VERTICES_PER_TRIANGLE * TRIANGLES_PER_QUAD) / VERTICES_PER_QUAD;
}

// This function convertes quad to triangle indices.
// When pXboxQuadIndexData is set, original quad indices are read from this buffer
// (this use-case is for when an indexed quad draw is to be emulated).
// When pXboxQuadIndexData is null, quad-emulating indices are generated
// (this use-case is for when a non-indexed quad draw is to be emulated).
// The number of indices to generate is specified through uNrOfTriangleIndices.
// Resulting triangle indices are written to pTriangleIndexData, which must
// be pre-allocated to fit the output data.
// (Note, this function is marked 'constexpr' to allow the compiler to optimize
// the case when pXboxQuadIndexData is null)
void CxbxConvertQuadListToTriangleListIndices(
	INDEX16* pXboxQuadIndexData,
	unsigned uNrOfTriangleIndices,
	INDEX16* pTriangleIndexData)
{
	assert(uNrOfTriangleIndices > 0);
	assert(pTriangleIndexData);

	unsigned i = 0;
	unsigned j = 0;
	while (i + (VERTICES_PER_TRIANGLE * TRIANGLES_PER_QUAD) <= uNrOfTriangleIndices) {
		if (CxbxGetClockWiseWindingOrder()) {
			// ABCD becomes ABD+BCD (split along the B-D diagonal, matching NV2A hardware)
			pTriangleIndexData[i + 0] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 0] : j + 0; // A
			pTriangleIndexData[i + 1] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 1] : j + 1; // B
			pTriangleIndexData[i + 2] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 3] : j + 3; // D
			i += VERTICES_PER_TRIANGLE;
			pTriangleIndexData[i + 0] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 1] : j + 1; // B
			pTriangleIndexData[i + 1] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 2] : j + 2; // C
			pTriangleIndexData[i + 2] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 3] : j + 3; // D
			i += VERTICES_PER_TRIANGLE;
		} else {
			// ABCD becomes ADB+BDC (split along the B-D diagonal, matching NV2A hardware)
			pTriangleIndexData[i + 0] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 0] : j + 0; // A
			pTriangleIndexData[i + 1] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 3] : j + 3; // D
			pTriangleIndexData[i + 2] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 1] : j + 1; // B
			i += VERTICES_PER_TRIANGLE;
			pTriangleIndexData[i + 0] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 1] : j + 1; // B
			pTriangleIndexData[i + 1] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 3] : j + 3; // D
			pTriangleIndexData[i + 2] = pXboxQuadIndexData ? pXboxQuadIndexData[j + 2] : j + 2; // C
			i += VERTICES_PER_TRIANGLE;
		}

		// Next quad, please :
		j += VERTICES_PER_QUAD;
	}
}

// Called from EMUPATCH(D3DDevice_DrawIndexedVerticesUP) when PrimitiveType == X_D3DPT_QUADLIST.
// This API receives the number of vertices to draw (VertexCount), the index data that references
// vertices and a single stream of vertex data. The number of vertices to draw indicates the number
// of indices that are going to be fetched. The vertex data is referenced up to the highest index
// number present in the index data.
// To emulate drawing indexed quads, g_pD3DDevice->DrawIndexedPrimitiveUP is called on host,
// whereby the quad indices are converted to triangle indices. This implies for every four
// quad indices, we have to generate (two times three is) six triangle indices. (Note, that
// vertex data undergoes it's own Xbox-to-host conversion, independent from these indices.)
INDEX16* CxbxCreateQuadListToTriangleListIndexData(INDEX16* pXboxQuadIndexData, unsigned QuadVertexCount)
{
	UINT NrOfTriangleIndices = QuadToTriangleVertexCount(QuadVertexCount);
	INDEX16* pQuadToTriangleIndexBuffer = (INDEX16*)malloc(NrOfTriangleIndices * sizeof(INDEX16));
	CxbxConvertQuadListToTriangleListIndices(pXboxQuadIndexData, NrOfTriangleIndices, pQuadToTriangleIndexBuffer);
	return pQuadToTriangleIndexBuffer;
}

void CxbxReleaseQuadListToTriangleListIndexData(void* pHostIndexData)
{
	free(pHostIndexData);
}

// ******************************************************************
// * Line loop closing-line helpers
// ******************************************************************

void CxbxDrawIndexedClosingLineUP(INDEX16 LowIndex, INDEX16 HighIndex, void *pHostVertexStreamZeroData, UINT uiHostVertexStreamZeroStride)
{
	// Close the line by copying the two endpoint vertices and drawing them as a LINELIST.
	uint8_t VertexData[512]; assert(512 >= 2 * uiHostVertexStreamZeroStride);
	uint8_t *FirstVertex = (uint8_t *)pHostVertexStreamZeroData + (LowIndex * uiHostVertexStreamZeroStride);
	uint8_t *SecondVertex = (uint8_t *)pHostVertexStreamZeroData + (HighIndex * uiHostVertexStreamZeroStride);

	memcpy(VertexData, FirstVertex, uiHostVertexStreamZeroStride);
	memcpy(VertexData + uiHostVertexStreamZeroStride, SecondVertex, uiHostVertexStreamZeroStride);

	static CxbxDynBuffer s_ClosingLineVB = { nullptr, 0, D3D11_BIND_VERTEX_BUFFER };
	UINT dataSize = 2 * uiHostVertexStreamZeroStride;
	ID3D11Buffer* pVB = s_ClosingLineVB.Update(VertexData, dataSize);
	if (pVB != nullptr) {
		UINT stride = uiHostVertexStreamZeroStride;
		UINT offset = 0;
		g_pD3DDeviceContext->IASetVertexBuffers(0, 1, &pVB, &stride, &offset);
		CxbxBindThickLineGS(NV097_SET_BEGIN_END_OP_LINES);
		g_pD3DDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
		g_pD3DDeviceContext->Draw(2, 0);
		CxbxUnbindThickLineGS(NV097_SET_BEGIN_END_OP_LINES);
	}

	g_dwPrimPerFrame++;
}
