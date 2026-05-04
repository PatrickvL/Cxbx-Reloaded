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
// ******************************************************************
// *  NV2A Hardware Tessellation for D3D11
// *
// *  Implements the PGRAPH SET_END_PATCH path: decodes forward-difference
// *  matrix data from strip curves collected by SET_CURVE_DATA, evaluates
// *  the tessellated surface, and renders via CxbxD3D11VertexFetchDraw.
// ******************************************************************

#include "PatchDraw.h"


#define LOG_PREFIX CXBXR_MODULE::D3D8

#include "RenderGlobals.h"
#include "Backend\Backend_D3D11.h" // CxbxD3D11VertexFetchDraw, CxbxUpdateNativeD3DResources
#include "core/hle/D3D8/XbVertexBuffer.h" // For CxbxDrawContext
#include "core/kernel/support/Emu.h"
#include "devices/video/nv2a_int.h" // For PGRAPHState, PatchState

#include <vector>
#include <cmath>
#include <cassert>

using namespace xbox;

// ---------------------------------------------------------------
// Float3 type
// ---------------------------------------------------------------
struct Float3 { float x, y, z; };

// ---------------------------------------------------------------
// Build a Float3 triangle list from a position grid (fast path)
// ---------------------------------------------------------------
static void BuildTriangleListFromGrid(
	const std::vector<Float3> &grid, UINT samplesU, UINT samplesV,
	std::vector<Float3> &outVertices)
{
	outVertices.clear();
	outVertices.reserve((samplesU - 1) * (samplesV - 1) * 6);

	for (UINT y = 0; y < samplesV - 1; y++) {
		for (UINT x = 0; x < samplesU - 1; x++) {
			const Float3 &p00 = grid[y * samplesU + x];
			const Float3 &p10 = grid[y * samplesU + x + 1];
			const Float3 &p01 = grid[(y + 1) * samplesU + x];
			const Float3 &p11 = grid[(y + 1) * samplesU + x + 1];

			outVertices.push_back(p00);
			outVertices.push_back(p10);
			outVertices.push_back(p01);
			outVertices.push_back(p10);
			outVertices.push_back(p11);
			outVertices.push_back(p01);
		}
	}
}

// ---------------------------------------------------------------
// NV2A Hardware Tessellation: D3D11_draw_patch
// Called from PGRAPH when SET_END_PATCH is received.
// Decodes forward-difference matrix data from strip curves and evaluates
// the tessellated surface, rendering via CxbxD3D11VertexFetchDraw.
// ---------------------------------------------------------------

extern void CxbxUpdateNativeD3DResources();

// Curve type constants (NV097_SET_BEGIN_END_CURVE_CMD values)
#define NV2A_CURVE_END_DATA         0
#define NV2A_CURVE_STRIP            1
#define NV2A_CURVE_LEFT_GUARD       2
#define NV2A_CURVE_RIGHT_GUARD      3

// Evaluate NV2A hardware tessellation patch and draw the resulting triangle grid.
//
// Known issue: a tiny gap remains between adjacent patches at tessellation levels >= 2.
// The gap is NOT caused by:
//   - FD accumulation error (exact endpoint computed via binomial sum)
//   - Wrong swatch handling (per-swatch drawing is implemented correctly)
//   - Transition curves (gap is identical with transitions disabled)
//   - Mismatched corner values (verified: adjacent patch corners match exactly)
//   - Guard curves (coefficients don't contain position FD data in the expected layout;
//     their first float4 doesn't match strip col0 values -- likely a different encoding)
// The gap likely requires understanding how the NV2A hardware uses guard curves
// (types 2/3) to provide exact edge values that override FD-stepped boundaries.
void D3D11_draw_patch(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;
	PatchState &patch = pg->patch;

	if (patch.curveCount == 0) {
		return;
	}

	// Decode patch0: per-attribute u-order (4 bits per hw attr, value = order-1, 0=disabled)
	// Find position attribute (first enabled = attr 0 typically)
	int posUOrder = 0;
	int posAttrOffset = 0; // float4 offset within a strip row to reach position coeffs
	{
		for (int a = 0; a < 8; a++) {
			int orderMinus1 = (patch.patch0 >> (a * 4)) & 0xF;
			if (orderMinus1 > 0) {
				if (posUOrder == 0) {
					posUOrder = orderMinus1 + 1;
					break;
				}
			}
		}
	}

	if (posUOrder < 2) {
		return; // No valid position attribute
	}

	// Decode patch3: position v-order, numCoeffs per strip row
	int posVOrder = ((patch.patch3 >> 6) & 0xF) + 1;
	int numCoeffsPerStrip = (patch.patch3 >> 24) & 0xFF;

	// Decode patch2: maxSwatch, partialWidth, partialHeight, nSwatchU/V
	int maxSwatch    = (patch.patch2 >> 16) & 0x1F;
	int partialWidth = (patch.patch2 >> 21) & 0x1F;
	int nSwatchU     = (patch.patch2 >> 8) & 0xFF;

	// The FD coefficients are scaled for a specific number of U-steps.
	// nSwatchU=0 means this is a partial-width swatch; use partialWidth.
	// Otherwise this is a full-width swatch; use maxSwatch.
	int numStepsU = (nSwatchU == 0) ? partialWidth : maxSwatch;
	if (numStepsU < 1) numStepsU = 8;

	// Collect strip curves (curveType == 1)
	std::vector<int> stripIndices;
	for (int c = 0; c < patch.curveCount; c++) {
		if (patch.curves[c].curveType == NV2A_CURVE_STRIP) {
			stripIndices.push_back(c);
		}
	}

	int numRows = (int)stripIndices.size();
	if (numRows < 1) {
		return;
	}

	// Grid dimensions: numOutU columns from FD stepping, numRows rows from strip curves.
	int numOutU = numStepsU + 1;

	// Transition curves extend the grid by one column or row to stitch adjacent
	// patches at different tessellation levels.
	// Type 5 (INNER_TRANSITION): extra column, runs in V direction
	// Type 4 (OUTER_TRANSITION): extra row, runs in U direction
	int transitionCurveU = -1;
	int transitionCurveV = -1;
	for (int c = 0; c < patch.curveCount; c++) {
		if (patch.curves[c].curveType == 5 && transitionCurveU < 0) transitionCurveU = c;
		if (patch.curves[c].curveType == 4 && transitionCurveV < 0) transitionCurveV = c;
	}

	bool hasExtraColU = (transitionCurveU >= 0) && (numOutU < numRows);
	bool hasExtraRowV = (transitionCurveV >= 0) && (numRows < numOutU);
	int finalOutU = numOutU + (hasExtraColU ? 1 : 0);
	int finalOutV = numRows + (hasExtraRowV ? 1 : 0);

	// Evaluate the FD grid from strip curves.
	// Each strip curve provides posUOrder forward-difference coefficients [f, df, d2f, ...]
	// that are stepped numStepsU times to produce numOutU output vertices per row.
	std::vector<Float3> grid(finalOutU * finalOutV);

	for (int row = 0; row < numRows; row++) {
		PatchCurve &curve = patch.curves[stripIndices[row]];
		const float *rowCoeffs = &patch.coefficients[curve.coeffStart * 4];

		float fd[16][3];
		for (int k = 0; k < posUOrder && k < 16; k++) {
			int idx = (posAttrOffset + k) * 4;
			fd[k][0] = rowCoeffs[idx + 0];
			fd[k][1] = rowCoeffs[idx + 1];
			fd[k][2] = rowCoeffs[idx + 2];
		}

		// Compute exact endpoint using binomial sum to avoid FD accumulation error.
		// exact(n) = sum_{k=0}^{order-1} C(n,k) * fd_initial[k]
		float exactEnd[3] = {0, 0, 0};
		{
			double binom = 1.0;
			for (int k = 0; k < posUOrder && k < 16; k++) {
				exactEnd[0] += (float)(binom * fd[k][0]);
				exactEnd[1] += (float)(binom * fd[k][1]);
				exactEnd[2] += (float)(binom * fd[k][2]);
				binom = binom * (double)(numStepsU - k) / (double)(k + 1);
			}
		}

		Float3 *outRow = &grid[row * finalOutU];
		outRow[0] = { fd[0][0], fd[0][1], fd[0][2] };

		for (int step = 1; step < numOutU; step++) {
			for (int i = 0; i < posUOrder - 1; i++) {
				fd[i][0] += fd[i + 1][0];
				fd[i][1] += fd[i + 1][1];
				fd[i][2] += fd[i + 1][2];
			}
			outRow[step] = { fd[0][0], fd[0][1], fd[0][2] };
		}

		// Replace last point with exact value to eliminate FD drift
		outRow[numStepsU] = { exactEnd[0], exactEnd[1], exactEnd[2] };
	}

	// U-transition: extra column (runs in V-direction)
	// The transition curve has same layout as a strip (numCoeffsPerStrip float4s per row-equivalent)
	// but it runs in V direction. Position coefficients are posVOrder entries at posAttrOffset.
	if (hasExtraColU) {
		PatchCurve &tCurve = patch.curves[transitionCurveU];
		const float *tCoeffs = &patch.coefficients[tCurve.coeffStart * 4];

		// For a V-running curve with same layout as strips, position is at posAttrOffset
		int tPosOffset = (tCurve.coeffCount == numCoeffsPerStrip) ? posAttrOffset : 0;
		int tOrder = posVOrder;

		float fd[16][3];
		for (int k = 0; k < tOrder && k < 16; k++) {
			int idx = (tPosOffset + k) * 4;
			fd[k][0] = tCoeffs[idx + 0];
			fd[k][1] = tCoeffs[idx + 1];
			fd[k][2] = tCoeffs[idx + 2];
		}

		grid[0 * finalOutU + numOutU] = { fd[0][0], fd[0][1], fd[0][2] };

		for (int row = 1; row < numRows; row++) {
			for (int i = 0; i < tOrder - 1; i++) {
				fd[i][0] += fd[i + 1][0];
				fd[i][1] += fd[i + 1][1];
				fd[i][2] += fd[i + 1][2];
			}
			grid[row * finalOutU + numOutU] = { fd[0][0], fd[0][1], fd[0][2] };
		}
	}

	// V-transition: extra row (runs in U-direction)
	// Same layout as a strip - position at posAttrOffset with posUOrder entries.
	if (hasExtraRowV) {
		PatchCurve &tCurve = patch.curves[transitionCurveV];
		const float *tCoeffs = &patch.coefficients[tCurve.coeffStart * 4];

		int extraRow = numRows;
		int tPosOffset = (tCurve.coeffCount == numCoeffsPerStrip) ? posAttrOffset : 0;

		float fd[16][3];
		for (int k = 0; k < posUOrder && k < 16; k++) {
			int idx = (tPosOffset + k) * 4;
			fd[k][0] = tCoeffs[idx + 0];
			fd[k][1] = tCoeffs[idx + 1];
			fd[k][2] = tCoeffs[idx + 2];
		}

		grid[extraRow * finalOutU + 0] = { fd[0][0], fd[0][1], fd[0][2] };

		for (int step = 1; step < numOutU; step++) {
			for (int i = 0; i < posUOrder - 1; i++) {
				fd[i][0] += fd[i + 1][0];
				fd[i][1] += fd[i + 1][1];
				fd[i][2] += fd[i + 1][2];
			}
			grid[extraRow * finalOutU + step] = { fd[0][0], fd[0][1], fd[0][2] };
		}

		// Corner point where both transitions meet
		if (hasExtraColU) {
			PatchCurve &uCurve = patch.curves[transitionCurveU];
			const float *uCoeffs = &patch.coefficients[uCurve.coeffStart * 4];
			int cPosOffset = (uCurve.coeffCount == numCoeffsPerStrip) ? posAttrOffset : 0;
			int cOrder = posVOrder;
			float fdC[16][3];
			for (int k = 0; k < cOrder && k < 16; k++) {
				int idx = (cPosOffset + k) * 4;
				fdC[k][0] = uCoeffs[idx + 0];
				fdC[k][1] = uCoeffs[idx + 1];
				fdC[k][2] = uCoeffs[idx + 2];
			}
			for (int s = 0; s < numRows; s++) {
				for (int i = 0; i < cOrder - 1; i++) {
					fdC[i][0] += fdC[i + 1][0];
					fdC[i][1] += fdC[i + 1][1];
					fdC[i][2] += fdC[i + 1][2];
				}
			}
			grid[extraRow * finalOutU + numOutU] = { fdC[0][0], fdC[0][1], fdC[0][2] };
		}
	}

	// Build triangle list from the grid (finalOutV x finalOutU)
	UINT samplesU = (UINT)finalOutU;
	UINT samplesV = (UINT)finalOutV;

	std::vector<Float3> triList;
	BuildTriangleListFromGrid(grid, samplesU, samplesV, triList);

	if (triList.empty()) {
		return;
	}

	// Route through VertexFetch as a UP (user-pointer) draw.
	VertexAttribute savedAttr0 = pg->vertex_attributes[0];

	pg->vertex_attributes[0].format = 2; // NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_F
	pg->vertex_attributes[0].count = 3;
	pg->vertex_attributes[0].stride = sizeof(Float3);
	pg->vertex_attributes[0].offset = 0;

	CxbxUpdateNativeD3DResources();

	CxbxDrawContext DrawContext = {};
	DrawContext.XboxPrimitiveType = X_D3DPT_TRIANGLELIST;
	DrawContext.dwVertexCount = (DWORD)triList.size();
	DrawContext.dwStartVertex = 0;
	DrawContext.pXboxIndexData = nullptr;
	DrawContext.dwBaseVertexIndex = 0;
	DrawContext.pXboxVertexStreamZeroData = triList.data();
	DrawContext.uiXboxVertexStreamZeroStride = sizeof(Float3);
	DrawContext.bNV2AInlineData = true;

	CxbxD3D11VertexFetchDraw(DrawContext);

	pg->vertex_attributes[0] = savedAttr0;
}
