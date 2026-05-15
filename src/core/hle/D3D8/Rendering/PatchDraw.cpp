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

#include <cmath>
#include <cassert>

using namespace xbox;

// ---------------------------------------------------------------
// Patch grid limits — NV2A supports up to order 10 (9 FD terms)
// and maxSwatch <= 17, giving at most 18 output points per axis.
// ---------------------------------------------------------------
static constexpr int MAX_FD_ORDER = 16;
static constexpr int MAX_GRID_DIM = 32; // generous upper bound per axis

// ---------------------------------------------------------------
// Float3 helper
// ---------------------------------------------------------------
struct Float3 { float x, y, z; };

static Float3 MakeFloat3(const float *p) { return { p[0], p[1], p[2] }; }

// ---------------------------------------------------------------
// Forward-difference curve: load coefficients from float4 buffer
// ---------------------------------------------------------------
struct FDCurve3 {
	float c[MAX_FD_ORDER][3];
	int   order;

	void Load(const float *float4Base, int attrOffset, int fdOrder) {
		order = (fdOrder > MAX_FD_ORDER) ? MAX_FD_ORDER : fdOrder;
		for (int k = 0; k < order; k++) {
			int idx = (attrOffset + k) * 4;
			c[k][0] = float4Base[idx + 0];
			c[k][1] = float4Base[idx + 1];
			c[k][2] = float4Base[idx + 2];
		}
	}

	// Single FD step: c[i] += c[i+1] for i = 0..order-2
	void Step() {
		for (int i = 0; i < order - 1; i++) {
			c[i][0] += c[i + 1][0];
			c[i][1] += c[i + 1][1];
			c[i][2] += c[i + 1][2];
		}
	}

	Float3 Value() const { return { c[0][0], c[0][1], c[0][2] }; }

	// Exact endpoint after n steps via binomial sum (eliminates FD drift)
	Float3 ExactEndpoint(int n) const {
		double binom = 1.0;
		double ex = 0, ey = 0, ez = 0;
		for (int k = 0; k < order; k++) {
			ex += binom * c[k][0];
			ey += binom * c[k][1];
			ez += binom * c[k][2];
			binom = binom * (double)(n - k) / (double)(k + 1);
		}
		return { (float)ex, (float)ey, (float)ez };
	}
};

// ---------------------------------------------------------------
// Patch register decode helpers (BEGIN_PATCH0-3 layout)
// ---------------------------------------------------------------

// patch0/patch1: per-attribute U-order, 4 bits each (value = order-1, 0=disabled)
static int DecodeAttrUOrder(uint32_t patch0, uint32_t patch1, int hwAttr) {
	uint32_t reg = (hwAttr < 8) ? patch0 : patch1;
	int shift = (hwAttr % 8) * 4;
	int orderMinus1 = (reg >> shift) & 0xF;
	return orderMinus1 ? (orderMinus1 + 1) : 0;
}

// patch2: swatch geometry
struct PatchGeometry {
	int maxSwatch;
	int partialWidth;
	int partialHeight;
	int nSwatchU;
	int nSwatchV;
};

static PatchGeometry DecodePatch2(uint32_t patch2) {
	return {
		(int)((patch2 >> 16) & 0x1F),
		(int)((patch2 >> 21) & 0x1F),
		(int)((patch2 >> 26) & 0x1F),
		(int)((patch2 >> 8) & 0xFF),
		(int)((patch2 >> 0) & 0xFF),
	};
}

// patch3: V-orders, numCoeffs, flags
static int DecodePosVOrder(uint32_t patch3)     { return ((patch3 >> 6) & 0xF) + 1; }
static int DecodeNormalVOrder(uint32_t patch3)   { return ((patch3 >> 10) & 0xF) + 1; }
static bool DecodeHasNormal(uint32_t patch3)     { return ((patch3 >> 10) & 0xF) != 0; }
static int DecodeNumCoeffsPerRow(uint32_t patch3){ return (patch3 >> 24) & 0xFF; }

// ---------------------------------------------------------------
// Curve type constants (NV097_SET_BEGIN_END_CURVE_CMD values)
// ---------------------------------------------------------------
#define NV2A_CURVE_END_DATA         0
#define NV2A_CURVE_STRIP            1
#define NV2A_CURVE_LEFT_GUARD       2
#define NV2A_CURVE_RIGHT_GUARD      3
#define NV2A_CURVE_OUTER_TRANSITION 4
#define NV2A_CURVE_INNER_TRANSITION 5

extern void CxbxUpdateNativeD3DResources();

// ---------------------------------------------------------------
// Evaluate NV2A hardware tessellation patch and draw the result.
//
// Guard curve data layout:
//   (1) Guard point normal   [1 float4, if normal enabled]
//   (2) Guard point position [1 float4, the V-edge endpoint]
//   (3) Guard curve position [posVOrder float4s, FD coefficients along V-edge]
//   (4) Guard curve normal   [normalVOrder float4s, if normal enabled]
// LEFT guard  = left column (col 0) of the swatch grid
// RIGHT guard = right column (last col) of the swatch grid
// Guard curve FD[0] = exact start point of the column edge.
// Guard point = exact end point of the column edge.
// Hardware uses these to replace FD-stepped edge values, eliminating
// inter-patch gaps.
// ---------------------------------------------------------------
void D3D11_draw_patch(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;
	PatchState &patch = pg->patch;

	if (patch.curveCount == 0)
		return;


	// --- Decode patch registers ---

	// Find position attribute (first enabled hw attr) and its U-order
	int posUOrder = 0;
	int posAttrOffset = 0; // float4 offset within a strip row
	for (int a = 0; a < 16; a++) {
		int order = DecodeAttrUOrder(patch.patch0, patch.patch1, a);
		if (order >= 2) {
			posUOrder = order;
			break;
		}
		if (order > 0) posAttrOffset += order; // skip lower-order attrs before position
	}
	if (posUOrder < 2)
		return;

	int posVOrder       = DecodePosVOrder(patch.patch3);
	int numCoeffsPerRow = DecodeNumCoeffsPerRow(patch.patch3);
	bool hasNormal      = DecodeHasNormal(patch.patch3);

	PatchGeometry geo = DecodePatch2(patch.patch2);
	int numStepsU = (geo.nSwatchU == 0) ? geo.partialWidth : geo.maxSwatch;
	if (numStepsU < 1) numStepsU = 8;

	// --- Collect curve indices by type ---

	int stripIndices[NV2A_PATCH_MAX_CURVES];
	int numRows = 0;
	int leftGuard = -1, rightGuard = -1;
	int transitionU = -1, transitionV = -1; // U=inner(col), V=outer(row)

	for (int c = 0; c < patch.curveCount; c++) {
		switch (patch.curves[c].curveType) {
		case NV2A_CURVE_STRIP:            stripIndices[numRows++] = c; break;
		case NV2A_CURVE_LEFT_GUARD:       if (leftGuard < 0)  leftGuard  = c; break;
		case NV2A_CURVE_RIGHT_GUARD:      if (rightGuard < 0) rightGuard = c; break;
		case NV2A_CURVE_INNER_TRANSITION: if (transitionU < 0) transitionU = c; break;
		case NV2A_CURVE_OUTER_TRANSITION: if (transitionV < 0) transitionV = c; break;
		}
	}
	// --- Handle transition-only swatches (0 strips, IT + OT) ---
	if (numRows < 1) {
		if (transitionU < 0 && transitionV < 0)
			return;

		// Transition-only swatch: draws a connecting strip between two adjacent
		// patches by evaluating the inner/outer transition curves.
		// IT and OT are V-direction FD curves that run in OPPOSITE directions.
		// Pair IT[k] with OT[N-k] to match spatially corresponding points.

		int numStepsV = (geo.nSwatchV == 0) ? geo.partialHeight : geo.maxSwatch;
		if (numStepsV < 1) numStepsV = 8;
		int numPointsV = numStepsV + 1;

		// Choose which curves to use (prefer IT+OT pair)
		int curveA = (transitionU >= 0) ? transitionU : transitionV;
		int curveB = (transitionV >= 0 && transitionV != curveA) ? transitionV : transitionU;
		if (curveA < 0 || curveB < 0 || curveA == curveB)
			return;

		// Evaluate both transition curves as V-direction FD curves
		Float3 edgeA[MAX_GRID_DIM], edgeB[MAX_GRID_DIM];

		{
			PatchCurve &tc = patch.curves[curveA];
			const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
			FDCurve3 fd;
			fd.Load(coeffs, 0, posVOrder);
			edgeA[0] = fd.Value();
			for (int s = 1; s < numPointsV; s++) {
				fd.Step();
				edgeA[s] = fd.Value();
			}
		}
		{
			PatchCurve &tc = patch.curves[curveB];
			const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
			FDCurve3 fd;
			fd.Load(coeffs, 0, posVOrder);
			edgeB[0] = fd.Value();
			for (int s = 1; s < numPointsV; s++) {
				fd.Step();
				edgeB[s] = fd.Value();
			}
		}

		// Build triangle list: pair edgeA[k] with edgeB[N-k] (reversed)
		Float3 triList[MAX_GRID_DIM * 6];
		int triCount = 0;

		for (int k = 0; k < numPointsV - 1; k++) {
			int rk  = (numPointsV - 1) - k;
			int rk1 = (numPointsV - 1) - (k + 1);

			Float3 a0 = edgeA[k];
			Float3 a1 = edgeA[k + 1];
			Float3 b0 = edgeB[rk];
			Float3 b1 = edgeB[rk1];

			triList[triCount++] = a0;
			triList[triCount++] = b0;
			triList[triCount++] = a1;
			triList[triCount++] = b0;
			triList[triCount++] = b1;
			triList[triCount++] = a1;
		}

		if (triCount == 0)
			return;

		// Draw transition strip
		VertexAttribute savedAttr0 = pg->vertex_attributes[0];
		pg->vertex_attributes[0].format = 2;
		pg->vertex_attributes[0].count = 3;
		pg->vertex_attributes[0].stride = sizeof(Float3);
		pg->vertex_attributes[0].offset = 0;

		CxbxUpdateNativeD3DResources();

		CxbxDrawContext DrawContext = {};
		DrawContext.XboxPrimitiveType = X_D3DPT_TRIANGLELIST;
		DrawContext.dwVertexCount = (DWORD)triCount;
		DrawContext.dwStartVertex = 0;
		DrawContext.pXboxIndexData = nullptr;
		DrawContext.dwBaseVertexIndex = 0;
		DrawContext.pXboxVertexStreamZeroData = triList;
		DrawContext.uiXboxVertexStreamZeroStride = sizeof(Float3);
		DrawContext.bNV2AInlineData = true;

		CxbxD3D11VertexFetchDraw(DrawContext);

		pg->vertex_attributes[0] = savedAttr0;
		return;
	}

	// --- Determine grid dimensions ---

	int numOutU = numStepsU + 1;

	// Detect missing boundary row: when we receive fewer strips than expected
	// for this swatch, the V=0 row was omitted (shared with an adjacent patch).
	// The guard curves span the FULL V-range including the omitted row.
	bool hasMissingFirstRow = (numRows < geo.partialHeight + 1) && (leftGuard >= 0) && (rightGuard >= 0);
	int totalRows = hasMissingFirstRow ? (numRows + 1) : numRows;

	bool hasExtraCol = (transitionU >= 0) && (numOutU < totalRows);
	bool hasExtraRow = (transitionV >= 0) && (totalRows < numOutU);

	int finalU = numOutU + (hasExtraCol ? 1 : 0);
	int finalV = totalRows + (hasExtraRow ? 1 : 0);

	assert(finalU <= MAX_GRID_DIM && finalV <= MAX_GRID_DIM);

	// --- Evaluate FD grid from strip curves ---

	Float3 grid[MAX_GRID_DIM * MAX_GRID_DIM];

	// If boundary row is missing, strips fill rows 1..numRows; row 0 is synthesized below
	int rowOffset = hasMissingFirstRow ? 1 : 0;

	for (int row = 0; row < numRows; row++) {
		PatchCurve &curve = patch.curves[stripIndices[row]];
		const float *coeffs = &patch.coefficients[curve.coeffStart * 4];

		FDCurve3 fd;
		fd.Load(coeffs, posAttrOffset, posUOrder);

		Float3 exact = fd.ExactEndpoint(numStepsU);
		Float3 *out = &grid[(row + rowOffset) * finalU];

		out[0] = fd.Value();
		for (int step = 1; step < numOutU; step++) {
			fd.Step();
			out[step] = fd.Value();
		}
		out[numStepsU] = exact; // replace last with drift-free value
	}

	// --- Synthesize missing first row from guard curves ---
	if (hasMissingFirstRow && leftGuard >= 0 && rightGuard >= 0) {
		// Guard curve layout: [normal_point?] [position_point] [position_FD(posVOrder)] [normal_FD]
		// FD[0] is the V=0 value for that column edge.
		auto GetGuardV0 = [&](int curveIdx) -> Float3 {
			PatchCurve &gc = patch.curves[curveIdx];
			const float *data = &patch.coefficients[gc.coeffStart * 4];
			int offset = hasNormal ? 1 : 0; // skip normal point
			offset++; // skip position endpoint
			// Position FD[0] = V=0 value
			return MakeFloat3(&data[offset * 4]);
		};

		Float3 v0Left  = GetGuardV0(leftGuard);
		Float3 v0Right = GetGuardV0(rightGuard);

		// First strip row (now at grid row 1) provides curvature reference
		Float3 strip0Left  = grid[1 * finalU + 0];
		Float3 strip0Right = grid[1 * finalU + numStepsU];

		for (int i = 0; i < numOutU; i++) {
			float t = (numStepsU > 0) ? (float)i / (float)numStepsU : 0.0f;

			// Linear blend between guard V=0 edge values
			Float3 base;
			base.x = v0Left.x + t * (v0Right.x - v0Left.x);
			base.y = v0Left.y + t * (v0Right.y - v0Left.y);
			base.z = v0Left.z + t * (v0Right.z - v0Left.z);

			// Add curvature from the first strip: deviation from its linear interpolant
			Float3 stripLinear;
			stripLinear.x = strip0Left.x + t * (strip0Right.x - strip0Left.x);
			stripLinear.y = strip0Left.y + t * (strip0Right.y - strip0Left.y);
			stripLinear.z = strip0Left.z + t * (strip0Right.z - strip0Left.z);

			Float3 curvature;
			curvature.x = grid[1 * finalU + i].x - stripLinear.x;
			curvature.y = grid[1 * finalU + i].y - stripLinear.y;
			curvature.z = grid[1 * finalU + i].z - stripLinear.z;

			grid[0 * finalU + i].x = base.x + curvature.x;
			grid[0 * finalU + i].y = base.y + curvature.y;
			grid[0 * finalU + i].z = base.z + curvature.z;
		}
	}

	// --- Apply guard curves to edge columns (watertight stitching) ---
	{
		auto ApplyGuardColumn = [&](int curveIdx, int col) {
			PatchCurve &gc = patch.curves[curveIdx];
			const float *data = &patch.coefficients[gc.coeffStart * 4];
			int offset = 0;
			if (hasNormal) offset++; // skip normal endpoint
			// Position endpoint (TOP of column = last row)
			Float3 endpoint = MakeFloat3(&data[offset * 4]);
			offset++;
			// Position FD curve (BOTTOM of column = first row, step per row)
			FDCurve3 guardFD;
			guardFD.Load(data, offset, posVOrder);

			// Apply guard values to each row except the last
			for (int row = 0; row < totalRows - 1; row++) {
				grid[row * finalU + col] = guardFD.Value();
				guardFD.Step();
			}
			// Last row: use exact endpoint (eliminates FD drift)
			grid[(totalRows - 1) * finalU + col] = endpoint;
		};

		if (leftGuard >= 0)
			ApplyGuardColumn(leftGuard, 0);
		if (rightGuard >= 0)
			ApplyGuardColumn(rightGuard, numOutU - 1);
	}

	// --- Transition curves (stitch adjacent patches at different tess levels) ---

	auto EvalTransitionColumn = [&](int curveIdx, int col) {
		PatchCurve &tc = patch.curves[curveIdx];
		const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
		int tOffset = (tc.coeffCount == numCoeffsPerRow) ? posAttrOffset : 0;

		FDCurve3 fd;
		fd.Load(coeffs, tOffset, posVOrder);

		grid[0 * finalU + col] = fd.Value();
		for (int row = 1; row < totalRows; row++) {
			fd.Step();
			grid[row * finalU + col] = fd.Value();
		}
	};

	auto EvalTransitionRow = [&](int curveIdx, int row) {
		PatchCurve &tc = patch.curves[curveIdx];
		const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
		int tOffset = (tc.coeffCount == numCoeffsPerRow) ? posAttrOffset : 0;

		FDCurve3 fd;
		fd.Load(coeffs, tOffset, posUOrder);

		grid[row * finalU + 0] = fd.Value();
		for (int step = 1; step < numOutU; step++) {
			fd.Step();
			grid[row * finalU + step] = fd.Value();
		}
	};

	if (hasExtraCol)
		EvalTransitionColumn(transitionU, numOutU);

	if (hasExtraRow) {
		EvalTransitionRow(transitionV, totalRows);

		// Corner where both transitions meet
		if (hasExtraCol) {
			PatchCurve &tc = patch.curves[transitionU];
			const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
			int tOffset = (tc.coeffCount == numCoeffsPerRow) ? posAttrOffset : 0;

			FDCurve3 fd;
			fd.Load(coeffs, tOffset, posVOrder);
			for (int s = 0; s < totalRows; s++)
				fd.Step();
			grid[totalRows * finalU + numOutU] = fd.Value();
		}
	}

	// --- Build triangle list and draw ---

	Float3 triList[MAX_GRID_DIM * MAX_GRID_DIM * 6];
	int triCount = 0;

	for (int y = 0; y < finalV - 1; y++) {
		for (int x = 0; x < finalU - 1; x++) {
			const Float3 &p00 = grid[y * finalU + x];
			const Float3 &p10 = grid[y * finalU + x + 1];
			const Float3 &p01 = grid[(y + 1) * finalU + x];
			const Float3 &p11 = grid[(y + 1) * finalU + x + 1];

			triList[triCount++] = p00;
			triList[triCount++] = p10;
			triList[triCount++] = p01;
			triList[triCount++] = p10;
			triList[triCount++] = p11;
			triList[triCount++] = p01;
		}
	}

	if (triCount == 0)
		return;

	// Route through VertexFetch as a UP (user-pointer) draw
	VertexAttribute savedAttr0 = pg->vertex_attributes[0];

	pg->vertex_attributes[0].format = 2; // NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_F
	pg->vertex_attributes[0].count = 3;
	pg->vertex_attributes[0].stride = sizeof(Float3);
	pg->vertex_attributes[0].offset = 0;

	CxbxUpdateNativeD3DResources();

	CxbxDrawContext DrawContext = {};
	DrawContext.XboxPrimitiveType = X_D3DPT_TRIANGLELIST;
	DrawContext.dwVertexCount = (DWORD)triCount;
	DrawContext.dwStartVertex = 0;
	DrawContext.pXboxIndexData = nullptr;
	DrawContext.dwBaseVertexIndex = 0;
	DrawContext.pXboxVertexStreamZeroData = triList;
	DrawContext.uiXboxVertexStreamZeroStride = sizeof(Float3);
	DrawContext.bNV2AInlineData = true;

	CxbxD3D11VertexFetchDraw(DrawContext);

	pg->vertex_attributes[0] = savedAttr0;
}
