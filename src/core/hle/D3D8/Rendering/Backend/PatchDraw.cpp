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

#include "../RenderGlobals.h"
#include "Backend_D3D11.h" // CxbxD3D11VertexFetchDraw, CxbxUpdateNativeD3DResources
#include "core/hle/D3D8/XbVertexBuffer.h" // For CxbxDrawContext
#include "core/kernel/support/Emu.h"
#include "devices/video/nv2a_int.h" // For PGRAPHState, PatchState

#include <cmath>
#include <cassert>
#include <cstring>

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

static Float3 Float3Sub(Float3 a, Float3 b) {
	return { a.x - b.x, a.y - b.y, a.z - b.z };
}

static Float3 Float3Cross(Float3 a, Float3 b) {
	return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}

static Float3 Float3Normalize(Float3 v) {
	float len = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
	if (len > 1e-8f) { float inv = 1.0f / len; return { v.x*inv, v.y*inv, v.z*inv }; }
	return { 0.0f, 0.0f, 1.0f };
}

// Compute per-vertex normals from a position grid using finite differences.
// Grid is rows × cols. Output normals array must be at least rows*cols.
static void ComputeGridNormals(const Float3 *grid, Float3 *normals, int rows, int cols) {
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			// Tangent in U (column) direction
			Float3 dU;
			if (x > 0 && x < cols - 1)
				dU = Float3Sub(grid[y*cols + x+1], grid[y*cols + x-1]);
			else if (x < cols - 1)
				dU = Float3Sub(grid[y*cols + x+1], grid[y*cols + x]);
			else
				dU = Float3Sub(grid[y*cols + x], grid[y*cols + x-1]);

			// Tangent in V (row) direction
			Float3 dV;
			if (y > 0 && y < rows - 1)
				dV = Float3Sub(grid[(y+1)*cols + x], grid[(y-1)*cols + x]);
			else if (y < rows - 1)
				dV = Float3Sub(grid[(y+1)*cols + x], grid[y*cols + x]);
			else
				dV = Float3Sub(grid[y*cols + x], grid[(y-1)*cols + x]);

			normals[y*cols + x] = Float3Normalize(Float3Cross(dU, dV));
		}
	}
}

// ---------------------------------------------------------------
// Forward-difference curve: load coefficients from float4 buffer
// ---------------------------------------------------------------
struct FDCurve4 {
	float c[MAX_FD_ORDER][4];
	int   order;

	void Load(const float *float4Base, int attrOffset, int fdOrder) {
		order = (fdOrder > MAX_FD_ORDER) ? MAX_FD_ORDER : fdOrder;
		for (int k = 0; k < order; k++) {
			int idx = (attrOffset + k) * 4;
			c[k][0] = float4Base[idx + 0];
			c[k][1] = float4Base[idx + 1];
			c[k][2] = float4Base[idx + 2];
			c[k][3] = float4Base[idx + 3];
		}
	}

	void Step() {
		for (int i = 0; i < order - 1; i++) {
			c[i][0] += c[i + 1][0];
			c[i][1] += c[i + 1][1];
			c[i][2] += c[i + 1][2];
			c[i][3] += c[i + 1][3];
		}
	}

	void Value(float *out) const {
		out[0] = c[0][0]; out[1] = c[0][1]; out[2] = c[0][2]; out[3] = c[0][3];
	}

	void ExactEndpoint(int n, float *out) const {
		double binom = 1.0;
		double r[4] = { 0, 0, 0, 0 };
		for (int k = 0; k < order; k++) {
			for (int j = 0; j < 4; j++) r[j] += binom * c[k][j];
			binom = binom * (double)(n - k) / (double)(k + 1);
		}
		for (int j = 0; j < 4; j++) out[j] = (float)r[j];
	}

	// Convenience: extract xyz as Float3
	Float3 Value3() const { return { c[0][0], c[0][1], c[0][2] }; }

	Float3 ExactEndpoint3(int n) const {
		float tmp[4];
		ExactEndpoint(n, tmp);
		return { tmp[0], tmp[1], tmp[2] };
	}
};

// ---------------------------------------------------------------
// Multi-attribute tessellation support
// ---------------------------------------------------------------

// Describes one tessellated vertex attribute
struct TessAttr {
	int hwIndex;        // NV2A attribute index (0-15)
	int uOrder;         // FD order in U direction (0 for auto attrs)
	int coeffOffset;    // float4 offset within strip row (for FD attrs)
	int outByteOffset;  // byte offset in packed output vertex
	int components;     // number of output float components (2, 3, or 4)
	bool isFD;          // true = evaluate from FD data
	bool isAutoNormal;  // true = derived from position cross product
};

static constexpr int MAX_TESS_ATTRS = 16;

// Per-attribute grid storage for non-position FD attributes (4 floats per grid point).
// Indexed by TessAttr array index; position uses its own Float3 grid.
static float s_fdGrids[MAX_TESS_ATTRS][MAX_GRID_DIM * MAX_GRID_DIM * 4];

// Output triangle list buffer (static to avoid stack overflow with large vertex strides)
static constexpr int MAX_TRI_BUF_BYTES = 512 * 1024;
static uint8_t s_triListBuf[MAX_TRI_BUF_BYTES];

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

	int posVOrder       = DecodePosVOrder(patch.patch3);
	int numCoeffsPerRow = DecodeNumCoeffsPerRow(patch.patch3);
	bool hasNormal      = DecodeHasNormal(patch.patch3);

	// Scan all 16 HW attributes for active tessellation targets.
	// Attributes with non-zero U-order in patch0/patch1 are FD-tessellated.
	// Attr 2 (normal) is auto-generated from position cross product when
	// hasNormal is set but attr 2 has no FD data of its own.
	TessAttr tessAttrs[MAX_TESS_ATTRS];
	int numTessAttrs = 0;
	int posIdx = -1;        // index into tessAttrs[] for position
	int normalIdx = -1;     // index into tessAttrs[] for auto-normal
	int posUOrder = 0;
	int posAttrOffset = 0;  // float4 offset of position within strip row

	{
		int coeffOff = 0;
		int byteOff = 0;
		for (int a = 0; a < 16; a++) {
			int order = DecodeAttrUOrder(patch.patch0, patch.patch1, a);
			if (order > 0) {
				// FD-tessellated attribute
				TessAttr &ta = tessAttrs[numTessAttrs];
				ta.hwIndex = a;
				ta.uOrder = order;
				ta.coeffOffset = coeffOff;
				ta.isFD = true;
				ta.isAutoNormal = false;

				if (posIdx < 0 && order >= 2) {
					// First attr with order >= 2 is position
					posIdx = numTessAttrs;
					posUOrder = order;
					posAttrOffset = coeffOff;
					ta.components = 3; // Position = xyz
				} else if (a == 2) {
					// Normal from FD data (driver-computed via TESSNORMAL)
					ta.components = 3;
				} else {
					ta.components = 4; // General attr = full float4
				}

				ta.outByteOffset = byteOff;
				byteOff += ta.components * (int)sizeof(float);
				coeffOff += order;
				numTessAttrs++;
			} else if (a == 2 && hasNormal) {
				// Auto-normal: derived from position cross product
				TessAttr &ta = tessAttrs[numTessAttrs];
				ta.hwIndex = 2;
				ta.uOrder = 0;
				ta.coeffOffset = 0;
				ta.components = 3;
				ta.isFD = false;
				ta.isAutoNormal = true;
				ta.outByteOffset = byteOff;
				byteOff += 3 * (int)sizeof(float);
				normalIdx = numTessAttrs;
				numTessAttrs++;
			}
		}
	}

	if (posIdx < 0 || posUOrder < 2)
		return;

	int outVertexStride = 0;
	for (int i = 0; i < numTessAttrs; i++)
		outVertexStride += tessAttrs[i].components * (int)sizeof(float);

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
			FDCurve4 fd;
			fd.Load(coeffs, 0, posVOrder);
			edgeA[0] = fd.Value3();
			for (int s = 1; s < numPointsV; s++) {
				fd.Step();
				edgeA[s] = fd.Value3();
			}
		}
		{
			PatchCurve &tc = patch.curves[curveB];
			const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
			FDCurve4 fd;
			fd.Load(coeffs, 0, posVOrder);
			edgeB[0] = fd.Value3();
			for (int s = 1; s < numPointsV; s++) {
				fd.Step();
				edgeB[s] = fd.Value3();
			}
		}

		// Transition strips always output position + auto-normal (cross product)
		int transStride = 6 * (int)sizeof(float); // pos(3) + normal(3)
		int transNormalOff = 3 * (int)sizeof(float);

		int maxTransVerts = (numPointsV - 1) * 6;
		int transBytes = maxTransVerts * transStride;
		assert(transBytes <= MAX_TRI_BUF_BYTES);
		if (transBytes > MAX_TRI_BUF_BYTES)
			return;

		int triCount = 0;

		for (int k = 0; k < numPointsV - 1; k++) {
			int rk  = (numPointsV - 1) - k;
			int rk1 = (numPointsV - 1) - (k + 1);

			Float3 a0 = edgeA[k];
			Float3 a1 = edgeA[k + 1];
			Float3 b0 = edgeB[rk];
			Float3 b1 = edgeB[rk1];

			// Compute face normal for this quad (both triangles share it)
			Float3 e1 = Float3Sub(b0, a0);
			Float3 e2 = Float3Sub(a1, a0);
			Float3 n = Float3Normalize(Float3Cross(e1, e2));

			Float3 verts[6] = { a0, b0, a1, b0, b1, a1 };
			for (int v = 0; v < 6; v++) {
				uint8_t *dst = &s_triListBuf[triCount * transStride];
				float *posDst = (float *)dst;
				posDst[0] = verts[v].x;
				posDst[1] = verts[v].y;
				posDst[2] = verts[v].z;
				float *nrmDst = (float *)(dst + transNormalOff);
				nrmDst[0] = n.x;
				nrmDst[1] = n.y;
				nrmDst[2] = n.z;
				triCount++;
			}
		}

		if (triCount == 0)
			return;

		// Draw transition strip: save ALL 16 attributes, disable unused ones
		int posHw = tessAttrs[posIdx].hwIndex;
		VertexAttribute savedTransAttrs[16];
		for (int i = 0; i < 16; i++) {
			savedTransAttrs[i] = pg->vertex_attributes[i];
			pg->vertex_attributes[i].count = 0; // disable
		}

		pg->vertex_attributes[posHw].format = 2;
		pg->vertex_attributes[posHw].size = 4;
		pg->vertex_attributes[posHw].count = 3;
		pg->vertex_attributes[posHw].stride = transStride;
		pg->vertex_attributes[posHw].offset = 0;

		pg->vertex_attributes[2].format = 2;
		pg->vertex_attributes[2].size = 4;
		pg->vertex_attributes[2].count = 3;
		pg->vertex_attributes[2].stride = transStride;
		pg->vertex_attributes[2].offset = 0;

		CxbxUpdateNativeD3DResources();

		CxbxDrawContext DrawContext = {};
		DrawContext.XboxPrimitiveType = X_D3DPT_TRIANGLELIST;
		DrawContext.dwVertexCount = (DWORD)triCount;
		DrawContext.dwStartVertex = 0;
		DrawContext.pXboxIndexData = nullptr;
		DrawContext.dwBaseVertexIndex = 0;
		DrawContext.pXboxVertexStreamZeroData = s_triListBuf;
		DrawContext.uiXboxVertexStreamZeroStride = transStride;
		DrawContext.bNV2AInlineData = true;

		CxbxD3D11VertexFetchDraw(DrawContext);

		for (int i = 0; i < 16; i++)
			pg->vertex_attributes[i] = savedTransAttrs[i];
		return;
	}

	// --- Determine grid dimensions ---

	int numOutU = numStepsU + 1;

	// Detect missing boundary row: when the driver omits the V=0 row (shared
	// with an adjacent swatch), the guard curve spans more V-rows than we
	// received strips.  We detect this using the guard's own data: the stored
	// endpoint is the exact position at the last grid row, and the FD curve
	// reconstructs it via ExactEndpoint(N) where N = totalRows - 1.
	// If ExactEndpoint(numRows - 1) matches the endpoint, the guard spans
	// exactly numRows (no missing row).  If it doesn't, the guard spans
	// numRows + 1 and the first row was omitted.
	bool hasMissingFirstRow = false;
	if (leftGuard >= 0 && numRows > 0) {
		PatchCurve &gc = patch.curves[leftGuard];
		const float *data = &patch.coefficients[gc.coeffStart * 4];
		int off = hasNormal ? 1 : 0;      // skip normal endpoint
		Float3 endpoint = MakeFloat3(&data[off * 4]);
		off++;                             // skip position endpoint
		FDCurve4 guardFD;
		guardFD.Load(data, off, posVOrder);

		Float3 ep = guardFD.ExactEndpoint3(numRows - 1);
		float dx = ep.x - endpoint.x;
		float dy = ep.y - endpoint.y;
		float dz = ep.z - endpoint.z;
		float distSq_nm1 = dx*dx + dy*dy + dz*dz;

		Float3 ep2 = guardFD.ExactEndpoint3(numRows);
		float dx2 = ep2.x - endpoint.x;
		float dy2 = ep2.y - endpoint.y;
		float dz2 = ep2.z - endpoint.z;
		float distSq_n = dx2*dx2 + dy2*dy2 + dz2*dz2;

		hasMissingFirstRow = (distSq_nm1 > 1e-3f);
	}
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

		FDCurve4 fd;
		fd.Load(coeffs, posAttrOffset, posUOrder);

		Float3 exact = fd.ExactEndpoint3(numStepsU);
		Float3 *out = &grid[(row + rowOffset) * finalU];

		out[0] = fd.Value3();
		for (int step = 1; step < numOutU; step++) {
			fd.Step();
			out[step] = fd.Value3();
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
			FDCurve4 guardFD;
			guardFD.Load(data, offset, posVOrder);

			// Apply guard values to each row except the last
			for (int row = 0; row < totalRows - 1; row++) {
				grid[row * finalU + col] = guardFD.Value3();
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

		FDCurve4 fd;
		fd.Load(coeffs, tOffset, posVOrder);

		grid[0 * finalU + col] = fd.Value3();
		for (int row = 1; row < totalRows; row++) {
			fd.Step();
			grid[row * finalU + col] = fd.Value3();
		}
	};

	auto EvalTransitionRow = [&](int curveIdx, int row) {
		PatchCurve &tc = patch.curves[curveIdx];
		const float *coeffs = &patch.coefficients[tc.coeffStart * 4];
		int tOffset = (tc.coeffCount == numCoeffsPerRow) ? posAttrOffset : 0;

		FDCurve4 fd;
		fd.Load(coeffs, tOffset, posUOrder);

		grid[row * finalU + 0] = fd.Value3();
		for (int step = 1; step < numOutU; step++) {
			fd.Step();
			grid[row * finalU + step] = fd.Value3();
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

			FDCurve4 fd;
			fd.Load(coeffs, tOffset, posVOrder);
			for (int s = 0; s < totalRows; s++)
				fd.Step();
			grid[totalRows * finalU + numOutU] = fd.Value3();
		}
	}

	// --- Evaluate extra FD attribute grids from strip curves ---
	// Position uses its own Float3 grid (already evaluated above with guards/transitions).
	// Other FD attributes are evaluated here from strip row data only.

	for (int ai = 0; ai < numTessAttrs; ai++) {
		TessAttr &ta = tessAttrs[ai];
		if (ai == posIdx || !ta.isFD)
			continue;

		for (int row = 0; row < numRows; row++) {
			PatchCurve &curve = patch.curves[stripIndices[row]];
			const float *coeffs = &patch.coefficients[curve.coeffStart * 4];

			FDCurve4 fd;
			fd.Load(coeffs, ta.coeffOffset, ta.uOrder);

			float exact[4];
			fd.ExactEndpoint(numStepsU, exact);

			float *out = &s_fdGrids[ai][(row + rowOffset) * finalU * 4];
			fd.Value(&out[0]);
			for (int step = 1; step < numOutU; step++) {
				fd.Step();
				fd.Value(&out[step * 4]);
			}
			memcpy(&out[numStepsU * 4], exact, 4 * sizeof(float));
		}

		// Handle missing first row: copy from first available strip row
		if (hasMissingFirstRow) {
			memcpy(&s_fdGrids[ai][0], &s_fdGrids[ai][1 * finalU * 4],
				finalU * 4 * sizeof(float));
		}

		// Fill transition columns/rows with nearest interior value
		if (hasExtraCol) {
			for (int y = 0; y < totalRows; y++)
				memcpy(&s_fdGrids[ai][(y * finalU + numOutU) * 4],
					&s_fdGrids[ai][(y * finalU + numOutU - 1) * 4], 4 * sizeof(float));
		}
		if (hasExtraRow) {
			memcpy(&s_fdGrids[ai][totalRows * finalU * 4],
				&s_fdGrids[ai][(totalRows - 1) * finalU * 4],
				finalU * 4 * sizeof(float));
		}
	}

	// --- Compute auto-normals from position grid (if needed) ---

	Float3 normals[MAX_GRID_DIM * MAX_GRID_DIM];
	if (normalIdx >= 0) {
		ComputeGridNormals(grid, normals, finalV, finalU);
	}

	// --- Build packed triangle list with all active attributes ---

	int maxVertices = (finalV - 1) * (finalU - 1) * 6;
	int bytesNeeded = maxVertices * outVertexStride;
	assert(bytesNeeded <= MAX_TRI_BUF_BYTES);
	if (bytesNeeded > MAX_TRI_BUF_BYTES)
		return;

	int triCount = 0;

	for (int y = 0; y < finalV - 1; y++) {
		for (int x = 0; x < finalU - 1; x++) {
			int gi[4] = {
				y * finalU + x,
				y * finalU + x + 1,
				(y + 1) * finalU + x,
				(y + 1) * finalU + x + 1
			};
			// Two triangles: (00, 10, 01), (10, 11, 01)
			int triVerts[6] = { gi[0], gi[1], gi[2], gi[1], gi[3], gi[2] };

			for (int v = 0; v < 6; v++) {
				uint8_t *dst = &s_triListBuf[triCount * outVertexStride];
				int gidx = triVerts[v];

				for (int ai = 0; ai < numTessAttrs; ai++) {
					TessAttr &ta = tessAttrs[ai];
					float *attrDst = (float *)(dst + ta.outByteOffset);

					if (ai == posIdx) {
						attrDst[0] = grid[gidx].x;
						attrDst[1] = grid[gidx].y;
						attrDst[2] = grid[gidx].z;
					} else if (ta.isAutoNormal) {
						attrDst[0] = normals[gidx].x;
						attrDst[1] = normals[gidx].y;
						attrDst[2] = normals[gidx].z;
					} else if (ta.isFD) {
						const float *src = &s_fdGrids[ai][gidx * 4];
						for (int c = 0; c < ta.components; c++)
							attrDst[c] = src[c];
					}
				}
				triCount++;
			}
		}
	}

	if (triCount == 0)
		return;

	// --- Set up vertex attributes for all active tessellation outputs ---

	// Save ALL 16 attributes and disable unused ones so the inline_array
	// packer doesn't include stale attributes from the game's prior draws.
	VertexAttribute savedAttrs[16];
	for (int i = 0; i < 16; i++) {
		savedAttrs[i] = pg->vertex_attributes[i];
		pg->vertex_attributes[i].count = 0; // disable
	}
	for (int ai = 0; ai < numTessAttrs; ai++) {
		int hw = tessAttrs[ai].hwIndex;
		pg->vertex_attributes[hw].format = 2; // NV097_SET_VERTEX_DATA_ARRAY_FORMAT_TYPE_F
		pg->vertex_attributes[hw].size = 4;   // sizeof(float)
		pg->vertex_attributes[hw].count = tessAttrs[ai].components;
		pg->vertex_attributes[hw].stride = outVertexStride;
		pg->vertex_attributes[hw].offset = 0;
	}

	CxbxUpdateNativeD3DResources();

	CxbxDrawContext DrawContext = {};
	DrawContext.XboxPrimitiveType = X_D3DPT_TRIANGLELIST;
	DrawContext.dwVertexCount = (DWORD)triCount;
	DrawContext.dwStartVertex = 0;
	DrawContext.pXboxIndexData = nullptr;
	DrawContext.dwBaseVertexIndex = 0;
	DrawContext.pXboxVertexStreamZeroData = s_triListBuf;
	DrawContext.uiXboxVertexStreamZeroStride = outVertexStride;
	DrawContext.bNV2AInlineData = true;

	CxbxD3D11VertexFetchDraw(DrawContext);

	for (int i = 0; i < 16; i++)
		pg->vertex_attributes[i] = savedAttrs[i];
}
