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
// *  CPU-side tessellation of rect and tri patches for D3D11
// *  (D3D11 has no DrawRectPatch/DrawTriPatch API)
// *
// *  The Xbox D3D8 API exposes DrawRectPatch and DrawTriPatch which
// *  tessellate B-spline/Bezier/Catmull-Rom patches on the GPU.
// *  D3D9 forwarded these to hardware. D3D11 removed this feature
// *  entirely, so we must perform the tessellation on the CPU and
// *  emit a triangle mesh for rendering.
// *
// *  FUTURE: Replace this CPU tessellation with D3D11 hardware
// *  tessellation (Hull Shader + Tessellator + Domain Shader):
// *   - Upload control points as a patch list primitive topology
// *     (e.g. D3D11_PRIMITIVE_TOPOLOGY_16_CONTROL_POINT_PATCHLIST)
// *   - Hull Shader: pass through control points, output tess factors
// *     from pNumSegs
// *   - Domain Shader: evaluate basis functions at each (u,v) sample
// *     point to produce position. AUTONORMAL can be computed
// *     analytically via normalize(cross(dP/du, dP/dv)) using the
// *     derivative basis functions (exact, not finite-difference).
// *     AUTOTEXCOORD is simply the (u,v) domain coordinates.
// *   - Each basis type (Bezier/B-spline/Catmull-Rom) x degree
// *     needs its own HS/DS pair (or a single pair with branching).
// *   - The existing VS logic would need to be split across VS+DS
// *     or incorporated into the DS, since the pipeline becomes
// *     VS -> HS -> Tessellator -> DS -> PS.
// *  This would be faster, more accurate, and architecturally cleaner.
// ******************************************************************

#include "PatchDraw.h"


#define LOG_PREFIX CXBXR_MODULE::D3D8

#include "RenderGlobals.h"
#include "Backend\Backend_D3D11.h" // FilterInputElementsByShaderSignature, CxbxD3D11VertexFetchDraw
#include "core/hle/D3D8/XbVertexShader.h"
#include "core/hle/D3D8/XbVertexBuffer.h" // For CxbxDrawContext
#include "core/kernel/support/Emu.h"
#include "devices/video/nv2a_int.h" // For PGRAPHState, PatchState
#include "common/AddressRanges.h"   // For CONTIGUOUS_MEMORY_BASE

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>

using namespace xbox;

// ---------------------------------------------------------------
// Utility: Evaluate 1D Bernstein basis for cubic Bezier (degree 3)
// ---------------------------------------------------------------
static inline float BernsteinBasis3(int i, float t)
{
	float s = 1.0f - t;
	switch (i) {
	case 0: return s * s * s;
	case 1: return 3.0f * s * s * t;
	case 2: return 3.0f * s * t * t;
	case 3: return t * t * t;
	default: return 0.0f;
	}
}

// ---------------------------------------------------------------
// Utility: Evaluate 1D Bernstein derivative for cubic Bezier
// ---------------------------------------------------------------
static inline float BernsteinDeriv3(int i, float t)
{
	float s = 1.0f - t;
	switch (i) {
	case 0: return -3.0f * s * s;
	case 1: return 3.0f * s * s - 6.0f * s * t;
	case 2: return 6.0f * s * t - 3.0f * t * t;
	case 3: return 3.0f * t * t;
	default: return 0.0f;
	}
}

// ---------------------------------------------------------------
// Utility: Evaluate 1D B-spline basis (uniform cubic, degree 3)
// The basis functions for a uniform cubic B-spline with knots at 0,1,2,3,4
// evaluated over the [0,1] parameter range (knot span [1,2])
// ---------------------------------------------------------------
static inline float BSplineBasis3(int i, float t)
{
	float s = 1.0f - t;
	switch (i) {
	case 0: return (s * s * s) / 6.0f;
	case 1: return (3.0f * t * t * t - 6.0f * t * t + 4.0f) / 6.0f;
	case 2: return (-3.0f * t * t * t + 3.0f * t * t + 3.0f * t + 1.0f) / 6.0f;
	case 3: return (t * t * t) / 6.0f;
	default: return 0.0f;
	}
}

// ---------------------------------------------------------------
// Utility: Evaluate 1D Catmull-Rom basis (cubic, degree 3)
// ---------------------------------------------------------------
static inline float CatmullRomBasis3(int i, float t)
{
	float t2 = t * t;
	float t3 = t2 * t;
	switch (i) {
	case 0: return (-t3 + 2.0f * t2 - t) / 2.0f;
	case 1: return (3.0f * t3 - 5.0f * t2 + 2.0f) / 2.0f;
	case 2: return (-3.0f * t3 + 4.0f * t2 + t) / 2.0f;
	case 3: return (t3 - t2) / 2.0f;
	default: return 0.0f;
	}
}

// ---------------------------------------------------------------
// Utility: Evaluate 1D linear basis (degree 1)
// ---------------------------------------------------------------
static inline float LinearBasis(int i, float t)
{
	switch (i) {
	case 0: return 1.0f - t;
	case 1: return t;
	default: return 0.0f;
	}
}

// ---------------------------------------------------------------
// Get the basis order (number of control points per row/col) for
// a given degree type
// ---------------------------------------------------------------
static int GetBasisOrder(X_D3DBASISTYPE Basis, X_D3DDEGREETYPE Degree)
{
	if (Degree == D3DDEGREE_LINEAR)
		return 2;
	// For cubic patches (Bezier, B-spline, Catmull-Rom)
	return 4;
}

// ---------------------------------------------------------------
// Evaluate 1D basis function for the given basis type
// ---------------------------------------------------------------
typedef float (*BasisFunc1D)(int i, float t);

static BasisFunc1D GetBasisFunction(X_D3DBASISTYPE Basis, X_D3DDEGREETYPE Degree)
{
	if (Degree == D3DDEGREE_LINEAR)
		return LinearBasis;

	switch (Basis) {
	case D3DBASIS_BEZIER:      return BernsteinBasis3;
	case D3DBASIS_BSPLINE:     return BSplineBasis3;
	case D3DBASIS_CATMULL_ROM: return CatmullRomBasis3;
	default:                   return BernsteinBasis3;
	}
}

// ---------------------------------------------------------------
// Float3 / Float4 types
// ---------------------------------------------------------------
struct Float3 { float x, y, z; };
struct Float4 { float x, y, z, w; };

static Float3 ReadVertexPosition(const uint8_t *pVertexData, UINT stride, UINT index)
{
	const float *p = reinterpret_cast<const float*>(pVertexData + index * stride);
	return { p[0], p[1], p[2] };
}

static Float3 Float3Cross(Float3 a, Float3 b)
{
	return { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}

static Float3 Float3Normalize(Float3 v)
{
	float len = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	if (len < 1e-12f) return { 0, 1, 0 };
	return { v.x/len, v.y/len, v.z/len };
}

// ---------------------------------------------------------------
// Evaluate a rectangular patch position grid
// Returns a 2D grid (totalSamplesU x totalSamplesV) of Float3.
// ---------------------------------------------------------------
static bool EvaluateRectPatchGrid(
	const uint8_t *pVertexData,
	UINT stride,
	UINT startU, UINT startV,
	UINT cpWidth, UINT cpHeight,
	UINT cpStride,
	X_D3DBASISTYPE basis,
	X_D3DDEGREETYPE degree,
	UINT tessU, UINT tessV,
	std::vector<Float3> &grid,
	UINT &outSamplesU, UINT &outSamplesV)
{
	int order = GetBasisOrder(basis, degree);
	BasisFunc1D basisFunc = GetBasisFunction(basis, degree);

	if (cpWidth < (UINT)order || cpHeight < (UINT)order)
		return false;

	// Number of sub-patches in each direction
	int patchesU = cpWidth - order + 1;
	int patchesV = cpHeight - order + 1;

	// Tessellation segments per sub-patch
	UINT segsPerPatchU = std::max(1u, tessU / (UINT)patchesU);
	UINT segsPerPatchV = std::max(1u, tessV / (UINT)patchesV);

	// Total grid of sample points
	outSamplesU = segsPerPatchU * patchesU + 1;
	outSamplesV = segsPerPatchV * patchesV + 1;

	// Evaluate all sample points
	grid.resize(outSamplesU * outSamplesV);

	for (int pv = 0; pv < patchesV; pv++) {
		for (int pu = 0; pu < patchesU; pu++) {
			for (UINT sv = 0; sv <= segsPerPatchV; sv++) {
				// Skip last row except for last sub-patch
				if (sv == segsPerPatchV && pv != patchesV - 1)
					continue;
				for (UINT su = 0; su <= segsPerPatchU; su++) {
					// Skip last column except for last sub-patch
					if (su == segsPerPatchU && pu != patchesU - 1)
						continue;

					float u = (float)su / (float)segsPerPatchU;
					float v = (float)sv / (float)segsPerPatchV;

					Float3 pos = { 0, 0, 0 };
					for (int j = 0; j < order; j++) {
						float bv = basisFunc(j, v);
						for (int i = 0; i < order; i++) {
							float bu = basisFunc(i, u);
							float w = bu * bv;
							UINT cpIndex = (startV + pv + j) * cpStride + (startU + pu + i);
							Float3 cp = ReadVertexPosition(pVertexData, stride, cpIndex);
							pos.x += w * cp.x;
							pos.y += w * cp.y;
							pos.z += w * cp.z;
						}
					}

					UINT gridX = pu * segsPerPatchU + su;
					UINT gridY = pv * segsPerPatchV + sv;
					grid[gridY * outSamplesU + gridX] = pos;
				}
			}
		}
	}

	return true;
}

// ---------------------------------------------------------------
// Compute surface normal at a grid point using finite differences
// ---------------------------------------------------------------
static Float3 ComputeGridNormal(const std::vector<Float3> &grid, UINT samplesU, UINT samplesV, UINT x, UINT y)
{
	// dP/du via central (or forward/backward at edges) differences
	Float3 dPdu;
	if (x > 0 && x < samplesU - 1) {
		const Float3 &left  = grid[y * samplesU + (x - 1)];
		const Float3 &right = grid[y * samplesU + (x + 1)];
		dPdu = { right.x - left.x, right.y - left.y, right.z - left.z };
	} else if (x == 0) {
		const Float3 &here  = grid[y * samplesU + x];
		const Float3 &right = grid[y * samplesU + (x + 1)];
		dPdu = { right.x - here.x, right.y - here.y, right.z - here.z };
	} else {
		const Float3 &left = grid[y * samplesU + (x - 1)];
		const Float3 &here = grid[y * samplesU + x];
		dPdu = { here.x - left.x, here.y - left.y, here.z - left.z };
	}
	// dP/dv via central (or forward/backward at edges) differences
	Float3 dPdv;
	if (y > 0 && y < samplesV - 1) {
		const Float3 &up   = grid[(y - 1) * samplesU + x];
		const Float3 &down = grid[(y + 1) * samplesU + x];
		dPdv = { down.x - up.x, down.y - up.y, down.z - up.z };
	} else if (y == 0) {
		const Float3 &here = grid[y * samplesU + x];
		const Float3 &down = grid[(y + 1) * samplesU + x];
		dPdv = { down.x - here.x, down.y - here.y, down.z - here.z };
	} else {
		const Float3 &up   = grid[(y - 1) * samplesU + x];
		const Float3 &here = grid[y * samplesU + x];
		dPdv = { here.x - up.x, here.y - up.y, here.z - up.z };
	}
	return Float3Normalize(Float3Cross(dPdu, dPdv));
}

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
// Tessellate a triangular patch into a triangle list
//
// Uses barycentric coordinates. For a cubic Bezier triangle,
// there are 10 control points indexed by (i,j,k) where i+j+k=3.
// We uniformly sample the barycentric domain.
// ---------------------------------------------------------------
static bool TessellateTriPatch(
	const uint8_t *pVertexData,
	UINT stride,
	UINT startVertex,
	UINT numVertices,
	X_D3DBASISTYPE basis,
	X_D3DDEGREETYPE degree,
	UINT tessSegs,
	std::vector<Float3> &outVertices)
{
	// For simplicity, support cubic Bezier triangle (10 control points)
	if (degree != D3DDEGREE_LINEAR && degree != D3DDEGREE_CUBIC)
		return false;

	UINT gridSize = tessSegs + 1;

	// Evaluate the patch at barycentric sample points
	std::vector<Float3> grid; // indexed by (row, col) within the triangular grid
	grid.resize(gridSize * (gridSize + 1) / 2);

	UINT idx = 0;
	for (UINT row = 0; row <= tessSegs; row++) {
		UINT colCount = tessSegs - row + 1;
		for (UINT col = 0; col < colCount; col++) {
			float u = (float)col / (float)tessSegs;
			float v = (float)row / (float)tessSegs;
			float w = 1.0f - u - v;

			if (degree == D3DDEGREE_LINEAR) {
				// Linear interpolation of 3 control points
				if (numVertices < 3) return false;
				Float3 p0 = ReadVertexPosition(pVertexData, stride, startVertex + 0);
				Float3 p1 = ReadVertexPosition(pVertexData, stride, startVertex + 1);
				Float3 p2 = ReadVertexPosition(pVertexData, stride, startVertex + 2);
				grid[idx] = {
					w * p0.x + u * p1.x + v * p2.x,
					w * p0.y + u * p1.y + v * p2.y,
					w * p0.z + u * p1.z + v * p2.z
				};
			} else {
				// Cubic Bezier triangle: 10 control points
				// Indexed by (i,j,k) where i+j+k = 3
				// Standard ordering: (3,0,0),(2,1,0),(1,2,0),(0,3,0),
				//                    (2,0,1),(1,1,1),(0,2,1),
				//                    (1,0,2),(0,1,2),
				//                    (0,0,3)
				if (numVertices < 10) return false;

				// Bernstein polynomial for triangular patch:
				// B(u,v,w) = sum over i+j+k=n of (n!/(i!j!k!)) * u^i * v^j * w^k * P_ijk
				float u2 = u * u, u3 = u2 * u;
				float v2 = v * v, v3 = v2 * v;
				float w2 = w * w, w3 = w2 * w;

				// Bernstein basis polynomials for degree 3
				float b300 = w3;
				float b210 = 3.0f * w2 * u;
				float b120 = 3.0f * w * u2;
				float b030 = u3;
				float b201 = 3.0f * w2 * v;
				float b111 = 6.0f * w * u * v;
				float b021 = 3.0f * u2 * v;
				float b102 = 3.0f * w * v2;
				float b012 = 3.0f * u * v2;
				float b003 = v3;

				// Read 10 control points in standard ordering
				Float3 cp[10];
				for (int c = 0; c < 10; c++)
					cp[c] = ReadVertexPosition(pVertexData, stride, startVertex + c);

				grid[idx] = {
					b300*cp[0].x + b210*cp[1].x + b120*cp[2].x + b030*cp[3].x +
					b201*cp[4].x + b111*cp[5].x + b021*cp[6].x +
					b102*cp[7].x + b012*cp[8].x + b003*cp[9].x,
					b300*cp[0].y + b210*cp[1].y + b120*cp[2].y + b030*cp[3].y +
					b201*cp[4].y + b111*cp[5].y + b021*cp[6].y +
					b102*cp[7].y + b012*cp[8].y + b003*cp[9].y,
					b300*cp[0].z + b210*cp[1].z + b120*cp[2].z + b030*cp[3].z +
					b201*cp[4].z + b111*cp[5].z + b021*cp[6].z +
					b102*cp[7].z + b012*cp[8].z + b003*cp[9].z
				};
			}
			idx++;
		}
	}

	// Build triangle list from the triangular grid
	outVertices.clear();
	outVertices.reserve(tessSegs * tessSegs * 3);

	// Helper to index into the triangular grid: row r, column c
	auto GridIndex = [tessSegs](UINT r, UINT c) -> UINT {
		// Sum of (tessSegs - i + 1) for i = 0..r-1 = r*(tessSegs+1) - r*(r-1)/2
		return r * (tessSegs + 1) - r * (r - 1) / 2 + c;
	};

	for (UINT row = 0; row < tessSegs; row++) {
		UINT colCount = tessSegs - row;
		for (UINT col = 0; col < colCount; col++) {
			// Upward triangle
			outVertices.push_back(grid[GridIndex(row, col)]);
			outVertices.push_back(grid[GridIndex(row, col + 1)]);
			outVertices.push_back(grid[GridIndex(row + 1, col)]);

			// Downward triangle (if not at the edge)
			if (col + 1 < colCount) {
				outVertices.push_back(grid[GridIndex(row, col + 1)]);
				outVertices.push_back(grid[GridIndex(row + 1, col + 1)]);
				outVertices.push_back(grid[GridIndex(row + 1, col)]);
			}
		}
	}

	return true;
}

// ---------------------------------------------------------------
// Build multi-attribute triangle list from a position grid.
// Outputs Float4 per register per vertex: position at register 0,
// AUTONORMAL and AUTOTEXCOORD at their respective registers,
// other declared registers at default (0,0,0,1).
//
// outData is a flat array of Float4; each vertex occupies
// numRegs consecutive Float4 values.
// ---------------------------------------------------------------
static void BuildMultiAttributeTriangleList(
	const std::vector<Float3> &grid, UINT samplesU, UINT samplesV,
	int autoNormalReg, int autoTexcoordReg,
	int numRegs, const int *regIndices,
	std::vector<Float4> &outData, UINT &outVertexCount)
{
	UINT numQuads = (samplesU - 1) * (samplesV - 1);
	outVertexCount = numQuads * 6;
	outData.resize((size_t)outVertexCount * numRegs);

	// Map register index to output slot for fast lookup
	int regSlot[X_VSH_MAX_ATTRIBUTES];
	for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) regSlot[i] = -1;
	for (int s = 0; s < numRegs; s++) regSlot[regIndices[s]] = s;

	auto EmitVertex = [&](UINT vIdx, UINT gx, UINT gy) {
		Float4 *pVert = &outData[(size_t)vIdx * numRegs];
		// Initialize all slots to default (0,0,0,1)
		for (int s = 0; s < numRegs; s++)
			pVert[s] = { 0.0f, 0.0f, 0.0f, 1.0f };

		// Position at register 0 (X_D3DVSDE_POSITION)
		const Float3 &pos = grid[gy * samplesU + gx];
		if (regSlot[0] >= 0)
			pVert[regSlot[0]] = { pos.x, pos.y, pos.z, 1.0f };

		// AUTONORMAL: surface normal from finite differences
		if (autoNormalReg >= 0 && regSlot[autoNormalReg] >= 0) {
			Float3 n = ComputeGridNormal(grid, samplesU, samplesV, gx, gy);
			pVert[regSlot[autoNormalReg]] = { n.x, n.y, n.z, 0.0f };
		}

		// AUTOTEXCOORD: parametric UV
		if (autoTexcoordReg >= 0 && regSlot[autoTexcoordReg] >= 0) {
			float u = (samplesU > 1) ? (float)gx / (float)(samplesU - 1) : 0.0f;
			float v = (samplesV > 1) ? (float)gy / (float)(samplesV - 1) : 0.0f;
			pVert[regSlot[autoTexcoordReg]] = { u, v, 0.0f, 1.0f };
		}
	};

	UINT vIdx = 0;
	for (UINT y = 0; y < samplesV - 1; y++) {
		for (UINT x = 0; x < samplesU - 1; x++) {
			EmitVertex(vIdx++, x,     y);
			EmitVertex(vIdx++, x + 1, y);
			EmitVertex(vIdx++, x,     y + 1);
			EmitVertex(vIdx++, x + 1, y);
			EmitVertex(vIdx++, x + 1, y + 1);
			EmitVertex(vIdx++, x,     y + 1);
		}
	}
}

// ---------------------------------------------------------------
// Cached tessellation input layout (one per VS bytecode pointer)
// ---------------------------------------------------------------
static ID3D11InputLayout *s_pTessInputLayout = nullptr;
static const void *s_TessInputLayoutVSPtr = nullptr;
static int s_TessInputLayoutNumRegs = 0;

static ID3D11InputLayout* GetOrCreateTessInputLayout(
	const int *regIndices, int numRegs, ID3DBlob *pVSBytecode)
{
	if (pVSBytecode == nullptr || numRegs == 0)
		return nullptr;

	// Invalidate cache if VS bytecode changed or register count changed
	if (pVSBytecode->GetBufferPointer() != s_TessInputLayoutVSPtr || numRegs != s_TessInputLayoutNumRegs) {
		if (s_pTessInputLayout) { s_pTessInputLayout->Release(); s_pTessInputLayout = nullptr; }

		std::vector<D3D11_INPUT_ELEMENT_DESC> elements(numRegs);
		for (int i = 0; i < numRegs; i++) {
			elements[i].SemanticName = "TEXCOORD";
			elements[i].SemanticIndex = (UINT)regIndices[i];
			elements[i].Format = DXGI_FORMAT_R32G32B32A32_FLOAT;
			elements[i].InputSlot = 0;
			elements[i].AlignedByteOffset = (UINT)(i * sizeof(Float4));
			elements[i].InputSlotClass = D3D11_INPUT_PER_VERTEX_DATA;
			elements[i].InstanceDataStepRate = 0;
		}

		// Filter out elements whose semantics the shader doesn't declare
		auto filtered = FilterInputElementsByShaderSignature(
			elements.data(), (UINT)elements.size(),
			pVSBytecode->GetBufferPointer(), pVSBytecode->GetBufferSize());

		HRESULT hr = E_FAIL;
		if (!filtered.empty()) {
			hr = g_pD3DDevice->CreateInputLayout(
				filtered.data(), (UINT)filtered.size(),
				pVSBytecode->GetBufferPointer(), pVSBytecode->GetBufferSize(),
				&s_pTessInputLayout);
		}
		if (FAILED(hr)) {
			EmuLog(LOG_LEVEL::WARNING, "GetOrCreateTessInputLayout: CreateInputLayout failed (0x%08X)", hr);
			s_pTessInputLayout = nullptr;
		}

		s_TessInputLayoutVSPtr = pVSBytecode->GetBufferPointer();
		s_TessInputLayoutNumRegs = numRegs;
	}
	return s_pTessInputLayout;
}

// ---------------------------------------------------------------
// Draw tessellated mesh (position-only fast path)
// ---------------------------------------------------------------
static HRESULT DrawTessellatedMesh(const std::vector<Float3> &vertices)
{
	if (vertices.empty())
		return E_FAIL;

	UINT vertexCount = (UINT)vertices.size();
	UINT dataSize = vertexCount * sizeof(Float3);

	static CxbxDynBuffer s_PatchVB = { nullptr, 0, D3D11_BIND_VERTEX_BUFFER };
	ID3D11Buffer *pVB = s_PatchVB.Update(vertices.data(), dataSize);
	if (pVB == nullptr)
		return E_FAIL;

	// Apply any pending state changes
	CxbxD3D11ApplyDirtyStates();

	UINT stride = sizeof(Float3);
	UINT offset = 0;
	g_pD3DDeviceContext->IASetVertexBuffers(0, 1, &pVB, &stride, &offset);
	g_pD3DDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	g_pD3DDeviceContext->Draw(vertexCount, 0);

	return S_OK;
}

// ---------------------------------------------------------------
// Draw tessellated mesh (multi-attribute path with auto-generated
// normals/texcoords and a dedicated input layout)
// ---------------------------------------------------------------
static HRESULT DrawTessellatedMeshMultiAttr(
	const std::vector<Float4> &vertexData, UINT vertexCount,
	int numRegs, const int *regIndices)
{
	if (vertexCount == 0)
		return E_FAIL;

	UINT dataSize = (UINT)(vertexData.size() * sizeof(Float4));

	static CxbxDynBuffer s_PatchVB = { nullptr, 0, D3D11_BIND_VERTEX_BUFFER };
	ID3D11Buffer *pVB = s_PatchVB.Update(vertexData.data(), dataSize);
	if (pVB == nullptr)
		return E_FAIL;

	// Apply pending state changes (shaders, constants, render targets, etc.)
	CxbxD3D11ApplyDirtyStates();

	// Override input layout with tessellation-specific one
	ID3DBlob *pVSBytecode = CxbxGetActiveVertexShaderBytecode();
	ID3D11InputLayout *pTessLayout = GetOrCreateTessInputLayout(regIndices, numRegs, pVSBytecode);
	if (pTessLayout != nullptr)
		g_pD3DDeviceContext->IASetInputLayout(pTessLayout);

	UINT stride = (UINT)(numRegs * sizeof(Float4));
	UINT offset = 0;
	g_pD3DDeviceContext->IASetVertexBuffers(0, 1, &pVB, &stride, &offset);
	g_pD3DDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	g_pD3DDeviceContext->Draw(vertexCount, 0);

	return S_OK;
}

// ---------------------------------------------------------------
// Determine whether auto-generation is needed, collect register
// list, and choose the appropriate draw path.
// ---------------------------------------------------------------
static HRESULT DrawRectPatchWithAutoGen(
	const std::vector<Float3> &grid, UINT samplesU, UINT samplesV)
{
	CxbxVertexDeclaration *pDecl = CxbxGetVertexDeclaration();
	bool needAutoGen = pDecl != nullptr &&
		(pDecl->autoNormalRegister >= 0 || pDecl->autoTexcoordRegister >= 0);

	if (!needAutoGen) {
		// Fast path: position-only output
		std::vector<Float3> triList;
		BuildTriangleListFromGrid(grid, samplesU, samplesV, triList);
		return DrawTessellatedMesh(triList);
	}

	// Multi-attribute path: collect all needed registers
	// Start with registers declared in the input layout
	int regIndices[X_VSH_MAX_ATTRIBUTES];
	int numRegs = 0;
	bool regIncluded[X_VSH_MAX_ATTRIBUTES] = {};

	if (pDecl->pD3D11InputElements != nullptr) {
		for (UINT i = 0; i < pDecl->D3D11InputElementCount; i++) {
			int reg = (int)pDecl->pD3D11InputElements[i].SemanticIndex;
			if (reg >= 0 && reg < X_VSH_MAX_ATTRIBUTES && !regIncluded[reg]) {
				regIndices[numRegs++] = reg;
				regIncluded[reg] = true;
			}
		}
	}

	// Add auto-generated registers (they're not in the input layout)
	if (pDecl->autoNormalRegister >= 0 && !regIncluded[pDecl->autoNormalRegister]) {
		regIndices[numRegs++] = pDecl->autoNormalRegister;
		regIncluded[pDecl->autoNormalRegister] = true;
	}
	if (pDecl->autoTexcoordRegister >= 0 && !regIncluded[pDecl->autoTexcoordRegister]) {
		regIndices[numRegs++] = pDecl->autoTexcoordRegister;
		regIncluded[pDecl->autoTexcoordRegister] = true;
	}

	// Sort by register index for consistent layout
	std::sort(regIndices, regIndices + numRegs);

	// Build multi-attribute triangle list
	std::vector<Float4> vertexData;
	UINT vertexCount = 0;
	BuildMultiAttributeTriangleList(
		grid, samplesU, samplesV,
		pDecl->autoNormalRegister, pDecl->autoTexcoordRegister,
		numRegs, regIndices,
		vertexData, vertexCount);

	return DrawTessellatedMeshMultiAttr(vertexData, vertexCount, numRegs, regIndices);
}

// ---------------------------------------------------------------
// Public API: Tessellate and draw a rectangular patch
// ---------------------------------------------------------------
HRESULT CxbxDrawRectPatchD3D11(
	UINT Handle,
	const float *pNumSegs,
	const xbox::X_D3DRECTPATCH_INFO *pRectPatchInfo)
{
	(void)Handle; // Not used in D3D11 path (no hardware tessellation cache)
	if (pRectPatchInfo == nullptr) {
		EmuLog(LOG_LEVEL::WARNING, "DrawRectPatch: pRectPatchInfo is NULL");
		return E_INVALIDARG;
	}

	// Get tessellation level
	float tessU = 8.0f, tessV = 8.0f;
	if (pNumSegs != nullptr) {
		tessU = pNumSegs[0];
		tessV = pNumSegs[1];
	}
	if (tessU < 1.0f) tessU = 1.0f;
	if (tessV < 1.0f) tessV = 1.0f;

	// Get vertex data from stream 0
	auto &stream = g_Xbox_SetStreamSource[0];
	if (stream.VertexBuffer == xbox::zeroptr) {
		EmuLog(LOG_LEVEL::WARNING, "DrawRectPatch: No vertex buffer bound to stream 0");
		return E_FAIL;
	}

	const uint8_t *pVertexData = (const uint8_t *)GetDataFromXboxResource(stream.VertexBuffer);
	if (pVertexData == nullptr) {
		EmuLog(LOG_LEVEL::WARNING, "DrawRectPatch: Could not get vertex data");
		return E_FAIL;
	}

	UINT stride = stream.Stride;
	if (stride == 0) stride = sizeof(Float3); // fallback

	std::vector<Float3> grid;
	UINT samplesU = 0, samplesV = 0;
	bool ok = EvaluateRectPatchGrid(
		pVertexData, stride,
		pRectPatchInfo->StartVertexOffsetWidth,
		pRectPatchInfo->StartVertexOffsetHeight,
		pRectPatchInfo->Width,
		pRectPatchInfo->Height,
		pRectPatchInfo->Stride,
		pRectPatchInfo->Basis,
		pRectPatchInfo->Degree,
		(UINT)tessU, (UINT)tessV,
		grid, samplesU, samplesV
	);

	if (!ok) {
		EmuLog(LOG_LEVEL::WARNING, "DrawRectPatch: Tessellation failed (Basis=%d, Degree=%d, %dx%d)",
			pRectPatchInfo->Basis, pRectPatchInfo->Degree,
			pRectPatchInfo->Width, pRectPatchInfo->Height);
		return E_FAIL;
	}

	return DrawRectPatchWithAutoGen(grid, samplesU, samplesV);
}

// ---------------------------------------------------------------
// Public API: Tessellate and draw a triangular patch
// ---------------------------------------------------------------
HRESULT CxbxDrawTriPatchD3D11(
	UINT Handle,
	const float *pNumSegs,
	const xbox::X_D3DTRIPATCH_INFO *pTriPatchInfo)
{
	(void)Handle; // Not used in D3D11 path (no hardware tessellation cache)
	if (pTriPatchInfo == nullptr) {
		EmuLog(LOG_LEVEL::WARNING, "DrawTriPatch: pTriPatchInfo is NULL");
		return E_INVALIDARG;
	}

	// Get tessellation level — for tri patches, all three edge segments
	// are specified but we use a uniform count
	float tessSegs = 8.0f;
	if (pNumSegs != nullptr) {
		tessSegs = std::max({ pNumSegs[0], pNumSegs[1], pNumSegs[2] });
	}
	if (tessSegs < 1.0f) tessSegs = 1.0f;

	// Get vertex data from stream 0
	auto &stream = g_Xbox_SetStreamSource[0];
	if (stream.VertexBuffer == xbox::zeroptr) {
		EmuLog(LOG_LEVEL::WARNING, "DrawTriPatch: No vertex buffer bound to stream 0");
		return E_FAIL;
	}

	const uint8_t *pVertexData = (const uint8_t *)GetDataFromXboxResource(stream.VertexBuffer);
	if (pVertexData == nullptr) {
		EmuLog(LOG_LEVEL::WARNING, "DrawTriPatch: Could not get vertex data");
		return E_FAIL;
	}

	UINT stride = stream.Stride;
	if (stride == 0) stride = sizeof(Float3);

	std::vector<Float3> tessellated;
	bool ok = TessellateTriPatch(
		pVertexData, stride,
		pTriPatchInfo->StartVertexOffset,
		pTriPatchInfo->NumVertices,
		pTriPatchInfo->Basis,
		pTriPatchInfo->Degree,
		(UINT)tessSegs,
		tessellated
	);

	if (!ok) {
		EmuLog(LOG_LEVEL::WARNING, "DrawTriPatch: Tessellation failed (Basis=%d, Degree=%d, NumVerts=%d)",
			pTriPatchInfo->Basis, pTriPatchInfo->Degree, pTriPatchInfo->NumVertices);
		return E_FAIL;
	}

	return DrawTessellatedMesh(tessellated);
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
//     their first float4 doesn't match strip col0 values — likely a different encoding)
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
	// Each strip curve provides posUOrder forward-difference coefficients [f, Δf, Δ²f, ...]
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
