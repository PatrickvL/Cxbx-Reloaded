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
// *  (c) 2002-2004 Aaron Robinson <caustik@caustik.com>
// *                Kingofc <kingofc@freenet.de>
// *
// *  All rights reserved
// *
// ******************************************************************
#define LOG_PREFIX CXBXR_MODULE::VTXSH

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h" // For CxbxPageTrackerGetMirrorSRV
#include "Backend_D3D11_Profiler.h"
#include "Shading\Shader.h" // For LoadPrecompiledCSO
#include "Shading\VertexShaderCache.h"

#include "core\hle\D3D8\XbVertexShader.h"
#include "core\hle\D3D8\XbPushBuffer.h" // For g_NV2A
#include "core\hle\D3D8\Rendering\NV2A_PGRAPH_Helpers.h"
#include "core\hle\D3D8\XbD3D8Logging.h" // For DEBUG_D3DRESULT
#include "devices\xbox.h"
#include "core\hle\D3D8\XbConvert.h" // For NV2A_VP_UPLOAD_INST
#include "devices\video\nv2a.h" // For D3DPUSH_DECODE
#include "common\Logging.h" // For LOG_INIT
#include "common\Settings.hpp" // for g_LibVersion_D3D8

#include "nv2a_vsh_emulator.h"

// Retained bytecode for FixedFunction vertex shader (needed for input layout creation)
static ID3DBlob* g_pD3D11FixedFunctionBytecode = nullptr;
static ID3DBlob* g_pD3D11JITVSBytecode = nullptr; // JIT-compiled VS bytecode for current draw
static ID3D11VertexShader* g_pD3D11JITCurrentVS = nullptr; // Currently active JIT VS

extern ID3D11VertexShader* CxbxCreateVertexShader(ID3DBlob* pCompiledShader, const char *shader_category)
{
	ID3D11VertexShader* pHostVertexShader = nullptr;

	if (g_pD3DDevice == nullptr) {
		EmuLog(LOG_LEVEL::WARNING, "Can't create %s vertex shader - no D3D device is set!", shader_category);
	}
	else {
		assert(pCompiledShader);

		HRESULT hRet;
		hRet = g_pD3DDevice->CreateVertexShader(
			(const void*)pCompiledShader->GetBufferPointer(),
			pCompiledShader->GetBufferSize(),
			nullptr,
			&pHostVertexShader
		);
		if (FAILED(hRet)) CxbxrAbort("Failed to create %s vertex shader", shader_category);
	}

	return pHostVertexShader;
}

ID3D11VertexShader* InitShader(const char* csoName, const char* label, ID3DBlob** ppRetainedBytecode = nullptr) {
	ID3D11VertexShader* shader = nullptr;

	ID3DBlob* pBlob = nullptr;
	LoadPrecompiledCSO(csoName, &pBlob);
	if (pBlob) {
		shader = CxbxCreateVertexShader(pBlob, label);
		if (ppRetainedBytecode) {
			*ppRetainedBytecode = pBlob; // Caller takes ownership
		} else {
			pBlob->Release();
		}
	}

	return shader;
}

// Upload NV2A XFPR (Transform Program RAM) and bind SRVs for the VS interpreter.
//
// Two StructuredBuffers feed the interpreter shader:
//   - g_PGRegs (t12): shared PGRAPH register array - already uploaded by
//     CxbxD3D11UploadRCInterpreterState(). We just bind it to the VS stage.
//   - g_XFPR (t5): pg->xf.xfpr[136][4] - the XFPR RAM mirror,
//     uploaded here.  On real NV2A hardware this is on-chip XF SRAM
//     behind the RDI interface, uploaded via NV097_SET_TRANSFORM_PROGRAM
//     with an auto-incrementing write pointer (CHEOPS_OFFSET.PROG_LD_PTR).
//
// The shader reads CHEOPS_PROGRAM_START from g_PGRegs to find the first
// active instruction slot and loops until FLD_FINAL.
void CxbxD3D11UploadVSInterpreterState(const xbox::dword_xt* /*pXboxMicrocode*/)
{
	if (!g_pD3D11XFPRBuf)
		return;

	// PGRAPH source: upload the entire program_data[] array (XFPR mirror).
	// The shader selects the active program via CHEOPS_PROGRAM_START.
	PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;

	// Skip XFPR upload if program data hasn't changed (dirty flag set by NV097_SET_TRANSFORM_PROGRAM)
	static uint32_t s_LastProgramGen = ~0u;
	if (pg->dirty[NV2A_DIRTY_PROGRAM] != s_LastProgramGen) {
		CxbxD3D11UpdateDynamicBuffer(g_pD3D11XFPRBuf,
			pg->xf.xfpr, sizeof(pg->xf.xfpr));
		s_LastProgramGen = pg->dirty[NV2A_DIRTY_PROGRAM];
	}

	// Bind VS interpreter SRVs once - pointers are stable for device lifetime
	static bool s_VSInterpreterSRVsBound = false;
	if (!s_VSInterpreterSRVsBound) {
		ID3D11ShaderResourceView* pMirrorSRV = CxbxPageTrackerGetMirrorSRV();
		g_pD3DDeviceContext->VSSetShaderResources(CXBX_D3D11_VS_PGREGS_SRV_SLOT, 1, &pMirrorSRV);
		g_pD3DDeviceContext->VSSetShaderResources(CXBX_D3D11_VS_XFPR_SRV_SLOT, 1, &g_pD3D11XFPRSRV);
		s_VSInterpreterSRVsBound = true;
	}
}

void CxbxUpdateHostVertexShader()
{
	// Vertex shaders are loaded once from embedded precompiled blobs (CSOs).
	// They persist for the lifetime of the D3D11 device; teardown is handled
	// by CxbxD3D11ReleaseBackendResources() on device release.
	static ID3D11VertexShader* fixedFunctionShader = nullptr;
	static bool shadersLoaded = false;

	if (!shadersLoaded) {
		shadersLoaded = true;
		CxbxSetVertexShader(nullptr);

		EmuLog(LOG_LEVEL::INFO, "Loading vertex shaders...");
		fixedFunctionShader = InitShader("CxbxFixedFunctionVS", "Fixed Function Vertex Shader", &g_pD3D11FixedFunctionBytecode);
		// VS interpreter is initialized lazily on first ShaderProgram draw via CxbxD3D11InitVSInterpreter()
	}

	// Select the active vertex shader based on current PGRAPH mode.
	// Called every draw; the actual shader objects are already loaded above.

	LOG_INIT; // Allows use of DEBUG_D3DRESULT

	PGRAPHState *pg = &g_NV2A->GetDeviceState()->pgraph;

	if (NV2AIsFixedFunctionMode(pg)) {
		HRESULT hRet = CxbxSetVertexShader(fixedFunctionShader);
		if (FAILED(hRet)) CxbxrAbort("Failed to set fixed-function shader");
	}
	else {
		// Read program tokens from PGRAPH program_data (the authoritative
		// source, written by the PFIFO puller from NV2A_VP_UPLOAD_INST).
		// The start address comes from CSV0_C CHEOPS_PROGRAM_START, which
		// the puller sets from NV097_SET_TRANSFORM_PROGRAM_START.

		// Skip JIT cache lookup when program hasn't changed since last draw
		static uint32_t s_lastProgramGen = UINT32_MAX;
		static ID3D11VertexShader* s_lastJITVS = nullptr;
		if (pg->dirty[NV2A_DIRTY_PROGRAM] == s_lastProgramGen && s_lastJITVS) {
			CxbxSetVertexShader(s_lastJITVS);
			return;
		}

		xbox::dword_xt *pTokens = nullptr;
		uint32_t startAddr = GET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_C)],
			NV_PGRAPH_CSV0_C_CHEOPS_PROGRAM_START);
		if (startAddr < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH) {
			pTokens = (xbox::dword_xt*)&pg->xf.xfpr[startAddr][0];
		}
		if (!pTokens) {
			LOG_TEST_CASE("PGRAPH program_data not available");
			return;
		}

		// Try JIT compilation first (10-100x faster than interpreter)
		{
			CXBX_PROFILE_SCOPE(PROF_VS_SHADER);
			ID3DBlob* pJITBytecode = nullptr;
			ID3D11VertexShader* pJITVS = g_VertexShaderCache.GetShader(
				pg->xf.xfpr, startAddr, g_pD3DDevice, &pJITBytecode);
			if (pJITVS) {
				InterlockedIncrement(&g_ProfileVSJITHits);
				// Release previous JIT bytecode ref
				if (g_pD3D11JITVSBytecode) { g_pD3D11JITVSBytecode->Release(); g_pD3D11JITVSBytecode = nullptr; }
				g_pD3D11JITVSBytecode = pJITBytecode;
				g_pD3D11JITCurrentVS = pJITVS;
				HRESULT hRet = CxbxSetVertexShader(pJITVS);
				DEBUG_D3DRESULT(hRet, "CxbxSetVertexShader(JIT)");

				s_lastProgramGen = pg->dirty[NV2A_DIRTY_PROGRAM];
				s_lastJITVS = pJITVS;
				return; // Skip interpreter path
			}
		}

		s_lastProgramGen = pg->dirty[NV2A_DIRTY_PROGRAM];
		s_lastJITVS = nullptr; // interpreter path — don't cache

		if (g_bUseVSInterpreter && CxbxD3D11InitVSInterpreter()) {
			InterlockedIncrement(&g_ProfileVSInterpreterHits);
			// Upload the raw NV2A microcode to the interpreter constant buffer
			CxbxD3D11UploadVSInterpreterState(pTokens);
			HRESULT hRet = CxbxSetVertexShader(g_pD3D11VSInterpreterVS);
			DEBUG_D3DRESULT(hRet, "CxbxSetVertexShader(VSInterpreter)");
		}
	}
}

ID3DBlob* CxbxGetActiveVertexShaderBytecode()
{
	if (NV2AIsFixedFunctionMode())
		return g_pD3D11FixedFunctionBytecode;
	// JIT-compiled VS provides its own bytecode
	if (g_pD3D11JITVSBytecode)
		return g_pD3D11JITVSBytecode;
	// VS interpreter provides its own bytecode for input layout creation
	if (g_bUseVSInterpreter && g_pD3D11VSInterpreterBytecode)
		return g_pD3D11VSInterpreterBytecode;
	return nullptr;
}

ID3DBlob* CxbxGetFixedFunctionVertexShaderBytecode()
{
	return g_pD3D11FixedFunctionBytecode;
}

void CxbxUpdateHostVertexDeclaration()
{
	// Titles can specify default values for registers via calls like SetVertexData4f
	// HLSL shaders need to know whether to use vertex data or default vertex shader values
	// Any register not in the vertex declaration should be set to the default value
	float vertexDefaultFlags[X_VSH_MAX_ATTRIBUTES];

	// PGRAPH-driven: read which attributes are active from NV2A state
	// Input layout is set to nullptr since the vertex pull CS handles all attribute fetching.
	g_pD3DDeviceContext->IASetInputLayout(nullptr);
	PGRAPHState* pg = (g_NV2A != nullptr) ? &g_NV2A->GetDeviceState()->pgraph : nullptr;
	for (int i = 0; i < X_VSH_MAX_ATTRIBUTES; i++) {
		bool active = pg && (pg->vertex_attributes[i].count > 0);
		vertexDefaultFlags[i] = active ? 0.0f : 1.0f;
	}

	// Only upload if the flags changed since last draw
	static float s_CachedVertexDefaultFlags[X_VSH_MAX_ATTRIBUTES] = {};
	static bool s_FirstDefaultFlagsCall = true;
	if (s_FirstDefaultFlagsCall || std::memcmp(vertexDefaultFlags, s_CachedVertexDefaultFlags, sizeof(vertexDefaultFlags)) != 0) {
		std::memcpy(s_CachedVertexDefaultFlags, vertexDefaultFlags, sizeof(vertexDefaultFlags));
		s_FirstDefaultFlagsCall = false;
		CxbxSetVertexShaderConstantF(CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_BASE, vertexDefaultFlags, CXBX_D3DVS_CONSTREG_VREGDEFAULTS_FLAG_SIZE);
	}
}

void D3D11_launch_transform_program(NV2AState *d, unsigned int program_start)
{
	PGRAPHState* pg = &(d->pgraph);

	// Cache the parsed program globally; only re-parse when program_data changes
	static Nv2aVshProgram s_CachedProgram = {};
	static bool s_CachedProgramValid = false;
	static uint32_t s_LastXFPRGen = ~0u;

	if (pg->dirty[NV2A_DIRTY_PROGRAM] != s_LastXFPRGen || !s_CachedProgramValid) {
		if (s_CachedProgramValid) {
			nv2a_vsh_program_destroy(&s_CachedProgram);
		}
		s_CachedProgram = {};
		Nv2aVshParseResult result = nv2a_vsh_parse_program(
			&s_CachedProgram,
			pg->xf.xfpr[0],
			NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
		if (result != NV2AVPR_SUCCESS) {
			LOG_TEST_CASE("nv2a_vsh_parse_program failed (cached full parse)");
			s_CachedProgramValid = false;
			return;
		}
		// Guard against buffer overflow: force is_final on the last slot so
		// the executor always terminates within the 136-entry allocation,
		// even if no instruction in the program sets the final bit.
		s_CachedProgram.steps[NV2A_MAX_TRANSFORM_PROGRAM_LENGTH - 1].is_final = true;
		s_CachedProgramValid = true;
		s_LastXFPRGen = pg->dirty[NV2A_DIRTY_PROGRAM];
	}

	// Create a view into the cached program starting at program_start
	// Execution stops naturally at the step with is_final==true
	Nv2aVshProgram program;
	program.steps = s_CachedProgram.steps + program_start;

	Nv2aVshCPUXVSSExecutionState state_linkage;
	Nv2aVshExecutionState state = nv2a_vsh_emu_initialize_xss_execution_state(
		&state_linkage, (float*)pg->xf.xfctx);
	memcpy(state_linkage.input_regs, pg->xf.vertex_state_shader_v0,
		sizeof(pg->xf.vertex_state_shader_v0));

	// Local execute-and-track loop: same as nv2a_vsh_emu_execute_track_context_writes
	// but writes directly to the uint32_t bitmap instead of a bool array.
	{
		Nv2aVshStep *step = program.steps;
		while (true) {
			nv2a_vsh_emu_apply(&state, step);
			// Track context register writes from MAC unit
			if (step->mac.opcode) {
				if (step->mac.outputs[0].type == NV2ART_CONTEXT) {
					unsigned idx = step->mac.outputs[0].index;
					pg->xf.xfctx_dirty[idx / 32] |= (1u << (idx % 32));
					pg->xf.xfctx_generation++;
				} else if (step->mac.outputs[1].type == NV2ART_CONTEXT) {
					unsigned idx = step->mac.outputs[1].index;
					pg->xf.xfctx_dirty[idx / 32] |= (1u << (idx % 32));
					pg->xf.xfctx_generation++;
				}
			}
			// Track context register writes from ILU unit
			if (step->ilu.opcode) {
				if (step->ilu.outputs[0].type == NV2ART_CONTEXT) {
					unsigned idx = step->ilu.outputs[0].index;
					pg->xf.xfctx_dirty[idx / 32] |= (1u << (idx % 32));
					pg->xf.xfctx_generation++;
				} else if (step->ilu.outputs[1].type == NV2ART_CONTEXT) {
					unsigned idx = step->ilu.outputs[1].index;
					pg->xf.xfctx_dirty[idx / 32] |= (1u << (idx % 32));
					pg->xf.xfctx_generation++;
				}
			}
			if (step->is_final) break;
			++step;
		}
	}
	// Do NOT call nv2a_vsh_program_destroy here - program.steps is a borrowed pointer
}
