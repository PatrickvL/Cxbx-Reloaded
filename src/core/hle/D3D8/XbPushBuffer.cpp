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
// *  (c) 2002-2003 Aaron Robinson <caustik@caustik.com>
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::PSHB

#include <assert.h> // For assert()

#include "core\kernel\support\Emu.h"
#include "core\hle\D3D8\XbD3D8Types.h" // For X_D3DFORMAT
#include "core\hle\D3D8\ResourceTracker.h"
#include "core\hle\D3D8\Rendering\RenderGlobals.h" // For g_Xbox_VertexShader_Handle
#include "core\hle\D3D8\XbPushBuffer.h"
#include "core\hle\D3D8\XbConvert.h"
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11.h" // For CxbxD3D11IABypassDraw
#include "devices/video/nv2a.h" // For g_NV2A, PGRAPHState
#include "devices/video/nv2a_int.h" // For NV** defines
#include "Logging.h"

// TODO: Find somewhere to put this that doesn't conflict with xbox::
extern void CxbxUpdateHostTextures();

const char *NV2AMethodToString(DWORD dwMethod); // forward

void EmuExecutePushBuffer
(
	xbox::X_D3DPushBuffer       *pPushBuffer,
	xbox::X_D3DFixup            *pFixup
)
{
	//Check whether Fixup exists or not. 
	if (pFixup != xbox::zeroptr) {
		LOG_TEST_CASE("PushBuffer has fixups");
		//Interpret address of PushBuffer Data and Fixup Data
		UINT8* pPushBufferData = (UINT8*)pPushBuffer->Data;
		UINT8* pFixupData = (UINT8*)(pFixup->Data + pFixup->Run);
		UINT32 SizeInBytes = 0;
		UINT32 OffsetInBytes = 0;

		while (TRUE) {
			SizeInBytes = *(UINT32*)pFixupData;
			//If SizeInBytes==0xFFFFFFFF, end of Fixup Data.
			if (SizeInBytes == 0xFFFFFFFF)
				break;
			pFixupData += 4;
			OffsetInBytes = *(UINT32*)pFixupData;
			pFixupData += 4;
			//fixup must not exceed the pushbuffer data range.
			if ((OffsetInBytes + SizeInBytes) <= pPushBuffer->Size)
			{
				memcpy(pPushBufferData + OffsetInBytes, pFixupData, SizeInBytes);
			}
			pFixupData += SizeInBytes;
			/*
			When IDirect3DDevice8::RunPushBuffer is called with a fix-up object specified,
			it will parse the fix-up data pointed to by pFixup and with a byte offset of Run.
			The fix-up data is encoded as follows.The first DWORD is the size, in bytes,
			of the push-buffer fix-up to be modified.The second DWORD is the offset, in bytes,
			from the start of the push-buffer where the fix-up is to be modified.
			The subsequent DWORDS are the data to be copied. This encoding repeats for every fix-up to be done,
			until it terminates with a size value of 0xffffffff.
			The offsets must be in an increasing order.
			*/
		}
	}
    pfifo_submit_pushbuffer(g_NV2A->GetDeviceState(), (void*)pPushBuffer->Data, pPushBuffer->Size);

    return;
}

void HLE_draw_arrays(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	for (unsigned int i = 0; i < pg->draw_arrays_length; i++) {
		CxbxDrawContext DrawContext = {};

		DrawContext.XboxPrimitiveType = (xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode;
		DrawContext.dwStartVertex = pg->gl_draw_arrays_start[i];
		DrawContext.dwVertexCount = pg->gl_draw_arrays_count[i];

		CxbxD3D11IABypassDraw(DrawContext);
		g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(
			DrawContext.XboxPrimitiveType, DrawContext.dwVertexCount);
	}
}

void HLE_draw_inline_buffer(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	if (pg->inline_buffer_length == 0)
		return;

	CxbxD3D11DrawInlineBuffer(pg);

	g_dwPrimPerFrame += ConvertXboxVertexCountToPrimitiveCount(
		(xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode, pg->inline_buffer_length);
}

void HLE_draw_inline_array(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	// Compute per-vertex stride from NV2A vertex attribute format registers.
	// Inline array data packs all enabled attributes contiguously per vertex,
	// unlike array-based draws which use the stride field from the format register.
	unsigned int nv2a_stride = 0;
	for (int i = 0; i < NV2A_VERTEXSHADER_ATTRIBUTES; i++) {
		if (pg->vertex_attributes[i].count != 0) { // count 0 = disabled (format 0 is valid: UB_D3D/D3DCOLOR)
			nv2a_stride += pg->vertex_attributes[i].count * pg->vertex_attributes[i].size;
		}
	}

	if (nv2a_stride == 0 || pg->inline_array_length == 0) {
		return;
	}

	UINT VertexCount = (pg->inline_array_length * sizeof(DWORD)) / nv2a_stride;
	if (VertexCount == 0) {
		return;
	}

	CxbxDrawContext DrawContext = {};
	DrawContext.XboxPrimitiveType = (xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode;
	DrawContext.dwVertexCount = VertexCount;
	DrawContext.pXboxVertexStreamZeroData = pg->inline_array;
	DrawContext.uiXboxVertexStreamZeroStride = nv2a_stride;
	DrawContext.bNV2AInlineData = true;

	CxbxDrawPrimitiveUP(DrawContext);
}

void HLE_draw_inline_elements(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	unsigned int uiIndexCount = pg->inline_elements_length;
	CxbxDrawContext DrawContext = {};

		DrawContext.XboxPrimitiveType = (xbox::X_D3DPRIMITIVETYPE)pg->primitive_mode;
	DrawContext.dwVertexCount = uiIndexCount;
	DrawContext.pXboxIndexData = d->pgraph.inline_elements;

	CxbxDrawIndexed(DrawContext);
}

void HLE_draw_state_update(NV2AState *d)
{
	PGRAPHState *pg = &d->pgraph;

	// Vertex attribute inline_value may have changed via NV2A push buffer
	// (SET_VERTEX_DATA4F/4UB/2S) since the last draw. Mark defaults dirty
	// so the IA bypass re-uploads them before the next draw.
	g_bD3D11IABypassDefaultsDirty = true;

	CxbxUpdateNativeD3DResources();

	LOG_INCOMPLETE(); // TODO : Read state from pgraph, convert to D3D
}

void HLE_draw_clear(NV2AState *d)
{
	// PGRAPHState *pg = &d->pgraph;

	CxbxUpdateNativeD3DResources();

	LOG_INCOMPLETE(); // TODO : Read state from pgraph, convert to D3D (call EMUPATCH(D3DDevice_Clear)?)
}

// Import pgraph_draw_* variables, declared in EmuNV2A_PGRAPH.cpp :
extern void(*pgraph_draw_arrays)(NV2AState *d);
extern void(*pgraph_draw_inline_buffer)(NV2AState *d);
extern void(*pgraph_draw_inline_array)(NV2AState *d);
extern void(*pgraph_draw_inline_elements)(NV2AState *d);
extern void(*pgraph_draw_state_update)(NV2AState *d);
extern void(*pgraph_draw_clear)(NV2AState *d);

void HLE_init_pgraph_plugins()
{
	/* attach HLE Direct3D render plugins */
	pgraph_draw_arrays = HLE_draw_arrays;
	pgraph_draw_inline_buffer = HLE_draw_inline_buffer;
	pgraph_draw_inline_array = HLE_draw_inline_array;
	pgraph_draw_inline_elements = HLE_draw_inline_elements;
	pgraph_draw_state_update = HLE_draw_state_update;
	pgraph_draw_clear = HLE_draw_clear;
}

extern void pgraph_handle_method(
	NV2AState *d,
	unsigned int subchannel,
	unsigned int method,
	uint32_t parameter);

uint32_t HLE_read_NV2A_pgraph_register(const int reg)
{
	NV2AState* dev = g_NV2A->GetDeviceState();
	PGRAPHState *pg = &(dev->pgraph);
	return pg->regs[RI(reg)];
}

float *HLE_get_NV2A_vertex_attribute_value_pointer(unsigned slot)
{
	NV2AState* dev = g_NV2A->GetDeviceState();
	PGRAPHState *pg = &(dev->pgraph);

	// See CASE_16(NV097_SET_VERTEX_DATA4UB, 4) in LLE pgraph_handle_method()
	VertexAttribute *vertex_attribute = &pg->vertex_attributes[slot];
	return vertex_attribute->inline_value;
}

uint32_t HLE_read_NV2A_vertex_program_slot(unsigned program_load, unsigned slot)
{
	NV2AState* dev = g_NV2A->GetDeviceState();
	PGRAPHState* pg = &(dev->pgraph);

	// See CASE_32(NV097_SET_TRANSFORM_PROGRAM, 4) in LLE pgraph_handle_method()
	assert(program_load < NV2A_MAX_TRANSFORM_PROGRAM_LENGTH);
	uint32_t value = pg->program_data[program_load][slot % 4];

	return value;
}

const char *NV2AMethodToString(DWORD dwMethod)
{
	switch (dwMethod) {

#define ENUM_RANGED_ToString_N(Name, Method, Pitch, N) \
	case Name(N): return #Name "((" #N ")*" #Pitch ")";

#define ENUM_RANGED_ToString_1(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 0)

#define ENUM_RANGED_ToString_2(Name, Method, Pitch) \
	ENUM_RANGED_ToString_1(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 1)

#define ENUM_RANGED_ToString_3(Name, Method, Pitch) \
	ENUM_RANGED_ToString_2(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 2)

#define ENUM_RANGED_ToString_4(Name, Method, Pitch) \
	ENUM_RANGED_ToString_3(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 3) 

#define ENUM_RANGED_ToString_6(Name, Method, Pitch) \
	ENUM_RANGED_ToString_4(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 4) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 5)

#define ENUM_RANGED_ToString_8(Name, Method, Pitch) \
	ENUM_RANGED_ToString_6(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 6) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 7)

#define ENUM_RANGED_ToString_10(Name, Method, Pitch) \
	ENUM_RANGED_ToString_8(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 8) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 9) \

#define ENUM_RANGED_ToString_16(Name, Method, Pitch) \
	ENUM_RANGED_ToString_10(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 10) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 11) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 12) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 13) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 14) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 15)

#define ENUM_RANGED_ToString_32(Name, Method, Pitch) \
	ENUM_RANGED_ToString_16(Name, Method, Pitch) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 16) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 17) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 18) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 19) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 20) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 21) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 22) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 23) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 24) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 25) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 26) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 27) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 28) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 29) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 30) \
	ENUM_RANGED_ToString_N(Name, Method, Pitch, 31)

#define ENUM_METHOD_ToString(Name, Method) case Method: return #Name;
#define ENUM_RANGED_ToString(Name, Method, Pitch, Repeat) ENUM_RANGED_ToString_##Repeat(Name, Method, Pitch)
#define ENUM_BITFLD_Ignore(Name, Value)
#define ENUM_VALUE_Ignore(Name, Value)

	ENUM_NV2A(ENUM_METHOD_ToString, ENUM_RANGED_ToString, ENUM_BITFLD_Ignore, ENUM_VALUE_Ignore)

	default:
		return "UNLABLED";
	}
}
