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
#include "../Backend\Backend_D3D11.h"

// D3DDevice_Begin, D3DDevice_SetVertexData2f, D3DDevice_SetVertexData2s,
// D3DDevice_SetVertexData4f_16__LTCG_edi1, D3DDevice_SetVertexData4f,
// D3DDevice_SetVertexData4ub, D3DDevice_SetVertexData4s,
// D3DDevice_SetVertexDataColor, D3DDevice_End — disabled.
// Xbox native Begin/End/SetVertexData pushes NV2A methods (NV097_SET_BEGIN_END,
// NV097_ARRAY_ELEMENT16 etc.) through the push buffer → PFIFO → PGRAPH.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_BeginPushBuffer, D3DDevice_BeginPushBuffer_0__LTCG_edi1,
// D3DDevice_EndPushBuffer — disabled.
// These only called the trampoline and toggled g_bRecordingPushBuffer
// which had zero readers. Xbox code runs unpatched.
// Patch disabled in Patches.cpp.

// D3DDevice_RunPushBuffer, D3DDevice_RunPushBuffer_4__LTCG_eax2 — disabled.
// Native RunPushBuffer applies fixups and pushes commands through
// the real GPU FIFO.  D3D_BlockOnTime drains the FIFO when the ring fills.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_DrawVertices, DrawVerticesUP, DrawVerticesUP_12__LTCG_ebx3,
// DrawIndexedVertices, DrawIndexedVerticesUP — disabled.
// Xbox native draw calls push NV2A methods (NV097_SET_BEGIN_END, NV097_DRAW_ARRAYS,
// NV097_ARRAY_ELEMENT16, etc.) through the push buffer → PFIFO → PGRAPH.
// Draws are now initiated from NV097 writes in the LLE path.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// CDevice_SetStateVB, CDevice_SetStateVB_8, CDevice_SetStateUP,
// CDevice_SetStateUP_4, CDevice_SetStateUP_0__LTCG_esi1 — disabled.
// These were unimplemented stubs (LOG_UNIMPLEMENTED). Xbox native code pushes
// all necessary NV2A state through the push buffer before draw calls.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetStipple — disabled (unimplemented stub, LOG_IGNORED).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetSwapCallback — disabled.
// Xbox native stores callback in device struct. Swap callback not invoked without this.
// TODO: Read callback pointer from Xbox device struct in Swap path if needed.
// Patch disabled in Patches.cpp.

// D3DDevice_PrimeVertexCache — disabled (unimplemented stub, LOG_UNIMPLEMENTED).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_DrawRectPatch, D3DDevice_DrawTriPatch — disabled.
// Xbox does CPU tessellation then submits vertices via push buffer.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
