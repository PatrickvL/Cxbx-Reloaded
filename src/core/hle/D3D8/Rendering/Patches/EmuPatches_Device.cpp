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

// D3DDevice_BeginPush / D3DDevice_EndPush — disabled.
// Xbox native BeginPush returns pointer into GPU ring buffer (CDevice+0x00 = pPut).
// PFIFO processes commands from ring buffer when DMA_PUT is updated by EndPush.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetIndices, D3DDevice_SetIndices_4__LTCG_ebx1 — disabled.
// Trampoline-only after DrawIndexedVertices was disabled (g_Xbox_BaseVertexIndex had no readers).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// All Direct3D_CreateDevice patches disabled — host D3D11 device created by CxbxInitHostD3DDevice(),
// native Xbox CreateDevice runs unpatched.
// See HostRender.cpp::CxbxInitHostD3DDevice() and Patches.cpp for details.

// Overload for logging
static void D3DDevice_SetIndices_4__LTCG_ebx1
(
   	xbox::X_D3DIndexBuffer   *pIndexData,
   	xbox::uint_xt             BaseVertexIndex
)
{
   	LOG_FUNC_BEGIN
   	   	LOG_FUNC_ARG(pIndexData)
   	   	LOG_FUNC_ARG(BaseVertexIndex)
   	   	LOG_FUNC_END;
}

// ******************************************************************
// * patch: D3DDevice_SetIndices_4__LTCG_ebx1
// LTCG specific D3DDevice_SetIndices function...
// D3DDevice_SetIndices_4__LTCG_ebx1 — removed (see comment above SetIndices).

// D3DDevice_Reset — disabled.
// Xbox D3DDevice_Reset only does Xbox D3D resource cleanup.
// Host rendering driven by NV2A state needs no reset — host resources
// are managed by PGRAPH RT cache and invalidated by dirty-page detection.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.


// ******************************************************************
// * patch: D3DDevice_GetDisplayFieldStatus (DISABLED)
// * Xbox native reads VBlank count from device struct and field from NV_PCRTC_RASTER bit 20.
// * PCRTC emulation now provides VERT_BLANK and FIELD bits; VBlank ISR increments the kernel counter.
// * Patch disabled in Patches.cpp — implementation moved to Direct3D9.cpp.unused-patches.
// ******************************************************************


// ******************************************************************
// * patch: D3DDevice_BeginVisibilityTest
// ******************************************************************
