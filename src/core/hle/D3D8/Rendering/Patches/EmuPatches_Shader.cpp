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

// D3DDevice_LoadVertexShader — disabled (trampoline-only after CxbxImpl removal).
// Patch disabled in Patches.cpp — Xbox code runs unpatched.

// D3DDevice_SelectVertexShader, D3DDevice_SelectVertexShader_0__LTCG_eax1_ebx2,
// D3DDevice_SelectVertexShader_4__LTCG_eax1 — disabled.
// Xbox native SelectVertexShader pushes NV097_SET_TRANSFORM_PROGRAM_START.
// The host reads the program start register from PGRAPH state.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
// Test-cases: Star Wars - Battlefront (LTCG_eax1_ebx2), Aggressive Inline (LTCG_eax1)

// D3DDevice_SetShaderConstantMode — disabled.
// g_Xbox_VertexShaderConstantMode has no render-thread readers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetVertexShaderConstant variants (8 patches) — disabled.
// Xbox native SetVertexShaderConstant generates NV097_SET_TRANSFORM_CONSTANT
// push buffer commands that flow through PFIFO → PGRAPH. No HLE interception needed.
// Patches disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetTexture_4__LTCG_eax2, D3DDevice_SetTexture_4__LTCG_eax1 — disabled.
// Xbox native SetTexture pushes NV097_SET_TEXTURE_OFFSET to the push buffer.
// Host texture lookup uses PGRAPH TEXOFFSET registers set by the push buffer.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
// Test-cases: NASCAR Heat 2002 (LTCG_eax2), Metal Wolf Chaos (LTCG_eax1)

// D3DDevice_SetPixelShader, D3DDevice_SetPixelShader_0__LTCG_eax1 — disabled.
// These only called CxbxImpl_SetPixelShader (after the trampoline) to write
// g_pXbox_PixelShader, which was only read by the COMBINECTL==0 HLE bridge
// fallback, now removed. Bodies moved to Direct3D9.cpp.unused-patches.

// D3DDevice_DrawVertices_4__LTCG_ecx2_eax3, D3DDevice_DrawVertices_8__LTCG_eax3 — disabled.
// LTCG variants of DrawVertices; patches disabled in Patches.cpp.

// D3DDevice_DeleteVertexShader, D3DDevice_DeleteVertexShader_0__LTCG_eax1 — disabled.
// Host shader cache uses PGRAPH program store; no per-shader cleanup needed.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetShaderConstantMode — disabled.
// g_Xbox_VertexShaderConstantMode has no render-thread readers.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetVertexShader — disabled.
// Getter reads g_Xbox_VertexShader_Handle; Xbox native reads from device struct.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_GetVertexShaderConstant — disabled.
// Getter reads HLE VS constant shadow; Xbox native reads from device constant table.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.

// D3DDevice_SetVertexShaderInput — disabled.
// Xbox native programs NV2A vertex attribute array registers via push buffer.
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
// Test-cases: PushBuffer XDK sample, Halo 2, Kung Fu Chaos, NBA LIVE 2005,
// Prince of Persia WW, Spyro A Hero's Tail

// D3DDevice_RunVertexStateShader — disabled.
// Xbox native uses NV097_SET_TRANSFORM_DATA + NV097_LAUNCH_TRANSFORM_PROGRAM,
// now handled by PGRAPH. Patch disabled in Patches.cpp.
// Implementation moved to Direct3D9.cpp.unused-patches.

// ******************************************************************
// D3DDevice_SetDepthClipPlanes — disabled (all cases are TODO stubs).
// Patch disabled in Patches.cpp — let Xbox code run unpatched.
