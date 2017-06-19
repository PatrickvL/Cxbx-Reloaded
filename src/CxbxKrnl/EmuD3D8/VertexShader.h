// ******************************************************************
// *
// *    .,-:::::    .,::      .::::::::.    .,::      .:
// *  ,;;;'````'    `;;;,  .,;;  ;;;'';;'   `;;;,  .,;;
// *  [[[             '[[,,[['   [[[__[[\.    '[[,,[['
// *  $$$              Y$$$P     $$""""Y$$     Y$$$P
// *  `88bo,__,o,    oP"``"Yo,  _88o,,od8P   oP"``"Yo,
// *    "YUMMMMMP",m"       "Mm,""YUMMMP" ,m"       "Mm,
// *
// *   Cxbx->Win32->CxbxKrnl->EmuD3D8->VertexShader.h
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
#ifndef VERTEXSHADER_H
#define VERTEXSHADER_H

#include "Cxbx.h"

// nv2a microcode header
typedef struct
{
    uint08 Type;
    uint08 Version;
    uint08 NumInst;
    uint08 Unknown0;
}
VSH_SHADER_HEADER;

#define VSH_INSTRUCTION_SIZE       4
#define VSH_INSTRUCTION_SIZE_BYTES (VSH_INSTRUCTION_SIZE * sizeof(DWORD))

// recompile xbox vertex shader declaration
extern DWORD EmuRecompileVshDeclaration
(
    DWORD                *pDeclaration,
    DWORD               **ppRecompiledDeclaration,
    DWORD                *pDeclarationSize,
    boolean               IsFixedFunction,
    CxbxVertexShaderDynamicPatch *pVertexDynamicPatch
);

// recompile xbox vertex shader function
extern HRESULT EmuRecompileVshFunction
(
    DWORD        *pFunction,
    LPD3DXBUFFER *ppRecompiled,
    DWORD        *pOriginalSize,
    boolean      bNoReservedConstants,
	boolean		 *pbUseDeclarationOnly
);

extern void FreeVertexDynamicPatch(CxbxVertexShader *pVertexShader);

// Checks for failed vertex shaders, and shaders that would need patching
extern boolean IsValidCurrentShader(void);
extern boolean VshHandleIsValidShader(DWORD Handle);

// TODO : FVF actually uses bit 16 up to 23 for texture sizes - instead, bit 1 indicates non-FVF shader handles!
// Dxbx note : On Xbox, a FVF is recognizable when the handle <= 0x0000FFFF
// (as all values above are allocated VertexShader addresses).
inline boolean VshHandleIsFVF(DWORD Handle) { return (Handle > NULL) && (Handle <= 0x0000FFFF); }
inline boolean VshHandleIsVertexShader(DWORD Handle) { return (Handle > 0x0000FFFF) ? TRUE : FALSE; }
inline X_D3DVertexShader *VshHandleGetVertexShader(DWORD Handle) { return VshHandleIsVertexShader(Handle) ? (X_D3DVertexShader *)Handle : nullptr; }

#ifdef _DEBUG_TRACK_VS
#define DbgVshPrintf if(g_bPrintfOn) printf
#else
inline void null_func_vsh(...) { }
#define DbgVshPrintf XTL::null_func_vsh
#endif

#endif
