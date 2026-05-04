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
#ifndef XBVERTEXSHADER_H
#define XBVERTEXSHADER_H

#include <d3dcompiler.h>

#include "core\hle\D3D8\XbD3D8Types.h" // for X_VSH_MAX_ATTRIBUTES

enum class VertexShaderMode {
	FixedFunction,
	// When titles use Xbox fixed function with pre-transformed vertices
	// it actually uses a special "passthrough" shader program
	Passthrough,
	ShaderProgram
};

extern VertexShaderMode g_Xbox_VertexShaderMode;

extern ID3DBlob* CxbxGetActiveVertexShaderBytecode();
extern ID3DBlob* CxbxGetFixedFunctionVertexShaderBytecode();

extern void CxbxrImpl_RunVertexStateShader(DWORD Address, CONST FLOAT* pData);
#endif
