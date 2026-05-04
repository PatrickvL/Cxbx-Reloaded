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
#ifndef XBVERTEXBUFFER_H
#define XBVERTEXBUFFER_H

#include "Cxbx.h"

#include "core\hle\D3D8\XbVertexShader.h"

typedef struct _CxbxDrawContext
{
    IN     xbox::X_D3DPRIMITIVETYPE    XboxPrimitiveType;
    IN     DWORD                 dwVertexCount;
    IN     DWORD                 dwStartVertex; // Only D3DDevice_DrawVertices sets this (potentially higher than default 0)
	IN	   PWORD				 pXboxIndexData; // Set by D3DDevice_DrawIndexedVertices, D3DDevice_DrawIndexedVerticesUP and D3D11_draw_inline_elements
	IN	   DWORD				 dwBaseVertexIndex; // Always 0 now (D3DDevice_SetIndices/DrawIndexedVertices patches disabled)
    // Data if Draw...UP call
    IN PVOID                     pXboxVertexStreamZeroData;
    IN UINT                      uiXboxVertexStreamZeroStride;
	IN bool                      bNV2AInlineData; // True when vertex data comes from NV2A inline_array (use PGRAPH attributes for layout)
}
CxbxDrawContext;

#endif
