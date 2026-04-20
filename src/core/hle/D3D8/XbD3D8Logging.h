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
// *  (c) 2017 Patrick van Logchem <pvanlogchem@gmail.com>
// *
// *  All rights reserved
// *
// ******************************************************************
#ifndef _EMU_D3D8_LOGGING_H
#define _EMU_D3D8_LOGGING_H

#include "Logging.h"
#include "XbD3D8Types.h"

extern const char* D3DErrorString(HRESULT hResult); // Implemented in RenderGlobals.cpp

#define DEBUG_D3DRESULT(hRet, message) \
	do { \
		LOG_CHECK_ENABLED(LOG_LEVEL::DEBUG) { \
			if (FAILED(hRet)) \
				if(g_bPrintfOn) \
					printf("%s%s : %s D3D error (0x%.08X: %s)\n", _logThreadPrefix.c_str(), _logFuncPrefix.c_str(), message, hRet, D3DErrorString(hRet)); \
		} \
	} while (0)

// Additional types, exclusively for logging (not really enums) :

namespace xbox {

enum X_D3DUSAGE : int;
enum X_D3DCOMMON_TYPE : int;
enum X_D3DRESOURCE_COMMON : int;
enum X_D3DRESOURCE_FORMAT : int;
enum X_D3DRESOURCE_SIZE : int;

} // end of namespace xbox

//
// Headers for rendering host D3D enum types :
//

// D3D11-specific logging headers
ENUM2STR_HEADER(DXGI_FORMAT)

namespace xbox {

//
// Headers for rendering Xbox D3D enum types :
//

// TODO: fill out these enumeration tables for convienance
//ENUM2STR_HEADER(X_D3DBLENDOP)
//ENUM2STR_HEADER(X_D3DBLEND)
//ENUM2STR_HEADER(X_D3DCMPFUNC)
//ENUM2STR_HEADER(X_D3DFILLMODE)
//ENUM2STR_HEADER(X_D3DSHADEMODE)
//ENUM2STR_HEADER(X_D3DSTENCILOP)
//ENUM2STR_HEADER(X_D3DTEXTURESTAGESTATETYPE)

ENUM2STR_HEADER(X_D3DCUBEMAP_FACES)
ENUM2STR_HEADER(X_D3DCULL)
//ENUM2STR_HEADER(X_D3DDEVTYPE)
ENUM2STR_HEADER(X_D3DFORMAT)
ENUM2STR_HEADER(X_D3DPRIMITIVETYPE)
LOGRENDER_HEADER(X_RECT)
ENUM2STR_HEADER(X_D3DRESOURCETYPE)
ENUM2STR_HEADER(X_D3DSET_DEPTH_CLIP_PLANES_FLAGS)
FLAGS2STR_HEADER(X_D3DUSAGE) // Not really an enum

//
// Cxbx D3D LOGRENDER(Type) implementations
//

ENUM2STR_HEADER(X_D3DCOMMON_TYPE)
LOGRENDER_HEADER(X_D3DRESOURCE_COMMON)
LOGRENDER_HEADER(X_D3DRESOURCE_FORMAT)
LOGRENDER_HEADER(X_D3DRESOURCE_SIZE)

//
// Xbox D3D LOGRENDER_HEADER(Type) declarations
//

LOGRENDER_HEADER(X_D3DVIEWPORT8)
LOGRENDER_HEADER(X_D3DDISPLAYMODE)
LOGRENDER_HEADER(X_D3DResource)
LOGRENDER_HEADER(X_D3DPixelContainer)
LOGRENDER_HEADER(X_D3DBOX)
LOGRENDER_HEADER(X_D3DLOCKED_BOX)
LOGRENDER_HEADER(X_D3DLOCKED_RECT)
LOGRENDER_HEADER(X_D3DRECT)

} // end of namespace xbox

#endif _EMU_D3D8_LOGGING_H
