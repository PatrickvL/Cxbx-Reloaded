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
// *  2020 PatrickvL
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::VTXSH // TODO : Introduce generic HLSL logging

#include <d3dcompiler.h>
#include "Shader.h"
#include "EmbeddedShaders.h"
#include "core\kernel\support\Emu.h" // EmuLog

// ============================================================================
// Precompiled CSO loading - uses bytecode embedded at compile time
// ============================================================================
bool LoadPrecompiledCSO(const char* csoName, ID3DBlob** ppBlob)
{
	*ppBlob = nullptr;

	const void* data = nullptr;
	size_t size = 0;
	if (!GetEmbeddedShaderData(csoName, &data, &size)) {
		EmuLog(LOG_LEVEL::WARNING, "LoadPrecompiledCSO: no embedded shader '%s'", csoName);
		return false;
	}

	HRESULT hr = D3DCreateBlob(size, ppBlob);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "LoadPrecompiledCSO: D3DCreateBlob failed (size=%zu)", size);
		return false;
	}

	memcpy((*ppBlob)->GetBufferPointer(), data, size);
	EmuLog(LOG_LEVEL::INFO, "LoadPrecompiledCSO: loaded embedded '%s' (%zu bytes)", csoName, size);
	return true;
}