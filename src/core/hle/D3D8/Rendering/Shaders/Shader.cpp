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
#include "core\kernel\support\Emu.h" // EmuLog

#include <string>

// ============================================================================
// Precompiled CSO loading â€” loads build-time compiled shaders from exe dir
// ============================================================================
bool LoadPrecompiledCSO(const char* csoName, ID3DBlob** ppBlob)
{
	*ppBlob = nullptr;

	// Build path: <exe_dir>/<csoName>.cso
	char exePath[MAX_PATH] = {};
	GetModuleFileNameA(nullptr, exePath, MAX_PATH);
	std::string path(exePath);
	auto lastSlash = path.find_last_of("\\/");
	if (lastSlash != std::string::npos)
		path = path.substr(0, lastSlash + 1);
	path += csoName;
	path += ".cso";

	FILE* fp = fopen(path.c_str(), "rb");
	if (!fp) {
		EmuLog(LOG_LEVEL::WARNING, "LoadPrecompiledCSO: file not found: %s", path.c_str());
		return false;
	}

	fseek(fp, 0, SEEK_END);
	long size = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	if (size < 8) {
		EmuLog(LOG_LEVEL::WARNING, "LoadPrecompiledCSO: file too small: %s (%ld bytes)", path.c_str(), size);
		fclose(fp);
		return false;
	}

	HRESULT hr = D3DCreateBlob(size, ppBlob);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "LoadPrecompiledCSO: D3DCreateBlob failed (size=%ld)", size);
		fclose(fp);
		return false;
	}

	size_t readBytes = fread((*ppBlob)->GetBufferPointer(), 1, size, fp);
	fclose(fp);

	if ((long)readBytes != size) {
		EmuLog(LOG_LEVEL::WARNING, "LoadPrecompiledCSO: partial read %s (%zu / %ld)", path.c_str(), readBytes, size);
		(*ppBlob)->Release();
		*ppBlob = nullptr;
		return false;
	}

	EmuLog(LOG_LEVEL::INFO, "LoadPrecompiledCSO: loaded %s (%ld bytes)", csoName, size);
	return true;
}
