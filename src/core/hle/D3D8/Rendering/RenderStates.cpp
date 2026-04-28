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
// *  (c) 2019 Luke Usher
// *
// *  All rights reserved
// *
// ******************************************************************
#define LOG_PREFIX CXBXR_MODULE::D3DST


#include "RenderStates.h"
#include "Logging.h"
#include "core/hle/D3D8/XbConvert.h"
#include "core/hle/D3D8/Rendering/RenderGlobals.h"

bool XboxRenderStateConverter::Init()
{
    // Get render state
    if (g_SymbolAddresses.find("D3D_g_RenderState") != g_SymbolAddresses.end()) {
        D3D__RenderState = (uint32_t*)g_SymbolAddresses["D3D_g_RenderState"];
    } else {
        EmuLog(LOG_LEVEL::ERROR2, "D3D_g_RenderState was not found!");
        return false;
    }

    // Build a mapping of Cxbx Render State indexes to indexes within the current XDK
    BuildRenderStateMappingTable();

    return true;
}

bool IsRenderStateAvailableInCurrentXboxD3D8Lib(RenderStateInfo& aRenderStateInfo)
{
	bool bIsRenderStateAvailable = (aRenderStateInfo.V <= g_LibVersion_D3D8);
	if (aRenderStateInfo.R > 0) { // Applies to xbox::X_D3DRS_MULTISAMPLETYPE
		// Note : X_D3DRS_MULTISAMPLETYPE seems the only render state that got
		// removed (from 4039 onwards), so we check that limitation here as well
		bIsRenderStateAvailable &= (g_LibVersion_D3D8 < aRenderStateInfo.R);
	}
	return bIsRenderStateAvailable;
}

void XboxRenderStateConverter::BuildRenderStateMappingTable()
{
    EmuLog(LOG_LEVEL::INFO, "Building Cxbx to XDK Render State Mapping Table");

    XboxRenderStateOffsets.fill(-1);

    int XboxIndex = 0;
    for (unsigned int RenderState = xbox::X_D3DRS_FIRST; RenderState <= xbox::X_D3DRS_LAST; RenderState++) {
        auto RenderStateInfo = GetDxbxRenderStateInfo(RenderState);
        if (IsRenderStateAvailableInCurrentXboxD3D8Lib(RenderStateInfo)) {
            XboxRenderStateOffsets[RenderState] = XboxIndex;
            EmuLog(LOG_LEVEL::INFO, "%s = %d", RenderStateInfo.S, XboxIndex);
            XboxIndex++;
            continue;
        }

        EmuLog(LOG_LEVEL::INFO, "%s Not Present", RenderStateInfo.S);
    }
}

bool XboxRenderStateConverter::XboxRenderStateExists(uint32_t State)
{
    if (XboxRenderStateOffsets[State] >= 0) {
        return true;
    }

    return false;
}

void XboxRenderStateConverter::SetXboxRenderState(uint32_t State, uint32_t Value)
{
    if (!XboxRenderStateExists(State)) {
        EmuLog(LOG_LEVEL::WARNING, "Attempt to write a Renderstate (%s) that does not exist in the current D3D8 XDK Version (%d)", GetDxbxRenderStateInfo(State).S, g_LibVersion_D3D8);
        return;
    }

    D3D__RenderState[XboxRenderStateOffsets[State]] = Value;
}

uint32_t XboxRenderStateConverter::GetXboxRenderState(uint32_t State)
{
    if (!XboxRenderStateExists(State)) {
        EmuLog(LOG_LEVEL::WARNING, "Attempt to read a Renderstate (%s) that does not exist in the current D3D8 XDK Version (%d)", GetDxbxRenderStateInfo(State).S, g_LibVersion_D3D8);
        return 0;
    }

    return D3D__RenderState[XboxRenderStateOffsets[State]];
}

float XboxRenderStateConverter::GetXboxRenderStateAsFloat(uint32_t State)
{
    if (!XboxRenderStateExists(State)) {
        EmuLog(LOG_LEVEL::WARNING, "Attempt to read a Renderstate (%s) that does not exist in the current D3D8 XDK Version (%d)", GetDxbxRenderStateInfo(State).S, g_LibVersion_D3D8);
        return 0;
    }

   	float f; std::memcpy(&f, &D3D__RenderState[XboxRenderStateOffsets[State]], sizeof(f)); return f;
}

void XboxRenderStateConverter::SetWireFrameMode(int wireframe)
{
    WireFrameMode = wireframe;
}
