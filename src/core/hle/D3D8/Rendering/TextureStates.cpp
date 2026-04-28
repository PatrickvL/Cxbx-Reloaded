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

#include "TextureStates.h"
#include "core\kernel\init\CxbxKrnl.h"
#include "core\kernel\support\Emu.h"
#include "Logging.h"
#include "EmuShared.h"
#include "core/hle/Intercept.hpp"

typedef struct {
    const char* S;          // String representation.
    bool IsSamplerState;    // True if the state maps to a Sampler State instead of Texture Stage
    DWORD PC;               // PC Index
} TextureStateInfo;

TextureStateInfo CxbxTextureStateInfo[] = {
    { "D3DTSS_ADDRESSU",                true,   D3DSAMP_ADDRESSU },
    { "D3DTSS_ADDRESSV",                true,   D3DSAMP_ADDRESSV },
    { "D3DTSS_ADDRESSW",                true,   D3DSAMP_ADDRESSW },
    { "D3DTSS_MAGFILTER",               true,   D3DSAMP_MAGFILTER },
    { "D3DTSS_MINFILTER",               true,   D3DSAMP_MINFILTER },
    { "D3DTSS_MIPFILTER",               true,   D3DSAMP_MIPFILTER },
    { "D3DTSS_MIPMAPLODBIAS",           true,   D3DSAMP_MIPMAPLODBIAS },
    { "D3DTSS_MAXMIPLEVEL",             true,   D3DSAMP_MAXMIPLEVEL },
    { "D3DTSS_MAXANISOTROPY",           true,   D3DSAMP_MAXANISOTROPY },
    { "D3DTSS_COLORKEYOP",              false,  0 },
    { "D3DTSS_COLORSIGN",               false,  0 },
    { "D3DTSS_ALPHAKILL",               false,  0 },
    { "D3DTSS_COLOROP",                 false,  D3DTSS_COLOROP },
    { "D3DTSS_COLORARG0",               false,  D3DTSS_COLORARG0 },
    { "D3DTSS_COLORARG1",               false,  D3DTSS_COLORARG1 },
    { "D3DTSS_COLORARG2",               false,  D3DTSS_COLORARG2 },
    { "D3DTSS_ALPHAOP",                 false,  D3DTSS_ALPHAOP },
    { "D3DTSS_ALPHAARG0",               false,  D3DTSS_ALPHAARG0 },
    { "D3DTSS_ALPHAARG1",               false,  D3DTSS_ALPHAARG1 },
    { "D3DTSS_ALPHAARG2",               false,  D3DTSS_ALPHAARG2 },
    { "D3DTSS_RESULTARG",               false,  D3DTSS_RESULTARG },
    { "D3DTSS_TEXTURETRANSFORMFLAGS",   false,  D3DTSS_TEXTURETRANSFORMFLAGS },
    { "D3DTSS_BUMPENVMAT00",            false,  D3DTSS_BUMPENVMAT00 },
    { "D3DTSS_BUMPENVMAT01",            false,  D3DTSS_BUMPENVMAT01 },
    { "D3DTSS_BUMPENVMAT11",            false,  D3DTSS_BUMPENVMAT11 },
    { "D3DTSS_BUMPENVMAT10",            false,  D3DTSS_BUMPENVMAT10 },
    { "D3DTSS_BUMPENVLSCALE",           false,  D3DTSS_BUMPENVLSCALE },
    { "D3DTSS_BUMPENVLOFFSET",          false,  D3DTSS_BUMPENVLOFFSET },
    { "D3DTSS_TEXCOORDINDEX",           false,  D3DTSS_TEXCOORDINDEX },
    { "D3DTSS_BORDERCOLOR",             true,   D3DSAMP_BORDERCOLOR },
    { "D3DTSS_COLORKEYCOLOR",           false,  0 },
};

bool XboxTextureStateConverter::Init()
{
    // Deferred states start at 0, this means that D3D_g_DeferredTextureState IS D3D__TextureState
    // No further works is required to derive the offset
    if (g_SymbolAddresses.find("D3D_g_DeferredTextureState") != g_SymbolAddresses.end()) {
        D3D__TextureState = (uint32_t*)g_SymbolAddresses["D3D_g_DeferredTextureState"];
    } else {
        EmuLog(LOG_LEVEL::ERROR2, "D3D_g_DeferredTextureState was not found!");
        return false;
    }

    // Build a mapping of Cxbx Texture State indexes to indexes within the current XDK
    BuildTextureStateMappingTable();

    return true;
}

void XboxTextureStateConverter::BuildTextureStateMappingTable()
{
    EmuLog(LOG_LEVEL::INFO, "Building Cxbx to XDK Texture State Mapping Table");
    for (int State = xbox::X_D3DTSS_FIRST; State <= xbox::X_D3DTSS_LAST; State++) {
        int index = State;

        // On early XDKs, we need to shuffle the values around a little
        // TODO: Verify which XDK version this change occurred at
        // Values range 0-9 (D3DTSS_COLOROP to D3DTSS_TEXTURETRANSFORMFLAGS) become 12-21
        // Values 10-21 (D3DTSS_ADDRESSU to D3DTSS_ALPHAKILL) become 0-11
        bool bOldOrder = g_LibVersion_D3D8 <= 3948; // Verfied old order in 3944, new order in 4039

        if (bOldOrder) {
            if (State <= 9) {
                index += 12;
            } else if (State <= 21) {
                index -= 10;
            }
        }

        EmuLog(LOG_LEVEL::INFO, "%s = %d", CxbxTextureStateInfo[State].S, index);
        XboxTextureStateOffsets[index] = State;
    }
}

// Normalize values which may have different mappings per XDK version
DWORD NormalizeValue(DWORD xboxState, DWORD value) {
    if (g_LibVersion_D3D8 <= 3948) {
        // D3DTOP verified old order in 3948, new order in 4039
        switch (xboxState) {
        case xbox::X_D3DTSS_COLOROP:
        case xbox::X_D3DTSS_ALPHAOP:
            switch (value) {
            case 13:
                return xbox::X_D3DTOP_BLENDTEXTUREALPHA;
            case 14:
                return xbox::X_D3DTOP_BLENDFACTORALPHA;
            case 15:
                return xbox::X_D3DTOP_BLENDTEXTUREALPHAPM;
            case 16:
                return xbox::X_D3DTOP_BLENDCURRENTALPHA;
            }
        }
    }

    return value;
}

uint32_t XboxTextureStateConverter::Get(int textureStage, DWORD xboxState) {
    if (textureStage < 0 || textureStage > 3)
        CxbxrAbort("Requested texture stage was out of range: %d", textureStage);
    if (xboxState < xbox::X_D3DTSS_FIRST || xboxState > xbox::X_D3DTSS_LAST)
        CxbxrAbort("Requested texture state was out of range: %d", xboxState);

    // Read the value of the current stage/state from the Xbox data structure
    DWORD rawValue = D3D__TextureState[(textureStage * xbox::X_D3DTS_STAGESIZE) + XboxTextureStateOffsets[xboxState]];

    return NormalizeValue(xboxState, rawValue);
}
