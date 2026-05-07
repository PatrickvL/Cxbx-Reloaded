#pragma once
// ---------------------------------------------------------------------------
// EmbeddedShaders.h — Precompiled HLSL shaders linked into the executable.
//
// Generated .h headers (from gen_cso_header.cmake) are included from the
// build directory via target_include_directories.  Each provides:
//   static const unsigned char <name>[];
//   static const size_t <name>_size;
// ---------------------------------------------------------------------------
#include <cstddef>
#include <cstring>

// Generated headers — one per precompiled shader (build/generated/shaders/)
#include "cso_CxbxRCInterpreterPS.h"
#include "cso_CxbxVSInterpreterVS.h"
#include "cso_CxbxFixedFunctionVS.h"
#include "cso_CxbxBlitVS.h"
#include "cso_CxbxBlitPS.h"
#include "cso_CxbxPointSpriteGS.h"
#include "cso_CxbxThickLineGS.h"
#include "cso_CxbxUnswizzleCS.h"
#include "cso_CxbxUnswizzleBGRA_CS.h"
#include "cso_CxbxIndexConvertCS.h"
#include "cso_CxbxPaletteExpandCS.h"
#include "cso_CxbxFormatConvertCS.h"
#include "cso_CxbxVertexConvertCS.h"

struct EmbeddedShaderEntry {
    const char*          name;
    const unsigned char* data;
    size_t               size;
};

inline bool GetEmbeddedShaderData(const char* name, const void** ppData, size_t* pSize)
{
    static const EmbeddedShaderEntry kShaders[] = {
        { "CxbxRCInterpreterPS",  cso_CxbxRCInterpreterPS,  cso_CxbxRCInterpreterPS_size },
        { "CxbxVSInterpreterVS",  cso_CxbxVSInterpreterVS,  cso_CxbxVSInterpreterVS_size },
        { "CxbxFixedFunctionVS",  cso_CxbxFixedFunctionVS,  cso_CxbxFixedFunctionVS_size },
        { "CxbxBlitVS",           cso_CxbxBlitVS,           cso_CxbxBlitVS_size },
        { "CxbxBlitPS",           cso_CxbxBlitPS,           cso_CxbxBlitPS_size },
        { "CxbxPointSpriteGS",   cso_CxbxPointSpriteGS,   cso_CxbxPointSpriteGS_size },
        { "CxbxThickLineGS",     cso_CxbxThickLineGS,     cso_CxbxThickLineGS_size },
        { "CxbxUnswizzleCS",     cso_CxbxUnswizzleCS,     cso_CxbxUnswizzleCS_size },
        { "CxbxUnswizzleBGRA_CS", cso_CxbxUnswizzleBGRA_CS, cso_CxbxUnswizzleBGRA_CS_size },
        { "CxbxIndexConvertCS",  cso_CxbxIndexConvertCS,  cso_CxbxIndexConvertCS_size },
        { "CxbxPaletteExpandCS", cso_CxbxPaletteExpandCS, cso_CxbxPaletteExpandCS_size },
        { "CxbxFormatConvertCS", cso_CxbxFormatConvertCS, cso_CxbxFormatConvertCS_size },
        { "CxbxVertexConvertCS", cso_CxbxVertexConvertCS, cso_CxbxVertexConvertCS_size },
    };

    for (const auto& s : kShaders) {
        if (strcmp(s.name, name) == 0) {
            *ppData = s.data;
            *pSize = s.size;
            return true;
        }
    }
    *ppData = nullptr;
    *pSize = 0;
    return false;
}
