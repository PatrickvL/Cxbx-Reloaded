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

#include <unordered_map>

#include "Cxbx.h"

#include "core\hle\D3D8\XbVertexShader.h" // for CxbxVertexShaderInfo

typedef struct _CxbxDrawContext
{
    IN     XTL::X_D3DPRIMITIVETYPE    XboxPrimitiveType;
    IN     DWORD                 dwVertexCount;
    IN     DWORD                 dwStartVertex; // Only D3DDevice_DrawVertices sets this (potentially higher than default 0)
	IN	   PWORD				 pXboxIndexData; // Set by D3DDevice_DrawIndexedVertices, D3DDevice_DrawIndexedVerticesUP and HLE_draw_inline_elements
	IN	   DWORD				 dwBaseVertexIndex; // Set to g_Xbox_BaseVertexIndex in D3DDevice_DrawIndexedVertices
	IN	   INDEX16				 LowIndex, HighIndex; // Set when pXboxIndexData is set
	IN	   size_t				 VerticesInBuffer; // Set by CxbxVertexBufferConverter::Apply
    // Data if Draw...UP call
    IN PVOID                     pXboxVertexStreamZeroData; // Set by D3DDevice_DrawVerticesUP, D3DDevice_DrawIndexedVerticesUP, HLE_draw_inline_array and EmuFlushIVB
    IN UINT                      uiXboxVertexStreamZeroStride;
	// Values to be used on host
	OUT PVOID                    pHostVertexStreamZeroData;
	OUT UINT                     uiHostVertexStreamZeroStride;
    OUT DWORD                    dwHostPrimitiveCount; // Set by CxbxVertexBufferConverter::Apply
}
CxbxDrawContext;

class CxbxPatchedStream
{
public:
    CxbxPatchedStream();
    ~CxbxPatchedStream();
    void Activate(CxbxDrawContext *pDrawContext, UINT uiStream) const;
    bool                    m_bIsCached = false; // When true, indicates this object's values have been written, false means it's not cached
    uint64_t                m_uiCachedVertexStreamInformationHash = 0;

	PVOID                   m_pCachedXboxVertexData = xbnullptr;
    UINT                    m_uiCachedXboxVertexDataSize = 0;
    uint64_t                m_uiCachedXboxVertexDataHash = 0;
    UINT                    m_uiCachedXboxVertexStride = 0;

    UINT                    m_uiCachedHostVertexStride = 0;
    bool                    m_bCachedUseHostVertexStreamZero = false;
    void                   *m_pCachedHostVertexStreamZeroData = nullptr;
    bool                    m_bCachedHostVertexStreamZeroDataIsAllocated = false;
    IDirect3DVertexBuffer  *m_pCachedHostVertexBuffer = nullptr;
};

class CxbxVertexBufferConverter
{
    public:
        CxbxVertexBufferConverter();
        void Apply(CxbxDrawContext *pPatchDesc);
        void PrintStats();
    private:
        // Stack tracking
        ULONG m_TotalCacheHits = 0;
        ULONG m_TotalCacheMisses = 0;

        UINT m_MaxCacheSize = 2000;                                        // Maximum number of entries in the cache
        UINT m_CacheElasticity = 200;                                      // Cache is allowed to grow this much more than maximum before being purged to maximum
        std::unordered_map<uint64_t, std::list<CxbxPatchedStream>::iterator> m_PatchedStreams;  // Stores references to patched streams for fast lookup
        std::list<CxbxPatchedStream> m_PatchedStreamUsageList;             // Linked list of vertex streams, least recently used is last in the list
        CxbxPatchedStream& GetPatchedStream(uint64_t);                     // Fetches (or inserts) a patched stream associated with the given key

        CxbxVertexShaderInfo *m_pVertexShaderInfo;

		bool m_bVshHandleIsFVF;
		UINT m_uiTextureCoordinatesByteOffsetInVertex;
		struct { int NrTexCoords; bool bTexIsLinear; int Width; int Height; int Depth; } m_FVFInfo_ActivePixelContainer[XTL::X_D3DTS_STAGECOUNT];
		bool m_bMustReset_FVF_XYZRHW;

		void PrepareStreamConversion();

		CxbxPatchedStream m_TempPatchedStream;
		CxbxPatchedStream &ConvertStream(CxbxDrawContext *pPatchDesc, UINT uiStream); // Prepares vertex data for the given stream number
};

// inline vertex buffer emulation
extern XTL::X_D3DPRIMITIVETYPE      g_InlineVertexBuffer_PrimitiveType;
extern DWORD                   g_InlineVertexBuffer_FVF;

extern struct _D3DIVB
{
    D3DXVECTOR3 Position;     // X_D3DVSDE_POSITION (*) > X_D3DFVF_XYZ or X_D3DFVF_XYZRHW
    FLOAT       Rhw;          // X_D3DVSDE_VERTEX (*)   > X_D3DFVF_XYZ or X_D3DFVF_XYZRHW
	FLOAT		Blend[4];	  // X_D3DVSDE_BLENDWEIGHT  > X_D3DFVF_XYZB1 (and 3 more up to X_D3DFVF_XYZB4)
    D3DXVECTOR3 Normal;       // X_D3DVSDE_NORMAL       > X_D3DFVF_NORMAL
	D3DCOLOR    Diffuse;      // X_D3DVSDE_DIFFUSE      > X_D3DFVF_DIFFUSE
	D3DCOLOR    Specular;     // X_D3DVSDE_SPECULAR     > X_D3DFVF_SPECULAR
	FLOAT       Fog;          // X_D3DVSDE_FOG          > X_D3DFVF_FOG unavailable; TODO : How to handle?
	D3DCOLOR    BackDiffuse;  // X_D3DVSDE_BACKDIFFUSE  > X_D3DFVF_BACKDIFFUSE unavailable; TODO : How to handle?
	D3DCOLOR    BackSpecular; // X_D3DVSDE_BACKSPECULAR > X_D3DFVF_BACKSPECULAR unavailable; TODO : How to handle?
    D3DXVECTOR4 TexCoord[4];  // X_D3DVSDE_TEXCOORD0    > X_D3DFVF_TEX1 (and 4 more up to X_D3DFVF_TEX4)

	// (*) X_D3DVSDE_POSITION and X_D3DVSDE_VERTEX both set Position, but Rhw seems optional,
	// hence, selection for X_D3DFVF_XYZ or X_D3DFVF_XYZRHW is rather fuzzy. We DO know that once
	// X_D3DFVF_NORMAL is given, X_D3DFVF_XYZRHW is forbidden (see D3DDevice_SetVertexData4f)
}
*g_InlineVertexBuffer_Table;

extern UINT g_InlineVertexBuffer_TableLength;
extern UINT g_InlineVertexBuffer_TableOffset;

extern VOID EmuFlushIVB();

extern VOID EmuUpdateActiveTexture();

extern DWORD g_dwPrimPerFrame;
 
#endif
