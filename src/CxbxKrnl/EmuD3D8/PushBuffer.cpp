// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
// ******************************************************************
// *
// *    .,-:::::    .,::      .::::::::.    .,::      .:
// *  ,;;;'````'    `;;;,  .,;;  ;;;'';;'   `;;;,  .,;;
// *  [[[             '[[,,[['   [[[__[[\.    '[[,,[['
// *  $$$              Y$$$P     $$""""Y$$     Y$$$P
// *  `88bo,__,o,    oP"``"Yo,  _88o,,od8P   oP"``"Yo,
// *    "YUMMMMMP",m"       "Mm,""YUMMMP" ,m"       "Mm,
// *
// *   Cxbx->Win32->CxbxKrnl->EmuD3D->PushBuffer.cpp
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
#define _CXBXKRNL_INTERNAL
#define _XBOXKRNL_DEFEXTRN_

#include "CxbxKrnl/Emu.h"
#include "CxbxKrnl/EmuXTL.h"
#include "CxbxKrnl/EmuD3D8Types.h" // For X_D3DFORMAT
#include "CxbxKrnl/ResourceTracker.h"
#include "CxbxKrnl/MemoryManager.h"
#include "State.h"

uint32  XTL::g_dwPrimaryPBCount = 0;
uint32 *XTL::g_pPrimaryPB = 0;

bool XTL::g_bStepPush = false;
bool XTL::g_bSkipPush = false;
bool XTL::g_bBrkPush  = false;

bool g_bPBSkipPusher = false;

static void DbgDumpMesh(XTL::INDEX16 *pIndexData, DWORD dwCount);

int XTL::DxbxFVF_GetTextureSize(DWORD dwFVF, int aTextureIndex)
// Determine the size (in bytes) of the texture format (indexed 0 .. 3).
// This is the reverse of the D3DFVF_TEXCOORDSIZE[0..3] macros.
{
	switch ((dwFVF >> ((aTextureIndex * 2) + 16)) & 3) {
	case D3DFVF_TEXTUREFORMAT1: return 1 * sizeof(FLOAT);
	case D3DFVF_TEXTUREFORMAT2: return 2 * sizeof(FLOAT);
	case D3DFVF_TEXTUREFORMAT3: return 3 * sizeof(FLOAT);
	case D3DFVF_TEXTUREFORMAT4: return 4 * sizeof(FLOAT);
	default:
		//assert(false || "DxbxFVF_GetTextureSize : Unhandled case");
		return 0;
	}
}

// Dxbx Note: This code is taken from EmuExecutePushBufferRaw and occured
// in EmuFlushIVB too, so it's generalize in this single implementation.
UINT XTL::DxbxFVFToVertexSizeInBytes(DWORD dwFVF, BOOL bIncludeTextures)
{
/*
	X_D3DFVF_POSITION_MASK    = $00E; // Dec  /2  #fl

	X_D3DFVF_XYZ              = $002; //  2 > 1 > 3
	X_D3DFVF_XYZRHW           = $004; //  4 > 2 > 4
	X_D3DFVF_XYZB1            = $006; //  6 > 3 > 4
	X_D3DFVF_XYZB2            = $008; //  8 > 4 > 5
	X_D3DFVF_XYZB3            = $00a; // 10 > 5 > 6
	X_D3DFVF_XYZB4            = $00c; // 12 > 6 > 7
*/
	// Divide the D3DFVF by two, this gives almost the number of floats needed for the format :
	UINT Result = (dwFVF & D3DFVF_POSITION_MASK) >> 1;
	if (Result >= (D3DFVF_XYZB1 >> 1)) {
		// Any format from D3DFVF_XYZB1 and above need 1 extra float :
		Result++;
	}
	else {
		// The other formats (XYZ and XYZRHW) need 2 extra floats :
		Result += 2;
	}

	// Express the size in bytes, instead of floats :
	Result *= sizeof(FLOAT);
	// D3DFVF_NORMAL cannot be combined with D3DFVF_XYZRHW :
	if ((dwFVF & D3DFVF_POSITION_MASK) != D3DFVF_XYZRHW) {
		if (dwFVF & D3DFVF_NORMAL) {
			Result += sizeof(FLOAT) * 3;
		}
	}

	if (dwFVF & D3DFVF_DIFFUSE) {
		Result += sizeof(XTL::D3DCOLOR);
	}

	if (dwFVF & D3DFVF_SPECULAR) {
		Result += sizeof(XTL::D3DCOLOR);
	}

	if (bIncludeTextures) {
		int NrTextures = ((dwFVF & D3DFVF_TEXCOUNT_MASK) >> D3DFVF_TEXCOUNT_SHIFT);
		while (NrTextures > 0) {
			NrTextures--;
			Result += DxbxFVF_GetTextureSize(dwFVF, NrTextures);
		}
	}

	return Result;
}

void XTL::EmuExecutePushBuffer
(
    X_D3DPushBuffer       *pPushBuffer,
    X_D3DFixup            *pFixup
)
{
    if (pFixup != NULL)
        CxbxKrnlCleanup("PushBuffer has fixups\n");

#ifdef _DEBUG_TRACK_PB
	DbgDumpPushBuffer((DWORD*)pPushBuffer->Data, pPushBuffer->Size);
#endif

    EmuExecutePushBufferRaw((DWORD*)pPushBuffer->Data);

    return;
}

#define NV2A_JMP_FLAG          0x00000001
#define NV2A_CALL_FLAG         0x00000002 // TODO : Should JMP & CALL be switched?
#define NV2A_ADDR_MASK         0xFFFFFFFC
#define NV2A_METHOD_MASK       0x00001FFC
#define NV2A_SUBCH_MASK        0x0000E000
#define NV2A_COUNT_MASK        0x0FFF0000 // 12 bits
#define NV2A_NOINCREMENT_FLAG  0x40000000
// Dxbx note : What do the other bits mean (mask $B0000000) ?

#define NV2A_METHOD_SHIFT 0 // Dxbx note : Not 2, because methods are actually DWORD offsets (and thus defined with increments of 4)
#define NV2A_SUBCH_SHIFT 12
#define NV2A_COUNT_SHIFT 18

#define NV2A_METHOD_MAX ((NV2A_METHOD_MASK | 3) >> NV2A_METHOD_SHIFT) // = 8191
#define NV2A_COUNT_MAX ((NV2A_COUNT_MASK >> NV2A_COUNT_SHIFT) - 1) // = 2047

void D3DPUSH_DECODE(const DWORD dwPushCommand, DWORD &dwMethod, DWORD &dwSubCh, DWORD &dwCount, BOOL &bNoInc)
{
	dwMethod = (dwPushCommand & NV2A_METHOD_MASK); // >> NV2A_METHOD_SHIFT;
	dwSubCh = (dwPushCommand & NV2A_SUBCH_MASK) >> NV2A_SUBCH_SHIFT;
	dwCount = (dwPushCommand & NV2A_COUNT_MASK) >> NV2A_COUNT_SHIFT;
	bNoInc = (dwPushCommand & NV2A_NOINCREMENT_FLAG) > 0;
}

// Globals and controller :
DWORD *pdwOrigPushData;
bool bShowPB = false;
DWORD* pdwPushArguments;

DWORD PrevMethod[2] = {};
XTL::INDEX16 *pIndexData = NULL;
PVOID pVertexData = NULL;

DWORD dwVertexShader = -1;
DWORD dwVertexStride = -1;

// cache of last 4 indices
XTL::INDEX16 pIBMem[4] = { 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF };

XTL::X_D3DPRIMITIVETYPE  XboxPrimitiveType = XTL::X_D3DPT_INVALID;

DWORD dwMethod;
DWORD dwSubCh;
DWORD dwCount;
BOOL bNoInc;

int HandledCount;
char *HandledBy;

DWORD XTL::NV2AInstance_Registers[8192] = {};

typedef void (*NV2ACallback_t)();

NV2ACallback_t NV2ACallbacks[8192] = {};

using namespace XTL;

void EmuNV2A_NOP() // 0x0100
{
	// HandledBy = "nop";
	HandledCount = dwCount;

	// TODO : NOP actually seems to do something on NV2A
}

void NVPB_SetBeginEnd() // 0x000017FC
{
	if (*pdwPushArguments == 0) {
#ifdef _DEBUG_TRACK_PB
		if (bShowPB) {
			printf("DONE)\n");
		}
#endif
		HandledBy = "DrawEnd()";
	}
	else {
#ifdef _DEBUG_TRACK_PB
		if (bShowPB) {
			printf("PrimitiveType := %d)\n", *pdwPushArguments);
		}
#endif

		XboxPrimitiveType = (X_D3DPRIMITIVETYPE)*pdwPushArguments;
		HandledBy = "DrawBegin()";
	}
}

void NVPB_InlineVertexArray() // 0x1818
{
	HandledBy = "DrawVertices()";
	HandledCount = dwCount;

	pVertexData = pdwPushArguments;

	// retrieve vertex shader
	g_pD3DDevice8->GetVertexShader(&dwVertexShader);
	if (dwVertexShader > 0xFFFF) {
		CxbxKrnlCleanup("Non-FVF Vertex Shaders not yet supported for PushBuffer emulation!");
		dwVertexShader = 0;
	}
	else if (dwVertexShader == NULL) {
		EmuWarning("FVF Vertex Shader is null");
		dwVertexShader = -1;
	}
	/*else if (dwVertexShader == 0x6) {
	dwVertexShader = (D3DFVF_XYZ|D3DFVF_DIFFUSE|D3DFVF_TEX1);
	}*/

	//	printf( "EmuExecutePushBufferRaw: FVF = 0x%.08X\n" );

	//
	// calculate stride
	//

	dwVertexStride = 0;
	if (VshHandleIsFVF(dwVertexShader)) {
		dwVertexStride = DxbxFVFToVertexSizeInBytes(dwVertexShader, /*bIncludeTextures=*/true);
	}

	/*
	// create cached vertex buffer only once, with maxed out size
	if (pVertexBuffer == 0) {
	HRESULT hRet = g_pD3DDevice8->CreateVertexBuffer(2047*sizeof(DWORD), D3DUSAGE_WRITEONLY, dwVertexShader, D3DPOOL_MANAGED, &pVertexBuffer);
	if (FAILED(hRet))
	CxbxKrnlCleanup("Unable to create vertex buffer cache for PushBuffer emulation (0x1818, dwCount : %d)", dwCount);
	}

	// copy vertex data
	{
	uint08 *pData = nullptr;
	HRESULT hRet = pVertexBuffer->Lock(0, dwCount * sizeof(DWORD), &pData, D3DLOCK_DISCARD);
	if (FAILED(hRet))
	CxbxKrnlCleanup("Unable to lock vertex buffer cache for PushBuffer emulation (0x1818, dwCount : %d)", dwCount);

	memcpy(pData, pVertexData, dwCount * sizeof(DWORD));
	pVertexBuffer->Unlock();
	}
	*/

#ifdef _DEBUG_TRACK_PB
	if (bShowPB) {
		printf("NVPB_InlineVertexArray(...)\n");
		printf("  dwCount : %d\n", dwCount);
		printf("  dwVertexShader : 0x%08X\n", dwVertexShader);
	}
#endif

	// render vertices
	if (dwVertexShader != -1) {
		UINT VertexCount = (dwCount * sizeof(DWORD)) / dwVertexStride;
		CxbxDrawContext DrawContext = {};
		DrawContext.XboxPrimitiveType = XboxPrimitiveType;
		DrawContext.dwVertexCount = VertexCount;
		DrawContext.pXboxVertexStreamZeroData = pVertexData;
		DrawContext.uiXboxVertexStreamZeroStride = dwVertexStride;
		DrawContext.hVertexShader = dwVertexShader;

		CxbxDrawPrimitiveUP(DrawContext);
	}
}

void NVPB_FixLoop() // 0x1808
{
	HandledBy = "DrawIndexedVertices()";
	HandledCount = dwCount;

	// Test case : Turok menu's
#ifdef _DEBUG_TRACK_PB
	if (bShowPB) {
		printf("  NVPB_FixLoop(%d)\n", dwCount);
		printf("\n");
		printf("  Index Array Data...\n");
		INDEX16 *pIndices = (INDEX16*)pdwPushArguments;
		for (uint s = 0; s < dwCount; s++) {
			if (s % 8 == 0)
				printf("\n  ");

			printf("  %.04X", *pIndices++);
		}

		printf("\n");
		printf("\n");
	}
#endif

	INDEX16 *pIndices = (INDEX16*)pdwPushArguments;
	for (uint mi = 0; mi < dwCount; mi++) {
		pIBMem[mi + 2] = pIndices[mi];
	}

	// perform rendering
	if (pIBMem[0] != 0xFFFF) {
		UINT uiIndexCount = dwCount + 2;
#ifdef _DEBUG_TRACK_PB
		if (!g_PBTrackDisable.exists(pdwOrigPushData))
#endif
			// render indexed vertices
		{
			if (!g_bPBSkipPusher) {
				if (IsValidCurrentShader()) {
					// TODO: This technically should be enabled
					XTL::DxbxUpdateDeferredStates(); // CxbxUpdateNativeD3DResources

					CxbxDrawContext DrawContext = {};
					DrawContext.XboxPrimitiveType = XboxPrimitiveType;
					DrawContext.dwVertexCount = EmuD3DIndexCountToVertexCount(XboxPrimitiveType, uiIndexCount);
					DrawContext.hVertexShader = g_CurrentVertexShader;

					CxbxDrawIndexed(DrawContext, pIBMem);
				}
			}
		}
	}
}

void NVPB_InlineIndexArray() // 0x1800
{
	HandledBy = "DrawIndices()";
	HandledCount = dwCount;

	// Test case : Turok menu's

	pIndexData = (INDEX16*)pdwPushArguments;
#ifdef _DEBUG_TRACK_PB
	if (bShowPB) {
		printf("  NVPB_InlineIndexArray(0x%.08X, %d)...\n", pIndexData, dwCount);
		printf("\n");
		printf("  Index Array Data...\n");
		INDEX16 *pIndices = pIndexData;
		for (uint s = 0; s < dwCount; s++) {
			if (s % 8 == 0)
				printf("\n  ");

			printf("  %.04X", *pIndices++);
		}

		printf("\n");

#if 0
		// retrieve stream data
		XTL::IDirect3DVertexBuffer8 *pActiveVB = nullptr;
		UINT  uiStride;

		// pActiveVB = CxbxUpdateVertexBuffer(Xbox_g_Stream[0].pVertexBuffer);
		// pActiveVB->AddRef(); // Avoid memory-curruption when this is Release()ed later
		// uiStride = Xbox_g_Stream[0].Stride;
		g_pD3DDevice8->GetStreamSource(0, &pActiveVB, &uiStride);
		// retrieve stream desc
		D3DVERTEXBUFFER_DESC VBDesc;
		pActiveVB->GetDesc(&VBDesc);
		// print out stream data
		{
			printf("\n");
			printf("  Vertex Stream Data (0x%.08X)...\n", pActiveVB);
			printf("\n");
			printf("  Format : %d\n", VBDesc.Format);
			printf("  Size   : %d bytes\n", VBDesc.Size);
			printf("  FVF    : 0x%.08X\n", VBDesc.FVF);
			printf("\n");
		}

		pActiveVB->Release(); // Was absent (thus leaked memory)
#endif

		DbgDumpMesh(pIndexData, dwCount);
	}
#endif

	// perform rendering
	{
		UINT dwIndexCount = dwCount * 2; // Each DWORD data in the pushbuffer carries 2 words

		// copy index data
		{
			// remember last 2 indices
			if (dwCount >= 2) { // TODO : Is 2 indices enough for all primitive types?
				pIBMem[0] = pIndexData[dwCount - 2];
				pIBMem[1] = pIndexData[dwCount - 1];
			}
			else {
				pIBMem[0] = 0xFFFF;
			}
		}

#ifdef _DEBUG_TRACK_PB
		if (!g_PBTrackDisable.exists(pdwOrigPushData))
#endif
			// render indexed vertices
		{
			if (!g_bPBSkipPusher) {
				if (IsValidCurrentShader()) {
					// TODO: This technically should be enabled
					XTL::DxbxUpdateDeferredStates(); // CxbxUpdateNativeD3DResources

					CxbxDrawContext DrawContext = {};
					DrawContext.XboxPrimitiveType = XboxPrimitiveType;
					DrawContext.dwVertexCount = EmuD3DIndexCountToVertexCount(XboxPrimitiveType, dwIndexCount);
					DrawContext.hVertexShader = g_CurrentVertexShader;

					CxbxDrawIndexed(DrawContext, pIndexData);
				}
			}
		}
	}
}

char *DxbxXboxMethodToString(DWORD dwMethod)
{
	return "UNLABLED"; // TODO
}

extern void XTL::EmuExecutePushBufferRaw
(
    DWORD                 *pdwPushData
)
{
	static bool NV2ACallbacks_Initialized = false;
	if (!NV2ACallbacks_Initialized) {
		NV2ACallbacks_Initialized = true;
		// Set handlers that do more than just store data in registers :
		NV2ACallbacks[NV2A_NOP / 4] = EmuNV2A_NOP; // 0x0100
		NV2ACallbacks[NV2A_VERTEX_BEGIN_END / 4] = NVPB_SetBeginEnd; // NV097_SET_BEGIN_END / 0x000017FC
		NV2ACallbacks[NV2A_VB_ELEMENT_U16 / 4] = NVPB_InlineIndexArray; // NV097_ARRAY_ELEMENT16 / 0x1800
		NV2ACallbacks[(NV2A_VB_ELEMENT_U16 + 8) / 4] = NVPB_FixLoop; // NV097_ARRAY_ELEMENT32 / 0x1808
		NV2ACallbacks[NV2A_VERTEX_DATA / 4] = NVPB_InlineVertexArray; // NV097_INLINE_ARRAY / 0x1818
	}

	// Test case : XDK Sample BeginPush
	if (g_bSkipPush) {
		return;
	}

	if (pdwPushData == NULL) {
		EmuWarning("pdwPushData is null");
		return;
	}

    pdwOrigPushData = pdwPushData;

    pIndexData = NULL;
    pVertexData = NULL;

    dwVertexShader = -1;
    dwVertexStride = -1;

    // cache of last 4 indices
	pIBMem[0] = 0xFFFF;
	pIBMem[1] = 0xFFFF;
	pIBMem[2] = 0xFFFF;
	pIBMem[3] = 0xFFFF;

    XboxPrimitiveType = X_D3DPT_INVALID;

    #ifdef _DEBUG_TRACK_PB
    bShowPB = false;

    g_PBTrackTotal.insert(pdwPushData);
    if (g_PBTrackShowOnce.remove(pdwPushData) != NULL) {
        printf("\n");
        printf("\n");
        printf("  PushBuffer@0x%.08X...\n", pdwPushData);
        printf("\n");
        bShowPB = true;
    }
    #endif

	DbgPrintf("  NV2A run from 0x%.08X\n", pdwPushData);

	while (*pdwPushData != 0) {
		char LogPrefixStr[200];
		int len = sprintf(LogPrefixStr, "  NV2A Get=$%.08X", pdwPushData);

		DWORD dwPushCommand = *pdwPushData;

		// Handle jumps and/or calls :
		if (((dwPushCommand & NV2A_JMP_FLAG) > 0)
			|| ((dwPushCommand & NV2A_CALL_FLAG) > 0)) {
			// Both 'jump' and 'call' just direct execution to the indicated address :
			//pdwPushData = (DWORD*)(dwPushCommand & NV2A_ADDR_MASK);
			DbgPrintf("%s *BREAK* Jump:0x%.08X\n", LogPrefixStr, pdwPushData);

			break;// continue;
		}

		// Decode push buffer contents (inverse of D3DPUSH_ENCODE) :
		D3DPUSH_DECODE(dwPushCommand, dwMethod, dwSubCh, dwCount, bNoInc);

		// Append a counter (variable part via %d, count already formatted) :
//		if (MayLog(lfUnit))
		{
			len += sprintf(LogPrefixStr + len, " %%2d/%2d:", dwCount); // intentional %% for StepNr
			if (dwSubCh > 0)
				len += sprintf(LogPrefixStr + len, " [SubCh:%d]", dwSubCh);

			if (bNoInc)
				len += sprintf(LogPrefixStr + len, " [NoInc]");
		}

		// Skip method DWORD, remember the address of the arguments and skip over the arguments already :
		pdwPushData++;
		pdwPushArguments = pdwPushData;
		pdwPushData += dwCount;

		// Initialize handled count & name to their default :
		HandledCount = 1;
		HandledBy = nullptr;

		// Simulate writes to the NV2A instance registers; We write all DWORDs before
		// executing them so that the callbacks can read this data via the named
		// NV2AInstance fields (see for example EmuNV2A_SetVertexShaderConstant) :
		// Note this only applies to 'inc' methods (no-inc methods read all data in-place).
		if ((dwSubCh == 0) && (!bNoInc)) {
			//assert((dwMethod + (dwCount * sizeof(DWORD))) <= sizeof(NV2AInstance_Registers));

			memcpy(&(NV2AInstance_Registers[dwMethod / 4]), pdwPushArguments, dwCount * sizeof(DWORD));
		}

		// Interpret GPU Instruction(s) :
		int StepNr = 1;
		while (dwCount > 0) {
			NV2ACallback_t NV2ACallback = nullptr;

			// Skip all commands not intended for channel 0 :
			if (dwSubCh > 0) {
				HandledCount = dwCount;
				HandledBy = "*CHANNEL IGNORED*";
			}
			else {
				// Assert(dwMethod < SizeOf(NV2AInstance));

				// For 'no inc' methods, write only the first register (in case it's referenced somewhere) :
				if (bNoInc)
					NV2AInstance_Registers[dwMethod / 4]  = *pdwPushArguments;

				// Retrieve the handler callback for this method (if any) :
				NV2ACallback = NV2ACallbacks[dwMethod / sizeof(DWORD)];
			}

//#ifdef DXBX_USE_D3D
			if (g_pD3DDevice8 == nullptr) {
				HandledBy = "*NO DEVICE*"; // Don't do anything if we have no device yet (should not occur anymore, but this helps spotting errors)
				NV2ACallback = nullptr;
			}
//#endif
/*#ifdef DXBX_USE_OPENGL
			if (g_EmuWindowsDC == 0){
				HandledBy = "*NO OGL DC*"; // Don't do anything if we have no device yet (should not occur anymore, but this helps spotting errors)
				NV2ACallback = nullptr;
			}
#endif*/
			// Before handling the method, display it's details :
			if (g_bPrintfOn) {
				printf("[0x%X] ", GetCurrentThreadId());
				printf(LogPrefixStr, StepNr);
				printf(" Method=%.04X Data=%.08X %s", dwMethod, *pdwPushArguments, DxbxXboxMethodToString(dwMethod));
				if (HandledBy != nullptr) {
					printf(HandledBy);
				}

				printf("\n");
			}

			HandledBy = nullptr;
			if (NV2ACallback != nullptr) {
				NV2ACallback();
			}

			// If there are more details, print them now :
			if (HandledBy != nullptr) {
				DbgPrintf("  NV2A > %s\n", HandledBy);
			}

			// Since some instructions use less arguments, we repeat this loop
			// for the next instruction so any leftover values are handled there :
			pdwPushArguments += HandledCount;
			dwCount -= HandledCount;
			StepNr += HandledCount;
			// Re-initialize handled count & name to their default, for the next command :
			HandledCount = 1;
			HandledBy = nullptr;

			// The no-increment flag applies to method only :
			if (!bNoInc) {
				dwMethod += 4; // 1 method further
				// Remember the last two methods, in case we need to differentiate contexts (using SeenRecentMethod):
				PrevMethod[1] = PrevMethod[0];
				PrevMethod[0] = dwMethod;
			}

			// Fake a read by the Nv2A, by moving the DMA 'Get' location
			// up to where the pushbuffer is executed, so that the BusyLoop
			// in CDevice.Init finishes cleanly :
			//g_NV2ADMAChannel.Get = pdwPushData;
			// TODO : We should probably set g_NV2ADMAChannel.Put to the same value first?

			// We trigger the DMA semaphore by setting GPU time to CPU time - 2 :
			//*/*D3DDevice.*/m_pGpuTime = /*D3DDevice.*/*m_pCpuTime -2;

			// TODO : We should register vblank counts somewhere?
		} // while dwCount > 0 try
    }

	// This line is to reset the GPU 'Get' pointer, so that busyloops will terminate :
	// g_NV2ADMAChannel.Get = pdwPushData;

	// Clear the handled pushbuffer commands :
	memset(pdwOrigPushData, 0, (intptr_t)pdwPushData - (intptr_t)pdwOrigPushData);

    #ifdef _DEBUG_TRACK_PB
    if (bShowPB) {
        printf("\n");
        printf("CxbxDbg> ");
        fflush(stdout);
    }
    #endif

    if (g_bStepPush) {
        g_pD3DDevice8->Present(0,0,0,0);
        Sleep(500);
    }
}

#ifdef _DEBUG_TRACK_PB
void DbgDumpMesh(XTL::INDEX16 *pIndexData, DWORD dwCount)
{
    if (!XTL::IsValidCurrentShader() || (dwCount == 0))
        return;

    // retrieve stream data
    char szFileName[128];
    sprintf(szFileName, "D:\\_cxbx\\mesh\\CxbxMesh-0x%.08X.x", pIndexData);
    FILE *dbgVertices = fopen(szFileName, "wt");

    BYTE *pVBData = (BYTE *)XTL::GetDataFromXboxResource(XTL::Xbox_g_Stream[0].pVertexBuffer);
    UINT  uiStride = XTL::Xbox_g_Stream[0].Stride;

    // print out stream data
    {
        XTL::INDEX16 maxIndex = 0;
		XTL::INDEX16 *pIndexCheck = pIndexData;
        for(uint chk=0;chk<dwCount;chk++) {
			XTL::INDEX16 x = *pIndexCheck++;
            if (x > maxIndex)
                maxIndex = x;
        }
#if 0
        if (maxIndex > ((VBDesc.Size/uiStride) - 1))
            maxIndex = (VBDesc.Size / uiStride) - 1;
#endif
        fprintf(dbgVertices, "xof 0303txt 0032\n");
        fprintf(dbgVertices, "\n");
        fprintf(dbgVertices, "//\n");
        fprintf(dbgVertices, "//  Vertex Stream Data (0x%.08X)...\n", pVBData);
        fprintf(dbgVertices, "//\n");
#if 0
		fprintf(dbgVertices, "//  Format : %d\n", VBDesc.Format);
        fprintf(dbgVertices, "//  Size   : %d bytes\n", VBDesc.Size);
        fprintf(dbgVertices, "//  FVF    : 0x%.08X\n", VBDesc.FVF);
#endif
        fprintf(dbgVertices, "//  iCount : %d\n", dwCount/2);
        fprintf(dbgVertices, "//\n");
		fprintf(dbgVertices, "\n");
        fprintf(dbgVertices, "Frame SCENE_ROOT {\n");
        fprintf(dbgVertices, "\n");
        fprintf(dbgVertices, "  FrameTransformMatrix {\n");
        fprintf(dbgVertices, "    1.000000,0.000000,0.000000,0.000000,\n");
        fprintf(dbgVertices, "    0.000000,1.000000,0.000000,0.000000,\n");
        fprintf(dbgVertices, "    0.000000,0.000000,1.000000,0.000000,\n");
        fprintf(dbgVertices, "    0.000000,0.000000,0.000000,1.000000;;\n");
        fprintf(dbgVertices, "  }\n");
        fprintf(dbgVertices, "\n");
        fprintf(dbgVertices, "  Frame Turok1 {\n");
        fprintf(dbgVertices, "\n");
        fprintf(dbgVertices, "    FrameTransformMatrix {\n");
        fprintf(dbgVertices, "      1.000000,0.000000,0.000000,0.000000,\n");
        fprintf(dbgVertices, "      0.000000,1.000000,0.000000,0.000000,\n");
        fprintf(dbgVertices, "      0.000000,0.000000,1.000000,0.000000,\n");
        fprintf(dbgVertices, "      0.000000,0.000000,0.000000,1.000000;;\n");
        fprintf(dbgVertices, "    }\n");
        fprintf(dbgVertices, "\n");
        fprintf(dbgVertices, "    Mesh {\n");
        fprintf(dbgVertices, "      %d;\n", maxIndex+1);

        uint max = maxIndex+1;
        for(uint v=0;v<max;v++)
        {
            fprintf(dbgVertices, "      %f;%f;%f;%s\n",
                *(FLOAT*)&pVBData[v*uiStride+0],
                *(FLOAT*)&pVBData[v*uiStride+4],
                *(FLOAT*)&pVBData[v*uiStride+8],
                (v < (max - 1)) ? "," : ";");
        }

        fprintf(dbgVertices, "      %d;\n", dwCount - 2);

		XTL::INDEX16 *pIndexValues = pIndexData;

        max = dwCount;

        DWORD a = *pIndexValues++;
        DWORD b = *pIndexValues++;
        DWORD c = *pIndexValues++;

        DWORD la = a,lb = b,lc = c;

        for(uint i=2;i<max;i++)
        {
            fprintf(dbgVertices, "      3;%d,%d,%d;%s\n",
                a,b,c, (i < (max - 1)) ? "," : ";");

            a = b;
            b = c;
            c = *pIndexValues++;

            la = a;
            lb = b;
            lc = c;
        }

        fprintf(dbgVertices, "    }\n");
        fprintf(dbgVertices, "  }\n");
        fprintf(dbgVertices, "}\n");

        fclose(dbgVertices);
    }
}

void XTL::DbgDumpPushBuffer(DWORD* PBData, DWORD dwSize)
{
	static int PbNumber = 0;	// Keep track of how many push buffers we've attemted to convert.
	DWORD dwVertexShader;
	char szPB[512];

	// Prevent dumping too many of these!
	if (PbNumber > 300) {
		return;
	}

	// Get a copy of the current vertex shader
	g_pD3DDevice8->GetVertexShader(&dwVertexShader);

	/*if (g_CurrentVertexShader != dwVertexShader) {
		printf( "g_CurrentVertexShader does not match FVF from GetVertexShader!\n"
					"g_CurrentVertexShader = 0x%.08X\n"
					"GetVertexShader = 0x%.08X\n" );
	}*/

	if (dwVertexShader > 0xFFFF) {
		EmuWarning("Cannot dump pushbuffer without an FVF (programmable shaders not supported)");
		return;
	}

	sprintf(szPB, "D:\\cxbx\\_pushbuffer\\pushbuffer%.03d.txt", PbNumber++);
	// Create a new file for this pushbuffer's data
	HANDLE hFile = CreateFile(szPB, GENERIC_WRITE, 0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, 0);
	if (hFile == INVALID_HANDLE_VALUE) {
		EmuWarning("Error creating pushbuffer file!");
	}

	DWORD dwBytesWritten;

	// Write pushbuffer data to the file.
	// TODO: Cache the 32-bit XXHash32::hash() of each pushbuffer to ensure that the same
	// pushbuffer is not written twice within a given emulation session.
	WriteFile(hFile, &g_CurrentVertexShader, sizeof(DWORD), &dwBytesWritten, nullptr);
	WriteFile(hFile, PBData, dwSize, &dwBytesWritten, nullptr);
	// Close handle
	CloseHandle(hFile);
}

#endif