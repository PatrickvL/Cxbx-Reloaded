// ******************************************************************
// *
// *    .,-:::::    .,::      .::::::::.    .,::      .:
// *  ,;;;'````'    `;;;,  .,;;  ;;;'';;'   `;;;,  .,;;
// *  [[[             '[[,,[['   [[[__[[\.    '[[,,[['
// *  $$$              Y$$$P     $$""""Y$$     Y$$$P
// *  `88bo,__,o,    oP"``"Yo,  _88o,,od8P   oP"``"Yo,
// *    "YUMMMMMP",m"       "Mm,""YUMMMP" ,m"       "Mm,
// *
// *   Cxbx->Win32->CxbxKrnl->D3D8.OOVPA.inl
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

#ifndef XONLINE_OOVPA_INL
#define XONLINE_OOVPA_INL

#include "OOVPA.h"

#include "HLEDataBase/D3D8.1.0.3925.inl"
#include "HLEDataBase/D3D8.1.0.4034.inl"
#include "HLEDataBase/D3D8.1.0.4134.inl"
#include "HLEDataBase/D3D8.1.0.4361.inl"
#include "HLEDataBase/D3D8.1.0.4432.inl"
#include "HLEDataBase/D3D8.1.0.4627.inl"
#include "HLEDataBase/D3D8.1.0.5028.inl"
#include "HLEDataBase/D3D8.1.0.5233.inl"
#include "HLEDataBase/D3D8.1.0.5344.inl"
#include "HLEDataBase/D3D8.1.0.5558.inl"
#include "HLEDataBase/D3D8.1.0.5788.inl"
#include "HLEDataBase/D3D8.1.0.5849.inl"

// ******************************************************************
// * D3D8_OOVPA
// ******************************************************************
OOVPATable D3D8_OOVPA[] = {

	REGISTER_OOVPAS(CMiniport_InitHardware, 3911, 5558), // PATCH 
	REGISTER_OOVPAS(CMiniport_CreateCtxDmaObject, 3911, 4361), // PATCH 
	REGISTER_OOVPAS(D3D_CMiniport_GetDisplayCapabilities, 3911, 3925, 5849), // DISABLED Was D3DDevice_Unknown1 PATCH
	REGISTER_OOVPAS(D3D_CheckDeviceFormat, 4034), // DISABLED
	REGISTER_OOVPAS(Direct3D_CreateDevice, 3911), // PATCH Was 3925, 4627
	REGISTER_OOVPAS(D3DDevice_IsBusy, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetDeviceCaps, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_BeginVisibilityTest, 3925, 4034, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_EndVisibilityTest, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVisibilityTestResult, 3925, 5344), // PATCH 
	REGISTER_OOVPAS(D3D_KickOffAndWaitForIdle, 3925, 4034, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetMaterial, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetBackMaterial, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetBackMaterial, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_LoadVertexShader, 3925, 4034, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SelectVertexShader, 3925, 4039, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_Release, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_BlockUntilVerticalBlank, 3925, 4039, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVerticalBlankCallback, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderTarget, 3925, 3948, 4039, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetStreamSource, 3925, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexShader, 3925, 4034, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreatePixelShader, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetPixelShader, 3925, 4039, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetIndices, 3925, 4034, 5558), // DISABLED
	REGISTER_OOVPAS(D3DDevice_SetViewport, 3925, 4034, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetTexture2, 3911, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTexture, 3911, 4034, 5344), // PATCH Was 3925, 5233
	REGISTER_OOVPAS(D3DDevice_Begin, 3925, 4039, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexData2f, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexData2s, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexData4f, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexData4ub, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexDataColor, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_End, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_Clear, 3925, 4034, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreatePalette, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_CreatePalette2, 5558), // DISABLED
	REGISTER_OOVPAS(D3DDevice_SetPalette, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreateTexture, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_CreateTexture2, 5788), // DISABLED
	REGISTER_OOVPAS(D3DDevice_CreateVolumeTexture, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_CreateCubeTexture, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_CreateIndexBuffer, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_SetVertexShaderConstant, 3925, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexShaderConstant1, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexShaderConstant4, 4627), // PATCH 
	REGISTER_OOVPA_ALIAS(D3DDevice_SetVertexShaderConstant1, 5558, D3DDevice_SetVertexShaderConstant1Fast),
	REGISTER_OOVPAS(D3DDevice_SetFlickerFilter, 3925, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetSoftDisplayFilter, 3925, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTextureState_TexCoordIndex, 3925, 4034, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTextureState_BorderColor, 3925, 4034, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_PSTextureModes, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_StencilFail, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_CullMode, 3925, 4034, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_Simple, 3925, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetTransform, 3925, 4039, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTransform, 3925, 4034, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_FogColor, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_FillMode, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_StencilEnable, 3925, 4034, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_Dxt1NoiseEnable, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_ZBias, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_ZEnable, 3925, 4034, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_Present, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_Swap, 4034, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetShaderConstantMode, 3925, 4039, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetBackBuffer, 3925, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetBackBuffer2, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetRenderTarget, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetDepthStencilSurface, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetDepthStencilSurface2, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreateVertexBuffer, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_CreateVertexBuffer2, 5558), // DISABLED
	REGISTER_OOVPAS(D3DVertexBuffer_Lock, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DResource_Register, 3925), // PATCH 
	REGISTER_OOVPAS(D3DResource_Release, 3925), // PATCH 
	REGISTER_OOVPAS(D3DResource_AddRef, 3925), // PATCH 
	REGISTER_OOVPAS(D3DResource_IsBusy, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DBaseTexture_GetLevelCount, 4361), // PATCH 
	REGISTER_OOVPAS(D3DSurface_LockRect, 3925), // PATCH 
	REGISTER_OOVPAS(D3DPalette_Lock, 3925), // PATCH 
	REGISTER_OOVPAS(D3DPalette_Lock2, 4627), // PATCH 
	REGISTER_OOVPAS(D3DTexture_LockRect, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DVolumeTexture_LockBox, 3925), // PATCH 
	REGISTER_OOVPAS(D3DCubeTexture_LockRect, 3925), // PATCH 
	REGISTER_OOVPAS(D3DCubeTexture_GetCubeMapSurface2, 4627), // PATCH 
	REGISTER_OOVPAS(D3DTexture_GetSurfaceLevel, 3925), // PATCH 
	REGISTER_OOVPAS(Lock2DSurface, 3925), // PATCH 
	REGISTER_OOVPAS(Get2DSurfaceDesc, 3925, 4034, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShaderSize, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetGammaRamp, 3925, 4928), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetMaterial, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_AddRef, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetViewport, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetGammaRamp, 3925, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetDisplayFieldStatus, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetCreationParameters, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_MultiSampleAntiAlias, 3925, 4034, 4627, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_VertexBlend, 3925, 4034), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_BackFillMode, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_TwoSidedLighting, 3925, 4034, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_NormalizeNormals, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_FrontFace, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_TextureFactor, 3925, 4034, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_LogicOp, 3925, 4034, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_EdgeAntiAlias, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_MultiSampleMask, 3925, 4034, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_MultiSampleMode, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_PersistDisplay, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_ShadowFunc, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_LineWidth, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_YuvEnable, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_OcclusionCullEnable, 3925, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_StencilCullEnable, 3925, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DrawVertices, 3925, 4034, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DrawVerticesUP, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DrawIndexedVertices, 3925, 4034, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DrawIndexedVerticesUP, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DrawRectPatch, 3911), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DrawTriPatch, 3911), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetDisplayMode, 3925, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTextureState_BumpEnv, 3925, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTextureState_ColorKeyColor, 3925, 4034, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexData4s, 3925), // PATCH 
	REGISTER_OOVPAS(D3D_SetPushBufferSize, 3925), // DISABLED
	REGISTER_OOVPAS(D3DResource_GetType, 3925), // DISABLED Probably not even necessary...
	REGISTER_OOVPAS(D3D_AllocContiguousMemory, 3925, 5788), // DISABLED Just calls MmAllocateContiguousMemory. Was PATCH
	REGISTER_OOVPAS(D3DDevice_SetRenderState_Deferred, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_GetLight, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetLight, 3925, 4034, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_LightEnable, 3911, 4034, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetLightEnable, 3911, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreateVertexShader, 3925), // PATCH 
	REGISTER_OOVPAS(D3DSurface_GetDesc, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_GetProjectionViewportMatrix, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetTile, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_ApplyStateBlock, 3925, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CaptureStateBlock, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DeleteStateBlock, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreateStateBlock, 3925, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DeletePixelShader, 3925, 5344), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetPixelShaderProgram, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_KickOff, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetModelView, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_FlushVertexCache, 3925, 4134), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetScissors, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexShaderInput, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_PrimeVertexCache, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetPixelShaderConstant, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_InsertCallback, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_BeginPushBuffer, 3925, 5788), // PATCH Not implemented yet.
	REGISTER_OOVPAS(D3DDevice_EndPushBuffer, 3925, 5788), // PATCH  Not implemented yet.
	REGISTER_OOVPAS(D3DDevice_SetRenderState_RopZCmpAlwaysRead, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_RopZRead, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_DoNotCullUncompressed, 3925), // PATCH 
	REGISTER_OOVPAS(XMETAL_StartPush, 3925), // DISABLED
	REGISTER_OOVPAS(D3D_SetFence, 3925, 5849), // XREF
	REGISTER_OOVPAS(D3DDevice_InsertFence, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_LoadVertexShaderProgram, 3925, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DeleteVertexShader, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_RunPushBuffer, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_Reset, 3925, 4134), // PATCH 
	REGISTER_OOVPAS(D3D_GetAdapterIdentifier, 3925), // DISABLED
	REGISTER_OOVPAS(D3D_GetDeviceCaps, 3925), // DISABLED
	REGISTER_OOVPAS(D3D_SetCommonDebugRegisters, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CreateImageSurface, 3925, 4034), // DISABLED
	REGISTER_OOVPAS(D3D_BlockOnTime, 3925, 5849), // PATCH Was XREF
	REGISTER_OOVPAS(D3D_BlockOnResource, 3925), // PATCH 
	REGISTER_OOVPAS(D3DResource_BlockUntilNotBusy, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_BlockOnFence, 3925, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_BeginStateBlock, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_EndStateBlock, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetTile, 3925, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SwitchTexture, 3925, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetModelView, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_CopyRects, 3925, 4034, 5233), // PATCH 
	REGISTER_OOVPAS(D3DVertexBuffer_GetDesc, 3925), // DISABLED
	REGISTER_OOVPAS(D3DDevice_GetShaderConstantMode, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShader, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShaderConstant, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShaderInput, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_RunVertexStateShader, 3925, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShaderType, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShaderDeclaration, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetVertexShaderFunction, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetPixelShader, 3925), // PATCH 
	REGISTER_OOVPAS(D3DDevice_IsFencePending, 3925, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_DeletePatch, 4034), // DISABLED (TODO)
	REGISTER_OOVPAS(D3DDevice_SetBackBufferScale, 5788), // PATCH 
	REGISTER_OOVPAS(MakeRequestedSpace, 5788), // PATCH NOT VERIFIED
	REGISTER_OOVPAS(D3DDevice_BeginPush, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_BeginStateBig, 5788), // XREF 
	REGISTER_OOVPAS(D3DDevice_CreateIndexBuffer2, 5558), // DISABLED
	REGISTER_OOVPAS(D3DDevice_EnableOverlay, 4361), // PATCH 
	REGISTER_OOVPAS(D3DDevice_EndPush, 4627), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetOverlayUpdateStatus, 5558), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetPushBufferOffset, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetRenderTarget2, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetScissors, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_GetViewportOffsetAndScale, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_LazySetStateVB, 5558), // XREF
	REGISTER_OOVPAS(D3DDevice_SetDepthClipPlanes, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_MultiSampleRenderTargetMode, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderState_SampleAlpha, 5028), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetRenderTargetFast, 5344), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetScreenSpaceOffset, 5849), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetSwapCallback, 5788), // PATCH 
	REGISTER_OOVPAS(D3DDevice_SetVertexShaderConstantNotInline, 5233), // PATCH 
	REGISTER_OOVPAS(D3DDevice_UpdateOverlay, 5558), // PATCH 
	REGISTER_OOVPAS(D3DTexture_GetSurfaceLevel2, 4627), // PATCH 
	REGISTER_OOVPAS(D3DVertexBuffer_Lock2, 5788), // PATCH 
	REGISTER_OOVPAS(D3D_ClearStateBlockFlags, 5788), // XREF
	REGISTER_OOVPAS(D3D_RecordStateBlock, 5788), // XREF
	REGISTER_OOVPAS(D3D_Unknown, 5788), // XREF
};

// ******************************************************************
// * D3D8_OOVPA_SIZE
// ******************************************************************
uint32 D3D8_OOVPA_SIZE = sizeof(D3D8_OOVPA);

#endif