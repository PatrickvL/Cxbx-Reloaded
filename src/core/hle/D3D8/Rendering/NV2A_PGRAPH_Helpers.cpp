// ******************************************************************
// *  NV2A PGRAPH Helper Functions — Implementation
// *
// *  Structured accessors for PGRAPH register state.
// *  All functions read directly from NV2A hardware state (no caching).
// *
// *  Primary overloads accept NV2AState*; convenience overloads (no arg)
// *  go through the g_NV2A global and are defined inline below.
// ******************************************************************
#include "EmuD3D8_common.h"
#include "NV2A_PGRAPH_Helpers.h"

// ==== Primary overloads (accept NV2AState*) ====

uint32_t NV2AGetTextureControlRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];
}

uint32_t NV2AGetTextureControl1Raw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXCTL1_0 + stage * 4)];
}

uint32_t NV2AGetTextureOffsetRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXOFFSET0 + stage * 4)];
}

uint32_t NV2AGetTextureAddressModeRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXADDRESS0 + stage * 4)];
}

uint32_t NV2AGetTextureFormatRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXFMT0 + stage * 4)];
}

uint32_t NV2AGetTextureImageRectRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXIMAGERECT0 + stage * 4)];
}

uint32_t NV2AGetTextureFilterRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXFILTER0 + stage * 4)];
}

uint32_t NV2AGetBorderColorRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_BORDERCOLOR0 + stage * 4)];
}

uint32_t NV2AGetTexturePaletteRaw(NV2AState* d, int stage)
{
	auto pg = &d->pgraph;
	return pg->regs[RI(NV_PGRAPH_TEXPALETTE0 + stage * 4)];
}

// Derived from Raw PGRAPH register values:

bool NV2AIsTextureEnabled(NV2AState* d, int stage)
{
	return (NV2AGetTextureControlRaw(d, stage) & NV_PGRAPH_TEXCTL0_0_ENABLE) != 0;
}

uint32_t NV2AGetTexturePitch(NV2AState* d, int stage)
{
	return GET_MASK(NV2AGetTextureControl1Raw(d, stage), NV_PGRAPH_TEXCTL1_0_IMAGE_PITCH);
}

uint32_t NV2AResolveTexturePhysicalAddress(NV2AState* d, int stage, uint32_t rawOffset)
{
	auto pg = &d->pgraph;
	uint32_t texFmt = NV2AGetTextureFormatRaw(d, stage);
	bool dmaSelect = (texFmt & NV_PGRAPH_TEXFMT0_CONTEXT_DMA) != 0;
	return pg->dma_base[dmaSelect] + rawOffset;
}

uint32_t NV2AResolvePalettePhysicalAddress(NV2AState* d, int stage, uint32_t rawOffset)
{
	auto pg = &d->pgraph;
	uint32_t texPalette = NV2AGetTexturePaletteRaw(d, stage);
	bool dmaSelect = (texPalette & NV_PGRAPH_TEXPALETTE0_CONTEXT_DMA) != 0;
	return pg->dma_base[dmaSelect] + rawOffset;
}

NV2ATextureAddress NV2AGetTextureAddress(NV2AState* d, int stage)
{
	NV2ATextureAddress result = {};

	uint32_t rawOffset = NV2AGetTextureOffsetRaw(d, stage);
	result.rawOffset = rawOffset;

	if (rawOffset == 0)
		return result;

	result.physicalAddress = NV2AResolveTexturePhysicalAddress(d, stage, rawOffset);
	result.dmaBase = result.physicalAddress - rawOffset;
	return result;
}

NV2ATextureFormat NV2AGetTextureFormat(NV2AState* d, int stage)
{
	NV2ATextureFormat result = {};

	uint32_t fmt = NV2AGetTextureFormatRaw(d, stage);
	result.raw = fmt;
	result.dmaSelect = (fmt & NV_PGRAPH_TEXFMT0_CONTEXT_DMA) != 0;
	result.cubemap = (fmt & NV_PGRAPH_TEXFMT0_CUBEMAPENABLE) != 0;
	result.borderFromColor = (fmt & NV_PGRAPH_TEXFMT0_BORDER_SOURCE) != 0;
	result.dimensionality = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_DIMENSIONALITY);
	result.colorFormat = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_COLOR);
	result.mipLevels = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_MIPMAP_LEVELS);
	result.logWidth = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_BASE_SIZE_U);
	result.logHeight = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_BASE_SIZE_V);
	result.logDepth = GET_MASK(fmt, NV_PGRAPH_TEXFMT0_BASE_SIZE_P);
	return result;
}

NV2ATextureImageRect NV2AGetTextureImageRect(NV2AState* d, int stage)
{
	NV2ATextureImageRect result = {};

	uint32_t rect = NV2AGetTextureImageRectRaw(d, stage);
	result.width = GET_MASK(rect, NV_PGRAPH_TEXIMAGERECT0_WIDTH);
	result.height = GET_MASK(rect, NV_PGRAPH_TEXIMAGERECT0_HEIGHT);
	return result;
}

NV2ATextureControl NV2AGetTextureControl(NV2AState* d, int stage)
{
	NV2ATextureControl result = {};

	uint32_t texctl0 = NV2AGetTextureControlRaw(d, stage);

	result.enabled = NV2AIsTextureEnabled(d, stage);
	result.minLodClamp = GET_MASK(texctl0, NV_PGRAPH_TEXCTL0_0_MIN_LOD_CLAMP);
	result.maxLodClamp = GET_MASK(texctl0, NV_PGRAPH_TEXCTL0_0_MAX_LOD_CLAMP);
	result.pitch = NV2AGetTexturePitch(d, stage);
	return result;
}

NV2APaletteState NV2AGetPaletteState(NV2AState* d, int stage)
{
	NV2APaletteState result = {};

	uint32_t texPalette = NV2AGetTexturePaletteRaw(d, stage);
	uint32_t palOffset = texPalette & NV_PGRAPH_TEXPALETTE0_OFFSET;
	if (palOffset == 0)
		return result;

	uint32_t lengthField = GET_MASK(texPalette, NV_PGRAPH_TEXPALETTE0_LENGTH);
	static const unsigned palSizes[] = { 256 * 4, 128 * 4, 64 * 4, 32 * 4 };

	result.physicalAddress = NV2AResolvePalettePhysicalAddress(d, stage, palOffset);
	result.data = (void*)(CONTIGUOUS_MEMORY_BASE + result.physicalAddress);
	result.size = palSizes[lengthField & 3];
	return result;
}

NV2ASurfaceState NV2AGetSurfaceState(PGRAPHState* pg)
{
	NV2ASurfaceState result = {};

	result.colorOffset = pg->regs[RI(NV_PGRAPH_BOFFSET3)];
	result.zetaOffset = pg->regs[RI(NV_PGRAPH_BOFFSET4)];

	uint32_t dmaPitch = pg->regs[RI(NV_PGRAPH_DMA_PITCH)];
	result.colorPitch = GET_MASK(dmaPitch, NV_PGRAPH_DMA_PITCH_COLOR);
	result.zetaPitch = GET_MASK(dmaPitch, NV_PGRAPH_DMA_PITCH_ZETA);

	uint32_t clipX = pg->regs[RI(NV_PGRAPH_SURFACECLIPX)];
	result.clipX = GET_MASK(clipX, NV_PGRAPH_SURFACECLIPX_X);
	result.clipWidth = GET_MASK(clipX, NV_PGRAPH_SURFACECLIPX_WIDTH);

	uint32_t clipY = pg->regs[RI(NV_PGRAPH_SURFACECLIPY)];
	result.clipY = GET_MASK(clipY, NV_PGRAPH_SURFACECLIPY_Y);
	result.clipHeight = GET_MASK(clipY, NV_PGRAPH_SURFACECLIPY_HEIGHT);

	uint32_t fmt = pg->regs[RI(NV_PGRAPH_SURFACEFORMAT)];
	result.colorFormat = GET_MASK(fmt, NV_PGRAPH_SURFACEFORMAT_COLOR);
	result.zetaFormat = GET_MASK(fmt, NV_PGRAPH_SURFACEFORMAT_ZETA);
	result.surfaceType = GET_MASK(fmt, NV_PGRAPH_SURFACEFORMAT_TYPE);
	result.antiAliasing = GET_MASK(fmt, NV_PGRAPH_SURFACEFORMAT_ANTI_ALIASING);
	result.logWidth = GET_MASK(fmt, NV_PGRAPH_SURFACEFORMAT_WIDTH);
	result.logHeight = GET_MASK(fmt, NV_PGRAPH_SURFACEFORMAT_HEIGHT);
	return result;
}

NV2ASurfaceState NV2AGetSurfaceState(NV2AState* d)
{
	return NV2AGetSurfaceState(&d->pgraph);
}

// ==== Vertex Shader Mode ====

bool NV2AIsFixedFunctionMode(PGRAPHState* pg)
{
	uint32_t mode = GET_MASK(pg->regs[RI(NV_PGRAPH_CSV0_D)], NV_PGRAPH_CSV0_D_MODE);
	return mode != NV097_SET_TRANSFORM_EXECUTION_MODE_MODE_PROGRAM;
}

bool NV2AIsFixedFunctionMode(NV2AState* d)
{
	return NV2AIsFixedFunctionMode(&d->pgraph);
}

bool NV2AIsFixedFunctionMode()
{
	return NV2AIsFixedFunctionMode(&g_NV2A->GetDeviceState()->pgraph);
}

// ==== Convenience overloads (go through g_NV2A global) ====

uint32_t NV2AGetTextureControlRaw(int stage) { return NV2AGetTextureControlRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTextureControl1Raw(int stage) { return NV2AGetTextureControl1Raw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTextureOffsetRaw(int stage) { return NV2AGetTextureOffsetRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTextureAddressModeRaw(int stage) { return NV2AGetTextureAddressModeRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTextureFormatRaw(int stage) { return NV2AGetTextureFormatRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTextureImageRectRaw(int stage) { return NV2AGetTextureImageRectRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTextureFilterRaw(int stage) { return NV2AGetTextureFilterRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetBorderColorRaw(int stage) { return NV2AGetBorderColorRaw(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTexturePaletteRaw(int stage) { return NV2AGetTexturePaletteRaw(g_NV2A->GetDeviceState(), stage); }

uint32_t NV2AResolveTexturePhysicalAddress(int stage, uint32_t rawOffset) { return NV2AResolveTexturePhysicalAddress(g_NV2A->GetDeviceState(), stage, rawOffset); }
uint32_t NV2AResolvePalettePhysicalAddress(int stage, uint32_t rawOffset) { return NV2AResolvePalettePhysicalAddress(g_NV2A->GetDeviceState(), stage, rawOffset); }

bool NV2AIsTextureEnabled(int stage) { return NV2AIsTextureEnabled(g_NV2A->GetDeviceState(), stage); }
uint32_t NV2AGetTexturePitch(int stage) { return NV2AGetTexturePitch(g_NV2A->GetDeviceState(), stage); }

NV2ATextureAddress NV2AGetTextureAddress(int stage) { return NV2AGetTextureAddress(g_NV2A->GetDeviceState(), stage); }
NV2ATextureFormat NV2AGetTextureFormat(int stage) { return NV2AGetTextureFormat(g_NV2A->GetDeviceState(), stage); }
NV2ATextureImageRect NV2AGetTextureImageRect(int stage) { return NV2AGetTextureImageRect(g_NV2A->GetDeviceState(), stage); }
NV2ATextureControl NV2AGetTextureControl(int stage) { return NV2AGetTextureControl(g_NV2A->GetDeviceState(), stage); }
NV2APaletteState NV2AGetPaletteState(int stage) { return NV2AGetPaletteState(g_NV2A->GetDeviceState(), stage); }
NV2ASurfaceState NV2AGetSurfaceState() { return NV2AGetSurfaceState(g_NV2A->GetDeviceState()); }
