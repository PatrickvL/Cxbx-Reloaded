// ******************************************************************
// *  NV2A PGRAPH Helper Functions
// *
// *  Structured accessors for PGRAPH register state.
// *  All functions read directly from NV2A hardware state (no caching).
// *  Functions that resolve physical addresses incorporate DMA context.
// *
// *  Each function has two overloads:
// *    - One accepting NV2AState* (for callers that already have the device state)
// *    - One with no device arg (goes through the g_NV2A global)
// ******************************************************************
#pragma once

#include "devices\video\nv2a_int.h" // For NV2AState, PGRAPHState, RI(), NV_PGRAPH_* constants

// ---- Texture Stage Scalars (inline — single register read) ----

inline uint32_t NV2AGetTextureControlRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXCTL0_0 + stage * 4)];
}

inline uint32_t NV2AGetTextureControl1Raw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXCTL1_0 + stage * 4)];
}

inline uint32_t NV2AGetTextureOffsetRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXOFFSET0 + stage * 4)];
}

inline uint32_t NV2AGetTextureAddressModeRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXADDRESS0 + stage * 4)];
}

inline uint32_t NV2AGetTextureFormatRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXFMT0 + stage * 4)];
}

inline uint32_t NV2AGetTextureImageRectRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXIMAGERECT0 + stage * 4)];
}

inline uint32_t NV2AGetTextureFilterRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXFILTER0 + stage * 4)];
}

inline uint32_t NV2AGetBorderColorRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_BORDERCOLOR0 + stage * 4)];
}

inline uint32_t NV2AGetTexturePaletteRaw(NV2AState* d, int stage) {
	return d->pgraph.regs[RI(NV_PGRAPH_TEXPALETTE0 + stage * 4)];
}

// Derived from Raw PGRAPH register values:

inline bool NV2AIsTextureEnabled(NV2AState* d, int stage) {
	return (NV2AGetTextureControlRaw(d, stage) & NV_PGRAPH_TEXCTL0_0_ENABLE) != 0;
}

inline uint32_t NV2AGetTexturePitch(NV2AState* d, int stage) {
	return GET_MASK(NV2AGetTextureControl1Raw(d, stage), NV_PGRAPH_TEXCTL1_0_IMAGE_PITCH);
}

// Convenience overloads (go through g_NV2A global; defined in .cpp)
uint32_t NV2AGetTextureControlRaw(int stage);
uint32_t NV2AGetTextureControl1Raw(int stage);
uint32_t NV2AGetTextureOffsetRaw(int stage);
uint32_t NV2AGetTextureAddressModeRaw(int stage);
uint32_t NV2AGetTextureFormatRaw(int stage);
uint32_t NV2AGetTextureImageRectRaw(int stage);
uint32_t NV2AGetTextureFilterRaw(int stage);
uint32_t NV2AGetBorderColorRaw(int stage);
uint32_t NV2AGetTexturePaletteRaw(int stage);
bool NV2AIsTextureEnabled(int stage);
uint32_t NV2AGetTexturePitch(int stage);

// ---- DMA Resolution ----

// Resolve a raw VRAM offset through the texture DMA context for a given stage.
// This adds dma_base[0 or 1] (per TEXFMT CONTEXT_DMA bit) to the offset.
uint32_t NV2AResolveTexturePhysicalAddress(NV2AState* d, int stage, uint32_t rawOffset);
uint32_t NV2AResolveTexturePhysicalAddress(int stage, uint32_t rawOffset);

// Resolve a raw VRAM offset through the palette DMA context for a given stage.
// This adds dma_base[0 or 1] (per TEXPALETTE CONTEXT_DMA bit) to the offset.
uint32_t NV2AResolvePalettePhysicalAddress(NV2AState* d, int stage, uint32_t rawOffset);
uint32_t NV2AResolvePalettePhysicalAddress(int stage, uint32_t rawOffset);

// Resolve a raw vertex buffer offset through the vertex DMA context.
// This adds dma_vertex_base[dmaSelect] to the offset, masked to 27-bit VRAM range.
inline uint32_t NV2AResolveVertexPhysicalAddress(NV2AState* d, bool dmaSelect, uint32_t rawOffset) {
	return (d->pgraph.dma_vertex_base[dmaSelect] + rawOffset) & 0x07FFFFFF;
}
uint32_t NV2AResolveVertexPhysicalAddress(bool dmaSelect, uint32_t rawOffset);

// ---- Texture Stage State ----

struct NV2ATextureFormat {
	uint32_t raw;             // Raw TEXFMT register value (= Xbox D3D Format field)
	uint32_t colorFormat;     // NV097_SET_TEXTURE_FORMAT_COLOR_* value
	uint32_t mipLevels;       // Number of mipmap levels
	uint32_t logWidth;        // log2(base width)
	uint32_t logHeight;       // log2(base height)
	uint32_t logDepth;        // log2(base depth)
	uint32_t dimensionality;  // 1, 2, or 3
	bool     cubemap;         // Cubemap enabled
	bool     borderFromColor; // Border source is color (vs texture data)
	bool     dmaSelect;       // CONTEXT_DMA bit (false=A, true=B)
};

// Read decoded texture format for a stage.
NV2ATextureFormat NV2AGetTextureFormat(NV2AState* d, int stage);
NV2ATextureFormat NV2AGetTextureFormat(int stage);

struct NV2ATextureImageRect {
	uint32_t width;           // Linear texture width (from TEXIMAGERECT)
	uint32_t height;          // Linear texture height
};

// Read linear texture dimensions (TEXIMAGERECT register).
NV2ATextureImageRect NV2AGetTextureImageRect(NV2AState* d, int stage);
NV2ATextureImageRect NV2AGetTextureImageRect(int stage);

struct NV2ATextureControl {
	bool     enabled;         // Texture stage enabled (TEXCTL0 bit 30)
	uint32_t pitch;           // Image pitch from TEXCTL1 (bits 16-31)
	uint32_t minLodClamp;     // Min LOD clamp (TEXCTL0 bits 18-29)
	uint32_t maxLodClamp;     // Max LOD clamp (TEXCTL0 bits 6-17)
};

// Read texture control state for a stage.
NV2ATextureControl NV2AGetTextureControl(NV2AState* d, int stage);
NV2ATextureControl NV2AGetTextureControl(int stage);

// ---- Palette Stage State ----

struct NV2APaletteState {
	void*    data;    // Pointer to palette data (CONTIGUOUS_MEMORY_BASE + physAddr), or nullptr
	unsigned size;    // Size in bytes (entries × 4)
	uint32_t physicalAddress; // DMA-resolved physical address
};

// Resolve palette data pointer and size for a texture stage.
// Reads TEXPALETTE0, resolves through DMA context. Returns empty if offset is zero.
// Caller is responsible for checking whether the texture format is palettized.
NV2APaletteState NV2AGetPaletteState(NV2AState* d, int stage);
NV2APaletteState NV2AGetPaletteState(int stage);

// ---- Surface State ----
// NV2A surfaces use the same PGRAPH register file as textures.
// NV2ASurfaceState is defined in nv2a_int.h (member of PGRAPHState).

// Read current surface configuration from PGRAPH.
NV2ASurfaceState NV2AGetSurfaceState(NV2AState* d);
NV2ASurfaceState NV2AGetSurfaceState(PGRAPHState* pg);
NV2ASurfaceState NV2AGetSurfaceState();

// ---- Vertex Shader Mode ----

// Returns true when CSV0_D_MODE == FIXED (hardware T&L pipeline).
// Returns false when CSV0_D_MODE == PROGRAM (user VS program in program_data).
bool NV2AIsFixedFunctionMode(PGRAPHState* pg);
bool NV2AIsFixedFunctionMode(NV2AState* d);
bool NV2AIsFixedFunctionMode(); // Uses g_NV2A global
