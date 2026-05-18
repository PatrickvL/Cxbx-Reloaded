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

#include <cstdint>

struct NV2AState; // Forward declaration
struct PGRAPHState; // Forward declaration

// ---- Texture Stage Scalars ----

uint32_t NV2AGetTextureControlRaw(NV2AState* d, int stage);
uint32_t NV2AGetTextureControlRaw(int stage);

// Read the raw TEXCTL1 register value for a stage (no decoding).
uint32_t NV2AGetTextureControl1Raw(NV2AState* d, int stage);
uint32_t NV2AGetTextureControl1Raw(int stage);

// Read the raw TEXOFFSET register value for a stage (no decoding).
uint32_t NV2AGetTextureOffsetRaw(NV2AState* d, int stage);
uint32_t NV2AGetTextureOffsetRaw(int stage);

// Read the raw TEXADDRESS register value for a stage (no decoding).
uint32_t NV2AGetTextureAddressModeRaw(NV2AState* d, int stage);
uint32_t NV2AGetTextureAddressModeRaw(int stage);

// Read the raw TEXFMT register value for a stage (no decoding).
uint32_t NV2AGetTextureFormatRaw(NV2AState* d, int stage);
uint32_t NV2AGetTextureFormatRaw(int stage);

// Read the raw TEXIMAGERECT register value for a stage (no decoding).
uint32_t NV2AGetTextureImageRectRaw(NV2AState* d, int stage);
uint32_t NV2AGetTextureImageRectRaw(int stage);

// Read the raw TEXFILTER register value for a stage (no decoding).
uint32_t NV2AGetTextureFilterRaw(NV2AState* d, int stage);
uint32_t NV2AGetTextureFilterRaw(int stage);

// Read the raw BORDERCOLOR register value for a stage (no decoding).
uint32_t NV2AGetBorderColorRaw(NV2AState* d, int stage);
uint32_t NV2AGetBorderColorRaw(int stage);

// Read the raw TEXPALETTE register value for a stage (no decoding).
uint32_t NV2AGetTexturePaletteRaw(NV2AState* d, int stage);
uint32_t NV2AGetTexturePaletteRaw(int stage);

// Check if a texture stage is enabled (TEXCTL0 bit 30).
bool NV2AIsTextureEnabled(NV2AState* d, int stage);
bool NV2AIsTextureEnabled(int stage);

// Read texture image pitch for a stage (from TEXCTL1 bits 16-31).
uint32_t NV2AGetTexturePitch(NV2AState* d, int stage);
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

// ---- Texture Stage State ----

struct NV2ATextureAddress {
	uint32_t physicalAddress; // DMA-resolved physical address in VRAM (0 if disabled/invalid)
	uint32_t rawOffset;       // Raw offset from PGRAPH register (before DMA resolution)
	uint32_t dmaBase;         // DMA context base address that was added
};

// Resolve the physical VRAM address of a texture's data for the given stage.
// Reads TEXOFFSET0 and resolves through DMA context A/B (selected by TEXFMT0 CONTEXT_DMA bit).
NV2ATextureAddress NV2AGetTextureAddress(NV2AState* d, int stage);
NV2ATextureAddress NV2AGetTextureAddress(int stage);

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
// NV097 SET_SURFACE_* methods write to these PGRAPH registers:
//   NV_PGRAPH_SURFACECLIPX   (0x19B4) - clip_x and clip_width
//   NV_PGRAPH_SURFACECLIPY   (0x19B8) - clip_y and clip_height
//   NV_PGRAPH_SURFACEFORMAT  (0x0714) - color/zeta format, type, AA, log width/height
//   NV_PGRAPH_DMA_PITCH      (0x0770) - color pitch (15:0), zeta pitch (31:16)
//   NV_PGRAPH_BOFFSET3       (0x082C) - color surface offset
//   NV_PGRAPH_BOFFSET4       (0x0830) - zeta surface offset
// All fields are decoded from register state; no side-struct is needed.

struct NV2ASurfaceState {
	uint32_t colorOffset;     // Raw color surface offset (relative to dma_color context)
	uint32_t zetaOffset;      // Raw zeta surface offset (relative to dma_zeta context)
	uint32_t colorPitch;      // Color surface pitch (bytes per row)
	uint32_t zetaPitch;       // Zeta surface pitch
	uint32_t clipX, clipY;    // Surface clip origin
	uint32_t clipWidth;       // Surface clip width
	uint32_t clipHeight;      // Surface clip height
	uint32_t antiAliasing;    // NV097_SET_SURFACE_FORMAT_ANTI_ALIASING_* value
	uint32_t colorFormat;     // Surface color format
	uint32_t zetaFormat;      // Surface zeta format
	uint32_t surfaceType;     // NV097_SET_SURFACE_FORMAT_TYPE_PITCH or _SWIZZLE
	uint32_t logWidth;        // log2(base width) for swizzle surfaces
	uint32_t logHeight;       // log2(base height) for swizzle surfaces
};

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
