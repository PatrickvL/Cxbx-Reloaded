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
// *  All rights reserved
// *
// ******************************************************************

// Backend_D3D11_PageTracker.cpp — Dirty page tracking for the contiguous memory mirror.
//
// CPU→GPU direction:
//   The 64 MiB contiguous region at 0x80000000 is allocated with VirtualAlloc +
//   MEM_WRITE_WATCH, providing zero-overhead write tracking via hardware PTE dirty
//   bits. At flush time, GetWriteWatch returns the list of modified pages, which
//   are then copied to the GPU mirror buffer. The write-watch is atomically reset.
//
// GPU→CPU direction:
//   When a render target is bound into contiguous memory, those pages are marked
//   GPU-dirty and set PAGE_NOACCESS. On CPU access (read or write), VEH restores
//   the page and triggers readback from the D3D11 RT.
//
// Tiled memory (0xF0000000):
//   Allocated as MEM_RESERVE + PAGE_NOACCESS (no file mapping alias). VEH commits
//   pages on demand, copying from the corresponding 0x80000000 address. At flush
//   time, any committed tiled pages are synced back to 0x80000000 and decommitted.
//
// The GPU mirror is a 64 MiB ByteAddressBuffer (DYNAMIC, SRV). The shader
// addresses it with: byteOffset = xboxPhysAddr & 0x07FFFFFF.

#include "Backend_D3D11_Internal.h"
#include "Backend_D3D11_PageTracker.h"
#include "common/AddressRanges.h"
#include "common/win32/WineEnv.h"
#include "devices/video/swizzle.h"

#include <cstring>
#include <immintrin.h>
#include <intrin.h>

// Wine may not reliably support MEM_WRITE_WATCH / GetWriteWatch.
// When running on Wine, always do a full memcpy instead.
static bool s_bWineFallback = false;

// AVX2 runtime detection — enables 256-bit bulk bitmap scan
static bool s_bHasAVX2 = [] {
	int cpuInfo[4] = {};
	__cpuid(cpuInfo, 0);
	if (cpuInfo[0] >= 7) {
		__cpuidex(cpuInfo, 7, 0);
		return (cpuInfo[1] & (1 << 5)) != 0; // EBX bit 5 = AVX2
	}
	return false;
}();

// ******************************************************************
// * Constants
// ******************************************************************
static constexpr uint32_t CONTIG_BASE     = CONTIGUOUS_MEMORY_BASE; // 0x80000000
static constexpr uint32_t CONTIG_SIZE     = XBOX_CONTIGUOUS_MEMORY_SIZE; // 64 MiB
static constexpr uint32_t TILED_BASE      = TILED_MEMORY_BASE; // 0xF0000000
static constexpr uint32_t TILED_SIZE      = TILED_MEMORY_SIZE; // 64 MiB
static constexpr uint32_t PAGE_SIZE_      = 4096;
static constexpr uint32_t PAGE_COUNT      = CONTIG_SIZE / PAGE_SIZE_; // 16384
static constexpr uint32_t BITMAP_DWORDS   = PAGE_COUNT / 32;         // 512

// PGRAPH register block appended after the 64 MiB RAM region.
// Offset 0x04000000 chosen to sit just past the 64 MiB RAM window.
static constexpr uint32_t GPU_PGRAPH_BASE = 0x04000000u;
static constexpr uint32_t GPU_PGRAPH_SIZE = 2048 * sizeof(uint32_t); // 8 KB

// PFB register block appended after PGRAPH.
// Offset 0x04002000 = GPU_PGRAPH_BASE + 8 KB (PGRAPH is 8 KB-aligned).
static constexpr uint32_t GPU_PFB_BASE = 0x04002000u;
static constexpr uint32_t GPU_PFB_SIZE = 1024 * sizeof(uint32_t); // 4 KB (NV_PFB_SIZE)

// PVIDEO register block appended after PFB.
// Offset 0x04003000 = GPU_PFB_BASE + 4 KB.
static constexpr uint32_t GPU_PVIDEO_BASE = 0x04003000u;
static constexpr uint32_t GPU_PVIDEO_SIZE = 1024 * sizeof(uint32_t); // 4 KB (NV_PVIDEO_SIZE)

// Total GPU buffer size: 64 MiB (RAM) + 8 KB (PGRAPH) + 4 KB (PFB) + 4 KB (PVIDEO)
static constexpr uint32_t GPU_BUFFER_SIZE = CONTIG_SIZE + GPU_PGRAPH_SIZE + GPU_PFB_SIZE + GPU_PVIDEO_SIZE;

// ******************************************************************
// * GPU-dirty bitmap (1 bit per 4 KB page) — set when RT writes here
// * Uses interlocked 32-bit ops for thread-safe VEH access.
// ******************************************************************
alignas(64) static volatile uint32_t s_GpuDirtyBitmap[BITMAP_DWORDS] = {};

// ******************************************************************
// * Texture-dirty bitmap — set when CPU-written pages are flushed to
// * the GPU mirror. Cleared per-texture after host upload completes.
// * This gates texture re-upload: if no texture-dirty pages overlap
// * a texture's address range, the host texture is still valid.
// * Applies to all texture types (swizzled, linear, compressed).
// * Uses interlocked 32-bit ops for thread-safe access.
// ******************************************************************
alignas(64) static volatile uint32_t s_TextureDirtyBitmap[BITMAP_DWORDS] = {};

// ******************************************************************
// * Tiled committed bitmap — tracks which 0xF0 pages are committed
// * (single-threaded access only, no atomics needed)
// ******************************************************************
static uint32_t s_TiledCommittedBitmap[BITMAP_DWORDS] = {};

// Quick flag: true when any tiled page is committed (avoids scanning bitmap)
static bool s_bHasTiledPages = false;

// Frame boundary flag: true for the first flush after Present.
// Only the first flush of a frame may use MAP_WRITE_DISCARD (which orphans
// the buffer). Subsequent mid-frame flushes use MAP_WRITE_NO_OVERWRITE so
// that earlier draw calls in the same frame still see their data.
static bool s_bFirstFlushOfFrame = true;

// ******************************************************************
// * Registered RT metadata for GPU→CPU readback
// ******************************************************************
struct RegisteredRT {
	uint32_t offset;        // VRAM byte offset from CONTIGUOUS_MEMORY_BASE
	uint32_t pitch;         // Xbox row pitch in bytes
	uint32_t width;         // Xbox width in pixels
	uint32_t height;        // Xbox height in pixels
	uint32_t bpp;           // Bytes per pixel (2 or 4)
	uint32_t surfaceType;   // NV097_SET_SURFACE_FORMAT_TYPE_PITCH or _SWIZZLE
	ID3D11Texture2D* pTexture;    // Host RT (NOT AddRef'd — owned by g_PgraphRTCache)
	ID3D11Texture2D* pStagingTex; // Cached staging texture (AddRef'd, we own it)
	bool needsReadback;     // True when RT has been rendered to since last readback
};

static constexpr uint32_t MAX_REGISTERED_RTS = 16;
static RegisteredRT s_RegisteredRTs[MAX_REGISTERED_RTS] = {};
static uint32_t s_NumRegisteredRTs = 0;

// Critical section for serializing D3D11 device context access between the
// puller thread (normal rendering) and the VEH readback path (CPU thread).
// The readback uses TryEnterCriticalSection — if the puller holds it, readback
// is skipped (graceful degradation to stale data, same as pre-fix behavior).
static CRITICAL_SECTION s_D3D11ContextLock;
static bool s_D3D11ContextLockInitialized = false;

// Thread ID that currently owns s_D3D11ContextLock.  Used by the VEH handler
// to detect same-thread re-entrancy: Windows critical sections are recursive,
// so TryEnterCriticalSection succeeds even when the *same* thread already holds
// the lock.  Re-entering D3D11 (CopyResource/Map while inside UpdateSubresource)
// corrupts driver state and crashes the Intel UMD.
static volatile DWORD s_D3D11ContextOwnerThread = 0;

// ******************************************************************
// * GPU mirror buffer (64 MiB ByteAddressBuffer — DEFAULT + SRV + UAV)
// ******************************************************************
static ID3D11Buffer*              s_pMirrorBuf = nullptr;
static ID3D11ShaderResourceView*  s_pMirrorSRV = nullptr;
static ID3D11ShaderResourceView*  s_pMirrorSRV_SNORM16x2 = nullptr; // R16G16_SNORM typed view
static ID3D11ShaderResourceView*  s_pMirrorSRV_UNORM8x4 = nullptr;  // R8G8B8A8_UNORM typed view
static ID3D11UnorderedAccessView* s_pMirrorUAV = nullptr;            // RWByteAddressBuffer for CS

// ******************************************************************
// * VEH handle
// ******************************************************************
static void* s_hVEH = nullptr;

// ******************************************************************
// * Static buffer for GetWriteWatch results (16384 pointers)
// ******************************************************************
static PVOID s_WriteWatchPages[PAGE_COUNT];

// ******************************************************************
// * Bitmap helpers — atomic 32-bit (for GpuDirty/TextureDirty)
// ******************************************************************
static inline void SetBitAtomic(volatile uint32_t* bitmap, uint32_t index)
{
	_InterlockedOr((volatile long*)&bitmap[index >> 5], 1L << (index & 31));
}

static inline void ClearBitAtomic(volatile uint32_t* bitmap, uint32_t index)
{
	_InterlockedAnd((volatile long*)&bitmap[index >> 5], ~(1L << (index & 31)));
}

static inline bool TestBitAtomic(const volatile uint32_t* bitmap, uint32_t index)
{
	return (bitmap[index >> 5] & (1u << (index & 31))) != 0;
}

// ******************************************************************
// * Bitmap helpers — non-atomic (for TiledCommitted, single-threaded)
// ******************************************************************
static inline void SetBit(uint32_t* bitmap, uint32_t index)
{
	bitmap[index >> 5] |= (1u << (index & 31));
}

static inline void ClearBit(uint32_t* bitmap, uint32_t index)
{
	bitmap[index >> 5] &= ~(1u << (index & 31));
}

static inline bool TestBit(const uint32_t* bitmap, uint32_t index)
{
	return (bitmap[index >> 5] & (1u << (index & 31))) != 0;
}

// ******************************************************************
// * Handle a fault (GPU-dirty pages or tiled redirect)
// * Called directly from the unified VEH in lleException.
// ******************************************************************

// Forward declaration (defined after CxbxPageTrackerRegisterRT)
static ID3D11Texture2D* CreateStagingForRT(ID3D11Texture2D* pTexture);

// Try to handle an access violation in the contiguous or tiled region.
// Returns true if the fault was handled (page committed/restored).
bool CxbxPageTrackerHandleFault(void* faultAddress, bool isWrite)
{
	uintptr_t addr = (uintptr_t)faultAddress;

	// --- Tiled memory redirect (0xF0000000 - 0xF3FFFFFF) ---
	if (addr >= TILED_BASE && addr < (TILED_BASE + TILED_SIZE)) {
		uint32_t offset = (uint32_t)(addr - TILED_BASE);
		uint32_t pageIdx = offset / PAGE_SIZE_;
		uint32_t pageOffset = pageIdx * PAGE_SIZE_;

		if (!TestBit(s_TiledCommittedBitmap, pageIdx)) {
			// Commit the page on demand
			LPVOID result = VirtualAlloc(
				(LPVOID)(TILED_BASE + pageOffset), PAGE_SIZE_,
				MEM_COMMIT, PAGE_READWRITE);
			if (result == nullptr) return false;

			// Copy current data from contiguous memory
			memcpy((void*)(TILED_BASE + pageOffset),
			       (void*)(CONTIG_BASE + pageOffset), PAGE_SIZE_);

			SetBit(s_TiledCommittedBitmap, pageIdx);
			s_bHasTiledPages = true;
		}
		return true;
	}

	// --- GPU-dirty page handling (0x80000000 - 0x83FFFFFF) ---
	if (addr >= CONTIG_BASE && addr < (CONTIG_BASE + CONTIG_SIZE)) {
		uint32_t offset = (uint32_t)(addr - CONTIG_BASE);
		uint32_t pageIdx = offset / PAGE_SIZE_;

		if (TestBitAtomic(s_GpuDirtyBitmap, pageIdx)) {
			// Perform readback from D3D11 render target into Xbox memory.
			// On CPU reads, copy the entire RT back to Xbox RAM so all pages
			// in the RT are restored at once (amortizes the GPU stall).
			// On CPU writes, skip readback — the CPU is overwriting the data.
			if (!isWrite && g_pD3DDeviceContext != nullptr &&
				s_D3D11ContextLockInitialized &&
				s_D3D11ContextOwnerThread != GetCurrentThreadId() &&
				TryEnterCriticalSection(&s_D3D11ContextLock)) {
				s_D3D11ContextOwnerThread = GetCurrentThreadId();
				// Find which registered RT covers this page
				const RegisteredRT* pRT = nullptr;
				for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
					uint32_t rtEnd = s_RegisteredRTs[i].offset +
						s_RegisteredRTs[i].pitch * s_RegisteredRTs[i].height;
					if (offset >= s_RegisteredRTs[i].offset && offset < rtEnd) {
						pRT = &s_RegisteredRTs[i];
						break;
					}
				}

				if (pRT && pRT->pTexture) {
					// FIRST: Clear GPU-dirty and restore access for the ENTIRE RT
					// before any memcpy. Otherwise the row-by-row copy would fault
					// on adjacent pages that are still PAGE_NOACCESS.
					uint32_t rtSize = pRT->pitch * pRT->height;
					uint32_t rtFirstPage = pRT->offset / PAGE_SIZE_;
					uint32_t rtLastPage = (pRT->offset + rtSize - 1) / PAGE_SIZE_;
					for (uint32_t p = rtFirstPage; p <= rtLastPage; p++) {
						if (TestBitAtomic(s_GpuDirtyBitmap, p)) {
							ClearBitAtomic(s_GpuDirtyBitmap, p);
							DWORD oldProtect;
							VirtualProtect((void*)(CONTIG_BASE + p * PAGE_SIZE_),
								PAGE_SIZE_, PAGE_READWRITE, &oldProtect);
						}
					}

					// Now perform the actual GPU→CPU readback via cached staging texture
					ID3D11Texture2D* pStaging = pRT->pStagingTex;
					if (!pStaging) {
						// Fallback: create on the fly if staging was not cached
						pStaging = CreateStagingForRT(pRT->pTexture);
					}
					if (pStaging) {
						D3D11_TEXTURE2D_DESC desc = {};
						pRT->pTexture->GetDesc(&desc);

						g_pD3DDeviceContext->CopyResource(pStaging, pRT->pTexture);

						D3D11_MAPPED_SUBRESOURCE mapped = {};
						HRESULT hr = g_pD3DDeviceContext->Map(pStaging, 0, D3D11_MAP_READ, 0, &mapped);
						if (SUCCEEDED(hr)) {
							// Copy from staging to Xbox RAM.
							// Host RT may have different pitch than Xbox RT.
							uint8_t* pDst = (uint8_t*)(CONTIG_BASE + pRT->offset);
							uint8_t* pSrc = (uint8_t*)mapped.pData;

							// If host is upscaled, only copy the top-left 1x region
							if (desc.Width > pRT->width || desc.Height > pRT->height) {
								// Can't directly copy upscaled data — skip readback
								// (would need a resolve/downscale pass)
							} else if (pRT->surfaceType == 0x2 /*SWIZZLE*/) {
								// NV2A swizzled surface: write pixels in Morton order.
								swizzle_rect(pSrc, pRT->width, pRT->height,
									pDst, mapped.RowPitch, pRT->bpp);
							} else {
								uint32_t xboxRowBytes = pRT->width * pRT->bpp;
								for (uint32_t row = 0; row < pRT->height; row++) {
									memcpy(pDst, pSrc, xboxRowBytes);
									pDst += pRT->pitch;
									pSrc += mapped.RowPitch;
								}
							}

							g_pD3DDeviceContext->Unmap(pStaging, 0);
						}
						// Only release if it was a fallback allocation
						if (pStaging != pRT->pStagingTex) {
							pStaging->Release();
						}
					}

					s_D3D11ContextOwnerThread = 0;
					LeaveCriticalSection(&s_D3D11ContextLock);
					return true;
				}

				s_D3D11ContextOwnerThread = 0;
				LeaveCriticalSection(&s_D3D11ContextLock);
			}

			// Fallback: no matching RT found, write access, or lock contended —
			// just restore the page (graceful degradation to stale data).
			ClearBitAtomic(s_GpuDirtyBitmap, pageIdx);

			DWORD oldProtect;
			VirtualProtect((void*)(CONTIG_BASE + pageIdx * PAGE_SIZE_),
				PAGE_SIZE_, PAGE_READWRITE, &oldProtect);
			return true;
		}

		// Safety net: the page is in the contiguous region but the GPU-dirty
		// bit is not set.  This can happen due to a race between MarkGPUDirty
		// (on the puller thread) and the VEH clearing the bit (on a game
		// thread).  Unconditionally restore PAGE_READWRITE — for committed
		// contiguous memory, PAGE_NOACCESS is only ever set by MarkGPUDirty,
		// so restoring access is always the correct recovery action.
		{
			DWORD oldProtect;
			VirtualProtect((void*)(CONTIG_BASE + pageIdx * PAGE_SIZE_),
				PAGE_SIZE_, PAGE_READWRITE, &oldProtect);
		}
		return true;
	}

	return false;
}

// ******************************************************************
// * VEH handler — GPU-dirty faults + tiled memory redirect
// ******************************************************************
// PageTrackerVEH removed — fault handling is now inlined in lleException
// to avoid paying double VEH dispatch overhead on every MMIO access.

// ******************************************************************
// * Sync committed tiled pages back to contiguous memory, then decommit
// ******************************************************************
static void SyncTiledPagesBack()
{
	if (!s_bHasTiledPages)
		return;

	for (uint32_t dw = 0; dw < BITMAP_DWORDS; dw++) {
		uint32_t bits = s_TiledCommittedBitmap[dw];
		if (bits == 0) continue;

		while (bits) {
			unsigned long pos;
			_BitScanForward(&pos, (unsigned long)bits);
			bits &= bits - 1;

			uint32_t pageIdx = dw * 32 + pos;
			uint32_t offset = pageIdx * PAGE_SIZE_;

			// Copy tiled page back to contiguous memory
			// (this write is automatically tracked by MEM_WRITE_WATCH)
			memcpy((void*)(CONTIG_BASE + offset),
			       (void*)(TILED_BASE + offset), PAGE_SIZE_);

			// Decommit — returns to MEM_RESERVE + PAGE_NOACCESS, faults again on next access
			VirtualFree((void*)(TILED_BASE + offset), PAGE_SIZE_, MEM_DECOMMIT);
		}

		s_TiledCommittedBitmap[dw] = 0;
	}

	s_bHasTiledPages = false;
}

// ******************************************************************
// * Public: Initialize page tracking
// ******************************************************************
void CxbxPageTrackerInit()
{
	memset((void*)s_GpuDirtyBitmap, 0, sizeof(s_GpuDirtyBitmap));
	memset(s_TiledCommittedBitmap, 0, sizeof(s_TiledCommittedBitmap));
	// All pages start texture-dirty so the first deswizzle for each texture is triggered
	memset((void*)s_TextureDirtyBitmap, 0xFF, sizeof(s_TextureDirtyBitmap));

	// Initialize D3D11 context lock for thread-safe readback from VEH
	if (!s_D3D11ContextLockInitialized) {
		InitializeCriticalSection(&s_D3D11ContextLock);
		s_D3D11ContextLockInitialized = true;
	}

	// Create GPU mirror buffer: 64 MiB RAM + 8 KB PGRAPH (DEFAULT, SRV + UAV)
	D3D11_BUFFER_DESC desc = {};
	desc.ByteWidth = GPU_BUFFER_SIZE;
	desc.Usage = D3D11_USAGE_DEFAULT;
	desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS;
	desc.CPUAccessFlags = 0;
	desc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS;

	HRESULT hr = g_pD3DDevice->CreateBuffer(&desc, nullptr, &s_pMirrorBuf);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to create mirror buffer (hr=0x%08X)", hr);
		return;
	}

	// Create raw buffer SRV (covers full buffer: RAM + PGRAPH)
	D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
	srvDesc.Format = DXGI_FORMAT_R32_TYPELESS;
	srvDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFEREX;
	srvDesc.BufferEx.FirstElement = 0;
	srvDesc.BufferEx.NumElements = GPU_BUFFER_SIZE / 4;
	srvDesc.BufferEx.Flags = D3D11_BUFFEREX_SRV_FLAG_RAW;

	hr = g_pD3DDevice->CreateShaderResourceView(s_pMirrorBuf, &srvDesc, &s_pMirrorSRV);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to create mirror SRV (hr=0x%08X)", hr);
		s_pMirrorBuf->Release();
		s_pMirrorBuf = nullptr;
		return;
	}

	// Create typed SRV views for hardware format decode (vertex fetch)
	{
		D3D11_SHADER_RESOURCE_VIEW_DESC typedDesc = {};
		typedDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
		typedDesc.Buffer.FirstElement = 0;
		typedDesc.Buffer.NumElements = CONTIG_SIZE / 4; // 4 bytes per element

		// R16G16_SNORM: each element = 4 bytes → 2 signed normalized shorts
		typedDesc.Format = DXGI_FORMAT_R16G16_SNORM;
		hr = g_pD3DDevice->CreateShaderResourceView(s_pMirrorBuf, &typedDesc, &s_pMirrorSRV_SNORM16x2);
		if (FAILED(hr))
			EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to create SNORM16x2 SRV (hr=0x%08X)", hr);

		// R8G8B8A8_UNORM: each element = 4 bytes → 4 unsigned normalized bytes
		typedDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
		hr = g_pD3DDevice->CreateShaderResourceView(s_pMirrorBuf, &typedDesc, &s_pMirrorSRV_UNORM8x4);
		if (FAILED(hr))
			EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to create UNORM8x4 SRV (hr=0x%08X)", hr);
	}

	// Create RWByteAddressBuffer UAV (covers full buffer: RAM + PGRAPH)
	{
		D3D11_UNORDERED_ACCESS_VIEW_DESC uavDesc = {};
		uavDesc.Format = DXGI_FORMAT_R32_TYPELESS;
		uavDesc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
		uavDesc.Buffer.FirstElement = 0;
		uavDesc.Buffer.NumElements = GPU_BUFFER_SIZE / 4;
		uavDesc.Buffer.Flags = D3D11_BUFFER_UAV_FLAG_RAW;
		hr = g_pD3DDevice->CreateUnorderedAccessView(s_pMirrorBuf, &uavDesc, &s_pMirrorUAV);
		if (FAILED(hr))
			EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to create mirror UAV (hr=0x%08X)", hr);
	}

	// Reset write-watch BEFORE the initial upload. This ensures that any writes
	// occurring concurrently (from the Xbox title thread) during the memcpy will
	// have their dirty bits preserved and picked up by the first flush. If we
	// reset AFTER the memcpy, writes that happen between the memcpy passing a page
	// and the reset would have their dirty bits cleared, leaving stale data in the
	// mirror with no way to detect it later.
	ResetWriteWatch((PVOID)CONTIG_BASE, CONTIG_SIZE);

	// Initial full upload of contiguous memory to GPU mirror.
	// Must use an explicit box — the buffer is larger than CONTIG_SIZE
	// (appended PGRAPH/PFB/PVIDEO), so a nullptr box would tell the
	// driver to read ByteWidth bytes from pSrcData, overrunning the
	// 64 MiB contiguous region.
	{
		D3D11_BOX box = { 0, 0, 0, CONTIG_SIZE, 1, 1 };
		g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box,
			(const void*)CONTIG_BASE, CONTIG_SIZE, 0);
	}

	// Detect Wine — GetWriteWatch may not reliably track dirty pages.
	// When running on Wine, always do a full upload on every flush.
	s_bWineFallback = isWineEnv();
	if (s_bWineFallback) {
		EmuLog(LOG_LEVEL::INFO, "PageTracker: Wine detected — using full-upload fallback (GetWriteWatch unreliable)");
	}

	// Verify that GetWriteWatch actually works on this allocation.
	// The contiguous region must have been allocated with MEM_WRITE_WATCH
	// (done in ReserveAddressRanges.cpp). If it wasn't, fall back to full upload.
	if (!s_bWineFallback) {
		PVOID testAddr;
		ULONG_PTR testCount = 1;
		ULONG testGranularity;
		UINT testResult = GetWriteWatch(0, (PVOID)CONTIG_BASE, CONTIG_SIZE,
			&testAddr, &testCount, &testGranularity);
		if (testResult != 0) {
			s_bWineFallback = true; // reuse Wine fallback path for full-upload
			EmuLog(LOG_LEVEL::WARNING, "PageTracker: GetWriteWatch failed (error %u) — using full-upload fallback", GetLastError());
		} else {
			EmuLog(LOG_LEVEL::INFO, "PageTracker: GetWriteWatch verified working (granularity=%u)", testGranularity);
		}
	}

	// VEH for page faults is now handled by the unified lleException handler.
	// No separate VEH registration needed here.

	EmuLog(LOG_LEVEL::INFO, "PageTracker: Initialized (%s, %u pages tracked)",
		s_bWineFallback ? "Wine full-upload" : "MEM_WRITE_WATCH", PAGE_COUNT);
}

// ******************************************************************
// * Public: Shutdown
// ******************************************************************
void CxbxPageTrackerShutdown()
{
	if (s_hVEH) {
		RemoveVectoredExceptionHandler(s_hVEH);
		s_hVEH = nullptr;
	}

	// Restore GPU-dirty pages to normal access
	for (uint32_t dw = 0; dw < BITMAP_DWORDS; dw++) {
		uint32_t bits = s_GpuDirtyBitmap[dw];
		if (bits == 0) continue;
		while (bits) {
			unsigned long pos;
			_BitScanForward(&pos, (unsigned long)bits);
			bits &= bits - 1;
			uint32_t offset = (dw * 32 + pos) * PAGE_SIZE_;
			DWORD oldProtect;
			VirtualProtect((void*)(CONTIG_BASE + offset), PAGE_SIZE_, PAGE_READWRITE, &oldProtect);
		}
	}

	// Decommit any committed tiled pages
	for (uint32_t dw = 0; dw < BITMAP_DWORDS; dw++) {
		uint32_t bits = s_TiledCommittedBitmap[dw];
		if (bits == 0) continue;
		while (bits) {
			unsigned long pos;
			_BitScanForward(&pos, (unsigned long)bits);
			bits &= bits - 1;
			uint32_t offset = (dw * 32 + pos) * PAGE_SIZE_;
			VirtualFree((void*)(TILED_BASE + offset), PAGE_SIZE_, MEM_DECOMMIT);
		}
	}

	// Release cached staging textures for all registered RTs
	for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
		if (s_RegisteredRTs[i].pStagingTex) {
			s_RegisteredRTs[i].pStagingTex->Release();
			s_RegisteredRTs[i].pStagingTex = nullptr;
		}
	}

	if (s_pMirrorUAV) { s_pMirrorUAV->Release(); s_pMirrorUAV = nullptr; }
	if (s_pMirrorSRV_UNORM8x4) { s_pMirrorSRV_UNORM8x4->Release(); s_pMirrorSRV_UNORM8x4 = nullptr; }
	if (s_pMirrorSRV_SNORM16x2) { s_pMirrorSRV_SNORM16x2->Release(); s_pMirrorSRV_SNORM16x2 = nullptr; }
	if (s_pMirrorSRV) { s_pMirrorSRV->Release(); s_pMirrorSRV = nullptr; }
	if (s_pMirrorBuf) { s_pMirrorBuf->Release(); s_pMirrorBuf = nullptr; }

	memset((void*)s_GpuDirtyBitmap, 0, sizeof(s_GpuDirtyBitmap));
	memset(s_TiledCommittedBitmap, 0, sizeof(s_TiledCommittedBitmap));
	memset((void*)s_TextureDirtyBitmap, 0, sizeof(s_TextureDirtyBitmap));
}

// ******************************************************************
// * Public: Flush CPU-dirty pages to GPU mirror
// ******************************************************************
uint32_t CxbxPageTrackerFlushToGPU()
{
	if (!s_pMirrorBuf)
		return 0;

	// Skip redundant flushes within the same frame. The Xbox CPU builds the
	// entire pushbuffer before kicking it, so between individual draws within
	// a single pushbuffer submission, no new CPU writes to VRAM are expected.
	// Flushing once at the start of each frame (first draw after Present) is
	// sufficient. This avoids 20k+ GetWriteWatch syscalls per frame.
	if (!s_bFirstFlushOfFrame)
		return 0;

	// Sync any committed tiled pages back to contiguous memory.
	// The memcpy writes to 0x80 are automatically tracked by write-watch.
	SyncTiledPagesBack();

	// Wine / broken-write-watch fallback: always do a full upload.
	if (s_bWineFallback) {
		memset((void*)s_TextureDirtyBitmap, 0xFF, sizeof(s_TextureDirtyBitmap));
		D3D11_BOX box = { 0, 0, 0, CONTIG_SIZE, 1, 1 };
		g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box,
			(const void*)CONTIG_BASE, CONTIG_SIZE, 0);
		s_bFirstFlushOfFrame = false;
		return PAGE_COUNT;
	}

	// Get dirty pages from write-watch (atomically reads and resets)
	ULONG_PTR count = PAGE_COUNT;
	ULONG granularity;
	UINT result = GetWriteWatch(
		WRITE_WATCH_FLAG_RESET,
		(PVOID)CONTIG_BASE, CONTIG_SIZE,
		s_WriteWatchPages, &count, &granularity);

	if (result != 0 || count == 0)
		return 0;

	// Mark all flushed pages as texture-dirty (for deswizzle gating).
	// This must happen regardless of the upload strategy below, so that
	// HostResourceRequiresUpdate can detect which textures need re-deswizzle.
	if (count > PAGE_COUNT / 4) {
		// Bulk dirty — mark all pages as texture-dirty
		memset((void*)s_TextureDirtyBitmap, 0xFF, sizeof(s_TextureDirtyBitmap));
	} else {
		for (ULONG_PTR i = 0; i < count; i++) {
			uint32_t pageIdx = (uint32_t)((uintptr_t)s_WriteWatchPages[i] - CONTIG_BASE) / PAGE_SIZE_;
			SetBitAtomic(s_TextureDirtyBitmap, pageIdx);
		}
	}

	// Many pages dirty: full-buffer upload is cheaper than many small UpdateSubresource calls.
	if (count > PAGE_COUNT / 4) {
		D3D11_BOX box = { 0, 0, 0, CONTIG_SIZE, 1, 1 };
		g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box,
			(const void*)CONTIG_BASE, CONTIG_SIZE, 0);
	} else {
		// Incremental update: one UpdateSubresource per coalesced page run.
		// GetWriteWatch guarantees ascending address order per MSDN — no sort needed.
		ULONG_PTR runStart = 0;
		while (runStart < count) {
			uintptr_t startAddr = (uintptr_t)s_WriteWatchPages[runStart];
			ULONG_PTR runEnd = runStart + 1;

			// Extend run while next page is contiguous
			while (runEnd < count &&
				(uintptr_t)s_WriteWatchPages[runEnd] == startAddr + (runEnd - runStart) * PAGE_SIZE_) {
				runEnd++;
			}

			uint32_t offset = (uint32_t)(startAddr - CONTIG_BASE);
			uint32_t runBytes = (uint32_t)(runEnd - runStart) * PAGE_SIZE_;
			D3D11_BOX box = { offset, 0, 0, offset + runBytes, 1, 1 };
			g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box,
				(const void*)startAddr, runBytes, 0);

			runStart = runEnd;
		}
	}
	s_bFirstFlushOfFrame = false;

	return (uint32_t)count;
}

// ******************************************************************
// * Public: Notify frame boundary (called from CxbxPresent)
// ******************************************************************
void CxbxPageTrackerOnPresent()
{
	s_bFirstFlushOfFrame = true;
}

// ******************************************************************
// * Public: Lock/unlock D3D11 context for puller thread
// ******************************************************************
void CxbxPageTrackerLockD3D11Context()
{
	if (s_D3D11ContextLockInitialized) {
		EnterCriticalSection(&s_D3D11ContextLock);
		s_D3D11ContextOwnerThread = GetCurrentThreadId();
	}
}

void CxbxPageTrackerUnlockD3D11Context()
{
	if (s_D3D11ContextLockInitialized) {
		s_D3D11ContextOwnerThread = 0;
		LeaveCriticalSection(&s_D3D11ContextLock);
	}
}

// ******************************************************************
// * Public: Check if any pages are CPU-dirty
// ******************************************************************
bool CxbxPageTrackerHasDirtyPages()
{
	// Quick check for committed tiled pages
	if (s_bHasTiledPages) return true;

	// Peek at write-watch (do not reset)
	// On Wine, always report dirty since we can't reliably check
	if (s_bWineFallback) return true;
	PVOID addr;
	ULONG_PTR count = 1;
	ULONG granularity;
	UINT result = GetWriteWatch(0, (PVOID)CONTIG_BASE, CONTIG_SIZE,
		&addr, &count, &granularity);
	return (result == 0 && count > 0);
}

// ******************************************************************
// * Public: Mark pages as GPU-dirty (render target bound here)
// ******************************************************************
void CxbxPageTrackerMarkGPUDirty(uint32_t startOffset, uint32_t size)
{
	if (startOffset >= CONTIG_SIZE) return;
	if (startOffset + size > CONTIG_SIZE) size = CONTIG_SIZE - startOffset;

	uint32_t firstPage = startOffset / PAGE_SIZE_;
	uint32_t lastPage = (startOffset + size - 1) / PAGE_SIZE_;

	for (uint32_t p = firstPage; p <= lastPage; p++) {
		SetBitAtomic(s_GpuDirtyBitmap, p);
		// Set PAGE_NOACCESS so CPU reads/writes fault for readback
		DWORD oldProtect;
		VirtualProtect((void*)(CONTIG_BASE + p * PAGE_SIZE_),
			PAGE_SIZE_, PAGE_NOACCESS, &oldProtect);
	}

	// Mark any registered RTs overlapping this range as needing readback.
	// This ensures that after MarkGPUDirty (RT rendered to), the next
	// FlushGPUDirtyToMirror will re-read the D3D11 texture contents.
	for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
		RegisteredRT& rt = s_RegisteredRTs[i];
		uint32_t rtEnd = rt.offset + rt.pitch * rt.height;
		if (startOffset < rtEnd && rt.offset < startOffset + size) {
			rt.needsReadback = true;
		}
	}
}

// ******************************************************************
// * Public: Register an RT for readback (called alongside MarkGPUDirty)
// ******************************************************************
// Helper: create a STAGING texture matching the given RT texture dimensions.
static ID3D11Texture2D* CreateStagingForRT(ID3D11Texture2D* pTexture)
{
	if (!pTexture) return nullptr;
	D3D11_TEXTURE2D_DESC desc = {};
	pTexture->GetDesc(&desc);
	desc.Usage = D3D11_USAGE_STAGING;
	desc.BindFlags = 0;
	desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
	desc.MiscFlags = 0;
	ID3D11Texture2D* pStaging = nullptr;
	HRESULT hr = g_pD3DDevice->CreateTexture2D(&desc, nullptr, &pStaging);
	return SUCCEEDED(hr) ? pStaging : nullptr;
}

void CxbxPageTrackerRegisterRT(uint32_t startOffset, uint32_t pitch,
	uint32_t width, uint32_t height, uint32_t bpp,
	uint32_t surfaceType, ID3D11Texture2D* pTexture)
{
	// Check if already registered at this offset — update in place
	for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
		if (s_RegisteredRTs[i].offset == startOffset) {
			// If RT dimensions/texture changed, recreate staging
			if (s_RegisteredRTs[i].pTexture != pTexture ||
				s_RegisteredRTs[i].width != width ||
				s_RegisteredRTs[i].height != height ||
				s_RegisteredRTs[i].bpp != bpp) {
				if (s_RegisteredRTs[i].pStagingTex) {
					s_RegisteredRTs[i].pStagingTex->Release();
				}
				s_RegisteredRTs[i] = { startOffset, pitch, width, height, bpp,
					surfaceType, pTexture, CreateStagingForRT(pTexture), true };
			} else {
				s_RegisteredRTs[i].pitch = pitch;
				s_RegisteredRTs[i].surfaceType = surfaceType;
				s_RegisteredRTs[i].needsReadback = true;
			}
			return;
		}
	}

	// Add new entry (evict oldest if full)
	if (s_NumRegisteredRTs >= MAX_REGISTERED_RTS) {
		// Release evicted entry's staging texture
		if (s_RegisteredRTs[0].pStagingTex) {
			s_RegisteredRTs[0].pStagingTex->Release();
		}
		// Shift down (FIFO eviction — oldest RT is least likely to be read back)
		memmove(&s_RegisteredRTs[0], &s_RegisteredRTs[1],
			(MAX_REGISTERED_RTS - 1) * sizeof(RegisteredRT));
		s_NumRegisteredRTs = MAX_REGISTERED_RTS - 1;
	}
	s_RegisteredRTs[s_NumRegisteredRTs++] = { startOffset, pitch, width, height, bpp,
		surfaceType, pTexture, CreateStagingForRT(pTexture), true };
}

// ******************************************************************
// * Public: Check if a page is GPU-dirty
// ******************************************************************
bool CxbxPageTrackerIsGPUDirty(uint32_t pageIndex)
{
	if (pageIndex >= PAGE_COUNT) return false;
	return TestBitAtomic(s_GpuDirtyBitmap, pageIndex);
}

// ******************************************************************
// * Public: Clear GPU-dirty flags after readback
// ******************************************************************
void CxbxPageTrackerClearGPUDirty(uint32_t startOffset, uint32_t size)
{
	if (startOffset >= CONTIG_SIZE) return;
	if (startOffset + size > CONTIG_SIZE) size = CONTIG_SIZE - startOffset;

	uint32_t firstPage = startOffset / PAGE_SIZE_;
	uint32_t lastPage = (startOffset + size - 1) / PAGE_SIZE_;

	for (uint32_t p = firstPage; p <= lastPage; p++) {
		ClearBitAtomic(s_GpuDirtyBitmap, p);
	}
}

// ******************************************************************
// * Public: Flush GPU-dirty pages overlapping a VB range to the mirror
// *
// * When a vertex buffer aliases render target memory, the D3D11 RT
// * content must be read back to Xbox RAM and then uploaded to the
// * GPU mirror buffer before the draw can read correct vertex data.
// ******************************************************************
bool CxbxPageTrackerFlushGPUDirtyToMirror(uint32_t startOffset, uint32_t size)
{
	if (!s_pMirrorBuf || !g_pD3DDeviceContext || size == 0 || startOffset >= CONTIG_SIZE)
		return false;
	if (startOffset + size > CONTIG_SIZE)
		size = CONTIG_SIZE - startOffset;

	uint32_t firstPage = startOffset / PAGE_SIZE_;
	uint32_t lastPage = (startOffset + size - 1) / PAGE_SIZE_;

	// Check if any GPU-dirty pages exist in this range.
	// If so, perform RT readback from D3D11 into Xbox RAM first.
	bool anyDirty = false;
	for (uint32_t p = firstPage; p <= lastPage; p++) {
		if (TestBitAtomic(s_GpuDirtyBitmap, p)) {
			anyDirty = true;
			break;
		}
	}

	// Always readback registered RTs that overlap this VB range, regardless of
	// dirty bits. The VEH fault handler may have cleared dirty bits after a
	// premature readback (before the D3D11 texture was fully rendered). By always
	// re-reading the D3D11 texture here, we guarantee the vertex shader sees the
	// latest rendered content.
	//
	// Swizzle policy: If ANY registered RT has surfaceType==SWIZZLE (0x2), the
	// game uses Morton-order index/vertex addressing (e.g. DisplacementMap XDK
	// sample). In that case ALL overlapping RTs are written in swizzled/Morton
	// layout regardless of their individual surfaceType — the IB expects all
	// attributes to share the same Morton addressing. If no registered RT is
	// SWIZZLE, each RT is written per its own surfaceType (PITCH = linear).
	//
	// Note: We check ALL registered RTs, not just those overlapping the current
	// VB range. The SWIZZLE RT might be at a different offset but still signals
	// that the entire draw uses Morton addressing for all attribute streams.
	{
		// First pass: determine if any registered RT is SWIZZLE
		bool bMortonDraw = false;
		for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
			if (s_RegisteredRTs[i].surfaceType == 0x2) {
				bMortonDraw = true;
				break;
			}
		}

		// Second pass: readback each overlapping RT
		for (uint32_t i = 0; i < s_NumRegisteredRTs; i++) {
			RegisteredRT& rt = s_RegisteredRTs[i];
			uint32_t rtSize = rt.pitch * rt.height;
			uint32_t rtEnd = rt.offset + rtSize;

			// Check overlap with VB range
			if (startOffset >= rtEnd || rt.offset >= startOffset + size)
				continue;

			if (!rt.pTexture || !g_pD3DDeviceContext)
				continue;

			// Ensure pages are accessible for the readback write
			uint32_t rtFirstPage = rt.offset / PAGE_SIZE_;
			uint32_t rtLastPage = (rtEnd - 1) / PAGE_SIZE_;
			for (uint32_t p = rtFirstPage; p <= rtLastPage; p++) {
				if (TestBitAtomic(s_GpuDirtyBitmap, p)) {
					ClearBitAtomic(s_GpuDirtyBitmap, p);
					DWORD oldProtect;
					VirtualProtect((void*)(CONTIG_BASE + p * PAGE_SIZE_),
						PAGE_SIZE_, PAGE_READWRITE, &oldProtect);
				}
			}

			// Readback RT from D3D11 → cached staging → Xbox RAM
			ID3D11Texture2D* pStaging = rt.pStagingTex;
			if (!pStaging) {
				// Fallback: create on the fly if staging was not cached
				pStaging = CreateStagingForRT(rt.pTexture);
			}
			if (pStaging) {
				D3D11_TEXTURE2D_DESC desc = {};
				rt.pTexture->GetDesc(&desc);

				g_pD3DDeviceContext->CopyResource(pStaging, rt.pTexture);

				D3D11_MAPPED_SUBRESOURCE mapped = {};
				HRESULT hr = g_pD3DDeviceContext->Map(pStaging, 0, D3D11_MAP_READ, 0, &mapped);
				if (SUCCEEDED(hr)) {
					uint8_t* pDst = (uint8_t*)(CONTIG_BASE + rt.offset);
					uint8_t* pSrc = (uint8_t*)mapped.pData;

					if (desc.Width >= rt.width && desc.Height >= rt.height) {
						if (bMortonDraw || rt.surfaceType == 0x2) {
							// Morton draw or individual swizzled surface:
							// write pixels in NV2A Morton/Z-order.
							swizzle_rect(pSrc, rt.width, rt.height,
								pDst, mapped.RowPitch, rt.bpp);
						} else {
							// Pitch-linear surface in a non-Morton draw:
							// copy row by row respecting Xbox pitch.
							uint32_t xboxRowBytes = rt.width * rt.bpp;
							for (uint32_t row = 0; row < rt.height; row++) {
								memcpy(pDst, pSrc, xboxRowBytes);
								pDst += rt.pitch;
								pSrc += mapped.RowPitch;
							}
						}
					}

					g_pD3DDeviceContext->Unmap(pStaging, 0);
				}
				// Only release if it was a fallback allocation
				if (pStaging != rt.pStagingTex) {
					pStaging->Release();
				}
			}

		}
	}

	// Always flush the requested VB pages to the GPU mirror.
	// Even if no GPU-dirty pages were found above, the VEH fault handler
	// may have already read back RT data into Xbox RAM (clearing the dirty
	// bits in the process). That data must still reach the GPU mirror for
	// the vertex fetch shader to see it. The normal FlushToGPU path skips
	// mid-frame flushes for performance, so this is the only opportunity.
	{
		uint32_t offset = firstPage * PAGE_SIZE_;
		uint32_t uploadSize = (lastPage - firstPage + 1) * PAGE_SIZE_;
		D3D11_BOX box = { offset, 0, 0, offset + uploadSize, 1, 1 };
		g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box,
			(const void*)(CONTIG_BASE + offset), uploadSize, 0);
	}

	return anyDirty;
}

// ******************************************************************
// * Public: Get the mirror buffer SRV
// ******************************************************************
ID3D11ShaderResourceView* CxbxPageTrackerGetMirrorSRV()
{
	return s_pMirrorSRV;
}

ID3D11ShaderResourceView* CxbxPageTrackerGetMirrorSRV_SNORM16x2()
{
	return s_pMirrorSRV_SNORM16x2;
}

ID3D11ShaderResourceView* CxbxPageTrackerGetMirrorSRV_UNORM8x4()
{
	return s_pMirrorSRV_UNORM8x4;
}

ID3D11UnorderedAccessView* CxbxPageTrackerGetMirrorUAV()
{
	return s_pMirrorUAV;
}

// ******************************************************************
// * Public: Upload PGRAPH regs[] to the appended region of the mirror buffer
// ******************************************************************
void CxbxPageTrackerUploadPGRAPH(const void* pRegs, uint32_t size)
{
	if (!s_pMirrorBuf || !g_pD3DDeviceContext || !pRegs || size == 0)
		return;
	if (size > GPU_PGRAPH_SIZE)
		size = GPU_PGRAPH_SIZE;
	D3D11_BOX box = { GPU_PGRAPH_BASE, 0, 0, GPU_PGRAPH_BASE + size, 1, 1 };
	g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box, pRegs, size, 0);
}

// ******************************************************************
// * Public: Upload PFB regs[] to the appended region of the mirror buffer
// ******************************************************************
void CxbxPageTrackerUploadPFB(const void* pRegs, uint32_t size)
{
	if (!s_pMirrorBuf || !g_pD3DDeviceContext || !pRegs || size == 0)
		return;
	if (size > GPU_PFB_SIZE)
		size = GPU_PFB_SIZE;
	D3D11_BOX box = { GPU_PFB_BASE, 0, 0, GPU_PFB_BASE + size, 1, 1 };
	g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box, pRegs, size, 0);
}

// ******************************************************************
// * Public: Upload PVIDEO regs[] to the appended region of the mirror buffer
// ******************************************************************
void CxbxPageTrackerUploadPVIDEO(const void* pRegs, uint32_t size)
{
	if (!s_pMirrorBuf || !g_pD3DDeviceContext || !pRegs || size == 0)
		return;
	if (size > GPU_PVIDEO_SIZE)
		size = GPU_PVIDEO_SIZE;
	D3D11_BOX box = { GPU_PVIDEO_BASE, 0, 0, GPU_PVIDEO_BASE + size, 1, 1 };
	g_pD3DDeviceContext->UpdateSubresource(s_pMirrorBuf, 0, &box, pRegs, size, 0);
}

// ******************************************************************
// * Public: Check if any pages in a texture's range are texture-dirty
// ******************************************************************
bool CxbxPageTrackerIsTextureDirty(uint32_t offset, uint32_t size)
{
	if (size == 0 || offset >= CONTIG_SIZE)
		return false;
	if (offset + size > CONTIG_SIZE)
		size = CONTIG_SIZE - offset;

	uint32_t firstPage = offset / PAGE_SIZE_;
	uint32_t lastPage = (offset + size - 1) / PAGE_SIZE_;

	// Fast path: check whole DWORDs in the middle
	uint32_t firstDW = firstPage >> 5;
	uint32_t lastDW = lastPage >> 5;

	if (firstDW == lastDW) {
		// All pages in a single DWORD — build a mask for [firstPage..lastPage]
		uint32_t lo = firstPage & 31;
		uint32_t hi = lastPage & 31;
		uint32_t mask = ((2u << hi) - 1) & ~((1u << lo) - 1);
		return (s_TextureDirtyBitmap[firstDW] & mask) != 0;
	}

	// Check partial first DWORD
	{
		uint32_t lo = firstPage & 31;
		uint32_t mask = ~((1u << lo) - 1); // bits [lo..31]
		if (s_TextureDirtyBitmap[firstDW] & mask)
			return true;
	}

	// Check full DWORDs in the middle — AVX2 path tests 256 pages (1 MiB) per iteration
	uint32_t midStart = firstDW + 1;
	uint32_t midEnd = lastDW;

	if (s_bHasAVX2 && midEnd - midStart >= 8) {
		// Advance to 8-DWORD (32-byte) alignment within the bitmap
		while (midStart < midEnd && (midStart & 7)) {
			if (s_TextureDirtyBitmap[midStart] != 0)
				return true;
			midStart++;
		}
		// Bulk AVX2 scan: 8 DWORDs (256 bits = 256 pages = 1 MiB) per test
		const uint32_t* base = const_cast<const uint32_t*>(s_TextureDirtyBitmap);
		while (midStart + 8 <= midEnd) {
			__m256i chunk = _mm256_load_si256((const __m256i*)&base[midStart]);
			if (!_mm256_testz_si256(chunk, chunk))
				return true;
			midStart += 8;
		}
	}

	// Scalar tail (or full scalar path when AVX2 unavailable)
	for (uint32_t dw = midStart; dw < midEnd; dw++) {
		if (s_TextureDirtyBitmap[dw] != 0)
			return true;
	}

	// Check partial last DWORD
	{
		uint32_t hi = lastPage & 31;
		uint32_t mask = (2u << hi) - 1; // bits [0..hi]
		if (s_TextureDirtyBitmap[lastDW] & mask)
			return true;
	}

	return false;
}

// ******************************************************************
// * Public: Clear texture-dirty bits after host texture upload
// ******************************************************************
void CxbxPageTrackerClearTextureDirty(uint32_t offset, uint32_t size)
{
	if (size == 0 || offset >= CONTIG_SIZE)
		return;
	if (offset + size > CONTIG_SIZE)
		size = CONTIG_SIZE - offset;

	uint32_t firstPage = offset / PAGE_SIZE_;
	uint32_t lastPage = (offset + size - 1) / PAGE_SIZE_;

	uint32_t firstDW = firstPage >> 5;
	uint32_t lastDW = lastPage >> 5;

	if (firstDW == lastDW) {
		uint32_t lo = firstPage & 31;
		uint32_t hi = lastPage & 31;
		uint32_t mask = ((2u << hi) - 1) & ~((1u << lo) - 1);
		_InterlockedAnd((volatile long*)&s_TextureDirtyBitmap[firstDW], (long)~mask);
		return;
	}

	// Clear partial first DWORD
	{
		uint32_t lo = firstPage & 31;
		uint32_t mask = ~((1u << lo) - 1);
		_InterlockedAnd((volatile long*)&s_TextureDirtyBitmap[firstDW], (long)~mask);
	}

	// Clear full DWORDs in the middle
	for (uint32_t dw = firstDW + 1; dw < lastDW; dw++) {
		s_TextureDirtyBitmap[dw] = 0;
	}

	// Clear partial last DWORD
	{
		uint32_t hi = lastPage & 31;
		uint32_t mask = (2u << hi) - 1;
		_InterlockedAnd((volatile long*)&s_TextureDirtyBitmap[lastDW], (long)~mask);
	}
}

// ******************************************************************
// * Sync committed tiled pages in a range back to contiguous memory
// * (without decommitting — pages stay committed for future writes)
// ******************************************************************
void CxbxSyncTiledRangeToContiguous(uint32_t startOffset, uint32_t size)
{
	if (!s_bHasTiledPages || size == 0 || startOffset >= CONTIG_SIZE)
		return;

	if (startOffset + size > CONTIG_SIZE)
		size = CONTIG_SIZE - startOffset;

	uint32_t firstPage = startOffset / PAGE_SIZE_;
	uint32_t lastPage = (startOffset + size - 1) / PAGE_SIZE_;

	for (uint32_t pageIdx = firstPage; pageIdx <= lastPage; pageIdx++) {
		uint32_t dw = pageIdx >> 5;
		uint32_t bit = pageIdx & 31;
		if (s_TiledCommittedBitmap[dw] & (1u << bit)) {
			uint32_t offset = pageIdx * PAGE_SIZE_;
			memcpy((void*)(CONTIG_BASE + offset),
			       (void*)(TILED_BASE + offset), PAGE_SIZE_);
		}
	}
}

