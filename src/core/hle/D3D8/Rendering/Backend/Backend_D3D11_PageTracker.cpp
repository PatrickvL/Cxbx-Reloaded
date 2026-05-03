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

#include <algorithm>
#include <cstring>

// Wine may not reliably support MEM_WRITE_WATCH / GetWriteWatch.
// When running on Wine, always do a full memcpy instead.
static bool s_bWineFallback = false;

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

// ******************************************************************
// * GPU-dirty bitmap (1 bit per 4 KB page) — set when RT writes here
// ******************************************************************
static uint32_t s_GpuDirtyBitmap[BITMAP_DWORDS] = {};

// ******************************************************************
// * Texture-dirty bitmap — set when CPU-written pages are flushed to
// * the GPU mirror. Cleared per-texture after host upload completes.
// * This gates texture re-upload: if no texture-dirty pages overlap
// * a texture's address range, the host texture is still valid.
// * Applies to all texture types (swizzled, linear, compressed).
// ******************************************************************
static uint32_t s_TextureDirtyBitmap[BITMAP_DWORDS] = {};

// ******************************************************************
// * Tiled committed bitmap — tracks which 0xF0 pages are committed
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
// * GPU mirror buffer (64 MiB ByteAddressBuffer + typed SRV views)
// ******************************************************************
static ID3D11Buffer*             s_pMirrorBuf = nullptr;
static ID3D11ShaderResourceView* s_pMirrorSRV = nullptr;
static ID3D11ShaderResourceView* s_pMirrorSRV_SNORM16x2 = nullptr; // R16G16_SNORM typed view
static ID3D11ShaderResourceView* s_pMirrorSRV_UNORM8x4 = nullptr;  // R8G8B8A8_UNORM typed view

// ******************************************************************
// * VEH handle
// ******************************************************************
static void* s_hVEH = nullptr;

// ******************************************************************
// * Static buffer for GetWriteWatch results (16384 pointers)
// ******************************************************************
static PVOID s_WriteWatchPages[PAGE_COUNT];

// ******************************************************************
// * Bitmap helpers
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
// * Called by PageTrackerVEH below.
// ******************************************************************

// Try to handle an access violation in the contiguous or tiled region.
// Returns true if the fault was handled (page committed/restored).
static bool CxbxPageTrackerHandleFault(void* faultAddress, bool isWrite)
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

		if (TestBit(s_GpuDirtyBitmap, pageIdx)) {
			// TODO: Trigger readback from D3D11 render target into Xbox memory.
			// For now, clear the flag and restore access.
			ClearBit(s_GpuDirtyBitmap, pageIdx);

			DWORD oldProtect;
			VirtualProtect((void*)(CONTIG_BASE + pageIdx * PAGE_SIZE_),
				PAGE_SIZE_, PAGE_READWRITE, &oldProtect);
			return true;
		}
		return false;
	}

	return false;
}

// ******************************************************************
// * VEH handler — GPU-dirty faults + tiled memory redirect
// ******************************************************************
static long WINAPI PageTrackerVEH(EXCEPTION_POINTERS* e)
{
	if (e->ExceptionRecord->ExceptionCode != EXCEPTION_ACCESS_VIOLATION)
		return EXCEPTION_CONTINUE_SEARCH;

	uintptr_t faultAddr = (uintptr_t)e->ExceptionRecord->ExceptionInformation[1];
	bool isWrite = (e->ExceptionRecord->ExceptionInformation[0] == 1);

	if (CxbxPageTrackerHandleFault((void*)faultAddr, isWrite))
		return EXCEPTION_CONTINUE_EXECUTION;

	return EXCEPTION_CONTINUE_SEARCH;
}

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
			_BitScanForward(&pos, bits);
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
	memset(s_GpuDirtyBitmap, 0, sizeof(s_GpuDirtyBitmap));
	memset(s_TiledCommittedBitmap, 0, sizeof(s_TiledCommittedBitmap));
	// All pages start texture-dirty so the first deswizzle for each texture is triggered
	memset(s_TextureDirtyBitmap, 0xFF, sizeof(s_TextureDirtyBitmap));

	// Create 64 MiB GPU mirror buffer (DYNAMIC ByteAddressBuffer with SRV)
	D3D11_BUFFER_DESC desc = {};
	desc.ByteWidth = CONTIG_SIZE;
	desc.Usage = D3D11_USAGE_DYNAMIC;
	desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
	desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	desc.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_ALLOW_RAW_VIEWS;

	HRESULT hr = g_pD3DDevice->CreateBuffer(&desc, nullptr, &s_pMirrorBuf);
	if (FAILED(hr)) {
		EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to create mirror buffer (hr=0x%08X)", hr);
		return;
	}

	// Create raw buffer SRV
	D3D11_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
	srvDesc.Format = DXGI_FORMAT_R32_TYPELESS;
	srvDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFEREX;
	srvDesc.BufferEx.FirstElement = 0;
	srvDesc.BufferEx.NumElements = CONTIG_SIZE / 4; // 16M elements of R32
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

	// Reset write-watch BEFORE the initial upload. This ensures that any writes
	// occurring concurrently (from the Xbox title thread) during the memcpy will
	// have their dirty bits preserved and picked up by the first flush. If we
	// reset AFTER the memcpy, writes that happen between the memcpy passing a page
	// and the reset would have their dirty bits cleared, leaving stale data in the
	// mirror with no way to detect it later.
	ResetWriteWatch((PVOID)CONTIG_BASE, CONTIG_SIZE);

	// Initial full upload of contiguous memory to GPU mirror
	{
		D3D11_MAPPED_SUBRESOURCE mapped = {};
		hr = g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
		if (SUCCEEDED(hr)) {
			memcpy(mapped.pData, (void*)CONTIG_BASE, CONTIG_SIZE);
			g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
		}
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

	// Register VEH for GPU-dirty page faults and tiled memory redirect
	s_hVEH = AddVectoredExceptionHandler(1, PageTrackerVEH);
	if (!s_hVEH) {
		EmuLog(LOG_LEVEL::WARNING, "PageTrackerInit: Failed to register VEH");
	}

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
			_BitScanForward(&pos, bits);
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
			_BitScanForward(&pos, bits);
			bits &= bits - 1;
			uint32_t offset = (dw * 32 + pos) * PAGE_SIZE_;
			VirtualFree((void*)(TILED_BASE + offset), PAGE_SIZE_, MEM_DECOMMIT);
		}
	}

	if (s_pMirrorSRV_UNORM8x4) { s_pMirrorSRV_UNORM8x4->Release(); s_pMirrorSRV_UNORM8x4 = nullptr; }
	if (s_pMirrorSRV_SNORM16x2) { s_pMirrorSRV_SNORM16x2->Release(); s_pMirrorSRV_SNORM16x2 = nullptr; }
	if (s_pMirrorSRV) { s_pMirrorSRV->Release(); s_pMirrorSRV = nullptr; }
	if (s_pMirrorBuf) { s_pMirrorBuf->Release(); s_pMirrorBuf = nullptr; }

	memset(s_GpuDirtyBitmap, 0, sizeof(s_GpuDirtyBitmap));
	memset(s_TiledCommittedBitmap, 0, sizeof(s_TiledCommittedBitmap));
	memset(s_TextureDirtyBitmap, 0, sizeof(s_TextureDirtyBitmap));
}

// ******************************************************************
// * Public: Flush CPU-dirty pages to GPU mirror
// ******************************************************************
uint32_t CxbxPageTrackerFlushToGPU()
{
	if (!s_pMirrorBuf)
		return 0;

	// Sync any committed tiled pages back to contiguous memory.
	// The memcpy writes to 0x80 are automatically tracked by write-watch.
	SyncTiledPagesBack();

	// Wine / broken-write-watch fallback: always do a full upload.
	// Use DISCARD only on the first flush of the frame to avoid orphaning
	// the buffer while earlier draws in this frame are still referencing it.
	if (s_bWineFallback) {
		memset(s_TextureDirtyBitmap, 0xFF, sizeof(s_TextureDirtyBitmap));
		D3D11_MAP mapType = s_bFirstFlushOfFrame
			? D3D11_MAP_WRITE_DISCARD : D3D11_MAP_WRITE_NO_OVERWRITE;
		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, mapType, 0, &mapped);
		if (SUCCEEDED(hr)) {
			memcpy(mapped.pData, (void*)CONTIG_BASE, CONTIG_SIZE);
			g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
		}
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
		memset(s_TextureDirtyBitmap, 0xFF, sizeof(s_TextureDirtyBitmap));
	} else {
		for (ULONG_PTR i = 0; i < count; i++) {
			uint32_t pageIdx = (uint32_t)((uintptr_t)s_WriteWatchPages[i] - CONTIG_BASE) / PAGE_SIZE_;
			SetBit(s_TextureDirtyBitmap, pageIdx);
		}
	}

	// Many pages dirty AND first flush of frame: safe to DISCARD + full memcpy.
	// DISCARD orphans the GPU buffer — any draw calls issued earlier in this frame
	// that haven't completed yet would read from the orphaned (old) buffer, which is
	// fine because DISCARD gives us a fresh allocation. But mid-frame DISCARD would
	// lose updates from earlier flushes, so we only allow it on the first flush.
	if (count > PAGE_COUNT / 4 && s_bFirstFlushOfFrame) {
		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped);
		if (SUCCEEDED(hr)) {
			memcpy(mapped.pData, (void*)CONTIG_BASE, CONTIG_SIZE);
			g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
		}
		s_bFirstFlushOfFrame = false;
	} else {
		// Incremental update: NO_OVERWRITE preserves prior flush data in the same frame.
		// Coalesce consecutive dirty pages into contiguous runs to minimize memcpy calls
		// and maximize throughput (large memcpy uses REP MOVSB / AVX at full bandwidth).

		// GetWriteWatch returns pages in ascending address order per MSDN, but sort
		// defensively to guarantee correctness of the coalescing algorithm.
		std::sort(s_WriteWatchPages, s_WriteWatchPages + count);

		D3D11_MAPPED_SUBRESOURCE mapped = {};
		HRESULT hr = g_pD3DDeviceContext->Map(s_pMirrorBuf, 0, D3D11_MAP_WRITE_NO_OVERWRITE, 0, &mapped);
		if (SUCCEEDED(hr)) {
			uint8_t* pDst = (uint8_t*)mapped.pData;

			// Walk the sorted page list and merge consecutive pages into runs
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
				memcpy(pDst + offset, (const void*)startAddr, runBytes);

				runStart = runEnd;
			}

			g_pD3DDeviceContext->Unmap(s_pMirrorBuf, 0);
		}
		s_bFirstFlushOfFrame = false;
	}

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
		SetBit(s_GpuDirtyBitmap, p);
		// Set PAGE_NOACCESS so CPU reads/writes fault for readback
		DWORD oldProtect;
		VirtualProtect((void*)(CONTIG_BASE + p * PAGE_SIZE_),
			PAGE_SIZE_, PAGE_NOACCESS, &oldProtect);
	}
}

// ******************************************************************
// * Public: Check if a page is GPU-dirty
// ******************************************************************
bool CxbxPageTrackerIsGPUDirty(uint32_t pageIndex)
{
	if (pageIndex >= PAGE_COUNT) return false;
	return TestBit(s_GpuDirtyBitmap, pageIndex);
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
		ClearBit(s_GpuDirtyBitmap, p);
	}
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

	// Check full DWORDs in the middle
	for (uint32_t dw = firstDW + 1; dw < lastDW; dw++) {
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
		s_TextureDirtyBitmap[firstDW] &= ~mask;
		return;
	}

	// Clear partial first DWORD
	{
		uint32_t lo = firstPage & 31;
		uint32_t mask = ~((1u << lo) - 1);
		s_TextureDirtyBitmap[firstDW] &= ~mask;
	}

	// Clear full DWORDs in the middle
	for (uint32_t dw = firstDW + 1; dw < lastDW; dw++) {
		s_TextureDirtyBitmap[dw] = 0;
	}

	// Clear partial last DWORD
	{
		uint32_t hi = lastPage & 31;
		uint32_t mask = (2u << hi) - 1;
		s_TextureDirtyBitmap[lastDW] &= ~mask;
	}
}

