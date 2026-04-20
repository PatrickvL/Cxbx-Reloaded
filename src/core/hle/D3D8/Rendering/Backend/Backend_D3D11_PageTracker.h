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
#ifndef BACKEND_D3D11_PAGE_TRACKER_H
#define BACKEND_D3D11_PAGE_TRACKER_H

#ifdef CXBX_USE_D3D11

// Backend_D3D11_PageTracker — Dirty page tracking for the 64 MiB contiguous memory mirror.
//
// Design:
//   CPU→GPU: MEM_WRITE_WATCH based dirty page tracking. The contiguous region
//            (0x80000000) is allocated with VirtualAlloc + MEM_WRITE_WATCH,
//            providing zero-overhead write tracking via hardware PTE dirty bits.
//            GetWriteWatch/ResetWriteWatch replaces the VEH write-fault approach.
//
//   GPU→CPU: Software bitmap. When a render target writes to contiguous memory,
//            those pages are marked GPU-dirty and set PAGE_NOACCESS. On CPU access,
//            VEH triggers readback from the D3D11 RT and restores access.
//
//   Tiled memory (0xF0000000): Allocated as MEM_RESERVE + PAGE_NOACCESS (no file
//            mapping alias). VEH commits pages on demand, copying from 0x80000000.
//            At flush time, committed tiled pages are synced back and decommitted.
//
// The 64 MiB contiguous region at CONTIGUOUS_MEMORY_BASE (0x80000000) maps to a
// GPU-side ByteAddressBuffer. Xbox physical addresses are masked to 27 bits
// (& 0x07FFFFFF) to index directly into this buffer.
//
// Page granularity: 4 KB (PAGE_SIZE). Bitmap size: 64 MiB / 4 KB = 16384 bits = 2 KB.

#include <cstdint>

// Number of pages in the contiguous region (64 MiB / 4 KB)
static constexpr uint32_t CXBX_CONTIG_PAGE_COUNT = (64 * 1024 * 1024) / 4096; // 16384

// ******************************************************************
// * Initialization / Shutdown
// ******************************************************************

// Initialize page tracking: allocate bitmaps, create the GPU mirror buffer (64 MiB
// ByteAddressBuffer), register VEH handler for write-fault tracking.
// Called once during D3D11 device initialization.
void CxbxPageTrackerInit();

// Release page tracking resources: unregister VEH, release GPU buffer, free bitmaps.
void CxbxPageTrackerShutdown();

// ******************************************************************
// * CPU→GPU: Write tracking (MEM_WRITE_WATCH)
// ******************************************************************

// Flush all CPU-dirty pages to the GPU mirror buffer.
// Uses GetWriteWatch to retrieve modified pages, copies them to the GPU
// mirror, and atomically resets the write-watch.
// Also syncs any committed tiled pages back to contiguous memory.
// Returns the number of pages flushed.
uint32_t CxbxPageTrackerFlushToGPU();

// Check if any pages are CPU-dirty (quick early-out for draw path).
bool CxbxPageTrackerHasDirtyPages();

// ******************************************************************
// * GPU→CPU: Render target tracking
// ******************************************************************

// Mark a range of contiguous memory as GPU-dirty (written by a render target).
// Called when SetRenderTarget binds an RT that maps into the contiguous region.
// startOffset: byte offset from CONTIGUOUS_MEMORY_BASE (masked to 27 bits).
// size: byte size of the render target surface.
void CxbxPageTrackerMarkGPUDirty(uint32_t startOffset, uint32_t size);

// Check if a specific contiguous page is GPU-dirty.
// Used by CPU read paths to trigger readback before accessing the data.
bool CxbxPageTrackerIsGPUDirty(uint32_t pageIndex);

// Clear GPU-dirty flags for a range after readback completes.
void CxbxPageTrackerClearGPUDirty(uint32_t startOffset, uint32_t size);

// ******************************************************************
// * GPU mirror buffer access (for draw path)
// ******************************************************************

// Get the SRV for the 64 MiB contiguous mirror ByteAddressBuffer (t0 binding).
struct ID3D11ShaderResourceView;
ID3D11ShaderResourceView* CxbxPageTrackerGetMirrorSRV();

#endif // CXBX_USE_D3D11
#endif // BACKEND_D3D11_PAGE_TRACKER_H
