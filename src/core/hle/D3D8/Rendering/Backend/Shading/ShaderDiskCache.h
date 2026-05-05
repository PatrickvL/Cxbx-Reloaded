#pragma once
// ShaderDiskCache.h — Persistent disk cache for JIT-compiled shader bytecode
//
// Saves compiled shader blobs to disk so subsequent runs skip D3DCompile.
// Per-game directory: <DataPath>\ShaderCache\<TitleID>-<GameName>\
// Files: <hash>.cso (raw DXBC bytecode)

#include <d3dcompiler.h>
#include <cstdint>

namespace ShaderDiskCache {

// Try to load a cached shader blob from disk.
// Returns a new-ref'd ID3DBlob on success, or nullptr on miss.
ID3DBlob* TryLoad(uint64_t hash);

// Queue a compiled shader blob for background save to disk.
// Takes ownership of nothing — copies the data internally.
void Save(uint64_t hash, ID3DBlob* pBlob);

// Flush pending writes and release resources. Call before process exit.
void Shutdown();

} // namespace ShaderDiskCache
