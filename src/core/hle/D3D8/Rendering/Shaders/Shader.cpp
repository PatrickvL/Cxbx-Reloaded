// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++ and C#: http://www.viva64.com
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
// *  2020 PatrickvL
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::VTXSH // TODO : Introduce generic HLSL logging

#include <d3dcompiler.h>
#include "Shader.h"
#include "common/FilePaths.hpp" // For szFilePath_CxbxReloaded_Exe
#include "core\kernel\init\CxbxKrnl.h" // LOG_TEST_CASE
#include "core\kernel\support\Emu.h" // EmuLog

#include <filesystem>
#include <fstream>
#include <array>
#include <thread>
#include <mutex>
#include <queue>
#include <atomic>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include "common\util\hasher.h" // For ComputeHash
//#include <sstream>

// Function pointer type matching D3DCompile's signature
typedef HRESULT(WINAPI *PFN_D3DCOMPILE)(
	LPCVOID pSrcData, SIZE_T SrcDataSize, LPCSTR pSourceName,
	const D3D_SHADER_MACRO *pDefines, ID3DInclude *pInclude,
	LPCSTR pEntrypoint, LPCSTR pTarget, UINT Flags1, UINT Flags2,
	ID3DBlob **ppCode, ID3DBlob **ppErrorMsgs);

// Dynamically resolve D3DCompile from a native d3dcompiler_47.dll placed
// next to the executable. This bypasses Wine's builtin d3dcompiler which
// cannot handle complex shaders like the VS/PS interpreter ubershaders.
// Falls back to the linked D3DCompile if the native DLL is not found.
static PFN_D3DCOMPILE GetD3DCompileFunc()
{
	static PFN_D3DCOMPILE s_func = nullptr;
	static bool s_tried = false;
	if (!s_tried) {
		s_tried = true;
		// Build an absolute path to d3dcompiler_47.dll next to our executable.
		// Using an absolute path forces Wine to load the native DLL instead of
		// its builtin, which cannot handle complex ubershaders.
		char exePath[MAX_PATH] = {};
		GetModuleFileNameA(nullptr, exePath, MAX_PATH);
		std::string dllPath(exePath);
		auto lastSlash = dllPath.find_last_of("\\/");
		if (lastSlash != std::string::npos)
			dllPath = dllPath.substr(0, lastSlash + 1);
		dllPath += "d3dcompiler_47.dll";
		HMODULE hMod = LoadLibraryA(dllPath.c_str());
		if (hMod) {
			s_func = (PFN_D3DCOMPILE)GetProcAddress(hMod, "D3DCompile");
			if (s_func) {
				EmuLog(LOG_LEVEL::INFO, "Loaded native D3DCompile from %s", dllPath.c_str());
			}
		}
		if (!s_func) {
			// Fall back to linked version
			s_func = &D3DCompile;
		}
	}
	return s_func;
}

ShaderSources g_ShaderSources;

// ============================================================================
// Disk-based shader bytecode cache
// ============================================================================

// Shader bytecode magic validation
// D3D9 SM1-3: vertex shaders start with 0xFFFExxxx, pixel shaders with 0xFFFFxxxx
// D3D10+ SM4+: DXBC container starts with 0x43425844 ("DXBC")
static bool IsValidShaderBytecode(uint32_t magic)
{
	return (magic >> 16) == 0xFFFE  // D3D9 vertex shader
	    || (magic >> 16) == 0xFFFF  // D3D9 pixel shader
	    || magic == 0x43425844;     // DXBC container (SM4+)
}

static std::string g_ShaderCacheDir;
static std::string g_SharedShaderCacheDir;
static std::atomic<int> g_CacheHits{0};
static std::atomic<int> g_CacheMisses{0};
static std::atomic<int> g_CacheSaves{0};
static std::atomic<int> g_CacheLoadErrors{0};

// Log file for shader cache (since emulation process may not have a console)
static FILE* g_ShaderCacheLogFile = nullptr;
static std::mutex g_LogMutex;

static void ShaderCacheLog(const char* fmt, ...)
{
	if (!g_ShaderCacheLogFile) return;
	std::lock_guard<std::mutex> lock(g_LogMutex);
	va_list args;
	va_start(args, fmt);
	vfprintf(g_ShaderCacheLogFile, fmt, args);
	va_end(args);
	fflush(g_ShaderCacheLogFile);
}

// Background save queue
static std::mutex g_SaveQueueMutex;
static std::queue<std::pair<std::string, std::vector<uint8_t>>> g_SaveQueue;
static std::thread g_SaveThread;
static std::atomic<bool> g_SaveThreadRunning{false};

static void ShaderCacheSaveWorker()
{
	while (g_SaveThreadRunning) {
		std::pair<std::string, std::vector<uint8_t>> item;
		bool hasItem = false;
		{
			std::lock_guard<std::mutex> lock(g_SaveQueueMutex);
			if (!g_SaveQueue.empty()) {
				item = std::move(g_SaveQueue.front());
				g_SaveQueue.pop();
				hasItem = true;
			}
		}

		if (!hasItem) {
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			continue;
		}

		FILE* fp = fopen(item.first.c_str(), "wb");
		if (fp) {
			fwrite(item.second.data(), 1, item.second.size(), fp);
			fclose(fp);
			g_CacheSaves++;
			ShaderCacheLog("SAVED %s (%zu bytes)\n", item.first.c_str(), item.second.size());
		} else {
			ShaderCacheLog("ERROR could not write %s (errno=%d)\n", item.first.c_str(), errno);
		}
	}
}

void ShaderCacheShutdown()
{
	// Signal the save thread to stop and drain whatever is left in the queue.
	g_SaveThreadRunning = false;

	// Drain anything still in the queue on this thread
	while (true) {
		std::pair<std::string, std::vector<uint8_t>> item;
		{
			std::lock_guard<std::mutex> lock(g_SaveQueueMutex);
			if (g_SaveQueue.empty()) break;
			item = std::move(g_SaveQueue.front());
			g_SaveQueue.pop();
		}
		FILE* fp = fopen(item.first.c_str(), "wb");
		if (fp) {
			fwrite(item.second.data(), 1, item.second.size(), fp);
			fclose(fp);
			ShaderCacheLog("SHUTDOWN-SAVE %s (%zu bytes)\n", item.first.c_str(), item.second.size());
		}
	}

	ShaderCacheLog("ShaderCache shutdown: hits=%d misses=%d saves=%d errors=%d\n",
		g_CacheHits.load(), g_CacheMisses.load(), g_CacheSaves.load(), g_CacheLoadErrors.load());

	if (g_ShaderCacheLogFile) {
		fclose(g_ShaderCacheLogFile);
		g_ShaderCacheLogFile = nullptr;
	}
}

// ============================================================================
// Shared in-memory shader caches
// ============================================================================

// Results from background async compiles
static std::mutex g_AsyncMutex;
static std::unordered_map<uint64_t, ID3DBlob*> g_AsyncResults;
static std::unordered_set<uint64_t> g_AsyncInFlight;

// Blob cache populated from disk loads and async compiles.
// Checked first on every EmuCompileShader call — pure hash-map lookup, no file I/O.
static std::unordered_map<uint64_t, ID3DBlob*> g_MemCache;

// Forward declarations for functions used by EnsureShaderCacheDir
static void PreloadShaderCache();
static void PreloadShaderCacheFrom(const std::string& dir);

// Returns true if the shared shader cache dir (00000000-All) is ready to use.
// This dir does not require a game certificate — only g_DataFilePath.
static bool EnsureSharedShaderCacheDir()
{
	if (!g_SharedShaderCacheDir.empty()) return true;
	if (g_DataFilePath.empty()) return false;

	g_SharedShaderCacheDir = g_DataFilePath + "\\ShaderCache\\00000000-All";
	std::error_code ec;
	if (!std::filesystem::exists(g_SharedShaderCacheDir)) {
		std::filesystem::create_directories(g_SharedShaderCacheDir, ec);
		if (ec) {
			g_SharedShaderCacheDir.clear();
			return false;
		}
	}

	// Create Dumped and Replacements subdirectories
	for (const char* sub : { "Dumped", "Replacements" }) {
		std::string subPath = g_SharedShaderCacheDir + "\\" + sub;
		if (!std::filesystem::exists(subPath)) {
			std::filesystem::create_directory(subPath, ec);
		}
	}

	// Preload all .cso files from the shared cache into memory
	PreloadShaderCacheFrom(g_SharedShaderCacheDir);

	ShaderCacheLog("Shared shader cache initialized: %s\n", g_SharedShaderCacheDir.c_str());
	return true;
}

// Returns true if cache dir is ready to use
static bool EnsureShaderCacheDir()
{
	if (!g_ShaderCacheDir.empty()) return true;

	// g_DataFilePath may not be set yet during early init
	if (g_DataFilePath.empty()) return false;

	// Need game certificate to create per-game directory
	if (!g_pCertificate) return false;

	// Build per-game cache dir: ShaderCache\<TitleID>-<GameName>
	// e.g. ShaderCache\4D530004-Halo
	char titleIdStr[16];
	snprintf(titleIdStr, sizeof(titleIdStr), "%08X", g_pCertificate->dwTitleId);

	// Get ASCII game title and sanitize for filesystem use
	std::string gameName;
	if (CxbxKrnl_Xbe && CxbxKrnl_Xbe->m_szAsciiTitle[0]) {
		gameName = CxbxKrnl_Xbe->m_szAsciiTitle;
		// Remove characters invalid in directory names
		for (char& c : gameName) {
			if (c == '\\' || c == '/' || c == ':' || c == '*' ||
				c == '?' || c == '"' || c == '<' || c == '>' || c == '|')
				c = '_';
		}
		// Trim trailing spaces
		while (!gameName.empty() && gameName.back() == ' ')
			gameName.pop_back();
	}

	std::string gameDir = std::string(titleIdStr);
	if (!gameName.empty()) {
		gameDir += "-" + gameName;
	}

	g_ShaderCacheDir = g_DataFilePath + "\\ShaderCache\\" + gameDir;
	std::error_code ec;
	if (!std::filesystem::exists(g_ShaderCacheDir)) {
		std::filesystem::create_directories(g_ShaderCacheDir, ec);
		if (ec) {
			// Failed to create — reset so we retry next time
			g_ShaderCacheDir.clear();
			return false;
		}
	}

	// Open log file in the per-game shader cache dir
	std::string logPath = g_ShaderCacheDir + "\\shader_cache.log";
	g_ShaderCacheLogFile = fopen(logPath.c_str(), "wt");
	ShaderCacheLog("ShaderCache initialized: %s\n", g_ShaderCacheDir.c_str());
	ShaderCacheLog("g_DataFilePath = %s\n", g_DataFilePath.c_str());
	ShaderCacheLog("TitleID = %s, GameName = %s\n", titleIdStr, gameName.c_str());

	// Start background save thread
	if (!g_SaveThreadRunning) {
		g_SaveThreadRunning = true;
		g_SaveThread = std::thread(ShaderCacheSaveWorker);
		g_SaveThread.detach();
	}

	// Create Dumped and Replacements subdirectories
	// Dumped: original HLSL sources written here for user inspection and patching
	// Replacements: user places modified HLSL here (same filename as Dumped) to override at runtime
	for (const char* sub : { "Dumped", "Replacements" }) {
		std::string subPath = g_ShaderCacheDir + "\\" + sub;
		if (!std::filesystem::exists(subPath)) {
			std::filesystem::create_directory(subPath, ec);
		}
	}

	// Also ensure the shared shader cache dir exists
	EnsureSharedShaderCacheDir();

	// Preload all .cso files into memory now that the cache dir is known
	PreloadShaderCacheFrom(g_ShaderCacheDir);

	return true;
}

static std::string GetShaderCachePath(uint64_t hash, const std::string& cacheDir)
{
	char filename[32];
	snprintf(filename, sizeof(filename), "%016llx.cso", hash);
	return cacheDir + "\\" + filename;
}

// Returns path of the HLSL file for the given hash+profile in 'Dumped' or 'Replacements'.
static std::string GetShaderHlslPath(uint64_t hash, const char* profile, const char* subdir, const std::string& cacheDir)
{
	char filename[80];
	snprintf(filename, sizeof(filename), "%016llx_%s.hlsl", hash, profile);
	return cacheDir + "\\" + subdir + "\\" + filename;
}

// Write the HLSL source to Dumped\ the first time a shader is seen.
static void DumpShaderSource(uint64_t hash, const char* profile, const std::string& hlsl, const std::string& cacheDir)
{
	if (cacheDir.empty()) return;
	std::string path = GetShaderHlslPath(hash, profile, "Dumped", cacheDir);
	if (std::filesystem::exists(path)) return; // already dumped
	std::ofstream f(path);
	if (f.good()) {
		f << hlsl;
		ShaderCacheLog("DUMPED %s\n", path.c_str());
	}
}

// Check whether the user placed a replacement HLSL in Replacements\.
static bool TryLoadReplacementShader(uint64_t hash, const char* profile, std::string& out_hlsl, const std::string& cacheDir)
{
	if (cacheDir.empty()) return false;
	std::string path = GetShaderHlslPath(hash, profile, "Replacements", cacheDir);
	std::ifstream f(path);
	if (!f.good()) return false;
	out_hlsl.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
	ShaderCacheLog("REPLACEMENT loaded %s\n", path.c_str());
	EmuLog(LOG_LEVEL::INFO, "ShaderCache: using replacement shader %s", path.c_str());
	return true;
}

static bool LoadCachedShader(uint64_t hash, ID3DBlob** ppBlob, const std::string& cacheDir)
{
	std::string path = GetShaderCachePath(hash, cacheDir);
	FILE* fp = fopen(path.c_str(), "rb");
	if (!fp) return false;

	fseek(fp, 0, SEEK_END);
	long size = ftell(fp);
	fseek(fp, 0, SEEK_SET);

	if (size < 8) {
		// Minimum: version token (4 bytes) + end token (4 bytes)
		ShaderCacheLog("REJECT %s (too small: %ld bytes)\n", path.c_str(), size);
		fclose(fp);
		g_CacheLoadErrors++;
		return false;
	}

	// Read the magic bytes first to validate
	uint32_t magic = 0;
	fread(&magic, 4, 1, fp);
	if (!IsValidShaderBytecode(magic)) {
		ShaderCacheLog("REJECT %s (bad magic: 0x%08X)\n", path.c_str(), magic);
		fclose(fp);
		g_CacheLoadErrors++;
		// Delete corrupt cache file
		std::error_code ec;
		std::filesystem::remove(path, ec);
		return false;
	}
	fseek(fp, 0, SEEK_SET);

	HRESULT hr = D3DCreateBlob(size, ppBlob);
	if (FAILED(hr)) {
		ShaderCacheLog("ERROR D3DCreateBlob failed for %s (size=%ld, hr=0x%08lX)\n", path.c_str(), size, hr);
		fclose(fp);
		g_CacheLoadErrors++;
		return false;
	}

	size_t readBytes = fread((*ppBlob)->GetBufferPointer(), 1, size, fp);
	fclose(fp);

	if ((long)readBytes != size) {
		ShaderCacheLog("ERROR partial read %s (%zu / %ld bytes)\n", path.c_str(), readBytes, size);
		(*ppBlob)->Release();
		*ppBlob = nullptr;
		g_CacheLoadErrors++;
		return false;
	}

	g_CacheHits++;
	ShaderCacheLog("HIT %016llx (%ld bytes) [hits=%d misses=%d]\n", hash, size, g_CacheHits.load(), g_CacheMisses.load());
	return true;
}

// Scan a shader cache directory and load every .cso into g_MemCache so that
// gameplay never needs a file open.
static void PreloadShaderCacheFrom(const std::string& dir)
{
	if (dir.empty()) return;
	std::error_code ec;
	int preloaded = 0;
	std::lock_guard<std::mutex> lock(g_AsyncMutex);
	for (auto& entry : std::filesystem::directory_iterator(dir, ec)) {
		if (entry.path().extension() != ".cso") continue;
		std::string stem = entry.path().stem().string();
		if (stem.size() != 16) continue;
		uint64_t hash = 0;
		try { hash = std::stoull(stem, nullptr, 16); }
		catch (...) { continue; }
		if (g_MemCache.count(hash)) continue;
		ID3DBlob* pBlob = nullptr;
		if (LoadCachedShader(hash, &pBlob, dir)) {
			g_MemCache[hash] = pBlob;
			preloaded++;
		}
	}
	ShaderCacheLog("Preloaded %d shader(s) into memory cache from %s\n", preloaded, dir.c_str());
	EmuLog(LOG_LEVEL::INFO, "ShaderCache: preloaded %d shader(s) into memory from %s", preloaded, dir.c_str());
}

// Legacy wrapper — preloads from the per-game cache directory.
static void PreloadShaderCache()
{
	PreloadShaderCacheFrom(g_ShaderCacheDir);
}

static void QueueSaveCachedShader(uint64_t hash, ID3DBlob* pBlob, const std::string& cacheDir)
{
	// Validate the blob has DXBC magic before saving
	if (pBlob->GetBufferSize() < 4) {
		ShaderCacheLog("SKIP save %016llx (blob too small: %zu bytes)\n", hash, pBlob->GetBufferSize());
		return;
	}

	uint32_t magic = *reinterpret_cast<const uint32_t*>(pBlob->GetBufferPointer());
	if (!IsValidShaderBytecode(magic)) {
		ShaderCacheLog("SKIP save %016llx (bad magic: 0x%08X)\n", hash, magic);
		return;
	}

	std::string path = GetShaderCachePath(hash, cacheDir);
	std::vector<uint8_t> data(
		static_cast<uint8_t*>(pBlob->GetBufferPointer()),
		static_cast<uint8_t*>(pBlob->GetBufferPointer()) + pBlob->GetBufferSize()
	);

	{
		std::lock_guard<std::mutex> lock(g_SaveQueueMutex);
		g_SaveQueue.push({path, std::move(data)});
	}
	ShaderCacheLog("QUEUED save %016llx (%zu bytes)\n", hash, pBlob->GetBufferSize());
}

// ============================================================================
// Async shader compilation (Dolphin-style)
// ============================================================================

// Precompiled fallback pixel shader (simple white output)
static ID3DBlob* g_FallbackPSBlob = nullptr;
static bool g_FallbacksInitialized = false;

static void EnsureFallbackShaders()
{
	if (g_FallbacksInitialized) return;
	g_FallbacksInitialized = true;

	// Minimal pixel shader: output white (only PS uses async fallback)
	const char* psSrc =
		"float4 main() : SV_Target0 { return float4(1,1,1,1); }\n";
	D3DCompile(psSrc, strlen(psSrc), nullptr, nullptr, nullptr,
		"main", "ps_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL0, 0, &g_FallbackPSBlob, nullptr);

	ShaderCacheLog("Fallback PS compiled: %p (%zu bytes)\n",
		g_FallbackPSBlob, g_FallbackPSBlob ? g_FallbackPSBlob->GetBufferSize() : 0);
}

// Background compile worker — runs real D3DCompile and stores result
static void AsyncCompileWorker(std::string hlsl_str, std::string profile,
	std::string sourceName, uint64_t hash, bool cacheReady, std::string cacheDir)
{
	auto tStart = std::chrono::high_resolution_clock::now();

	ID3DBlob* pResult = nullptr;
	ID3DBlob* pErrors = nullptr;
	// Match the synchronous path: O1, all shaders use DX11 native Texture2D.Sample()
	UINT flags1 = D3DCOMPILE_OPTIMIZATION_LEVEL1;

	auto pfnD3DCompile = GetD3DCompileFunc();

	HRESULT hr = pfnD3DCompile(
		hlsl_str.c_str(), hlsl_str.length(),
		sourceName.empty() ? nullptr : sourceName.c_str(),
		nullptr, D3D_COMPILE_STANDARD_FILE_INCLUDE,
		"main", profile.c_str(), flags1, 0, &pResult, &pErrors);

	if (FAILED(hr)) {
		if (pErrors) { pErrors->Release(); pErrors = nullptr; }
		flags1 = D3DCOMPILE_OPTIMIZATION_LEVEL0;
		hr = pfnD3DCompile(
			hlsl_str.c_str(), hlsl_str.length(),
			sourceName.empty() ? nullptr : sourceName.c_str(),
			nullptr, D3D_COMPILE_STANDARD_FILE_INCLUDE,
			"main", profile.c_str(), flags1, 0, &pResult, &pErrors);
	}

	if (pErrors) { pErrors->Release(); pErrors = nullptr; }

	auto tEnd = std::chrono::high_resolution_clock::now();
	double ms = std::chrono::duration<double, std::milli>(tEnd - tStart).count();

	{
		std::lock_guard<std::mutex> lock(g_AsyncMutex);
		g_AsyncInFlight.erase(hash);
		if (!FAILED(hr) && pResult) {
			g_AsyncResults[hash] = pResult; // takes ownership
			pResult->AddRef();
			g_MemCache[hash] = pResult;     // also in fast in-memory lookup cache
			if (cacheReady) {
				QueueSaveCachedShader(hash, pResult, cacheDir);
			}
			ShaderCacheLog("ASYNC COMPILED %016llx in %.2f ms (profile=%s, blob=%zu bytes)\n",
				hash, ms, profile.c_str(), pResult->GetBufferSize());
		} else {
			ShaderCacheLog("ASYNC COMPILE FAILED %016llx in %.2f ms (profile=%s, hr=0x%08lX)\n",
				hash, ms, profile.c_str(), hr);
		}
	}
}

std::string DebugPrependLineNumbers(std::string shaderString) {
	std::stringstream shader(shaderString);
	auto debugShader = std::stringstream();

	int i = 1;
	for (std::string line; std::getline(shader, line); ) {
		auto lineNumber = std::to_string(i++);
		auto paddedLineNumber = lineNumber.insert(0, 3 - lineNumber.size(), ' ');
		debugShader << "/* " << paddedLineNumber << " */ " << line << "\n";
	}

	return debugShader.str();
}

extern HRESULT EmuCompileShader
(
	std::string hlsl_str,
	const char* shader_profile,
	ID3DBlob** ppHostShader,
	const char* pSourceName,
	bool asyncAllowed,
	bool useSharedCache
)
{
	// Compute a cache key from the HLSL source + shader profile
	std::string cacheInput = hlsl_str + "|" + shader_profile;
	uint64_t cacheHash = ComputeHash(cacheInput.c_str(), cacheInput.size());

	// 0. Fast path: check in-memory cache (no file I/O or mutex beyond the map lookup)
	{
		std::lock_guard<std::mutex> lock(g_AsyncMutex);
		auto memIt = g_MemCache.find(cacheHash);
		if (memIt != g_MemCache.end()) {
			memIt->second->AddRef();
			*ppHostShader = memIt->second;
			return S_OK;
		}
	}

	// 1. Try loading from disk cache (only if cache dir is ready)
	// For shared (generic) shaders, use the 00000000-All directory;
	// for per-game shaders, use the title-specific directory.
	bool cacheReady;
	std::string effectiveCacheDir;
	if (useSharedCache) {
		cacheReady = EnsureSharedShaderCacheDir();
		effectiveCacheDir = g_SharedShaderCacheDir;
		// Also ensure the per-game dir is set up (starts save thread, etc.)
		EnsureShaderCacheDir();
	} else {
		cacheReady = EnsureShaderCacheDir();
		effectiveCacheDir = g_ShaderCacheDir;
	}
	if (cacheReady) {
		auto t0 = std::chrono::high_resolution_clock::now();
		if (LoadCachedShader(cacheHash, ppHostShader, effectiveCacheDir)) {
			auto t1 = std::chrono::high_resolution_clock::now();
			double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
			ShaderCacheLog("LOAD took %.2f ms (profile=%s)\n", ms, shader_profile);
			// Add to memcache so subsequent hits skip disk I/O
			(*ppHostShader)->AddRef();
			{
				std::lock_guard<std::mutex> lock(g_AsyncMutex);
				g_MemCache[cacheHash] = *ppHostShader;
			}
			return S_OK;
		}
	}

	// 1.5. Dump original HLSL + check for user-provided replacement
	DumpShaderSource(cacheHash, shader_profile, hlsl_str, effectiveCacheDir);
	std::string replacementHlsl;
	if (TryLoadReplacementShader(cacheHash, shader_profile, replacementHlsl, effectiveCacheDir)) {
		hlsl_str = std::move(replacementHlsl);
	}

	// 2. Async pixel shader path (only for PS, not VS)
	if (asyncAllowed && shader_profile[0] == 'p') {
		// Ensure fallback shader is ready (outside mutex — D3DCompile is slow first time)
		EnsureFallbackShaders();

		std::unique_lock<std::mutex> lock(g_AsyncMutex);

		// Check in-memory async results (completed background compiles)
		auto it = g_AsyncResults.find(cacheHash);
		if (it != g_AsyncResults.end()) {
			// Async compilation finished — return the real shader
			it->second->AddRef();
			*ppHostShader = it->second;
			ShaderCacheLog("ASYNC HIT %016llx (profile=%s)\n", cacheHash, shader_profile);
			return S_OK;
		}

		if (g_AsyncInFlight.count(cacheHash)) {
			// Already compiling in background — return fallback
			if (g_FallbackPSBlob) {
				g_FallbackPSBlob->AddRef();
				*ppHostShader = g_FallbackPSBlob;
				ShaderCacheLog("ASYNC PENDING %016llx -> fallback (profile=%s)\n", cacheHash, shader_profile);
				return S_FALSE;
			}
			// If no fallback available, fall through to synchronous compile
		} else {
			// Start background compile — release lock before thread creation
			g_AsyncInFlight.insert(cacheHash);
			lock.unlock();

			std::string profileStr = shader_profile;
			std::string sourceStr = pSourceName ? pSourceName : "";
			std::thread(AsyncCompileWorker, hlsl_str, profileStr, sourceStr, cacheHash, cacheReady, effectiveCacheDir).detach();

			// Return fallback shader
			if (g_FallbackPSBlob) {
				g_FallbackPSBlob->AddRef();
				*ppHostShader = g_FallbackPSBlob;
				ShaderCacheLog("ASYNC START %016llx -> fallback (profile=%s)\n", cacheHash, shader_profile);
				return S_FALSE;
			}
			// If no fallback (shouldn't happen), fall through to synchronous
			std::lock_guard<std::mutex> reLock(g_AsyncMutex);
			g_AsyncInFlight.erase(cacheHash);
		}
	} else if (asyncAllowed) {
		// VS path: check async results in case a previous async compile for this hash finished
		std::lock_guard<std::mutex> lock(g_AsyncMutex);
		auto it = g_AsyncResults.find(cacheHash);
		if (it != g_AsyncResults.end()) {
			it->second->AddRef();
			*ppHostShader = it->second;
			ShaderCacheLog("ASYNC HIT %016llx (profile=%s)\n", cacheHash, shader_profile);
			return S_OK;
		}
	}

	// 3. Synchronous compilation (for vertex shaders or when async fallback isn't available)
	ID3DBlob* pErrors = nullptr;
	ID3DBlob* pErrorsCompatibility = nullptr;
	HRESULT             hRet = 0;

	*ppHostShader = nullptr;

	EmuLog(LOG_LEVEL::DEBUG, "--- HLSL conversion ---");
	EmuLog(LOG_LEVEL::DEBUG, DebugPrependLineNumbers(hlsl_str).c_str());
	EmuLog(LOG_LEVEL::DEBUG, "-----------------------");


	auto tCompileStart = std::chrono::high_resolution_clock::now();
	UINT flags1 = D3DCOMPILE_OPTIMIZATION_LEVEL3;

	// Use O1 for all Cxbx shaders:
	// - Vertex shaders include FetchAllAttributes() which unrolls a
	//   20-format switch × 16 attributes; at O3 the HLSL compiler spends
	//   exponential time optimising this.
	// - The RC interpreter pixel shader (718 lines with dynamic loops/switches)
	//   can take minutes to compile at O3 under Wine/vkd3d-shader.
	// O1 compiles in seconds and produces correct code for both.
	flags1 = (flags1 & ~D3DCOMPILE_OPTIMIZATION_LEVEL3) | D3DCOMPILE_OPTIMIZATION_LEVEL1;

	auto pfnD3DCompile = GetD3DCompileFunc();

	hRet = pfnD3DCompile(
		hlsl_str.c_str(),
		hlsl_str.length(),
		pSourceName,
		nullptr, // pDefines
		D3D_COMPILE_STANDARD_FILE_INCLUDE,
		"main",
		shader_profile,
		flags1,
		0,
		ppHostShader,
		&pErrors
	);
	if (FAILED(hRet)) {
		if (pErrors) {
			EmuLog(LOG_LEVEL::WARNING, "Shader compile failed: %s", (char*)(pErrors->GetBufferPointer()));
			pErrors->Release();
			pErrors = nullptr;
		} else {
			EmuLog(LOG_LEVEL::WARNING, "Shader compile failed. Recompiling in compatibility mode");
		}
		// Retry at O0. Avoid AVOID_FLOW_CONTROL — it flattens
		// the VS/PS interpreter loops and produces wrong results.
		flags1 = D3DCOMPILE_OPTIMIZATION_LEVEL0;
		hRet = pfnD3DCompile(
			hlsl_str.c_str(),
			hlsl_str.length(),
			pSourceName,
			nullptr, // pDefines
			D3D_COMPILE_STANDARD_FILE_INCLUDE, // pInclude // TODO precompile x_* HLSL functions?
			"main", // shader entry poiint
			shader_profile,
			flags1, // flags1
			0, // flags2
			ppHostShader, // out
			&pErrorsCompatibility // ppErrorMsgs out
		);

		if (FAILED(hRet)) {
			LOG_TEST_CASE("Couldn't assemble recompiled shader");
			//EmuLog(LOG_LEVEL::WARNING, "Couldn't assemble recompiled shader");
		}
	}

	auto tCompileEnd = std::chrono::high_resolution_clock::now();
	double compileMs = std::chrono::duration<double, std::milli>(tCompileEnd - tCompileStart).count();
	HRESULT compileResult = hRet; // Preserve the actual compile result

	// Determine the log level
	auto hlslErrorLogLevel = FAILED(hRet) ? LOG_LEVEL::ERROR2 : LOG_LEVEL::DEBUG;
	if (pErrors) {
		// Log errors from the initial compilation
		EmuLog(hlslErrorLogLevel, "%s", (char*)(pErrors->GetBufferPointer()));
		pErrors->Release();
		pErrors = nullptr;
	}

	// Failure to recompile in compatibility mode ignored for now
	if (pErrorsCompatibility != nullptr) {
		pErrorsCompatibility->Release();
		pErrorsCompatibility = nullptr;
	}

	LOG_CHECK_ENABLED(LOG_LEVEL::DEBUG) {
		if (g_bPrintfOn) {
			if (!FAILED(compileResult)) {
				// Log disassembly — use a separate HRESULT so we don't clobber compileResult
				HRESULT hDisasm = D3DDisassemble(
					(*ppHostShader)->GetBufferPointer(),
					(*ppHostShader)->GetBufferSize(),
					D3D_DISASM_ENABLE_DEFAULT_VALUE_PRINTS | D3D_DISASM_ENABLE_INSTRUCTION_NUMBERING,
					NULL,
					&pErrors
				);
				if (pErrors) {
					EmuLog(hlslErrorLogLevel, "%s", (char*)(pErrors->GetBufferPointer()));
					pErrors->Release();
				}
			}
		}
	}

	// Save successfully compiled shader to disk cache (async — queued to background thread)
	// Also store in memcache so subsequent lookups never touch disk
	if (!FAILED(compileResult) && *ppHostShader) {
		{
			(*ppHostShader)->AddRef();
			std::lock_guard<std::mutex> lock(g_AsyncMutex);
			g_MemCache[cacheHash] = *ppHostShader;
		}
		if (cacheReady) {
			g_CacheMisses++;
			ShaderCacheLog("SYNC MISS %016llx compile took %.2f ms (profile=%s, blob=%zu bytes) [hits=%d misses=%d]\n",
				cacheHash, compileMs, shader_profile,
				(*ppHostShader)->GetBufferSize(),
				g_CacheHits.load(), g_CacheMisses.load());
			QueueSaveCachedShader(cacheHash, *ppHostShader, effectiveCacheDir);
		}
	} else if (FAILED(compileResult)) {
		ShaderCacheLog("COMPILE FAILED hash=%016llx profile=%s hr=0x%08lX\n",
			cacheHash, shader_profile, compileResult);
	} else if (!cacheReady) {
		ShaderCacheLog("SKIP (cache not ready, g_DataFilePath='%s') hash=%016llx\n",
			g_DataFilePath.c_str(), cacheHash);
	}

	return compileResult;
}

std::ifstream OpenWithRetry(const std::string& path) {
	auto fstream = std::ifstream(path);
	int failures = 0;
	while (fstream.fail()) {
		Sleep(50);
		fstream = std::ifstream(path);

		if (failures++ > 10) {
			// crash?
			CxbxrAbort("Error opening shader file: %s", path);
			break;
		}
	}

	return fstream;
}

int ShaderSources::Update() {
	int versionOnDisk = shaderVersionOnDisk;
	if (shaderVersionLoadedFromDisk != versionOnDisk) {
		bool isInitialLoad = (shaderVersionLoadedFromDisk < 0);
		LoadShadersFromDisk();
		shaderVersionLoadedFromDisk = versionOnDisk;

		// Invalidate disk shader cache only on hot-reload (file watcher fired).
		// Skip on initial cold load: the compiled .cso cache is still valid
		// because EmuCompileShader hashes the full HLSL source + profile,
		// so stale entries from a previous template version are simply never
		// matched and sit harmlessly on disk.
		if (!isInitialLoad) {
			if (!g_ShaderCacheDir.empty() && std::filesystem::exists(g_ShaderCacheDir)) {
				std::error_code ec;
				std::filesystem::remove_all(g_ShaderCacheDir, ec);
				std::filesystem::create_directories(g_ShaderCacheDir, ec);
			}
			// Also invalidate the shared shader cache
			if (!g_SharedShaderCacheDir.empty() && std::filesystem::exists(g_SharedShaderCacheDir)) {
				std::error_code ec;
				std::filesystem::remove_all(g_SharedShaderCacheDir, ec);
				std::filesystem::create_directories(g_SharedShaderCacheDir, ec);
			}
		}
	}

	return shaderVersionLoadedFromDisk;
}

void ShaderSources::LoadShadersFromDisk() {
	const auto hlslDir = std::filesystem::path(szFilePath_CxbxReloaded_Exe)
		.parent_path()
		.append("hlsl");

	// Fixed Function Pixel Shader
	{
		auto dir = hlslDir;
		this->fixedFunctionPixelShaderPath = dir.append("CxbxFixedFunctionPixelShader.hlsl").string();
		std::stringstream tmp;
		tmp << OpenWithRetry(this->fixedFunctionPixelShaderPath).rdbuf();
		this->fixedFunctionPixelShaderHlsl = tmp.str();
	}

	// Vertex Shader Template
	{
		std::stringstream tmp;
		auto dir = hlslDir;
		dir.append("CxbxVertexShaderTemplate.hlsl");
		this->vertexShaderTemplatePath = dir.string();
		tmp << OpenWithRetry(dir.string()).rdbuf();
		std::string hlsl = tmp.str();

		const std::string insertionPoint = "// <XBOX SHADER PROGRAM GOES HERE>\n";
		auto index = hlsl.find(insertionPoint);

		if (index == std::string::npos) {
			// Handle broken shaders
			this->vertexShaderTemplateHlsl[0] = hlsl;
			this->vertexShaderTemplateHlsl[1] = "";
		}
		else
		{
			this->vertexShaderTemplateHlsl[0] = hlsl.substr(0, index);
			this->vertexShaderTemplateHlsl[1] = hlsl.substr(index + insertionPoint.length());
		}
	}

	// Fixed Function Vertex Shader
	{
		auto dir = hlslDir;
		this->fixedFunctionVertexShaderPath = dir.append("CxbxFixedFunctionVertexShader.hlsl").string();
		std::stringstream tmp;
		tmp << OpenWithRetry(this->fixedFunctionVertexShaderPath).rdbuf();
		this->fixedFunctionVertexShaderHlsl = tmp.str();
	}

	// Passthrough Vertex Shader
	{
		auto dir = hlslDir;
		this->vertexShaderPassthroughPath = dir.append("CxbxVertexShaderPassthrough.hlsl").string();
		std::stringstream tmp;
		tmp << OpenWithRetry(this->vertexShaderPassthroughPath).rdbuf();
		this->vertexShaderPassthroughHlsl = tmp.str();
	}

	// Register Combiner Interpreter (PS ubershader)
	{
		auto dir = hlslDir;
		this->registerCombinerInterpreterPath = dir.append("CxbxRegisterCombinerInterpreter.hlsl").string();
		std::stringstream tmp;
		tmp << OpenWithRetry(this->registerCombinerInterpreterPath).rdbuf();
		this->registerCombinerInterpreterHlsl = tmp.str();
	}

	// Vertex Shader Interpreter (VS ubershader)
	{
		auto dir = hlslDir;
		this->vertexShaderInterpreterPath = dir.append("CxbxVertexShaderInterpreter.hlsl").string();
		std::stringstream tmp;
		tmp << OpenWithRetry(this->vertexShaderInterpreterPath).rdbuf();
		this->vertexShaderInterpreterHlsl = tmp.str();
	}
}

void ShaderSources::InitShaderHotloading() {
	static std::jthread fsWatcherThread;

	if (fsWatcherThread.joinable()) {
		EmuLog(LOG_LEVEL::ERROR2, "Ignoring request to start shader file watcher - it has already been started.");
		return;
	}

	EmuLog(LOG_LEVEL::DEBUG, "Starting shader file watcher...");

	fsWatcherThread = std::jthread([]{
		// Determine the filename and directory for the fixed function shader
		char cxbxExePath[MAX_PATH];
		GetModuleFileName(GetModuleHandle(nullptr), cxbxExePath, MAX_PATH);
		auto hlslDir = std::filesystem::path(cxbxExePath).parent_path().append("hlsl/");

		HANDLE changeHandle = FindFirstChangeNotification(hlslDir.string().c_str(), false, FILE_NOTIFY_CHANGE_LAST_WRITE);

		if (changeHandle == INVALID_HANDLE_VALUE) {
			DWORD errorCode = GetLastError();
			EmuLog(LOG_LEVEL::ERROR2, "Error initializing shader file watcher: %d", errorCode);

			return 1;
		}

		while (true) {
			if (FindNextChangeNotification(changeHandle)) {
				WaitForSingleObject(changeHandle, INFINITE);

				// Wait for changes to stop..
				// Will usually be at least two - one for the file and one for the directory
				while (true) {
					FindNextChangeNotification(changeHandle);
					if (WaitForSingleObject(changeHandle, 100) == WAIT_TIMEOUT) {
						break;
					}
				}

				EmuLog(LOG_LEVEL::DEBUG, "Change detected in shader folder");

				g_ShaderSources.shaderVersionOnDisk++;
			}
			else {
				EmuLog(LOG_LEVEL::ERROR2, "Shader filewatcher failed to get the next notification");
				break;
			}
		}

		EmuLog(LOG_LEVEL::DEBUG, "Shader file watcher exiting...");

		// until there is a way to disable hotloading
		// this is always an error
		FindCloseChangeNotification(changeHandle);
		return 1;
	});
}
