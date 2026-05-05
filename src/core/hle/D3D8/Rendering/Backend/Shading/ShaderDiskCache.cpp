// ShaderDiskCache.cpp — Persistent disk cache for JIT-compiled shader bytecode

#define LOG_PREFIX CXBXR_MODULE::VTXSH

#include "ShaderDiskCache.h"
#include <string>
#include "common/FilePaths.hpp"
#include "core/kernel/init/CxbxKrnl.h"
#include "core/kernel/support/Emu.h"
#include "common/xbe/Xbe.h"

#include <filesystem>
#include <mutex>
#include <queue>
#include <thread>
#include <atomic>
#include <vector>
#include <cstdio>
#include <cstring>

namespace ShaderDiskCache {

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
static std::string g_CacheDir;
static std::mutex g_InitMutex;

static std::mutex g_SaveQueueMutex;
static std::queue<std::pair<std::string, std::vector<uint8_t>>> g_SaveQueue;
static std::thread g_SaveThread;
static std::atomic<bool> g_SaveThreadRunning{false};

// ---------------------------------------------------------------------------
// Background save worker
// ---------------------------------------------------------------------------
static void SaveWorker()
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
            std::this_thread::yield();
            continue;
        }

        FILE* fp = fopen(item.first.c_str(), "wb");
        if (fp) {
            fwrite(item.second.data(), 1, item.second.size(), fp);
            fclose(fp);
        }
    }
}

// ---------------------------------------------------------------------------
// Ensure cache directory exists (lazy init)
// ---------------------------------------------------------------------------
static bool EnsureCacheDir()
{
    if (!g_CacheDir.empty()) return true;

    std::lock_guard<std::mutex> lock(g_InitMutex);
    if (!g_CacheDir.empty()) return true; // double-check

    if (g_DataFilePath.empty()) return false;
    if (!g_pCertificate) return false;

    // Build per-game directory name: <TitleID>-<GameName>
    char titleIdStr[16];
    snprintf(titleIdStr, sizeof(titleIdStr), "%08X", g_pCertificate->dwTitleId);

    std::string gameName;
    if (CxbxKrnl_Xbe && CxbxKrnl_Xbe->m_szAsciiTitle[0]) {
        gameName = CxbxKrnl_Xbe->m_szAsciiTitle;
        // Sanitize for filesystem
        for (char& c : gameName) {
            if (c == '\\' || c == '/' || c == ':' || c == '*' ||
                c == '?' || c == '"' || c == '<' || c == '>' || c == '|')
                c = '_';
        }
        // Trim trailing spaces/dots
        while (!gameName.empty() && (gameName.back() == ' ' || gameName.back() == '.'))
            gameName.pop_back();
    }

    std::string dirName = std::string(titleIdStr);
    if (!gameName.empty()) dirName += "-" + gameName;

    g_CacheDir = g_DataFilePath + "\\ShaderCache\\" + dirName;

    std::error_code ec;
    std::filesystem::create_directories(g_CacheDir, ec);
    if (ec) {
        g_CacheDir.clear();
        return false;
    }

    // Start background save thread
    if (!g_SaveThreadRunning) {
        g_SaveThreadRunning = true;
        g_SaveThread = std::thread(SaveWorker);
        g_SaveThread.detach();
    }

    EmuLog(LOG_LEVEL::INFO, "ShaderDiskCache: %s", g_CacheDir.c_str());
    return true;
}

static std::string GetCachePath(uint64_t hash)
{
    char filename[32];
    snprintf(filename, sizeof(filename), "%016llx.cso", hash);
    return g_CacheDir + "\\" + filename;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------
ID3DBlob* TryLoad(uint64_t hash)
{
    if (!EnsureCacheDir()) return nullptr;

    std::string path = GetCachePath(hash);
    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) return nullptr;

    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    if (size < 16) { // minimum valid DXBC size
        fclose(fp);
        return nullptr;
    }
    fseek(fp, 0, SEEK_SET);

    ID3DBlob* pBlob = nullptr;
    HRESULT hr = D3DCreateBlob(size, &pBlob);
    if (FAILED(hr)) {
        fclose(fp);
        return nullptr;
    }

    size_t readBytes = fread(pBlob->GetBufferPointer(), 1, size, fp);
    fclose(fp);

    if ((long)readBytes != size) {
        pBlob->Release();
        return nullptr;
    }

    return pBlob;
}

void Save(uint64_t hash, ID3DBlob* pBlob)
{
    if (!pBlob) return;
    if (!EnsureCacheDir()) return;

    std::string path = GetCachePath(hash);

    // Don't overwrite existing cache files
    if (std::filesystem::exists(path)) return;

    // Copy bytecode into save queue
    std::vector<uint8_t> data(pBlob->GetBufferSize());
    memcpy(data.data(), pBlob->GetBufferPointer(), pBlob->GetBufferSize());

    {
        std::lock_guard<std::mutex> lock(g_SaveQueueMutex);
        g_SaveQueue.emplace(std::move(path), std::move(data));
    }
}

void Shutdown()
{
    g_SaveThreadRunning = false;

    // Drain remaining items on this thread
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
        }
    }
}

} // namespace ShaderDiskCache
