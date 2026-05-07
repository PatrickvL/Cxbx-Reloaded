// Backend_D3D11_Profiler.h — Lightweight whole-system CPU timing profiler.
// Accumulates QPC ticks across named phases, dumps breakdown once/sec via OutputDebugString.
// Covers: MMIO/VEH, PFIFO, page tracking, render state, draw calls, present.
// Zero overhead when disabled (g_bCxbxProfilerEnabled = false).
//
// Per-frame statistics: tracks min/max/avg/stddev of frame intervals and per-phase
// times to identify stutter-causing outliers.
//
// Thread safety: render-thread phases are single-threaded (no atomics needed).
// Cross-thread counters (MMIO, PFIFO) use InterlockedAdd64 for lock-free accumulation.
#pragma once

#include <windows.h>
#include <cstdio>
#include <cmath>
#include <algorithm>

// ============================================================
// Phase IDs — add new phases here and in the name table below.
// ============================================================
enum CxbxProfilePhase {
    // --- Render thread phases (single-threaded, no atomics) ---
    PROF_PFIFO_FLUSH = 0,   // Inline pfifo_flush_to_pgraph on render thread
    PROF_PGRAPH_LOCK_WAIT,  // Time waiting to acquire pgraph_lock
    PROF_VS_SHADER,         // Vertex shader selection (JIT attempt + fallback)
    PROF_VS_JIT_COMPILE,    // VS JIT D3DCompile (cache miss only)
    PROF_VS_CONSTANTS,      // Upload VS constants (cbuffer map)
    PROF_TEXTURES,          // Texture lookup + SRV creation + scaling
    PROF_PIPELINE_STATE,    // Blend/depth/rasterizer state objects
    PROF_SAMPLERS,          // Sampler state creation/cache
    PROF_PIXEL_SHADER,      // Pixel shader update (JIT attempt + fallback)
    PROF_PS_JIT_COMPILE,    // PS JIT D3DCompile (cache miss only)
    PROF_RENDER_TARGET,     // RT resolve from PGRAPH surface offsets
    PROF_VIEWPORT,          // Viewport/scissor from PGRAPH
    PROF_PAGE_FLUSH,        // GetWriteWatch + memcpy dirty pages to GPU mirror
    PROF_DRAW_CALL,         // ID3D11DeviceContext::Draw
    PROF_PRESENT_BLIT,      // CxbxBltSurface (backbuffer blit)
    PROF_PRESENT_OVERLAY,   // YUY2→ARGB overlay conversion
    PROF_PRESENT_SWAP,      // IDXGISwapChain::Present call
    PROF_GPU_FLUSH,         // ID3D11DeviceContext::Flush before Present
    PROF_BLOCKONTTIME,      // D3D_BlockOnTime → pfifo_flush_to_pgraph

    // --- Cross-thread phases (use atomic accumulators) ---
    PROF_MMIO,              // VEH MMIO handler (page fault + decode + dispatch)
    PROF_PFIFO_PUSHER,      // Pusher thread: DMA parse + method dispatch
    PROF_VBLANK_JITTER,     // VBlank scheduling jitter (deviation from ideal period)
    PROF_SYSEVENTS_LOOP,    // system_events loop iteration time

    PROF_PHASE_COUNT
};

// First cross-thread phase index (everything >= this uses atomics)
#define PROF_CROSS_THREAD_START PROF_MMIO

inline const char* g_ProfilePhaseNames[PROF_PHASE_COUNT] = {
    "PfifoFlush",
    "PgraphLockWait",
    "VSShader",
    "VSJitCompile",
    "VSConstants",
    "Textures",
    "PipelineState",
    "Samplers",
    "PixelShader",
    "PSJitCompile",
    "RenderTarget",
    "Viewport",
    "PageFlush",
    "DrawCall",
    "PresentBlit",
    "PresentOverlay",
    "PresentSwap",
    "GpuFlush",
    "BlockOnTime",
    "MMIO",
    "PfifoPusher",
    "VBlankJitter",
    "SysEventsLoop"
};

// ============================================================
// Per-frame statistics accumulator (for frame interval + key phases)
// ============================================================
struct CxbxFrameStats {
    double sum;
    double sumSq;
    double minVal;
    double maxVal;
    UINT   count;

    void reset() {
        sum = sumSq = 0.0;
        minVal = 1e30;
        maxVal = 0.0;
        count = 0;
    }

    void record(double ms) {
        sum += ms;
        sumSq += ms * ms;
        if (ms < minVal) minVal = ms;
        if (ms > maxVal) maxVal = ms;
        count++;
    }

    double avg() const { return count > 0 ? sum / count : 0.0; }
    double stddev() const {
        if (count < 2) return 0.0;
        double mean = avg();
        return sqrt(sumSq / count - mean * mean);
    }
};

// ============================================================
// Global state
// ============================================================
inline bool     g_bCxbxProfilerEnabled = true;
inline UINT     g_ProfileDrawCount = 0;
inline UINT     g_ProfileFrameCount = 0;
inline LONGLONG g_ProfileAccum[PROF_PHASE_COUNT] = {};

// Per-frame distribution tracking
inline CxbxFrameStats g_ProfileFrameIntervalStats;     // Frame-to-frame interval
inline CxbxFrameStats g_ProfilePhaseFrameStats[PROF_PHASE_COUNT]; // Per-phase per-frame stats
inline LONGLONG       g_ProfileFramePhaseAccum[PROF_PHASE_COUNT] = {}; // Current frame's per-phase ticks
inline LARGE_INTEGER  g_ProfileLastFrameQPC = {};       // QPC at last frame tick

// Cross-thread counters (hit counts for high-frequency events)
inline volatile LONG g_ProfileMMIOCount = 0;     // Number of MMIO accesses
inline volatile LONG g_ProfilePusherMethods = 0; // Number of PFIFO methods dispatched
inline volatile LONG g_ProfileCSDispatchCount = 0; // Number of CS dispatches

// JIT vs interpreter hit counters
inline volatile LONG g_ProfileVSJITHits = 0;       // VS JIT cache hits
inline volatile LONG g_ProfileVSJITCompiles = 0;   // VS JIT cache misses (new compiles)
inline volatile LONG g_ProfileVSInterpreterHits = 0; // VS interpreter fallback
inline volatile LONG g_ProfilePSJITHits = 0;       // PS JIT cache hits
inline volatile LONG g_ProfilePSJITCompiles = 0;   // PS JIT cache misses (new compiles)
inline volatile LONG g_ProfilePSInterpreterHits = 0; // PS interpreter fallback

inline LARGE_INTEGER g_ProfileLastReport = {};

// ============================================================
// RAII scoped timer — accumulates into the given phase bucket.
// Uses InterlockedAdd for cross-thread phases, plain add for render-thread phases.
// Also records per-frame stats for render-thread phases.
// ============================================================
struct CxbxProfileScope {
    CxbxProfilePhase phase;
    LARGE_INTEGER start;

    __forceinline CxbxProfileScope(CxbxProfilePhase p) : phase(p) {
        if (g_bCxbxProfilerEnabled) QueryPerformanceCounter(&start);
        else start.QuadPart = 0;
    }

    __forceinline ~CxbxProfileScope() {
        if (!g_bCxbxProfilerEnabled || start.QuadPart == 0) return;
        LARGE_INTEGER end;
        QueryPerformanceCounter(&end);
        LONGLONG delta = end.QuadPart - start.QuadPart;
        if (phase >= PROF_CROSS_THREAD_START) {
            InterlockedAdd64(&g_ProfileAccum[phase], delta);
        } else {
            g_ProfileAccum[phase] += delta;
            g_ProfileFramePhaseAccum[phase] += delta;
        }
    }
};

// Convenience macro
#define CXBX_PROFILE_SCOPE(phase) CxbxProfileScope _prof_##phase(phase)

// ============================================================
// Call once per frame (in present path) to check if 1s elapsed
// and dump the breakdown with per-frame statistics.
// ============================================================
inline void CxbxProfilerFrameTick()
{
    if (!g_bCxbxProfilerEnabled) return;

    LARGE_INTEGER now, freq;
    QueryPerformanceCounter(&now);
    QueryPerformanceFrequency(&freq);

    // Record frame interval
    if (g_ProfileLastFrameQPC.QuadPart != 0) {
        double frameMs = (double)(now.QuadPart - g_ProfileLastFrameQPC.QuadPart) * 1000.0 / (double)freq.QuadPart;
        g_ProfileFrameIntervalStats.record(frameMs);
    }
    g_ProfileLastFrameQPC = now;

    // Record per-phase stats for this frame (render-thread phases only)
    for (int i = 0; i < PROF_CROSS_THREAD_START; i++) {
        if (g_ProfileFramePhaseAccum[i] > 0) {
            double ms = (double)g_ProfileFramePhaseAccum[i] * 1000.0 / (double)freq.QuadPart;
            g_ProfilePhaseFrameStats[i].record(ms);
        }
        g_ProfileFramePhaseAccum[i] = 0;
    }

    g_ProfileFrameCount++;

    if (g_ProfileLastReport.QuadPart == 0) {
        g_ProfileLastReport = now;
        g_ProfileFrameIntervalStats.reset();
        for (int i = 0; i < PROF_PHASE_COUNT; i++)
            g_ProfilePhaseFrameStats[i].reset();
        return;
    }

    double elapsed = (double)(now.QuadPart - g_ProfileLastReport.QuadPart) / (double)freq.QuadPart;
    if (elapsed < 1.0) return;

    // ---- Dump aggregate breakdown ----
    char buf[8192];
    int pos = sprintf_s(buf, "[CXBX-PROF] %u frames, %u draws, %ld MMIOs, %ld methods, %ld CS | VS(jit=%ld comp=%ld interp=%ld) PS(jit=%ld comp=%ld interp=%ld)\n",
        g_ProfileFrameCount, g_ProfileDrawCount,
        (long)g_ProfileMMIOCount, (long)g_ProfilePusherMethods, (long)g_ProfileCSDispatchCount,
        (long)g_ProfileVSJITHits, (long)g_ProfileVSJITCompiles, (long)g_ProfileVSInterpreterHits,
        (long)g_ProfilePSJITHits, (long)g_ProfilePSJITCompiles, (long)g_ProfilePSInterpreterHits);

    // Frame interval stats
    {
        auto& s = g_ProfileFrameIntervalStats;
        if (s.count > 0) {
            pos += sprintf_s(buf + pos, sizeof(buf) - pos,
                "  FrameInterval: avg=%.2fms min=%.2fms max=%.2fms stddev=%.2fms (%u frames)\n",
                s.avg(), s.minVal, s.maxVal, s.stddev(), s.count);
        }
    }

    // Per-phase breakdown with per-frame stats
    pos += sprintf_s(buf + pos, sizeof(buf) - pos,
        "  %-16s %8s %8s %8s %8s %8s %8s\n",
        "Phase", "Total", "Avg/F", "Min/F", "Max/F", "StdDev", "Count");

    for (int i = 0; i < PROF_PHASE_COUNT; i++) {
        LONGLONG ticks = (i >= PROF_CROSS_THREAD_START)
            ? InterlockedExchange64(&g_ProfileAccum[i], 0)
            : g_ProfileAccum[i];
        if (ticks > 0) {
            double totalMs = (double)ticks * 1000.0 / (double)freq.QuadPart;
            auto& s = g_ProfilePhaseFrameStats[i];
            if (s.count > 0) {
                pos += sprintf_s(buf + pos, sizeof(buf) - pos,
                    "  %-16s %7.2fms %7.2fms %7.2fms %7.2fms %7.2fms %5u\n",
                    g_ProfilePhaseNames[i], totalMs,
                    s.avg(), s.minVal, s.maxVal, s.stddev(), s.count);
            } else {
                pos += sprintf_s(buf + pos, sizeof(buf) - pos,
                    "  %-16s %7.2fms\n",
                    g_ProfilePhaseNames[i], totalMs);
            }
        }
    }

    pos += sprintf_s(buf + pos, sizeof(buf) - pos, "\n");
    printf("%s", buf);
    fflush(stdout);
    OutputDebugStringA(buf);

    // Write to a log file next to the executable (survives stdout redirection)
    {
        static HANDLE s_hLog = INVALID_HANDLE_VALUE;
        if (s_hLog == INVALID_HANDLE_VALUE) {
            char exePath[MAX_PATH] = {};
            GetModuleFileNameA(nullptr, exePath, MAX_PATH);
            char* lastSlash = strrchr(exePath, '\\');
            if (lastSlash) *(lastSlash + 1) = '\0';
            strcat_s(exePath, "CxbxProfiler.log");
            s_hLog = CreateFileA(exePath,
                GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS,
                FILE_ATTRIBUTE_NORMAL, NULL);
        }
        if (s_hLog != INVALID_HANDLE_VALUE) {
            DWORD written;
            WriteFile(s_hLog, buf, (DWORD)pos, &written, NULL);
            FlushFileBuffers(s_hLog);
        }
    }

    // Reset render-thread accumulators (cross-thread ones reset via InterlockedExchange above)
    for (int i = 0; i < PROF_CROSS_THREAD_START; i++)
        g_ProfileAccum[i] = 0;
    g_ProfileFrameCount = 0;
    g_ProfileDrawCount = 0;
    InterlockedExchange(&g_ProfileMMIOCount, 0);
    InterlockedExchange(&g_ProfilePusherMethods, 0);
    InterlockedExchange(&g_ProfileCSDispatchCount, 0);
    InterlockedExchange(&g_ProfileVSJITHits, 0);
    InterlockedExchange(&g_ProfileVSJITCompiles, 0);
    InterlockedExchange(&g_ProfileVSInterpreterHits, 0);
    InterlockedExchange(&g_ProfilePSJITHits, 0);
    InterlockedExchange(&g_ProfilePSJITCompiles, 0);
    InterlockedExchange(&g_ProfilePSInterpreterHits, 0);
    g_ProfileLastReport = now;

    // Reset per-frame stats
    g_ProfileFrameIntervalStats.reset();
    for (int i = 0; i < PROF_PHASE_COUNT; i++)
        g_ProfilePhaseFrameStats[i].reset();
}
