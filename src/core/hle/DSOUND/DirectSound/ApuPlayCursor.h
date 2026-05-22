#pragma once

#include <windows.h>
#include <cstdint>

// APU play cursor emulation for games that poll voice descriptor CBO directly.
// On real Xbox hardware, the APU Voice Processor continuously advances the
// Current Buffer Offset (CBO) in voice descriptors as audio data is consumed.
// HLE DirectSound replaces the entire audio path, so no voice descriptors exist
// and CBO is never updated — games that bypass GetCurrentPosition and read CBO
// directly will spin forever (or block on an event that's never signaled by APU
// interrupts).  This structure tracks a detected play cursor and advances it
// from dsound_worker at a rate matching typical Xbox audio output.
struct ApuPlayCursorState {
    volatile DWORD* pCursor;    // Xbox virtual address of the play cursor (contiguous memory)
    DWORD           bufSize;    // Total buffer size (cursor wraps or stops at this value)
    int64_t         lastQPC;    // QPC timestamp of last advancement
    DWORD           rate;       // Advancement rate in bytes per second
    void*           pEvent;     // Xbox KEVENT to signal when cursor advances (may be nullptr)
};

extern ApuPlayCursorState g_ApuPlayCursor;
