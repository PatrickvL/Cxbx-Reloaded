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
// *  (c) 2018      ergo720
// *
// *  All rights reserved
// *
// ******************************************************************

#include <core\kernel\exports\xboxkrnl.h>
#include "devices\Xbox.h"

#include <windows.h>
#include <tlhelp32.h>
#include <array>
#include "Timer.h"
#include "common\util\CxbxUtil.h"
#include "core\kernel\support\EmuFS.h"
#include "core\kernel\exports\EmuKrnlPs.hpp"
#include "core\kernel\exports\EmuKrnl.h"
#include "devices\Xbox.h"
#include "devices\usb\OHCI.h"
#include "core\hle\DSOUND\DirectSound\DirectSoundGlobal.hpp"
#include "core\hle\D3D8\Rendering\Backend\Backend_D3D11_Profiler.h"


std::atomic_uint64_t HostLastQPC; // last absolute host QPC reading
static uint64_t pit_last; // QPC ticks (relative to start) when PIT last fired
static int64_t PIT_PERIOD_QPC; // 1ms in QPC ticks, set by timer_init()
// The frequency of the high resolution clock of the host, and the start time
int64_t HostQPCFrequency, HostQPCStartTime;

// High-resolution waitable timer handle, created once in timer_init().
// CREATE_WAITABLE_TIMER_HIGH_RESOLUTION (Win10 1803+) uses a dedicated
// kernel timer queue with ~0.5ms resolution and zero CPU burn.
static HANDLE g_hPreciseTimer = NULL;


void timer_init()
{
	QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER *>(&HostQPCFrequency));
	QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER *>(&HostQPCStartTime));
	HostLastQPC = HostQPCStartTime;
	PIT_PERIOD_QPC = HostQPCFrequency / 1000; // 1ms in QPC ticks
	pit_last = 0; // get_now() returns 0 at start

	// Create high-resolution waitable timer (Win10 1803+, build 17134).
	g_hPreciseTimer = CreateWaitableTimerEx(NULL, NULL,
		CREATE_WAITABLE_TIMER_HIGH_RESOLUTION, TIMER_ALL_ACCESS);
	// If CREATE_WAITABLE_TIMER_HIGH_RESOLUTION fails, fall back to a
	// regular waitable timer — with timeBeginPeriod(1) already called
	// in CxbxrKrnlInit, even regular waitable timers get ~1ms precision,
	// sufficient for our needs.
	if (g_hPreciseTimer == NULL) {
		g_hPreciseTimer = CreateWaitableTimerEx(NULL, NULL, 0, TIMER_ALL_ACCESS);
	}

	// Synchronize xbox system time with host time
	LARGE_INTEGER HostSystemTime;
	GetSystemTimeAsFileTime((LPFILETIME)&HostSystemTime);
	xbox::KeSystemTime.High2Time = HostSystemTime.u.HighPart;
	xbox::KeSystemTime.LowPart = HostSystemTime.u.LowPart;
	xbox::KeSystemTime.High1Time = HostSystemTime.u.HighPart;
}

// Precise sleep until an absolute QPC deadline, with zero busy-wait.
// Returns the final QPC value at wake-up, so callers can use it as
// the next anchor without a redundant QueryPerformanceCounter call.
int64_t SleepPrecise(int64_t targetQPC)
{
	LARGE_INTEGER now;
	QueryPerformanceCounter(&now);

	// Early-out: target already passed
	if (now.QuadPart >= targetQPC)
		return now.QuadPart;

	if (g_hPreciseTimer != NULL) {
		// Convert QPC delta to 100ns units (negative = relative deadline).
		// Multiply-before-divide preserves precision for short intervals.
		LARGE_INTEGER dueTime;
		dueTime.QuadPart = -((targetQPC - now.QuadPart) * 10000000LL / HostQPCFrequency);

		// Clamp to at least -1 (100ns) to avoid zero (which means "already signaled").
		if (dueTime.QuadPart == 0)
			dueTime.QuadPart = -1;

		// High-resolution waitable timer: kernel wakes us at ~0.5ms precision
		// with no CPU burn. The slight undershoot (<0.5ms) is acceptable since
		// callers re-read QPC after wake-up.
		SetWaitableTimerEx(g_hPreciseTimer, &dueTime, 0, NULL, NULL, NULL, 0);
		WaitForSingleObject(g_hPreciseTimer, INFINITE);
	} else {
		// Fallback for older Windows: convert to ms and use Sleep().
		// timeBeginPeriod(1) is already called in CxbxrKrnlInit, so
		// Sleep(1) actually sleeps ~1ms. Sleep(0) yields timeslice.
		Sleep((DWORD)((targetQPC - now.QuadPart) * 1000 / HostQPCFrequency));
	}

	QueryPerformanceCounter(&now);
	return now.QuadPart;
}

// ── Emulated Xbox clock ──────────────────────────────────────────

// Read the host QPC and return elapsed ticks since timer_init().
// All subsystems (NV2A, OHCI, DSound, PIT) operate in QPC ticks
// to avoid integer truncation from µs conversion.
uint64_t get_now()
{
	LARGE_INTEGER now;
	QueryPerformanceCounter(&now);
	HostLastQPC = now.QuadPart;
	return now.QuadPart - HostQPCStartTime;
}

// ── PIT (Programmable Interval Timer) ────────────────────────────

// Dispatch the PIT clock ISR if overdue, return absolute next deadline.
static uint64_t pit_tick(uint64_t now)
{
	uint64_t next = pit_last + PIT_PERIOD_QPC;
	if (now >= next) {
		uint64_t elapsed_us = (now - pit_last) * 1000000 / HostQPCFrequency;
		xbox::KiClockIsr(elapsed_us);
		pit_last = now;
		return now + PIT_PERIOD_QPC;
	}
	return next;
}

// ── Non-periodic event dispatch ──────────────────────────────────

// Detect game threads stuck polling an APU play cursor that HLE DirectSound
// doesn't advance.  Scans all threads looking for a tight spin loop reading
// from contiguous memory (0x80000000+).  When found, initializes
// g_ApuPlayCursor so dsound_worker can advance it continuously.
//
// Safety: requires the thread to be at the SAME EIP reading the SAME address
// with an UNCHANGED value across two consecutive 1-second scans.  This prevents
// false positives from threads that briefly pass through spin-waits during
// normal operation (movie playback, vsync waits, etc.).
static void detect_apu_play_cursor_spin()
{
	// Already detected — nothing more to do.
	if (g_ApuPlayCursor.pCursor != nullptr)
		return;

	static DWORD s_lastCheck = 0;
	DWORD now = GetTickCount();
	if (now - s_lastCheck < 1000) return;
	s_lastCheck = now;

	// State from previous scan for two-scan confirmation.
	static DWORD s_candidateTid = 0;
	static DWORD s_candidateEip = 0;
	static DWORD s_candidateEdx = 0;   // cursor pointer address
	static DWORD s_candidateEsi = 0;   // buffer size (from ESI register)
	static DWORD s_candidateVal = 0;   // *cursor value at last scan

	DWORD myTid = GetCurrentThreadId();
	DWORD pid = GetCurrentProcessId();
	HANDLE snap = CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD, 0);
	if (snap == INVALID_HANDLE_VALUE) return;

	bool found = false;
	THREADENTRY32 te;
	te.dwSize = sizeof(te);
	if (Thread32First(snap, &te)) {
		do {
			if (te.th32OwnerProcessID != pid) continue;
			if (te.th32ThreadID == myTid) continue;

			HANDLE hThread = OpenThread(THREAD_SUSPEND_RESUME | THREAD_GET_CONTEXT, FALSE, te.th32ThreadID);
			if (!hThread) continue;

			SuspendThread(hThread);
			CONTEXT ctx = {};
			ctx.ContextFlags = CONTEXT_CONTROL | CONTEXT_INTEGER;
			if (GetThreadContext(hThread, &ctx)) {
				// Thread must be in game code range with EDX pointing to contiguous memory
				if (ctx.Eip >= 0x00100000 && ctx.Eip < 0x00400000 &&
				    ctx.Edx >= 0x80000000 && ctx.Edx < 0x84000000) {
					DWORD curVal = *reinterpret_cast<volatile DWORD*>(ctx.Edx);

					// Check if this matches our previous candidate: same thread, same EIP,
					// same cursor address, and cursor value UNCHANGED = genuinely stuck.
					if (te.th32ThreadID == s_candidateTid &&
					    ctx.Eip == s_candidateEip &&
					    ctx.Edx == s_candidateEdx &&
					    curVal == s_candidateVal) {
						// Confirmed stuck — activate cursor advancement.
						LARGE_INTEGER qpc;
						QueryPerformanceCounter(&qpc);
						g_ApuPlayCursor.pCursor = reinterpret_cast<volatile DWORD*>(ctx.Edx);
						g_ApuPlayCursor.bufSize = s_candidateEsi;
						// 48 kHz * 2 bytes (16-bit mono) = 96,000 bytes/sec
						g_ApuPlayCursor.rate = 96000;
						g_ApuPlayCursor.lastQPC = qpc.QuadPart;
						EmuLogEx(CXBXR_MODULE::DSOUND, LOG_LEVEL::INFO,
							"APU play cursor detected: addr=0x%08X bufSize=%u val=%u",
							ctx.Edx, s_candidateEsi, curVal);
						found = true;
					} else {
						// Record as candidate for next scan's confirmation.
						s_candidateTid = te.th32ThreadID;
						s_candidateEip = ctx.Eip;
						s_candidateEdx = ctx.Edx;
						s_candidateEsi = ctx.Esi;
						s_candidateVal = curVal;
					}
				}
			}
			ResumeThread(hThread);
			CloseHandle(hThread);
			if (found) break;
		} while (Thread32Next(snap, &te));
	}
	CloseHandle(snap);
}

// Proactive memory scan for APU play cursor pointer.
// Scans the XBE .data section for DWORDs pointing into Xbox contiguous memory
// (0x80000000-0x84000000) adjacent to a plausible buffer-size DWORD.
// Unlike detect_apu_play_cursor_spin(), this does NOT suspend threads — it
// reads game memory directly from the system_events thread, safe from any
// deadlock or timing interference.  Uses two-scan confirmation (2s apart) to
// prevent false positives.
static void detect_apu_play_cursor_scan()
{
	if (g_ApuPlayCursor.pCursor != nullptr)
		return;

	static DWORD s_lastCheck = 0;
	DWORD now = GetTickCount();
	if (now - s_lastCheck < 2000) return;
	s_lastCheck = now;

	static uintptr_t s_candidateAddr = 0;
	static DWORD     s_candidateCbo = 0;

	if (s_candidateAddr != 0) {
		DWORD* p = reinterpret_cast<DWORD*>(s_candidateAddr);
		DWORD val = *p;
		if (val >= 0x80000000 && val < 0x84000000 && (val & 0x3) == 0) {
			DWORD cbo = *reinterpret_cast<volatile DWORD*>(val);
			if (cbo == s_candidateCbo) {
				LARGE_INTEGER qpc;
				QueryPerformanceCounter(&qpc);
				g_ApuPlayCursor.pCursor = reinterpret_cast<volatile DWORD*>(val);
				g_ApuPlayCursor.bufSize = 0x200000;
				g_ApuPlayCursor.rate = 96000;
				g_ApuPlayCursor.lastQPC = qpc.QuadPart;
				g_ApuPlayCursor.pEvent = nullptr;
				EmuLogEx(CXBXR_MODULE::DSOUND, LOG_LEVEL::INFO,
					"APU play cursor detected via memory scan: addr=0x%08X cbo=%u",
					val, cbo);
				return;
			}
		}
		s_candidateAddr = 0;
	}

	DWORD* const SCAN_START = reinterpret_cast<DWORD*>(0x1A0000);
	for (DWORD* p = SCAN_START; p < reinterpret_cast<DWORD*>(0x300000); p++) {
		DWORD val = *p;
		if (val >= 0x80000000 && val < 0x84000000 && (val & 0x3) == 0) {
			DWORD cbo = *reinterpret_cast<volatile DWORD*>(val);
			if (cbo < 0x100000) {
				s_candidateAddr = reinterpret_cast<uintptr_t>(p);
				s_candidateCbo = cbo;
				return;
			}
		}
	}
}

static void dispatch_non_periodic_events()
{
	dsound_worker();

	// Advance APU voice CBO for games that read voice descriptors directly
	// from contiguous memory (bypassing HLE DirectSound GetCurrentPosition).
	// This is the hardware-accurate equivalent of the real APU advancing CBO
	// as audio samples are consumed by the DMA engine.
	if (g_APU) {
		g_APU->AdvanceVoiceCursors();
	}

	for (int i = 0; i < MAX_BUS_INTERRUPT_LEVEL; i++) {
		// Skip IRQ 3 (GPU/NV2A) — delivered explicitly by
		// nv2a_vblank_interrupt and PGRAPH INTR_ERROR mechanism.
		if (i == 3) continue;

		if (g_bEnableAllInterrupts && HalSystemInterrupts[i].IsPending() && EmuInterruptList[i] && EmuInterruptList[i]->Connected) {
			HalSystemInterrupts[i].Trigger(EmuInterruptList[i]);
		}
	}
}

// ── Periodic event dispatch + deadline ────────────────────────────

// Tick all periodic subsystems — each dispatches if overdue and
// returns the absolute QPC time of its next deadline. Returns
// the earliest deadline (relative to HostQPCStartTime).
static uint64_t dispatch_periodic_events(uint64_t now)
{
	std::array<uint64_t, 5> deadlines = {
		pit_tick(now),
		g_NV2A->vblank_tick(now),
		g_NV2A->ptimer_tick(now),
		g_USB0->m_HostController->OHCI_tick(now),
		dsound_tick(now)
	};
	return *std::min_element(deadlines.begin(), deadlines.end());
}

// ── System events thread ─────────────────────────────────────────

xbox::void_xt NTAPI system_events(xbox::PVOID arg)
{
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);

	// Run at DPC level to prevent this thread from executing APCs/DPCs
	xbox::KeRaiseIrqlToDpcLevel();

	while (true) {
		LARGE_INTEGER loop_start;
		if (g_bCxbxProfilerEnabled) QueryPerformanceCounter(&loop_start);

		// 1. Read the host clock (QPC ticks since start)
		const uint64_t now = get_now();

		// 2. Dispatch all events and find earliest next deadline
		dispatch_non_periodic_events();
		const uint64_t next_deadline = dispatch_periodic_events(now);

		// 3. Sleep until the absolute deadline (skip if no subsystem is active)
		if (next_deadline != UINT64_MAX) {
			SleepPrecise((int64_t)next_deadline + HostQPCStartTime);
		}

		if (g_bCxbxProfilerEnabled) {
			LARGE_INTEGER loop_end;
			QueryPerformanceCounter(&loop_end);
			InterlockedAdd64(&g_ProfileAccum[PROF_SYSEVENTS_LOOP], loop_end.QuadPart - loop_start.QuadPart);
		}
	}
}

int64_t Timer_GetScaledPerformanceCounter(int64_t Period)
{
	LARGE_INTEGER currentQPC;
	QueryPerformanceCounter(&currentQPC);

	// Scale frequency with overflow avoidance, like in std::chrono
	// https://github.com/microsoft/STL/blob/6d2f8b0ed88ea6cba26cc2151f47f678442c1663/stl/inc/chrono#L703
	const int64_t currentTime = currentQPC.QuadPart - HostQPCStartTime;
	const int64_t whole = (currentTime / HostQPCFrequency) * Period;
	const int64_t part  = (currentTime % HostQPCFrequency) * Period / HostQPCFrequency;

	return whole + part;
}

