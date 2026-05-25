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

#include <windows.h>
#include <tlhelp32.h>
#include <array>
#include "Timer.h"
#include "common\util\CxbxUtil.h"
#include "core\kernel\support\EmuFS.h"
#include "core\kernel\exports\EmuKrnlPs.hpp"
#include "core\kernel\exports\EmuKrnl.h"
#include "devices\Xbox.h"
#include "devices\audio\AC97Device.h"
#include "devices\audio\APUDevice.h"
#include "devices\usb\OHCI.h"
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
// All subsystems (NV2A, OHCI, PIT) operate in QPC ticks
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

static void dispatch_non_periodic_events()
{
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
	// AC97 ServiceAudio is called from apu_tick (if present) or directly
	if (g_AC97 != nullptr) {
		g_AC97->ServiceAudio();
	}

	std::array<uint64_t, 5> deadlines = {
		pit_tick(now),
		g_NV2A->vblank_tick(now),
		g_NV2A->ptimer_tick(now),
		g_USB0->m_HostController->OHCI_tick(now),
		g_APU != nullptr ? g_APU->apu_tick(now) : UINT64_MAX
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
