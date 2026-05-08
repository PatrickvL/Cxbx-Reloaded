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
#include <chrono>
#include <thread>
#include <vector>
#include <mutex>
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


static std::atomic_uint64_t last_qpc; // last time when QPC was called
static std::atomic_uint64_t exec_time; // total execution time in us since the emulation started
static uint64_t pit_last; // last time when the pit time was updated
static uint64_t pit_last_qpc; // last QPC time of the pit
// The frequency of the high resolution clock of the host, and the start time
int64_t HostQPCFrequency, HostQPCStartTime;


void timer_init()
{
	QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER *>(&HostQPCFrequency));
	QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER *>(&HostQPCStartTime));
	pit_last_qpc = last_qpc = HostQPCStartTime;
	pit_last = get_now();

	// Synchronize xbox system time with host time
	LARGE_INTEGER HostSystemTime;
	GetSystemTimeAsFileTime((LPFILETIME)&HostSystemTime);
	xbox::KeSystemTime.High2Time = HostSystemTime.u.HighPart;
	xbox::KeSystemTime.LowPart = HostSystemTime.u.LowPart;
	xbox::KeSystemTime.High1Time = HostSystemTime.u.HighPart;
}

// More precise sleep, but with increased CPU usage.
// Takes an absolute QPC target — no conversion, no drift.
void SleepPrecise(int64_t targetQPC)
{
	// Adaptive sleep strategy — every phase self-calibrates to never overshoot:
	// 1. Sleep() for the bulk, with margin based on worst-case Sleep() overshoot
	// 2. SwitchToThread() yielding, exits when remaining < 2x average yield duration
	// 3. Final tight spin for sub-yield precision

	// Max-tracked Sleep overshoot with slow decay (not EMA — avoids overshooting
	// on outlier spikes). Yield duration uses EMA since overshooting a single
	// yield just means one extra spin iteration, not a missed deadline.
	// Atomic for thread safety (called from system_events + PGRAPH puller).
	// Cache-line aligned to prevent false sharing between the two atomics
	// and with any adjacent static data.
	alignas(64) static std::atomic<int64_t> s_maxSleepOvershoot{HostQPCFrequency * 2 / 1000}; // init ~2ms
	alignas(64) static std::atomic<int64_t> s_avgYieldTicks{HostQPCFrequency / 1000};          // init ~1ms
	const int64_t kMaxYieldThreshold = HostQPCFrequency * 5 / 1000; // cap yield exit at ~5ms

	LARGE_INTEGER now;
	QueryPerformanceCounter(&now);

	// Early-out: target already passed
	if (now.QuadPart >= targetQPC)
		return;

	// Phase 1: Sleep() for the bulk, with margin based on worst-case overshoot
	int64_t avgYield = s_avgYieldTicks.load(std::memory_order_relaxed);
	int64_t yieldExit = avgYield * 2;
	if (yieldExit > kMaxYieldThreshold)
		yieldExit = kMaxYieldThreshold;
	int64_t remaining = targetQPC - now.QuadPart;
	int64_t sleepMargin = s_maxSleepOvershoot.load(std::memory_order_relaxed) + yieldExit;
	if (remaining > sleepMargin) {
		DWORD sleepMs = (DWORD)((remaining - sleepMargin) * 1000 / HostQPCFrequency);
		if (sleepMs > 0) {
			LARGE_INTEGER before = now;
			Sleep(sleepMs);
			QueryPerformanceCounter(&now);
			// Track worst-case overshoot with slow decay
			int64_t requestedTicks = (int64_t)sleepMs * HostQPCFrequency / 1000;
			int64_t overshoot = (now.QuadPart - before.QuadPart) - requestedTicks;
			if (overshoot < 0) overshoot = 0;
			int64_t prev = s_maxSleepOvershoot.load(std::memory_order_relaxed);
			if (overshoot > prev) {
				s_maxSleepOvershoot.store(overshoot, std::memory_order_relaxed);
			} else {
				// Slow decay: shrink by 1/64 per sample so it adapts down over time
				// Floor at 0.5ms to prevent near-zero margin after long stable periods
				int64_t decayed = prev - (prev >> 6);
				int64_t floor = HostQPCFrequency / 2000; // 0.5ms
				if (decayed < floor) decayed = floor;
				s_maxSleepOvershoot.store(decayed, std::memory_order_relaxed);
			}
		}
	}

	// Phase 2: Adaptive yield via SwitchToThread(), exit when remaining < 2x avg yield
	// Reuses post-yield QPC as next iteration's timestamp (no redundant QPC call).
	while (true) {
		remaining = targetQPC - now.QuadPart;
		if (remaining <= yieldExit)
			break;
		SwitchToThread();
		LARGE_INTEGER prev = now;
		QueryPerformanceCounter(&now);
		int64_t yieldTicks = now.QuadPart - prev.QuadPart;
		avgYield += (yieldTicks - avgYield) >> 3; // EMA 1/8
	}
	s_avgYieldTicks.store(avgYield, std::memory_order_relaxed);

	// Phase 3: Final tight spin — compare QPC directly, no clock domain crossing
	while (now.QuadPart < targetQPC) {
		QueryPerformanceCounter(&now);
	}
}

// NOTE: the pit device is not implemented right now, so we put this here
static uint64_t pit_next(uint64_t now)
{
	constexpr uint64_t pit_period = 1000;
	uint64_t next = pit_last + pit_period;

	if (now >= next) {
		xbox::KiClockIsr(now - pit_last);
		pit_last = get_now();
		return pit_period;
	}

	return pit_last + pit_period - now; // time remaining until next clock interrupt
}

static void update_non_periodic_events()
{
	// update dsound
	dsound_worker();

	// check for hw interrupts
	for (int i = 0; i < MAX_BUS_INTERRUPT_LEVEL; i++) {
		// Skip IRQ 3 (GPU/NV2A) - it's delivered explicitly by nv2a_vblank_interrupt
		// and the PGRAPH INTR_ERROR mechanism. Triggering it here races with the DPC
		// that re-enables PMC_INTR_EN_0, causing an ISR/DPC ping-pong deadlock.
		if (i == 3) continue;

		// If the interrupt is pending and connected, process it
		if (g_bEnableAllInterrupts && HalSystemInterrupts[i].IsPending() && EmuInterruptList[i] && EmuInterruptList[i]->Connected) {
			HalSystemInterrupts[i].Trigger(EmuInterruptList[i]);
		}
	}
}

uint64_t get_now()
{
	LARGE_INTEGER now;
	QueryPerformanceCounter(&now);
	uint64_t elapsed_us = now.QuadPart - last_qpc;
	last_qpc = now.QuadPart;
	elapsed_us *= 1000000;
	elapsed_us /= HostQPCFrequency;
	exec_time += elapsed_us;
	return exec_time;
}

static uint64_t get_next(uint64_t now)
{
	std::array<uint64_t, 5> next = {
		pit_next(now),
		g_NV2A->vblank_next(now),
		g_NV2A->ptimer_next(now),
		g_USB0->m_HostController->OHCI_next(now),
		dsound_next(now)
	};
	return *std::min_element(next.begin(), next.end());
}

xbox::void_xt NTAPI system_events(xbox::PVOID arg)
{
	// Testing shows that, if this thread has the same priority of the other xbox threads, it can take tens, even hundreds of ms to complete a single loop.
	// So we increase its priority to above normal, so that it scheduled more often
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);

	// Always run this thread at dpc level to prevent it from ever executing APCs/DPCs
	xbox::KeRaiseIrqlToDpcLevel();

	// Persistent anchor: after SleepPrecise spins to targetQPC, we use
	// that exact value as the base for the next iteration — no fresh QPC
	// read in between, so no accumulating drift.
	LARGE_INTEGER qpc;
	QueryPerformanceCounter(&qpc);
	int64_t wall_anchor = qpc.QuadPart;

	while (true) {
		LARGE_INTEGER loop_start;
		if (g_bCxbxProfilerEnabled) QueryPerformanceCounter(&loop_start);

		const uint64_t last_time = get_now();
		const uint64_t nearest_next = get_next(last_time);

		// Process non-periodic events once at the start of each cycle
		update_non_periodic_events();

		// Wait precisely for the next periodic event deadline.
		// Target is anchored to the previous SleepPrecise wake-up, not "now".
		if (nearest_next > 0) {
			int64_t targetQPC = wall_anchor
				+ (int64_t)nearest_next * HostQPCFrequency / 1000000;
			SleepPrecise(targetQPC);
			wall_anchor = targetQPC; // exact wake-up becomes next anchor
		}

		// Process non-periodic events again after waking (handles any
		// that arrived during the sleep)
		update_non_periodic_events();

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

