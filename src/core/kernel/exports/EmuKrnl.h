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
// *  You should have received a copy of the GNU General Public License
// *  along with this program; see the file COPYING.
// *  If not, write to the Free Software Foundation, Inc.,
// *  59 Temple Place - Suite 330, Bostom, MA 02111-1307, USA.
// *
// *  (c) 2002-2003 Aaron Robinson <caustik@caustik.com>
// *
// *  All rights reserved
// *
// ******************************************************************
#ifndef EMUKRNL_H
#define EMUKRNL_H

#include "core\kernel\init\CxbxKrnl.h"
#include "core\kernel\support\Emu.h"
#include "core\kernel\support\EmuFS.h"
#include "devices\audio\APUDevice.h"

extern class APUDevice* g_APU;
#include "EmuKrnlKi.h"
#include "core\hle\DSOUND\DirectSound\ApuPlayCursor.h"
#include <future>
#include <cstdio>

// CONTAINING_RECORD macro
// Gets the value of structure member (field - num1),given the type(MYSTRUCT, in this code) and the List_Entry head(temp, in this code)
// See https://stackoverflow.com/questions/8240273/a-portable-way-to-calculate-pointer-to-the-whole-structure-using-pointer-to-a-fi
//#define CONTAINING_RECORD(ptr, type, field) \
//	(((type) *)((char *)(ptr) - offsetof((type), member)))

#define OBJECT_TO_OBJECT_HEADER(Object) \
    CONTAINING_RECORD(Object, OBJECT_HEADER, Body)

void InitializeListHead(xbox::PLIST_ENTRY pListHead);
bool IsListEmpty(xbox::PLIST_ENTRY pListHead);
void InsertHeadList(xbox::PLIST_ENTRY pListHead, xbox::PLIST_ENTRY pEntry);
void InsertTailList(xbox::PLIST_ENTRY pListHead, xbox::PLIST_ENTRY pEntry);
//#define RemoveEntryList(e) do { PLIST_ENTRY f = (e)->Flink, b = (e)->Blink; f->Blink = b; b->Flink = f; (e)->Flink = (e)->Blink = NULL; } while (0)

xbox::boolean_xt RemoveEntryList(xbox::PLIST_ENTRY pEntry);
xbox::PLIST_ENTRY RemoveHeadList(xbox::PLIST_ENTRY pListHead);
xbox::PLIST_ENTRY RemoveTailList(xbox::PLIST_ENTRY pListHead);

extern xbox::LAUNCH_DATA_PAGE DefaultLaunchDataPage;
extern xbox::PKINTERRUPT EmuInterruptList[MAX_BUS_INTERRUPT_LEVEL + 1];
extern xbox::PKINTERRUPT EmuInterruptChained[MAX_BUS_INTERRUPT_LEVEL + 1];
// Indicates to disable/enable all interrupts when cli and sti instructions are executed
inline std::atomic_bool g_bEnableAllInterrupts = true;

class HalSystemInterrupt {
public:
	void Assert(bool state) {
		// If the interrupt was marked as Asserted, and was previously not, set the pending flag too!
		if (m_Asserted == 0 && state == 1) {
			m_Pending = true;
		}

		m_Asserted = state;
	};

	void Enable() {
		m_Enabled = true;
	}

	void Disable() {
		m_Enabled = false;
	}

	bool IsEnabled() {
		return m_Enabled;
	}

	bool IsPending() {
		return m_Asserted && m_Pending;
	}

	void SetInterruptMode(xbox::KINTERRUPT_MODE InterruptMode) {
		m_InterruptMode = InterruptMode;
	}

	void Trigger(xbox::PKINTERRUPT Interrupt) {
		// If interrupt was level sensitive, we clear the pending flag, preventing the interrupt from being triggered 
		// until it is deasserted then asserted again. Latched interrupts are triggered until the line is Deasserted!
		if (m_InterruptMode == xbox::KINTERRUPT_MODE::LevelSensitive) {
			m_Pending = false;
		}

		xbox::boolean_xt(__stdcall *ServiceRoutine)(xbox::PKINTERRUPT, void*) = (xbox::boolean_xt(__stdcall *)(xbox::PKINTERRUPT, void*))Interrupt->ServiceRoutine;
		xbox::boolean_xt result = ServiceRoutine(Interrupt, Interrupt->ServiceContext);
	}
private:
	bool m_Asserted = false;
	bool m_Enabled = false;
	xbox::KINTERRUPT_MODE m_InterruptMode;
	bool m_Pending = false;
};

extern HalSystemInterrupt HalSystemInterrupts[MAX_BUS_INTERRUPT_LEVEL + 1];

bool DisableInterrupts();
void RestoreInterruptMode(bool value);
void CallSoftwareInterrupt(const xbox::KIRQL SoftwareIrql);
bool AddWaitObject(xbox::PKTHREAD kThread, xbox::PLARGE_INTEGER Timeout);

template<typename T>
std::optional<xbox::ntstatus_xt> SatisfyWait(T &&Lambda, xbox::PKTHREAD kThread, xbox::boolean_xt Alertable, xbox::char_xt WaitMode)
{
	if (const auto ret = Lambda(kThread)) {
		return ret;
	}

	xbox::KiApcListMtx.lock();
	bool EmptyKernel = IsListEmpty(&kThread->ApcState.ApcListHead[xbox::KernelMode]);
	bool EmptyUser = IsListEmpty(&kThread->ApcState.ApcListHead[xbox::UserMode]);
	xbox::KiApcListMtx.unlock();

	if (EmptyKernel == false) {
		xbox::KiExecuteKernelApc();
	}

	if ((EmptyUser == false) &&
		(Alertable == TRUE) &&
		(WaitMode == xbox::UserMode)) {
		xbox::KiExecuteUserApc();
		xbox::KiUnwaitThreadAndLock(kThread, X_STATUS_USER_APC, 0);
		return kThread->WaitStatus;
	}

	return std::nullopt;
}

template<bool host_wait, typename T>
xbox::ntstatus_xt WaitApc(T &&Lambda, xbox::PLARGE_INTEGER Timeout, xbox::boolean_xt Alertable, xbox::char_xt WaitMode, xbox::PKTHREAD kThread)
{
	// NOTE1: kThread->Alerted is currently never set. When the alerted mechanism is implemented, the alerts should
	// also interrupt the wait.

	xbox::ntstatus_xt status;
	if (Timeout == nullptr) {
		// No timout specified, so this is an infinite wait until an alert, a user apc or the object(s) become(s) signalled
		HANDLE hWake = CxbxGetThreadWakeEvent(kThread);
		int stallIters = 0;
		while (true) {
			if (const auto ret = SatisfyWait(Lambda, kThread, Alertable, WaitMode)) {
				status = *ret;
				break;
			}

			// Block on the per-thread wake event instead of polling.
			// The event is signaled by KiUnwaitThread (when the waited
			// object becomes signaled or a timeout fires) and by
			// KiInsertQueueApc (when an APC is queued to this thread).
			// Use alertable wait so host I/O completion APCs still work.
			// Use a timeout instead of INFINITE to handle the case where
			// game code directly writes to Header.SignalState (bypassing
			// KeSetEvent), which would not trigger KiWaitTest/KiUnwaitThread.
			if (hWake) {
				WaitForSingleObjectEx(hWake, 10, TRUE);
			} else {
				SleepEx(1, TRUE);
			}

			// Stall diagnostic: if we've been waiting >3 seconds, scan for the
			// APU play cursor pointer.
			if (++stallIters == 300) {
				auto* objPtr = (unsigned char*)kThread->WaitBlockList->Object;
				int objType = objPtr[0];
				int signalState = *(int*)(objPtr + 4);
				int waiters = !(((uintptr_t*)(objPtr + 8))[0] == (uintptr_t)(objPtr + 8));
				fprintf(stderr, "[WAIT-STALL] tid=0x%X obj=0x%p type=%d signalState=%d waiters=%d state=%d\n",
					GetCurrentThreadId(), kThread->WaitBlockList->Object,
					objType, signalState, waiters, (int)kThread->State);
				fflush(stderr);

				uintptr_t eventAddr = (uintptr_t)kThread->WaitBlockList->Object;
				if (eventAddr >= 0x00011000 && eventAddr < 0x00300000) {
					DWORD* expectedP = reinterpret_cast<DWORD*>(eventAddr - 0x2530);
					DWORD expectedVal = *expectedP;
					if (expectedVal >= 0x80000000 && expectedVal < 0x84000000 && (expectedVal & 0x3) == 0) {
						DWORD expectedCbo = *reinterpret_cast<volatile DWORD*>(expectedVal);
						if (expectedCbo < 0x100000) {
							LARGE_INTEGER qpc;
							QueryPerformanceCounter(&qpc);
							g_ApuPlayCursor.pCursor = reinterpret_cast<volatile DWORD*>(expectedVal);
							g_ApuPlayCursor.bufSize = 0x200000;
							g_ApuPlayCursor.rate = 96000;
							g_ApuPlayCursor.lastQPC = qpc.QuadPart;
							g_ApuPlayCursor.pEvent = kThread->WaitBlockList->Object;
							fprintf(stderr, "  [APU-CURSOR] Activated at XDK offset: cursor=0x%08X cbo=%u event=0x%p\n",
								expectedVal, expectedCbo, g_ApuPlayCursor.pEvent);
							fflush(stderr);

							if (g_APU) {
								g_APU->SetFallbackVoiceBase(
									reinterpret_cast<uint8_t*>(expectedVal - 0x58));
							}

							// Diagnostic: directly write to CBO and verify it persists
							volatile uint32_t* pCbo = reinterpret_cast<volatile uint32_t*>(expectedVal);
							uint32_t cboBefore = *pCbo;
							*pCbo = cboBefore + 48000; // advance by ~1s of mono 48kHz audio
							uint32_t cboAfter = *pCbo;
							fprintf(stderr, "  [APU-WRITE] Direct CBO write: before=%u wrote=%u after=%u\n",
								cboBefore, cboBefore + 48000, cboAfter);
							fflush(stderr);
						}
					}
					if (g_ApuPlayCursor.pCursor == nullptr) {
					// Fallback scan
					DWORD* scanStart = reinterpret_cast<DWORD*>(
						(eventAddr > 0x4000) ? (eventAddr - 0x4000) : 0x00011000);
					DWORD* scanEnd = reinterpret_cast<DWORD*>(eventAddr);
					for (DWORD* p = scanEnd - 1; p >= scanStart; p--) {
						DWORD val = *p;
						if (val >= 0x80000000 && val < 0x84000000 && (val & 0x3) == 0) {
							DWORD cbo = *reinterpret_cast<volatile DWORD*>(val);
							if (cbo < 0x100000) {
								LARGE_INTEGER qpc;
								QueryPerformanceCounter(&qpc);
								g_ApuPlayCursor.pCursor = reinterpret_cast<volatile DWORD*>(val);
								g_ApuPlayCursor.bufSize = 0x200000;
								g_ApuPlayCursor.rate = 96000;
								g_ApuPlayCursor.lastQPC = qpc.QuadPart;
								g_ApuPlayCursor.pEvent = kThread->WaitBlockList->Object;
								fprintf(stderr, "  [APU-CURSOR] Activated: cursor=0x%08X cbo=%u event=0x%p\n",
									val, cbo, g_ApuPlayCursor.pEvent);
								fflush(stderr);
								break;
							}
						}
					}
					} /* g_ApuPlayCursor.pCursor == nullptr */
				}
			}
		}
	}
	else if (Timeout->QuadPart == 0) {
		assert(host_wait);
		// A zero timeout means that we only have to check the conditions once and then return immediately if they are not satisfied
		if (const auto ret = SatisfyWait(Lambda, kThread, Alertable, WaitMode)) {
			status = *ret;
		}
		else {
			// If the wait failed, then always remove the wait block. Note that this can only happen with host waits, since guest waits never call us at all
			// when Timeout->QuadPart == 0. Test case: Halo 2 (sporadically when playing the intro video)
			xbox::KiUnwaitThreadAndLock(kThread, X_STATUS_TIMEOUT, 0);
			status = kThread->WaitStatus;
		}
	}
	else {
		// A non-zero timeout means we have to check the conditions until we reach the requested time.
		// The kernel timer (set up by the caller) will fire KiTimerExpiration → KiUnwaitThread
		// which signals our wake event, so we can block efficiently here.
		HANDLE hWake = CxbxGetThreadWakeEvent(kThread);
		int finiteStallIters = 0;
		while (true) {
			if (const auto ret = SatisfyWait(Lambda, kThread, Alertable, WaitMode)) {
				status = *ret;
				break;
			}

			if (host_wait && (kThread->State == xbox::Ready)) {
				status = kThread->WaitStatus;
				break;
			}

			if (hWake) {
				WaitForSingleObjectEx(hWake, 10, TRUE);
			} else {
				SleepEx(1, TRUE);
			}

			if (++finiteStallIters == 300) {
				auto* objPtr2 = (unsigned char*)kThread->WaitBlockList->Object;
				fprintf(stderr, "[WAIT-STALL-FIN] tid=0x%X obj=0x%p type=%d signalState=%d timeout=%lld state=%d\n",
					GetCurrentThreadId(), kThread->WaitBlockList->Object,
					objPtr2[0], *(int*)(objPtr2 + 4), Timeout->QuadPart, (int)kThread->State);
				fflush(stderr);

				uintptr_t eventAddr = (uintptr_t)kThread->WaitBlockList->Object;
				if (eventAddr >= 0x00011000 && eventAddr < 0x00300000) {
					DWORD* expectedP = reinterpret_cast<DWORD*>(eventAddr - 0x2530);
					DWORD expectedVal = *expectedP;
					if (expectedVal >= 0x80000000 && expectedVal < 0x84000000 && (expectedVal & 0x3) == 0) {
						DWORD expectedCbo = *reinterpret_cast<volatile DWORD*>(expectedVal);
						if (expectedCbo < 0x100000) {
							LARGE_INTEGER qpc;
							QueryPerformanceCounter(&qpc);
							g_ApuPlayCursor.pCursor = reinterpret_cast<volatile DWORD*>(expectedVal);
							g_ApuPlayCursor.bufSize = 0x200000;
							g_ApuPlayCursor.rate = 96000;
							g_ApuPlayCursor.lastQPC = qpc.QuadPart;
							g_ApuPlayCursor.pEvent = kThread->WaitBlockList->Object;
							fprintf(stderr, "  [APU-CURSOR] Activated at XDK offset: cursor=0x%08X cbo=%u event=0x%p\n",
								expectedVal, expectedCbo, g_ApuPlayCursor.pEvent);
							fflush(stderr);
						}
					}
					if (g_ApuPlayCursor.pCursor == nullptr) {
					DWORD* scanStart = reinterpret_cast<DWORD*>(
						(eventAddr > 0x4000) ? (eventAddr - 0x4000) : 0x00011000);
					DWORD* scanEnd = reinterpret_cast<DWORD*>(eventAddr);
					for (DWORD* p = scanEnd - 1; p >= scanStart; p--) {
						DWORD val = *p;
						if (val >= 0x80000000 && val < 0x84000000 && (val & 0x3) == 0) {
							DWORD cbo = *reinterpret_cast<volatile DWORD*>(val);
							if (cbo < 0x100000) {
								LARGE_INTEGER qpc;
								QueryPerformanceCounter(&qpc);
								g_ApuPlayCursor.pCursor = reinterpret_cast<volatile DWORD*>(val);
								g_ApuPlayCursor.bufSize = 0x200000;
								g_ApuPlayCursor.rate = 96000;
								g_ApuPlayCursor.lastQPC = qpc.QuadPart;
								g_ApuPlayCursor.pEvent = kThread->WaitBlockList->Object;
								fprintf(stderr, "  [APU-CURSOR] Activated: cursor=0x%08X cbo=%u event=0x%p\n",
									val, cbo, g_ApuPlayCursor.pEvent);
								fflush(stderr);
								break;
							}
						}
					}
					} /* g_ApuPlayCursor.pCursor == nullptr */
				}
			}
		}
	}

	if constexpr (host_wait) {
		kThread->State = xbox::Running;
	}
	return status;
}

#endif
