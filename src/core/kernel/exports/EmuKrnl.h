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
#include "EmuKrnlKi.h"
#include <future>
#include <cstdio>

#define OBJECT_TO_OBJECT_HEADER(Object) \
    CONTAINING_RECORD(Object, OBJECT_HEADER, Body)

void InitializeListHead(xbox::PLIST_ENTRY pListHead);
bool IsListEmpty(xbox::PLIST_ENTRY pListHead);
void InsertHeadList(xbox::PLIST_ENTRY pListHead, xbox::PLIST_ENTRY pEntry);
void InsertTailList(xbox::PLIST_ENTRY pListHead, xbox::PLIST_ENTRY pEntry);

xbox::boolean_xt RemoveEntryList(xbox::PLIST_ENTRY pEntry);
xbox::PLIST_ENTRY RemoveHeadList(xbox::PLIST_ENTRY pListHead);
xbox::PLIST_ENTRY RemoveTailList(xbox::PLIST_ENTRY pListHead);

extern xbox::LAUNCH_DATA_PAGE DefaultLaunchDataPage;
extern xbox::PKINTERRUPT EmuInterruptList[MAX_BUS_INTERRUPT_LEVEL + 1];
extern xbox::PKINTERRUPT EmuInterruptChained[MAX_BUS_INTERRUPT_LEVEL + 1];
inline std::atomic_bool g_bEnableAllInterrupts = true;

class HalSystemInterrupt {
public:
	void Assert(bool state) {
		if (m_Asserted == 0 && state == 1) m_Pending = true;
		m_Asserted = state;
	};
	void Enable() { m_Enabled = true; }
	void Disable() { m_Enabled = false; }
	bool IsEnabled() { return m_Enabled; }
	bool IsPending() { return m_Asserted && m_Pending; }
	void SetInterruptMode(xbox::KINTERRUPT_MODE InterruptMode) { m_InterruptMode = InterruptMode; }
	void Trigger(xbox::PKINTERRUPT Interrupt) {
		if (m_InterruptMode == xbox::KINTERRUPT_MODE::LevelSensitive) m_Pending = false;

		// On real Xbox, the ISR runs at its assigned device IRQL (above
		// DISPATCH_LEVEL). This prevents KeInsertQueueDpc from dispatching
		// DPCs inline during the ISR — DPCs are queued and fire later when
		// IRQL drops below DISPATCH_LEVEL. Without this, the ISR's
		// KeInsertQueueDpc calls see IRQL < DISPATCH_LEVEL and dispatch
		// immediately, causing tight loops for DPCs that re-queue themselves
		// or for ISRs that fire repeatedly.
		volatile xbox::KPCR* Pcr = EmuKeGetPcr();
		xbox::KIRQL OldIrql = (xbox::KIRQL)Pcr->Irql;
		xbox::KIRQL IsrIrql = (xbox::KIRQL)Interrupt->Irql;
		if (IsrIrql < DISPATCH_LEVEL) IsrIrql = DISPATCH_LEVEL;
		if (IsrIrql > OldIrql) Pcr->Irql = IsrIrql;

		auto ServiceRoutine = (xbox::boolean_xt(__stdcall*)(xbox::PKINTERRUPT, void*))Interrupt->ServiceRoutine;
		ServiceRoutine(Interrupt, Interrupt->ServiceContext);

		Pcr->Irql = OldIrql;
	}
private:
	bool m_Asserted = false, m_Enabled = false, m_Pending = false;
	xbox::KINTERRUPT_MODE m_InterruptMode;
};

extern HalSystemInterrupt HalSystemInterrupts[MAX_BUS_INTERRUPT_LEVEL + 1];

bool DisableInterrupts();
void RestoreInterruptMode(bool value);
void CallSoftwareInterrupt(const xbox::KIRQL SoftwareIrql);
bool AddWaitObject(xbox::PKTHREAD kThread, xbox::PLARGE_INTEGER Timeout);

template<typename T>
std::optional<xbox::ntstatus_xt> SatisfyWait(T &&Lambda, xbox::PKTHREAD kThread,
	xbox::boolean_xt Alertable, xbox::char_xt WaitMode)
{
	if (const auto ret = Lambda(kThread)) return ret;
	xbox::KiApcListMtx.lock();
	bool EmptyKernel = IsListEmpty(&kThread->ApcState.ApcListHead[xbox::KernelMode]);
	bool EmptyUser = IsListEmpty(&kThread->ApcState.ApcListHead[xbox::UserMode]);
	xbox::KiApcListMtx.unlock();
	if (!EmptyKernel) xbox::KiExecuteKernelApc();
	if (!EmptyUser && Alertable == TRUE && WaitMode == xbox::UserMode) {
		xbox::KiExecuteUserApc();
		xbox::KiUnwaitThreadAndLock(kThread, X_STATUS_USER_APC, 0);
		return kThread->WaitStatus;
	}
	return std::nullopt;
}

template<bool host_wait, typename T>
xbox::ntstatus_xt WaitApc(T &&Lambda, xbox::PLARGE_INTEGER Timeout,
	xbox::boolean_xt Alertable, xbox::char_xt WaitMode, xbox::PKTHREAD kThread)
{
	xbox::ntstatus_xt status;
	if (Timeout == nullptr) {
		HANDLE hWake = CxbxGetThreadWakeEvent(kThread);
		int stallIters = 0;
		while (true) {
			if (const auto ret = SatisfyWait(Lambda, kThread, Alertable, WaitMode)) {
				status = *ret; break;
			}
			if (hWake) WaitForSingleObjectEx(hWake, 10, TRUE);
			else SleepEx(1, TRUE);
			if (++stallIters == 300) {
				auto* objPtr = (unsigned char*)kThread->WaitBlockList->Object;
				fprintf(stderr, "[WAIT-STALL] tid=0x%X obj=0x%p type=%d signalState=%d state=%d\n",
					GetCurrentThreadId(), kThread->WaitBlockList->Object,
					objPtr[0], *(int*)(objPtr + 4), (int)kThread->State);
				fflush(stderr);
			}
		}
	}
	else if (Timeout->QuadPart == 0) {
		assert(host_wait);
		if (const auto ret = SatisfyWait(Lambda, kThread, Alertable, WaitMode)) status = *ret;
		else {
			xbox::KiUnwaitThreadAndLock(kThread, X_STATUS_TIMEOUT, 0);
			status = kThread->WaitStatus;
		}
	}
	else {
		HANDLE hWake = CxbxGetThreadWakeEvent(kThread);
		int finiteStallIters = 0;
		while (true) {
			if (const auto ret = SatisfyWait(Lambda, kThread, Alertable, WaitMode)) {
				status = *ret; break;
			}
			if (host_wait && kThread->State == xbox::Ready) {
				status = kThread->WaitStatus; break;
			}
			if (hWake) WaitForSingleObjectEx(hWake, 10, TRUE);
			else SleepEx(1, TRUE);
			if (++finiteStallIters == 300) {
				auto* objPtr2 = (unsigned char*)kThread->WaitBlockList->Object;
				fprintf(stderr, "[WAIT-STALL-FIN] tid=0x%X obj=0x%p type=%d timeout=%lld state=%d\n",
					GetCurrentThreadId(), kThread->WaitBlockList->Object,
					objPtr2[0], Timeout->QuadPart, (int)kThread->State);
				fflush(stderr);
			}
		}
	}
	if constexpr (host_wait) kThread->State = xbox::Running;
	return status;
}

#endif
