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
// *  All rights reserved
// *
// ******************************************************************
#pragma once

#include "core/kernel/common/types.h"
#include <mutex>

// ******************************************************************
// * ETIMER - Executive Timer object (wraps KTIMER for Nt*Timer API)
// ******************************************************************
// Source: ReactOS, simplified for Xbox compatibility layer (no wake
// timers, no per-thread active timer list)
typedef struct _ETIMER {
	xbox::KTIMER KeTimer;          // Embedded kernel timer
	xbox::KAPC TimerApc;           // APC to queue when timer fires
	xbox::KDPC TimerDpc;           // DPC that queues the APC
	std::mutex Lock;               // Protects ApcAssociated state
	xbox::boolean_xt ApcAssociated;// TRUE when TimerApc is active
	xbox::long_xt Period;          // Periodic interval (ms), 0 = one-shot
} ETIMER, *PETIMER;

namespace xbox
{

// Executive timer helper functions (implemented in EmuKrnlEx.cpp)
void_xt NTAPI ExpDeleteTimer(IN PVOID ObjectBody);
void_xt NTAPI ExpTimerDpcRoutine(IN PKDPC Dpc, IN PVOID DeferredContext, IN PVOID SystemArgument1, IN PVOID SystemArgument2);
void_xt NTAPI ExpTimerApcKernelRoutine(IN PKAPC Apc, IN PKNORMAL_ROUTINE *NormalRoutine, IN PVOID *NormalContext, IN PVOID *SystemArgument1, IN PVOID *SystemArgument2);

// Executive mutant helper functions (implemented in EmuKrnlEx.cpp)
void_xt NTAPI ExpDeleteMutant(IN PVOID ObjectBody);

};
