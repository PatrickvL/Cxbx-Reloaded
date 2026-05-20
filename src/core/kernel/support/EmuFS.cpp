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

#define LOG_PREFIX CXBXR_MODULE::FS


#include <core\kernel\exports\xboxkrnl.h>
#include "core\kernel\exports\EmuKrnl.h" // For InitializeListHead(), etc.
#include "core\kernel\exports\EmuKrnlKe.h"
#include "core\kernel\exports\EmuKrnlKi.h"
#include "core\kernel\support\EmuFS.h" // For fs_instruction_t
#include "core\kernel\support\NativeHandle.h"
#include "core\kernel\init\CxbxKrnl.h"
#include "Logging.h"

#include <windows.h>
#include <cstdio>
#include <vector>

// Global DpcRoutineActive flag (single-CPU emulation, defined in EmuKrnlKe.cpp)
extern volatile xbox::ulong_xt g_DpcRoutineActive;

#ifdef RtlZeroMemory
#undef RtlZeroMemory
#endif

#define TLS_ALIGNMENT_OFFSET 12

// NT_TIB (Thread Information Block) offsets - see https://www.microsoft.com/msj/archive/S2CE.aspx
#define TIB_ExceptionList         offsetof(NT_TIB, ExceptionList)         // = 0x00/0
#define TIB_StackBase             offsetof(NT_TIB, StackBase)             // = 0x04/4
#define TIB_StackLimit            offsetof(NT_TIB, StackLimit)            // = 0x08/8
#define TIB_SubSystemTib          offsetof(NT_TIB, SubSystemTib)          // = 0x0C/12
#define TIB_FiberData             offsetof(NT_TIB, FiberData)             // = 0x10/16
#define TIB_ArbitraryUserPointer  offsetof(NT_TIB, ArbitraryUserPointer)  // = 0x14/20
#define TIB_Self                  offsetof(NT_TIB, Self)                  // = 0x18/24

// KPCR (Kernel Processor Control Region) offsets
#define KPCR_NtTib                offsetof(KPCR, NtTib)                   // = 0x00/0
#define KPCR_SelfPcr              offsetof(KPCR, SelfPcr)                 // = 0x1C/28
#define KPCR_Prcb                 offsetof(KPCR, Prcb)                    // = 0x20/32
#define KPCR_Irql                 offsetof(KPCR, Irql)                    // = 0x24/36
#define KPCR_PrcbData             offsetof(KPCR, PrcbData)                // = 0x28/40

// KPRCB (Kernel PRocesor Control Block) offsets
#define KPRCB_CurrentThread       offsetof(KPRCB, CurrentThread)          // = 0x00, KPCR : 0x28/40
#define KPRCB_NextThread          offsetof(KPRCB, NextThread)             // = 0x04, KPCR : 0x2C/44
#define KPRCB_IdleThread          offsetof(KPRCB, IdleThread)             // = 0x08, KPCR : 0x30/48
#define KPRCB_DpcListHead         offsetof(KPRCB, DpcListHead)            // = 0x28, KPCR : 0x50/80
#define KPRCB_DpcRoutineActive    offsetof(KPRCB, DpcRoutineActive)       // = 0x30, KPCR : 0x58/88

// KTHREAD (Kernel Thread) offsets
#define KTHREAD_Header            offsetof(KTHREAD, Header)               // = 0x0/0
#define KTHREAD_MutantListHead    offsetof(KTHREAD, MutantListHead)       // = 0x10/16
#define KTHREAD_KernelTime        offsetof(KTHREAD, KernelTime)           // = 0x18/24
#define KTHREAD_StackBase         offsetof(KTHREAD, StackBase)            // = 0x1C/28
#define KTHREAD_StackLimit        offsetof(KTHREAD, StackLimit)           // = 0x20/32
#define KTHREAD_KernelStack       offsetof(KTHREAD, KernelStack)          // = 0x24/36
#define KTHREAD_TlsData           offsetof(KTHREAD, TlsData)              // = 0x28/40
#define KTHREAD_State             offsetof(KTHREAD, State)                // = 0x2C/44
#define KTHREAD_Alerted           offsetof(KTHREAD, Alerted)              // = 0x2D/45
#define KTHREAD_Alertable         offsetof(KTHREAD, Alertable)            // = 0x2F/47
#define KTHREAD_NpxState          offsetof(KTHREAD, NpxState)             // = 0x30/48
// = 0x31/49 */ CHAR Saturation;
// = 0x32/50 */ CHAR Priority;
// = 0x33/51 */ CHAR Padding;
// = 0x34/52 */ KAPC_STATE ApcState;
// = 0x4C/76 */ ULONG ContextSwitches;
// = 0x50/80 */ ULONG WaitStatus;
// = 0x54/84 */ CHAR WaitIrql;
// = 0x55/85 */ CHAR WaitMode;
// = 0x56/86 */ CHAR WaitNext;
// = 0x57/87 */ CHAR WaitReason;
// = 0x58/88 */ PVOID WaitBlockList;
// = 0x5C/92 */ LIST_ENTRY WaitListEntry;
// = 0x64/100 */ ULONG WaitTime;
// = 0x68/104 */ ULONG KernelApcDisable;
// = 0x6C/108 */ ULONG Quantum;
// = 0x70/112 */ CHAR BasePriority;
// = 0x71/113 */ CHAR DecrementCount;
// = 0x72/114 */ CHAR PriorityDecrement;
// = 0x73/115 */ CHAR DisableBoost;
// = 0x74/116 */ CHAR NpxIrql;
// = 0x75/117 */ CHAR SuspendCount;
// = 0x76/118 */ CHAR Preempted;
// = 0x77/119 */ CHAR HasTerminated;
// = 0x78/120 */ PVOID Queue;
// = 0x7C/124 */ LIST_ENTRY QueueListEntry;
// = 0x88/136 */ UCHAR rsvd1[4];
// = 0x88/136 */ KTIMER Timer;
// = 0xB0/176 */ KWAIT_BLOCK TimerWaitBlock;
// = 0xC8/200 */ KAPC SuspendApc;
// = 0xF0/240 */ KSEMAPHORE SuspendSemaphore;
// = 0x104/260 */ LIST_ENTRY ThreadListEntry;
// = 0x10C/268 */ UCHAR _padding[4];

template void EmuGenerateFS<true>(xbox::PETHREAD Ethread, unsigned XboxStackBaseReserved, unsigned XboxStackSizeReserved);
template void EmuGenerateFS<false>(xbox::PETHREAD Ethread, unsigned XboxStackBaseReserved, unsigned XboxStackSizeReserved);

NT_TIB *GetNtTib()
{
	return (NT_TIB *)__readfsdword(TIB_LinearSelfAddress);
}

// ******************************************************************
// * FS Thunks — Trailing-offset-byte scheme
// *
// * Each patched instruction is overwritten as:
// *   CALL rel32 (5 bytes) + KPCR offset byte (1 byte) + NOP padding
// *
// * The thunk reads the offset byte from [return_address], increments
// * the return address past it, then uses it as [KPCR_base + offset].
// * One thunk per register per direction handles all valid KPCR offsets,
// * reducing ~50 naked functions to 12 regular + 5 special = 17 total.
// ******************************************************************

// Valid KPCR offsets that Xbox code may access via fs:[]
static const uint8_t kValidKpcrOffsets[] = {
	0x00, // KPCR.NtTib              — ExceptionList (start of NT_TIB)
	0x04, // KPCR.NtTib.StackBase    — Thread stack base
	0x20, // KPCR.Prcb              — Pointer to KPRCB
	0x28, // KPCR.PrcbData          — Inline KPRCB (CurrentThread, etc.)
	0x58, // KPCR.PrcbData.DpcRoutineActive
};

void EmuKeSetPcr(xbox::KPCR *Pcr)
{
	// Store the Xbox KPCR pointer in FS (See EmuKeGetPcr())
	// 
	// Note : Cxbx currently doesn't do preemptive thread switching,
	// which implies that thread-state management is done by Windows.
	//
	// Xbox executable code expects thread-specific state data to
	// be available via the FS segment register. To emulate this,
	// Cxbx uses the user data-slot feature of Windows threads.
	// 
	// Cxbx puts a pointer to a thread-specific copy of an entire
	// Kernel Processor Control Region (KPCR) into this data-slot.
	//
	// In the Xbox there's only be KPCR (as it's a per-processor-
	// structure, and the Xbox has only one processor).
	//
	// Since Cxbx doesn't control thread-switches (yet), each thread
	// must have a thread-specific copy of the KPCR, to contain all
	// thread-specific data that can be reached via this structure
	// (like the NT_TIB structure and ETHREAD CurrentThread pointer).
	// 
	// For this to work, Cxbx patches all executable code accessing
	// the FS segment register, so that the KPCR is accessed via
	// the user data-slot of each Windows thread Cxbx uses for an
	// Xbox thread.
	//
	__writefsdword(TIB_ArbitraryDataSlot, (DWORD)Pcr);
}

void EmuKeFreePcr()
{
	using namespace xbox;
	PVOID Pcr = (PVOID)EmuKeGetPcr();
	ulong_xt Size = zero;
	ntstatus_xt Status = NtFreeVirtualMemory(&Pcr, &Size, XBOX_MEM_RELEASE); // free pcr
	assert(Status == X_STATUS_SUCCESS);
	__writefsdword(TIB_ArbitraryDataSlot, NULL);
}

__declspec(naked) void EmuFS_RefreshKPCR()
{
	// Backup all registers, call EmuKeGetPcr and then restore all registers
	// EmuKeGetPcr makes sure a valid KPCR exists for the current thread
	// and creates it if missing, we backup and restore all registers
	// to keep it safe to call in our patches
	// This function can be later expanded to do nice things 
	// like setup the per-thread KPCR values for us too!
	__asm {
		pushfd
		pushad
		call EmuKeGetPcr
		popad
		popfd
		ret
	}
}

// ******************************************************************
// * Special thunks (exact-match, no trailing offset byte)
// ******************************************************************

__declspec(naked) void EmuFS_CmpEsiFs00()
{
	// cmp esi, large fs:0 — compare esi with KPCR[0] (ExceptionList)
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		mov eax, fs : [TIB_ArbitraryDataSlot]
		cmp esi, [eax]
		pop eax
		ret
	}
}

__declspec(naked) void EmuFS_MovzxEaxBytePtrFs24()
{
	// movzx eax, byte ptr fs:[24h] — inlined KeGetCurrentIrql()
	__asm
	{
		call EmuFS_RefreshKPCR
		mov eax, fs : [TIB_ArbitraryDataSlot]
		movzx eax, byte ptr[eax + 24h]
		ret
	}
}

__declspec(naked) void EmuFS_WriteEsp()
{
	// mov fs:[offset], esp — stores adjusted ESP to KPCR[offset]
	// Uses trailing-offset-byte scheme but also needs flags preservation
	// and ESP adjustment to compensate for CALL + pushfd overhead.
	__asm
	{
		pushfd
		call EmuFS_RefreshKPCR
		push eax
		push ecx
		// Stack: [ecx, eax, flags, ret_addr, <caller>]
		//         +0   +4   +8     +12
		mov ecx, [esp + 12]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 12]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov[eax + ecx], esp
		add dword ptr[eax + ecx], 16  // adjust: ecx(4) + eax(4) + flags(4) + ret_addr(4)
		pop ecx
		pop eax
		popfd
		ret
	}
}

__declspec(naked) void EmuFS_PushDwordPtrFs00()
{
	// push large dword ptr fs:0 — push KPCR[0] onto caller's stack
	// Uses two xchg instructions to shuffle the stack without any
	// off-stack storage, eliminating the static-variable data race.
	__asm
	{
		call EmuFS_RefreshKPCR
		// Stack: [ret_addr, <caller>]
		push eax                      // Stack: [eax_saved, ret_addr, <caller>]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov eax, [eax]                // eax = KPCR[0] (value to push)
		xchg eax, [esp + 4]          // Stack: [eax_saved, value, <caller>], eax = ret_addr
		xchg eax, [esp]              // Stack: [ret_addr, value, <caller>], eax = eax_saved
		ret                           // Stack: [value, <caller>] — push simulated
	}
}

__declspec(naked) void EmuFS_PopDwordPtrFs00()
{
	// pop large dword ptr fs:0 — pop caller's stack into KPCR[0]
	// Uses ret 4 to discard the consumed stack slot, matching the
	// ESP increment a real pop would have performed.
	__asm
	{
		call EmuFS_RefreshKPCR
		// Stack: [ret_addr, popped_value, <caller>]
		push eax
		push ebx
		// Stack: [ebx, eax, ret_addr, popped_value, <caller>]
		//         +0   +4   +8        +12
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov ebx, [esp + 12]           // ebx = popped_value
		mov[eax], ebx                 // KPCR[0] = popped_value
		pop ebx
		pop eax
		// Stack: [ret_addr, popped_value, <caller>]
		ret 4                         // return and discard popped_value slot
		// Stack: [<caller>] — pop simulated
	}
}

// ******************************************************************
// * Regular read thunks (trailing-offset-byte scheme)
// *
// * The thunk reads the KPCR offset byte from [return_address],
// * increments the return address past it, then loads KPCR[offset]
// * into the target register.
// ******************************************************************

__declspec(naked) void EmuFS_ReadEax()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push ecx
		mov ecx, [esp + 4]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 4]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov eax, [eax + ecx]
		pop ecx
		ret
	}
}

__declspec(naked) void EmuFS_ReadEbx()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push ecx
		mov ecx, [esp + 4]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 4]
		mov ebx, fs : [TIB_ArbitraryDataSlot]
		mov ebx, [ebx + ecx]
		pop ecx
		ret
	}
}

__declspec(naked) void EmuFS_ReadEcx()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		mov eax, [esp + 4]
		movzx eax, byte ptr[eax]
		inc dword ptr[esp + 4]
		mov ecx, fs : [TIB_ArbitraryDataSlot]
		mov ecx, [ecx + eax]
		pop eax
		ret
	}
}

__declspec(naked) void EmuFS_ReadEdx()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push ecx
		mov ecx, [esp + 4]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 4]
		mov edx, fs : [TIB_ArbitraryDataSlot]
		mov edx, [edx + ecx]
		pop ecx
		ret
	}
}

__declspec(naked) void EmuFS_ReadEsi()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push ecx
		mov ecx, [esp + 4]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 4]
		mov esi, fs : [TIB_ArbitraryDataSlot]
		mov esi, [esi + ecx]
		pop ecx
		ret
	}
}

__declspec(naked) void EmuFS_ReadEdi()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push ecx
		mov ecx, [esp + 4]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 4]
		mov edi, fs : [TIB_ArbitraryDataSlot]
		mov edi, [edi + ecx]
		pop ecx
		ret
	}
}

// ******************************************************************
// * Regular write thunks (trailing-offset-byte scheme)
// *
// * The thunk reads the KPCR offset byte from [return_address],
// * increments the return address past it, then stores the source
// * register into KPCR[offset]. Writes to offset 0x58 also sync
// * the global g_DpcRoutineActive flag.
// ******************************************************************

__declspec(naked) void EmuFS_WriteEax()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push ecx
		push ebx
		mov ecx, [esp + 8]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 8]
		mov ebx, fs : [TIB_ArbitraryDataSlot]
		mov[ebx + ecx], eax
		cmp ecx, 58h
		jne no_dpc_eax
		mov[g_DpcRoutineActive], eax
	no_dpc_eax:
		pop ebx
		pop ecx
		ret
	}
}

__declspec(naked) void EmuFS_WriteEbx()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		push ecx
		mov ecx, [esp + 8]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 8]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov[eax + ecx], ebx
		cmp ecx, 58h
		jne no_dpc_ebx
		mov[g_DpcRoutineActive], ebx
	no_dpc_ebx:
		pop ecx
		pop eax
		ret
	}
}

__declspec(naked) void EmuFS_WriteEcx()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		push ebx
		mov ebx, [esp + 8]
		movzx ebx, byte ptr[ebx]
		inc dword ptr[esp + 8]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov[eax + ebx], ecx
		cmp ebx, 58h
		jne no_dpc_ecx
		mov[g_DpcRoutineActive], ecx
	no_dpc_ecx:
		pop ebx
		pop eax
		ret
	}
}

__declspec(naked) void EmuFS_WriteEdx()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		push ecx
		mov ecx, [esp + 8]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 8]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov[eax + ecx], edx
		cmp ecx, 58h
		jne no_dpc_edx
		mov[g_DpcRoutineActive], edx
	no_dpc_edx:
		pop ecx
		pop eax
		ret
	}
}

__declspec(naked) void EmuFS_WriteEsi()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		push ecx
		mov ecx, [esp + 8]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 8]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov[eax + ecx], esi
		cmp ecx, 58h
		jne no_dpc_esi
		mov[g_DpcRoutineActive], esi
	no_dpc_esi:
		pop ecx
		pop eax
		ret
	}
}

__declspec(naked) void EmuFS_WriteEdi()
{
	__asm
	{
		call EmuFS_RefreshKPCR
		push eax
		push ecx
		mov ecx, [esp + 8]
		movzx ecx, byte ptr[ecx]
		inc dword ptr[esp + 8]
		mov eax, fs : [TIB_ArbitraryDataSlot]
		mov[eax + ecx], edi
		cmp ecx, 58h
		jne no_dpc_edi
		mov[g_DpcRoutineActive], edi
	no_dpc_edi:
		pop ecx
		pop eax
		ret
	}
}

// ******************************************************************
// * Patch table and matching
// ******************************************************************

struct fs_patch_t {
	std::vector<uint8_t> data;   // byte pattern to match
	void* functionPtr;           // thunk to redirect to
	int offsetBytePos;           // position of wildcard KPCR offset byte (-1 = exact match)
};

// initialize fs segment selector emulation
void EmuInitFS()
{
	// Build the patch table. Entries with offsetBytePos >= 0 use the
	// trailing-offset-byte scheme: the byte at that position is a wildcard
	// during matching (validated against kValidKpcrOffsets), and emitted as
	// the trailing byte after the CALL opcode. Entries are ordered by size
	// (longest first) to avoid false-positive partial matches.
	std::vector<fs_patch_t> patches;

	// 8-byte entries (exact match)
	patches.push_back({ { 0x64, 0x0F, 0xB6, 0x05, 0x24, 0x00, 0x00, 0x00 }, (void*)&EmuFS_MovzxEaxBytePtrFs24, -1 });  // movzx eax, byte ptr fs:[24h]

	// 7-byte entries — reads (wildcard offset at position 3)
	patches.push_back({ { 0x64, 0x8B, 0x1D, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_ReadEbx, 3 });  // mov ebx, fs:[??]
	patches.push_back({ { 0x64, 0x8B, 0x0D, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_ReadEcx, 3 });  // mov ecx, fs:[??]
	patches.push_back({ { 0x64, 0x8B, 0x15, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_ReadEdx, 3 });  // mov edx, fs:[??]
	patches.push_back({ { 0x64, 0x8B, 0x35, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_ReadEsi, 3 });  // mov esi, fs:[??]
	patches.push_back({ { 0x64, 0x8B, 0x3D, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_ReadEdi, 3 });  // mov edi, fs:[??]

	// 7-byte entries — writes (wildcard offset at position 3)
	patches.push_back({ { 0x64, 0x89, 0x1D, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEbx, 3 }); // mov fs:[??], ebx
	patches.push_back({ { 0x64, 0x89, 0x0D, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEcx, 3 }); // mov fs:[??], ecx
	patches.push_back({ { 0x64, 0x89, 0x15, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEdx, 3 }); // mov fs:[??], edx
	patches.push_back({ { 0x64, 0x89, 0x35, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEsi, 3 }); // mov fs:[??], esi
	patches.push_back({ { 0x64, 0x89, 0x3D, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEdi, 3 }); // mov fs:[??], edi
	patches.push_back({ { 0x64, 0x89, 0x25, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEsp, 3 }); // mov fs:[??], esp

	// 7-byte entries — specials (exact match)
	patches.push_back({ { 0x64, 0x3B, 0x35, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_CmpEsiFs00, -1 });      // cmp esi, fs:[0]
	patches.push_back({ { 0x64, 0x8F, 0x05, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_PopDwordPtrFs00, -1 });  // pop dword ptr fs:[0]
	patches.push_back({ { 0x64, 0xFF, 0x35, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_PushDwordPtrFs00, -1 }); // push dword ptr fs:[0]

	// 6-byte entries — eax short-form reads/writes (wildcard offset at position 2)
	patches.push_back({ { 0x64, 0xA1, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_ReadEax, 2 });   // mov eax, fs:[??]
	patches.push_back({ { 0x64, 0xA3, 0x00, 0x00, 0x00, 0x00 }, (void*)&EmuFS_WriteEax, 2 });  // mov fs:[??], eax

	EmuLogEx(CXBXR_MODULE::INIT, LOG_LEVEL::DEBUG, "Patching FS Register Accesses (%d patterns)\n", patches.size());
	long numberOfPatches = patches.size();

	// Iterate through each CODE section
	for (uint32_t sectionIndex = 0; sectionIndex < CxbxKrnl_Xbe->m_Header.dwSections; sectionIndex++) {
		if (!CxbxKrnl_Xbe->m_SectionHeader[sectionIndex].dwFlags.bExecutable) {
			continue;
		}

		EmuLogEx(CXBXR_MODULE::INIT, LOG_LEVEL::DEBUG, "Searching for FS instructions in section %s\n",
			CxbxKrnl_Xbe->m_szSectionName[sectionIndex]);
		xbox::addr_xt startAddr = CxbxKrnl_Xbe->m_SectionHeader[sectionIndex].dwVirtualAddr;
		xbox::addr_xt endAddr = startAddr + CxbxKrnl_Xbe->m_SectionHeader[sectionIndex].dwSizeofRaw;

		for (xbox::addr_xt addr = startAddr; addr < endAddr; addr++)
		{
			for (int i = 0; i < numberOfPatches; i++)
			{
				long sizeOfData = patches[i].data.size();
				if (addr + sizeOfData >= endAddr) {
					continue;
				}

				// Match the pattern, treating offsetBytePos as a wildcard
				bool match = true;
				uint8_t offsetByte = 0;

				for (int j = 0; j < sizeOfData; j++) {
					if (j == patches[i].offsetBytePos) {
						// Wildcard position — capture and validate against whitelist
						offsetByte = *(uint8_t*)(addr + j);
						bool valid = false;
						for (auto off : kValidKpcrOffsets) {
							if (offsetByte == off) { valid = true; break; }
						}
						if (!valid) { match = false; break; }
					} else {
						if (*(uint8_t*)(addr + j) != patches[i].data[j]) {
							match = false; break;
						}
					}
				}

				if (match) {
					EmuLogEx(CXBXR_MODULE::INIT, LOG_LEVEL::DEBUG, "Patching FS instruction at 0x%.8X\n", addr);

					// Write CALL rel32 (5 bytes)
					*(uint8_t*)addr = OPCODE_CALL_E8;
					*(uint32_t*)(addr + 1) = (uint32_t)patches[i].functionPtr - addr - 5;

					if (patches[i].offsetBytePos >= 0) {
						// Trailing-offset-byte scheme: write captured offset after CALL
						*(uint8_t*)(addr + 5) = offsetByte;
						// NOP-fill the remainder
						int remaining = sizeOfData - 6;
						if (remaining > 0) {
							memset((void*)(addr + 6), OPCODE_NOP_90, remaining);
						}
					} else {
						// Exact-match entry: NOP-fill after CALL
						int remaining = sizeOfData - 5;
						memset((void*)(addr + 5), OPCODE_NOP_90, remaining);
					}

					addr += sizeOfData - 1;
					break;
				}
			}
		}
	}

	EmuLogEx(CXBXR_MODULE::INIT, LOG_LEVEL::DEBUG, "Done patching FS Register Accesses\n");
}

// Get Xbox's TIB StackBase address from thread's StackBase.
xbox::PVOID EmuGetTIBStackBase(xbox::PVOID ThreadStackBase) {
	xbox::addr_xt StackBaseAddr = reinterpret_cast<xbox::addr_xt>(ThreadStackBase);
	StackBaseAddr -= sizeof(xbox::FX_SAVE_AREA);
	return reinterpret_cast<xbox::PVOID>(StackBaseAddr);
}

// generate stack size reserved for xbox threads to write on.
xbox::dword_xt EmuGenerateStackSize(xbox::addr_xt& espBaseAddress, xbox::ulong_xt TlsDataSize) {
	using namespace xbox;
	dword_xt StackSize = espBaseAddress & 15; // Fix 16 byte alignment
	espBaseAddress -= StackSize;
	StackSize += sizeof(FX_SAVE_AREA);
	TlsDataSize = ALIGN_UP(TlsDataSize, ulong_xt);
	StackSize += TlsDataSize; // (optional)
	StackSize += sizeof(KSTART_FRAME);
	StackSize += sizeof(KSWITCHFRAME);
	return StackSize;
}

// generate fs segment selector
template<bool IsHostThread>
void EmuGenerateFS(xbox::PETHREAD Ethread, unsigned Host2XbStackBaseReserved, unsigned Host2XbStackSizeReserved)
{
	xbox::PVOID base;
	xbox::ulong_xt size;

	// Allocate the xbox KPCR structure
	base = xbox::zeroptr;
	size = sizeof(xbox::KPCR);
	xbox::NtAllocateVirtualMemory(&base, 0, &size, XBOX_MEM_RESERVE | XBOX_MEM_COMMIT, XBOX_PAGE_READWRITE);
	xbox::KPCR *NewPcr = (xbox::KPCR*)base;
	xbox::RtlZeroMemory(NewPcr, sizeof(xbox::KPCR));
	xbox::NT_TIB *XbTib = &(NewPcr->NtTib);
	xbox::PKPRCB Prcb = &(NewPcr->PrcbData);
	// Note : As explained above (at EmuKeSetPcr), Cxbx cannot allocate one NT_TIB and KPRCB
	// structure per thread, since Cxbx currently doesn't do thread-switching.
	// Thus, the only way to give each thread it's own PrcbData.CurrentThread, is to put the
	// KPCR pointer in the TIB_ArbitraryDataSlot, which is read by the above EmuFS_* patches.
	//
	// Once we simulate thread switching ourselves, we can update PrcbData.CurrentThread
	// and simplify this initialization, by using only one KPCR for the single Xbox processor.
	//
	// One way to do our own (preemprive) thread-switching would be to use this technique :
	// http://www.eran.io/implementing-a-preemptive-kernel-within-a-single-windows-thread/
	// See https://github.com/Cxbx-Reloaded/Cxbx-Reloaded/issues/146 for more info.

	// Copy the Nt TIB over to the emulated TIB :
	NT_TIB* hTib = GetNtTib();
	{
		memcpy(XbTib, hTib, sizeof(NT_TIB));
		// Fixup the TIB self pointer :
		NewPcr->NtTib.Self = XbTib;

		// NOTE: The actual issue was TlsData was not within Host's stack which is now implemented.
		//       But instead of direct Host's stack, (which should not be tampered from Host's kernel stack block!)
		//       we allocated through inline asm to reserve xbox stack dynamically in order to have xbox's kernel stack reside in
		//       host's stack (in permitted function's stack usage).
		// Write the Xbox stack base to the Host, allows ConvertThreadToFiber to work correctly
		// Test case:
		// * DoA2
		// * DoA3
		// NOTE: This is disabled due to cause of corruption to host's TIB and
		//       silent crash for xbox threads creation.
		// Test case:
		// * Direct3DCreate9Ex call from inside xbox thread
		// * PCSTProxy (used from PsCreateSystemThreadEx export function)
		//__writefsdword(TIB_StackBase, (DWORD)NewPcr->NtTib.StackBase);
		//__writefsdword(TIB_StackLimit, (DWORD)NewPcr->NtTib.StackLimit);
	}

	// Set flat address of this PCR
	NewPcr->SelfPcr = NewPcr;
	// Set pointer to Prcb
	NewPcr->Prcb = Prcb;

	// Initialize the prcb :
	{
		// TODO : Once we do our own thread-switching (as mentioned above),
		// we can also start using Prcb->DpcListHead instead of DpcQueue :
		InitializeListHead(&(Prcb->DpcListHead));
		Prcb->DpcRoutineActive = FALSE;

		NewPcr->Irql = PASSIVE_LEVEL; // See KeLowerIrql;
	}

	if constexpr (IsHostThread) {
		// This only happens for the kernel initialization thread of cxbxr
		// Another thing to note, we do not insert into xbox's system as it will not be used in running xbox environment.
		// Instead it will be sleeping until title/user make a decision what to do next.
		assert(Ethread == xbox::zeroptr);

		base = xbox::zeroptr;
		size = sizeof(xbox::ETHREAD);
		xbox::NtAllocateVirtualMemory(&base, 0, &size, XBOX_MEM_RESERVE | XBOX_MEM_COMMIT, XBOX_PAGE_READWRITE);
		Ethread = (xbox::PETHREAD)base;
		xbox::RtlZeroMemory(Ethread, sizeof(xbox::ETHREAD)); // Clear, to prevent side-effects on random contents
		// Initialize the IRP tracking list for this thread
		InitializeListHead(&Ethread->IrpList);
		// Emulate kernel stack size as we can't use exact size.
		xbox::ulong_xt KernelStackSize = Host2XbStackBaseReserved - reinterpret_cast<xbox::ulong_xt>(hTib->StackLimit);
		// Since the cxbxr's kernel initialization occur there, we do not create a new thread
		// and therefore doesn't need to set any additional System/Start details set in the xbox's kernel stack.
		xbox::KeInitializeThread<IsHostThread>(
			&Ethread->Tcb,
			(xbox::PVOID)Host2XbStackBaseReserved,
			KernelStackSize,
			xbox::zero,
			xbox::zeroptr,  // Unused (SystemRoutine)
			xbox::zeroptr,  // Unused (StartRoutine)
			xbox::zeroptr,  // Unused (StartContext)
			xbox::zeroptr); // Unused (&KiUniqueProcess)
	}
#ifndef ENABLE_KTHREAD_SWITCHING
	else {
		// Otherwise, xbox::PsCreateSystemThreadEx is called and xbox::KeInitializeThread is already called from it.
		// But we need to carry the reserved part onto host's stack to able align with xbox and host sharing the same stack in a new thread.
		// Since we are using direct execution than in virtualization environment.
		// Tcb.StackBase always point at the beginning of kernel stack (DOWN).
		xbox::addr_xt xStackBase = reinterpret_cast<xbox::addr_xt>(Ethread->Tcb.StackBase);
		xbox::addr_xt xStackLimit = reinterpret_cast<xbox::addr_xt>(Ethread->Tcb.StackLimit);
		xbox::addr_xt xTlsData = reinterpret_cast<xbox::addr_xt>(Ethread->Tcb.TlsData);
		xbox::addr_xt xKernelStack = reinterpret_cast<xbox::addr_xt>(Ethread->Tcb.KernelStack);
		xbox::dword_xt xKernelStackSize = xStackBase - xKernelStack;
		assert(xKernelStackSize <= Host2XbStackSizeReserved);
		PVOID hKernelStack = reinterpret_cast<PVOID>(Host2XbStackBaseReserved - xKernelStackSize);
		std::memcpy(hKernelStack, Ethread->Tcb.KernelStack, xKernelStackSize);
		// Update TlsData address if used
		if (Ethread->Tcb.TlsData) {
			Ethread->Tcb.TlsData = reinterpret_cast<xbox::PVOID>(Host2XbStackBaseReserved - (xStackBase - xTlsData));
		}
		// Set stacks addresses
		Ethread->Tcb.StackBase = reinterpret_cast<xbox::PVOID>(Host2XbStackBaseReserved);
		Ethread->Tcb.StackLimit = hTib->StackLimit; // Always point to host's StackLimit.
		Ethread->Tcb.KernelStack = hKernelStack;
		// We can safely delete kernel stack as there is no virtualization environment implemented.
		xbox::MmDeleteKernelStack(reinterpret_cast<xbox::PVOID>(xStackBase), reinterpret_cast<xbox::PVOID>(xStackLimit));
	}
#endif

	// Initialize TlsData in order to avoid xbox's APIs attempt to try access null pointer from CxbxrCreateThread.
	// As far as I have seen, Xapi's CreateThread's startup function is the only one that does the tls' initialization.
	// It will repeat same process yet will not cause performance impact.
	// NOTE: PsCreateSystemThread's startup function does not do tls' initialization.
	if (Ethread->Tcb.TlsData) {
		Xbe::TLS* XbeTls = (Xbe::TLS*)CxbxKrnl_Xbe->m_Header.dwTLSAddr;
		uint32_t RawTlsDataSize = XbeTls->dwDataEndAddr - XbeTls->dwDataStartAddr;
		// First index is a pointer to the array of tls datas.
		xbox::addr_xt* TlsData = reinterpret_cast<xbox::addr_xt*>(Ethread->Tcb.TlsData);
		*TlsData = reinterpret_cast<xbox::addr_xt>(Ethread->Tcb.TlsData) + sizeof(xbox::addr_xt);
		// Set the actual tls data from xbe.
		TlsData += 1;
		std::memcpy(TlsData, reinterpret_cast<xbox::PVOID>(XbeTls->dwDataStartAddr), RawTlsDataSize);

		if (XbeTls->dwSizeofZeroFill) {
			std::memset(reinterpret_cast<xbox::PBYTE>(TlsData) + RawTlsDataSize, 0, XbeTls->dwSizeofZeroFill);
		}
	}

	// Initialize a fake PrcbData.CurrentThread 
	{
		// TODO: Do we need NtTib's overwrite in ENABLE_KTHREAD_SWITCHING usage?
		// Set the stack details over to NtTib's structure.
		NewPcr->NtTib.StackBase = EmuGetTIBStackBase(Ethread->Tcb.StackBase);
		NewPcr->NtTib.StackLimit = Ethread->Tcb.StackLimit;
		// Set PrcbData.CurrentThread
		Prcb->CurrentThread = (xbox::PKTHREAD)Ethread;
	}

	// Create a host wake event for this thread's dispatcher waits
	CxbxRegisterThreadWakeEvent((xbox::PKTHREAD)Ethread);

	// Make the KPCR struct available to EmuKeGetPcr()
	EmuKeSetPcr(NewPcr);

	EmuLog(LOG_LEVEL::DEBUG, "Installed KPCR in TIB_ArbitraryDataSlot (with Ethread->Tcb.TlsData = 0x%.8X)", Ethread->Tcb.TlsData);

	_controlfp(_PC_53, _MCW_PC); // Set Precision control to 53 bits (verified setting)
	_controlfp(_RC_NEAR, _MCW_RC); // Set Rounding control to near (unsure about this)
}
