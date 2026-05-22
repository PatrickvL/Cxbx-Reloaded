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
// *  (c) 2018 Luke Usher <luke.usher@outlook.coM>
// *
// *  All rights reserved
// *
// ******************************************************************

#include <cstdio>

#include "APUDevice.h"
#include "../../common/AddressRanges.h"

#include <Windows.h>

extern uint32_t GetAPUTime();

// TODO: Everything :P
// TODO: Audio Processing/Thread

#define APU_VP_BASE 0x20000
#define APU_VP_SIZE 0x10000

#define APU_GP_BASE 0x30000
#define APU_GP_SIZE 0x10000

#define APU_EP_BASE 0x50000
#define APU_EP_SIZE 0x10000

void APUDevice::Init()
{
	PCIBarRegister r;
	r.Raw.type = PCI_BAR_TYPE_IO;
	r.IO.address = 0xD000 >> 4;
	RegisterBAR(0, 256, r.value);

	r.Raw.type = PCI_BAR_TYPE_IO;
	r.IO.address = 0xD200 >> 4;
	RegisterBAR(0, 128, r.value);

	r.Raw.type = PCI_BAR_TYPE_MEMORY;
	r.Memory.address = APU_BASE >> 4;
	RegisterBAR(2, APU_SIZE, r.value);

	m_DeviceId = 0x01B0;
	m_VendorId = PCI_VENDOR_ID_NVIDIA;
}
	
void APUDevice::Reset()
{

}

uint32_t APUDevice::IORead(int barIndex, uint32_t addr, unsigned size)
{
	return 0;
}

void APUDevice::IOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
}

uint32_t APUDevice::MMIORead(int barIndex, uint32_t addr, unsigned size)
{
	// DIAG: Log all APU reads
	static int readCount = 0;
	if (readCount++ < 50) {
		fprintf(stderr, "[APU-MMIO] Read addr=0x%05X size=%u\n", addr, size);
	}

	if (addr >= APU_VP_BASE && addr < APU_VP_BASE + APU_VP_SIZE) {
		return VPRead(addr - APU_VP_BASE, size);
	}

	if (addr >= APU_GP_BASE && addr < APU_GP_BASE + APU_GP_SIZE) {
		return GPRead(addr - APU_GP_BASE, size);
	}

	if (addr >= APU_EP_BASE && addr < APU_EP_BASE + APU_EP_SIZE) {
		return EPRead(addr - APU_EP_BASE, size);
	}

	switch (addr) {
		case 0x200C: return GetAPUTime();	
	}

	return 0;
}

void APUDevice::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	// DIAG: Log all APU writes
	static int writeCount = 0;
	if (writeCount++ < 100) {
		fprintf(stderr, "[APU-MMIO] Write addr=0x%05X value=0x%08X size=%u\n", addr, value, size);
	}

	if (addr >= APU_VP_BASE && addr < APU_VP_BASE + APU_VP_SIZE) {
		VPWrite(addr - APU_VP_BASE, value, size);
		return;
	}

	if (addr >= APU_GP_BASE && addr < APU_GP_BASE + APU_GP_SIZE) {
		GPWrite(addr - APU_GP_BASE, value, size);
		return;
	}

	if (addr >= APU_EP_BASE && addr < APU_EP_BASE + APU_EP_SIZE) {
		EPWrite(addr - APU_EP_BASE, value, size);
		return;
	}

	// Unhandled APU register write
}


uint32_t APUDevice::GPRead(uint32_t addr, unsigned size)
{
	return 0;
}

void APUDevice::GPWrite(uint32_t addr, uint32_t value, unsigned size)
{
}


uint32_t APUDevice::VPRead(uint32_t addr, unsigned size)
{
	switch (addr) {
		case 0x10: return 0x80; // HACK: Pretend the FIFO is always empty, bypasses hangs when APU isn't fully implemented
		case NV_PAPU_VPVADDR_OFF: return m_vpvaddr;
	}

	return 0;
}

void APUDevice::VPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	switch (addr) {
		case NV_PAPU_VPVADDR_OFF:
			m_vpvaddr = value;
			fprintf(stderr, "[APU] VPVADDR set to 0x%08X (VA=0x%08X)\n", value, CONTIGUOUS_MEMORY_BASE + value);
			return;
	}
}

void APUDevice::AdvanceVoiceCursors()
{
	if (m_vpvaddr == 0) {
		static int logCount = 0;
		if (logCount++ < 5) fprintf(stderr, "[APU] AdvanceVoiceCursors called but VPVADDR=0\n");
		return;
	}

	DWORD now = GetTickCount();
	if (m_lastTickMs == 0) { m_lastTickMs = now; return; }
	DWORD elapsed = now - m_lastTickMs;
	if (elapsed == 0) return;
	m_lastTickMs = now;

	// Voice descriptor table lives in contiguous (physical) memory.
	// In Cxbx, CONTIGUOUS_MEMORY_BASE (0x80000000) + physical offset is the virtual address.
	uint8_t* voiceTable = (uint8_t*)(CONTIGUOUS_MEMORY_BASE + m_vpvaddr);

	for (int v = 0; v < NV_PAVS_MAX_VOICES; v++) {
		uint8_t* vd = voiceTable + v * NV_PAVS_VOICE_SIZE;
		volatile uint32_t* pState = (volatile uint32_t*)(vd + NV_PAVS_VOICE_PAR_STATE_OFF);
		volatile uint32_t* pOffset = (volatile uint32_t*)(vd + NV_PAVS_VOICE_PAR_OFFSET_OFF);
		volatile uint32_t* pFmt = (volatile uint32_t*)(vd + NV_PAVS_VOICE_CFG_FMT_OFF);

		// Check if voice has valid format (non-zero means initialized)
		uint32_t fmt = *pFmt;
		if (fmt == 0) continue;

		// Read current CBO
		uint32_t offset_reg = *pOffset;
		uint32_t cbo = offset_reg & NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK;

		// Derive sample rate from format register.
		// NV_PAVS_VOICE_CFG_FMT contains sample rate and format info.
		// For a rough approximation, assume 48000 Hz stereo 16-bit (most common Xbox format).
		// bytes_per_ms = sampleRate * channels * bytesPerSample / 1000
		// = 48000 * 2 * 2 / 1000 = 192 bytes/ms
		uint32_t bytesPerMs = 192;
		uint32_t advance = elapsed * bytesPerMs;

		// Advance CBO
		uint32_t newCbo = cbo + advance;

		// Write back (preserve upper bits of the register)
		*pOffset = (offset_reg & ~NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK) | (newCbo & NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK);
	}
}


uint32_t APUDevice::EPRead(uint32_t addr, unsigned size)
{
	return 0;
}

void APUDevice::EPWrite(uint32_t addr, uint32_t value, unsigned size)
{
}
