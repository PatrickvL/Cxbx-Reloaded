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
#include <vector>
#include <algorithm>

#include "APUDevice.h"
#include "common\util\CxbxUtil.h"

extern uint32_t GetAPUTime();

// Cache of detected voice descriptor bases for the fallback path (VPVADDR=0).
static std::vector<uint8_t*> s_fallbackVoices;
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

	StartVoiceProcessingThread();
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
    DWORD now = GetTickCount();
    if (m_lastTickMs == 0) { m_lastTickMs = now; return; }
    DWORD elapsedMs = now - m_lastTickMs;
    if (elapsedMs == 0) return;
    m_lastTickMs = now;

    auto decodeBytesPerMs = [](uint32_t fmt) -> uint32_t {
        uint32_t codec  =  fmt & 0xF;
        uint32_t chans  = (fmt >> 4) & 0xF;
        uint32_t rate   = (fmt >> 8) & 0xFFF;
        if (chans == 0) chans = 2;
        if (rate == 0)  rate = 48000;
        return (rate * chans * ((codec == 0) ? 1 : 2)) / 1000;
    };

    // Full VPVADDR mode: process all 256 voices.
    if (m_vpvaddr != 0) {
        uint8_t* table = (uint8_t*)(CONTIGUOUS_MEMORY_BASE + m_vpvaddr);
        for (int v = 0; v < NV_PAVS_MAX_VOICES; v++) {
            uint8_t* vd = table + v * NV_PAVS_VOICE_SIZE;
            volatile uint32_t* pFmt    = (volatile uint32_t*)(vd + NV_PAVS_VOICE_CFG_FMT_OFF);
            volatile uint32_t* pOffset = (volatile uint32_t*)(vd + NV_PAVS_VOICE_PAR_OFFSET_OFF);
            uint32_t fmt = *pFmt;
            if (fmt == 0) continue;
            uint32_t offReg = *pOffset;
            uint32_t cbo = offReg & NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK;
            uint32_t newCbo = cbo + elapsedMs * decodeBytesPerMs(fmt);
            *pOffset = (offReg & ~NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK)
                     | (newCbo & NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK);
        }
        return;
    }

    // Advance CBO for all registered fallback voices.
    for (uint8_t* vd : m_fallbackVoices) {
        volatile uint32_t* pFmt    = (volatile uint32_t*)(vd + NV_PAVS_VOICE_CFG_FMT_OFF);
        volatile uint32_t* pOffset = (volatile uint32_t*)(vd + NV_PAVS_VOICE_PAR_OFFSET_OFF);
        uint32_t fmt = *pFmt;
        uint32_t bytesPerMs = (fmt != 0) ? decodeBytesPerMs(fmt) : 192;
        uint32_t offReg = *pOffset;
        uint32_t cbo = offReg & NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK;
        uint32_t newCbo = cbo + elapsedMs * bytesPerMs;
        *pOffset = (offReg & ~NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK)
                 | (newCbo & NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK);
    }
}

void APUDevice::StartVoiceProcessingThread()
{
    m_exit.store(false);
    m_thread = std::thread([this]() {
        while (!m_exit.load()) {
            AdvanceVoiceCursors();
            Sleep(1); // ~1ms polling interval, matching real APU tick granularity
        }
    });
}

void APUDevice::StopVoiceProcessingThread()
{
    m_exit.store(true);
    if (m_thread.joinable()) {
        m_thread.join();
    }
}


uint32_t APUDevice::EPRead(uint32_t addr, unsigned size)
{
	return 0;
}

void APUDevice::EPWrite(uint32_t addr, uint32_t value, unsigned size)
{
}
