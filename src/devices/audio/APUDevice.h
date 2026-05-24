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
// *  (c) 2018 Luke Usher <luke.usher@outlook.com>
// *
// *  All rights reserved
// *
// ******************************************************************

#ifndef _APU_H_
#define _APU_H_

#include "../PCIDevice.h"
#include <vector>
#include <thread>
#include <atomic>

// VP register offsets (relative to VP base 0x20000)
#define NV_PAPU_VPVADDR_OFF   0x002C  // Voice descriptor table physical address
#define NV_PAPU_VPSGEADDR_OFF 0x0030  // SGE table physical address

// Voice descriptor field offsets
#define NV_PAVS_VOICE_CFG_FMT_OFF      0x04
#define NV_PAVS_VOICE_PAR_STATE_OFF    0x54
#define NV_PAVS_VOICE_PAR_OFFSET_OFF   0x58  // Contains CBO (Current Buffer Offset)

// Voice descriptor size
#define NV_PAVS_VOICE_SIZE    0x80
#define NV_PAVS_MAX_VOICES    256

// Masks
#define NV_PAVS_VOICE_PAR_OFFSET_CBO_MASK 0x00FFFFFF

class APUDevice : public PCIDevice {
public:
	using PCIDevice::PCIDevice;

	// PCI Functions
	void Init();
	void Reset();

	uint32_t IORead(int barIndex, uint32_t addr, unsigned size = sizeof(uint8_t));
	void IOWrite(int barIndex, uint32_t addr, uint32_t data, unsigned size = sizeof(uint8_t));

	uint32_t MMIORead(int barIndex, uint32_t addr, unsigned size);
	void MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size);

	// Voice table base address (physical) — set by game writing NV_PAPU_VPVADDR
	uint32_t GetVPVADDR() const { return m_vpvaddr; }

	// Fallback voice descriptor base for games that bypass CMcpxAPU::Initialize
	// and write voice descriptors directly to contiguous memory (VPVADDR = 0).
	void SetFallbackVoiceBase(uint8_t* voiceBase) {
		// Add to list of tracked voices (don't add duplicates)
		for (auto* v : m_fallbackVoices) { if (v == voiceBase) return; }
		m_fallbackVoices.push_back(voiceBase);
	}

	// Advance CBO for all active voices based on elapsed time.
	// Called periodically from the APU worker thread.
	void AdvanceVoiceCursors();

	// Start the APU voice processing thread (like pfifo_puller for NV2A)
	void StartVoiceProcessingThread();

	// Stop the APU voice processing thread
	void StopVoiceProcessingThread();

private:
	uint32_t GPRead(uint32_t addr, unsigned size);
	void GPWrite(uint32_t addr, uint32_t value, unsigned size);
	uint32_t EPRead(uint32_t addr, unsigned size);
	void EPWrite(uint32_t addr, uint32_t value, unsigned size);
	uint32_t VPRead(uint32_t addr, unsigned size);
	void VPWrite(uint32_t addr, uint32_t value, unsigned size);

	uint32_t m_vpvaddr = 0;           // Voice descriptor table physical address
	std::vector<uint8_t*> m_fallbackVoices; // Active voices when VPVADDR=0
	uint32_t m_lastTickMs = 0;        // Last tick time for CBO advancement
	std::atomic<bool> m_exit{false};  // Thread stop signal
	std::thread m_thread;              // Voice processing thread
};

#endif
