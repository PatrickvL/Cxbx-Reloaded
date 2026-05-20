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

#include <array>
#include <cstddef>
#include <cstdint>

#include "../PCIDevice.h"
class APUDevice : public PCIDevice {
public:
	using PCIDevice::PCIDevice;

	static constexpr size_t MAX_VOICE_HANDLES = 0xFFFF;

	// PCI Functions
	void Init();
	void Reset();

	uint32_t IORead(int barIndex, uint32_t addr, unsigned size = sizeof(uint8_t));
	void IOWrite(int barIndex, uint32_t addr, uint32_t data, unsigned size = sizeof(uint8_t));

	uint32_t MMIORead(int barIndex, uint32_t addr, unsigned size);
	void MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size);
	void SynchronizeAudio();
private:
	struct SSLData {
		uint32_t base[2]{};
		uint8_t count[2]{};
		uint32_t ssl_index = 0;
		uint32_t ssl_seg = 0;
	};

	struct PlaybackState {
		uint32_t offset = 0;
		double fraction = 0.0;
		bool valid = false;
	};

	struct LowPassFilterState {
		float high = 0.0f;
		float band = 0.0f;
		float low = 0.0f;
	};

	uint32_t GPRead(uint32_t addr, unsigned size);
	void GPWrite(uint32_t addr, uint32_t value, unsigned size);
	uint32_t EPRead(uint32_t addr, unsigned size);
	void EPWrite(uint32_t addr, uint32_t value, unsigned size);
	uint32_t VPRead(uint32_t addr, unsigned size);
	void VPWrite(uint32_t addr, uint32_t value, unsigned size);
	void ConsumeVPMethod(uint32_t addr, uint32_t value, unsigned size);
	void UpdateVPFifo();
	void RefreshVPStatus();
	void RefreshInterruptStatus();
	void RenderBasicAudioChunk(size_t frameCount);
	void WriteOutputBuffers(const int32_t* mixBins, size_t frameCount);
	void RenderBasicVoiceList(uint32_t topRegister, int32_t* mixBins, size_t frameCount);
	void RenderBasicVoice(uint32_t voiceHandle, int32_t* mixBins, size_t frameCount);
	void InitializeVoiceEnvelopes(uint32_t voiceHandle, uint32_t voiceOnValue);
	void BeginVoiceRelease(uint32_t voiceHandle);
	float StepVoiceEnvelope(uint32_t voiceHandle, uint32_t reg0, uint32_t regA,
		uint32_t rrReg, uint32_t rrMask, uint32_t levelRegister, uint32_t levelMask,
		uint32_t countMask, uint32_t stateMask);
	bool ReadGuestWord(uint32_t guestAddress, uint32_t& value) const;
	bool ReadGuestBytes(uint32_t guestAddress, void* dest, size_t size) const;
	bool WriteGuestBytes(uint32_t guestAddress, const void* src, size_t size);
	bool WriteGuestWord(uint32_t guestAddress, uint32_t value);
	bool WriteGuestWordMasked(uint32_t guestAddress, uint32_t mask, uint32_t value);
	bool ReadVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t& value) const;
	bool WriteVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t value);
	bool WriteVPScatterGatherEntry(uint32_t handle, uint32_t value);
	void WriteNotifierStatus(uint32_t voiceHandle, uint32_t notifier, uint8_t status);
	bool ResolveVoiceAddress(uint32_t linearAddress, uint32_t& guestAddress) const;
	bool ReadVoiceBufferBytes(uint32_t linearAddress, void* dest, size_t size) const;
	bool WriteGuestCircularBuffer(uint32_t guestAddress, uint32_t length, uint32_t& cursor,
		const void* src, size_t size);

	uint32_t ReadRegister(uint32_t addr, unsigned size) const;
	void WriteRegister(uint32_t addr, uint32_t value, unsigned size);
	void SetRegister32(uint32_t addr, uint32_t value);
	uint32_t GetRegister32(uint32_t addr) const;

	std::array<uint8_t, APU_SIZE> m_Registers{};
	uint32_t m_VPFifoLevel = 0;
	uint32_t m_VPFifoLastUpdate = 0;
	uint32_t m_LastAudioUpdate = 0;
	uint32_t m_VPInputSgeHandle = 0;
	uint32_t m_VPOutputSgeHandle = 0;
	uint32_t m_VPSSLBasePage = 0;
	std::array<uint32_t, 4> m_VPOutBufferCursor{};
	std::array<SSLData, MAX_VOICE_HANDLES> m_VPSSLData{};
	std::array<PlaybackState, MAX_VOICE_HANDLES> m_VPPlaybackState{};
	std::array<std::array<LowPassFilterState, 2>, MAX_VOICE_HANDLES> m_VPLowPassState{};
	bool m_LoggedXADPCMDecodeFailure = false;
};

#endif
