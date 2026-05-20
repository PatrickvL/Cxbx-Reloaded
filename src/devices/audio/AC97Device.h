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

#ifndef _AC97_H_
#define _AC97_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include <AL/al.h>
#include <AL/alc.h>

#include "../PCIDevice.h"
class AC97Device : public PCIDevice {
	public:
		using PCIDevice::PCIDevice;

		// PCI Functions
		void Init();
		void Reset();

		uint32_t IORead(int barIndex, uint32_t addr, unsigned size = sizeof(uint8_t));
		void IOWrite(int barIndex, uint32_t addr, uint32_t data, unsigned size = sizeof(uint8_t));

		uint32_t MMIORead(int barIndex, uint32_t addr, unsigned size);
		void MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size);
		void Begin3DVoiceFrameBatch(size_t frameCount);
		void Submit3DVoiceFrames(uint32_t voiceHandle, uint32_t hrtfEntryIndex, bool stereo,
			const std::array<uint8_t, 4>& hrtfSubmix, uint8_t hrtfHeadroom,
			const int16_t* samples, size_t frameCount);
		void SubmitPCMFrames(const int16_t* samples, size_t frameCount);
	private:
		enum class PrimeResult : uint8_t {
			Ready,
			EndOfList,
			DescriptorError,
		};

		struct SpatialVoiceState {
			bool active = false;
			bool stereo = false;
			uint32_t hrtfEntryIndex = 0xFFFFFFFF;
			std::array<uint8_t, 4> hrtfSubmix{};
			uint8_t hrtfHeadroom = 0;
			std::vector<int16_t> samples{};
		};

		bool EnsureOutputDevice();
		void ResetOutputStream();
		uint32_t ReadRegister(uint32_t addr, unsigned size) const;
		void WriteRegister(uint32_t addr, uint32_t value, unsigned size);
		uint16_t ReadRegister16(uint32_t addr) const;
		void WriteRegister16(uint32_t addr, uint16_t value);
		void UpdateGlobalStatus();
		void UpdateBusMasterChannels();
		void ResetBusMasterChannel(uint32_t channelBase);
		PrimeResult PrimeBusMasterChannel(uint32_t channelBase);
		void UpdateBusMasterStatus(uint32_t channelBase);
		uint32_t GetBusMasterSampleRate(uint32_t channelBase) const;
		uint8_t GetPrefetchedIndexValue(uint32_t channelBase, uint8_t currentIndex, uint8_t lastValidIndex) const;
		bool ReadGuest32(uint32_t guestAddress, uint32_t& value) const;
		bool IsDescriptorErrorAcknowledged(uint32_t channelBase) const;

		std::array<uint8_t, 0x180> m_Registers{};
		std::array<uint32_t, 3> m_ChannelLastUpdate{};
		std::array<uint32_t, 3> m_ChannelSampleRemainder{};
		std::array<bool, 3> m_ChannelAdvanceOnRestart{};
		std::array<bool, 3> m_ChannelQueuedAfterHalt{};
		std::array<bool, 3> m_ChannelDescriptorError{};
		std::vector<int16_t> m_OutputScratch{};
		std::unordered_map<uint32_t, SpatialVoiceState> m_Pending3DVoices{};
		std::vector<ALuint> m_FreeOutputBuffers{};
		std::unordered_map<ALuint, size_t> m_OutputBufferIndex{};
		ALCdevice* m_OutputDevice = nullptr;
		ALCcontext* m_OutputContext = nullptr;
		ALuint m_OutputSource = 0;
		std::array<ALuint, 16> m_OutputBuffers{};
		std::array<uint32_t, 16> m_OutputBufferBytes{};
		uint32_t m_QueuedAudioBytes = 0;
		bool m_OutputDeviceFailed = false;
		bool m_LoggedQueueFull = false;
};

#endif
