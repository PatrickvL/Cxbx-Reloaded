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
#include <vector>

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
		void SubmitPCMFrames(const int16_t* samples, size_t frameCount);
	private:
		bool EnsureOutputDevice();
		uint32_t ReadRegister(uint32_t addr, unsigned size) const;
		void WriteRegister(uint32_t addr, uint32_t value, unsigned size);
		uint16_t ReadRegister16(uint32_t addr) const;
		void WriteRegister16(uint32_t addr, uint16_t value);
		void UpdateBusMasterChannels();
		void ResetBusMasterChannel(uint32_t channelBase);
		bool PrimeBusMasterChannel(uint32_t channelBase);
		void UpdateBusMasterStatus(uint32_t channelBase);
		uint32_t GetBusMasterSampleRate(uint32_t channelBase) const;
		bool ReadGuest32(uint32_t guestAddress, uint32_t& value) const;

		std::array<uint8_t, 0x180> m_Registers{};
		std::array<uint32_t, 3> m_ChannelLastUpdate{};
		std::array<uint32_t, 3> m_ChannelSampleRemainder{};
		std::vector<int16_t> m_OutputScratch{};
		uint32_t m_OutputDevice = 0;
		bool m_OutputDeviceFailed = false;
		bool m_LoggedQueueFull = false;
};

#endif
