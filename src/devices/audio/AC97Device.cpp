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

#include "AC97Device.h"

#include <cstring>

namespace {

constexpr uint32_t AC97_NAM_SIZE = 0x100;
constexpr uint32_t AC97_NABM_SIZE = 0x80;
constexpr uint32_t AC97_MMIO_SIZE = AC97_NAM_SIZE + AC97_NABM_SIZE;

constexpr uint32_t AC97_Reset = 0x00;
constexpr uint32_t AC97_Powerdown_Ctrl_Stat = 0x26;
constexpr uint32_t AC97_Extended_Audio_ID = 0x28;
constexpr uint32_t AC97_Extended_Audio_Ctrl_Stat = 0x2A;
constexpr uint32_t AC97_PCM_Front_DAC_Rate = 0x2C;
constexpr uint32_t AC97_PCM_Surround_DAC_Rate = 0x2E;
constexpr uint32_t AC97_PCM_LFE_DAC_Rate = 0x30;
constexpr uint32_t AC97_PCM_LR_ADC_Rate = 0x32;
constexpr uint32_t AC97_MIC_ADC_Rate = 0x34;
constexpr uint32_t AC97_Vendor_ID1 = 0x7C;
constexpr uint32_t AC97_Vendor_ID2 = 0x7E;

constexpr uint32_t NABM_PI_BASE = 0x00;
constexpr uint32_t NABM_PO_BASE = 0x10;
constexpr uint32_t NABM_MC_BASE = 0x20;
constexpr uint32_t NABM_GLOB_CNT = 0x2C;
constexpr uint32_t NABM_GLOB_STA = 0x30;

constexpr uint32_t BM_CIV = 0x04;
constexpr uint32_t BM_LVI = 0x05;
constexpr uint32_t BM_SR = 0x06;
constexpr uint32_t BM_PICB = 0x08;
constexpr uint32_t BM_PIV = 0x0A;
constexpr uint32_t BM_CR = 0x0B;

constexpr uint16_t SR_FIFOE = 1 << 4;
constexpr uint16_t SR_BCIS = 1 << 3;
constexpr uint16_t SR_LVBCI = 1 << 2;
constexpr uint16_t SR_CELV = 1 << 1;
constexpr uint16_t SR_DCH = 1 << 0;
constexpr uint16_t SR_WCLEAR_MASK = SR_FIFOE | SR_BCIS | SR_LVBCI;

constexpr uint8_t CR_RR = 1 << 1;
constexpr uint8_t CR_RPBM = 1 << 0;

constexpr uint16_t AC97_EXT_AUDIO_ID_VRA = 1 << 0;
constexpr uint16_t AC97_EXT_AUDIO_ID_VRM = 1 << 3;
constexpr uint16_t AC97_POWER_READY = 0x000F;
constexpr uint16_t AC97_RATE_48KHZ = 48000;
constexpr uint16_t AC97_VENDOR_SIGMATEL_1 = 0x8384;
constexpr uint16_t AC97_VENDOR_SIGMATEL_2 = 0x7608;

uint32_t ReadLE(const uint8_t* data, uint32_t addr, unsigned size)
{
	uint32_t value = 0;
	for (unsigned i = 0; i < size; ++i) {
		value |= static_cast<uint32_t>(data[addr + i]) << (i * 8);
	}
	return value;
}

void WriteLE(uint8_t* data, uint32_t addr, uint32_t value, unsigned size)
{
	for (unsigned i = 0; i < size; ++i) {
		data[addr + i] = static_cast<uint8_t>((value >> (i * 8)) & 0xFF);
	}
}

}

void AC97Device::Init()
{
	PCIBarRegister r;

	r.Raw.type = PCI_BAR_TYPE_IO;
	r.IO.address = 0xD000;
	RegisterBAR(0, AC97_NAM_SIZE, r.value);

	r.Raw.type = PCI_BAR_TYPE_IO;
	r.IO.address = 0xD200;
	RegisterBAR(1, AC97_NABM_SIZE, r.value);

	r.Raw.type = PCI_BAR_TYPE_MEMORY;
	r.Memory.address = AC97_BASE >> 4;
	RegisterBAR(2, AC97_SIZE, r.value);

	m_DeviceId = 0x01B1;
	m_VendorId = PCI_VENDOR_ID_NVIDIA;

	Reset();
}

void AC97Device::Reset()
{
	std::memset(m_Registers.data(), 0, m_Registers.size());

	WriteRegister16(AC97_Powerdown_Ctrl_Stat, AC97_POWER_READY);
	WriteRegister16(AC97_Extended_Audio_ID, AC97_EXT_AUDIO_ID_VRA | AC97_EXT_AUDIO_ID_VRM);
	WriteRegister16(AC97_Extended_Audio_Ctrl_Stat, AC97_EXT_AUDIO_ID_VRA);
	WriteRegister16(AC97_PCM_Front_DAC_Rate, AC97_RATE_48KHZ);
	WriteRegister16(AC97_PCM_Surround_DAC_Rate, AC97_RATE_48KHZ);
	WriteRegister16(AC97_PCM_LFE_DAC_Rate, AC97_RATE_48KHZ);
	WriteRegister16(AC97_PCM_LR_ADC_Rate, AC97_RATE_48KHZ);
	WriteRegister16(AC97_MIC_ADC_Rate, AC97_RATE_48KHZ);
	WriteRegister16(AC97_Vendor_ID1, AC97_VENDOR_SIGMATEL_1);
	WriteRegister16(AC97_Vendor_ID2, AC97_VENDOR_SIGMATEL_2);

	ResetBusMasterChannel(NABM_PI_BASE);
	ResetBusMasterChannel(NABM_PO_BASE);
	ResetBusMasterChannel(NABM_MC_BASE);
	WriteRegister(AC97_NAM_SIZE + NABM_GLOB_CNT, 0, sizeof(uint32_t));
	WriteRegister(AC97_NAM_SIZE + NABM_GLOB_STA, 0, sizeof(uint32_t));
}

uint32_t AC97Device::IORead(int barIndex, uint32_t addr, unsigned size)
{
	switch (barIndex) {
	case 0:
		return ReadRegister(addr, size);
	case 1:
		return ReadRegister(AC97_NAM_SIZE + addr, size);
	default:
		return 0;
	}
}

void AC97Device::IOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	switch (barIndex) {
	case 0:
		if (addr == AC97_Reset && size >= sizeof(uint16_t)) {
			Reset();
			return;
		}
		WriteRegister(addr, value, size);
		return;
	case 1:
		addr += AC97_NAM_SIZE;
		if (addr == AC97_NAM_SIZE + NABM_PI_BASE + BM_SR ||
			addr == AC97_NAM_SIZE + NABM_PO_BASE + BM_SR ||
			addr == AC97_NAM_SIZE + NABM_MC_BASE + BM_SR) {
			const uint16_t current = ReadRegister16(addr);
			WriteRegister16(addr, current & ~(static_cast<uint16_t>(value) & SR_WCLEAR_MASK));
			return;
		}
		if (addr == AC97_NAM_SIZE + NABM_PI_BASE + BM_CR ||
			addr == AC97_NAM_SIZE + NABM_PO_BASE + BM_CR ||
			addr == AC97_NAM_SIZE + NABM_MC_BASE + BM_CR) {
			WriteRegister(addr, value, size);
			const uint32_t channelBase = (addr - AC97_NAM_SIZE) & ~0xF;
			if (value & CR_RR) {
				ResetBusMasterChannel(channelBase);
			} else {
				UpdateBusMasterStatus(channelBase);
			}
			return;
		}
		WriteRegister(addr, value, size);
		return;
	default:
		return;
	}
}

uint32_t AC97Device::MMIORead(int barIndex, uint32_t addr, unsigned size)
{
	(void)barIndex;

	if (addr < AC97_MMIO_SIZE) {
		return ReadRegister(addr, size);
	}

	return 0;
}

void AC97Device::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	(void)barIndex;

	if (addr < AC97_MMIO_SIZE) {
		IOWrite(addr < AC97_NAM_SIZE ? 0 : 1, addr < AC97_NAM_SIZE ? addr : addr - AC97_NAM_SIZE, value, size);
	}
}

uint32_t AC97Device::ReadRegister(uint32_t addr, unsigned size) const
{
	if (size == 0 || addr + size > m_Registers.size()) {
		return 0;
	}

	return ReadLE(m_Registers.data(), addr, size);
}

void AC97Device::WriteRegister(uint32_t addr, uint32_t value, unsigned size)
{
	if (size == 0 || addr + size > m_Registers.size()) {
		return;
	}

	WriteLE(m_Registers.data(), addr, value, size);
}

uint16_t AC97Device::ReadRegister16(uint32_t addr) const
{
	return static_cast<uint16_t>(ReadRegister(addr, sizeof(uint16_t)));
}

void AC97Device::WriteRegister16(uint32_t addr, uint16_t value)
{
	WriteRegister(addr, value, sizeof(uint16_t));
}

void AC97Device::ResetBusMasterChannel(uint32_t channelBase)
{
	WriteRegister(AC97_NAM_SIZE + channelBase + BM_CIV, 0, sizeof(uint8_t));
	WriteRegister(AC97_NAM_SIZE + channelBase + BM_LVI, 0, sizeof(uint8_t));
	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_SR, SR_DCH);
	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_PICB, 0);
	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_PIV, 0);
	WriteRegister(AC97_NAM_SIZE + channelBase + BM_CR, 0, sizeof(uint8_t));
}

void AC97Device::UpdateBusMasterStatus(uint32_t channelBase)
{
	const uint32_t crAddr = AC97_NAM_SIZE + channelBase + BM_CR;
	const uint32_t srAddr = AC97_NAM_SIZE + channelBase + BM_SR;
	const uint8_t control = static_cast<uint8_t>(ReadRegister(crAddr, sizeof(uint8_t)));
	uint16_t status = ReadRegister16(srAddr) & SR_WCLEAR_MASK;
	if (control & CR_RPBM) {
		status &= ~SR_DCH;
	} else {
		status |= SR_DCH;
	}

	if (ReadRegister(AC97_NAM_SIZE + channelBase + BM_LVI, sizeof(uint8_t)) == 0) {
		status |= SR_CELV;
	} else {
		status &= ~SR_CELV;
	}

	WriteRegister16(srAddr, status);
}
