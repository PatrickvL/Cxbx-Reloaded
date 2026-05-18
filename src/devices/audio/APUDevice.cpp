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

#include "APUDevice.h"
#include "APUTimer.h"
#include "common/AddressRanges.h"

#include <cstring>

namespace {

constexpr uint32_t APU_VP_BASE = 0x20000;
constexpr uint32_t APU_VP_SIZE = 0x10000;
constexpr uint32_t APU_VP_FREE = 0x10;
constexpr uint32_t APU_VP_FIFO_CAPACITY = 0x80;
constexpr uint32_t APU_VP_STATUS_EMPTY = 0x80;
constexpr uint32_t APU_VP_VOICE_MAX_HANDLE = 0xFFFF;

constexpr uint32_t APU_GP_BASE = 0x30000;
constexpr uint32_t APU_GP_SIZE = 0x10000;

constexpr uint32_t APU_EP_BASE = 0x50000;
constexpr uint32_t APU_EP_SIZE = 0x10000;

constexpr uint32_t NV_PAPU_ISTS = 0x00001000;
constexpr uint32_t NV_PAPU_IEN = 0x00001004;
constexpr uint32_t NV_PAPU_FECTL = 0x00001100;
constexpr uint32_t NV_PAPU_FECV = 0x00001110;
constexpr uint32_t NV_PAPU_FEAV = 0x00001118;
constexpr uint32_t NV_PAPU_FENADDR = 0x0000115C;
constexpr uint32_t NV_PAPU_FEDECMETH = 0x00001300;
constexpr uint32_t NV_PAPU_FEDECPARAM = 0x00001304;
constexpr uint32_t NV_PAPU_FEMEMADDR = 0x00001324;
constexpr uint32_t NV_PAPU_FEMEMDATA = 0x00001334;
constexpr uint32_t NV_PAPU_FETFORCE0 = 0x00001500;
constexpr uint32_t NV_PAPU_FETFORCE1 = 0x00001504;
constexpr uint32_t NV_PAPU_SECTL = 0x00002000;
constexpr uint32_t NV_PAPU_XGSCNT = 0x0000200C;
constexpr uint32_t NV_PAPU_VPVADDR = 0x0000202C;
constexpr uint32_t NV_PAPU_VPSGEADDR = 0x00002030;
constexpr uint32_t NV_PAPU_VPSSLADDR = 0x00002034;
constexpr uint32_t NV_PAPU_GPSADDR = 0x00002040;
constexpr uint32_t NV_PAPU_GPFADDR = 0x00002044;
constexpr uint32_t NV_PAPU_EPSADDR = 0x00002048;
constexpr uint32_t NV_PAPU_EPFADDR = 0x0000204C;
constexpr uint32_t NV_PAPU_TVL2D = 0x00002054;
constexpr uint32_t NV_PAPU_TVL3D = 0x00002060;
constexpr uint32_t NV_PAPU_TVLMP = 0x0000206C;
constexpr uint32_t NV_PAPU_GPSMAXSGE = 0x000020D4;
constexpr uint32_t NV_PAPU_GPFMAXSGE = 0x000020D8;
constexpr uint32_t NV_PAPU_EPSMAXSGE = 0x000020DC;
constexpr uint32_t NV_PAPU_EPFMAXSGE = 0x000020E0;

constexpr uint32_t NV_PAPU_FEAV_VALUE = 0x0000FFFF;
constexpr uint32_t NV_PAPU_FEAV_LST = 0x00030000;

constexpr uint32_t NV1BA0_PIO_SET_ANTECEDENT_VOICE = 0x00000120;
constexpr uint32_t NV1BA0_PIO_VOICE_ON = 0x00000124;
constexpr uint32_t NV1BA0_PIO_VOICE_OFF = 0x00000128;
constexpr uint32_t NV1BA0_PIO_VOICE_PAUSE = 0x00000140;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_VOICE = 0x000002F8;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_VBIN = 0x00000300;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_FMT = 0x00000304;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENV0 = 0x00000308;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENVA = 0x0000030C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENV1 = 0x00000310;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENVF = 0x00000314;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_MISC = 0x00000318;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_VOLA = 0x00000360;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_VOLB = 0x00000364;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_VOLC = 0x00000368;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_LFO_ENV = 0x0000036C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_PITCH = 0x0000037C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_BUF_BASE = 0x000003A0;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_BUF_LBO = 0x000003A4;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_BUF_CBO = 0x000003D8;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_BUF_EBO = 0x000003DC;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE = 0x00000804;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE_OFFSET = 0x00000808;
constexpr uint32_t NV1BA0_PIO_SET_OUTBUF_BA = 0x00001000;
constexpr uint32_t NV1BA0_PIO_SET_OUTBUF_LEN = 0x00001004;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE = 0x00001800;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_OFFSET = 0x00001808;

constexpr uint32_t NV1BA0_PIO_VOICE_ON_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_ON_ENVF = 0x0F000000;
constexpr uint32_t NV1BA0_PIO_VOICE_ON_ENVA = 0xF0000000;
constexpr uint32_t NV1BA0_PIO_VOICE_OFF_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_PAUSE_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_PAUSE_ACTION = 1 << 18;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_PITCH_STEP = 0xFFFF0000;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE_HANDLE = 0xFFFFFFFF;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE_OFFSET_PARAMETER = 0xFFFFF000;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_HANDLE = 0xFFFFFFFF;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_OFFSET_PARAMETER = 0xFFFFF000;

constexpr uint32_t NV_PAVS_SIZE = 0x00000080;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN = 0x00000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT = 0x00000004;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV0 = 0x00000008;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVA = 0x0000000C;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV1 = 0x00000010;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVF = 0x00000014;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC = 0x00000018;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSL_START = 0x00000020;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSH_SAMPLE = 0x00000024;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE = 0x00000054;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET = 0x00000058;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT = 0x0000005C;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLA = 0x00000060;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLB = 0x00000064;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLC = 0x00000068;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_ENV = 0x0000006C;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK = 0x0000007C;

constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_PAUSED = 1 << 18;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_NEW_VOICE = 1 << 20;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE = 1 << 21;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR = 0x0F000000;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EACUR = 0xF0000000;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSL_START_BA = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET_CBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT_EBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH = 0xFFFF0000;

constexpr uint32_t APU_VOICE_LIST_INHERIT = 0;

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

uint32_t ReadRegisterFragment(uint32_t value, uint32_t byteOffset, unsigned size)
{
	const uint32_t shift = byteOffset * 8;
	if (size >= sizeof(uint32_t)) {
		return value;
	}

	const uint32_t mask = (1u << (size * 8)) - 1;
	return (value >> shift) & mask;
}

uint32_t Ctz32(uint32_t value)
{
	uint32_t shift = 0;
	while (((value >> shift) & 1u) == 0u && shift < 32) {
		++shift;
	}
	return shift;
}

bool IsGuestRangeAccessible(uint32_t guestAddress, uint32_t size)
{
	return size > 0 && guestAddress <= PHYSICAL_MAP_SIZE && size <= (PHYSICAL_MAP_SIZE - guestAddress);
}

}

// TODO: Everything :P
// TODO: Audio Processing/Thread

void APUDevice::Init()
{
	PCIBarRegister r;
	r.Raw.type = PCI_BAR_TYPE_MEMORY;
	r.Memory.address = APU_BASE >> 4;
	RegisterBAR(0, APU_SIZE, r.value);

	m_DeviceId = 0x01B0;
	m_VendorId = PCI_VENDOR_ID_NVIDIA;

	Reset();
}

void APUDevice::Reset()
{
	std::memset(m_Registers.data(), 0, m_Registers.size());
	m_VPFifoLevel = 0;
	m_VPFifoLastUpdate = GetAPUTime();
	m_VPInputSgeHandle = 0;
	m_VPOutputSgeHandle = 0;

	SetRegister32(NV_PAPU_ISTS, 0);
	SetRegister32(NV_PAPU_IEN, 0);
	SetRegister32(NV_PAPU_FECTL, 0);
	SetRegister32(NV_PAPU_FECV, 0);
	SetRegister32(NV_PAPU_FEAV, 0);
	SetRegister32(NV_PAPU_FENADDR, 0);
	SetRegister32(NV_PAPU_FEDECMETH, 0);
	SetRegister32(NV_PAPU_FEDECPARAM, 0);
	SetRegister32(NV_PAPU_FEMEMADDR, 0);
	SetRegister32(NV_PAPU_FEMEMDATA, 0);
	SetRegister32(NV_PAPU_FETFORCE0, 0);
	SetRegister32(NV_PAPU_FETFORCE1, 0);
	SetRegister32(NV_PAPU_SECTL, 0);
	SetRegister32(NV_PAPU_VPVADDR, 0);
	SetRegister32(NV_PAPU_VPSGEADDR, 0);
	SetRegister32(NV_PAPU_VPSSLADDR, 0);
	SetRegister32(NV_PAPU_GPSADDR, 0);
	SetRegister32(NV_PAPU_GPFADDR, 0);
	SetRegister32(NV_PAPU_EPSADDR, 0);
	SetRegister32(NV_PAPU_EPFADDR, 0);
	SetRegister32(NV_PAPU_GPSMAXSGE, 0);
	SetRegister32(NV_PAPU_GPFMAXSGE, 0);
	SetRegister32(NV_PAPU_EPSMAXSGE, 0);
	SetRegister32(NV_PAPU_EPFMAXSGE, 0);
	RefreshVPStatus();
}

uint32_t APUDevice::IORead(int barIndex, uint32_t addr, unsigned size)
{
	(void)barIndex;
	(void)addr;
	(void)size;
	return 0;
}

void APUDevice::IOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	(void)barIndex;
	(void)addr;
	(void)value;
	(void)size;
}

uint32_t APUDevice::MMIORead(int barIndex, uint32_t addr, unsigned size)
{
	(void)barIndex;

	if (addr >= APU_VP_BASE && addr < APU_VP_BASE + APU_VP_SIZE) {
		return VPRead(addr - APU_VP_BASE, size);
	}

	if (addr >= APU_GP_BASE && addr < APU_GP_BASE + APU_GP_SIZE) {
		return GPRead(addr - APU_GP_BASE, size);
	}

	if (addr >= APU_EP_BASE && addr < APU_EP_BASE + APU_EP_SIZE) {
		return EPRead(addr - APU_EP_BASE, size);
	}

	if (addr >= NV_PAPU_XGSCNT && addr < NV_PAPU_XGSCNT + sizeof(uint32_t)) {
		return ReadRegisterFragment(GetAPUTime(), addr - NV_PAPU_XGSCNT, size);
	}

	return ReadRegister(addr, size);
}

void APUDevice::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	(void)barIndex;

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

	if (addr >= NV_PAPU_ISTS && addr < NV_PAPU_ISTS + sizeof(uint32_t)) {
		const uint32_t clearMask = value << ((addr - NV_PAPU_ISTS) * 8);
		SetRegister32(NV_PAPU_ISTS, GetRegister32(NV_PAPU_ISTS) & ~clearMask);
		return;
	}

	if (addr >= NV_PAPU_FEMEMDATA && addr < NV_PAPU_FEMEMDATA + sizeof(uint32_t)) {
		WriteRegister(addr, value, size);
		WriteGuestWord(GetRegister32(NV_PAPU_FEMEMADDR), GetRegister32(NV_PAPU_FEMEMDATA));
		return;
	}

	WriteRegister(addr, value, size);
}


uint32_t APUDevice::GPRead(uint32_t addr, unsigned size)
{
	return ReadRegister(APU_GP_BASE + addr, size);
}

void APUDevice::GPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	WriteRegister(APU_GP_BASE + addr, value, size);
}


uint32_t APUDevice::VPRead(uint32_t addr, unsigned size)
{
	UpdateVPFifo();

	if (addr >= APU_VP_FREE && addr < APU_VP_FREE + sizeof(uint32_t)) {
		return ReadRegisterFragment(
			GetRegister32(APU_VP_BASE + APU_VP_FREE),
			addr - APU_VP_FREE,
			size);
	}

	return ReadRegister(APU_VP_BASE + addr, size);
}

void APUDevice::VPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	UpdateVPFifo();

	if (addr >= APU_VP_FREE && addr < APU_VP_FREE + sizeof(uint32_t)) {
		return;
	}

	WriteRegister(APU_VP_BASE + addr, value, size);
	ConsumeVPMethod(addr, value, size);
	if (m_VPFifoLevel < APU_VP_FIFO_CAPACITY) {
		++m_VPFifoLevel;
	}
	RefreshVPStatus();
}


uint32_t APUDevice::EPRead(uint32_t addr, unsigned size)
{
	return ReadRegister(APU_EP_BASE + addr, size);
}

void APUDevice::EPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	WriteRegister(APU_EP_BASE + addr, value, size);
}

uint32_t APUDevice::ReadRegister(uint32_t addr, unsigned size) const
{
	if (size == 0 || addr + size > m_Registers.size()) {
		return 0;
	}

	return ReadLE(m_Registers.data(), addr, size);
}

void APUDevice::WriteRegister(uint32_t addr, uint32_t value, unsigned size)
{
	if (size == 0 || addr + size > m_Registers.size()) {
		return;
	}

	WriteLE(m_Registers.data(), addr, value, size);
}

void APUDevice::SetRegister32(uint32_t addr, uint32_t value)
{
	WriteRegister(addr, value, sizeof(uint32_t));
}

uint32_t APUDevice::GetRegister32(uint32_t addr) const
{
	return ReadRegister(addr, sizeof(uint32_t));
}

void APUDevice::ConsumeVPMethod(uint32_t addr, uint32_t value, unsigned size)
{
	if (size != sizeof(uint32_t)) {
		return;
	}

	SetRegister32(NV_PAPU_FEDECMETH, addr);
	SetRegister32(NV_PAPU_FEDECPARAM, value);

	const auto currentVoice = [this]() {
		return GetRegister32(NV_PAPU_FECV);
	};

	switch (addr) {
	case NV1BA0_PIO_SET_ANTECEDENT_VOICE:
		SetRegister32(NV_PAPU_FEAV, value);
		return;
	case NV1BA0_PIO_SET_CURRENT_VOICE:
		SetRegister32(NV_PAPU_FECV, value);
		return;
	case NV1BA0_PIO_VOICE_ON: {
		const uint32_t selectedHandle = value & NV1BA0_PIO_VOICE_ON_HANDLE;
		if (selectedHandle >= APU_VP_VOICE_MAX_HANDLE) {
			return;
		}

		const uint32_t feav = GetRegister32(NV_PAPU_FEAV);
		const uint32_t list = (feav & NV_PAPU_FEAV_LST) >> Ctz32(NV_PAPU_FEAV_LST);
		if (list != APU_VOICE_LIST_INHERIT) {
			uint32_t topRegister = 0;
			switch (list) {
			case 1: topRegister = NV_PAPU_TVL2D; break;
			case 2: topRegister = NV_PAPU_TVL3D; break;
			case 3: topRegister = NV_PAPU_TVLMP; break;
			default: break;
			}
			if (topRegister != 0) {
				WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
					NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE,
					GetRegister32(topRegister));
				SetRegister32(topRegister, selectedHandle);
			}
		} else {
			const uint32_t antecedentVoice = feav & NV_PAPU_FEAV_VALUE;
			if (antecedentVoice < APU_VP_VOICE_MAX_HANDLE) {
				uint32_t nextHandle = 0;
				if (ReadVoiceMask(antecedentVoice, NV_PAVS_VOICE_TAR_PITCH_LINK,
					NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle)) {
					WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
						NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle);
					WriteVoiceMask(antecedentVoice, NV_PAVS_VOICE_TAR_PITCH_LINK,
						NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, selectedHandle);
				}
			}
		}

		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_EACUR,
			(value & NV1BA0_PIO_VOICE_ON_ENVA) >> Ctz32(NV1BA0_PIO_VOICE_ON_ENVA));
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_EFCUR,
			(value & NV1BA0_PIO_VOICE_ON_ENVF) >> Ctz32(NV1BA0_PIO_VOICE_ON_ENVF));
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_PAUSED, 0);
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 1);
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 1);
		return;
	}
	case NV1BA0_PIO_VOICE_OFF: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_OFF_HANDLE;
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 0);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
		return;
	}
	case NV1BA0_PIO_VOICE_PAUSE: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_PAUSE_HANDLE;
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_PAUSED,
			(value & NV1BA0_PIO_VOICE_PAUSE_ACTION) != 0 ? 1u : 0u);
		return;
	}
	case NV1BA0_PIO_SET_VOICE_CFG_VBIN:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_VBIN, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_FMT:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_FMT, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_ENV0:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_ENV0, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_ENVA:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_ENVA, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_ENV1:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_ENV1, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_ENVF:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_ENVF, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_MISC:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_MISC, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_VOLA:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_VOLA, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_VOLB:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_VOLB, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_VOLC:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_VOLC, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_LFO_ENV:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_LFO_ENV, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_PITCH:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_PITCH_LINK,
			NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH,
			(value & NV1BA0_PIO_SET_VOICE_TAR_PITCH_STEP) >> Ctz32(NV1BA0_PIO_SET_VOICE_TAR_PITCH_STEP));
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_BASE:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CUR_PSL_START,
			NV_PAVS_VOICE_CUR_PSL_START_BA, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_LBO:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CUR_PSH_SAMPLE,
			NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO, value);
		return;
	case NV1BA0_PIO_SET_VOICE_BUF_CBO:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_PAR_OFFSET,
			NV_PAVS_VOICE_PAR_OFFSET_CBO, value);
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_EBO:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_PAR_NEXT,
			NV_PAVS_VOICE_PAR_NEXT_EBO, value);
		return;
	case NV1BA0_PIO_SET_CURRENT_INBUF_SGE:
		m_VPInputSgeHandle = value & NV1BA0_PIO_SET_CURRENT_INBUF_SGE_HANDLE;
		return;
	case NV1BA0_PIO_SET_CURRENT_INBUF_SGE_OFFSET:
		WriteVPScatterGatherEntry(m_VPInputSgeHandle, value & NV1BA0_PIO_SET_CURRENT_INBUF_SGE_OFFSET_PARAMETER);
		return;
	case NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE:
		m_VPOutputSgeHandle = value & NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_HANDLE;
		return;
	case NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_OFFSET:
		WriteVPScatterGatherEntry(m_VPOutputSgeHandle, value & NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_OFFSET_PARAMETER);
		return;
	default:
		if ((addr >= NV1BA0_PIO_SET_OUTBUF_BA && addr < NV1BA0_PIO_SET_OUTBUF_BA + 0x20 && ((addr - NV1BA0_PIO_SET_OUTBUF_BA) % 8) == 0) ||
			(addr >= NV1BA0_PIO_SET_OUTBUF_LEN && addr < NV1BA0_PIO_SET_OUTBUF_LEN + 0x20 && ((addr - NV1BA0_PIO_SET_OUTBUF_LEN) % 8) == 0)) {
			return;
		}
		return;
	}
}

bool APUDevice::ReadGuestWord(uint32_t guestAddress, uint32_t& value) const
{
	if (!IsGuestRangeAccessible(guestAddress, sizeof(uint32_t))) {
		return false;
	}

	std::memcpy(&value, reinterpret_cast<const void*>(static_cast<uintptr_t>(CONTIGUOUS_MEMORY_BASE + guestAddress)), sizeof(value));
	return true;
}

bool APUDevice::WriteGuestWord(uint32_t guestAddress, uint32_t value)
{
	if (!IsGuestRangeAccessible(guestAddress, sizeof(uint32_t))) {
		return false;
	}

	std::memcpy(reinterpret_cast<void*>(static_cast<uintptr_t>(CONTIGUOUS_MEMORY_BASE + guestAddress)), &value, sizeof(value));
	return true;
}

bool APUDevice::WriteGuestWordMasked(uint32_t guestAddress, uint32_t mask, uint32_t value)
{
	if (mask == 0) {
		return true;
	}

	uint32_t current = 0;
	if (!ReadGuestWord(guestAddress, current)) {
		return false;
	}

	const uint32_t shift = mask == 0xFFFFFFFF ? 0 : Ctz32(mask);
	current &= ~mask;
	current |= (value << shift) & mask;
	return WriteGuestWord(guestAddress, current);
}

bool APUDevice::ReadVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t& value) const
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return false;
	}

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		return false;
	}

	const uint32_t voiceBase = voiceTableBase + voiceHandle * NV_PAVS_SIZE + offset;
	uint32_t current = 0;
	if (!ReadGuestWord(voiceBase, current)) {
		return false;
	}

	const uint32_t shift = mask == 0xFFFFFFFF ? 0 : Ctz32(mask);
	value = mask == 0xFFFFFFFF ? current : ((current & mask) >> shift);
	return true;
}

bool APUDevice::WriteVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t value)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return false;
	}

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		return false;
	}

	const uint32_t voiceBase = voiceTableBase + voiceHandle * NV_PAVS_SIZE + offset;
	return WriteGuestWordMasked(voiceBase, mask, value);
}

bool APUDevice::WriteVPScatterGatherEntry(uint32_t handle, uint32_t value)
{
	const uint32_t sgeTableBase = GetRegister32(NV_PAPU_VPSGEADDR);
	if (sgeTableBase == 0) {
		return false;
	}

	const uint32_t sgeBase = sgeTableBase + handle * 8;
	return WriteGuestWord(sgeBase, value);
}

void APUDevice::UpdateVPFifo()
{
	const uint32_t now = GetAPUTime();
	const uint32_t elapsed = now - m_VPFifoLastUpdate;
	if (elapsed > 0) {
		if (elapsed >= m_VPFifoLevel) {
			m_VPFifoLevel = 0;
		} else {
			m_VPFifoLevel -= elapsed;
		}
		m_VPFifoLastUpdate = now;
		RefreshVPStatus();
	}
}

void APUDevice::RefreshVPStatus()
{
	uint32_t status = 0;
	if (m_VPFifoLevel == 0) {
		status |= APU_VP_STATUS_EMPTY;
	}
	SetRegister32(APU_VP_BASE + APU_VP_FREE, status);
}
