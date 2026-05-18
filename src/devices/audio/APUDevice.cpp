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

#include <cstring>

namespace {

constexpr uint32_t APU_VP_BASE = 0x20000;
constexpr uint32_t APU_VP_SIZE = 0x10000;
constexpr uint32_t APU_VP_FREE = 0x10;
constexpr uint32_t APU_VP_FIFO_CAPACITY = 0x80;
constexpr uint32_t APU_VP_STATUS_EMPTY = 0x80;

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
constexpr uint32_t NV_PAPU_GPSMAXSGE = 0x000020D4;
constexpr uint32_t NV_PAPU_GPFMAXSGE = 0x000020D8;
constexpr uint32_t NV_PAPU_EPSMAXSGE = 0x000020DC;
constexpr uint32_t NV_PAPU_EPFMAXSGE = 0x000020E0;

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
