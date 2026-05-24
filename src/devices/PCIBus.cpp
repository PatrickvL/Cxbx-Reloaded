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
// *  (c) 2018 Patrick van Logchem <pvanlogchem@gmail.com>
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::X86

#include "PCIBus.h"
#include <cstdio>
#include "common\Logging.h"
#include "audio\AudioDiagnostics.h"

void PCIBus::ConnectDevice(uint32_t deviceId, PCIDevice *pDevice)
{
	if (m_Devices.find(deviceId) != m_Devices.end()) {
		printf("PCIBus: Attempting to connect two devices to the same device address\n");
		return;
	}

	m_Devices[deviceId] = pDevice;
	pDevice->Init();
}

void PCIBus::IOWriteConfigAddress(uint32_t pData) 
{
	memcpy(&m_configAddressRegister, &pData, sizeof(PCIConfigAddressRegister));
}

uint32_t PCIBus::IOReadConfigData()
{
	auto it = m_Devices.find(PCI_DEVID(m_configAddressRegister.busNumber, PCI_DEVFN(m_configAddressRegister.deviceNumber, m_configAddressRegister.functionNumber)));
	if (it != m_Devices.end()) {
		return it->second->ReadConfigRegister(m_configAddressRegister.registerNumber & PCI_CONFIG_REGISTER_MASK);
	}

	printf("PCIBus::IOReadConfigData: Invalid Device Read (Bus: %d\t Slot: %d\t Function: %d)\n", m_configAddressRegister.busNumber,
		m_configAddressRegister.deviceNumber,
		m_configAddressRegister.functionNumber);

	// Unpopulated PCI slots return 0xFFFFFFFF
	return 0xFFFFFFFF;
}

void PCIBus::IOWriteConfigData(uint32_t pData) {
	auto it = m_Devices.find(PCI_DEVID(m_configAddressRegister.busNumber, PCI_DEVFN(m_configAddressRegister.deviceNumber, m_configAddressRegister.functionNumber)));
	if (it != m_Devices.end()) {
		it->second->WriteConfigRegister(m_configAddressRegister.registerNumber & PCI_CONFIG_REGISTER_MASK, pData);
		return;
	}

	printf("PCIBus::IOWriteConfigData: Invalid Device Write (Bus: %d\t Slot: %d\t Function: %d)\n", m_configAddressRegister.busNumber,
		m_configAddressRegister.deviceNumber,
		m_configAddressRegister.functionNumber);
}

bool PCIBus::IORead(uint32_t addr, uint32_t* data, unsigned size)
{
	switch (addr) {
	case PORT_PCI_CONFIG_DATA: // 0xCFC
	case PORT_PCI_CONFIG_DATA + 1:
	case PORT_PCI_CONFIG_DATA + 2:
	case PORT_PCI_CONFIG_DATA + 3:
	{
		uint32_t configData = IOReadConfigData();
		int byteOffset = addr - PORT_PCI_CONFIG_DATA;
		switch (size) {
		case sizeof(uint32_t):
			*data = configData;
			break;
		case sizeof(uint16_t):
			*data = (configData >> (byteOffset * 8)) & 0xFFFF;
			break;
		case sizeof(uint8_t):
			*data = (configData >> (byteOffset * 8)) & 0xFF;
			break;
		}
		return true;
	}
	default:
		for (auto it = m_Devices.begin(); it != m_Devices.end(); ++it) {
			PCIBar bar;
			if (it->second->GetIOBar(addr, &bar)) {
				*data = it->second->IORead(bar.index, addr - bar.reg.IO.address, size);
				return true;
			}
		}
	}

	return false;
}

bool PCIBus::IOWrite(uint32_t addr, uint32_t value, unsigned size)
{
	switch (addr) {
	case PORT_PCI_CONFIG_ADDRESS: // 0xCF8
		if (size == sizeof(uint32_t)) {
			IOWriteConfigAddress(value);
			return true;
		} // TODO : else log wrong size-access?
		break;
	case PORT_PCI_CONFIG_DATA: // 0xCFC
		if (size == sizeof(uint32_t)) {
			IOWriteConfigData(value);
			return true; // TODO : Should IOWriteConfigData() success/failure be returned?
		} // TODO : else log wrong size-access?
		break;
	default:
		for (auto it = m_Devices.begin(); it != m_Devices.end(); ++it) {
			PCIBar bar;
			if (it->second->GetIOBar(addr, &bar)) {
				it->second->IOWrite(bar.index, addr - bar.reg.IO.address, value, size);
				return true;
			}
		}
	}

	return false;
}

bool PCIBus::MMIORead(uint32_t addr, uint32_t* data, unsigned size)
{
	for (auto it = m_Devices.begin(); it != m_Devices.end(); ++it) {
		PCIBar bar;
		if (it->second->GetMMIOBar(addr, &bar)) {
			*data = it->second->MMIORead(bar.index, addr - (bar.reg.Memory.address << 4), size);
			return true;
		}
	}

	return false;
}

bool PCIBus::MMIOWrite(uint32_t addr, uint32_t value, unsigned size)
{
	bool traceAPUWrite = false;
	const char* apuTraceType = nullptr;
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		if (IsAPUVPBaseRegisterTrace(addr)) {
			traceAPUWrite = true;
			apuTraceType = "VP base-register";
		} else if (IsAPUVPMethodTrace(addr)) {
			traceAPUWrite = true;
			apuTraceType = "VP method-window";
		}
	}

	for (auto it = m_Devices.begin(); it != m_Devices.end(); ++it) {
		PCIBar bar;
		if (it->second->GetMMIOBar(addr, &bar)) {
			if (traceAPUWrite) {
				EmuLog(LOG_LEVEL::INFO,
					"PCI MMIO %s write routed addr=0x%08x device=0x%08x bar=%d barBase=0x%08x offset=0x%08x value=0x%08x size=%u",
					apuTraceType,
					addr,
					it->second->ReadConfigRegister(PCI_CONFIG_DEVICE),
					bar.index,
					(bar.reg.Memory.address << 4),
					addr - (bar.reg.Memory.address << 4),
					value,
					size);
			}
			it->second->MMIOWrite(bar.index, addr - (bar.reg.Memory.address << 4), value, size);
			return true;
		}
	}

	if (traceAPUWrite) {
		EmuLog(LOG_LEVEL::INFO,
			"PCI MMIO %s write was not routed addr=0x%08x value=0x%08x size=%u",
			apuTraceType,
			addr,
			value,
			size);
	}

	return false;
}

void PCIBus::Reset()
{
	for (auto it = m_Devices.begin(); it != m_Devices.end(); ++it) {
		it->second->Reset();
	}
}
