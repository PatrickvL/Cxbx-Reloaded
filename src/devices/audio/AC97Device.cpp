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
#include "APUDevice.h"
#include "APUTimer.h"
#include "common/AddressRanges.h"
#include "core/kernel/support/Emu.h"

#include "SDL.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#define LOG_PREFIX CXBXR_MODULE::MCPX

namespace {

constexpr uint32_t AC97_NAM_SIZE = 0x100;
constexpr uint32_t AC97_NABM_SIZE = 0x80;
constexpr uint32_t AC97_MMIO_SIZE = AC97_NAM_SIZE + AC97_NABM_SIZE;

constexpr uint32_t AC97_Reset = 0x00;
constexpr uint32_t AC97_Master_Volume = 0x02;
constexpr uint32_t AC97_Powerdown_Ctrl_Stat = 0x26;
constexpr uint32_t AC97_PCM_Out_Volume = 0x18;
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

constexpr uint32_t BM_BDBAR = 0x00;
constexpr uint32_t BM_CIV = 0x04;
constexpr uint32_t BM_LVI = 0x05;
constexpr uint32_t BM_SR = 0x06;
constexpr uint32_t BM_PICB = 0x08;
constexpr uint32_t BM_PIV = 0x0A;
constexpr uint32_t BM_CR = 0x0B;

constexpr uint32_t AC97_DESCRIPTOR_COUNT = 32;
constexpr uint32_t AC97_DESCRIPTOR_STRIDE = 8;
constexpr uint32_t AC97_DESCRIPTOR_LENGTH_MASK = 0x0000FFFF;
constexpr uint32_t AC97_DESCRIPTOR_IOC = 0x80000000;

constexpr uint32_t GLOB_CNT_WRST = 1 << 0;
constexpr uint32_t GLOB_CNT_CRST = 1 << 1;
constexpr uint32_t GLOB_CNT_MASK = GLOB_CNT_WRST | GLOB_CNT_CRST;
constexpr uint32_t GLOB_STA_PI_INT = 1 << 8;
constexpr uint32_t GLOB_STA_PO_INT = 1 << 9;
constexpr uint32_t GLOB_STA_MC_INT = 1 << 10;
constexpr uint32_t GLOB_STA_CHANNEL_INT_MASK = GLOB_STA_PI_INT | GLOB_STA_PO_INT | GLOB_STA_MC_INT;
constexpr uint32_t GLOB_STA_RDY = 1 << 15;

constexpr uint16_t SR_FIFOE = 1 << 4;
constexpr uint16_t SR_BCIS = 1 << 3;
constexpr uint16_t SR_LVBCI = 1 << 2;
constexpr uint16_t SR_CELV = 1 << 1;
constexpr uint16_t SR_DCH = 1 << 0;
constexpr uint16_t SR_WCLEAR_MASK = SR_FIFOE | SR_BCIS | SR_LVBCI;

constexpr uint8_t CR_FEIE = 1 << 4;
constexpr uint8_t CR_LVBIE = 1 << 3;
constexpr uint8_t CR_IOCE = 1 << 2;
constexpr uint8_t CR_RR = 1 << 1;
constexpr uint8_t CR_RPBM = 1 << 0;
constexpr uint8_t CR_VALID_MASK = 0x1F;

constexpr uint16_t AC97_EXT_AUDIO_ID_VRA = 1 << 0;
constexpr uint16_t AC97_EXT_AUDIO_ID_VRM = 1 << 3;
constexpr uint16_t AC97_EXT_AUDIO_CTRL_MASK = AC97_EXT_AUDIO_ID_VRA | AC97_EXT_AUDIO_ID_VRM;
constexpr uint16_t AC97_POWER_READY = 0x000F;
constexpr uint16_t AC97_MIN_RATE = 8000;
constexpr uint16_t AC97_RATE_48KHZ = 48000;
constexpr uint16_t AC97_VENDOR_SIGMATEL_1 = 0x8384;
constexpr uint16_t AC97_VENDOR_SIGMATEL_2 = 0x7608;
constexpr uint32_t AC97_OUTPUT_CHANNELS = 2;
constexpr uint32_t AC97_OUTPUT_BYTES_PER_FRAME = sizeof(int16_t) * AC97_OUTPUT_CHANNELS;
constexpr uint32_t AC97_MAX_QUEUED_AUDIO_BYTES = APU_TIMER_FREQUENCY * AC97_OUTPUT_BYTES_PER_FRAME / 2;
constexpr uint16_t AC97_VOLUME_MUTE = 0x8000;
constexpr uint16_t AC97_VOLUME_LEFT_MASK = 0x1F00;
constexpr uint16_t AC97_VOLUME_RIGHT_MASK = 0x001F;
constexpr uint32_t AC97_VOLUME_LEFT_SHIFT = 8;
constexpr float AC97_VOLUME_STEP_DB = 1.5f;

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

size_t ChannelIndex(uint32_t channelBase)
{
	switch (channelBase) {
	case NABM_PI_BASE:
		return 0;
	case NABM_PO_BASE:
		return 1;
	case NABM_MC_BASE:
	default:
		return 2;
	}
}

// PIV reports the next descriptor the controller can prefetch. When CIV already
// points at the last valid descriptor, hardware has no later valid entry to
// prefetch, so PIV remains aligned with the current descriptor index.
uint8_t GetPrefetchedIndexValue(uint8_t currentIndex, uint8_t lastValidIndex)
{
	return currentIndex == lastValidIndex
		? currentIndex
		: static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
}

bool HasBusMasterInterrupt(uint16_t status, uint8_t control)
{
	return ((status & SR_FIFOE) != 0 && (control & CR_FEIE) != 0) ||
		((status & SR_BCIS) != 0 && (control & CR_IOCE) != 0) ||
		((status & SR_LVBCI) != 0 && (control & CR_LVBIE) != 0);
}

bool IsGuestRangeAccessible(uint32_t guestAddress, uint32_t size)
{
	return size > 0 && guestAddress <= PHYSICAL_MAP_SIZE && size <= (PHYSICAL_MAP_SIZE - guestAddress);
}

float DecodeOutputAttenuation(uint16_t volumeRegister, bool leftChannel)
{
	if ((volumeRegister & AC97_VOLUME_MUTE) != 0) {
		return 0.0f;
	}

	const uint32_t attenuation = leftChannel
		? ((volumeRegister & AC97_VOLUME_LEFT_MASK) >> AC97_VOLUME_LEFT_SHIFT)
		: (volumeRegister & AC97_VOLUME_RIGHT_MASK);
	return std::pow(10.0f, -(static_cast<float>(attenuation) * AC97_VOLUME_STEP_DB) / 20.0f);
}

int16_t ClampToInt16(int32_t sample)
{
	return static_cast<int16_t>(std::clamp(sample, static_cast<int32_t>(INT16_MIN), static_cast<int32_t>(INT16_MAX)));
}

int32_t ScaleSample(int16_t sample, float gain)
{
	const float scaled = static_cast<float>(sample) * gain;
	return static_cast<int32_t>(scaled >= 0.0f ? (scaled + 0.5f) : (scaled - 0.5f));
}

uint16_t ClampSampleRateRegister(uint16_t value)
{
	if (value == 0) {
		return AC97_RATE_48KHZ;
	}
	return std::clamp(value, AC97_MIN_RATE, AC97_RATE_48KHZ);
}

}

extern APUDevice* g_APU;

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
	m_ChannelLastUpdate.fill(GetAPUTime());
	m_ChannelSampleRemainder.fill(0);
	m_ChannelAdvanceOnRestart.fill(false);
	m_LoggedQueueFull = false;
	if (m_OutputDevice != 0) {
		SDL_ClearQueuedAudio(static_cast<SDL_AudioDeviceID>(m_OutputDevice));
	}

	ResetBusMasterChannel(NABM_PI_BASE);
	ResetBusMasterChannel(NABM_PO_BASE);
	ResetBusMasterChannel(NABM_MC_BASE);
	WriteRegister(AC97_NAM_SIZE + NABM_GLOB_CNT, 0, sizeof(uint32_t));
	UpdateGlobalStatus();
}

bool AC97Device::EnsureOutputDevice()
{
	if (m_OutputDevice != 0) {
		return true;
	}

	if (m_OutputDeviceFailed) {
		return false;
	}

	if ((SDL_WasInit(SDL_INIT_AUDIO) & SDL_INIT_AUDIO) == 0 && SDL_InitSubSystem(SDL_INIT_AUDIO) != 0) {
		m_OutputDeviceFailed = true;
		return false;
	}

	SDL_AudioSpec desired{};
	desired.freq = APU_TIMER_FREQUENCY;
	desired.format = AUDIO_S16SYS;
	desired.channels = static_cast<Uint8>(AC97_OUTPUT_CHANNELS);
	desired.samples = 1024;

	SDL_AudioSpec obtained{};
	const SDL_AudioDeviceID device = SDL_OpenAudioDevice(nullptr, 0, &desired, &obtained, 0);
	if (device == 0) {
		EmuLog(LOG_LEVEL::WARNING, "Failed to open audio device: %s", SDL_GetError());
		m_OutputDeviceFailed = true;
		return false;
	}

	m_OutputDevice = static_cast<uint32_t>(device);
	SDL_PauseAudioDevice(device, 0);
	return true;
}

void AC97Device::SubmitPCMFrames(const int16_t* samples, size_t frameCount)
{
	if (samples == nullptr || frameCount == 0 || !EnsureOutputDevice()) {
		return;
	}

	const SDL_AudioDeviceID device = static_cast<SDL_AudioDeviceID>(m_OutputDevice);
	if (SDL_GetQueuedAudioSize(device) >= AC97_MAX_QUEUED_AUDIO_BYTES) {
		if (!m_LoggedQueueFull) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 output queue full, dropping PCM frames");
			m_LoggedQueueFull = true;
		}
		return;
	}
	m_LoggedQueueFull = false;

	const uint16_t masterVolume = ReadRegister16(AC97_Master_Volume);
	const uint16_t pcmOutVolume = ReadRegister16(AC97_PCM_Out_Volume);
	const float leftGain = DecodeOutputAttenuation(masterVolume, true) * DecodeOutputAttenuation(pcmOutVolume, true);
	const float rightGain = DecodeOutputAttenuation(masterVolume, false) * DecodeOutputAttenuation(pcmOutVolume, false);
	if (leftGain == 1.0f && rightGain == 1.0f) {
		SDL_QueueAudio(device, samples, static_cast<Uint32>(frameCount * AC97_OUTPUT_BYTES_PER_FRAME));
		return;
	}

	m_OutputScratch.resize(frameCount * AC97_OUTPUT_CHANNELS);
	for (size_t frame = 0; frame < frameCount; ++frame) {
		const size_t sampleIndex = frame * AC97_OUTPUT_CHANNELS;
		m_OutputScratch[sampleIndex] = ClampToInt16(ScaleSample(samples[sampleIndex], leftGain));
		m_OutputScratch[sampleIndex + 1] = ClampToInt16(ScaleSample(samples[sampleIndex + 1], rightGain));
	}
	SDL_QueueAudio(device, m_OutputScratch.data(), static_cast<Uint32>(m_OutputScratch.size() * sizeof(int16_t)));
}

uint32_t AC97Device::IORead(int barIndex, uint32_t addr, unsigned size)
{
	if (barIndex == 1) {
		UpdateBusMasterChannels();
	}

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
	if (barIndex == 1) {
		UpdateBusMasterChannels();
	}

	switch (barIndex) {
	case 0:
		if (addr == AC97_Reset && size >= sizeof(uint16_t)) {
			Reset();
			return;
		}
		if (size >= sizeof(uint16_t)) {
			const uint16_t value16 = static_cast<uint16_t>(value);
			switch (addr) {
			case AC97_Extended_Audio_ID:
			case AC97_Vendor_ID1:
			case AC97_Vendor_ID2:
				return;
			case AC97_Extended_Audio_Ctrl_Stat: {
				const uint16_t current = ReadRegister16(AC97_Extended_Audio_Ctrl_Stat);
				const uint16_t next = static_cast<uint16_t>((current & ~AC97_EXT_AUDIO_CTRL_MASK) | (value16 & AC97_EXT_AUDIO_CTRL_MASK));
				WriteRegister16(AC97_Extended_Audio_Ctrl_Stat, next);
				if ((next & AC97_EXT_AUDIO_ID_VRA) == 0) {
					WriteRegister16(AC97_PCM_Front_DAC_Rate, AC97_RATE_48KHZ);
					WriteRegister16(AC97_PCM_Surround_DAC_Rate, AC97_RATE_48KHZ);
					WriteRegister16(AC97_PCM_LFE_DAC_Rate, AC97_RATE_48KHZ);
					WriteRegister16(AC97_PCM_LR_ADC_Rate, AC97_RATE_48KHZ);
				}
				if ((next & AC97_EXT_AUDIO_ID_VRM) == 0) {
					WriteRegister16(AC97_MIC_ADC_Rate, AC97_RATE_48KHZ);
				}
				return;
			}
			case AC97_PCM_Front_DAC_Rate:
			case AC97_PCM_Surround_DAC_Rate:
			case AC97_PCM_LFE_DAC_Rate:
			case AC97_PCM_LR_ADC_Rate:
				if ((ReadRegister16(AC97_Extended_Audio_Ctrl_Stat) & AC97_EXT_AUDIO_ID_VRA) == 0) {
					WriteRegister16(addr, AC97_RATE_48KHZ);
					return;
				}
				WriteRegister16(addr, ClampSampleRateRegister(value16));
				return;
			case AC97_MIC_ADC_Rate:
				if ((ReadRegister16(AC97_Extended_Audio_Ctrl_Stat) & AC97_EXT_AUDIO_ID_VRM) == 0) {
					WriteRegister16(addr, AC97_RATE_48KHZ);
					return;
				}
				WriteRegister16(addr, ClampSampleRateRegister(value16));
				return;
			default:
				break;
			}
		}
		WriteRegister(addr, value, size);
		return;
	case 1:
		addr += AC97_NAM_SIZE;
		if (addr == AC97_NAM_SIZE + NABM_GLOB_CNT && size >= sizeof(uint32_t)) {
			const uint32_t control = value & GLOB_CNT_MASK;
			if ((control & (GLOB_CNT_WRST | GLOB_CNT_CRST)) != 0) {
				Reset();
				return;
			}
			WriteRegister(addr, control, sizeof(uint32_t));
			UpdateGlobalStatus();
			return;
		}
		if (addr == AC97_NAM_SIZE + NABM_GLOB_STA) {
			UpdateGlobalStatus();
			return;
		}
		{
			const uint32_t channelAddr = addr - AC97_NAM_SIZE;
			const uint32_t channelBase = channelAddr & ~0xFu;
			if (channelBase == NABM_PI_BASE || channelBase == NABM_PO_BASE || channelBase == NABM_MC_BASE) {
				switch (channelAddr & 0x0Fu) {
				case BM_BDBAR:
					if (size >= sizeof(uint32_t)) {
						WriteRegister(addr, value & ~0x7u, sizeof(uint32_t));
						m_ChannelAdvanceOnRestart[ChannelIndex(channelBase)] = false;
						UpdateBusMasterStatus(channelBase);
					}
					return;
				case BM_CIV:
				case BM_PICB:
				case BM_PIV:
					return;
				case BM_LVI:
				{
					const uint32_t baseAddr = AC97_NAM_SIZE + channelBase;
					const uint8_t previousLastValid = static_cast<uint8_t>(ReadRegister(baseAddr + BM_LVI, sizeof(uint8_t)) & 0x1F);
					const uint8_t newLastValid = static_cast<uint8_t>(value & 0x1F);
					WriteRegister(addr, newLastValid, sizeof(uint8_t));

					const uint8_t control = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CR, sizeof(uint8_t)));
					const uint16_t status = ReadRegister16(baseAddr + BM_SR);
					const uint16_t remaining = ReadRegister16(baseAddr + BM_PICB);
					const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CIV, sizeof(uint8_t)) & 0x1F);
					if ((control & CR_RPBM) != 0 &&
						(status & SR_DCH) != 0 &&
						remaining == 0 &&
						currentIndex == previousLastValid &&
						newLastValid != currentIndex) {
						const uint8_t nextIndex = static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
						WriteRegister(baseAddr + BM_CIV, nextIndex, sizeof(uint8_t));
						m_ChannelAdvanceOnRestart[ChannelIndex(channelBase)] = false;
					}

					UpdateBusMasterStatus(channelBase);
					return;
				}
				case BM_SR: {
					const uint16_t current = ReadRegister16(addr);
					WriteRegister16(addr, current & ~(static_cast<uint16_t>(value) & SR_WCLEAR_MASK));
					UpdateGlobalStatus();
					return;
				}
				case BM_CR: {
					const size_t channelIndex = ChannelIndex(channelBase);
					const uint32_t baseAddr = AC97_NAM_SIZE + channelBase;
					const uint8_t previousControl = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CR, sizeof(uint8_t)));
					const uint8_t control = static_cast<uint8_t>(value) & CR_VALID_MASK;
					WriteRegister(addr, control, sizeof(uint8_t));
					if ((control & CR_RR) != 0) {
						ResetBusMasterChannel(channelBase);
					} else {
						if ((previousControl & CR_RPBM) == 0 &&
							(control & CR_RPBM) != 0 &&
							m_ChannelAdvanceOnRestart[channelIndex]) {
							const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CIV, sizeof(uint8_t)) & 0x1F);
							const uint8_t lastValidIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_LVI, sizeof(uint8_t)) & 0x1F);
							// If software stopped the engine after it halted on the previous LVI
							// and then extended LVI, restart from the newly queued descriptor.
							if (ReadRegister16(baseAddr + BM_PICB) == 0 && currentIndex != lastValidIndex) {
								const uint8_t nextIndex = static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
								WriteRegister(baseAddr + BM_CIV, nextIndex, sizeof(uint8_t));
								m_ChannelAdvanceOnRestart[channelIndex] = false;
							}
						}
						UpdateBusMasterStatus(channelBase);
					}
					return;
				}
				default:
					break;
				}
			}
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
		if (addr >= AC97_NAM_SIZE) {
			UpdateBusMasterChannels();
		}
		return ReadRegister(addr, size);
	}

	return 0;
}

void AC97Device::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	(void)barIndex;

	if (addr < AC97_MMIO_SIZE) {
		if (addr >= AC97_NAM_SIZE) {
			UpdateBusMasterChannels();
		}
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

void AC97Device::UpdateGlobalStatus()
{
	uint32_t status = ReadRegister(AC97_NAM_SIZE + NABM_GLOB_STA, sizeof(uint32_t)) & ~GLOB_STA_CHANNEL_INT_MASK;

	if ((ReadRegister16(AC97_Powerdown_Ctrl_Stat) & AC97_POWER_READY) == AC97_POWER_READY) {
		status |= GLOB_STA_RDY;
	} else {
		status &= ~GLOB_STA_RDY;
	}

	const uint16_t piStatus = ReadRegister16(AC97_NAM_SIZE + NABM_PI_BASE + BM_SR);
	const uint16_t poStatus = ReadRegister16(AC97_NAM_SIZE + NABM_PO_BASE + BM_SR);
	const uint16_t mcStatus = ReadRegister16(AC97_NAM_SIZE + NABM_MC_BASE + BM_SR);
	const uint8_t piControl = static_cast<uint8_t>(ReadRegister(AC97_NAM_SIZE + NABM_PI_BASE + BM_CR, sizeof(uint8_t)));
	const uint8_t poControl = static_cast<uint8_t>(ReadRegister(AC97_NAM_SIZE + NABM_PO_BASE + BM_CR, sizeof(uint8_t)));
	const uint8_t mcControl = static_cast<uint8_t>(ReadRegister(AC97_NAM_SIZE + NABM_MC_BASE + BM_CR, sizeof(uint8_t)));
	if (HasBusMasterInterrupt(piStatus, piControl)) {
		status |= GLOB_STA_PI_INT;
	}
	if (HasBusMasterInterrupt(poStatus, poControl)) {
		status |= GLOB_STA_PO_INT;
	}
	if (HasBusMasterInterrupt(mcStatus, mcControl)) {
		status |= GLOB_STA_MC_INT;
	}

	WriteRegister(AC97_NAM_SIZE + NABM_GLOB_STA, status, sizeof(uint32_t));
}

void AC97Device::UpdateBusMasterChannels()
{
	if (g_APU != nullptr) {
		g_APU->SynchronizeAudio();
	}
	UpdateBusMasterStatus(NABM_PI_BASE);
	UpdateBusMasterStatus(NABM_PO_BASE);
	UpdateBusMasterStatus(NABM_MC_BASE);
	UpdateGlobalStatus();
}

void AC97Device::ResetBusMasterChannel(uint32_t channelBase)
{
	const size_t channelIndex = ChannelIndex(channelBase);
	m_ChannelLastUpdate[channelIndex] = GetAPUTime();
	m_ChannelSampleRemainder[channelIndex] = 0;
	m_ChannelAdvanceOnRestart[channelIndex] = false;

	WriteRegister(AC97_NAM_SIZE + channelBase + BM_BDBAR, 0, sizeof(uint32_t));
	WriteRegister(AC97_NAM_SIZE + channelBase + BM_CIV, 0, sizeof(uint8_t));
	WriteRegister(AC97_NAM_SIZE + channelBase + BM_LVI, 0, sizeof(uint8_t));
	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_SR, SR_DCH | SR_CELV);
	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_PICB, 0);
	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_PIV, 0);
	WriteRegister(AC97_NAM_SIZE + channelBase + BM_CR, 0, sizeof(uint8_t));
	UpdateGlobalStatus();
}

AC97Device::PrimeResult AC97Device::PrimeBusMasterChannel(uint32_t channelBase)
{
	const uint32_t baseAddr = AC97_NAM_SIZE + channelBase;
	const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CIV, sizeof(uint8_t)) & 0x1F);
	const uint8_t lastValidIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_LVI, sizeof(uint8_t)) & 0x1F);
	const uint32_t descriptorBase = ReadRegister(baseAddr + BM_BDBAR, sizeof(uint32_t)) & ~0x7u;
	if (descriptorBase == 0) {
		return PrimeResult::DescriptorError;
	}

	uint32_t descriptorControl = 0;
	if (!ReadGuest32(descriptorBase + currentIndex * AC97_DESCRIPTOR_STRIDE + 4, descriptorControl)) {
		return PrimeResult::DescriptorError;
	}

	const uint16_t descriptorLength = static_cast<uint16_t>(descriptorControl & AC97_DESCRIPTOR_LENGTH_MASK);
	if (descriptorLength == 0) {
		return currentIndex == lastValidIndex ? PrimeResult::EndOfList : PrimeResult::DescriptorError;
	}

	WriteRegister16(baseAddr + BM_PICB, descriptorLength);
	WriteRegister16(baseAddr + BM_PIV, GetPrefetchedIndexValue(currentIndex, lastValidIndex));
	return PrimeResult::Ready;
}

void AC97Device::UpdateBusMasterStatus(uint32_t channelBase)
{
	const uint32_t crAddr = AC97_NAM_SIZE + channelBase + BM_CR;
	const uint32_t srAddr = AC97_NAM_SIZE + channelBase + BM_SR;
	const uint32_t civAddr = AC97_NAM_SIZE + channelBase + BM_CIV;
	const uint32_t lviAddr = AC97_NAM_SIZE + channelBase + BM_LVI;
	const uint32_t picbAddr = AC97_NAM_SIZE + channelBase + BM_PICB;
	const uint32_t bdbarAddr = AC97_NAM_SIZE + channelBase + BM_BDBAR;
	const uint8_t control = static_cast<uint8_t>(ReadRegister(crAddr, sizeof(uint8_t)));
	uint16_t status = ReadRegister16(srAddr) & SR_WCLEAR_MASK;
	const size_t channelIndex = ChannelIndex(channelBase);
	const uint32_t now = GetAPUTime();
	const uint32_t elapsed = now - m_ChannelLastUpdate[channelIndex];
	m_ChannelLastUpdate[channelIndex] = now;

	const auto primeChannel = [&]() {
		switch (PrimeBusMasterChannel(channelBase)) {
		case PrimeResult::Ready:
			return true;
		case PrimeResult::EndOfList:
			m_ChannelAdvanceOnRestart[channelIndex] = true;
			status |= SR_LVBCI;
			return false;
		case PrimeResult::DescriptorError:
			m_ChannelAdvanceOnRestart[channelIndex] = false;
			status |= SR_FIFOE;
			return false;
		}
		return false;
	};

	if (control & CR_RPBM) {
		uint64_t samplesToConsume = m_ChannelSampleRemainder[channelIndex];
		samplesToConsume += static_cast<uint64_t>(elapsed) * GetBusMasterSampleRate(channelBase);
		m_ChannelSampleRemainder[channelIndex] = static_cast<uint32_t>(samplesToConsume % APU_TIMER_FREQUENCY);
		samplesToConsume /= APU_TIMER_FREQUENCY;

		if (ReadRegister16(picbAddr) == 0 && !primeChannel()) {
			samplesToConsume = 0;
		}

		while (samplesToConsume > 0) {
			uint16_t remaining = ReadRegister16(picbAddr);
			if (remaining == 0 && !primeChannel()) {
				break;
			}

			remaining = ReadRegister16(picbAddr);
			if (remaining == 0) {
				break;
			}

			const uint16_t consumed = static_cast<uint16_t>(samplesToConsume > remaining ? remaining : samplesToConsume);
			remaining = static_cast<uint16_t>(remaining - consumed);
			WriteRegister16(picbAddr, remaining);
			samplesToConsume -= consumed;

			if (remaining != 0) {
				continue;
			}

			const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(civAddr, sizeof(uint8_t)) & 0x1F);
			const uint8_t lastValidIndex = static_cast<uint8_t>(ReadRegister(lviAddr, sizeof(uint8_t)) & 0x1F);
			const uint32_t descriptorBase = ReadRegister(bdbarAddr, sizeof(uint32_t)) & ~0x7u;
			uint32_t descriptorControl = 0;
			if (descriptorBase != 0 &&
				ReadGuest32(descriptorBase + currentIndex * AC97_DESCRIPTOR_STRIDE + 4, descriptorControl) &&
				(descriptorControl & AC97_DESCRIPTOR_IOC) != 0) {
				status |= SR_BCIS;
			}

			if (currentIndex == lastValidIndex) {
				m_ChannelAdvanceOnRestart[channelIndex] = true;
				status |= SR_LVBCI;
				break;
			}

			const uint8_t nextIndex = static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
			WriteRegister(civAddr, nextIndex, sizeof(uint8_t));
			m_ChannelAdvanceOnRestart[channelIndex] = false;
			if (!primeChannel()) {
				break;
			}
		}
	}

	if ((control & CR_RPBM) == 0 || ReadRegister16(picbAddr) == 0) {
		status |= SR_DCH;
	} else {
		status &= ~SR_DCH;
	}

	const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(civAddr, sizeof(uint8_t)) & 0x1F);
	const uint8_t lastValidIndex = static_cast<uint8_t>(ReadRegister(lviAddr, sizeof(uint8_t)) & 0x1F);
	// CELV reflects descriptor position, not whether DMA is currently running.
	if (currentIndex == lastValidIndex) {
		status |= SR_CELV;
	} else {
		status &= ~SR_CELV;
	}

	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_PIV, GetPrefetchedIndexValue(currentIndex, lastValidIndex));
	WriteRegister16(srAddr, status);
	UpdateGlobalStatus();
}

uint32_t AC97Device::GetBusMasterSampleRate(uint32_t channelBase) const
{
	switch (channelBase) {
	case NABM_PI_BASE:
		return ReadRegister16(AC97_PCM_LR_ADC_Rate);
	case NABM_MC_BASE:
		return ReadRegister16(AC97_MIC_ADC_Rate);
	case NABM_PO_BASE:
	default:
		return ReadRegister16(AC97_PCM_Front_DAC_Rate);
	}
}

bool AC97Device::ReadGuest32(uint32_t guestAddress, uint32_t& value) const
{
	if (!IsGuestRangeAccessible(guestAddress, sizeof(uint32_t))) {
		return false;
	}

	std::memcpy(&value, reinterpret_cast<const void*>(static_cast<uintptr_t>(CONTIGUOUS_MEMORY_BASE + guestAddress)), sizeof(value));
	return true;
}
