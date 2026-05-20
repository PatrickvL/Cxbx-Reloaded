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

#include <AL/al.h>
#include <AL/alc.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>

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
constexpr uint32_t AC97_STREAM_BUFFER_BYTES = 1024;
constexpr uint32_t AC97_STREAM_BUFFER_FRAMES = AC97_STREAM_BUFFER_BYTES / AC97_OUTPUT_BYTES_PER_FRAME;
constexpr uint32_t AC97_STARTUP_BUFFER_CHUNKS = 8;
constexpr uint32_t AC97_STARTUP_BUFFER_FRAMES = AC97_STREAM_BUFFER_FRAMES * AC97_STARTUP_BUFFER_CHUNKS;
constexpr uint32_t AC97_MAX_BUFFERED_AUDIO_BYTES = AC97_MAX_QUEUED_AUDIO_BYTES * 2;
constexpr uint16_t AC97_VOLUME_MUTE = 0x8000;
constexpr uint16_t AC97_VOLUME_LEFT_MASK = 0x1F00;
constexpr uint16_t AC97_VOLUME_RIGHT_MASK = 0x001F;
constexpr uint32_t AC97_VOLUME_LEFT_SHIFT = 8;
constexpr float AC97_VOLUME_STEP_DB = 1.5f;
constexpr const char* AC97_OPENAL_DEVICE_ENV = "CXBXR_OPENAL_DEVICE";
constexpr size_t AC97_OPENAL_DEVICE_NAME_LIMIT = 4096;
constexpr size_t AC97_OPENAL_DEVICE_COUNT_LIMIT = 256;

const char* GetRequestedOpenALDevice()
{
	static std::once_flag requestedDeviceOnce;
	static std::string requestedDevice;
	std::call_once(requestedDeviceOnce, []() {
		const char* environmentValue = std::getenv(AC97_OPENAL_DEVICE_ENV);
		if (environmentValue != nullptr && environmentValue[0] != '\0') {
			requestedDevice = environmentValue;
		}
	});
	return requestedDevice.empty() ? nullptr : requestedDevice.c_str();
}

std::vector<std::string> GetAvailableOpenALDevices()
{
	const ALCchar* deviceList = nullptr;
	if (alcIsExtensionPresent(nullptr, "ALC_ENUMERATE_ALL_EXT") == ALC_TRUE) {
		deviceList = alcGetString(nullptr, ALC_ALL_DEVICES_SPECIFIER);
	} else if (alcIsExtensionPresent(nullptr, "ALC_ENUMERATION_EXT") == ALC_TRUE) {
		deviceList = alcGetString(nullptr, ALC_DEVICE_SPECIFIER);
	}

	std::vector<std::string> devices;
	if (deviceList == nullptr) {
		return devices;
	}

	for (const ALCchar* current = deviceList; *current != '\0' && devices.size() < AC97_OPENAL_DEVICE_COUNT_LIMIT;) {
		const size_t nameLength = strnlen(current, AC97_OPENAL_DEVICE_NAME_LIMIT);
		if (nameLength == AC97_OPENAL_DEVICE_NAME_LIMIT) {
			EmuLog(LOG_LEVEL::WARNING, "OpenAL device enumeration stopped after an unexpectedly long device name");
			break;
		}
		devices.emplace_back(current, nameLength);
		current += nameLength + 1;
	}
	return devices;
}

std::string JoinOpenALDeviceNames(const std::vector<std::string>& devices)
{
	if (devices.empty()) {
		return "<not reported>";
	}

	std::string joinedNames;
	for (size_t i = 0; i < devices.size(); ++i) {
		if (i != 0) {
			joinedNames += ", ";
		}
		joinedNames += devices[i];
	}
	return joinedNames;
}

const char* GetOpenedOpenALDeviceName(ALCdevice* device)
{
	if (device == nullptr) {
		return nullptr;
	}

	const ALCchar* deviceName = nullptr;
	if (alcIsExtensionPresent(device, "ALC_ENUMERATE_ALL_EXT") == ALC_TRUE) {
		deviceName = alcGetString(device, ALC_ALL_DEVICES_SPECIFIER);
	}
	if (deviceName == nullptr || deviceName[0] == '\0') {
		deviceName = alcGetString(device, ALC_DEVICE_SPECIFIER);
	}
	return (deviceName != nullptr && deviceName[0] != '\0') ? deviceName : nullptr;
}

const char* GetOpenALSourceStateName(ALint state)
{
	switch (state) {
	case AL_INITIAL:
		return "AL_INITIAL";
	case AL_PLAYING:
		return "AL_PLAYING";
	case AL_PAUSED:
		return "AL_PAUSED";
	case AL_STOPPED:
		return "AL_STOPPED";
	default:
		return "AL_UNKNOWN_STATE";
	}
}

const char* GetOpenALErrorName(ALenum error)
{
	switch (error) {
	case AL_NO_ERROR:
		return "AL_NO_ERROR";
	case AL_INVALID_NAME:
		return "AL_INVALID_NAME";
	case AL_INVALID_ENUM:
		return "AL_INVALID_ENUM";
	case AL_INVALID_VALUE:
		return "AL_INVALID_VALUE";
	case AL_INVALID_OPERATION:
		return "AL_INVALID_OPERATION";
	case AL_OUT_OF_MEMORY:
		return "AL_OUT_OF_MEMORY";
	default:
		return "AL_UNKNOWN_ERROR";
	}
}

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

uint8_t GetNextDescriptorIndex(uint8_t currentIndex)
{
	return static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
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
	m_ChannelQueuedAfterHalt.fill(false);
	m_ChannelDescriptorError.fill(false);
	m_Pending3DVoices.clear();
	m_LoggedQueueFull = false;
	m_LoggedPlaybackStartFailure = false;
	ResetOutputStream();

	ResetBusMasterChannel(NABM_PI_BASE);
	ResetBusMasterChannel(NABM_PO_BASE);
	ResetBusMasterChannel(NABM_MC_BASE);
	WriteRegister(AC97_NAM_SIZE + NABM_GLOB_CNT, 0, sizeof(uint32_t));
	UpdateGlobalStatus();
}

bool AC97Device::EnsureOutputDevice()
{
	if (m_OutputContext != nullptr) {
		if (alcGetCurrentContext() != m_OutputContext && !alcMakeContextCurrent(m_OutputContext)) {
			return false;
		}
		return true;
	}

	if (m_OutputDeviceFailed) {
		return false;
	}

	const char* requestedDevice = GetRequestedOpenALDevice();
	const std::vector<std::string> availableDevices = GetAvailableOpenALDevices();
	const std::string availableDeviceList = JoinOpenALDeviceNames(availableDevices);
	bool openedRequestedDevice = false;

	if (requestedDevice != nullptr) {
		m_OutputDevice = alcOpenDevice(requestedDevice);
		openedRequestedDevice = (m_OutputDevice != nullptr);
		if (m_OutputDevice == nullptr) {
			EmuLog(LOG_LEVEL::WARNING,
				"Failed to open requested OpenAL device '%s' from %s; falling back to the default device. Available devices: %s",
				requestedDevice, AC97_OPENAL_DEVICE_ENV, availableDeviceList.c_str());
		}
	}

	if (m_OutputDevice == nullptr) {
		m_OutputDevice = alcOpenDevice(nullptr);
	}
	if (m_OutputDevice == nullptr) {
		if (requestedDevice != nullptr) {
			EmuLog(LOG_LEVEL::WARNING,
				"Failed to open both the requested OpenAL device '%s' and the default OpenAL device. Available devices: %s",
				requestedDevice, availableDeviceList.c_str());
		} else {
			EmuLog(LOG_LEVEL::WARNING, "Failed to open OpenAL device. Available devices: %s", availableDeviceList.c_str());
		}
		m_OutputDeviceFailed = true;
		return false;
	}

	const char* openedDeviceName = GetOpenedOpenALDeviceName(m_OutputDevice);
	if (openedDeviceName != nullptr) {
		if (openedRequestedDevice) {
			EmuLog(LOG_LEVEL::INFO, "Opened OpenAL output device '%s' (requested via %s)", openedDeviceName, AC97_OPENAL_DEVICE_ENV);
		} else {
			const char* logMessage = (requestedDevice != nullptr)
				? "Opened default OpenAL output device '%s' after requested-device fallback"
				: "Opened default OpenAL output device '%s'";
			EmuLog(LOG_LEVEL::INFO, logMessage, openedDeviceName);
		}
	}

	m_OutputContext = alcCreateContext(m_OutputDevice, nullptr);
	if (m_OutputContext == nullptr || !alcMakeContextCurrent(m_OutputContext)) {
		EmuLog(LOG_LEVEL::WARNING, "Failed to create OpenAL context");
		if (m_OutputContext != nullptr) {
			alcDestroyContext(m_OutputContext);
			m_OutputContext = nullptr;
		}
		alcCloseDevice(m_OutputDevice);
		m_OutputDevice = nullptr;
		m_OutputDeviceFailed = true;
		return false;
	}

	alGenSources(1, &m_OutputSource);
	if (alGetError() != AL_NO_ERROR || m_OutputSource == 0) {
		EmuLog(LOG_LEVEL::WARNING, "Failed to create OpenAL source");
		alcMakeContextCurrent(nullptr);
		alcDestroyContext(m_OutputContext);
		alcCloseDevice(m_OutputDevice);
		m_OutputContext = nullptr;
		m_OutputDevice = nullptr;
		m_OutputDeviceFailed = true;
		return false;
	}

	alGenBuffers(static_cast<ALsizei>(m_OutputBuffers.size()), m_OutputBuffers.data());
	if (alGetError() != AL_NO_ERROR) {
		EmuLog(LOG_LEVEL::WARNING, "Failed to create OpenAL stream buffers");
		alDeleteSources(1, &m_OutputSource);
		m_OutputSource = 0;
		alcMakeContextCurrent(nullptr);
		alcDestroyContext(m_OutputContext);
		alcCloseDevice(m_OutputDevice);
		m_OutputContext = nullptr;
		m_OutputDevice = nullptr;
		m_OutputDeviceFailed = true;
		return false;
	}

	m_OutputBufferIndex.clear();
	for (size_t i = 0; i < m_OutputBuffers.size(); ++i) {
		m_OutputBufferIndex.emplace(m_OutputBuffers[i], i);
	}
	m_FreeOutputBuffers.assign(m_OutputBuffers.begin(), m_OutputBuffers.end());
	m_OutputBufferBytes.fill(0);
	m_QueuedAudioBytes = 0;
	m_StagedOutputFrames.clear();
	m_LastOutputSourceState = -1;
	m_LoggedPlaybackStartFailure = false;
	alSourcef(m_OutputSource, AL_GAIN, 1.0f);
	return true;
}

bool AC97Device::QueryOutputSourceSnapshot(ALint& state, ALint& queued, ALint& processed, ALenum& error) const
{
	error = AL_NO_ERROR;

	alGetSourcei(m_OutputSource, AL_SOURCE_STATE, &state);
	error = alGetError();
	if (error != AL_NO_ERROR) {
		return false;
	}

	alGetSourcei(m_OutputSource, AL_BUFFERS_QUEUED, &queued);
	error = alGetError();
	if (error != AL_NO_ERROR) {
		return false;
	}

	alGetSourcei(m_OutputSource, AL_BUFFERS_PROCESSED, &processed);
	error = alGetError();
	return error == AL_NO_ERROR;
}

void AC97Device::LogOutputOperationFailure(const char* operation, ALenum sourceError) const
{
	ALint state = AL_INITIAL;
	ALint queued = 0;
	ALint processed = 0;
	ALenum snapshotError = AL_NO_ERROR;
	if (QueryOutputSourceSnapshot(state, queued, processed, snapshotError)) {
		EmuLog(LOG_LEVEL::WARNING,
			"AC97 OpenAL %s failed: %s (0x%04x), state=%s queuedBuffers=%d processedBuffers=%d queuedBytes=%u stagedFrames=%zu",
			operation, GetOpenALErrorName(sourceError), static_cast<unsigned>(sourceError),
			GetOpenALSourceStateName(state), queued, processed,
			m_QueuedAudioBytes, m_StagedOutputFrames.size() / AC97_OUTPUT_CHANNELS);
	} else {
		EmuLog(LOG_LEVEL::WARNING,
			"AC97 OpenAL %s failed: %s (0x%04x), snapshot query also failed with %s (0x%04x), queuedBytes=%u stagedFrames=%zu",
			operation, GetOpenALErrorName(sourceError), static_cast<unsigned>(sourceError),
			GetOpenALErrorName(snapshotError), static_cast<unsigned>(snapshotError),
			m_QueuedAudioBytes, m_StagedOutputFrames.size() / AC97_OUTPUT_CHANNELS);
	}
}

void AC97Device::ResetOutputStream()
{
	if (m_OutputContext == nullptr || m_OutputSource == 0) {
		return;
	}

	if (alcGetCurrentContext() != m_OutputContext && !alcMakeContextCurrent(m_OutputContext)) {
		return;
	}

	alSourceStop(m_OutputSource);

	ALint queued = 0;
	alGetSourcei(m_OutputSource, AL_BUFFERS_QUEUED, &queued);
	while (queued > 0) {
		ALuint buffer = 0;
		alSourceUnqueueBuffers(m_OutputSource, 1, &buffer);
		if (alGetError() != AL_NO_ERROR) {
			break;
		}
		--queued;
	}

	m_FreeOutputBuffers.assign(m_OutputBuffers.begin(), m_OutputBuffers.end());
	m_OutputBufferBytes.fill(0);
	m_QueuedAudioBytes = 0;
	m_StagedOutputFrames.clear();
	m_LastOutputSourceState = -1;
	m_LoggedPlaybackStartFailure = false;
}

void AC97Device::Begin3DVoiceFrameBatch()
{
	for (auto& [voiceHandle, voiceState] : m_Pending3DVoices) {
		(void)voiceHandle;
		voiceState.active = false;
	}
}

void AC97Device::Submit3DVoiceFrames(uint32_t voiceHandle, uint32_t hrtfEntryIndex, bool sourceStereo,
	const std::array<uint8_t, 4>& hrtfSubmix, const std::array<uint32_t, 4>& hrtfSubmixVolumes,
	uint8_t hrtfHeadroom,
	const int16_t* stereoSamples, size_t frameCount)
{
	if (stereoSamples == nullptr || frameCount == 0) {
		return;
	}

	auto& voiceState = m_Pending3DVoices[voiceHandle];
	voiceState.active = true;
	voiceState.sourceStereo = sourceStereo;
	voiceState.hrtfEntryIndex = hrtfEntryIndex;
	voiceState.hrtfSubmix = hrtfSubmix;
	voiceState.hrtfSubmixVolumes = hrtfSubmixVolumes;
	voiceState.hrtfHeadroom = hrtfHeadroom;
	voiceState.samples.resize(frameCount * AC97_OUTPUT_CHANNELS);
	std::copy_n(stereoSamples, frameCount * AC97_OUTPUT_CHANNELS, voiceState.samples.begin());
}

void AC97Device::SubmitPCMFrames(const int16_t* samples, size_t frameCount)
{
	if (samples == nullptr || frameCount == 0 || !EnsureOutputDevice()) {
		return;
	}

	ALint processed = 0;
	alGetSourcei(m_OutputSource, AL_BUFFERS_PROCESSED, &processed);
	ALenum sourceError = alGetError();
	if (sourceError != AL_NO_ERROR) {
		EmuLog(LOG_LEVEL::WARNING,
			"AC97 OpenAL failed to query processed buffers: %s (0x%04x)",
			GetOpenALErrorName(sourceError), static_cast<unsigned>(sourceError));
		return;
	}
	while (processed > 0) {
		ALuint buffer = 0;
		alSourceUnqueueBuffers(m_OutputSource, 1, &buffer);
		if (alGetError() != AL_NO_ERROR) {
			break;
		}
		const auto bufferIndexIt = m_OutputBufferIndex.find(buffer);
		if (bufferIndexIt == m_OutputBufferIndex.end()) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 OpenAL unqueued unknown buffer, resetting stream");
			ResetOutputStream();
			return;
		}
		const size_t bufferIndex = bufferIndexIt->second;
		if (m_QueuedAudioBytes < m_OutputBufferBytes[bufferIndex]) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 OpenAL queue byte accounting underflow, resetting counter");
			m_QueuedAudioBytes = 0;
		} else {
			m_QueuedAudioBytes -= m_OutputBufferBytes[bufferIndex];
		}
		m_OutputBufferBytes[bufferIndex] = 0;
		m_FreeOutputBuffers.push_back(buffer);
		--processed;
	}

	const uint32_t pendingBytes = static_cast<uint32_t>(m_StagedOutputFrames.size() * sizeof(int16_t));
	const uint32_t incomingBytes = static_cast<uint32_t>(frameCount * AC97_OUTPUT_BYTES_PER_FRAME);
	const uint64_t totalBufferedBytes = static_cast<uint64_t>(pendingBytes) + incomingBytes + m_QueuedAudioBytes;
	if (totalBufferedBytes > AC97_MAX_BUFFERED_AUDIO_BYTES) {
		if (!m_LoggedQueueFull) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 OpenAL staged audio backlog full, dropping PCM frames");
			m_LoggedQueueFull = true;
		}
		return;
	}

	const uint16_t masterVolume = ReadRegister16(AC97_Master_Volume);
	const uint16_t pcmOutVolume = ReadRegister16(AC97_PCM_Out_Volume);
	const float leftGain = DecodeOutputAttenuation(masterVolume, true) * DecodeOutputAttenuation(pcmOutVolume, true);
	const float rightGain = DecodeOutputAttenuation(masterVolume, false) * DecodeOutputAttenuation(pcmOutVolume, false);

	const int16_t* output = samples;
	if (leftGain != 1.0f || rightGain != 1.0f) {
		m_OutputScratch.resize(frameCount * AC97_OUTPUT_CHANNELS);
		for (size_t frame = 0; frame < frameCount; ++frame) {
			const size_t sampleIndex = frame * AC97_OUTPUT_CHANNELS;
			m_OutputScratch[sampleIndex] = ClampToInt16(ScaleSample(samples[sampleIndex], leftGain));
			m_OutputScratch[sampleIndex + 1] = ClampToInt16(ScaleSample(samples[sampleIndex + 1], rightGain));
		}
		output = m_OutputScratch.data();
	}

	m_StagedOutputFrames.insert(m_StagedOutputFrames.end(), output, output + (frameCount * AC97_OUTPUT_CHANNELS));

	while (!m_FreeOutputBuffers.empty() && !m_StagedOutputFrames.empty()) {
		const size_t stagedFrameCount = m_StagedOutputFrames.size() / AC97_OUTPUT_CHANNELS;
		const size_t minimumFramesToQueue = (m_QueuedAudioBytes == 0)
			? static_cast<size_t>(AC97_STARTUP_BUFFER_FRAMES)
			: static_cast<size_t>(AC97_STREAM_BUFFER_FRAMES);
		if (stagedFrameCount < minimumFramesToQueue) {
			break;
		}

		const size_t framesToQueue = std::min(stagedFrameCount, static_cast<size_t>(AC97_STREAM_BUFFER_FRAMES));
		const uint32_t queuedBytes = static_cast<uint32_t>(framesToQueue * AC97_OUTPUT_BYTES_PER_FRAME);
		if (static_cast<uint64_t>(m_QueuedAudioBytes) + queuedBytes > AC97_MAX_QUEUED_AUDIO_BYTES) {
			break;
		}

		const ALuint buffer = m_FreeOutputBuffers.back();
		m_FreeOutputBuffers.pop_back();
		const auto bufferIndexIt = m_OutputBufferIndex.find(buffer);
		if (bufferIndexIt == m_OutputBufferIndex.end()) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 OpenAL selected unknown free buffer, dropping PCM frames");
			m_FreeOutputBuffers.push_back(buffer);
			return;
		}
		const size_t bufferIndex = bufferIndexIt->second;
		alBufferData(buffer, AL_FORMAT_STEREO16, m_StagedOutputFrames.data(),
			static_cast<ALsizei>(queuedBytes),
			static_cast<ALsizei>(APU_TIMER_FREQUENCY));
		sourceError = alGetError();
		if (sourceError != AL_NO_ERROR) {
			LogOutputOperationFailure("alBufferData", sourceError);
			m_OutputBufferBytes[bufferIndex] = 0;
			m_FreeOutputBuffers.push_back(buffer);
			return;
		}

		alSourceQueueBuffers(m_OutputSource, 1, &buffer);
		sourceError = alGetError();
		if (sourceError != AL_NO_ERROR) {
			LogOutputOperationFailure("alSourceQueueBuffers", sourceError);
			m_OutputBufferBytes[bufferIndex] = 0;
			m_FreeOutputBuffers.push_back(buffer);
			return;
		}

		m_OutputBufferBytes[bufferIndex] = queuedBytes;
		m_QueuedAudioBytes += queuedBytes;
		m_StagedOutputFrames.erase(m_StagedOutputFrames.begin(),
			m_StagedOutputFrames.begin() + (framesToQueue * AC97_OUTPUT_CHANNELS));
	}

	if (m_StagedOutputFrames.empty()) {
		if (!m_FreeOutputBuffers.empty() && m_QueuedAudioBytes < AC97_MAX_QUEUED_AUDIO_BYTES) {
			m_LoggedQueueFull = false;
		}
	} else if (m_QueuedAudioBytes >= AC97_MAX_QUEUED_AUDIO_BYTES) {
		if (!m_LoggedQueueFull) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 OpenAL queue full, dropping PCM frames");
			m_LoggedQueueFull = true;
		}
	} else if (m_FreeOutputBuffers.empty()) {
		if (!m_LoggedQueueFull) {
			EmuLog(LOG_LEVEL::WARNING, "AC97 OpenAL buffer pool exhausted, dropping PCM frames");
			m_LoggedQueueFull = true;
		}
	} else {
		m_LoggedQueueFull = false;
	}

	ALint state = 0;
	ALint queued = 0;
	ALint processedBuffers = 0;
	if (!QueryOutputSourceSnapshot(state, queued, processedBuffers, sourceError)) {
		EmuLog(LOG_LEVEL::WARNING,
			"AC97 OpenAL failed to query source snapshot: %s (0x%04x)",
			GetOpenALErrorName(sourceError), static_cast<unsigned>(sourceError));
		return;
	}
	if (state != m_LastOutputSourceState) {
		EmuLog(LOG_LEVEL::INFO,
			"AC97 OpenAL source state -> %s, queuedBuffers=%d processedBuffers=%d queuedBytes=%u freeBuffers=%zu stagedFrames=%zu",
			GetOpenALSourceStateName(state), queued, processedBuffers,
			m_QueuedAudioBytes, m_FreeOutputBuffers.size(),
			m_StagedOutputFrames.size() / AC97_OUTPUT_CHANNELS);
		m_LastOutputSourceState = state;
	}
	if (state != AL_PLAYING && m_QueuedAudioBytes != 0) {
		EmuLog(LOG_LEVEL::INFO,
			"AC97 OpenAL requesting playback start from %s, queuedBuffers=%d processedBuffers=%d queuedBytes=%u",
			GetOpenALSourceStateName(state), queued, processedBuffers, m_QueuedAudioBytes);
		alSourcePlay(m_OutputSource);
		sourceError = alGetError();
		ALint postState = state;
		ALint postQueued = queued;
		ALint postProcessed = processedBuffers;
		ALenum postQueryError = AL_NO_ERROR;
		if (!QueryOutputSourceSnapshot(postState, postQueued, postProcessed, postQueryError)) {
			EmuLog(LOG_LEVEL::WARNING,
				"AC97 OpenAL failed to query source snapshot after alSourcePlay: %s (0x%04x)",
				GetOpenALErrorName(postQueryError), static_cast<unsigned>(postQueryError));
			m_LoggedPlaybackStartFailure = true;
			return;
		}
		if (postState != m_LastOutputSourceState) {
			EmuLog(LOG_LEVEL::INFO,
				"AC97 OpenAL source state -> %s, queuedBuffers=%d processedBuffers=%d queuedBytes=%u freeBuffers=%zu stagedFrames=%zu",
				GetOpenALSourceStateName(postState), postQueued, postProcessed,
				m_QueuedAudioBytes, m_FreeOutputBuffers.size(),
				m_StagedOutputFrames.size() / AC97_OUTPUT_CHANNELS);
			m_LastOutputSourceState = postState;
		}
		if (sourceError != AL_NO_ERROR || postState != AL_PLAYING) {
			if (!m_LoggedPlaybackStartFailure) {
				EmuLog(LOG_LEVEL::WARNING,
					"AC97 OpenAL source did not report AL_PLAYING after alSourcePlay: playResult=%s (0x%04x), postPlayState=%s queuedBuffers=%d processedBuffers=%d queuedBytes=%u",
					GetOpenALErrorName(sourceError), static_cast<unsigned>(sourceError),
					GetOpenALSourceStateName(postState), postQueued, postProcessed, m_QueuedAudioBytes);
				m_LoggedPlaybackStartFailure = true;
			}
		} else {
			EmuLog(LOG_LEVEL::INFO,
				"AC97 OpenAL source entered AL_PLAYING, queuedBuffers=%d processedBuffers=%d queuedBytes=%u",
				postQueued, postProcessed, m_QueuedAudioBytes);
			m_LoggedPlaybackStartFailure = false;
		}
	} else if (state == AL_PLAYING) {
		m_LoggedPlaybackStartFailure = false;
	}
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
						const size_t channelIndex = ChannelIndex(channelBase);
						const uint8_t control = static_cast<uint8_t>(ReadRegister(AC97_NAM_SIZE + channelBase + BM_CR, sizeof(uint8_t)));
						m_ChannelAdvanceOnRestart[channelIndex] = false;
						m_ChannelQueuedAfterHalt[channelIndex] = false;
						if ((control & CR_RPBM) == 0 &&
							IsDescriptorErrorAcknowledged(channelBase)) {
							m_ChannelDescriptorError[channelIndex] = false;
						}
						UpdateBusMasterStatus(channelBase);
					}
					return;
				case BM_CIV:
				case BM_PICB:
				case BM_PIV:
					return;
				case BM_LVI:
				{
					const size_t channelIndex = ChannelIndex(channelBase);
					const uint32_t baseAddr = AC97_NAM_SIZE + channelBase;
					const uint8_t previousLastValid = static_cast<uint8_t>(ReadRegister(baseAddr + BM_LVI, sizeof(uint8_t)) & 0x1F);
					const uint8_t newLastValid = static_cast<uint8_t>(value & 0x1F);
					const uint8_t control = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CR, sizeof(uint8_t)));
					WriteRegister(addr, newLastValid, sizeof(uint8_t));
					if ((control & CR_RPBM) == 0 &&
						IsDescriptorErrorAcknowledged(channelBase)) {
						m_ChannelDescriptorError[channelIndex] = false;
					}

					const uint16_t status = ReadRegister16(baseAddr + BM_SR);
					const uint16_t remaining = ReadRegister16(baseAddr + BM_PICB);
					const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CIV, sizeof(uint8_t)) & 0x1F);
					const bool haltedAtEndOfList =
						remaining == 0 &&
						currentIndex == previousLastValid &&
						m_ChannelAdvanceOnRestart[channelIndex];
					if (haltedAtEndOfList && newLastValid != previousLastValid) {
						m_ChannelQueuedAfterHalt[channelIndex] = true;
					}
					if ((control & CR_RPBM) != 0 &&
						(status & SR_DCH) != 0 &&
						haltedAtEndOfList &&
						m_ChannelQueuedAfterHalt[channelIndex] &&
						newLastValid != currentIndex) {
						const uint8_t nextIndex = static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
						WriteRegister(baseAddr + BM_CIV, nextIndex, sizeof(uint8_t));
						WriteRegister16(baseAddr + BM_SR, status & ~(SR_DCH | SR_CELV));
						m_ChannelAdvanceOnRestart[channelIndex] = false;
						m_ChannelQueuedAfterHalt[channelIndex] = false;
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
							IsDescriptorErrorAcknowledged(channelBase)) {
							m_ChannelDescriptorError[channelIndex] = false;
						}
						if ((previousControl & CR_RPBM) == 0 &&
							(control & CR_RPBM) != 0 &&
							m_ChannelAdvanceOnRestart[channelIndex]) {
							const uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CIV, sizeof(uint8_t)) & 0x1F);
							const uint8_t lastValidIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_LVI, sizeof(uint8_t)) & 0x1F);
							// If software stopped the engine after it halted on the previous LVI
							// and then extended LVI, restart from the newly queued descriptor.
							if (ReadRegister16(baseAddr + BM_PICB) == 0 &&
								(currentIndex != lastValidIndex || m_ChannelQueuedAfterHalt[channelIndex])) {
								const uint8_t nextIndex = static_cast<uint8_t>((currentIndex + 1) & (AC97_DESCRIPTOR_COUNT - 1));
								WriteRegister(baseAddr + BM_CIV, nextIndex, sizeof(uint8_t));
								m_ChannelAdvanceOnRestart[channelIndex] = false;
								m_ChannelQueuedAfterHalt[channelIndex] = false;
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

bool AC97Device::IsDescriptorErrorAcknowledged(uint32_t channelBase) const
{
	return (ReadRegister16(AC97_NAM_SIZE + channelBase + BM_SR) & SR_FIFOE) == 0;
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
	m_ChannelQueuedAfterHalt[channelIndex] = false;
	m_ChannelDescriptorError[channelIndex] = false;

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
	uint8_t currentIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_CIV, sizeof(uint8_t)) & 0x1F);
	const uint8_t lastValidIndex = static_cast<uint8_t>(ReadRegister(baseAddr + BM_LVI, sizeof(uint8_t)) & 0x1F);
	const uint32_t descriptorBase = ReadRegister(baseAddr + BM_BDBAR, sizeof(uint32_t)) & ~0x7u;
	if (descriptorBase == 0) {
		return PrimeResult::DescriptorError;
	}

	for (uint32_t descriptorCount = 0; descriptorCount < AC97_DESCRIPTOR_COUNT; ++descriptorCount) {
		uint32_t descriptorControl = 0;
		if (!ReadGuest32(descriptorBase + currentIndex * AC97_DESCRIPTOR_STRIDE + 4, descriptorControl)) {
			return PrimeResult::DescriptorError;
		}

		const uint16_t descriptorLength = static_cast<uint16_t>(descriptorControl & AC97_DESCRIPTOR_LENGTH_MASK);
		if (descriptorLength != 0) {
			WriteRegister(baseAddr + BM_CIV, currentIndex, sizeof(uint8_t));
			WriteRegister16(baseAddr + BM_PICB, descriptorLength);
			WriteRegister16(baseAddr + BM_PIV, GetPrefetchedIndexValue(channelBase, currentIndex, lastValidIndex));
			return PrimeResult::Ready;
		}

		if (currentIndex == lastValidIndex) {
			WriteRegister(baseAddr + BM_CIV, currentIndex, sizeof(uint8_t));
			WriteRegister16(baseAddr + BM_PIV, GetPrefetchedIndexValue(channelBase, currentIndex, lastValidIndex));
			return PrimeResult::EndOfList;
		}

		currentIndex = GetNextDescriptorIndex(currentIndex);
	}

	return PrimeResult::DescriptorError;
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
	const uint8_t initialCurrentIndex = static_cast<uint8_t>(ReadRegister(civAddr, sizeof(uint8_t)) & 0x1F);
	const uint8_t initialLastValidIndex = static_cast<uint8_t>(ReadRegister(lviAddr, sizeof(uint8_t)) & 0x1F);
	const bool haltedAtEndOfList =
		ReadRegister16(picbAddr) == 0 &&
		m_ChannelAdvanceOnRestart[channelIndex] &&
		!m_ChannelQueuedAfterHalt[channelIndex] &&
		initialCurrentIndex == initialLastValidIndex;
	const uint32_t now = GetAPUTime();
	const uint32_t elapsed = now - m_ChannelLastUpdate[channelIndex];
	m_ChannelLastUpdate[channelIndex] = now;

	const auto primeChannel = [&]() {
		switch (PrimeBusMasterChannel(channelBase)) {
		case PrimeResult::Ready:
			m_ChannelQueuedAfterHalt[channelIndex] = false;
			m_ChannelDescriptorError[channelIndex] = false;
			return true;
		case PrimeResult::EndOfList:
			m_ChannelAdvanceOnRestart[channelIndex] = true;
			m_ChannelQueuedAfterHalt[channelIndex] = false;
			m_ChannelDescriptorError[channelIndex] = false;
			status |= SR_LVBCI;
			return false;
		case PrimeResult::DescriptorError:
			m_ChannelAdvanceOnRestart[channelIndex] = false;
			m_ChannelQueuedAfterHalt[channelIndex] = false;
			m_ChannelDescriptorError[channelIndex] = true;
			status |= SR_FIFOE;
			return false;
		default:
			return false;
		}
	};

	if (control & CR_RPBM) {
		uint64_t samplesToConsume = m_ChannelSampleRemainder[channelIndex];
		samplesToConsume += static_cast<uint64_t>(elapsed) * GetBusMasterSampleRate(channelBase);
		m_ChannelSampleRemainder[channelIndex] = static_cast<uint32_t>(samplesToConsume % APU_TIMER_FREQUENCY);
		samplesToConsume /= APU_TIMER_FREQUENCY;

		if (haltedAtEndOfList) {
			samplesToConsume = 0;
		} else if (ReadRegister16(picbAddr) == 0 &&
			(m_ChannelDescriptorError[channelIndex] || !primeChannel())) {
			samplesToConsume = 0;
		}

		while (samplesToConsume > 0) {
			uint16_t remaining = ReadRegister16(picbAddr);
			if (remaining == 0 &&
				(m_ChannelDescriptorError[channelIndex] || !primeChannel())) {
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
				m_ChannelQueuedAfterHalt[channelIndex] = false;
				status |= SR_LVBCI;
				break;
			}

			const uint8_t nextIndex = GetNextDescriptorIndex(currentIndex);
			WriteRegister(civAddr, nextIndex, sizeof(uint8_t));
			m_ChannelAdvanceOnRestart[channelIndex] = false;
			m_ChannelQueuedAfterHalt[channelIndex] = false;
			m_ChannelDescriptorError[channelIndex] = false;
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

	WriteRegister16(AC97_NAM_SIZE + channelBase + BM_PIV, GetPrefetchedIndexValue(channelBase, currentIndex, lastValidIndex));
	WriteRegister16(srAddr, status);
	UpdateGlobalStatus();
}

uint8_t AC97Device::GetPrefetchedIndexValue(uint32_t channelBase, uint8_t currentIndex, uint8_t lastValidIndex) const
{
	// PIV reports the next descriptor the controller can actually prefetch. If
	// later queued entries are zero-length, hardware skips them and exposes the
	// next non-empty descriptor instead. If no later valid entry remains, PIV
	// stays aligned with CIV.
	if (currentIndex == lastValidIndex) {
		return currentIndex;
	}

	const uint32_t descriptorBase = ReadRegister(AC97_NAM_SIZE + channelBase + BM_BDBAR, sizeof(uint32_t)) & ~0x7u;
	if (descriptorBase == 0) {
		return currentIndex;
	}

	uint8_t nextIndex = GetNextDescriptorIndex(currentIndex);
	for (uint32_t descriptorCount = 0; descriptorCount < AC97_DESCRIPTOR_COUNT - 1; ++descriptorCount) {
		uint32_t descriptorControl = 0;
		if (!ReadGuest32(descriptorBase + nextIndex * AC97_DESCRIPTOR_STRIDE + 4, descriptorControl)) {
			return currentIndex;
		}

		if ((descriptorControl & AC97_DESCRIPTOR_LENGTH_MASK) != 0) {
			return nextIndex;
		}

		if (nextIndex == lastValidIndex) {
			return currentIndex;
		}

		nextIndex = GetNextDescriptorIndex(nextIndex);
	}

	return currentIndex;
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
