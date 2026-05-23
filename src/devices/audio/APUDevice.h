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
#include <mutex>
#include <vector>

#include "../PCIDevice.h"
class APUDevice : public PCIDevice {
public:
	using PCIDevice::PCIDevice;

	static constexpr size_t MAX_VOICE_HANDLES = 0xFFFF;
	static constexpr size_t MAX_HRTF_VOICES = 64;
	static constexpr size_t HRTF_FILTER_TAPS = 31;
	static constexpr size_t HRTF_FILTER_DELAY_SAMPLES = 42;
	static constexpr size_t HRTF_FILTER_BUFFER_LENGTH = HRTF_FILTER_TAPS + HRTF_FILTER_DELAY_SAMPLES;
	static constexpr size_t MAX_RECENT_FE_METHODS = 16;

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
		size_t previewDecodeFailures = 0;
		bool valid = false;
	};

	struct LowPassFilterState {
		float high = 0.0f;
		float band = 0.0f;
		float low = 0.0f;
	};

	struct HRTFEntryState {
		std::array<std::array<int8_t, HRTF_FILTER_TAPS>, 2> coeffs{};
		int16_t itd = 0;
	};

	struct HRTFFilterChannelState {
		std::array<float, HRTF_FILTER_BUFFER_LENGTH> buf{};
		std::array<float, HRTF_FILTER_TAPS> hrir_coeff_cur{};
		std::array<float, HRTF_FILTER_TAPS> hrir_coeff_tar{};
	};

	struct HRTFFilterState {
		size_t buf_pos = 0;
		std::array<HRTFFilterChannelState, 2> ch{};
		float itd_cur = 0.0f;
		float itd_tar = 0.0f;
	};

	struct RecentFEMethodDiagnostic {
		uint32_t sequence = 0;
		uint32_t addr = 0;
		uint32_t value = 0;
		uint32_t currentVoice = 0;
		uint32_t targetVoice = 0;
		uint32_t feav = 0;
		uint32_t vpvaddr = 0;
		uint32_t vpsgeaddr = 0;
		uint32_t vpssladdr = 0;
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
	void ApplySubmixHeadroom(int32_t* mixBins, size_t frameCount);
	void WriteOutputBuffers(const int32_t* mixBins, size_t frameCount);
	size_t RenderBasicVoiceList(uint32_t topRegister, int32_t* mixBins, size_t frameCount);
	struct BasicVoiceDiagnosticSummary;
	void RenderBasicVoice(uint32_t voiceHandle, int32_t* mixBins, size_t frameCount,
		BasicVoiceDiagnosticSummary* diagnostics = nullptr);
	void LogVoiceTableDiagnostics() const;
	void LogRecentVoiceStateDiagnostics() const;
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
	bool ReadScatterGatherBytes(uint32_t sgeBase, uint32_t maxSge, uint32_t addr, void* dest, size_t size) const;
	bool WriteScatterGatherBytes(uint32_t sgeBase, uint32_t maxSge, uint32_t addr, const void* src, size_t size);
	uint32_t ReadScratchWindowWithDMA(uint32_t sgeBaseRegister, uint32_t maxSgeRegister,
		const uint8_t* data, size_t length, uint32_t addr, unsigned size) const;
	void WriteScratchWindowWithDMA(uint32_t sgeBaseRegister, uint32_t maxSgeRegister,
		uint8_t* data, size_t length, uint32_t addr, uint32_t value, unsigned size);
	uint32_t RefreshFEMemDataRegister(uint32_t fallbackValue);
	bool ResolveOptionalGuestTableBase(uint32_t registerAddress, uint32_t fallbackGuestAddress, uint32_t& guestBase) const;
	void SignalNotifierInterrupt();
	bool ReadVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t& value) const;
	bool WriteVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t value);
	bool WriteVPScatterGatherEntry(uint32_t handle, uint32_t value);
	void WriteNotifierValue(uint32_t voiceHandle, uint32_t notifier, uint32_t value);
	void WriteNotifierStatus(uint32_t voiceHandle, uint32_t notifier, uint8_t status);
	void NotifyVoiceCompletion(uint32_t voiceHandle, uint8_t status);
	uint32_t GetVoicePlaybackOffset(uint32_t voiceHandle) const;
	uint32_t GetVoiceNextHandle(uint32_t voiceHandle) const;
	void SetVoiceNextHandle(uint32_t voiceHandle, uint32_t nextHandle);
	void UnlinkVoiceFromList(uint32_t topRegister, uint32_t voiceHandle);
	void UnlinkVoiceFromLists(uint32_t voiceHandle);
	bool IsVoiceLocked(uint32_t voiceHandle) const;
	void SetVoiceLocked(uint32_t voiceHandle, bool locked);
	bool IsVoiceActiveHinted(uint32_t voiceHandle) const;
	void SetVoiceActiveHint(uint32_t voiceHandle, bool active);
	bool ResolveVoiceAddress(uint32_t linearAddress, uint32_t& guestAddress) const;
	bool ReadVoiceBufferBytes(uint32_t linearAddress, void* dest, size_t size) const;
	bool ReadGuestCircularBuffer(uint32_t guestAddress, uint32_t length, uint32_t& cursor,
		void* dest, size_t size) const;
	bool WriteGuestCircularBuffer(uint32_t guestAddress, uint32_t length, uint32_t& cursor,
		const void* src, size_t size);
	bool HasGuestDspExecution() const;
	bool HasGuestVPOutputBufferPlaybackPath() const;
	bool MixGuestVPOutputBuffers(int16_t* output, size_t frameCount, std::array<uint32_t, 4>* slotPeak);
	bool SubmitGuestVPOutputBuffersToAC97(size_t frameCount, std::array<uint32_t, 4>* slotPeak);
	uint32_t ReadMemoryWindow(const uint8_t* data, size_t length, uint32_t addr, unsigned size) const;
	void WriteMemoryWindow(uint8_t* data, size_t length, uint32_t addr, uint32_t value, unsigned size);
	void WriteHRTFCoefficient(uint32_t entryIndex, size_t channel, size_t coefficientIndex, int8_t value);
	void ClearHRTFFilterState(uint32_t voiceHandle);
	void SetHRTFFilterTarget(uint32_t voiceHandle, const HRTFEntryState& entry);
	void ProcessHRTFSample(uint32_t voiceHandle, float& sampleLeft, float& sampleRight);
	void RecordRecentFEMethod(uint32_t addr, uint32_t value, uint32_t currentVoiceValue);
	void LogRecentFEMethodDiagnostics() const;

	uint32_t ReadRegister(uint32_t addr, unsigned size) const;
	void WriteRegister(uint32_t addr, uint32_t value, unsigned size);
	void SetRegister32(uint32_t addr, uint32_t value);
	uint32_t GetRegister32(uint32_t addr) const;

	std::array<uint8_t, APU_SIZE> m_Registers{};
	mutable std::mutex m_AudioUpdateMutex{};
	uint32_t m_VPFifoLevel = 0;
	uint32_t m_VPFifoLastUpdate = 0;
	uint32_t m_LastAudioUpdate = 0;
	uint32_t m_VPInputSgeHandle = 0;
	uint32_t m_VPOutputSgeHandle = 0;
	uint32_t m_VPNotifyContextDMA = 0;
	uint32_t m_VPCurrentSSLContextDMA = 0;
	uint32_t m_VPSSLBasePage = 0;
	uint32_t m_VPCurrentHRTFEntry = 0;
	std::array<uint8_t, 0x1000 * sizeof(uint32_t)> m_GPXMem{};
	std::array<uint8_t, 0x400 * sizeof(uint32_t)> m_GPMixBuf{};
	std::array<uint8_t, 0x800 * sizeof(uint32_t)> m_GPYMem{};
	std::array<uint8_t, 0x1000 * sizeof(uint32_t)> m_GPPMem{};
	std::array<uint8_t, 0x0C00 * sizeof(uint32_t)> m_EPXMem{};
	std::array<uint8_t, 0x0100 * sizeof(uint32_t)> m_EPYMem{};
	std::array<uint8_t, 0x1000 * sizeof(uint32_t)> m_EPPMem{};
	std::array<HRTFEntryState, 128> m_VPHRTFEntries{};
	std::array<uint8_t, 4> m_VPHRTFSubmix{};
	uint8_t m_VPHRTFHeadroom = 0;
	std::array<uint8_t, 32> m_VPSubmixHeadroom{};
	std::array<uint64_t, (MAX_VOICE_HANDLES + 63) / 64> m_VPVoiceLocked{};
	std::array<uint64_t, (MAX_VOICE_HANDLES + 63) / 64> m_VPActiveVoiceHints{};
	std::array<uint32_t, 4> m_VPOutBufferCursor{};
	std::array<uint32_t, 4> m_VPOutBufferPlaybackCursor{};
	std::array<SSLData, MAX_VOICE_HANDLES> m_VPSSLData{};
	std::array<PlaybackState, MAX_VOICE_HANDLES> m_VPPlaybackState{};
	std::array<std::array<LowPassFilterState, 2>, MAX_VOICE_HANDLES> m_VPLowPassState{};
	std::array<HRTFFilterState, MAX_HRTF_VOICES> m_VPHRTFFilterState{};
	std::array<RecentFEMethodDiagnostic, MAX_RECENT_FE_METHODS> m_RecentFEMethods{};
	size_t m_RecentFEMethodCount = 0;
	size_t m_RecentFEMethodNext = 0;
	uint32_t m_RecentFEMethodSequence = 0;
	std::vector<int16_t> m_VP3DVoiceCaptureScratch{};
	bool m_LoggedXADPCMDecodeFailure = false;
	bool m_LoggedEmptyVoiceTableDiagnostics = false;
	mutable bool m_LoggedVoiceTableReadFailure = false;
	bool m_LoggedVoiceTableWriteFailure = false;
	bool m_LoggedScatterGatherWriteFailure = false;
	bool m_LoggedVoiceListInsertFailure = false;
	bool m_LoggedMissingVoiceTableDuringRender = false;
	bool m_LoggedAC97Missing = false;
	bool m_LoggedStreamingSSLFailure = false;
	bool m_EnableHostSpatialHandoff = true;
	bool m_LoggedVPOutputBufferReadFailure = false;
	bool m_LoggedFallbackActiveVoiceRender = false;
	size_t m_ChunkCaptured3DVoiceCount = 0;
	size_t m_ChunkSubmittedHostSpatialVoiceCount = 0;
};

#endif
