#ifndef CXBXR_AUDIO_DIAGNOSTICS_H
#define CXBXR_AUDIO_DIAGNOSTICS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

namespace audio_diagnostics {

#ifndef CXBXR_ENABLE_AUDIO_DIAGNOSTIC_LOGGING
#define CXBXR_ENABLE_AUDIO_DIAGNOSTIC_LOGGING 1
#endif

inline constexpr bool kEnableDiagnosticLogging = CXBXR_ENABLE_AUDIO_DIAGNOSTIC_LOGGING != 0;
inline constexpr uint32_t kInt16AbsoluteMax = static_cast<uint32_t>(INT16_MAX) + 1;

inline uint32_t PeakAbsoluteSampleAmplitude(const int16_t* samples, size_t sampleCount)
{
	uint32_t peak = 0;
	for (size_t i = 0; i < sampleCount; ++i) {
		const int32_t signedSample = static_cast<int32_t>(samples[i]);
		const uint32_t magnitude = static_cast<uint32_t>(signedSample < 0 ? -signedSample : signedSample);
		peak = std::max(peak, magnitude);
	}
	return std::min<uint32_t>(peak, kInt16AbsoluteMax);
}

}

#endif
