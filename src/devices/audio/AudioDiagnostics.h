#ifndef CXBXR_AUDIO_DIAGNOSTICS_H
#define CXBXR_AUDIO_DIAGNOSTICS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

namespace audio_diagnostics {

inline constexpr bool kEnableDiagnosticLogging = true;

inline uint16_t PeakAbsoluteSampleAmplitude(const int16_t* samples, size_t sampleCount)
{
	uint32_t peak = 0;
	for (size_t i = 0; i < sampleCount; ++i) {
		const uint32_t magnitude = static_cast<uint32_t>(std::abs(static_cast<int32_t>(samples[i])));
		peak = std::max(peak, magnitude);
	}
	return static_cast<uint16_t>(std::min<uint32_t>(peak, static_cast<uint32_t>(INT16_MAX) + 1));
}

}

#endif
