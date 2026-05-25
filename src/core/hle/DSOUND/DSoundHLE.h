// ******************************************************************
// *
// *  This file is part of the Cxbx project.
// *
// *  (c) 2018-2026 Cxbx-Reloaded contributors
// *
// *  All rights reserved
// *
// ******************************************************************
// Minimal DirectSound HLE bridge — patches DSoundDoWork so the
// game's per-frame audio tick drives our APU/AC97 LLE pipeline.

#pragma once

#include "core\hle\XAPI\Xapi.h"
#include "core\hle\Patches.hpp"
#include "devices\Xbox.h"

// Declared here because our include chain doesn't pull in the Xbox
// kernel headers that define it; the actual value is 0x0 (no bits
// checked), see Intercept.cpp line 83.
extern bool bLLE_APU;

namespace xbox {

void_xt WINAPI EMUPATCH(DirectSoundDoWork)()
{
	// Drive the LLE audio pipeline at the game's frame rate so
	// voice-rendered and DMA-driven audio is consumed.  Do NOT
	// forward to the native implementation — it completes
	// instantly when there is no audio work and creates the
	// infinite NtSetEvent / KeInsertQueueDpc spin.
	if (g_APU != nullptr) {
		g_APU->SynchronizeAudio();
	}
	if (g_AC97 != nullptr) {
		g_AC97->ServiceAudio();
	}
}

} // namespace xbox
