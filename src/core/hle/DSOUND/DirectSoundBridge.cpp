// ******************************************************************
// *
// *  This file is part of the Cxbx project.
// *
// *  (c) 2018-2026 Cxbx-Reloaded contributors
// *
// *  All rights reserved
// *
// ******************************************************************
// Minimal DirectSound HLE bridge that patches the Xbox DSound entry
// points so that game audio triggers the LLE APU/AC97 pipeline at the
// cadence the game expects, rather than relying solely on the VBlank
// DPC timer.

#define LOG_PREFIX CXBXR_MODULE::DSOUND

#include "core\kernel\exports\xboxkrnl.h"
#include "core\kernel\init\CxbxKrnl.h"
#include "common\Logging.h"
#include "core\kernel\support\Emu.h"
#include "core\hle\XAPI\Xapi.h"
#include "core\hle\Patches.hpp"
#include "devices\Xbox.h"

namespace xbox {

// ******************************************************************
// * patch: DirectSoundDoWork
// ******************************************************************
// The game calls DirectSoundDoWork every frame to drive its audio
// pipeline (voice dispatch, buffer completion, stream processing).
// We intercept the call so that our LLE APU/AC97 backend synchronises
// with every game frame instead of waiting for the VBlank timer.
void_xt WINAPI EMUPATCH(DirectSoundDoWork)()
{
	// Flush APU voice rendering and AC97 DMA before the native
	// DoWork inspects voice positions and notifier state.
	if (g_APU != nullptr) {
		g_APU->SynchronizeAudio();
	}
	if (g_AC97 != nullptr) {
		g_AC97->ServiceAudio();
	}

	// Forward to the original Xbox kernel DirectSoundDoWork so that
	// guest-side DSound bookkeeping (voice chains, buffer-completion,
	// stream-packet dispatch) still executes.
	void_xt (WINAPI * nativeDoWork)() = nullptr;
	nativeDoWork = reinterpret_cast<decltype(nativeDoWork)>(
		GetPatchedFunctionTrampoline("DirectSoundDoWork"));
	if (nativeDoWork != nullptr) {
		nativeDoWork();
	}
}

} // namespace xbox
