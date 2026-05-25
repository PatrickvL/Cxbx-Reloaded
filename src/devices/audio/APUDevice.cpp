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
#include "AC97Device.h"
#include "AudioDiagnostics.h"
#include "APUTimer.h"
#include "dsp/dsp.h"
#include "dsp/dsp_state.h"
#include "common/AddressRanges.h"
#include "common/audio/XADPCM.h"
#include "core/kernel/exports/EmuKrnl.h"
#include "core/kernel/support/Emu.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <vector>

#define LOG_PREFIX CXBXR_MODULE::MCPX

namespace {

constexpr uint32_t APU_VP_BASE = 0x20000;
constexpr uint32_t APU_VP_SIZE = 0x10000;
constexpr uint32_t APU_VP_FREE = 0x10;
constexpr uint32_t APU_VP_FIFO_CAPACITY = 0x80;
constexpr uint32_t APU_VP_VOICE_MAX_HANDLE = 0xFFFF;

constexpr uint32_t APU_GP_BASE = 0x30000;
constexpr uint32_t APU_GP_SIZE = 0x10000;
constexpr uint32_t NV_PAPU_GPXMEM = 0x00000000;
constexpr uint32_t NV_PAPU_GPMIXBUF = 0x00005000;
constexpr uint32_t NV_PAPU_GPYMEM = 0x00006000;
constexpr uint32_t NV_PAPU_GPPMEM = 0x0000A000;
constexpr uint32_t NV_PAPU_GPRST = 0x0000FFFC;

constexpr uint32_t APU_EP_BASE = 0x50000;
constexpr uint32_t APU_EP_SIZE = 0x10000;
constexpr uint32_t NV_PAPU_EPXMEM = 0x00000000;
constexpr uint32_t NV_PAPU_EPYMEM = 0x00006000;
constexpr uint32_t NV_PAPU_EPPMEM = 0x0000A000;
constexpr uint32_t NV_PAPU_EPRST = 0x0000FFFC;

constexpr uint32_t NV_PAPU_ISTS = 0x00001000;
constexpr uint32_t NV_PAPU_ISTS_GINTSTS = 1 << 0;
constexpr uint32_t NV_PAPU_ISTS_FETINTSTS = 1 << 4;
constexpr uint32_t NV_PAPU_ISTS_FENINTSTS = 1 << 5;
constexpr uint32_t NV_PAPU_ISTS_FEVINTSTS = 1 << 6;
constexpr uint32_t NV_PAPU_IEN = 0x00001004;
constexpr uint32_t NV_PAPU_FECTL = 0x00001100;
constexpr uint32_t NV_PAPU_FECTL_FEMETHMODE = 0x000000E0;
constexpr uint32_t NV_PAPU_FECTL_FEMETHMODE_TRAPPED = 0x000000E0;
constexpr uint32_t NV_PAPU_FECTL_FETRAPREASON = 0x00000F00;
constexpr uint32_t NV_PAPU_FECTL_FETRAPREASON_REQUESTED = 0x00000F00;
constexpr uint32_t NV_PAPU_FECV = 0x00001110;
constexpr uint32_t NV_PAPU_FEAV = 0x00001118;
constexpr uint32_t NV_PAPU_FENADDR = 0x0000115C;
constexpr uint32_t NV_PAPU_FEDECMETH = 0x00001300;
constexpr uint32_t NV_PAPU_FEDECPARAM = 0x00001304;
constexpr uint32_t NV_PAPU_FEMEMADDR = 0x00001324;
constexpr uint32_t NV_PAPU_FEMEMDATA = 0x00001334;
constexpr uint32_t NV_PAPU_FETFORCE0 = 0x00001500;
constexpr uint32_t NV_PAPU_FETFORCE1 = 0x00001504;
constexpr uint32_t NV_PAPU_FETFORCE1_SE2FE_IDLE_VOICE = 1 << 15;
constexpr uint32_t NV_PAPU_SECTL = 0x00002000;
constexpr uint32_t NV_PAPU_SECTL_XCNTMODE = 0x00000018;
constexpr uint32_t NV_PAPU_SECTL_XCNTMODE_OFF = 0;
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

constexpr uint32_t NV_PAPU_GPRST_GPRST = 1 << 0;
constexpr uint32_t NV_PAPU_GPRST_GPDSPRST = 1 << 1;
constexpr uint32_t NV_PAPU_EPRST_EPRST = 1 << 0;
constexpr uint32_t NV_PAPU_EPRST_EPDSPRST = NV_PAPU_GPRST_GPDSPRST;

constexpr uint32_t NV_PAPU_GPOFBASE0 = 0x00003024;
constexpr uint32_t NV_PAPU_GPOFEND0 = 0x00003028;
constexpr uint32_t NV_PAPU_GPOFCUR0 = 0x0000302C;
constexpr uint32_t NV_PAPU_GPIFBASE0 = 0x00003064;
constexpr uint32_t NV_PAPU_GPIFEND0 = 0x00003068;
constexpr uint32_t NV_PAPU_GPIFCUR0 = 0x0000306C;
constexpr uint32_t NV_PAPU_EPOFBASE0 = 0x00004024;
constexpr uint32_t NV_PAPU_EPOFEND0 = 0x00004028;
constexpr uint32_t NV_PAPU_EPOFCUR0 = 0x0000402C;
constexpr uint32_t NV_PAPU_EPIFBASE0 = 0x00004064;
constexpr uint32_t NV_PAPU_EPIFEND0 = 0x00004068;
constexpr uint32_t NV_PAPU_EPIFCUR0 = 0x0000406C;
constexpr uint32_t NV_PAPU_FIFO_VALUE = 0x00FFFFFF;

constexpr size_t APU_DSP_FRAME_SAMPLES = 32;
constexpr size_t APU_DSP_FRAME_STEREO_SAMPLES = APU_DSP_FRAME_SAMPLES * 2;
constexpr size_t APU_DSP_EP_FRAME_DIVIDER = 8;
constexpr size_t APU_DSP_GP_OUTPUT_FIFO_COUNT = 4;
constexpr size_t APU_DSP_GP_INPUT_FIFO_COUNT = 2;
constexpr size_t APU_DSP_EP_OUTPUT_FIFO_COUNT = 4;
constexpr size_t APU_DSP_EP_INPUT_FIFO_COUNT = 2;
constexpr uint32_t APU_DSP_GP_MIXBUF_BASE = 0x001400;

constexpr uint32_t NV_PAPU_FEAV_VALUE = 0x0000FFFF;
constexpr uint32_t NV_PAPU_FEAV_LST = 0x00030000;

constexpr uint32_t NV1BA0_PIO_SET_ANTECEDENT_VOICE = 0x00000120;
constexpr uint32_t NV1BA0_PIO_VOICE_ON = 0x00000124;
constexpr uint32_t NV1BA0_PIO_VOICE_OFF = 0x00000128;
constexpr uint32_t NV1BA0_PIO_VOICE_RELEASE = 0x0000012C;
constexpr uint32_t NV1BA0_PIO_GET_VOICE_POSITION = 0x00000130;
constexpr uint32_t NV1BA0_PIO_VOICE_PAUSE = 0x00000140;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY = 0x00000160;
constexpr uint32_t NV1BA0_PIO_SET_CONTEXT_DMA_NOTIFY = 0x00000180;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_SSL_CONTEXT_DMA = 0x0000018C;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_SSL = 0x00000190;
constexpr uint32_t NV1BA0_PIO_SET_SUBMIX_HEADROOM = 0x00000200;
constexpr uint32_t NV1BA0_PIO_SET_HRTF_HEADROOM = 0x00000280;
constexpr uint32_t NV1BA0_PIO_SET_HRTF_SUBMIXES = 0x000002C0;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_VOICE = 0x000002F8;
constexpr uint32_t NV1BA0_PIO_VOICE_LOCK = 0x000002FC;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_VBIN = 0x00000300;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_FMT = 0x00000304;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENV0 = 0x00000308;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENVA = 0x0000030C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENV1 = 0x00000310;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_ENVF = 0x00000314;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_MISC = 0x00000318;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_HRTF = 0x0000031C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_SSL_A = 0x00000320;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_SSL_B = 0x0000035C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_VOLA = 0x00000360;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_VOLB = 0x00000364;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_VOLC = 0x00000368;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_LFO_ENV = 0x0000036C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_LFO_MOD = 0x00000370;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_FCA = 0x00000374;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_FCB = 0x00000378;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_PITCH = 0x0000037C;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_BUF_BASE = 0x000003A0;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_BUF_LBO = 0x000003A4;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_BUF_CBO = 0x000003D8;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_CFG_BUF_EBO = 0x000003DC;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_METHOD_FIRST = 0x00000300;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_METHOD_LAST = 0x000003DC;
constexpr uint32_t NV1BA0_PIO_SET_HRIR = 0x00000400;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_X = 0x0000043C;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE = 0x00000804;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE_OFFSET = 0x00000808;
constexpr uint32_t NV1BA0_PIO_SET_OUTBUF_BA = 0x00001000;
constexpr uint32_t NV1BA0_PIO_SET_OUTBUF_LEN = 0x00001004;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE = 0x00001800;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_OFFSET = 0x00001808;
constexpr uint32_t SE2FE_IDLE_VOICE = 0x00008000;

constexpr uint32_t NV1BA0_PIO_VOICE_ON_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_ON_ENVF = 0x0F000000;
constexpr uint32_t NV1BA0_PIO_VOICE_ON_ENVA = 0xF0000000;
constexpr uint32_t NV1BA0_PIO_VOICE_OFF_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_RELEASE_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_GET_VOICE_POSITION_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_PAUSE_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_VOICE_PAUSE_ACTION = 1 << 18;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_SSL_BASE_PAGE = 0x003FFFC0;
constexpr uint32_t NV1BA0_PIO_SET_SUBMIX_HEADROOM_AMOUNT = 0x00000007;
constexpr uint32_t NV1BA0_PIO_SET_HRTF_HEADROOM_AMOUNT = 0x00000007;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_HRTF_HANDLE = 0x0000FFFF;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_SSL_A_COUNT = 0x000000FF;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_SSL_A_BASE = 0xFFFFFF00;
constexpr uint32_t NV1BA0_PIO_SET_VOICE_TAR_PITCH_STEP = 0xFFFF0000;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_LEFT0 = 0x000000FF;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_RIGHT0 = 0x0000FF00;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_LEFT1 = 0x00FF0000;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_RIGHT1 = 0xFF000000;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_X_LEFT30 = 0x000000FF;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_X_RIGHT30 = 0x0000FF00;
constexpr uint32_t NV1BA0_PIO_SET_HRIR_X_ITD = 0xFFFF0000;
constexpr uint32_t NV1BA0_PIO_SET_SSL_SEGMENT_OFFSET = 0x00000600;
constexpr uint32_t NV1BA0_PIO_SET_SSL_SEGMENT_LENGTH = 0x00000604;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE_HANDLE = 0xFFFFFFFF;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_INBUF_SGE_OFFSET_PARAMETER = 0xFFFFF000;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_HANDLE = 0xFFFFFFFF;
constexpr uint32_t NV1BA0_PIO_SET_CURRENT_OUTBUF_SGE_OFFSET_PARAMETER = 0xFFFFF000;
constexpr uint32_t NV1BA0_PIO_SET_OUTBUF_BA_ADDRESS = 0x007FFF00;
constexpr uint32_t NV1BA0_PIO_SET_OUTBUF_LEN_VALUE = 0x007FFF00;

const char* GetAPURegisterTraceName(uint32_t addr)
{
	switch (addr) {
	case NV_PAPU_VPVADDR: return "NV_PAPU_VPVADDR";
	case NV_PAPU_VPSGEADDR: return "NV_PAPU_VPSGEADDR";
	case NV_PAPU_VPSSLADDR: return "NV_PAPU_VPSSLADDR";
	case NV_PAPU_TVL2D: return "NV_PAPU_TVL2D";
	case NV_PAPU_TVL3D: return "NV_PAPU_TVL3D";
	case NV_PAPU_TVLMP: return "NV_PAPU_TVLMP";
	default: return nullptr;
	}
}

const char* GetAPUVPMethodTraceName(uint32_t addr)
{
	switch (addr) {
	case NV1BA0_PIO_SET_ANTECEDENT_VOICE: return "NV1BA0_PIO_SET_ANTECEDENT_VOICE";
	case NV1BA0_PIO_VOICE_ON: return "NV1BA0_PIO_VOICE_ON";
	case NV1BA0_PIO_VOICE_OFF: return "NV1BA0_PIO_VOICE_OFF";
	case NV1BA0_PIO_VOICE_RELEASE: return "NV1BA0_PIO_VOICE_RELEASE";
	case NV1BA0_PIO_GET_VOICE_POSITION: return "NV1BA0_PIO_GET_VOICE_POSITION";
	case NV1BA0_PIO_VOICE_PAUSE: return "NV1BA0_PIO_VOICE_PAUSE";
	case NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY: return "NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY";
	case NV1BA0_PIO_SET_CONTEXT_DMA_NOTIFY: return "NV1BA0_PIO_SET_CONTEXT_DMA_NOTIFY";
	case NV1BA0_PIO_SET_CURRENT_SSL_CONTEXT_DMA: return "NV1BA0_PIO_SET_CURRENT_SSL_CONTEXT_DMA";
	case NV1BA0_PIO_SET_CURRENT_SSL: return "NV1BA0_PIO_SET_CURRENT_SSL";
	case NV1BA0_PIO_SET_SUBMIX_HEADROOM: return "NV1BA0_PIO_SET_SUBMIX_HEADROOM";
	case NV1BA0_PIO_SET_HRTF_HEADROOM: return "NV1BA0_PIO_SET_HRTF_HEADROOM";
	case NV1BA0_PIO_SET_HRTF_SUBMIXES: return "NV1BA0_PIO_SET_HRTF_SUBMIXES";
	case NV1BA0_PIO_SET_CURRENT_VOICE: return "NV1BA0_PIO_SET_CURRENT_VOICE";
	case NV1BA0_PIO_VOICE_LOCK: return "NV1BA0_PIO_VOICE_LOCK";
	case NV1BA0_PIO_SET_VOICE_LFO_MOD: return "NV1BA0_PIO_SET_VOICE_LFO_MOD";
	default:
		if (addr >= NV1BA0_PIO_SET_VOICE_METHOD_FIRST && addr <= NV1BA0_PIO_SET_VOICE_METHOD_LAST) {
			return "NV1BA0_PIO_SET_VOICE_*";
		}
		return nullptr;
	}
}

constexpr uint32_t NV_PAVS_SIZE = 0x00000080;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN = 0x00000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN_V0BIN = 0x0000001F;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN_V1BIN = 0x000003E0;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN_V2BIN = 0x00007C00;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN_V3BIN = 0x001F0000;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN_V4BIN = 0x03E00000;
constexpr uint32_t NV_PAVS_VOICE_CFG_VBIN_V5BIN = 0x7C000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT = 0x00000004;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_V6BIN = 0x0000001F;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_V7BIN = 0x000003E0;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_SAMPLES_PER_BLOCK = 0x001F0000;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_MULTIPASS_BIN = 0x001F0000;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_MULTIPASS = 1 << 21;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_PERSIST = 1 << 23;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_DATA_TYPE = 1 << 24;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_LOOP = 1 << 25;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_CLEAR_MIX = 1 << 26;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_STEREO = 1 << 27;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE = 0x30000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE = 0xC0000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_U8 = 0;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S16 = 1;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S24 = 2;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S32 = 3;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B8 = 0;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B16 = 1;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_ADPCM = 2;
constexpr uint32_t NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B32 = 3;
constexpr uint32_t NV_PAVS_VOICE_CFG_HRTF_TARGET_HANDLE = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV0 = 0x00000008;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV0_EA_ATTACKRATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV0_EA_DELAYTIME = 0x00FFF000;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVA = 0x0000000C;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVA_EA_DECAYRATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVA_EA_HOLDTIME = 0x00FFF000;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVA_EA_SUSTAINLEVEL = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV1 = 0x00000010;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV1_EF_ATTACKRATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENV1_EF_DELAYTIME = 0x00FFF000;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVF = 0x00000014;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVF_EF_DECAYRATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVF_EF_HOLDTIME = 0x00FFF000;
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVF_EF_SUSTAINLEVEL = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC = 0x00000018;
constexpr uint32_t NV_PAVS_VOICE_CFG_HRTF_TARGET = 0x0000001C;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSL_START = 0x00000020;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSH_SAMPLE = 0x00000024;
constexpr uint32_t NV_PAVS_VOICE_CUR_ECNT = 0x00000034;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE = 0x00000054;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET = 0x00000058;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT = 0x0000005C;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLA = 0x00000060;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLB = 0x00000064;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLC = 0x00000068;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_ENV = 0x0000006C;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK = 0x0000007C;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLA_VOLUME6_B3_0 = 0x0000000F;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLA_VOLUME0 = 0x0000FFF0;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLA_VOLUME7_B3_0 = 0x000F0000;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLA_VOLUME1 = 0xFFF00000;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLB_VOLUME6_B7_4 = 0x0000000F;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLB_VOLUME2 = 0x0000FFF0;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLB_VOLUME7_B7_4 = 0x000F0000;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLB_VOLUME3 = 0xFFF00000;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLC_VOLUME6_B11_8 = 0x0000000F;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLC_VOLUME4 = 0x0000FFF0;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLC_VOLUME7_B11_8 = 0x000F0000;
constexpr uint32_t NV_PAVS_VOICE_TAR_VOLC_VOLUME5 = 0xFFF00000;

constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_PAUSED = 1 << 18;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_NEW_VOICE = 1 << 20;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE = 1 << 21;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR = 0x0F000000;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF = 0;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_DELAY = 1;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK = 2;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD = 3;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY = 4;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN = 5;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_RELEASE = 6;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EFCUR_FORCE_RELEASE = 7;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EACUR = 0xF0000000;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_EACUR_OFF = 0;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSL_START_BA = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_ECNT_EACOUNT = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_ECNT_EFCOUNT = 0xFFFF0000;
constexpr uint32_t NV_PAVS_VOICE_PAR_LFO = 0x00000050;
constexpr uint32_t NV_PAVS_VOICE_PAR_LFO_LFOALVL = 0x00007FFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_LFO_LFOADR = 0x00008000;
constexpr uint32_t NV_PAVS_VOICE_PAR_LFO_LFOFLVL = 0x7FFF0000;
constexpr uint32_t NV_PAVS_VOICE_PAR_LFO_LFOFDR = 0x80000000;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_LFOA_DELAYMODE = 1 << 16;
constexpr uint32_t NV_PAVS_VOICE_PAR_STATE_LFOF_DELAYMODE = 1 << 17;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET_CBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET_EALVL = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT_EBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT_EFLVL = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH = 0xFFFF0000;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_ENV_EA_RELEASERATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_ENV_LFOADLT = 0x003FF000;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_ENV_LFOFDLT = 0xFFC00000;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_MOD = 0x00000070;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_MOD_LFOAAM = 0x000000FF;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_MOD_LFOAFM = 0x0000FF00;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_MOD_LFOAFC = 0x00FF0000;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_MOD_LFOFFM = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC_LFOA_DELAYMODE = 0x00004000;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC_LFOF_DELAYMODE = 0x00008000;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC_FMODE = 0x00030000;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCA = 0x00000074;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCB = 0x00000078;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCA_FC0 = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCA_FC1 = 0xFFFF0000;

constexpr uint32_t APU_VOICE_LIST_INHERIT = 0;
constexpr uint32_t APU_SGE_PAGE_SIZE = 0x1000;
constexpr size_t APU_AUDIO_CHUNK_FRAMES = 256;
constexpr float APU_VOLUME_DECIBEL_DIVISOR = 64.0f * -20.0f;
// MCPX notifier records are 16-byte entries. The hardware stores two generic
// notifier entries before the per-voice audio notifiers, so
// MCPX_HW_NOTIFIER_BASE_OFFSET means 2 entries * 16 bytes = 32 bytes rather
// than a raw byte offset literal.
constexpr uint32_t MCPX_HW_NOTIFIER_ENTRY_SIZE = 16;
constexpr uint32_t MCPX_HW_NOTIFIER_BASE_OFFSET = 2;
// MCPX exposes four notifier slots per voice in the audio notifier table.
constexpr uint32_t MCPX_HW_NOTIFIER_COUNT = 4;
constexpr uint32_t MCPX_HW_NOTIFIER_SSLA_DONE = 0;
constexpr uint32_t MCPX_HW_NOTIFIER_SSLB_DONE = 1;
constexpr uint32_t MCPX_HW_NOTIFIER_VOICE_POSITION = 2;
// NV1BA0 reports successful notifier completion with status 0x01 in observed
// guest/hardware traces; using 0xFF caused guest polling loops to wait
// indefinitely because it does not match the expected success code.
constexpr uint8_t NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS = 0x01;
constexpr uint8_t NV1BA0_NOTIFICATION_STATUS_DONE_ERROR = 0x80;
constexpr size_t APU_MAX_CONSECUTIVE_PREVIEW_DECODE_FAILURES = 8;
constexpr double APU_PITCH_STEP_EXPONENT = 4096.0;
constexpr size_t APU_XADPCM_PCM_SAMPLES_PER_BLOCK = XBOX_ADPCM_DSTSIZE / sizeof(int16_t);
constexpr size_t APU_XADPCM_MAX_CHANNELS = 2;
constexpr size_t APU_XADPCM_MAX_SOURCE_BLOCK_BYTES = XBOX_ADPCM_SRCSIZE * APU_XADPCM_MAX_CHANNELS;
constexpr size_t APU_XADPCM_MAX_DECODED_SAMPLES = APU_XADPCM_PCM_SAMPLES_PER_BLOCK * APU_XADPCM_MAX_CHANNELS;
constexpr size_t APU_MIXBIN_COUNT = 32;
constexpr uint32_t APU_FIRST_NON_STEREO_BIN = 2;
constexpr uint32_t APU_MAX_3D_VOICES = static_cast<uint32_t>(APUDevice::MAX_HRTF_VOICES);
constexpr size_t APU_HRTF_SUBMIX_COUNT = 4;
static_assert(APU_HRTF_SUBMIX_COUNT <= 8, "HRTF submix handoff expects no more than eight voice volumes");
constexpr size_t APU_HRTF_ENTRY_COUNT = 128;
constexpr size_t APU_HRTF_COEFFICIENT_COUNT = APUDevice::HRTF_FILTER_TAPS;
constexpr uint32_t APU_INVALID_HRTF_ENTRY_INDEX = 0xFFFF;
constexpr uint32_t APU_IRQ = 5;
constexpr float APU_HRTF_ITD_SCALE = 512.0f;
constexpr float APU_HRTF_PAN_ITD_WEIGHT = 0.75f;
constexpr float APU_HRTF_PAN_MAGNITUDE_WEIGHT = 0.25f;
constexpr float APU_HRTF_PARAM_SMOOTH_ALPHA = 0.01f;
constexpr float APU_HRTF_NORMALIZATION_EPSILON = 0.000001f;
constexpr float APU_HRTF_MAX_DELAY_SAMPLES_FLOAT = static_cast<float>(APUDevice::HRTF_FILTER_DELAY_SAMPLES);
constexpr float APU_HRTF_ITD_NORMALIZER = APU_HRTF_ITD_SCALE * APU_HRTF_MAX_DELAY_SAMPLES_FLOAT;
// Scale normalized floating-point samples to signed 16-bit PCM amplitude.
constexpr float APU_SAMPLE_SCALE_FACTOR = 32767.0f;
constexpr size_t APU_DIAGNOSTIC_MAX_VOICES_TO_LOG = 4;
constexpr size_t APU_DIAGNOSTIC_MAX_ACTIVE_VOICES_TO_LOG = 4;

uint32_t GetFEMethodTargetVoiceOrDefault(uint32_t addr, uint32_t value, uint32_t currentVoiceValue)
{
	switch (addr) {
	case NV1BA0_PIO_SET_CURRENT_VOICE:
		m_CurrentVoice = value & NV1BA0_PIO_VOICE_ON_HANDLE;
		SetRegister32(NV_PAPU_FECV, value & NV1BA0_PIO_VOICE_ON_HANDLE);
		return;
	}

	const int32_t maxLevel = static_cast<int32_t>(APU_LFO_LEVEL_MAX);
	const int32_t period = maxLevel * 2;
	int32_t value = static_cast<int32_t>(std::min<uint32_t>(level, APU_LFO_LEVEL_MAX));
	value += descending ? -static_cast<int32_t>(delta) : static_cast<int32_t>(delta);
	// Reflect off the 0..0x7FFF bounds and flip the direction bit so the guest-
	// visible PAR_LFO state advances as a continuous triangle waveform; modulo
	// arithmetic keeps large deltas from iterating across multiple bounces.
	value %= period;
	if (value < 0) {
		value += period;
	}
	if (value >= maxLevel) {
		value = period - value;
		descending = true;
	} else {
		descending = false;
	}

	level = static_cast<uint32_t>(value);
}

bool IsVoiceLFODelayActive(uint32_t voiceState, uint32_t delayMask, uint32_t envelopeStateMask)
{
	return (voiceState & delayMask) != 0 &&
		GetMaskedValue(voiceState, envelopeStateMask) == NV_PAVS_VOICE_PAR_STATE_EFCUR_DELAY;
}

float DecodeSignedLFOAmount(uint32_t value)
{
	return static_cast<float>(static_cast<int8_t>(value & 0xFF)) / APU_LFO_MODULATION_NORMALIZER;
}

uint32_t ExtractLFOField(uint32_t value, uint32_t mask)
{
	return (value & mask) >> Ctz32(mask);
}

int16_t ClampToInt16(int32_t value)
{
	if (value > 32767) {
		return 32767;
	}
	if (value < -32768) {
		return -32768;
	}
	return static_cast<int16_t>(value);
}

int16_t ConvertFloatSampleToInt16(float sample)
{
	return ClampToInt16(static_cast<int32_t>(std::lrint(
		static_cast<double>(ClampUnitSample(sample)) * 32767.0)));
}

uint32_t ConvertInt16ToDSP24(int16_t sample)
{
	return static_cast<uint32_t>(static_cast<int32_t>(sample) << 8) & 0x00FFFFFF;
}

int16_t ConvertDSP24ToInt16(uint32_t sample)
{
	const int32_t signedSample = static_cast<int32_t>(sample << 8) >> 8;
	return ClampToInt16(signedSample >> 8);
}

}

struct APUDevice::BasicVoiceDiagnosticSummary {
	uint32_t voiceHandle = 0;
	uint32_t bins[8]{};
	uint32_t volumes[8]{};
	uint8_t headroom[8]{};
	uint32_t startOffset = 0;
	uint32_t offsetAdvance = 0;
	size_t framesRendered = 0;
	double pitchStep = 0.0;
	double maxPitchStep = 0.0;
	float maxEnvelopeGain = 0.0f;
	float maxFilterEnvelopeGain = 0.0f;
	float maxAmplitudeLFOModulation = 1.0f;
	uint32_t decodedPeak = 0;
	uint32_t mixedPeak = 0;
	bool visited = false;
	bool active = false;
	bool paused = false;
	bool mixed = false;
	bool decodedNonZero = false;
	bool stereoContribution = false;
	bool nonStereoContribution = false;
};

extern AC97Device* g_AC97;

// Basic VP playback and guest-visible buffer plumbing exist here, but full
// GP/EP DSP execution and threaded audio scheduling are still incomplete.

void APUDevice::GPDspScratchRW(void* opaque, uint8_t* ptr, uint32_t addr, size_t len, bool dir)
{
	auto* apu = static_cast<APUDevice*>(opaque);
	if (apu != nullptr) {
		apu->TransferDSPScratch(true, ptr, addr, len, dir);
	}
}

void APUDevice::EPDspScratchRW(void* opaque, uint8_t* ptr, uint32_t addr, size_t len, bool dir)
{
	auto* apu = static_cast<APUDevice*>(opaque);
	if (apu != nullptr) {
		apu->TransferDSPScratch(false, ptr, addr, len, dir);
	}
}

void APUDevice::GPDspFifoRW(void* opaque, uint8_t* ptr, unsigned index, size_t len, bool dir)
{
	auto* apu = static_cast<APUDevice*>(opaque);
	if (apu != nullptr) {
		apu->TransferDSPFifo(true, ptr, index, len, dir);
	}
}

void APUDevice::EPDspFifoRW(void* opaque, uint8_t* ptr, unsigned index, size_t len, bool dir)
{
	auto* apu = static_cast<APUDevice*>(opaque);
	if (apu != nullptr) {
		apu->TransferDSPFifo(false, ptr, index, len, dir);
	}
}

void APUDevice::InitializeDSP()
{
	if (m_GPDsp == nullptr) {
		m_GPDsp = dsp_init(this, &APUDevice::GPDspScratchRW, &APUDevice::GPDspFifoRW);
	}
	if (m_EPDsp == nullptr) {
		m_EPDsp = dsp_init(this, &APUDevice::EPDspScratchRW, &APUDevice::EPDspFifoRW);
	}
}

void APUDevice::ResetDSPState()
{
	InitializeDSP();
	if (m_GPDsp != nullptr) {
		dsp_reset(m_GPDsp);
	}
	if (m_EPDsp != nullptr) {
		dsp_reset(m_EPDsp);
	}
	m_DSPFrameDivider = 0;
	m_DSPOutputScratch.clear();
	m_LoggedDSPOutputCaptureFailure = false;
}

bool APUDevice::IsGPDSPEnabled() const
{
	const uint32_t reset = GetRegister32(APU_GP_BASE + NV_PAPU_GPRST);
	return (reset & NV_PAPU_GPRST_GPRST) != 0 &&
		(reset & NV_PAPU_GPRST_GPDSPRST) != 0 &&
		m_GPDsp != nullptr;
}

bool APUDevice::IsEPDSPEnabled() const
{
	const uint32_t reset = GetRegister32(APU_EP_BASE + NV_PAPU_EPRST);
	return (reset & NV_PAPU_EPRST_EPRST) != 0 &&
		(reset & NV_PAPU_EPRST_EPDSPRST) != 0 &&
		m_EPDsp != nullptr;
}

bool APUDevice::IsAnyDSPEnabled() const
{
	return IsGPDSPEnabled() || IsEPDSPEnabled();
}

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
	m_LastAudioUpdate = m_VPFifoLastUpdate;
	m_XGSCounter = 0;
	m_VPInputSgeHandle = 0;
	m_VPOutputSgeHandle = 0;
	m_VPNotifyContextDMA = 0;
	m_VPCurrentSSLContextDMA = 0;
	m_VPSSLBasePage = 0;
	m_VPCurrentHRTFEntry = 0;
	m_GPXMem.fill(0);
	m_GPMixBuf.fill(0);
	m_GPYMem.fill(0);
	m_GPPMem.fill(0);
	m_EPXMem.fill(0);
	m_EPYMem.fill(0);
	m_EPPMem.fill(0);
	m_VPHRTFEntries.fill(HRTFEntryState{});
	m_VPHRTFSubmix.fill(0);
	m_VPHRTFHeadroom = 0;
	m_VPSubmixHeadroom.fill(0);
	m_VPVoiceLocked.fill(0);
	// Initialize the FE/VP fallback state consistently on every reset by
	// clearing the voice hints and shadow table before the guest provides
	// NV_PAPU_VPVADDR.
	m_VPActiveVoiceHints.fill(0);
	m_VPVoiceTableShadow.fill(0);
	m_VPOutBufferCursor.fill(0);
	m_VPOutBufferPlaybackCursor.fill(0);
	m_VPOutBufferQueuedBytes.fill(0);
	m_VPSSLData.fill(APUDevice::SSLData{});
	m_VPPlaybackState.fill(APUDevice::PlaybackState{});
	m_VPNotifierEnvelopeState.fill(NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF);
	m_VPHRTFFilterState.fill(APUDevice::HRTFFilterState{});
	m_RecentFEMethods.fill(APUDevice::RecentFEMethodDiagnostic{});
	m_RecentFEMethodCount = 0;
	m_RecentFEMethodNext = 0;
	m_RecentFEMethodSequence = 0;
	m_LoggedXADPCMDecodeFailure = false;
	m_LoggedEmptyVoiceTableDiagnostics = false;
	m_LoggedVoiceTableReadFailure = false;
	m_LoggedVoiceTableWriteFailure = false;
	m_LoggedScatterGatherWriteFailure = false;
	m_LoggedMissingVoiceTableDuringRender = false;
	m_EnableHostSpatialHandoff = true;
	m_LoggedVPOutputBufferReadFailure = false;
	m_LoggedVPOutputBufferUnderrun = false;
	m_LoggedVPOutputBufferOverrun = false;
	ResetDSPState();
	m_LoggedFallbackActiveVoiceRender = false;
	m_ChunkCaptured3DVoiceCount = 0;
	m_ChunkSubmittedHostSpatialVoiceCount = 0;

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
	SetRegister32(NV_PAPU_SECTL, 0x00000008);
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
	SetRegister32(NV_PAPU_TVL2D, APU_VP_VOICE_MAX_HANDLE);
	SetRegister32(NV_PAPU_TVL3D, APU_VP_VOICE_MAX_HANDLE);
	SetRegister32(NV_PAPU_TVLMP, APU_VP_VOICE_MAX_HANDLE);
	RefreshVPStatus();
	RefreshInterruptStatus();
}

uint32_t APUDevice::IORead(int barIndex, uint32_t addr, unsigned size)
{
	return MMIORead(barIndex, addr, size);
}

void APUDevice::IOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	MMIOWrite(barIndex, addr, value, size);
}

uint32_t APUDevice::MMIORead(int barIndex, uint32_t addr, unsigned size)
{
	(void)barIndex;
	SynchronizeAudio();

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
		// XGSCNT reflects guest-visible rendered sample progress rather than raw
		// host time; SynchronizeAudio only advances this counter while XCNTMODE
		// allows audio progress, so MMIO reads expose the frozen value directly.
		return ReadRegisterFragment(m_XGSCounter, addr - NV_PAPU_XGSCNT, size);
	}

	if (addr >= NV_PAPU_FEMEMDATA && addr < NV_PAPU_FEMEMDATA + sizeof(uint32_t)) {
		const uint32_t currentValue = RefreshFEMemDataRegister(0);
		return ReadRegisterFragment(currentValue, addr - NV_PAPU_FEMEMDATA, size);
	}

	if (addr >= NV_PAPU_ISTS && addr < NV_PAPU_ISTS + sizeof(uint32_t)) {
		RefreshInterruptStatus();
		return ReadRegisterFragment(GetRegister32(NV_PAPU_ISTS), addr - NV_PAPU_ISTS, size);
	}

	return ReadRegister(addr, size);
}

void APUDevice::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	(void)barIndex;
	SynchronizeAudio();

	// Diagnostic: log every distinct register write so we can confirm
	// the guest is actually configuring the APU.
	{
		static std::once_flag once;
		static uint32_t lastLogFrame;
		static uint32_t frameCounter;
		std::call_once(once, []{ lastLogFrame = GetAPUTime(); });
		++frameCounter;
		if (GetAPUTime() - lastLogFrame >= 48000) { // every ~1 sec
			lastLogFrame = GetAPUTime();
			EmuLog(LOG_LEVEL::INFO,
				"APU diag: MMIO writes/sec=%u VPVADDR=0x%08X SECTL=0x%08X TVL2D=0x%04X TVL3D=0x%04X activeVoiceHints=%u",
				frameCounter,
				GetRegister32(NV_PAPU_VPVADDR),
				GetRegister32(NV_PAPU_SECTL),
				GetRegister32(NV_PAPU_TVL2D),
				GetRegister32(NV_PAPU_TVL3D),
				m_VPActiveVoiceHints[0] != 0 ? 1 : 0);
			frameCounter = 0;
		}
	}

	// Diagnostic: log the first N VP-method writes with address/value
	if (addr >= APU_VP_BASE && addr < APU_VP_BASE + APU_VP_SIZE) {
		static uint32_t vpDiagIdx;
		if (vpDiagIdx < 100) {
			++vpDiagIdx;
			const uint32_t vpOffset = addr - APU_VP_BASE;
			EmuLog(LOG_LEVEL::INFO,
				"APU diag: VP write #%u VPoff=0x%04X value=0x%08X size=%u %s",
				vpDiagIdx, vpOffset, value, size,
				(vpOffset == 0x124) ? "*** VOICE_ON ***" :
				(vpOffset == 0x10C) ? "CLEAR_VOICES" :
				(vpOffset == 0x2F8) ? "SET_CURRENT_VOICE" :
				(vpOffset == 0x304) ? "SET_VOICE_CFG_FMT" :
				(vpOffset == 0x318) ? "SET_VOICE_CFG_VBIN" :
				(vpOffset == 0x31C) ? "SET_VOICE_CFG_VBOUT" : "");
		}
		// Always log VOICE_ON even beyond the diag limit
		if ((addr - APU_VP_BASE) == 0x124) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU diag: VOICE_ON detected! handle=0x%04X VPVADDR=0x%08X",
				value & 0xFFFF,
				GetRegister32(NV_PAPU_VPVADDR));
		}
	}

	// Diagnostic: also log first GP FIFO writes with more detail
	if (addr >= APU_GP_BASE && addr < APU_GP_BASE + APU_GP_SIZE) {
		static uint32_t gpDiagIdx;
		if (gpDiagIdx < 25) {
			++gpDiagIdx;
			const uint32_t gpOffset = addr - APU_GP_BASE;
			EmuLog(LOG_LEVEL::INFO,
				"APU diag: GP write #%u GPoff=0x%05X value=0x%08X size=%u",
				gpDiagIdx, gpOffset, value, size);
		}
	}

	// Diagnostic: log first general-region writes (non-VP, non-GP, non-EP)
	if (!(addr >= APU_VP_BASE && addr < APU_VP_BASE + APU_VP_SIZE) &&
	    !(addr >= APU_GP_BASE && addr < APU_GP_BASE + APU_GP_SIZE) &&
	    !(addr >= APU_EP_BASE && addr < APU_EP_BASE + APU_EP_SIZE)) {
		static uint32_t genDiagIdx;
		if (genDiagIdx < 15) {
			++genDiagIdx;
			EmuLog(LOG_LEVEL::INFO,
				"APU diag: general write #%u addr=0x%05X value=0x%08X size=%u",
				genDiagIdx, addr, value, size);
		}
	}

	if (addr >= APU_VP_BASE && addr < APU_VP_BASE + APU_VP_SIZE) {
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			const uint32_t methodAddr = addr - APU_VP_BASE;
			if (const char* name = GetAPUVPMethodTraceName(methodAddr)) {
				EmuLog(LOG_LEVEL::INFO,
					"APU MMIO VP method write %s method=0x%08x value=0x%08x size=%u",
					name,
					methodAddr,
					value,
					size);
			}
		}
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
		RefreshInterruptStatus();
		return;
	}

	if ((addr >= NV_PAPU_IEN && addr < NV_PAPU_IEN + sizeof(uint32_t)) ||
		(addr >= NV_PAPU_FECTL && addr < NV_PAPU_FECTL + sizeof(uint32_t))) {
		WriteRegister(addr, value, size);
		RefreshInterruptStatus();
		return;
	}

	if (addr >= NV_PAPU_XGSCNT && addr < NV_PAPU_XGSCNT + sizeof(uint32_t)) {
		// Treat XGSCNT as read-only: guests can probe/reset it without clobbering
		// the live counter value returned by MMIO reads.
		return;
	}

	if ((addr >= NV_PAPU_FETFORCE0 && addr < NV_PAPU_FETFORCE0 + sizeof(uint32_t)) ||
		(addr >= NV_PAPU_FETFORCE1 && addr < NV_PAPU_FETFORCE1 + sizeof(uint32_t))) {
		WriteRegister(addr, value, size);
		const bool wroteFETFORCE0 =
			addr >= NV_PAPU_FETFORCE0 && addr < NV_PAPU_FETFORCE0 + sizeof(uint32_t);
		uint32_t fectl = GetRegister32(NV_PAPU_FECTL);
		const bool idleVoiceTrapPending =
			(fectl & NV_PAPU_FECTL_FEMETHMODE) == NV_PAPU_FECTL_FEMETHMODE_TRAPPED &&
			(fectl & NV_PAPU_FECTL_FETRAPREASON) == NV_PAPU_FECTL_FETRAPREASON_REQUESTED;
		if (idleVoiceTrapPending) {
			const bool clearIdleVoiceTrap =
				wroteFETFORCE0 ||
				((GetRegister32(NV_PAPU_FETFORCE1) & NV_PAPU_FETFORCE1_SE2FE_IDLE_VOICE) == 0);
			if (clearIdleVoiceTrap) {
				fectl &= ~(NV_PAPU_FECTL_FEMETHMODE | NV_PAPU_FECTL_FETRAPREASON);
				SetRegister32(NV_PAPU_FECTL, fectl);
			}
		}
		RefreshInterruptStatus();
		return;
	}

	if (addr >= NV_PAPU_FEMEMDATA && addr < NV_PAPU_FEMEMDATA + sizeof(uint32_t)) {
		WriteRegister(addr, value, size);
		WriteGuestWord(GetRegister32(NV_PAPU_FEMEMADDR), GetRegister32(NV_PAPU_FEMEMDATA));
		RefreshFEMemDataRegister(GetRegister32(NV_PAPU_FEMEMDATA));
		return;
	}

	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		if (size == sizeof(uint32_t)) {
			if (const char* name = GetAPURegisterTraceName(addr)) {
				EmuLog(LOG_LEVEL::INFO,
					"APU MMIO register write %s old=0x%08x new=0x%08x",
					name,
					GetRegister32(addr),
					value);
				switch (addr) {
				case NV_PAPU_VPVADDR:
					m_LoggedVoiceTableReadFailure = false;
					m_LoggedVoiceTableWriteFailure = false;
					break;
				case NV_PAPU_VPSGEADDR:
					m_LoggedScatterGatherWriteFailure = false;
					break;
				default:
					break;
				}
			}
		}
	}

	WriteRegister(addr, value, size);
	const uint32_t registerBase = addr & ~0x3u;
	if (IsAPUWordAddressRegister(registerBase)) {
		SetRegister32(registerBase, NormalizeAPUWordAddress(GetRegister32(registerBase)));
		if (registerBase == NV_PAPU_FEMEMADDR) {
			RefreshFEMemDataRegister(GetRegister32(NV_PAPU_FEMEMDATA));
		}
	} else if (IsAPUMaxSgeRegister(registerBase)) {
		SetRegister32(registerBase, NormalizeAPUMaxSge(GetRegister32(registerBase)));
	}
}


uint32_t APUDevice::GPRead(uint32_t addr, unsigned size)
{
	if (addr >= NV_PAPU_GPXMEM && addr < NV_PAPU_GPXMEM + m_GPXMem.size()) {
		const uint32_t wordAddr = (addr - NV_PAPU_GPXMEM) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_GPDsp, 'X', wordAddr), addr & 0x3u, size);
	}
	if (addr >= NV_PAPU_GPMIXBUF && addr < NV_PAPU_GPMIXBUF + m_GPMixBuf.size()) {
		const uint32_t wordAddr = APU_DSP_GP_MIXBUF_BASE + (addr - NV_PAPU_GPMIXBUF) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_GPDsp, 'X', wordAddr), addr & 0x3u, size);
	}
	if (addr >= NV_PAPU_GPYMEM && addr < NV_PAPU_GPYMEM + m_GPYMem.size()) {
		const uint32_t wordAddr = (addr - NV_PAPU_GPYMEM) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_GPDsp, 'Y', wordAddr), addr & 0x3u, size);
	}
	if (addr >= NV_PAPU_GPPMEM && addr < NV_PAPU_GPPMEM + m_GPPMem.size()) {
		const uint32_t wordAddr = (addr - NV_PAPU_GPPMEM) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_GPDsp, 'P', wordAddr), addr & 0x3u, size);
	}
	return ReadRegister(APU_GP_BASE + addr, size);
}

void APUDevice::GPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	if (addr >= NV_PAPU_GPXMEM && addr < NV_PAPU_GPXMEM + m_GPXMem.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_GPXMEM) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_GPXMEM) & 0x3u;
		const uint32_t wordAddr = wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_GPDsp, 'X', wordAddr);
		dsp_write_memory(m_GPDsp, 'X', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	if (addr >= NV_PAPU_GPMIXBUF && addr < NV_PAPU_GPMIXBUF + m_GPMixBuf.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_GPMIXBUF) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_GPMIXBUF) & 0x3u;
		const uint32_t wordAddr = APU_DSP_GP_MIXBUF_BASE + wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_GPDsp, 'X', wordAddr);
		dsp_write_memory(m_GPDsp, 'X', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	if (addr >= NV_PAPU_GPYMEM && addr < NV_PAPU_GPYMEM + m_GPYMem.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_GPYMEM) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_GPYMEM) & 0x3u;
		const uint32_t wordAddr = wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_GPDsp, 'Y', wordAddr);
		dsp_write_memory(m_GPDsp, 'Y', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	if (addr >= NV_PAPU_GPPMEM && addr < NV_PAPU_GPPMEM + m_GPPMem.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_GPPMEM) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_GPPMEM) & 0x3u;
		const uint32_t wordAddr = wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_GPDsp, 'P', wordAddr);
		dsp_write_memory(m_GPDsp, 'P', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	const uint32_t oldValue = ReadRegister(APU_GP_BASE + addr, sizeof(uint32_t));
	WriteRegister(APU_GP_BASE + addr, value, size);
	if (addr == NV_PAPU_GPRST && size == sizeof(uint32_t)) {
		const bool wasEnabled = (oldValue & (NV_PAPU_GPRST_GPRST | NV_PAPU_GPRST_GPDSPRST)) ==
			(NV_PAPU_GPRST_GPRST | NV_PAPU_GPRST_GPDSPRST);
		const bool isEnabled = (value & (NV_PAPU_GPRST_GPRST | NV_PAPU_GPRST_GPDSPRST)) ==
			(NV_PAPU_GPRST_GPRST | NV_PAPU_GPRST_GPDSPRST);
		if (!isEnabled) {
			dsp_reset(m_GPDsp);
			m_DSPFrameDivider = 0;
		} else if (!wasEnabled) {
			dsp_bootstrap(m_GPDsp);
		}
	}
}

bool APUDevice::TransferDSPScratch(bool gp, uint8_t* ptr, uint32_t addr, size_t len, bool dir)
{
	if (ptr == nullptr || len == 0) {
		return true;
	}

	const uint32_t sgeBase = GetRegister32(gp ? NV_PAPU_GPSADDR : NV_PAPU_EPSADDR);
	const uint32_t maxSge = GetRegister32(gp ? NV_PAPU_GPSMAXSGE : NV_PAPU_EPSMAXSGE);
	if (sgeBase == 0) {
		if (!dir) {
			std::memset(ptr, 0, len);
		}
		return false;
	}

	return dir
		? WriteScatterGatherBytes(sgeBase, maxSge, addr, ptr, len)
		: ReadScatterGatherBytes(sgeBase, maxSge, addr, ptr, len);
}

uint32_t APUDevice::TransferDSPCircularScatterGather(uint32_t sgeBase, uint32_t maxSge, uint8_t* ptr,
	uint32_t base, uint32_t end, uint32_t cur, size_t len, bool dir)
{
	if (ptr == nullptr || len == 0 || end <= base) {
		return base;
	}
	if (cur >= end) {
		cur = base + ((cur - base) % (end - base));
	} else if (cur < base) {
		cur = base;
	}

	size_t remaining = len;
	uint8_t* bytes = ptr;
	while (remaining != 0) {
		const uint32_t chunk = std::min<uint32_t>(end - cur, static_cast<uint32_t>(remaining));
		const bool ok = dir
			? WriteScatterGatherBytes(sgeBase, maxSge, cur, bytes, chunk)
			: ReadScatterGatherBytes(sgeBase, maxSge, cur, bytes, chunk);
		if (!ok) {
			if (!dir) {
				std::memset(bytes, 0, remaining);
			}
			break;
		}
		bytes += chunk;
		remaining -= chunk;
		cur += chunk;
		if (cur >= end) {
			cur = base;
		}
	}

	return cur;
}

void APUDevice::CaptureEPFifoOutput(uint8_t* ptr, size_t len)
{
	if (ptr == nullptr || len == 0 || (len % sizeof(int16_t)) != 0) {
		if (!m_LoggedDSPOutputCaptureFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU DSP EP output capture rejected payload ptr=%p len=%zu aligned=%d",
				ptr,
				len,
				(len % sizeof(int16_t)) == 0 ? 1 : 0);
			m_LoggedDSPOutputCaptureFailure = true;
		}
		return;
	}

	const size_t sampleCount = len / sizeof(int16_t);
	m_DSPOutputScratch.resize(sampleCount);
	std::memcpy(m_DSPOutputScratch.data(), ptr, len);
	m_LoggedDSPOutputCaptureFailure = false;
}

void APUDevice::TransferDSPFifo(bool gp, uint8_t* ptr, unsigned index, size_t len, bool dir)
{
	if (ptr == nullptr || len == 0) {
		return;
	}

	const uint32_t sgeBase = GetRegister32(gp ? NV_PAPU_GPFADDR : NV_PAPU_EPFADDR);
	const uint32_t maxSge = GetRegister32(gp ? NV_PAPU_GPFMAXSGE : NV_PAPU_EPFMAXSGE);
	if (sgeBase == 0) {
		if (!dir) {
			std::memset(ptr, 0, len);
		}
		return;
	}

	uint32_t baseRegister = 0;
	uint32_t endRegister = 0;
	uint32_t curRegister = 0;
	if (gp) {
		if (dir) {
			if (index >= APU_DSP_GP_OUTPUT_FIFO_COUNT) {
				return;
			}
			baseRegister = NV_PAPU_GPOFBASE0 + static_cast<uint32_t>(index) * 0x10;
			endRegister = NV_PAPU_GPOFEND0 + static_cast<uint32_t>(index) * 0x10;
			curRegister = NV_PAPU_GPOFCUR0 + static_cast<uint32_t>(index) * 0x10;
		} else {
			if (index >= APU_DSP_GP_INPUT_FIFO_COUNT) {
				std::memset(ptr, 0, len);
				return;
			}
			baseRegister = NV_PAPU_GPIFBASE0 + static_cast<uint32_t>(index) * 0x10;
			endRegister = NV_PAPU_GPIFEND0 + static_cast<uint32_t>(index) * 0x10;
			curRegister = NV_PAPU_GPIFCUR0 + static_cast<uint32_t>(index) * 0x10;
		}
	} else {
		if (dir) {
			if (index >= APU_DSP_EP_OUTPUT_FIFO_COUNT) {
				return;
			}
			baseRegister = NV_PAPU_EPOFBASE0 + static_cast<uint32_t>(index) * 0x10;
			endRegister = NV_PAPU_EPOFEND0 + static_cast<uint32_t>(index) * 0x10;
			curRegister = NV_PAPU_EPOFCUR0 + static_cast<uint32_t>(index) * 0x10;
		} else {
			if (index >= APU_DSP_EP_INPUT_FIFO_COUNT) {
				std::memset(ptr, 0, len);
				return;
			}
			baseRegister = NV_PAPU_EPIFBASE0 + static_cast<uint32_t>(index) * 0x10;
			endRegister = NV_PAPU_EPIFEND0 + static_cast<uint32_t>(index) * 0x10;
			curRegister = NV_PAPU_EPIFCUR0 + static_cast<uint32_t>(index) * 0x10;
		}
	}

	const uint32_t base = GetRegister32(baseRegister) & NV_PAPU_FIFO_VALUE;
	const uint32_t end = GetRegister32(endRegister) & NV_PAPU_FIFO_VALUE;
	if (base == 0 || end <= base) {
		if (!dir) {
			std::memset(ptr, 0, len);
		}
		return;
	}

	if (!gp && dir && index == 0) {
		CaptureEPFifoOutput(ptr, len);
	}

	const uint32_t cur = GetRegister32(curRegister) & NV_PAPU_FIFO_VALUE;
	const uint32_t next = TransferDSPCircularScatterGather(sgeBase, maxSge, ptr, base, end, cur, len, dir);
	SetRegister32(curRegister, next & NV_PAPU_FIFO_VALUE);
}

bool APUDevice::ProcessDSPAudio(int16_t* output, const int32_t* mixBins, size_t frameCount)
{
	if (output == nullptr || mixBins == nullptr || frameCount == 0 || !IsAnyDSPEnabled() ||
		(frameCount % APU_DSP_FRAME_SAMPLES) != 0) {
		return false;
	}

	m_DSPOutputScratch.clear();
	const uint64_t maxCycleBudget = std::max<size_t>(frameCount, APU_DSP_FRAME_SAMPLES) * 4096ull;
	const auto runDSP = [maxCycleBudget](DSPState* dsp) {
		dsp_start_frame(dsp);
		dsp->core.is_idle = false;
		dsp->core.cycle_count = 0;
		uint64_t scheduledCycles = 0;
		while (!dsp->core.is_idle && scheduledCycles < maxCycleBudget) {
			dsp_run(dsp, 1000);
			scheduledCycles += 1000;
		}
	};

	for (size_t frameBase = 0; frameBase < frameCount; frameBase += APU_DSP_FRAME_SAMPLES) {
		if (IsGPDSPEnabled()) {
			for (size_t mixbin = 0; mixbin < std::min<size_t>(APU_MIXBIN_COUNT, 32); ++mixbin) {
				for (size_t sample = 0; sample < APU_DSP_FRAME_SAMPLES; ++sample) {
					const int16_t value = ClampToInt16(mixBins[mixbin * frameCount + frameBase + sample]);
					dsp_write_memory(m_GPDsp, 'X',
						APU_DSP_GP_MIXBUF_BASE + static_cast<uint32_t>(mixbin * APU_DSP_FRAME_SAMPLES + sample),
						ConvertInt16ToDSP24(value));
				}
			}
			runDSP(m_GPDsp);
		}

		if (IsEPDSPEnabled()) {
			++m_DSPFrameDivider;
			if ((m_DSPFrameDivider % APU_DSP_EP_FRAME_DIVIDER) == 0) {
				runDSP(m_EPDsp);
			}
		}
	}

	if (m_DSPOutputScratch.size() != frameCount * 2) {
		if (!m_LoggedDSPOutputCaptureFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU DSP frame produced %zu samples, expected %zu for %zu stereo frames; falling back to non-DSP host mix",
				m_DSPOutputScratch.size(),
				frameCount * 2,
				frameCount);
			m_LoggedDSPOutputCaptureFailure = true;
		}
		return false;
	}
	m_LoggedDSPOutputCaptureFailure = false;
	std::copy(m_DSPOutputScratch.begin(), m_DSPOutputScratch.end(), output);
	return true;
}


uint32_t APUDevice::VPRead(uint32_t addr, unsigned size)
{
	UpdateVPFifo();

	if (addr >= APU_VP_FREE && addr < APU_VP_FREE + sizeof(uint32_t)) {
		return ReadRegisterFragment(
			GetVPFifoFreeSlots(),
			addr - APU_VP_FREE,
			size);
	}

	if (addr >= NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY && addr < NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY + sizeof(uint32_t)) {
		return ReadRegisterFragment(m_VPCurrentHRTFEntry, addr - NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY, size);
	}

	if (addr >= NV1BA0_PIO_VOICE_LOCK && addr < NV1BA0_PIO_VOICE_LOCK + sizeof(uint32_t)) {
		const uint32_t currentVoice = GetRegister32(NV_PAPU_FECV) & APU_VP_VOICE_MAX_HANDLE;
		const uint32_t lockValue = IsVoiceLocked(currentVoice) ? 1u : 0u;
		return ReadRegisterFragment(lockValue, addr - NV1BA0_PIO_VOICE_LOCK, size);
	}

	if (addr >= NV1BA0_PIO_SET_CURRENT_SSL && addr < NV1BA0_PIO_SET_CURRENT_SSL + sizeof(uint32_t)) {
		return ReadRegisterFragment(m_VPSSLBasePage, addr - NV1BA0_PIO_SET_CURRENT_SSL, size);
	}

	if (addr >= NV1BA0_PIO_SET_HRTF_SUBMIXES && addr < NV1BA0_PIO_SET_HRTF_SUBMIXES + sizeof(uint32_t)) {
		const uint32_t submixValue = static_cast<uint32_t>(m_VPHRTFSubmix[0])
			| (static_cast<uint32_t>(m_VPHRTFSubmix[1]) << 8)
			| (static_cast<uint32_t>(m_VPHRTFSubmix[2]) << 16)
			| (static_cast<uint32_t>(m_VPHRTFSubmix[3]) << 24);
		return ReadRegisterFragment(submixValue, addr - NV1BA0_PIO_SET_HRTF_SUBMIXES, size);
	}

	if (addr >= NV1BA0_PIO_SET_HRTF_HEADROOM && addr < NV1BA0_PIO_SET_HRTF_HEADROOM + sizeof(uint32_t)) {
		return ReadRegisterFragment(m_VPHRTFHeadroom, addr - NV1BA0_PIO_SET_HRTF_HEADROOM, size);
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
	if (addr >= NV_PAPU_EPXMEM && addr < NV_PAPU_EPXMEM + m_EPXMem.size()) {
		const uint32_t wordAddr = (addr - NV_PAPU_EPXMEM) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_EPDsp, 'X', wordAddr), addr & 0x3u, size);
	}
	if (addr >= NV_PAPU_EPYMEM && addr < NV_PAPU_EPYMEM + m_EPYMem.size()) {
		const uint32_t wordAddr = (addr - NV_PAPU_EPYMEM) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_EPDsp, 'Y', wordAddr), addr & 0x3u, size);
	}
	if (addr >= NV_PAPU_EPPMEM && addr < NV_PAPU_EPPMEM + m_EPPMem.size()) {
		const uint32_t wordAddr = (addr - NV_PAPU_EPPMEM) / sizeof(uint32_t);
		return ReadRegisterFragment(dsp_read_memory(m_EPDsp, 'P', wordAddr), addr & 0x3u, size);
	}
	return ReadRegister(APU_EP_BASE + addr, size);
}

void APUDevice::EPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	if (addr >= NV_PAPU_EPXMEM && addr < NV_PAPU_EPXMEM + m_EPXMem.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_EPXMEM) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_EPXMEM) & 0x3u;
		const uint32_t wordAddr = wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_EPDsp, 'X', wordAddr);
		dsp_write_memory(m_EPDsp, 'X', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	if (addr >= NV_PAPU_EPYMEM && addr < NV_PAPU_EPYMEM + m_EPYMem.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_EPYMEM) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_EPYMEM) & 0x3u;
		const uint32_t wordAddr = wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_EPDsp, 'Y', wordAddr);
		dsp_write_memory(m_EPDsp, 'Y', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	if (addr >= NV_PAPU_EPPMEM && addr < NV_PAPU_EPPMEM + m_EPPMem.size()) {
		const uint32_t wordOffset = (addr - NV_PAPU_EPPMEM) & ~0x3u;
		const uint32_t fragmentOffset = (addr - NV_PAPU_EPPMEM) & 0x3u;
		const uint32_t wordAddr = wordOffset / sizeof(uint32_t);
		const uint32_t current = dsp_read_memory(m_EPDsp, 'P', wordAddr);
		dsp_write_memory(m_EPDsp, 'P', wordAddr, WriteRegisterFragment(current, value, fragmentOffset, size));
		return;
	}
	const uint32_t oldValue = ReadRegister(APU_EP_BASE + addr, sizeof(uint32_t));
	WriteRegister(APU_EP_BASE + addr, value, size);
	if (addr == NV_PAPU_EPRST && size == sizeof(uint32_t)) {
		const bool wasEnabled = (oldValue & (NV_PAPU_EPRST_EPRST | NV_PAPU_EPRST_EPDSPRST)) ==
			(NV_PAPU_EPRST_EPRST | NV_PAPU_EPRST_EPDSPRST);
		const bool isEnabled = (value & (NV_PAPU_EPRST_EPRST | NV_PAPU_EPRST_EPDSPRST)) ==
			(NV_PAPU_EPRST_EPRST | NV_PAPU_EPRST_EPDSPRST);
		if (!isEnabled) {
			dsp_reset(m_EPDsp);
			m_DSPFrameDivider = 0;
		} else if (!wasEnabled) {
			dsp_bootstrap(m_EPDsp);
			m_DSPFrameDivider = 0;
		}
	}
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
	const uint32_t currentVoiceValue = currentVoice();

	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		RecordRecentFEMethod(addr, value, currentVoiceValue);
	}

	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		if (addr >= NV1BA0_PIO_SET_VOICE_METHOD_FIRST && addr <= NV1BA0_PIO_SET_VOICE_METHOD_LAST) {
			EmuLog(LOG_LEVEL::INFO,
				"APU SET_VOICE_* method=0x%08x value=0x%08x voice=0x%04x vpvaddr=0x%08x vpsgeaddr=0x%08x",
				addr,
				value,
				currentVoiceValue & APU_VP_VOICE_MAX_HANDLE,
				GetRegister32(NV_PAPU_VPVADDR),
				GetRegister32(NV_PAPU_VPSGEADDR));
		}
	}

	switch (addr) {
	case NV1BA0_PIO_SET_ANTECEDENT_VOICE:
		SetRegister32(NV_PAPU_FEAV, value);
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			EmuLog(LOG_LEVEL::INFO,
				"APU SET_ANTECEDENT_VOICE value=0x%08x list=%u antecedent=0x%04x",
				value,
				(value & NV_PAPU_FEAV_LST) >> Ctz32(NV_PAPU_FEAV_LST),
				value & NV_PAPU_FEAV_VALUE);
		}
		m_LoggedEmptyVoiceTableDiagnostics = false;
		return;
	case NV1BA0_PIO_SET_CURRENT_VOICE:
		SetRegister32(NV_PAPU_FECV, value);
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			EmuLog(LOG_LEVEL::INFO,
				"APU SET_CURRENT_VOICE value=0x%08x voice=0x%04x vpvaddr=0x%08x vpsgeaddr=0x%08x vpssladdr=0x%08x",
				value,
				value & APU_VP_VOICE_MAX_HANDLE,
				GetRegister32(NV_PAPU_VPVADDR),
				GetRegister32(NV_PAPU_VPSGEADDR),
				GetRegister32(NV_PAPU_VPSSLADDR));
		}
		return;
	case NV1BA0_PIO_VOICE_ON: {
		const uint32_t selectedHandle = value & NV1BA0_PIO_VOICE_ON_HANDLE;
		if (selectedHandle >= APU_VP_VOICE_MAX_HANDLE) {
			return;
		}

		{
			static bool once;
			if (!once) {
				once = true;
				EmuLog(LOG_LEVEL::INFO,
					"APU diag: first VOICE_ON handle=%u SECTL=0x%08X VPVADDR=0x%08X",
					selectedHandle,
					GetRegister32(NV_PAPU_SECTL),
					GetRegister32(NV_PAPU_VPVADDR));
			}
		}

		UnlinkVoiceFromLists(selectedHandle);

		const uint32_t feav = GetRegister32(NV_PAPU_FEAV);
		const uint32_t list = (feav & NV_PAPU_FEAV_LST) >> Ctz32(NV_PAPU_FEAV_LST);
		const uint32_t antecedentVoice = feav & NV_PAPU_FEAV_VALUE;
		uint32_t topRegister = 0;
		uint32_t topBefore = APU_VP_VOICE_MAX_HANDLE;
		bool inserted = false;
		if (list != APU_VOICE_LIST_INHERIT) {
			switch (list) {
			case 1: topRegister = NV_PAPU_TVL2D; break;
			case 2: topRegister = NV_PAPU_TVL3D; break;
			case 3: topRegister = NV_PAPU_TVLMP; break;
			default: break;
			}
			if (topRegister != 0) {
				topBefore = GetRegister32(topRegister);
				WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
					NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE,
					topBefore);
				SetRegister32(topRegister, selectedHandle);
				inserted = true;
			}
		} else {
			if (antecedentVoice < APU_VP_VOICE_MAX_HANDLE) {
				uint32_t nextHandle = 0;
				if (ReadVoiceMask(antecedentVoice, NV_PAVS_VOICE_TAR_PITCH_LINK,
					NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle)) {
					WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
						NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle);
					WriteVoiceMask(antecedentVoice, NV_PAVS_VOICE_TAR_PITCH_LINK,
						NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, selectedHandle);
					inserted = true;
				}
			}
		}
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			EmuLog(LOG_LEVEL::INFO,
				"APU VOICE_ON handle=%u feav=0x%08x list=%u antecedent=0x%04x topRegister=0x%08x topBefore=0x%08x inserted=%d vpvaddr=0x%08x vpsgeaddr=0x%08x vpssladdr=0x%08x",
				selectedHandle,
				feav,
				list,
				antecedentVoice,
				topRegister,
				topBefore,
				inserted ? 1 : 0,
				GetRegister32(NV_PAPU_VPVADDR),
				GetRegister32(NV_PAPU_VPSGEADDR),
				GetRegister32(NV_PAPU_VPSSLADDR));
		}
		if (!inserted) {
			if (!m_LoggedVoiceListInsertFailure) {
				EmuLog(LOG_LEVEL::WARNING,
					"APU VOICE_ON failed to link handle=%u list=%u antecedent=0x%04x topRegister=0x%08x vpvaddr=0x%08x vpsgeaddr=0x%08x vpssladdr=0x%08x",
					selectedHandle,
					list,
					antecedentVoice,
					topRegister,
					GetRegister32(NV_PAPU_VPVADDR),
					GetRegister32(NV_PAPU_VPSGEADDR),
					GetRegister32(NV_PAPU_VPSSLADDR));
				m_LoggedVoiceListInsertFailure = true;
			}
		}
		if (inserted) {
			m_LoggedVoiceListInsertFailure = false;
		}

		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_PAUSED, 0);
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 1);
		SetVoiceActiveHint(selectedHandle, true);
		SetVoiceLocked(selectedHandle, false);
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_OFFSET,
			NV_PAVS_VOICE_PAR_OFFSET_CBO, 0);
		m_VPSSLData[selectedHandle].ssl_index = 0;
		m_VPSSLData[selectedHandle].ssl_seg = 0;
		m_VPSSLData[selectedHandle].persistCompleted = {};
		m_VPPlaybackState[selectedHandle] = PlaybackState{};
		uint32_t notifierBase = 0;
		if (ResolveOptionalGuestTableBase(NV_PAPU_FENADDR, m_VPNotifyContextDMA, notifierBase)) {
			WriteNotifierValue(selectedHandle, MCPX_HW_NOTIFIER_VOICE_POSITION, 0);
			WriteNotifierStatus(selectedHandle, MCPX_HW_NOTIFIER_VOICE_POSITION,
				NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS);
		}
		ClearHRTFFilterState(selectedHandle);
		InitializeVoiceEnvelopes(selectedHandle, value);
		m_LoggedEmptyVoiceTableDiagnostics = false;
		return;
	}
	case NV1BA0_PIO_VOICE_OFF: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_OFF_HANDLE;
		ClearStoppedVoiceState(voiceHandle);
		UnlinkVoiceFromLists(voiceHandle);
		// Sample the guest-visible offset register before clearing the local playback cache so
		// the completion notifier reflects the position software last programmed/observed.
		WriteNotifierValue(voiceHandle, MCPX_HW_NOTIFIER_VOICE_POSITION, GetVoicePlaybackOffset(voiceHandle));
		NotifyVoiceCompletion(voiceHandle, NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS);
		if (voiceHandle < m_VPPlaybackState.size()) {
			m_VPPlaybackState[voiceHandle] = PlaybackState{};
		}
		ClearHRTFFilterState(voiceHandle);
		m_LoggedEmptyVoiceTableDiagnostics = false;
		return;
	}
	case NV1BA0_PIO_VOICE_PAUSE: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_PAUSE_HANDLE;
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_PAUSED,
			(value & NV1BA0_PIO_VOICE_PAUSE_ACTION) != 0 ? 1u : 0u);
		return;
	}
	case NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY:
		m_VPCurrentHRTFEntry = value & NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY_HANDLE;
		return;
	case NV1BA0_PIO_VOICE_LOCK: {
		const uint32_t vh = currentVoice();
		const bool locking = (value & 1u) != 0;
		const bool wasLocked = IsVoiceLocked(vh);
		SetVoiceLocked(vh, locking);
		// When a voice is unlocked after being configured, and the
		// game never explicitly calls VOICE_ON, activate it so the
		// Stream Engine renders it.  Matches Bink/DirectSound
		// patterns that use lock/configure/unlock to finalise voices.
		if (!locking && wasLocked && vh < APU_VP_VOICE_MAX_HANDLE) {
			uint32_t state = 0;
			if (ReadVoiceMask(vh, NV_PAVS_VOICE_PAR_STATE,
				NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, state) &&
				(state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0) {
				// Activate the voice via the same path as VOICE_ON
				UnlinkVoiceFromLists(vh);
				WriteVoiceMask(vh, NV_PAVS_VOICE_PAR_STATE,
					NV_PAVS_VOICE_PAR_STATE_PAUSED, 0);
				WriteVoiceMask(vh, NV_PAVS_VOICE_PAR_STATE,
					NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 1);
				SetVoiceActiveHint(vh, true);
				SetVoiceLocked(vh, false);
				WriteVoiceMask(vh, NV_PAVS_VOICE_PAR_OFFSET,
					NV_PAVS_VOICE_PAR_OFFSET_CBO, 0);
				EmuLog(LOG_LEVEL::INFO,
					"APU diag: auto-activated voice %u on VOICE_LOCK unlock", vh);
			}
		}
		return;
	}
	case NV1BA0_PIO_VOICE_RELEASE: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_RELEASE_HANDLE;
		if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
			return;
		}
		BeginVoiceRelease(voiceHandle);
		return;
	}
	case NV1BA0_PIO_GET_VOICE_POSITION:
		WriteNotifierValue(value & NV1BA0_PIO_GET_VOICE_POSITION_HANDLE,
			MCPX_HW_NOTIFIER_VOICE_POSITION,
			GetVoicePlaybackOffset(value & NV1BA0_PIO_GET_VOICE_POSITION_HANDLE));
		WriteNotifierStatus(value & NV1BA0_PIO_GET_VOICE_POSITION_HANDLE,
			MCPX_HW_NOTIFIER_VOICE_POSITION,
			NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS);
		return;
	case NV1BA0_PIO_SET_CONTEXT_DMA_NOTIFY:
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			EmuLog(LOG_LEVEL::INFO,
				"APU SET_CONTEXT_DMA_NOTIFY value=0x%08x voice=0x%04x",
				value,
				currentVoice() & APU_VP_VOICE_MAX_HANDLE);
		}
		m_VPNotifyContextDMA = value;
		WriteRegister(APU_VP_BASE + addr, value, sizeof(uint32_t));
		return;
	case NV1BA0_PIO_SET_CURRENT_SSL_CONTEXT_DMA:
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			EmuLog(LOG_LEVEL::INFO,
				"APU SET_CURRENT_SSL_CONTEXT_DMA value=0x%08x voice=0x%04x",
				value,
				currentVoice() & APU_VP_VOICE_MAX_HANDLE);
		}
		m_VPCurrentSSLContextDMA = value;
		WriteRegister(APU_VP_BASE + addr, value, sizeof(uint32_t));
		return;
	case NV1BA0_PIO_SET_CURRENT_SSL:
		m_VPSSLBasePage = value & NV1BA0_PIO_SET_CURRENT_SSL_BASE_PAGE;
		return;
	case NV1BA0_PIO_SET_HRTF_SUBMIXES:
		m_VPHRTFSubmix[0] = static_cast<uint8_t>((value >> 0) & 0x1Fu);
		m_VPHRTFSubmix[1] = static_cast<uint8_t>((value >> 8) & 0x1Fu);
		m_VPHRTFSubmix[2] = static_cast<uint8_t>((value >> 16) & 0x1Fu);
		m_VPHRTFSubmix[3] = static_cast<uint8_t>((value >> 24) & 0x1Fu);
		return;
	case NV1BA0_PIO_SET_HRTF_HEADROOM:
		m_VPHRTFHeadroom = static_cast<uint8_t>(value & NV1BA0_PIO_SET_HRTF_HEADROOM_AMOUNT);
		return;
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
	case NV1BA0_PIO_SET_VOICE_TAR_HRTF:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CFG_HRTF_TARGET,
			NV_PAVS_VOICE_CFG_HRTF_TARGET_HANDLE,
			value & NV1BA0_PIO_SET_VOICE_TAR_HRTF_HANDLE);
		return;
	case NV1BA0_PIO_SET_VOICE_SSL_A:
		if (currentVoice() < m_VPSSLData.size()) {
			m_VPSSLData[currentVoice()].base[0] = (value & NV1BA0_PIO_SET_VOICE_SSL_A_BASE) >> Ctz32(NV1BA0_PIO_SET_VOICE_SSL_A_BASE);
			m_VPSSLData[currentVoice()].count[0] = static_cast<uint8_t>((value & NV1BA0_PIO_SET_VOICE_SSL_A_COUNT) >> Ctz32(NV1BA0_PIO_SET_VOICE_SSL_A_COUNT));
			m_VPSSLData[currentVoice()].persistCompleted[0] = false;
		}
		return;
	case NV1BA0_PIO_SET_VOICE_SSL_B:
		if (currentVoice() < m_VPSSLData.size()) {
			m_VPSSLData[currentVoice()].base[1] = (value & NV1BA0_PIO_SET_VOICE_SSL_A_BASE) >> Ctz32(NV1BA0_PIO_SET_VOICE_SSL_A_BASE);
			m_VPSSLData[currentVoice()].count[1] = static_cast<uint8_t>((value & NV1BA0_PIO_SET_VOICE_SSL_A_COUNT) >> Ctz32(NV1BA0_PIO_SET_VOICE_SSL_A_COUNT));
			m_VPSSLData[currentVoice()].persistCompleted[1] = false;
		}
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
	case NV1BA0_PIO_SET_VOICE_LFO_MOD:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_LFO_MOD, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_FCA:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_FCA, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_FCB:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_FCB, 0xFFFFFFFF, value);
		return;
	case NV1BA0_PIO_SET_VOICE_TAR_PITCH:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_TAR_PITCH_LINK,
			NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH,
			(value & NV1BA0_PIO_SET_VOICE_TAR_PITCH_STEP) >> Ctz32(NV1BA0_PIO_SET_VOICE_TAR_PITCH_STEP));
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_BASE:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CUR_PSL_START,
			NV_PAVS_VOICE_CUR_PSL_START_BA, value);
		if (currentVoice() < m_VPPlaybackState.size()) {
			m_VPPlaybackState[currentVoice()].valid = false;
		}
		ClearHRTFFilterState(currentVoice());
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_LBO:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_CUR_PSH_SAMPLE,
			NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO, value);
		if (currentVoice() < m_VPPlaybackState.size()) {
			m_VPPlaybackState[currentVoice()].valid = false;
		}
		ClearHRTFFilterState(currentVoice());
		return;
	case NV1BA0_PIO_SET_VOICE_BUF_CBO:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_PAR_OFFSET,
			NV_PAVS_VOICE_PAR_OFFSET_CBO, value);
		if (currentVoice() < m_VPPlaybackState.size()) {
			m_VPPlaybackState[currentVoice()] = PlaybackState{};
		}
		ClearHRTFFilterState(currentVoice());
		return;
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_EBO:
		WriteVoiceMask(currentVoice(), NV_PAVS_VOICE_PAR_NEXT,
			NV_PAVS_VOICE_PAR_NEXT_EBO, value);
		if (currentVoice() < m_VPPlaybackState.size()) {
			m_VPPlaybackState[currentVoice()].valid = false;
		}
		ClearHRTFFilterState(currentVoice());
		return;
	case NV1BA0_PIO_SET_HRIR:
	case NV1BA0_PIO_SET_HRIR + 0x04:
	case NV1BA0_PIO_SET_HRIR + 0x08:
	case NV1BA0_PIO_SET_HRIR + 0x0C:
	case NV1BA0_PIO_SET_HRIR + 0x10:
	case NV1BA0_PIO_SET_HRIR + 0x14:
	case NV1BA0_PIO_SET_HRIR + 0x18:
	case NV1BA0_PIO_SET_HRIR + 0x1C:
	case NV1BA0_PIO_SET_HRIR + 0x20:
	case NV1BA0_PIO_SET_HRIR + 0x24:
	case NV1BA0_PIO_SET_HRIR + 0x28:
	case NV1BA0_PIO_SET_HRIR + 0x2C:
	case NV1BA0_PIO_SET_HRIR + 0x30:
	case NV1BA0_PIO_SET_HRIR + 0x34:
	case NV1BA0_PIO_SET_HRIR + 0x38: {
		if (m_VPCurrentHRTFEntry >= APU_HRTF_ENTRY_COUNT) {
			return;
		}

		const size_t slot = (addr - NV1BA0_PIO_SET_HRIR) / sizeof(uint32_t);
		const size_t coefficientIndex = slot * 2;
		WriteHRTFCoefficient(m_VPCurrentHRTFEntry, 0, coefficientIndex,
			static_cast<int8_t>((value & NV1BA0_PIO_SET_HRIR_LEFT0) >> Ctz32(NV1BA0_PIO_SET_HRIR_LEFT0)));
		WriteHRTFCoefficient(m_VPCurrentHRTFEntry, 1, coefficientIndex,
			static_cast<int8_t>((value & NV1BA0_PIO_SET_HRIR_RIGHT0) >> Ctz32(NV1BA0_PIO_SET_HRIR_RIGHT0)));
		WriteHRTFCoefficient(m_VPCurrentHRTFEntry, 0, coefficientIndex + 1,
			static_cast<int8_t>((value & NV1BA0_PIO_SET_HRIR_LEFT1) >> Ctz32(NV1BA0_PIO_SET_HRIR_LEFT1)));
		WriteHRTFCoefficient(m_VPCurrentHRTFEntry, 1, coefficientIndex + 1,
			static_cast<int8_t>((value & NV1BA0_PIO_SET_HRIR_RIGHT1) >> Ctz32(NV1BA0_PIO_SET_HRIR_RIGHT1)));
		return;
	}
	case NV1BA0_PIO_SET_HRIR_X:
		if (m_VPCurrentHRTFEntry >= APU_HRTF_ENTRY_COUNT) {
			return;
		}

		WriteHRTFCoefficient(m_VPCurrentHRTFEntry, 0, APU_HRTF_COEFFICIENT_COUNT - 1,
			static_cast<int8_t>((value & NV1BA0_PIO_SET_HRIR_X_LEFT30) >> Ctz32(NV1BA0_PIO_SET_HRIR_X_LEFT30)));
		WriteHRTFCoefficient(m_VPCurrentHRTFEntry, 1, APU_HRTF_COEFFICIENT_COUNT - 1,
			static_cast<int8_t>((value & NV1BA0_PIO_SET_HRIR_X_RIGHT30) >> Ctz32(NV1BA0_PIO_SET_HRIR_X_RIGHT30)));
		m_VPHRTFEntries[m_VPCurrentHRTFEntry].itd =
			static_cast<int16_t>((value & NV1BA0_PIO_SET_HRIR_X_ITD) >> Ctz32(NV1BA0_PIO_SET_HRIR_X_ITD));
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
	case SE2FE_IDLE_VOICE:
		SetRegister32(NV_PAPU_FECV, value & APU_VP_VOICE_MAX_HANDLE);
		if ((GetRegister32(NV_PAPU_FETFORCE1) & NV_PAPU_FETFORCE1_SE2FE_IDLE_VOICE) != 0) {
			uint32_t fectl = GetRegister32(NV_PAPU_FECTL);
			fectl &= ~(NV_PAPU_FECTL_FEMETHMODE | NV_PAPU_FECTL_FETRAPREASON);
			fectl |= NV_PAPU_FECTL_FEMETHMODE_TRAPPED | NV_PAPU_FECTL_FETRAPREASON_REQUESTED;
			SetRegister32(NV_PAPU_FECTL, fectl);
			RefreshInterruptStatus();
		}
		return;
	default:
		if (addr >= NV1BA0_PIO_SET_VOICE_METHOD_FIRST && addr < NV1BA0_PIO_SET_VOICE_CFG_BUF_BASE &&
			((addr - NV1BA0_PIO_SET_VOICE_METHOD_FIRST) % sizeof(uint32_t)) == 0) {
			const uint32_t rawVoiceOffset = addr - NV1BA0_PIO_SET_VOICE_METHOD_FIRST;
			if (rawVoiceOffset < NV_PAVS_SIZE) {
				const bool affectsPlaybackCursor =
					rawVoiceOffset == NV_PAVS_VOICE_CUR_PSL_START ||
					rawVoiceOffset == NV_PAVS_VOICE_CUR_PSH_SAMPLE ||
					rawVoiceOffset == NV_PAVS_VOICE_PAR_OFFSET ||
					rawVoiceOffset == NV_PAVS_VOICE_PAR_NEXT;
				WriteVoiceMask(currentVoice(), rawVoiceOffset, 0xFFFFFFFF, value);
				if (currentVoice() < m_VPPlaybackState.size()) {
					if (rawVoiceOffset == NV_PAVS_VOICE_PAR_OFFSET) {
						m_VPPlaybackState[currentVoice()] = PlaybackState{};
					} else if (affectsPlaybackCursor) {
						m_VPPlaybackState[currentVoice()].valid = false;
					}
				}
				if (affectsPlaybackCursor) {
					ClearHRTFFilterState(currentVoice());
				}
				return;
			}
		}
		if (addr >= NV1BA0_PIO_SET_SUBMIX_HEADROOM &&
			addr < NV1BA0_PIO_SET_SUBMIX_HEADROOM + sizeof(uint32_t) * APU_MIXBIN_COUNT &&
			((addr - NV1BA0_PIO_SET_SUBMIX_HEADROOM) % sizeof(uint32_t)) == 0) {
			const size_t slot = (addr - NV1BA0_PIO_SET_SUBMIX_HEADROOM) / sizeof(uint32_t);
			if (slot < m_VPSubmixHeadroom.size()) {
				m_VPSubmixHeadroom[slot] = static_cast<uint8_t>(value & NV1BA0_PIO_SET_SUBMIX_HEADROOM_AMOUNT);
			}
			return;
		}
		if (addr >= NV1BA0_PIO_SET_SSL_SEGMENT_OFFSET && addr < 0x00000800) {
			uint32_t sslTableBase = 0;
			if (ResolveOptionalGuestTableBase(NV_PAPU_VPSSLADDR, m_VPCurrentSSLContextDMA, sslTableBase)) {
				WriteGuestWord(sslTableBase + m_VPSSLBasePage * 8 + (addr - NV1BA0_PIO_SET_SSL_SEGMENT_OFFSET), value);
			}
			return;
		}
		if (addr >= NV1BA0_PIO_SET_OUTBUF_BA && addr < NV1BA0_PIO_SET_OUTBUF_BA + 0x20 && ((addr - NV1BA0_PIO_SET_OUTBUF_BA) % 8) == 0) {
			const size_t slot = (addr - NV1BA0_PIO_SET_OUTBUF_BA) / 8;
			WriteRegister(APU_VP_BASE + addr, value & NV1BA0_PIO_SET_OUTBUF_BA_ADDRESS, sizeof(uint32_t));
			if (slot < m_VPOutBufferCursor.size()) {
				m_VPOutBufferCursor[slot] = 0;
			}
			if (slot < m_VPOutBufferPlaybackCursor.size()) {
				m_VPOutBufferPlaybackCursor[slot] = 0;
			}
			if (slot < m_VPOutBufferQueuedBytes.size()) {
				m_VPOutBufferQueuedBytes[slot] = 0;
			}
			return;
		}
		if (addr >= NV1BA0_PIO_SET_OUTBUF_LEN && addr < NV1BA0_PIO_SET_OUTBUF_LEN + 0x20 && ((addr - NV1BA0_PIO_SET_OUTBUF_LEN) % 8) == 0) {
			const size_t slot = (addr - NV1BA0_PIO_SET_OUTBUF_LEN) / 8;
			WriteRegister(APU_VP_BASE + addr, value & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE, sizeof(uint32_t));
			if (slot < m_VPOutBufferCursor.size()) {
				m_VPOutBufferCursor[slot] = 0;
			}
			if (slot < m_VPOutBufferPlaybackCursor.size()) {
				m_VPOutBufferPlaybackCursor[slot] = 0;
			}
			if (slot < m_VPOutBufferQueuedBytes.size()) {
				m_VPOutBufferQueuedBytes[slot] = 0;
			}
			return;
		}
		return;
	}
}

bool APUDevice::ReadGuestWord(uint32_t guestAddress, uint32_t& value) const
{
	if (!ReadGuestBytes(guestAddress, &value, sizeof(value))) {
		return false;
	}
	return true;
}

bool APUDevice::ReadGuestBytes(uint32_t guestAddress, void* dest, size_t size) const
{
	uintptr_t hostAddress = 0;
	if (dest == nullptr || !ResolveGuestMemoryPointer(guestAddress, size, hostAddress)) {
		return false;
	}

	std::memcpy(dest, reinterpret_cast<const void*>(hostAddress), size);
	return true;
}

bool APUDevice::WriteGuestWord(uint32_t guestAddress, uint32_t value)
{
	return WriteGuestBytes(guestAddress, &value, sizeof(value));
}

bool APUDevice::WriteGuestBytes(uint32_t guestAddress, const void* src, size_t size)
{
	uintptr_t hostAddress = 0;
	if (src == nullptr || !ResolveGuestMemoryPointer(guestAddress, size, hostAddress)) {
		return false;
	}

	std::memcpy(reinterpret_cast<void*>(hostAddress), src, size);
	return true;
}

bool APUDevice::ResolveOptionalGuestTableBase(uint32_t registerAddress, uint32_t fallbackGuestAddress, uint32_t& guestBase) const
{
	const uint32_t registerBase = GetRegister32(registerAddress);
	if (registerBase != 0 && ResolveGuestMemoryPointer(registerBase, 1)) {
		guestBase = registerBase;
		return true;
	}

	if (fallbackGuestAddress != 0 && ResolveGuestMemoryPointer(fallbackGuestAddress, 1)) {
		guestBase = fallbackGuestAddress;
		return true;
	}

	guestBase = 0;
	return false;
}

void APUDevice::SignalNotifierInterrupt()
{
	SetRegister32(NV_PAPU_ISTS, GetRegister32(NV_PAPU_ISTS) | NV_PAPU_ISTS_FEVINTSTS | NV_PAPU_ISTS_FENINTSTS);
	RefreshInterruptStatus();
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

bool APUDevice::ReadScatterGatherBytes(uint32_t sgeBase, uint32_t maxSge, uint32_t addr, void* dest, size_t size) const
{
	if (dest == nullptr || size == 0) {
		return false;
	}

	auto* bytes = reinterpret_cast<uint8_t*>(dest);
	uint32_t currentAddress = addr;
	size_t remaining = size;
	while (remaining > 0) {
		const uint32_t pageEntry = currentAddress / APU_SGE_PAGE_SIZE;
		if (pageEntry > maxSge) {
			return false;
		}

		uint32_t pageBase = 0;
		if (!ReadGuestWord(sgeBase + pageEntry * 8, pageBase)) {
			return false;
		}

		const uint32_t offsetInPage = currentAddress % APU_SGE_PAGE_SIZE;
		const size_t chunk = std::min<size_t>(remaining, APU_SGE_PAGE_SIZE - offsetInPage);
		if (!ReadGuestBytes(pageBase + offsetInPage, bytes, chunk)) {
			return false;
		}

		bytes += chunk;
		currentAddress += static_cast<uint32_t>(chunk);
		remaining -= chunk;
	}

	return true;
}

bool APUDevice::WriteScatterGatherBytes(uint32_t sgeBase, uint32_t maxSge, uint32_t addr, const void* src, size_t size)
{
	if (src == nullptr || size == 0) {
		return false;
	}

	const auto* bytes = reinterpret_cast<const uint8_t*>(src);
	uint32_t currentAddress = addr;
	size_t remaining = size;
	while (remaining > 0) {
		const uint32_t pageEntry = currentAddress / APU_SGE_PAGE_SIZE;
		if (pageEntry > maxSge) {
			return false;
		}

		uint32_t pageBase = 0;
		if (!ReadGuestWord(sgeBase + pageEntry * 8, pageBase)) {
			return false;
		}

		const uint32_t offsetInPage = currentAddress % APU_SGE_PAGE_SIZE;
		const size_t chunk = std::min<size_t>(remaining, APU_SGE_PAGE_SIZE - offsetInPage);
		if (!WriteGuestBytes(pageBase + offsetInPage, bytes, chunk)) {
			return false;
		}

		bytes += chunk;
		currentAddress += static_cast<uint32_t>(chunk);
		remaining -= chunk;
	}

	return true;
}

uint32_t APUDevice::ReadScratchWindowWithDMA(uint32_t sgeBaseRegister, uint32_t maxSgeRegister,
	const uint8_t* data, size_t length, uint32_t addr, unsigned size) const
{
	if (size == 0 || addr + size > length) {
		return 0;
	}

	const uint32_t sgeBase = GetRegister32(sgeBaseRegister);
	if (sgeBase != 0) {
		std::array<uint8_t, sizeof(uint32_t)> scratch{};
		if (ReadScatterGatherBytes(sgeBase, GetRegister32(maxSgeRegister), addr, scratch.data(), size)) {
			return ReadLE(scratch.data(), 0, size);
		}
	}

	return ReadMemoryWindow(data, length, addr, size);
}

void APUDevice::WriteScratchWindowWithDMA(uint32_t sgeBaseRegister, uint32_t maxSgeRegister,
	uint8_t* data, size_t length, uint32_t addr, uint32_t value, unsigned size)
{
	if (size == 0 || addr + size > length) {
		return;
	}

	const uint32_t sgeBase = GetRegister32(sgeBaseRegister);
	if (sgeBase != 0) {
		std::array<uint8_t, sizeof(uint32_t)> scratch{};
		WriteLE(scratch.data(), 0, value, size);
		if (WriteScatterGatherBytes(sgeBase, GetRegister32(maxSgeRegister), addr, scratch.data(), size)) {
			return;
		}
	}

	WriteMemoryWindow(data, length, addr, value, size);
}

uint32_t APUDevice::RefreshFEMemDataRegister(uint32_t fallbackValue)
{
	uint32_t value = fallbackValue;
	ReadGuestWord(GetRegister32(NV_PAPU_FEMEMADDR), value);
	SetRegister32(NV_PAPU_FEMEMDATA, value);
	return value;
}

bool APUDevice::ReadVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t& value) const
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return false;
	}
	if (!IsVoiceEntryOffsetWithinBounds(offset)) {
		return false;
	}

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		const size_t shadowOffset = static_cast<size_t>(voiceHandle) * NV_PAVS_SIZE + offset;
		const uint32_t current = ReadMemoryWindow(
			m_VPVoiceTableShadow.data(),
			m_VPVoiceTableShadow.size(),
			static_cast<uint32_t>(shadowOffset),
			sizeof(uint32_t));
		value = GetMaskedValue(current, mask);
		m_LoggedVoiceTableReadFailure = false;
		return true;
	}

	const uint32_t voiceBase = voiceTableBase + voiceHandle * NV_PAVS_SIZE + offset;
	uint32_t current = 0;
	if (!ReadGuestWord(voiceBase, current)) {
		if (!m_LoggedVoiceTableReadFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU ReadVoiceMask failed voiceTableBase=0x%08x voiceBase=0x%08x handle=%u offset=0x%08x mask=0x%08x",
				voiceTableBase,
				voiceBase,
				voiceHandle,
				offset,
				mask);
			m_LoggedVoiceTableReadFailure = true;
		}
		return false;
	}

	value = GetMaskedValue(current, mask);
	// A successful read means the voice table is reachable again, so allow a future
	// access regression to emit a fresh one-shot diagnostic.
	m_LoggedVoiceTableReadFailure = false;
	return true;
}

bool APUDevice::WriteVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t value)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return false;
	}
	if (!IsVoiceEntryOffsetWithinBounds(offset)) {
		return false;
	}

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		const size_t shadowOffset = static_cast<size_t>(voiceHandle) * NV_PAVS_SIZE + offset;
		const uint32_t current = ReadMemoryWindow(
			m_VPVoiceTableShadow.data(),
			m_VPVoiceTableShadow.size(),
			static_cast<uint32_t>(shadowOffset),
			sizeof(uint32_t));
		const uint32_t mergedValue = MergeMaskedValue(current, mask, value);
		WriteMemoryWindow(
			m_VPVoiceTableShadow.data(),
			m_VPVoiceTableShadow.size(),
			static_cast<uint32_t>(shadowOffset),
			mergedValue,
			sizeof(uint32_t));
		m_LoggedVoiceTableWriteFailure = false;
		return true;
	}

	const uint32_t voiceBase = voiceTableBase + voiceHandle * NV_PAVS_SIZE + offset;
	const bool success = WriteGuestWordMasked(voiceBase, mask, value);
	if (!success) {
		if (!m_LoggedVoiceTableWriteFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU WriteVoiceMask failed voiceTableBase=0x%08x voiceBase=0x%08x handle=%u offset=0x%08x mask=0x%08x value=0x%08x",
				voiceTableBase,
				voiceBase,
				voiceHandle,
				offset,
				mask,
				value);
			m_LoggedVoiceTableWriteFailure = true;
		}
	} else {
		m_LoggedVoiceTableWriteFailure = false;
	}
	return success;
}

bool APUDevice::WriteVPScatterGatherEntry(uint32_t handle, uint32_t value)
{
	const uint32_t sgeTableBase = GetRegister32(NV_PAPU_VPSGEADDR);
	if (sgeTableBase == 0) {
		if (!m_LoggedScatterGatherWriteFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU WriteVPScatterGatherEntry blocked sgeTableBase=0x00000000 handle=%u value=0x%08x",
				handle,
				value);
			m_LoggedScatterGatherWriteFailure = true;
		}
		return false;
	}

	const uint32_t sgeBase = sgeTableBase + handle * 8;
	const bool success = WriteGuestWord(sgeBase, value);
	if (!success) {
		if (!m_LoggedScatterGatherWriteFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU WriteVPScatterGatherEntry failed sgeTableBase=0x%08x sgeBase=0x%08x handle=%u value=0x%08x",
				sgeTableBase,
				sgeBase,
				handle,
				value);
			m_LoggedScatterGatherWriteFailure = true;
		}
	} else {
		m_LoggedScatterGatherWriteFailure = false;
	}
	return success;
}

void APUDevice::WriteNotifierValue(uint32_t voiceHandle, uint32_t notifier, uint32_t value)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE || notifier >= MCPX_HW_NOTIFIER_COUNT) {
		return;
	}

	uint32_t notifierBase = 0;
	if (!ResolveOptionalGuestTableBase(NV_PAPU_FENADDR, m_VPNotifyContextDMA, notifierBase)) {
		SignalNotifierInterrupt();
		return;
	}

	const uint32_t offset = MCPX_HW_NOTIFIER_ENTRY_SIZE *
		(MCPX_HW_NOTIFIER_BASE_OFFSET + voiceHandle * MCPX_HW_NOTIFIER_COUNT + notifier);
	WriteGuestWord(notifierBase + offset, value);
}

void APUDevice::WriteNotifierEnvelopeState(uint32_t voiceHandle, uint32_t notifier, uint8_t envState)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE || notifier >= MCPX_HW_NOTIFIER_COUNT) {
		return;
	}

	uint32_t notifierBase = 0;
	if (!ResolveOptionalGuestTableBase(NV_PAPU_FENADDR, m_VPNotifyContextDMA, notifierBase)) {
		return;
	}

	const uint32_t offset = MCPX_HW_NOTIFIER_ENTRY_SIZE *
		(MCPX_HW_NOTIFIER_BASE_OFFSET + voiceHandle * MCPX_HW_NOTIFIER_COUNT + notifier);
	WriteGuestBytes(notifierBase + offset + 14, &envState, sizeof(envState));
}

void APUDevice::SetVoiceNotifierEnvelopeState(uint32_t voiceHandle, uint8_t envState, bool force)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE || voiceHandle >= m_VPNotifierEnvelopeState.size()) {
		return;
	}

	if (!force && m_VPNotifierEnvelopeState[voiceHandle] == envState) {
		return;
	}

	m_VPNotifierEnvelopeState[voiceHandle] = envState;
	for (uint32_t notifier = 0; notifier < MCPX_HW_NOTIFIER_COUNT; ++notifier) {
		WriteNotifierEnvelopeState(voiceHandle, notifier, envState);
	}
}

uint8_t APUDevice::GetVoiceNotifierEnvelopeState(uint32_t voiceHandle) const
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE || voiceHandle >= m_VPNotifierEnvelopeState.size()) {
		return static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF);
	}

	return m_VPNotifierEnvelopeState[voiceHandle];
}

void APUDevice::WriteNotifierStatus(uint32_t voiceHandle, uint32_t notifier, uint8_t status)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE || notifier >= MCPX_HW_NOTIFIER_COUNT) {
		return;
	}

	uint32_t notifierBase = 0;
	if (!ResolveOptionalGuestTableBase(NV_PAPU_FENADDR, m_VPNotifyContextDMA, notifierBase)) {
		SignalNotifierInterrupt();
		return;
	}

	const uint32_t offset = MCPX_HW_NOTIFIER_ENTRY_SIZE *
		(MCPX_HW_NOTIFIER_BASE_OFFSET + voiceHandle * MCPX_HW_NOTIFIER_COUNT + notifier);
	const uint8_t envState = GetVoiceNotifierEnvelopeState(voiceHandle);
	WriteGuestBytes(notifierBase + offset + 14, &envState, sizeof(envState));
	WriteGuestBytes(notifierBase + offset + 15, &status, sizeof(status));

	SignalNotifierInterrupt();
}

void APUDevice::NotifyVoiceCompletion(uint32_t voiceHandle, uint8_t status)
{
	uint32_t voiceDataType = 0;
	if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_FMT, NV_PAVS_VOICE_CFG_FMT_DATA_TYPE, voiceDataType)) {
		return;
	}

	uint32_t notifier = MCPX_HW_NOTIFIER_SSLA_DONE;
	if (voiceDataType != 0 && voiceHandle < m_VPSSLData.size() &&
		m_VPSSLData[voiceHandle].ssl_index == 1) {
		notifier = MCPX_HW_NOTIFIER_SSLB_DONE;
	}

	WriteNotifierValue(voiceHandle, notifier, GetVoicePlaybackOffset(voiceHandle));
	WriteNotifierStatus(voiceHandle, notifier, status);
}

uint32_t APUDevice::GetVoicePlaybackOffset(uint32_t voiceHandle) const
{
	if (voiceHandle < m_VPPlaybackState.size()) {
		const auto& playbackState = m_VPPlaybackState[voiceHandle];
		if (playbackState.valid) {
			return playbackState.offset;
		}
	}

	uint32_t currentOffset = 0;
	// ReadVoiceMask already emits one-shot diagnostics on failure; keep the notifier
	// payload deterministic by falling back to zero when the guest voice table cannot
	// currently be read.
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_CBO, currentOffset);
	return currentOffset;
}

uint32_t APUDevice::GetVoiceNextHandle(uint32_t voiceHandle) const
{
	uint32_t nextHandle = APU_VP_VOICE_MAX_HANDLE;
	if (voiceHandle < APU_VP_VOICE_MAX_HANDLE) {
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
			NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle);
	}
	return nextHandle;
}

void APUDevice::SetVoiceNextHandle(uint32_t voiceHandle, uint32_t nextHandle)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return;
	}

	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
		NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE,
		std::min(nextHandle, APU_VP_VOICE_MAX_HANDLE));
}

void APUDevice::UnlinkVoiceFromList(uint32_t topRegister, uint32_t voiceHandle)
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return;
	}

	const uint32_t head = GetRegister32(topRegister);
	if (head == voiceHandle) {
		SetRegister32(topRegister, GetVoiceNextHandle(voiceHandle));
		return;
	}

	uint32_t previousHandle = head;
	for (size_t visited = 0; visited < APU_VP_VOICE_MAX_HANDLE && previousHandle < APU_VP_VOICE_MAX_HANDLE; ++visited) {
		const uint32_t nextHandle = GetVoiceNextHandle(previousHandle);
		if (nextHandle == voiceHandle) {
			SetVoiceNextHandle(previousHandle, GetVoiceNextHandle(voiceHandle));
			return;
		}
		if (nextHandle == previousHandle) {
			break;
		}
		previousHandle = nextHandle;
	}
}

void APUDevice::UnlinkVoiceFromLists(uint32_t voiceHandle)
{
	UnlinkVoiceFromList(NV_PAPU_TVL2D, voiceHandle);
	UnlinkVoiceFromList(NV_PAPU_TVL3D, voiceHandle);
	UnlinkVoiceFromList(NV_PAPU_TVLMP, voiceHandle);
	SetVoiceNextHandle(voiceHandle, APU_VP_VOICE_MAX_HANDLE);
}

void APUDevice::ClearStoppedVoiceState(uint32_t voiceHandle)
{
	SetVoiceActiveHint(voiceHandle, false);
	SetVoiceLocked(voiceHandle, false);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_PAUSED, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_LFOA_DELAYMODE, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_LFOF_DELAYMODE, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_EACUR,
		NV_PAVS_VOICE_PAR_STATE_EACUR_OFF);
	SetVoiceNotifierEnvelopeState(voiceHandle,
		static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF), true);
}

bool APUDevice::IsVoiceLocked(uint32_t voiceHandle) const
{
	if (voiceHandle >= MAX_VOICE_HANDLES) {
		return false;
	}

	const uint64_t mask = uint64_t{1} << (voiceHandle % 64);
	return (m_VPVoiceLocked[voiceHandle / 64] & mask) != 0;
}

void APUDevice::SetVoiceLocked(uint32_t voiceHandle, bool locked)
{
	if (voiceHandle >= MAX_VOICE_HANDLES) {
		return;
	}

	const uint64_t mask = uint64_t{1} << (voiceHandle % 64);
	if (locked) {
		m_VPVoiceLocked[voiceHandle / 64] |= mask;
	} else {
		m_VPVoiceLocked[voiceHandle / 64] &= ~mask;
	}
}

bool APUDevice::IsVoiceActiveHinted(uint32_t voiceHandle) const
{
	if (voiceHandle >= MAX_VOICE_HANDLES) {
		return false;
	}

	const uint64_t mask = uint64_t{1} << (voiceHandle % 64);
	return (m_VPActiveVoiceHints[voiceHandle / 64] & mask) != 0;
}

void APUDevice::SetVoiceActiveHint(uint32_t voiceHandle, bool active)
{
	if (voiceHandle >= MAX_VOICE_HANDLES) {
		return;
	}

	const uint64_t mask = uint64_t{1} << (voiceHandle % 64);
	if (active) {
		m_VPActiveVoiceHints[voiceHandle / 64] |= mask;
	} else {
		m_VPActiveVoiceHints[voiceHandle / 64] &= ~mask;
	}
}

bool APUDevice::ResolveVoiceAddress(uint32_t linearAddress, uint32_t& guestAddress) const
{
	if (ResolveGuestMemoryPointer(linearAddress, 1)) {
		guestAddress = linearAddress;
		return true;
	}

	const uint32_t sgeTableBase = GetRegister32(NV_PAPU_VPSGEADDR);
	if (sgeTableBase == 0) {
		guestAddress = linearAddress;
		return IsGuestRangeAccessible(guestAddress, 1);
	}

	const uint32_t entry = linearAddress / APU_SGE_PAGE_SIZE;
	uint32_t pageBase = 0;
	if (!ReadGuestWord(sgeTableBase + entry * 8, pageBase)) {
		return false;
	}

	guestAddress = pageBase + (linearAddress & (APU_SGE_PAGE_SIZE - 1));
	return IsGuestRangeAccessible(guestAddress, 1);
}

bool APUDevice::ReadVoiceBufferBytes(uint32_t linearAddress, void* dest, size_t size) const
{
	auto* out = static_cast<uint8_t*>(dest);
	if (out == nullptr) {
		return false;
	}

	for (size_t i = 0; i < size; ++i) {
		uint32_t guestAddress = 0;
		if (!ResolveVoiceAddress(linearAddress + static_cast<uint32_t>(i), guestAddress) ||
			!ReadGuestBytes(guestAddress, out + i, 1)) {
			return false;
		}
	}

	return true;
}

bool APUDevice::WriteGuestCircularBuffer(uint32_t guestAddress, uint32_t length, uint32_t& cursor,
	const void* src, size_t size)
{
	if (src == nullptr || length == 0) {
		return false;
	}

	const auto* bytes = static_cast<const uint8_t*>(src);
	size_t remaining = size;
	cursor %= length;
	while (remaining > 0) {
		const uint32_t chunkLength = std::min<uint32_t>(length - cursor, static_cast<uint32_t>(remaining));
		if (!WriteGuestBytes(guestAddress + cursor, bytes, chunkLength)) {
			return false;
		}

		bytes += chunkLength;
		remaining -= chunkLength;
		cursor = (cursor + chunkLength) % length;
	}

	return true;
}

bool APUDevice::ReadGuestCircularBuffer(uint32_t guestAddress, uint32_t length, uint32_t& cursor,
	void* dest, size_t size) const
{
	if (dest == nullptr || length == 0) {
		return false;
	}

	auto* bytes = static_cast<uint8_t*>(dest);
	size_t remaining = size;
	cursor %= length;
	while (remaining > 0) {
		const uint32_t chunkLength = std::min<uint32_t>(length - cursor, static_cast<uint32_t>(remaining));
		if (!ReadGuestBytes(guestAddress + cursor, bytes, chunkLength)) {
			return false;
		}

		bytes += chunkLength;
		remaining -= chunkLength;
		cursor = (cursor + chunkLength) % length;
	}

	return true;
}

bool APUDevice::HasGuestVPOutputBufferPlaybackPath() const
{
	for (size_t slot = 0; slot < APU_HRTF_SUBMIX_COUNT; ++slot) {
		const uint32_t bin = m_VPHRTFSubmix[slot];
		if (bin < APU_FIRST_NON_STEREO_BIN || bin >= APU_MIXBIN_COUNT || slot >= m_VPOutBufferPlaybackCursor.size()) {
			continue;
		}

		const uint32_t outBufferBaseRegister =
			GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_BA + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferLengthRegister =
			GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_LEN + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferBase = outBufferBaseRegister & NV1BA0_PIO_SET_OUTBUF_BA_ADDRESS;
		const uint32_t outBufferLength = outBufferLengthRegister & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE;
		if (outBufferBase == 0 || outBufferLength < sizeof(int16_t)) {
			continue;
		}

		// The VP render path already writes the mixed submixes into the programmed
		// guest output buffers, so once the guest enables those buffers we can route
		// playback through the guest-facing buffer bridge instead of relying on the
		// temporary host-side spatial handoff path.
		return true;
	}

	return false;
}

bool APUDevice::ConsumeGuestVPOutputBuffer(size_t slot, uint32_t guestAddress, uint32_t length,
	int16_t* dest, size_t frameCount, uint32_t* peak)
{
	const uint32_t usableLength = length & ~uint32_t(sizeof(int16_t) - 1);
	if (dest == nullptr || frameCount == 0 || slot >= m_VPOutBufferPlaybackCursor.size() ||
		slot >= m_VPOutBufferQueuedBytes.size() || usableLength < sizeof(int16_t)) {
		return false;
	}

	std::fill_n(dest, frameCount, 0);
	if (peak != nullptr) {
		*peak = 0;
	}

	const size_t requestedBytes = frameCount * sizeof(int16_t);
	const uint32_t availableBytes = std::min<uint32_t>(m_VPOutBufferQueuedBytes[slot], usableLength);
	if (availableBytes == 0) {
		return false;
	}

	const size_t readBytes = std::min<size_t>(requestedBytes, availableBytes);
	if (!ReadGuestCircularBuffer(guestAddress, usableLength, m_VPOutBufferPlaybackCursor[slot], dest, readBytes)) {
		if (!m_LoggedVPOutputBufferReadFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU guest VP output-buffer playback failed for slot=%zu base=0x%08x length=%u queued=%u",
				slot, guestAddress, usableLength, availableBytes);
			m_LoggedVPOutputBufferReadFailure = true;
		}
		return false;
	}

	m_VPOutBufferQueuedBytes[slot] -= static_cast<uint32_t>(readBytes);
	if (readBytes < requestedBytes) {
		if (!m_LoggedVPOutputBufferUnderrun) {
			EmuLog(LOG_LEVEL::INFO,
				"APU guest VP output-buffer underrun slot=%zu requested=%zu available=%u length=%u; zero-filling remainder",
				slot, requestedBytes, availableBytes, usableLength);
			m_LoggedVPOutputBufferUnderrun = true;
		}
	} else {
		m_LoggedVPOutputBufferUnderrun = false;
	}

	m_LoggedVPOutputBufferReadFailure = false;
	if (peak != nullptr) {
		*peak = audio_diagnostics::PeakAbsoluteSampleAmplitude(dest, frameCount);
	}
	return true;
}

bool APUDevice::MixGuestVPOutputBuffers(int16_t* output, size_t frameCount, std::array<uint32_t, 4>* slotPeak)
{
	if (output == nullptr || frameCount == 0) {
		return false;
	}

	constexpr std::array<size_t, APU_HRTF_SUBMIX_COUNT> kHRTFOutputChannelMapping{ 0, 1, 0, 1 };
	bool mixedAnySubmix = false;

	if (slotPeak != nullptr) {
		slotPeak->fill(0);
	}

	std::vector<int16_t> slotSamples(frameCount);
	for (size_t slot = 0; slot < APU_HRTF_SUBMIX_COUNT; ++slot) {
		const uint32_t bin = m_VPHRTFSubmix[slot];
		if (bin < APU_FIRST_NON_STEREO_BIN || bin >= APU_MIXBIN_COUNT || slot >= m_VPOutBufferPlaybackCursor.size()) {
			continue;
		}

		const uint32_t outBufferBaseRegister =
			GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_BA + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferLengthRegister =
			GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_LEN + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferBase = outBufferBaseRegister & NV1BA0_PIO_SET_OUTBUF_BA_ADDRESS;
		const uint32_t outBufferLength = outBufferLengthRegister & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE;
		if (outBufferBase == 0 || outBufferLength < sizeof(int16_t)) {
			continue;
		}

		uint32_t peakValue = 0;
		if (!ConsumeGuestVPOutputBuffer(slot, outBufferBase, outBufferLength, slotSamples.data(), frameCount, &peakValue)) {
			continue;
		}

		const size_t outputChannel = kHRTFOutputChannelMapping[slot];
		for (size_t frame = 0; frame < frameCount; ++frame) {
			const size_t outputIndex = frame * 2 + outputChannel;
			output[outputIndex] = ClampToInt16(static_cast<int32_t>(output[outputIndex]) + slotSamples[frame]);
		}
		mixedAnySubmix = true;
		if (slotPeak != nullptr) {
			(*slotPeak)[slot] = peakValue;
		}
	}

	return mixedAnySubmix;
}

bool APUDevice::SubmitGuestVPOutputBuffersToAC97(size_t frameCount, std::array<uint32_t, 4>* slotPeak)
{
	if (g_AC97 == nullptr || frameCount == 0) {
		return false;
	}

	bool stagedAnySubmix = false;
	if (slotPeak != nullptr) {
		slotPeak->fill(0);
	}

	std::vector<int16_t> slotSamples(frameCount);
	for (size_t slot = 0; slot < APU_HRTF_SUBMIX_COUNT; ++slot) {
		const uint32_t bin = m_VPHRTFSubmix[slot];
		if (bin < APU_FIRST_NON_STEREO_BIN || bin >= APU_MIXBIN_COUNT || slot >= m_VPOutBufferPlaybackCursor.size()) {
			continue;
		}

		const uint32_t outBufferBaseRegister =
			GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_BA + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferLengthRegister =
			GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_LEN + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferBase = outBufferBaseRegister & NV1BA0_PIO_SET_OUTBUF_BA_ADDRESS;
		const uint32_t outBufferLength = outBufferLengthRegister & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE;
		if (outBufferBase == 0 || outBufferLength < sizeof(int16_t)) {
			continue;
		}

		uint32_t peakValue = 0;
		if (!ConsumeGuestVPOutputBuffer(slot, outBufferBase, outBufferLength, slotSamples.data(), frameCount, &peakValue)) {
			continue;
		}

		g_AC97->SubmitGuestSpatialSubmixFrames(static_cast<uint32_t>(slot),
			static_cast<uint8_t>(bin),
			slotSamples.data(),
			frameCount);
		stagedAnySubmix = true;
		if (slotPeak != nullptr) {
			(*slotPeak)[slot] = peakValue;
		}
	}

	return stagedAnySubmix;
}

uint32_t APUDevice::ReadMemoryWindow(const uint8_t* data, size_t length, uint32_t addr, unsigned size) const
{
	if (data == nullptr || size == 0 || addr + size > length) {
		return 0;
	}

	return ReadLE(data, addr, size);
}

void APUDevice::WriteMemoryWindow(uint8_t* data, size_t length, uint32_t addr, uint32_t value, unsigned size)
{
	if (data == nullptr || size == 0 || addr + size > length) {
		return;
	}

	WriteLE(data, addr, value, size);
}

void APUDevice::WriteHRTFCoefficient(uint32_t entryIndex, size_t channel, size_t coefficientIndex, int8_t value)
{
	if (entryIndex >= m_VPHRTFEntries.size() || channel >= m_VPHRTFEntries[entryIndex].coeffs.size() ||
		coefficientIndex >= m_VPHRTFEntries[entryIndex].coeffs[channel].size()) {
		return;
	}

	m_VPHRTFEntries[entryIndex].coeffs[channel][coefficientIndex] = value;
}

void APUDevice::ClearHRTFFilterState(uint32_t voiceHandle)
{
	if (voiceHandle >= m_VPHRTFFilterState.size()) {
		return;
	}

	m_VPHRTFFilterState[voiceHandle] = HRTFFilterState{};
}

void APUDevice::SetHRTFFilterTarget(uint32_t voiceHandle, const HRTFEntryState& entry)
{
	if (voiceHandle >= m_VPHRTFFilterState.size()) {
		return;
	}

	auto& filter = m_VPHRTFFilterState[voiceHandle];
	filter.itd_tar = std::clamp(static_cast<float>(entry.itd) / APU_HRTF_ITD_SCALE,
		-APU_HRTF_MAX_DELAY_SAMPLES_FLOAT,
		APU_HRTF_MAX_DELAY_SAMPLES_FLOAT);

	for (size_t channel = 0; channel < filter.ch.size(); ++channel) {
		float sum = 0.0f;
		for (size_t coefficientIndex = 0; coefficientIndex < entry.coeffs[channel].size(); ++coefficientIndex) {
			const float coefficient = static_cast<float>(entry.coeffs[channel][coefficientIndex]) / 127.0f;
			filter.ch[channel].hrir_coeff_tar[coefficientIndex] = coefficient;
			sum += std::fabs(coefficient);
		}
		if (sum > APU_HRTF_NORMALIZATION_EPSILON &&
			std::fabs(sum - 1.0f) > APU_HRTF_NORMALIZATION_EPSILON) {
			for (float& coefficient : filter.ch[channel].hrir_coeff_tar) {
				coefficient /= sum;
			}
		}
	}
}

void APUDevice::ProcessHRTFSample(uint32_t voiceHandle, float& sampleLeft, float& sampleRight)
{
	if (voiceHandle >= m_VPHRTFFilterState.size()) {
		return;
	}

	auto& filter = m_VPHRTFFilterState[voiceHandle];
	for (size_t channel = 0; channel < filter.ch.size(); ++channel) {
		auto& currentCoefficients = filter.ch[channel].hrir_coeff_cur;
		const auto& targetCoefficients = filter.ch[channel].hrir_coeff_tar;
		for (size_t coefficientIndex = 0; coefficientIndex < currentCoefficients.size(); ++coefficientIndex) {
			currentCoefficients[coefficientIndex] += APU_HRTF_PARAM_SMOOTH_ALPHA *
				(targetCoefficients[coefficientIndex] - currentCoefficients[coefficientIndex]);
		}
	}
	filter.itd_cur += APU_HRTF_PARAM_SMOOTH_ALPHA * (filter.itd_tar - filter.itd_cur);

	const float inputSamples[2]{ sampleLeft, sampleRight };
	float outputSamples[2]{};
	const int bufferLength = static_cast<int>(HRTF_FILTER_BUFFER_LENGTH);
	const int bufferPosition = static_cast<int>(filter.buf_pos);
	for (size_t channel = 0; channel < filter.ch.size(); ++channel) {
		auto& channelState = filter.ch[channel];
		channelState.buf[filter.buf_pos] = inputSamples[channel];

		float delay = 0.0f;
		const float delayMagnitude = std::fabs(filter.itd_cur);
		if ((filter.itd_cur >= 0.0f && channel == 0) || (filter.itd_cur < 0.0f && channel == 1)) {
			delay = delayMagnitude;
		}
		const int delayInteger = static_cast<int>(delay);
		const float delayFraction = delay - static_cast<float>(delayInteger);

		float convolutionSum = 0.0f;
		for (size_t coefficientIndex = 0; coefficientIndex < channelState.hrir_coeff_cur.size(); ++coefficientIndex) {
			const int tapOffset = (delayInteger + static_cast<int>(coefficientIndex)) % bufferLength;
			const int index1 = (bufferPosition - tapOffset + bufferLength) % bufferLength;
			float delayedSample = channelState.buf[index1];
			if (delayFraction > 0.0f) {
				const int index2 = (index1 - 1 + bufferLength) % bufferLength;
				delayedSample = delayedSample * (1.0f - delayFraction) + channelState.buf[index2] * delayFraction;
			}
			convolutionSum += channelState.hrir_coeff_cur[coefficientIndex] * delayedSample;
		}

		outputSamples[channel] = convolutionSum;
	}

	sampleLeft = outputSamples[0];
	sampleRight = outputSamples[1];
	++filter.buf_pos;
	if (filter.buf_pos >= HRTF_FILTER_BUFFER_LENGTH) {
		filter.buf_pos = 0;
	}
}

void APUDevice::InitializeVoiceEnvelopes(uint32_t voiceHandle, uint32_t voiceOnValue)
{
	const auto initializeEnvelope = [&](uint32_t startState, uint32_t reg0, uint32_t regA,
		uint32_t releaseRegister, uint32_t releaseMask, uint32_t levelRegister,
		uint32_t levelMask, uint32_t countMask, uint32_t stateMask) {
		const auto& fieldConfig = GetEnvelopeFieldConfig(reg0, regA);
		uint32_t count = 0;
		uint32_t level = 0xFF;

		switch (startState) {
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF:
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_DELAY:
			ReadVoiceMask(voiceHandle, reg0, fieldConfig.delayTimeMask, count);
			count *= 16;
			level = 0;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK:
			level = 0;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD:
			ReadVoiceMask(voiceHandle, regA, fieldConfig.holdTimeMask, count);
			count *= 16;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY:
			ReadVoiceMask(voiceHandle, regA, fieldConfig.decayRateMask, count);
			count *= 16;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN:
			ReadVoiceMask(voiceHandle, regA, fieldConfig.sustainLevelMask, level);
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_RELEASE:
			ReadVoiceMask(voiceHandle, releaseRegister, releaseMask, count);
			count *= 16;
			level = 0;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_FORCE_RELEASE:
			level = 0;
			break;
		default:
			break;
		}

		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask, startState);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, level);
	};

	initializeEnvelope(
		(voiceOnValue & NV1BA0_PIO_VOICE_ON_ENVA) >> Ctz32(NV1BA0_PIO_VOICE_ON_ENVA),
		NV_PAVS_VOICE_CFG_ENV0, NV_PAVS_VOICE_CFG_ENVA,
		NV_PAVS_VOICE_TAR_LFO_ENV, NV_PAVS_VOICE_TAR_LFO_ENV_EA_RELEASERATE,
		NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_EALVL,
		NV_PAVS_VOICE_CUR_ECNT_EACOUNT, NV_PAVS_VOICE_PAR_STATE_EACUR);
	initializeEnvelope(
		(voiceOnValue & NV1BA0_PIO_VOICE_ON_ENVF) >> Ctz32(NV1BA0_PIO_VOICE_ON_ENVF),
		NV_PAVS_VOICE_CFG_ENV1, NV_PAVS_VOICE_CFG_ENVF,
		NV_PAVS_VOICE_CFG_MISC, NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE,
		NV_PAVS_VOICE_PAR_NEXT, NV_PAVS_VOICE_PAR_NEXT_EFLVL,
		NV_PAVS_VOICE_CUR_ECNT_EFCOUNT, NV_PAVS_VOICE_PAR_STATE_EFCUR);

	uint32_t misc = 0;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_MISC, 0xFFFFFFFF, misc);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOALVL, APU_LFO_LEVEL_CENTER);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOADR, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFLVL, APU_LFO_LEVEL_CENTER);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFDR, 0);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE,
		NV_PAVS_VOICE_PAR_STATE_LFOA_DELAYMODE,
		(misc & NV_PAVS_VOICE_CFG_MISC_LFOA_DELAYMODE) != 0 ? 1u : 0u);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE,
		NV_PAVS_VOICE_PAR_STATE_LFOF_DELAYMODE,
		(misc & NV_PAVS_VOICE_CFG_MISC_LFOF_DELAYMODE) != 0 ? 1u : 0u);
	uint32_t envelopeState = NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_EACUR, envelopeState);
	SetVoiceNotifierEnvelopeState(voiceHandle, static_cast<uint8_t>(envelopeState), true);
}

void APUDevice::BeginVoiceRelease(uint32_t voiceHandle)
{
	uint32_t releaseRate = 0;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_LFO_ENV,
		NV_PAVS_VOICE_TAR_LFO_ENV_EA_RELEASERATE, releaseRate);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT,
		NV_PAVS_VOICE_CUR_ECNT_EACOUNT, releaseRate * 16);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE,
		NV_PAVS_VOICE_PAR_STATE_EACUR, NV_PAVS_VOICE_PAR_STATE_EFCUR_RELEASE);

	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_MISC,
		NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE, releaseRate);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT,
		NV_PAVS_VOICE_CUR_ECNT_EFCOUNT, releaseRate * 16);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE,
		NV_PAVS_VOICE_PAR_STATE_EFCUR, NV_PAVS_VOICE_PAR_STATE_EFCUR_RELEASE);
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE,
		NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
	SetVoiceNotifierEnvelopeState(voiceHandle,
		static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_RELEASE), true);
}

void APUDevice::AdvancePausedVoiceState(uint32_t voiceHandle, size_t frameCount,
	BasicVoiceDiagnosticSummary* diagnostics)
{
	const bool captureVoiceDiagnostics = diagnostics != nullptr;
	uint32_t lfoEnv = 0;
	uint32_t lfoMod = 0;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_LFO_ENV, 0xFFFFFFFF, lfoEnv);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_LFO_MOD, 0xFFFFFFFF, lfoMod);
	const uint32_t lfoADelta = ExtractLFOField(lfoEnv, NV_PAVS_VOICE_TAR_LFO_ENV_LFOADLT);
	const uint32_t lfoFDelta = ExtractLFOField(lfoEnv, NV_PAVS_VOICE_TAR_LFO_ENV_LFOFDLT);
	const float lfoAmplitudeAmount = DecodeSignedLFOAmount(
		ExtractLFOField(lfoMod, NV_PAVS_VOICE_TAR_LFO_MOD_LFOAAM));

	for (size_t frame = 0; frame < frameCount; ++frame) {
		uint32_t state = 0;
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state) ||
			(state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0 ||
			(state & NV_PAVS_VOICE_PAR_STATE_PAUSED) == 0) {
			break;
		}

		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
		const float envelopeGain = StepVoiceEnvelope(
			voiceHandle,
			NV_PAVS_VOICE_CFG_ENV0, NV_PAVS_VOICE_CFG_ENVA,
			NV_PAVS_VOICE_TAR_LFO_ENV, NV_PAVS_VOICE_TAR_LFO_ENV_EA_RELEASERATE,
			NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_EALVL,
			NV_PAVS_VOICE_CUR_ECNT_EACOUNT, NV_PAVS_VOICE_PAR_STATE_EACUR);
		const float filterEnvelopeGain = StepVoiceEnvelope(
			voiceHandle,
			NV_PAVS_VOICE_CFG_ENV1, NV_PAVS_VOICE_CFG_ENVF,
			NV_PAVS_VOICE_CFG_MISC, NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE,
			NV_PAVS_VOICE_PAR_NEXT, NV_PAVS_VOICE_PAR_NEXT_EFLVL,
			NV_PAVS_VOICE_CUR_ECNT_EFCOUNT, NV_PAVS_VOICE_PAR_STATE_EFCUR);
		if (captureVoiceDiagnostics) {
			diagnostics->maxEnvelopeGain = std::max(diagnostics->maxEnvelopeGain, envelopeGain);
			diagnostics->maxFilterEnvelopeGain = std::max(diagnostics->maxFilterEnvelopeGain, filterEnvelopeGain);
		}
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state) ||
			(state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0 ||
			(state & NV_PAVS_VOICE_PAR_STATE_PAUSED) == 0) {
			break;
		}

		uint32_t lfoALevel = APU_LFO_LEVEL_CENTER;
		uint32_t lfoAReverse = 0;
		uint32_t lfoFLevel = APU_LFO_LEVEL_CENTER;
		uint32_t lfoFReverse = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOALVL, lfoALevel);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOADR, lfoAReverse);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFLVL, lfoFLevel);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFDR, lfoFReverse);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state);
		bool lfoADescending = lfoAReverse != 0;
		bool lfoFDescending = lfoFReverse != 0;
		StepVoiceLFOLevel(
			IsVoiceLFODelayActive(state, NV_PAVS_VOICE_PAR_STATE_LFOA_DELAYMODE, NV_PAVS_VOICE_PAR_STATE_EACUR)
				? 0u
				: lfoADelta,
			lfoALevel, lfoADescending);
		StepVoiceLFOLevel(
			IsVoiceLFODelayActive(state, NV_PAVS_VOICE_PAR_STATE_LFOF_DELAYMODE, NV_PAVS_VOICE_PAR_STATE_EFCUR)
				? 0u
				: lfoFDelta,
			lfoFLevel, lfoFDescending);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOALVL, lfoALevel);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOADR, lfoADescending ? 1u : 0u);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFLVL, lfoFLevel);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFDR, lfoFDescending ? 1u : 0u);
		if (captureVoiceDiagnostics) {
			const float lfoAValue = NormalizeVoiceLFOModulationLevel(lfoALevel);
			const float amplitudeLFOModulation = std::clamp(1.0f + lfoAValue * lfoAmplitudeAmount, 0.0f, 2.0f);
			diagnostics->maxAmplitudeLFOModulation = std::max(
				diagnostics->maxAmplitudeLFOModulation, amplitudeLFOModulation);
		}
	}
}

float APUDevice::StepVoiceEnvelope(uint32_t voiceHandle, uint32_t reg0, uint32_t regA,
	uint32_t rrReg, uint32_t rrMask, uint32_t levelRegister, uint32_t levelMask,
	uint32_t countMask, uint32_t stateMask)
{
	const auto& fieldConfig = GetEnvelopeFieldConfig(reg0, regA);
	uint32_t currentState = 0;
	if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask, currentState)) {
		return 1.0f;
	}

	const bool amplitudeEnvelope = countMask == NV_PAVS_VOICE_CUR_ECNT_EACOUNT;
	const auto stopVoice = [&]() {
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 0);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
		if (amplitudeEnvelope) {
			SetVoiceNotifierEnvelopeState(voiceHandle,
				static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF), true);
		}
	};

	switch (currentState) {
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF:
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, 0);
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0xFF);
		return 1.0f;
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_DELAY: {
		uint32_t count = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0);
		if (count == 0) {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK);
			if (amplitudeEnvelope) {
				SetVoiceNotifierEnvelopeState(voiceHandle,
					static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK),true);
			}
		} else {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count - 1);
		}
		return 0.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK: {
		uint32_t count = 0;
		uint32_t attackRate = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		ReadVoiceMask(voiceHandle, reg0, fieldConfig.attackRateMask, attackRate);

		const uint32_t attackSpan = attackRate * 16;
		uint32_t level = 0xFF;
		if (attackRate != 0 && attackSpan != 0) {
			level = std::min<uint32_t>(0xFF, static_cast<uint32_t>((count * 0xFFu) / attackSpan));
		}
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, level);

		if (attackRate == 0 || count >= attackSpan) {
			uint32_t holdTime = 0;
			ReadVoiceMask(voiceHandle, regA, fieldConfig.holdTimeMask, holdTime);
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD);
			if (amplitudeEnvelope) {
				SetVoiceNotifierEnvelopeState(voiceHandle,
					static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD),true);
			}
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, holdTime * 16);
			WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0xFF);
			return 1.0f;
		}

		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count + 1);
		return static_cast<float>(level) / 255.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD: {
		uint32_t count = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0xFF);
		if (count == 0) {
			uint32_t decayRate = 0;
			ReadVoiceMask(voiceHandle, regA, fieldConfig.decayRateMask, decayRate);
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY);
			if (amplitudeEnvelope) {
				SetVoiceNotifierEnvelopeState(voiceHandle,
					static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY),true);
			}
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, decayRate * 16);
		} else {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count - 1);
		}
		return 1.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY: {
		uint32_t count = 0;
		uint32_t decayRate = 0;
		uint32_t sustainLevel = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		ReadVoiceMask(voiceHandle, regA, fieldConfig.decayRateMask, decayRate);
		ReadVoiceMask(voiceHandle, regA, fieldConfig.sustainLevelMask, sustainLevel);

		if (decayRate == 0 || count == 0) {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN);
			if (amplitudeEnvelope) {
				SetVoiceNotifierEnvelopeState(voiceHandle,
					static_cast<uint8_t>(NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN),true);
			}
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, 0);
			WriteVoiceMask(voiceHandle, levelRegister, levelMask, sustainLevel);
			return static_cast<float>(sustainLevel) / 255.0f;
		}

		const uint32_t decaySpan = std::max<uint32_t>(1, decayRate * 16);
		const float progress = std::clamp(
			static_cast<float>(decaySpan - std::min<uint32_t>(count, decaySpan)) / static_cast<float>(decaySpan),
			0.0f, 1.0f);
		const uint32_t level = static_cast<uint32_t>(std::clamp(
			255.0f + (static_cast<float>(sustainLevel) - 255.0f) * progress,
			static_cast<float>(sustainLevel), 255.0f));

		WriteVoiceMask(voiceHandle, levelRegister, levelMask, level);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count - 1);
		return static_cast<float>(level) / 255.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN: {
		uint32_t sustainLevel = 0;
		ReadVoiceMask(voiceHandle, regA, fieldConfig.sustainLevelMask, sustainLevel);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, 0);
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, sustainLevel);
		return static_cast<float>(sustainLevel) / 255.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_RELEASE: {
		uint32_t count = 0;
		uint32_t level = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		ReadVoiceMask(voiceHandle, levelRegister, levelMask, level);

		if (count == 0) {
			WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0);
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_FORCE_RELEASE);
			if (amplitudeEnvelope) {
				stopVoice();
			}
			return 0.0f;
		}

		uint32_t releaseRate = 0;
		ReadVoiceMask(voiceHandle, rrReg, rrMask, releaseRate);
		if (releaseRate == 0) {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, 0);
			WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0);
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_FORCE_RELEASE);
			if (amplitudeEnvelope) {
				stopVoice();
			}
			return 0.0f;
		}

		const uint32_t nextCount = count - 1;
		const uint32_t nextLevel = count > 0 ? static_cast<uint32_t>((static_cast<uint64_t>(level) * nextCount) / count) : 0;
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, nextCount);
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, nextLevel);
		return static_cast<float>(nextLevel) / 255.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_FORCE_RELEASE:
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, 0);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, 0);
		if (amplitudeEnvelope) {
			stopVoice();
		}
		return 0.0f;
	default:
		return 1.0f;
	}
}

void APUDevice::SynchronizeAudio()
{
    std::lock_guard<std::mutex> lock(m_AudioUpdateMutex);

    const uint32_t now = GetAPUTime();
    const bool counterOff =
        ((GetRegister32(NV_PAPU_SECTL) & NV_PAPU_SECTL_XCNTMODE) >> Ctz32(NV_PAPU_SECTL_XCNTMODE)) ==
        NV_PAPU_SECTL_XCNTMODE_OFF;

    // Periodically scan the guest voice table for voices that were
    // activated via direct memory write (games can set ACTIVE_VOICE
    // without calling the VP VOICE_ON method).  When a voice is found
    // active in guest memory but not tracked in our shadow hints,
    // perform the activation side-effects that ConsumeVPMethod(VOICE_ON)
    // would normally handle.
    {
        static uint32_t lastVoiceScan;
        if (now - lastVoiceScan >= 4800) { // ~10 Hz
            lastVoiceScan = now;
            const uint32_t vpvaddr = GetRegister32(NV_PAPU_VPVADDR);
            if (vpvaddr != 0) {
                uint32_t newlyActivated = 0;
                for (uint32_t vh = 0; vh < MAX_VOICE_HANDLES; ++vh) {
                    uint32_t state = 0;
                    if (!ReadVoiceMask(vh, NV_PAVS_VOICE_PAR_STATE,
                        NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, state)) {
                        break;
                    }
                    const bool active = (state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) != 0;
                    const bool tracked = IsVoiceActiveHinted(vh);
                    if (active && !tracked) {
                        UnlinkVoiceFromLists(vh);
                        const uint32_t feav = GetRegister32(NV_PAPU_FEAV);
                        const uint32_t list = (feav & NV_PAPU_FEAV_LST) >> Ctz32(NV_PAPU_FEAV_LST);
                        WriteVoiceMask(vh, NV_PAVS_VOICE_PAR_STATE,
                            NV_PAVS_VOICE_PAR_STATE_PAUSED, 0);
                        SetVoiceActiveHint(vh, true);
                        SetVoiceLocked(vh, false);
                        WriteVoiceMask(vh, NV_PAVS_VOICE_PAR_OFFSET,
                            NV_PAVS_VOICE_PAR_OFFSET_CBO, 0);
                        ++newlyActivated;
                    }
                }
                if (newlyActivated > 0) {
                    EmuLog(LOG_LEVEL::INFO,
                        "APU diag: guest-memory voice scan activated %u voices not tracked in VP hints",
                        newlyActivated);
                }
            }
        }
    }

    if (counterOff) {
        const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
        bool hasActiveVoices = voiceTableBase != 0;
        if (!hasActiveVoices) {
            for (size_t i = 0; i < m_VPActiveVoiceHints.size(); ++i) {
                if (m_VPActiveVoiceHints[i] != 0) {
                    hasActiveVoices = true;
                    break;
                }
            }
        }
        if (!hasActiveVoices) {
            m_LastAudioUpdate = now;
            return;
        }
    }

    uint32_t remaining = now - m_LastAudioUpdate;
    if (remaining > 0) {
        static uint32_t lastRenderLog;
        if (now - lastRenderLog >= 48000) { // log at most once per second
            lastRenderLog = now;
            EmuLog(LOG_LEVEL::INFO,
                "APU diag: rendering %u frames (counter=%s VPVADDR=0x%08X activeHints=%u)",
                remaining,
                counterOff ? "off/bypass" : "on",
                GetRegister32(NV_PAPU_VPVADDR),
                m_VPActiveVoiceHints[0] != 0 || m_VPActiveVoiceHints[1] != 0 ? 1 : 0);
        }
    }
    while (remaining > 0) {
        const size_t chunk = std::min<size_t>(remaining, APU_AUDIO_CHUNK_FRAMES);
        RenderBasicAudioChunk(chunk);
        m_LastAudioUpdate += static_cast<uint32_t>(chunk);
        m_XGSCounter += static_cast<uint32_t>(chunk);
        remaining -= static_cast<uint32_t>(chunk);
    }
}

void APUDevice::RenderBasicAudioChunk(size_t frameCount)
{
	if (frameCount == 0) {
		return;
	}

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		if (!m_LoggedMissingVoiceTableDuringRender) {
			EmuLog(LOG_LEVEL::INFO,
				"APU using internal shadow voice table (NV_PAPU_VPVADDR not yet initialized by guest)");
			m_LoggedMissingVoiceTableDuringRender = true;
		}
	} else {
		m_LoggedMissingVoiceTableDuringRender = false;
	}

	m_EnableHostSpatialHandoff = !HasGuestVPOutputBufferPlaybackPath() && !IsAnyDSPEnabled();
	m_ChunkCaptured3DVoiceCount = 0;
	m_ChunkSubmittedHostSpatialVoiceCount = 0;
	if (g_AC97 != nullptr) {
		g_AC97->Begin3DVoiceFrameBatch();
	}

	std::vector<int32_t> mixBins(frameCount * APU_MIXBIN_COUNT, 0);
	const size_t visited2D = RenderBasicVoiceList(NV_PAPU_TVL2D, mixBins.data(), frameCount);
	const size_t visited3D = RenderBasicVoiceList(NV_PAPU_TVL3D, mixBins.data(), frameCount);
	const size_t visitedMP = RenderBasicVoiceList(NV_PAPU_TVLMP, mixBins.data(), frameCount);
	size_t fallbackVisited = 0;
	if ((visited2D + visited3D + visitedMP) == 0 &&
		GetRegister32(NV_PAPU_TVL2D) >= APU_VP_VOICE_MAX_HANDLE &&
		GetRegister32(NV_PAPU_TVL3D) >= APU_VP_VOICE_MAX_HANDLE &&
		GetRegister32(NV_PAPU_TVLMP) >= APU_VP_VOICE_MAX_HANDLE) {
		std::unordered_set<uint32_t> renderedVoiceHandles;
		const auto renderFallbackVoice = [&](uint32_t voiceHandle) {
			if (voiceHandle >= MAX_VOICE_HANDLES || IsVoiceLocked(voiceHandle) ||
				!renderedVoiceHandles.insert(voiceHandle).second) {
				return;
			}

			++fallbackVisited;
			RenderBasicVoice(voiceHandle, mixBins.data(), frameCount);
		};
		for (size_t wordIndex = 0; wordIndex < m_VPActiveVoiceHints.size(); ++wordIndex) {
			uint64_t pendingVoices = m_VPActiveVoiceHints[wordIndex];
			while (pendingVoices != 0) {
				const uint32_t bitIndex = CountTrailingZeros64(pendingVoices);
				const uint32_t voiceHandle = static_cast<uint32_t>(wordIndex * 64 + bitIndex);
				pendingVoices &= ~(uint64_t{1} << bitIndex);
				renderFallbackVoice(voiceHandle);
			}
		}
		renderFallbackVoice(GetRegister32(NV_PAPU_FECV) & NV1BA0_PIO_VOICE_ON_HANDLE);
		const size_t recentMethodCount = std::min(m_RecentFEMethodCount, m_RecentFEMethods.size());
		for (size_t i = 0; i < recentMethodCount; ++i) {
			const size_t recentIndex =
				GetRecentFEMethodIndex(m_RecentFEMethodNext, m_RecentFEMethods.size(), i);
			const auto& event = m_RecentFEMethods[recentIndex];
			renderFallbackVoice(event.targetVoice);
			renderFallbackVoice(event.currentVoice);
		}
	}
	const bool hasVoiceActivity = (visited2D + visited3D + visitedMP + fallbackVisited) != 0;
	if (fallbackVisited != 0) {
		if (!m_LoggedFallbackActiveVoiceRender) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU rendered %zu fallback candidate voices outside the guest TVL lists; list heads were all idle",
				fallbackVisited);
			m_LoggedFallbackActiveVoiceRender = true;
		}
	}
	// FE still expects the SE2FE idle-voice path to progress even before the guest
	// publishes NV_PAPU_VPVADDR; ConsumeVPMethod remains safe here because the
	// shadow voice table backs the FE/VP state touched by that path.
	if (!hasVoiceActivity &&
		(GetRegister32(NV_PAPU_FETFORCE1) & NV_PAPU_FETFORCE1_SE2FE_IDLE_VOICE) != 0 &&
		(GetRegister32(NV_PAPU_FECTL) & NV_PAPU_FECTL_FEMETHMODE) != NV_PAPU_FECTL_FEMETHMODE_TRAPPED) {
		ConsumeVPMethod(SE2FE_IDLE_VOICE, APU_VP_VOICE_MAX_HANDLE, sizeof(uint32_t));
	}
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		if (!hasVoiceActivity) {
			if (!m_LoggedEmptyVoiceTableDiagnostics) {
				LogVoiceTableDiagnostics();
				LogRecentVoiceStateDiagnostics();
				m_LoggedEmptyVoiceTableDiagnostics = true;
			}
		} else {
			m_LoggedEmptyVoiceTableDiagnostics = false;
		}
	}

	uint32_t preHeadroomStereoPeak[2]{};
	uint32_t preHeadroomDominantPeak = 0;
	size_t preHeadroomDominantBin = 0;
	bool shouldLogChunkDiagnostics = hasVoiceActivity;
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		preHeadroomStereoPeak[0] = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 0);
		preHeadroomStereoPeak[1] = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 1);
		for (size_t slot = 2; slot < APU_MIXBIN_COUNT; ++slot) {
			const uint32_t peak = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, slot);
			if (peak > preHeadroomDominantPeak) {
				preHeadroomDominantPeak = peak;
				preHeadroomDominantBin = slot;
			}
		}
		shouldLogChunkDiagnostics = shouldLogChunkDiagnostics ||
			preHeadroomStereoPeak[0] != 0 ||
			preHeadroomStereoPeak[1] != 0 ||
			preHeadroomDominantPeak != 0;
	}
	ApplySubmixHeadroom(mixBins.data(), frameCount);
	WriteOutputBuffers(mixBins.data(), frameCount);

	// Diagnostic: log mix-bin peaks and voice activity periodically
	{
		static uint32_t lastChunkLog;
		const uint32_t now = GetAPUTime();
		if (hasVoiceActivity && now - lastChunkLog >= 48000) {
			lastChunkLog = now;
			const uint32_t peak0 = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 0);
			const uint32_t peak1 = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 1);
			EmuLog(LOG_LEVEL::INFO,
				"APU diag: chunk frames=%zu voices(2D=%zu 3D=%zu MP=%zu fb=%zu) stereoPeak=[%u, %u] VPVADDR=0x%08X SECTL=0x%08X",
				frameCount, visited2D, visited3D, visitedMP, fallbackVisited,
				peak0, peak1,
				GetRegister32(NV_PAPU_VPVADDR),
				GetRegister32(NV_PAPU_SECTL));
		}
	}

	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		uint32_t postHeadroomStereoPeak[2]{
			PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 0),
			PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 1)
		};
		uint32_t postHeadroomDominantPeak = 0;
		size_t postHeadroomDominantBin = 0;
		for (size_t slot = 2; slot < APU_MIXBIN_COUNT; ++slot) {
			const uint32_t peak = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, slot);
			if (peak > postHeadroomDominantPeak) {
				postHeadroomDominantPeak = peak;
				postHeadroomDominantBin = slot;
			}
		}
		shouldLogChunkDiagnostics = shouldLogChunkDiagnostics ||
			postHeadroomStereoPeak[0] != 0 ||
			postHeadroomStereoPeak[1] != 0 ||
			postHeadroomDominantPeak != 0;

		std::array<uint32_t, 4> outBufferPeak{};
		for (size_t slot = 0; slot < outBufferPeak.size() && slot < APU_MIXBIN_COUNT; ++slot) {
			size_t sourceBin = slot;
			if (slot < APU_HRTF_SUBMIX_COUNT && m_VPHRTFSubmix[slot] < APU_MIXBIN_COUNT) {
				sourceBin = m_VPHRTFSubmix[slot];
			}
			outBufferPeak[slot] = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, sourceBin);
			shouldLogChunkDiagnostics = shouldLogChunkDiagnostics || outBufferPeak[slot] != 0;
		}

		if (shouldLogChunkDiagnostics) {
			EmuLog(LOG_LEVEL::INFO,
				"APU chunk diagnostics frames=%zu preStereo=[%u,%u] preDominantBin=%zu:%u postStereo=[%u,%u] postDominantBin=%zu:%u stereoHeadroom=[%u,%u] outBuffers=[%u,%u,%u,%u]",
				frameCount,
				preHeadroomStereoPeak[0],
				preHeadroomStereoPeak[1],
				preHeadroomDominantBin,
				preHeadroomDominantPeak,
				postHeadroomStereoPeak[0],
				postHeadroomStereoPeak[1],
				postHeadroomDominantBin,
				postHeadroomDominantPeak,
				static_cast<unsigned>(m_VPSubmixHeadroom[0]),
				static_cast<unsigned>(m_VPSubmixHeadroom[1]),
				outBufferPeak[0],
				outBufferPeak[1],
				outBufferPeak[2],
				outBufferPeak[3]);
		}
	}
	if (g_AC97 == nullptr) {
		if ((shouldLogChunkDiagnostics || hasVoiceActivity) && !m_LoggedAC97Missing) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU rendered %zu frames but AC97 is not connected, so no host audio can be submitted",
				frameCount);
			m_LoggedAC97Missing = true;
		}
		return;
	}
	m_LoggedAC97Missing = false;

	std::vector<int16_t> output(frameCount * 2);
	const bool dspOutputActive = ProcessDSPAudio(output.data(), mixBins.data(), frameCount);
	if (!dspOutputActive) {
		for (size_t frame = 0; frame < frameCount; ++frame) {
			// mixBins are stored slot-major: all frames for bin 0, then all frames for bin 1, etc.
			output[frame * 2] = static_cast<int16_t>(std::clamp<int64_t>(mixBins[frame], INT16_MIN, INT16_MAX));
			output[frame * 2 + 1] = static_cast<int16_t>(std::clamp<int64_t>(mixBins[frameCount + frame], INT16_MIN, INT16_MAX));
		}
	}
	const bool stereoBinsSilent = dspOutputActive
		? audio_diagnostics::PeakAbsoluteSampleAmplitude(output.data(), output.size()) == 0
		: (PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 0) == 0 &&
			PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, 1) == 0);
	const bool guestVPOutputPlaybackConfigured = !m_EnableHostSpatialHandoff && !dspOutputActive;
	bool guestVPOutputPlaybackActive = false;
	bool guestVPOutputStereoMixed = false;
	std::array<uint32_t, 4> guestVPOutputPeak{};
	const bool hostSpatialSubmitted = m_ChunkSubmittedHostSpatialVoiceCount != 0;
	bool useHRTFStereoFallback = false;
	bool useNonStereoBinFallback = false;
	bool blendFallbackIntoStereo = false;
	constexpr std::array<size_t, APU_HRTF_SUBMIX_COUNT> kHRTFFallbackChannelMapping{ 0, 1, 0, 1 };
	if (!guestVPOutputPlaybackConfigured && !hostSpatialSubmitted && !dspOutputActive) {
		blendFallbackIntoStereo = !stereoBinsSilent;
		for (size_t slot = 0; slot < APU_HRTF_SUBMIX_COUNT; ++slot) {
			const uint32_t bin = m_VPHRTFSubmix[slot];
			if (bin >= APU_FIRST_NON_STEREO_BIN && bin < APU_MIXBIN_COUNT &&
				PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, bin) != 0) {
				useHRTFStereoFallback = true;
				break;
			}
		}
		if (!useHRTFStereoFallback) {
			for (size_t bin = APU_FIRST_NON_STEREO_BIN; bin < APU_MIXBIN_COUNT; ++bin) {
				if (PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, bin) != 0) {
					useNonStereoBinFallback = true;
					break;
				}
			}
		}
	}

	if (guestVPOutputPlaybackConfigured) {
		guestVPOutputPlaybackActive = SubmitGuestVPOutputBuffersToAC97(frameCount, &guestVPOutputPeak);
		if (!guestVPOutputPlaybackActive) {
			guestVPOutputStereoMixed = MixGuestVPOutputBuffers(output.data(), frameCount, &guestVPOutputPeak);
		}
	}

	if (!guestVPOutputPlaybackActive && !guestVPOutputStereoMixed && !dspOutputActive) {
		for (size_t frame = 0; frame < frameCount; ++frame) {
			int64_t left = output[frame * 2];
			int64_t right = output[frame * 2 + 1];
			if (useHRTFStereoFallback) {
				for (size_t slot = 0; slot < APU_HRTF_SUBMIX_COUNT; ++slot) {
					const uint32_t bin = m_VPHRTFSubmix[slot];
					if (bin < APU_FIRST_NON_STEREO_BIN || bin >= APU_MIXBIN_COUNT) {
						continue;
					}

					const size_t binBase = static_cast<size_t>(bin) * frameCount;
					int32_t contribution = mixBins[binBase + frame];
					if (blendFallbackIntoStereo) {
						contribution /= 2;
					}
					// Fold the four global HRTF submix slots back to host stereo as L,R,L,R
					// until the dedicated OpenAL 3D handoff consumes them directly. If stereo
					// bins already contain signal, blend the fallback at half strength so
					// routed-only content stays audible without doubling fully duplicated dry
					// paths as aggressively.
					if (kHRTFFallbackChannelMapping[slot] == 0) {
						left += contribution;
					} else {
						right += contribution;
					}
				}
			} else if (useNonStereoBinFallback) {
				for (size_t bin = 2; bin < APU_MIXBIN_COUNT; ++bin) {
					const size_t binBase = bin * frameCount;
					int32_t contribution = mixBins[binBase + frame];
					if (blendFallbackIntoStereo) {
						contribution /= 2;
					}
					// Without the DSP/output-buffer stages, some voices only reach non-stereo
					// mixbins. Fold them back to host stereo by bin parity so their audio stays
					// audible until the full guest routing path is implemented, even when a
					// small amount of direct stereo signal is already present. When stereo bins
					// are already active, blend the fallback at half strength to reduce obvious
					// double-mixing of paths that intentionally target both stereo and effect
					// bins. This mirrors the voice-mixing convention above where even-numbered
					// routes originate from the left sample and odd-numbered routes originate
					// from the right sample.
					if ((bin & 1u) == 0) {
						left += contribution;
					} else {
						right += contribution;
					}
				}
			}

			output[frame * 2] = static_cast<int16_t>(std::clamp<int64_t>(left, INT16_MIN, INT16_MAX));
			output[frame * 2 + 1] = static_cast<int16_t>(std::clamp<int64_t>(right, INT16_MIN, INT16_MAX));
		}
	}

	uint32_t stereoPeak = 0;
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		const char* arbitrationWinner = "stereo-direct";
		if (dspOutputActive) {
			arbitrationWinner = "gp-ep-dsp";
		} else if (guestVPOutputPlaybackActive) {
			arbitrationWinner = "guest-vp-spatial";
		} else if (guestVPOutputStereoMixed) {
			arbitrationWinner = "guest-vp-stereo";
		} else if (hostSpatialSubmitted) {
			arbitrationWinner = "host-spatial";
		} else if (useHRTFStereoFallback) {
			arbitrationWinner = "stereo-hrtf-fallback";
		} else if (useNonStereoBinFallback) {
			arbitrationWinner = "stereo-nonstereo-fallback";
		}
		stereoPeak = audio_diagnostics::PeakAbsoluteSampleAmplitude(output.data(), output.size());
		if (hasVoiceActivity || stereoPeak != 0) {
			EmuLog(LOG_LEVEL::INFO,
				"APU playback arbitration frames=%zu peak=%u captured3DVoices=%zu hostSpatialVoices=%zu guestVPConfigured=%d guestVPActive=%d guestVPStereo=%d dspActive=%d hrtfFallback=%d nonStereoFallback=%d winner=%s bins=[%u,%u,%u,%u] guestVPPeaks=[%u,%u,%u,%u]",
				frameCount,
				static_cast<unsigned>(stereoPeak),
				m_ChunkCaptured3DVoiceCount,
				m_ChunkSubmittedHostSpatialVoiceCount,
				guestVPOutputPlaybackConfigured ? 1 : 0,
				guestVPOutputPlaybackActive ? 1 : 0,
				guestVPOutputStereoMixed ? 1 : 0,
				dspOutputActive ? 1 : 0,
				useHRTFStereoFallback ? 1 : 0,
				useNonStereoBinFallback ? 1 : 0,
				arbitrationWinner,
				static_cast<unsigned>(m_VPHRTFSubmix[0]),
				static_cast<unsigned>(m_VPHRTFSubmix[1]),
				static_cast<unsigned>(m_VPHRTFSubmix[2]),
				static_cast<unsigned>(m_VPHRTFSubmix[3]),
				guestVPOutputPeak[0],
				guestVPOutputPeak[1],
				guestVPOutputPeak[2],
				guestVPOutputPeak[3]);
		}
	}

	// Diagnostic: log output peak right before submitting to AC97
	{
		static uint32_t lastOutputLog;
		const uint32_t now = GetAPUTime();
		if (now - lastOutputLog >= 48000) {
			lastOutputLog = now;
			const uint32_t peak = audio_diagnostics::PeakAbsoluteSampleAmplitude(output.data(), frameCount * 2);
			EmuLog(LOG_LEVEL::INFO,
				"APU diag: SubmitPCMFrames frames=%zu outputPeak=%u dsp=%d guestVP=%d hostSpatial=%zu",
				frameCount, peak,
				dspOutputActive ? 1 : 0,
				!m_EnableHostSpatialHandoff ? 1 : 0,
				m_ChunkSubmittedHostSpatialVoiceCount);
		}
	}

	g_AC97->SubmitPCMFrames(output.data(), frameCount);
}

void APUDevice::ApplySubmixHeadroom(int32_t* mixBins, size_t frameCount)
{
	if (mixBins == nullptr || frameCount == 0) {
		return;
	}

	const size_t maxSlot = std::min(APU_MIXBIN_COUNT, m_VPSubmixHeadroom.size());
	for (size_t slot = 0; slot < maxSlot; ++slot) {
		const uint8_t headroom = m_VPSubmixHeadroom[slot];
		if (headroom == 0) {
			continue;
		}

		int32_t* slotMix = mixBins + slot * frameCount;
		for (size_t frame = 0; frame < frameCount; ++frame) {
			slotMix[frame] = ApplyHeadroomToMixSample(slotMix[frame], headroom);
		}
	}
}

void APUDevice::WriteOutputBuffers(const int32_t* mixBins, size_t frameCount)
{
	if (mixBins == nullptr || frameCount == 0) {
		return;
	}

	std::vector<int16_t> output(frameCount);
	for (size_t slot = 0; slot < m_VPOutBufferCursor.size() && slot < APU_MIXBIN_COUNT; ++slot) {
		const uint32_t outBufferBaseRegister = GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_BA + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferLengthRegister = GetRegister32(APU_VP_BASE + NV1BA0_PIO_SET_OUTBUF_LEN + static_cast<uint32_t>(slot) * 8);
		const uint32_t outBufferBase = outBufferBaseRegister & NV1BA0_PIO_SET_OUTBUF_BA_ADDRESS;
		const uint32_t outBufferLength = (outBufferLengthRegister & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE) &
			~uint32_t(sizeof(int16_t) - 1);
		if (outBufferBase == 0 || outBufferLength < sizeof(int16_t)) {
			m_VPOutBufferCursor[slot] = 0;
			if (slot < m_VPOutBufferPlaybackCursor.size()) {
				m_VPOutBufferPlaybackCursor[slot] = 0;
			}
			if (slot < m_VPOutBufferQueuedBytes.size()) {
				m_VPOutBufferQueuedBytes[slot] = 0;
			}
			continue;
		}

		size_t sourceBin = slot;
		if (slot < APU_HRTF_SUBMIX_COUNT && m_VPHRTFSubmix[slot] < APU_MIXBIN_COUNT) {
			sourceBin = m_VPHRTFSubmix[slot];
		}
		for (size_t frame = 0; frame < frameCount; ++frame) {
			output[frame] = ClampToInt16(mixBins[sourceBin * frameCount + frame]);
		}

		const uint32_t writtenBytes = static_cast<uint32_t>(output.size() * sizeof(output[0]));
		if (!WriteGuestCircularBuffer(outBufferBase, outBufferLength, m_VPOutBufferCursor[slot],
			output.data(), writtenBytes)) {
			continue;
		}

		if (slot < m_VPOutBufferQueuedBytes.size() && slot < m_VPOutBufferPlaybackCursor.size()) {
			const uint64_t queuedAfterWrite = static_cast<uint64_t>(m_VPOutBufferQueuedBytes[slot]) + writtenBytes;
			if (queuedAfterWrite > outBufferLength) {
				const uint32_t droppedBytes = static_cast<uint32_t>(queuedAfterWrite - outBufferLength);
				m_VPOutBufferPlaybackCursor[slot] =
					(m_VPOutBufferPlaybackCursor[slot] + droppedBytes) % outBufferLength;
				m_VPOutBufferQueuedBytes[slot] = outBufferLength;
				if (!m_LoggedVPOutputBufferOverrun) {
					EmuLog(LOG_LEVEL::INFO,
						"APU guest VP output-buffer overrun slot=%zu dropped=%u length=%u; keeping newest audio",
						slot, droppedBytes, outBufferLength);
					m_LoggedVPOutputBufferOverrun = true;
				}
			} else {
				m_VPOutBufferQueuedBytes[slot] = static_cast<uint32_t>(queuedAfterWrite);
				m_LoggedVPOutputBufferOverrun = false;
			}
		}
	}
}

size_t APUDevice::RenderBasicVoiceList(uint32_t topRegister, int32_t* mixBins, size_t frameCount)
{
	uint32_t voiceHandle = GetRegister32(topRegister);
	const uint32_t listHead = voiceHandle;
	size_t visitedVoiceCount = 0;
	size_t activeVoiceCount = 0;
	size_t mixedVoiceCount = 0;
	size_t decodedNonZeroCount = 0;
	size_t stereoContributionCount = 0;
	size_t nonStereoOnlyCount = 0;
	std::vector<BasicVoiceDiagnosticSummary> interestingVoices;
	const auto requiresDiagnosticAttention = [](const BasicVoiceDiagnosticSummary& diagnostics) {
		// Capture the main failure modes we are tracing:
		// 1) decode produced signal but none of it reached stereo bins,
		// 2) the active voice stayed silent all the way through decode,
		// 3) frames rendered without any offset advance despite a positive pitch step.
		return diagnostics.active && !diagnostics.paused &&
			((diagnostics.decodedNonZero && !diagnostics.stereoContribution) ||
			 !diagnostics.decodedNonZero ||
			 (diagnostics.framesRendered != 0 && diagnostics.offsetAdvance == 0 && diagnostics.pitchStep > 0.0));
	};
	for (size_t visited = 0; visited < 16384 && voiceHandle < APU_VP_VOICE_MAX_HANDLE; ++visited) {
		uint32_t nextHandle = APU_VP_VOICE_MAX_HANDLE;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
			NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle);
		BasicVoiceDiagnosticSummary diagnostics;
		uint32_t state = 0;
		const bool hasState = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state);
		const bool active = hasState && (state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) != 0;
		if (!active) {
			SetVoiceActiveHint(voiceHandle, false);
			ConsumeVPMethod(SE2FE_IDLE_VOICE, voiceHandle, sizeof(uint32_t));
			UnlinkVoiceFromLists(voiceHandle);
		} else if (!IsVoiceLocked(voiceHandle)) {
			RenderBasicVoice(voiceHandle, mixBins, frameCount, &diagnostics);
		}
		if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
			++visitedVoiceCount;
			if (diagnostics.active) {
				++activeVoiceCount;
			}
			if (diagnostics.decodedNonZero) {
				++decodedNonZeroCount;
			}
			if (diagnostics.mixed) {
				++mixedVoiceCount;
			}
			if (diagnostics.stereoContribution) {
				++stereoContributionCount;
			}
			if (diagnostics.nonStereoContribution && !diagnostics.stereoContribution) {
				++nonStereoOnlyCount;
			}
			if (interestingVoices.size() < APU_DIAGNOSTIC_MAX_VOICES_TO_LOG &&
				requiresDiagnosticAttention(diagnostics)) {
				interestingVoices.push_back(diagnostics);
			}
		}
		if (nextHandle == voiceHandle) {
			break;
		}
		voiceHandle = nextHandle;
	}
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		if (visitedVoiceCount != 0 || !interestingVoices.empty()) {
			EmuLog(LOG_LEVEL::INFO,
				"APU voice list diagnostics top=0x%08x head=0x%08x visited=%zu active=%zu decodedNonZero=%zu mixed=%zu stereoContrib=%zu nonStereoOnly=%zu frames=%zu",
				topRegister,
				listHead,
				visitedVoiceCount,
				activeVoiceCount,
				decodedNonZeroCount,
				mixedVoiceCount,
				stereoContributionCount,
				nonStereoOnlyCount,
				frameCount);
			for (const auto& diagnostics : interestingVoices) {
				EmuLog(LOG_LEVEL::INFO,
					"APU voice diag handle=%u paused=%d offsets=%u+%u frames=%zu pitchStep=%.6f/%.6f decodedPeak=%u mixedPeak=%u envMax=%.3f filterEnvMax=%.3f lfoGainMax=%.3f stereo=%d other=%d bins=[%u,%u,%u,%u,%u,%u,%u,%u] volumes=[%u,%u,%u,%u,%u,%u,%u,%u] headroom=[%u,%u,%u,%u,%u,%u,%u,%u]",
					diagnostics.voiceHandle,
					diagnostics.paused ? 1 : 0,
					diagnostics.startOffset,
					diagnostics.offsetAdvance,
					diagnostics.framesRendered,
					diagnostics.pitchStep,
					diagnostics.maxPitchStep,
					diagnostics.decodedPeak,
					diagnostics.mixedPeak,
					diagnostics.maxEnvelopeGain,
					diagnostics.maxFilterEnvelopeGain,
					diagnostics.maxAmplitudeLFOModulation,
					diagnostics.stereoContribution ? 1 : 0,
					diagnostics.nonStereoContribution ? 1 : 0,
					diagnostics.bins[0], diagnostics.bins[1], diagnostics.bins[2], diagnostics.bins[3],
					diagnostics.bins[4], diagnostics.bins[5], diagnostics.bins[6], diagnostics.bins[7],
					diagnostics.volumes[0], diagnostics.volumes[1], diagnostics.volumes[2], diagnostics.volumes[3],
					diagnostics.volumes[4], diagnostics.volumes[5], diagnostics.volumes[6], diagnostics.volumes[7],
					static_cast<unsigned>(diagnostics.headroom[0]), static_cast<unsigned>(diagnostics.headroom[1]),
					static_cast<unsigned>(diagnostics.headroom[2]), static_cast<unsigned>(diagnostics.headroom[3]),
					static_cast<unsigned>(diagnostics.headroom[4]), static_cast<unsigned>(diagnostics.headroom[5]),
					static_cast<unsigned>(diagnostics.headroom[6]), static_cast<unsigned>(diagnostics.headroom[7]));
			}
		}
	}
	return visitedVoiceCount;
}

void APUDevice::RecordRecentFEMethod(uint32_t addr, uint32_t value, uint32_t currentVoiceValue)
{
	if constexpr (!audio_diagnostics::kEnableDiagnosticLogging) {
		return;
	}

	RecentFEMethodDiagnostic event{};
	event.sequence = ++m_RecentFEMethodSequence;
	event.addr = addr;
	event.value = value;
	event.currentVoice = currentVoiceValue & APU_VP_VOICE_MAX_HANDLE;
	event.targetVoice = GetFEMethodTargetVoiceOrDefault(addr, value, currentVoiceValue);
	event.feav = GetRegister32(NV_PAPU_FEAV);
	event.vpvaddr = GetRegister32(NV_PAPU_VPVADDR);
	event.vpsgeaddr = GetRegister32(NV_PAPU_VPSGEADDR);
	event.vpssladdr = GetRegister32(NV_PAPU_VPSSLADDR);

	m_RecentFEMethods[m_RecentFEMethodNext] = event;
	m_RecentFEMethodNext = (m_RecentFEMethodNext + 1) % m_RecentFEMethods.size();
	if (m_RecentFEMethodCount < m_RecentFEMethods.size()) {
		++m_RecentFEMethodCount;
	}
}

void APUDevice::LogRecentFEMethodDiagnostics() const
{
	if constexpr (!audio_diagnostics::kEnableDiagnosticLogging) {
		return;
	}

	if (m_RecentFEMethodCount == 0) {
		EmuLog(LOG_LEVEL::INFO, "APU recent FE method diagnostics [none]");
		return;
	}

	size_t startIndex = 0;
	if (m_RecentFEMethodCount == m_RecentFEMethods.size()) {
		startIndex = m_RecentFEMethodNext;
	}
	for (size_t i = 0; i < m_RecentFEMethodCount; ++i) {
		const auto& event = m_RecentFEMethods[(startIndex + i) % m_RecentFEMethods.size()];
		const char* name = GetAPUVPMethodTraceName(event.addr);
		const uint32_t list = (event.feav & NV_PAPU_FEAV_LST) >> Ctz32(NV_PAPU_FEAV_LST);
		const uint32_t antecedentVoice = event.feav & NV_PAPU_FEAV_VALUE;
		EmuLog(LOG_LEVEL::INFO,
			"APU recent FE method[%zu/%zu] seq=%u name=%s addr=0x%08x value=0x%08x currentVoice=0x%04x targetVoice=0x%04x list=%u antecedent=0x%04x vpvaddr=0x%08x vpsgeaddr=0x%08x vpssladdr=0x%08x",
			i + 1,
			m_RecentFEMethodCount,
			event.sequence,
			name != nullptr ? name : "UNKNOWN",
			event.addr,
			event.value,
			event.currentVoice,
			event.targetVoice,
			list,
			antecedentVoice,
			event.vpvaddr,
			event.vpsgeaddr,
			event.vpssladdr);
	}
}

void APUDevice::LogVoiceTableDiagnostics() const
{
	if constexpr (!audio_diagnostics::kEnableDiagnosticLogging) {
		return;
	}

	const auto logFEVPRegisterDiagnostics = [this]() {
		const uint32_t feav = GetRegister32(NV_PAPU_FEAV);
		const uint32_t list = (feav & NV_PAPU_FEAV_LST) >> Ctz32(NV_PAPU_FEAV_LST);
		const uint32_t antecedentVoice = feav & NV_PAPU_FEAV_VALUE;
		const uint32_t decodedMethod = GetRegister32(NV_PAPU_FEDECMETH);
		const char* decodedMethodName = GetAPUVPMethodTraceName(decodedMethod);
		EmuLog(LOG_LEVEL::INFO,
			"APU FE/VP register diagnostics currentVoice=0x%04x targetVoice=0x%04x list=%u antecedent=0x%04x lastMethod=%s addr=0x%08x value=0x%08x vpFree=0x%08x vpvaddr=0x%08x vpsgeaddr=0x%08x vpssladdr=0x%08x top=[0x%08x,0x%08x,0x%08x]",
			GetRegister32(NV_PAPU_FECV) & APU_VP_VOICE_MAX_HANDLE,
			decodedMethod == NV1BA0_PIO_VOICE_ON
				? (GetRegister32(NV_PAPU_FEDECPARAM) & NV1BA0_PIO_VOICE_ON_HANDLE)
				: (GetRegister32(NV_PAPU_FECV) & APU_VP_VOICE_MAX_HANDLE),
			list,
			antecedentVoice,
			decodedMethodName != nullptr ? decodedMethodName : "UNKNOWN",
			decodedMethod,
			GetRegister32(NV_PAPU_FEDECPARAM),
			GetRegister32(APU_VP_BASE + APU_VP_FREE),
			GetRegister32(NV_PAPU_VPVADDR),
			GetRegister32(NV_PAPU_VPSGEADDR),
			GetRegister32(NV_PAPU_VPSSLADDR),
			GetRegister32(NV_PAPU_TVL2D),
			GetRegister32(NV_PAPU_TVL3D),
			GetRegister32(NV_PAPU_TVLMP));
	};

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	const char* voiceTableSource = voiceTableBase == 0 ? "shadow" : "guest";

	size_t activeVoiceCount = 0;
	size_t pausedVoiceCount = 0;
	size_t newVoiceCount = 0;
	std::array<uint32_t, APU_DIAGNOSTIC_MAX_ACTIVE_VOICES_TO_LOG> activeHandles;
	activeHandles.fill(APU_VP_VOICE_MAX_HANDLE);
	size_t loggedActiveHandles = 0;
	for (uint32_t voiceHandle = 0; voiceHandle < APU_VP_VOICE_MAX_HANDLE; ++voiceHandle) {
		uint32_t state = 0;
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state) ||
			(state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0) {
			continue;
		}
		++activeVoiceCount;
		if ((state & NV_PAVS_VOICE_PAR_STATE_PAUSED) != 0) {
			++pausedVoiceCount;
		}
		if ((state & NV_PAVS_VOICE_PAR_STATE_NEW_VOICE) != 0) {
			++newVoiceCount;
		}
		if (loggedActiveHandles < activeHandles.size()) {
			activeHandles[loggedActiveHandles++] = voiceHandle;
		}
	}

	if (loggedActiveHandles == 0) {
		EmuLog(LOG_LEVEL::INFO,
			"APU voice table diagnostics source=%s voiceTableBase=0x%08x active=%zu paused=%zu new=%zu handles=[none]",
			voiceTableSource,
			voiceTableBase,
			activeVoiceCount,
			pausedVoiceCount,
			newVoiceCount);
		logFEVPRegisterDiagnostics();
		LogRecentFEMethodDiagnostics();
		return;
	}

	EmuLog(LOG_LEVEL::INFO,
		"APU voice table diagnostics source=%s voiceTableBase=0x%08x active=%zu paused=%zu new=%zu handles=[%u,%u,%u,%u]",
		voiceTableSource,
		voiceTableBase,
		activeVoiceCount,
		pausedVoiceCount,
		newVoiceCount,
		activeHandles[0],
		loggedActiveHandles > 1 ? activeHandles[1] : APU_VP_VOICE_MAX_HANDLE,
		loggedActiveHandles > 2 ? activeHandles[2] : APU_VP_VOICE_MAX_HANDLE,
		loggedActiveHandles > 3 ? activeHandles[3] : APU_VP_VOICE_MAX_HANDLE);
	logFEVPRegisterDiagnostics();
}

void APUDevice::LogRecentVoiceStateDiagnostics() const
{
	if constexpr (!audio_diagnostics::kEnableDiagnosticLogging) {
		return;
	}

	struct VoiceDiagnosticCandidate {
		uint32_t handle = APU_VP_VOICE_MAX_HANDLE;
		const char* reason = nullptr;
	};

	std::array<VoiceDiagnosticCandidate, 4> candidates{};
	size_t candidateCount = 0;
	const auto addCandidate = [&](uint32_t handle, const char* reason) {
		if (handle >= APU_VP_VOICE_MAX_HANDLE) {
			return;
		}
		for (size_t i = 0; i < candidateCount; ++i) {
			if (candidates[i].handle == handle) {
				return;
			}
		}
		if (candidateCount < candidates.size()) {
			candidates[candidateCount++] = VoiceDiagnosticCandidate{ handle, reason };
		}
	};

	addCandidate(GetRegister32(NV_PAPU_FECV) & NV1BA0_PIO_VOICE_ON_HANDLE, "fecv");
	const size_t recentMethodCount = std::min(m_RecentFEMethodCount, m_RecentFEMethods.size());
	for (size_t i = 0; i < recentMethodCount && candidateCount < candidates.size(); ++i) {
		const size_t recentIndex =
			GetRecentFEMethodIndex(m_RecentFEMethodNext, m_RecentFEMethods.size(), i);
		const auto& event = m_RecentFEMethods[recentIndex];
		addCandidate(event.currentVoice, "recent-current");
		addCandidate(event.targetVoice, "recent-target");
	}

	if (candidateCount == 0) {
		EmuLog(LOG_LEVEL::INFO, "APU recent voice state diagnostics [none]");
		return;
	}

	for (size_t candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex) {
		const uint32_t voiceHandle = candidates[candidateIndex].handle;
		uint32_t state = 0;
		uint32_t format = 0;
		uint32_t baseAddress = 0;
		uint32_t currentOffset = 0;
		uint32_t endOffset = 0;
		uint32_t loopOffset = 0;
		uint32_t bin0 = 0;
		uint32_t bin1 = 0;
		uint32_t volume0 = 0x0FFF;
		uint32_t volume1 = 0x0FFF;
		const bool hasState = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state);
		const bool hasFormat = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_FMT, 0xFFFFFFFF, format);
		const bool hasBase = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_PSL_START,
			NV_PAVS_VOICE_CUR_PSL_START_BA, baseAddress);
		const bool hasCurrent = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET,
			NV_PAVS_VOICE_PAR_OFFSET_CBO, currentOffset);
		const bool hasEnd = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_NEXT,
			NV_PAVS_VOICE_PAR_NEXT_EBO, endOffset);
		const bool hasLoop = ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_PSH_SAMPLE,
			NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO, loopOffset);
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V0BIN, bin0)) {
			bin0 = 0;
		}
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V1BIN, bin1)) {
			bin1 = 1;
		}
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLA, NV_PAVS_VOICE_TAR_VOLA_VOLUME0, volume0);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLA, NV_PAVS_VOICE_TAR_VOLA_VOLUME1, volume1);

		const bool streaming = hasFormat && (format & NV_PAVS_VOICE_CFG_FMT_DATA_TYPE) != 0;
		const bool stereo = hasFormat && (format & NV_PAVS_VOICE_CFG_FMT_STEREO) != 0;
		const uint32_t channels = stereo ? 2u : 1u;
		const uint32_t sampleSize = hasFormat
			? ((format & NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE))
			: 0u;
		const uint32_t containerSizeMode = hasFormat
			? ((format & NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE))
			: 0u;
		const uint32_t samplesPerBlock = hasFormat
			? (((format & NV_PAVS_VOICE_CFG_FMT_SAMPLES_PER_BLOCK) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_SAMPLES_PER_BLOCK)) + 1u)
			: 0u;
		const bool adpcm = containerSizeMode == NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_ADPCM;

		uint32_t containerSize = 0;
		switch (containerSizeMode) {
		case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B8:
			containerSize = 1;
			break;
		case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B16:
			containerSize = 2;
			break;
		case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_ADPCM:
			containerSize = XBOX_ADPCM_SRCSIZE;
			break;
		case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B32:
			containerSize = 4;
			break;
		default:
			break;
		}

		uint32_t previewAddress = 0;
		bool previewAddressValid = hasBase && hasCurrent && containerSize != 0;
		if (previewAddressValid) {
			if (adpcm && samplesPerBlock != 0) {
				const uint32_t bytesPerBlock = containerSize * channels;
				previewAddress = baseAddress + (currentOffset / samplesPerBlock) * bytesPerBlock;
			} else if (streaming) {
				previewAddress = baseAddress + currentOffset * containerSize * channels;
			} else {
				previewAddress = baseAddress + currentOffset;
			}
		}

		std::array<uint8_t, 8> previewBytes{};
		bool previewRead = false;
		if (previewAddressValid) {
			previewRead = ReadVoiceBufferBytes(previewAddress, previewBytes.data(), previewBytes.size());
		}

		int32_t previewSampleLeft = 0;
		int32_t previewSampleRight = 0;
		bool decodedPreview = false;
		if (previewAddressValid && adpcm && samplesPerBlock != 0) {
			const uint32_t bytesPerBlock = containerSize * channels;
			std::array<uint8_t, APU_XADPCM_MAX_SOURCE_BLOCK_BYTES> encodedBlock{};
			std::array<int16_t, APU_XADPCM_MAX_DECODED_SAMPLES> decodedSamples{};
			if (ReadVoiceBufferBytes(previewAddress, encodedBlock.data(), bytesPerBlock)) {
				const int decodedBytes = TXboxAdpcmDecoder_Decode_Memory(
					encodedBlock.data(),
					static_cast<int>(bytesPerBlock),
					reinterpret_cast<uint8_t*>(decodedSamples.data()),
					static_cast<int>(channels));
				if (decodedBytes > 0) {
					const uint32_t sampleIndex = (currentOffset % samplesPerBlock) * channels;
					if (sampleIndex < decodedSamples.size() &&
						(channels == 1 || (sampleIndex + 1) < decodedSamples.size())) {
						previewSampleLeft = decodedSamples[sampleIndex];
						previewSampleRight = channels > 1 ? decodedSamples[sampleIndex + 1] : decodedSamples[sampleIndex];
						decodedPreview = true;
					}
				}
			}
		} else if (previewRead) {
			switch (sampleSize) {
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_U8:
				previewSampleLeft = static_cast<int32_t>(previewBytes[0]) - 128;
				previewSampleRight = channels > 1 ? static_cast<int32_t>(previewBytes[1]) - 128 : previewSampleLeft;
				decodedPreview = true;
				break;
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S16:
				previewSampleLeft = static_cast<int32_t>(static_cast<int16_t>(
					static_cast<uint16_t>(previewBytes[0]) | (static_cast<uint16_t>(previewBytes[1]) << 8)));
				if (channels > 1) {
					previewSampleRight = static_cast<int32_t>(static_cast<int16_t>(
						static_cast<uint16_t>(previewBytes[2]) | (static_cast<uint16_t>(previewBytes[3]) << 8)));
				} else {
					previewSampleRight = previewSampleLeft;
				}
				decodedPreview = true;
				break;
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S24: {
				const uint32_t packedLeft = static_cast<uint32_t>(previewBytes[0]) |
					(static_cast<uint32_t>(previewBytes[1]) << 8) |
					(static_cast<uint32_t>(previewBytes[2]) << 16);
				previewSampleLeft = (static_cast<int32_t>(packedLeft << 8)) >> 8;
				if (channels > 1) {
					const uint32_t packedRight = static_cast<uint32_t>(previewBytes[3]) |
						(static_cast<uint32_t>(previewBytes[4]) << 8) |
						(static_cast<uint32_t>(previewBytes[5]) << 16);
					previewSampleRight = (static_cast<int32_t>(packedRight << 8)) >> 8;
				} else {
					previewSampleRight = previewSampleLeft;
				}
				decodedPreview = true;
				break;
			}
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S32:
				previewSampleLeft = static_cast<int32_t>(
					static_cast<uint32_t>(previewBytes[0]) |
					(static_cast<uint32_t>(previewBytes[1]) << 8) |
					(static_cast<uint32_t>(previewBytes[2]) << 16) |
					(static_cast<uint32_t>(previewBytes[3]) << 24));
				if (channels > 1) {
					previewSampleRight = static_cast<int32_t>(
						static_cast<uint32_t>(previewBytes[4]) |
						(static_cast<uint32_t>(previewBytes[5]) << 8) |
						(static_cast<uint32_t>(previewBytes[6]) << 16) |
						(static_cast<uint32_t>(previewBytes[7]) << 24));
				} else {
					previewSampleRight = previewSampleLeft;
				}
				decodedPreview = true;
				break;
			default:
				break;
			}
		}

		EmuLog(LOG_LEVEL::INFO,
			"APU recent voice state diag reason=%s handle=%u active=%d paused=%d new=%d fmt=0x%08x base=0x%08x cbo=0x%08x ebo=0x%08x lbo=0x%08x bins=[%u,%u] volumes=[%u,%u] previewAddr=0x%08x previewRead=%d previewBytes=[%02x,%02x,%02x,%02x,%02x,%02x,%02x,%02x] previewSample=[%d,%d]",
			candidates[candidateIndex].reason != nullptr ? candidates[candidateIndex].reason : "unknown",
			voiceHandle,
			hasState && (state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) != 0 ? 1 : 0,
			hasState && (state & NV_PAVS_VOICE_PAR_STATE_PAUSED) != 0 ? 1 : 0,
			hasState && (state & NV_PAVS_VOICE_PAR_STATE_NEW_VOICE) != 0 ? 1 : 0,
			hasFormat ? format : 0,
			hasBase ? baseAddress : 0,
			hasCurrent ? currentOffset : 0,
			hasEnd ? endOffset : 0,
			hasLoop ? loopOffset : 0,
			bin0,
			bin1,
			volume0,
			volume1,
			previewAddressValid ? previewAddress : 0,
			previewRead ? 1 : 0,
			static_cast<unsigned>(previewBytes[0]), static_cast<unsigned>(previewBytes[1]),
			static_cast<unsigned>(previewBytes[2]), static_cast<unsigned>(previewBytes[3]),
			static_cast<unsigned>(previewBytes[4]), static_cast<unsigned>(previewBytes[5]),
			static_cast<unsigned>(previewBytes[6]), static_cast<unsigned>(previewBytes[7]),
			decodedPreview ? previewSampleLeft : 0,
			decodedPreview ? previewSampleRight : 0);
	}
}

void APUDevice::RenderBasicVoice(uint32_t voiceHandle, int32_t* mixBins, size_t frameCount,
	BasicVoiceDiagnosticSummary* diagnostics)
{
	const bool captureVoiceDiagnostics = diagnostics != nullptr;
	if (captureVoiceDiagnostics) {
		*diagnostics = {};
		diagnostics->voiceHandle = voiceHandle;
		diagnostics->visited = true;
	}

	uint32_t state = 0;
	if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, state) ||
		(state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0) {
		if ((state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0) {
			SetVoiceActiveHint(voiceHandle, false);
		}
		return;
	}
	if ((state & NV_PAVS_VOICE_PAR_STATE_PAUSED) != 0) {
		if (captureVoiceDiagnostics) {
			diagnostics->active = true;
			diagnostics->paused = true;
		}
		AdvancePausedVoiceState(voiceHandle, frameCount, diagnostics);
		return;
	}

	uint32_t format = 0;
	if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_FMT, 0xFFFFFFFF, format)) {
		return;
	}
	if (captureVoiceDiagnostics) {
		diagnostics->active = true;
	}

	const bool streaming = (format & NV_PAVS_VOICE_CFG_FMT_DATA_TYPE) != 0;
	const bool multipass = (format & NV_PAVS_VOICE_CFG_FMT_MULTIPASS) != 0;
	const bool persist = (format & NV_PAVS_VOICE_CFG_FMT_PERSIST) != 0;
	const bool clearMix = (format & NV_PAVS_VOICE_CFG_FMT_CLEAR_MIX) != 0;
	const uint32_t multipassBin = (format & NV_PAVS_VOICE_CFG_FMT_MULTIPASS_BIN) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_MULTIPASS_BIN);
	const uint32_t samplesPerBlock = ((format & NV_PAVS_VOICE_CFG_FMT_SAMPLES_PER_BLOCK) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_SAMPLES_PER_BLOCK)) + 1;

	const bool loop = (format & NV_PAVS_VOICE_CFG_FMT_LOOP) != 0;
	const bool stereo = (format & NV_PAVS_VOICE_CFG_FMT_STEREO) != 0;
	const uint32_t channels = stereo ? 2u : 1u;
	const uint32_t sampleSize = (format & NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE);
	const uint32_t containerSizeMode = (format & NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE) >> Ctz32(NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE);
	const bool adpcm = containerSizeMode == NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_ADPCM;
	if (multipass) {
		if (multipassBin >= APU_MIXBIN_COUNT) {
			return;
		}
	} else if (adpcm) {
		if (samplesPerBlock != APU_XADPCM_PCM_SAMPLES_PER_BLOCK) {
			return;
		}
	}

	uint32_t containerSize = 0;
	switch (containerSizeMode) {
	case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B8:
		containerSize = 1;
		break;
	case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B16:
		containerSize = 2;
		break;
	case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_ADPCM:
		containerSize = XBOX_ADPCM_SRCSIZE;
		break;
	case NV_PAVS_VOICE_CFG_FMT_CONTAINER_SIZE_B32:
		containerSize = 4;
		break;
	default:
		return;
	}

	uint32_t bins[8]{};
	bins[0] = 0;
	bins[1] = 1;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V0BIN, bins[0]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V1BIN, bins[1]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V2BIN, bins[2]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V3BIN, bins[3]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V4BIN, bins[4]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_VBIN, NV_PAVS_VOICE_CFG_VBIN_V5BIN, bins[5]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_FMT, NV_PAVS_VOICE_CFG_FMT_V6BIN, bins[6]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_FMT, NV_PAVS_VOICE_CFG_FMT_V7BIN, bins[7]);
	if (voiceHandle < APU_MAX_3D_VOICES) {
		// MCPX 3D voices override bins 0-3 with the global HRTF submix destinations.
		for (size_t binIndex = 0; binIndex < APU_HRTF_SUBMIX_COUNT; ++binIndex) {
			bins[binIndex] = m_VPHRTFSubmix[binIndex];
		}
	}
	if (bins[0] == 0 && bins[1] == 0) {
		bins[1] = 1;
	}

	uint32_t volumes[8]{};
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLA, NV_PAVS_VOICE_TAR_VOLA_VOLUME0, volumes[0]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLA, NV_PAVS_VOICE_TAR_VOLA_VOLUME1, volumes[1]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLB, NV_PAVS_VOICE_TAR_VOLB_VOLUME2, volumes[2]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLB, NV_PAVS_VOICE_TAR_VOLB_VOLUME3, volumes[3]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLC, NV_PAVS_VOICE_TAR_VOLC_VOLUME4, volumes[4]);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLC, NV_PAVS_VOICE_TAR_VOLC_VOLUME5, volumes[5]);
	uint32_t volume6 = 0;
	uint32_t volume7 = 0;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLA, NV_PAVS_VOICE_TAR_VOLA_VOLUME6_B3_0, volume6);
	volumes[6] |= volume6;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLB, NV_PAVS_VOICE_TAR_VOLB_VOLUME6_B7_4, volume6);
	volumes[6] |= volume6 << 4;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLC, NV_PAVS_VOICE_TAR_VOLC_VOLUME6_B11_8, volume6);
	volumes[6] |= volume6 << 8;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLA, NV_PAVS_VOICE_TAR_VOLA_VOLUME7_B3_0, volume7);
	volumes[7] |= volume7;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLB, NV_PAVS_VOICE_TAR_VOLB_VOLUME7_B7_4, volume7);
	volumes[7] |= volume7 << 4;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_VOLC, NV_PAVS_VOICE_TAR_VOLC_VOLUME7_B11_8, volume7);
	volumes[7] |= volume7 << 8;
	std::array<uint32_t, APU_HRTF_SUBMIX_COUNT> hrtfSubmixVolumes{};
	// MCPX 3D voices route HRTF output through volumes[0..3], which feed the four global HRTF submix slots.
	std::copy_n(volumes, APU_HRTF_SUBMIX_COUNT, hrtfSubmixVolumes.begin());
	if (captureVoiceDiagnostics) {
		std::copy_n(bins, 8, diagnostics->bins);
		std::copy_n(volumes, 8, diagnostics->volumes);
	}

	uint32_t baseAddress = 0;
	uint32_t currentOffset = 0;
	uint32_t endOffset = 0;
	uint32_t loopOffset = 0;
	uint32_t pitch = 0;
	uint32_t filterMode = 0;
	uint32_t lfoEnv = 0;
	uint32_t lfoMod = 0;
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_PSL_START, NV_PAVS_VOICE_CUR_PSL_START_BA, baseAddress);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_CBO, currentOffset);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_NEXT, NV_PAVS_VOICE_PAR_NEXT_EBO, endOffset);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_PSH_SAMPLE, NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO, loopOffset);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_PITCH_LINK, NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH, pitch);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_LFO_ENV, 0xFFFFFFFF, lfoEnv);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_LFO_MOD, 0xFFFFFFFF, lfoMod);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_MISC, NV_PAVS_VOICE_CFG_MISC_FMODE, filterMode);

	if (!multipass && !streaming && !loop && currentOffset > endOffset) {
		return;
	}
	if (!multipass && !streaming && loop && loopOffset > endOffset) {
		return;
	}

	const double pitchStep = DecodePitchStep(pitch);
	if (captureVoiceDiagnostics) {
		diagnostics->startOffset = currentOffset;
		diagnostics->pitchStep = pitchStep;
	}
	const uint32_t lfoADelta = ExtractLFOField(lfoEnv, NV_PAVS_VOICE_TAR_LFO_ENV_LFOADLT);
	const uint32_t lfoFDelta = ExtractLFOField(lfoEnv, NV_PAVS_VOICE_TAR_LFO_ENV_LFOFDLT);
	const float lfoAmplitudeAmount = DecodeSignedLFOAmount(
		ExtractLFOField(lfoMod, NV_PAVS_VOICE_TAR_LFO_MOD_LFOAAM));
	const float lfoAmplitudePitchAmount = DecodeSignedLFOAmount(
		ExtractLFOField(lfoMod, NV_PAVS_VOICE_TAR_LFO_MOD_LFOAFM));
	const float lfoAmplitudeCutoffAmount = DecodeSignedLFOAmount(
		ExtractLFOField(lfoMod, NV_PAVS_VOICE_TAR_LFO_MOD_LFOAFC));
	const float lfoPitchAmount = DecodeSignedLFOAmount(
		ExtractLFOField(lfoMod, NV_PAVS_VOICE_TAR_LFO_MOD_LFOFFM));
	const uint32_t bytesPerFrame = adpcm ? 0u : containerSize * channels;
	// The guest programs CBO/EBO/LBO in byte units for non-streaming PCM voices.
	// Streaming SSL segments and ADPCM blocks still advance in decoded-sample units.
	const uint32_t offsetStep = (!multipass && !streaming && !adpcm) ? bytesPerFrame : 1u;
	const uint32_t bytesPerBlock = containerSize * channels;
	auto& playbackState = m_VPPlaybackState[voiceHandle];
	if (!playbackState.valid || playbackState.offset != currentOffset) {
		playbackState.offset = currentOffset;
		playbackState.fraction = 0.0;
		playbackState.previewDecodeFailures = 0;
		playbackState.valid = true;
		m_VPLowPassState[voiceHandle] = {};
		ClearHRTFFilterState(voiceHandle);
	}
	uint32_t cachedADPCMBlockIndex = 0;
	uint32_t cachedADPCMBaseAddress = 0;
	uint8_t cachedADPCMChannels = 0;
	bool cachedADPCMValid = false;
	std::array<int16_t, APU_XADPCM_MAX_DECODED_SAMPLES> cachedADPCMSamples{};

	auto sslData = m_VPSSLData[voiceHandle];
	bool shouldSkipStateWrites = false;
	auto terminateVoiceWithStatus = [&](uint8_t completionStatus) {
		shouldSkipStateWrites = true;
		playbackState.previewDecodeFailures = 0;
		ClearStoppedVoiceState(voiceHandle);
		WriteNotifierValue(voiceHandle, MCPX_HW_NOTIFIER_VOICE_POSITION, currentOffset);
		NotifyVoiceCompletion(voiceHandle, completionStatus);
		UnlinkVoiceFromLists(voiceHandle);
		playbackState = PlaybackState{};
		m_VPLowPassState[voiceHandle] = {};
		ClearHRTFFilterState(voiceHandle);
	};
	auto stopVoice = [&]() {
		terminateVoiceWithStatus(NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS);
	};
	auto failVoice = [&]() {
		terminateVoiceWithStatus(NV1BA0_NOTIFICATION_STATUS_DONE_ERROR);
	};
	std::vector<int32_t> multipassSource;
	if (multipass && clearMix) {
		multipassSource.assign(mixBins + multipassBin * frameCount, mixBins + (multipassBin + 1) * frameCount);
		std::fill_n(mixBins + multipassBin * frameCount, frameCount, 0);
	}
	auto recordDecodedPeak = [&](float sampleLeft, float sampleRight) {
		if (!captureVoiceDiagnostics) {
			return;
		}
		const float peakSample = std::max(std::fabs(sampleLeft), std::fabs(sampleRight));
		if (peakSample <= 0.0f) {
			return;
		}
		diagnostics->decodedNonZero = true;
		diagnostics->decodedPeak = std::max(diagnostics->decodedPeak,
			static_cast<uint32_t>(peakSample * APU_SAMPLE_SCALE_FACTOR));
	};
	auto mixSamples = [&](float sampleLeft, float sampleRight, float envelopeGain, size_t frame) {
		const float channelSamples[2]{ sampleLeft, sampleRight };
		for (size_t binIndex = 0; binIndex < 8; ++binIndex) {
			if (bins[binIndex] >= APU_MIXBIN_COUNT) {
				continue;
			}
			uint8_t headroom = 0;
			if (voiceHandle < APU_MAX_3D_VOICES && binIndex < APU_HRTF_SUBMIX_COUNT) {
				headroom = m_VPHRTFHeadroom;
			} else if (bins[binIndex] < m_VPSubmixHeadroom.size()) {
				headroom = m_VPSubmixHeadroom[bins[binIndex]];
			}
			if (captureVoiceDiagnostics) {
				diagnostics->headroom[binIndex] = headroom;
			}
			// Guest submix/HRTF headroom is applied after voices have accumulated into
			// their destination mixbins (or on the AC97/OpenAL handoff path for host
			// spatial playback). Attenuating here as well squares the headroom amount.
			const float gain = AttenuateVoiceVolume(volumes[binIndex]) * envelopeGain;
			if (gain == 0.0f) {
				continue;
			}
			const float sample = channelSamples[binIndex % channels];
			const int32_t contribution = static_cast<int32_t>(sample * gain * APU_SAMPLE_SCALE_FACTOR);
			mixBins[bins[binIndex] * frameCount + frame] += contribution;
			if (captureVoiceDiagnostics && contribution != 0) {
				const uint32_t magnitude = AbsoluteMixMagnitude(contribution);
				diagnostics->mixed = true;
				diagnostics->mixedPeak = std::max(diagnostics->mixedPeak, magnitude);
				if (bins[binIndex] <= 1) {
					diagnostics->stereoContribution = true;
				} else {
					diagnostics->nonStereoContribution = true;
				}
			}
		}
	};
	const bool lowPassEnabled = voiceHandle < APU_MAX_3D_VOICES
		? filterMode == 1
		: (stereo ? filterMode == 1 : (filterMode & 1u) != 0);
	float lowPassCutoff[2]{};
	float lowPassResonance[2]{};
	if (lowPassEnabled) {
		for (uint32_t channel = 0; channel < channels; ++channel) {
			const uint32_t registerOffset = channel == 0 ? NV_PAVS_VOICE_TAR_FCA : NV_PAVS_VOICE_TAR_FCB;
			uint32_t cutoff = 0;
			uint32_t resonance = 0;
			ReadVoiceMask(voiceHandle, registerOffset, NV_PAVS_VOICE_TAR_FCA_FC0, cutoff);
			ReadVoiceMask(voiceHandle, registerOffset, NV_PAVS_VOICE_TAR_FCA_FC1, resonance);
			lowPassCutoff[channel] = std::clamp(std::pow(2.0f, static_cast<float>(static_cast<int16_t>(cutoff)) / 4096.0f),
				APU_FILTER_MIN_FREQUENCY, 1.0f);
			lowPassResonance[channel] = std::clamp(static_cast<float>(resonance) / APU_FILTER_Q_NORMALIZER,
				APU_FILTER_MIN_Q, 1.0f);
		}
		if (channels == 1) {
			lowPassCutoff[1] = lowPassCutoff[0];
			lowPassResonance[1] = lowPassResonance[0];
		}
	}
	auto applyLowPass = [&](float& sampleLeft, float& sampleRight, float envGain, float cutoffLFOOctaves) {
		if (!lowPassEnabled) {
			return;
		}
		const float clampedFilterEnvelopeGain = std::clamp(
			envGain, APU_FILTER_ENV_MIN_GAIN, APU_FILTER_ENV_MAX_GAIN);
		float modulatedCutoff[2]{};
		for (size_t channel = 0; channel < std::size(modulatedCutoff); ++channel) {
			modulatedCutoff[channel] = APU_FILTER_MIN_FREQUENCY +
				(lowPassCutoff[channel] - APU_FILTER_MIN_FREQUENCY) * clampedFilterEnvelopeGain;
			modulatedCutoff[channel] = std::clamp(
				static_cast<float>(modulatedCutoff[channel] * std::exp2(cutoffLFOOctaves)),
				APU_FILTER_MIN_FREQUENCY,
				1.0f);
		}
		auto& filterState = m_VPLowPassState[voiceHandle];
		sampleLeft = ClampUnitSample(RunLowPassFilter(filterState[0].high, filterState[0].band, filterState[0].low,
			modulatedCutoff[0], lowPassResonance[0], sampleLeft));
		sampleRight = ClampUnitSample(RunLowPassFilter(filterState[1].high, filterState[1].band, filterState[1].low,
			modulatedCutoff[1], lowPassResonance[1], sampleRight));
	};
	uint32_t hrtfEntryIndex = APU_INVALID_HRTF_ENTRY_INDEX;
	const bool hrtfEnabled = voiceHandle < APU_MAX_3D_VOICES &&
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_HRTF_TARGET,
			NV_PAVS_VOICE_CFG_HRTF_TARGET_HANDLE, hrtfEntryIndex) &&
		hrtfEntryIndex < m_VPHRTFEntries.size();
	const float guestHRTFPan = [&]() {
		if (!hrtfEnabled) {
			return 0.0f;
		}

		const auto& hrtfEntry = m_VPHRTFEntries[hrtfEntryIndex];
		float leftMagnitude = 0.0f;
		float rightMagnitude = 0.0f;
		for (const int8_t coefficient : hrtfEntry.coeffs[0]) {
			leftMagnitude += std::fabs(static_cast<float>(coefficient));
		}
		for (const int8_t coefficient : hrtfEntry.coeffs[1]) {
			rightMagnitude += std::fabs(static_cast<float>(coefficient));
		}
		const float normalizedItd = std::clamp(
			static_cast<float>(hrtfEntry.itd) / APU_HRTF_ITD_NORMALIZER,
			-1.0f,
			1.0f);
		const float magnitudeSum = leftMagnitude + rightMagnitude;
		const float magnitudeBalance = magnitudeSum > APU_HRTF_NORMALIZATION_EPSILON
			? ((rightMagnitude - leftMagnitude) / magnitudeSum)
			: 0.0f;
		return std::clamp(
			normalizedItd * APU_HRTF_PAN_ITD_WEIGHT +
			magnitudeBalance * APU_HRTF_PAN_MAGNITUDE_WEIGHT,
			-1.0f,
			1.0f);
	}();
	const bool capture3DHandoff = hrtfEnabled && g_AC97 != nullptr;
	const bool submit3DToAC97 = capture3DHandoff && m_EnableHostSpatialHandoff;
	if (capture3DHandoff) {
		m_VP3DVoiceCaptureScratch.assign(frameCount * 2, 0);
		++m_ChunkCaptured3DVoiceCount;
	}
	auto storeCaptured3DSample = [&](size_t frame, float sampleLeft, float sampleRight, float envelopeGain) {
		const size_t sampleIndex = frame * 2;
		if (sampleIndex >= m_VP3DVoiceCaptureScratch.size() || (sampleIndex + 1) >= m_VP3DVoiceCaptureScratch.size()) {
			return;
		}
		m_VP3DVoiceCaptureScratch[sampleIndex] = ConvertFloatSampleToInt16(sampleLeft * envelopeGain);
		m_VP3DVoiceCaptureScratch[sampleIndex + 1] = ConvertFloatSampleToInt16(sampleRight * envelopeGain);
	};
	if (hrtfEnabled) {
		SetHRTFFilterTarget(voiceHandle, m_VPHRTFEntries[hrtfEntryIndex]);
	}
	auto holdPersistentStream = [&](SSLData& voiceSSLData, uint32_t exhaustedIndex) {
		currentOffset = 0;
		voiceSSLData.ssl_seg = 0;
		playbackState = PlaybackState{};
		m_VPLowPassState[voiceHandle] = {};
		ClearHRTFFilterState(voiceHandle);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_CBO, 0);
		if (exhaustedIndex < voiceSSLData.persistCompleted.size() &&
			!voiceSSLData.persistCompleted[exhaustedIndex]) {
			WriteNotifierStatus(voiceHandle,
				exhaustedIndex == 0 ? MCPX_HW_NOTIFIER_SSLA_DONE : MCPX_HW_NOTIFIER_SSLB_DONE,
				NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS);
			voiceSSLData.persistCompleted[exhaustedIndex] = true;
		}
	};
	auto readSampleBytes = [&](uint32_t sampleAddress, void* dest, size_t size) {
		// Streaming voices can surface VP-linear offsets, raw physical offsets, or
		// KSEG0/physical-map addresses depending on how the title programmed the SSL
		// tables, so route them through the same voice-buffer translation path as
		// non-streaming voices.
		return ReadVoiceBufferBytes(sampleAddress, dest, size);
	};
	auto loadStreamingSegment = [&](SSLData& voiceSSLData, uint32_t& segmentBaseAddress, uint32_t& segmentEndOffset, uint32_t& segmentCurrentOffset, bool commit) -> bool {
		for (size_t attempts = 0; attempts < 4; ++attempts) {
			if (voiceSSLData.ssl_index > 1) {
				voiceSSLData.ssl_index = 0;
			}
			const uint32_t sslIndex = voiceSSLData.ssl_index;
			if (voiceSSLData.count[sslIndex] == 0) {
				if (commit) {
					if (persist) {
						holdPersistentStream(voiceSSLData, sslIndex);
					} else {
						voiceSSLData.ssl_index = 0;
						stopVoice();
					}
				}
				return false;
			}
			if (voiceSSLData.ssl_seg >= voiceSSLData.count[sslIndex]) {
				if (commit) {
					WriteNotifierStatus(voiceHandle,
						sslIndex == 0 ? MCPX_HW_NOTIFIER_SSLA_DONE : MCPX_HW_NOTIFIER_SSLB_DONE,
						NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS);
				}
				voiceSSLData.ssl_index = 1 - voiceSSLData.ssl_index;
				voiceSSLData.ssl_seg = 0;
				segmentCurrentOffset = 0;
				continue;
			}

			uint32_t sslTableBase = 0;
			if (!ResolveOptionalGuestTableBase(NV_PAPU_VPSSLADDR, m_VPCurrentSSLContextDMA, sslTableBase)) {
				if (!m_LoggedStreamingSSLFailure) {
					EmuLog(LOG_LEVEL::WARNING,
						"APU streaming voice %u needs an SSL table, but neither NV_PAPU_VPSSLADDR nor the current SSL context DMA resolved to guest memory",
						voiceHandle);
					m_LoggedStreamingSSLFailure = true;
				}
				if (commit) {
					failVoice();
				}
				return false;
			}
			m_LoggedStreamingSSLFailure = false;

			const uint32_t segmentPage = voiceSSLData.base[sslIndex] + static_cast<uint32_t>(voiceSSLData.ssl_seg);
			uint32_t segmentOffset = 0;
			uint32_t segmentLength = 0;
			if (!ReadGuestWord(sslTableBase + segmentPage * 8, segmentOffset) ||
				!ReadGuestWord(sslTableBase + segmentPage * 8 + 4, segmentLength)) {
				if (commit) {
					failVoice();
				}
				return false;
			}

			const uint32_t segmentSamples = segmentLength & 0xFFFF;
			const uint32_t segmentContainerSizeMode = (segmentLength >> 16) & 0x3;
			const uint32_t segmentSamplesPerBlock = ((segmentLength >> 18) & 0x1F) + 1;
			const bool segmentStereo = ((segmentLength >> 23) & 1u) != 0;
			if (segmentSamples == 0) {
				++voiceSSLData.ssl_seg;
				segmentCurrentOffset = 0;
				continue;
			}
			if (segmentContainerSizeMode != containerSizeMode ||
				segmentSamplesPerBlock != samplesPerBlock ||
				segmentStereo != stereo) {
				if (commit) {
					failVoice();
				}
				return false;
			}

			segmentBaseAddress = segmentOffset;
			segmentEndOffset = segmentSamples - 1;
			if (sslIndex < voiceSSLData.persistCompleted.size()) {
				voiceSSLData.persistCompleted[sslIndex] = false;
			}
			return true;
		}

		if (commit) {
			stopVoice();
		}
		return false;
	};
	auto advancePlaybackPosition = [&](SSLData& voiceSSLData, uint32_t& segmentBaseAddress, uint32_t& segmentEndOffset, uint32_t& segmentCurrentOffset, bool commit) -> bool {
		while (segmentCurrentOffset > segmentEndOffset) {
			if (streaming) {
				++voiceSSLData.ssl_seg;
				segmentCurrentOffset = 0;
				if (!loadStreamingSegment(voiceSSLData, segmentBaseAddress, segmentEndOffset, segmentCurrentOffset, commit)) {
					return false;
				}
			} else if (loop) {
				segmentCurrentOffset = loopOffset;
			} else {
				if (commit) {
					stopVoice();
				}
				return false;
			}
		}
		return true;
	};
	auto decodeFrame = [&](uint32_t segmentBaseAddress, uint32_t segmentCurrentOffset, float& sampleLeft, float& sampleRight) -> bool {
		sampleLeft = 0.0f;
		sampleRight = 0.0f;
		if (multipass) {
			const int32_t mixedSample = multipassSource.empty()
				? mixBins[multipassBin * frameCount + segmentCurrentOffset]
				: multipassSource[segmentCurrentOffset];
			sampleLeft = static_cast<float>(mixedSample) / APU_SAMPLE_SCALE_FACTOR;
			sampleRight = sampleLeft;
			recordDecodedPeak(sampleLeft, sampleRight);
			return true;
		}
		if (adpcm) {
			const uint32_t blockIndex = segmentCurrentOffset / samplesPerBlock;
			const uint32_t sampleIndexInBlock = segmentCurrentOffset % samplesPerBlock;
			std::array<uint8_t, APU_XADPCM_MAX_SOURCE_BLOCK_BYTES> encodedBlock{};
			if (!cachedADPCMValid ||
				cachedADPCMBlockIndex != blockIndex ||
				cachedADPCMBaseAddress != segmentBaseAddress ||
				cachedADPCMChannels != channels) {
				const uint32_t blockAddress = segmentBaseAddress + blockIndex * bytesPerBlock;
				if (!readSampleBytes(blockAddress, encodedBlock.data(), bytesPerBlock)) {
					return false;
				}

				const uint32_t expectedDecodedBytes = samplesPerBlock * channels * sizeof(int16_t);
				const int decodedBytes = TXboxAdpcmDecoder_Decode_Memory(
					encodedBlock.data(),
					static_cast<int>(bytesPerBlock),
					reinterpret_cast<uint8_t*>(cachedADPCMSamples.data()),
					static_cast<int>(channels));
				if (decodedBytes != static_cast<int>(expectedDecodedBytes)) {
					if (!m_LoggedXADPCMDecodeFailure) {
						EmuLog(LOG_LEVEL::WARNING,
							"APU XADPCM decode failed for voice %u block %u (decoded %d bytes, expected %u)",
							voiceHandle, blockIndex, decodedBytes, expectedDecodedBytes);
						m_LoggedXADPCMDecodeFailure = true;
					}
					return false;
				}

				cachedADPCMBlockIndex = blockIndex;
				cachedADPCMBaseAddress = segmentBaseAddress;
				cachedADPCMChannels = static_cast<uint8_t>(channels);
				cachedADPCMValid = true;
			}

			const uint32_t sampleIndex = sampleIndexInBlock * channels;
			const int16_t leftSample = cachedADPCMSamples[sampleIndex];
			sampleLeft = ConvertSigned16(leftSample);
			if (channels > 1) {
				const int16_t rightSample = cachedADPCMSamples[sampleIndex + 1];
				sampleRight = ConvertSigned16(rightSample);
			} else {
				sampleRight = sampleLeft;
			}
			recordDecodedPeak(sampleLeft, sampleRight);
			return true;
		}

		const uint32_t linearAddress = segmentBaseAddress +
			(streaming ? (segmentCurrentOffset * bytesPerFrame) : segmentCurrentOffset);
		for (uint32_t channel = 0; channel < channels; ++channel) {
			const uint32_t sampleAddress = linearAddress + channel * containerSize;
			float sample = 0.0f;
			switch (sampleSize) {
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_U8: {
				uint8_t raw = 0;
				if (!readSampleBytes(sampleAddress, &raw, sizeof(raw))) {
					return false;
				}
				sample = ConvertUnsigned8(raw);
				break;
			}
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S16: {
				int16_t raw = 0;
				if (!readSampleBytes(sampleAddress, &raw, sizeof(raw))) {
					return false;
				}
				sample = ConvertSigned16(raw);
				break;
			}
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S24: {
				uint8_t raw[4]{};
				if (!readSampleBytes(sampleAddress, raw, 3)) {
					return false;
				}
				const uint32_t packed = static_cast<uint32_t>(raw[0]) |
					(static_cast<uint32_t>(raw[1]) << 8) |
					(static_cast<uint32_t>(raw[2]) << 16);
				sample = ConvertSigned24(packed);
				break;
			}
			case NV_PAVS_VOICE_CFG_FMT_SAMPLE_SIZE_S32: {
				int32_t raw = 0;
				if (!readSampleBytes(sampleAddress, &raw, sizeof(raw))) {
					return false;
				}
				sample = ConvertSigned32(raw);
				break;
			}
			default:
				return false;
			}

			if (channel == 0) {
				sampleLeft = sample;
				sampleRight = sample;
			} else {
				sampleRight = sample;
			}
		}
		recordDecodedPeak(sampleLeft, sampleRight);
		return true;
	};
	if (!multipass && streaming && !loadStreamingSegment(sslData, baseAddress, endOffset, currentOffset, true)) {
		if (!shouldSkipStateWrites) {
			m_VPSSLData[voiceHandle] = sslData;
			playbackState.offset = currentOffset;
		}
		return;
	}

	for (size_t frame = 0; frame < frameCount; ++frame) {
		if (!multipass && !advancePlaybackPosition(sslData, baseAddress, endOffset, currentOffset, true)) {
			break;
		}

		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
		const float envelopeGain = StepVoiceEnvelope(
			voiceHandle,
			NV_PAVS_VOICE_CFG_ENV0, NV_PAVS_VOICE_CFG_ENVA,
			NV_PAVS_VOICE_TAR_LFO_ENV, NV_PAVS_VOICE_TAR_LFO_ENV_EA_RELEASERATE,
			NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_EALVL,
			NV_PAVS_VOICE_CUR_ECNT_EACOUNT, NV_PAVS_VOICE_PAR_STATE_EACUR);
		if (captureVoiceDiagnostics) {
			diagnostics->maxEnvelopeGain = std::max(diagnostics->maxEnvelopeGain, envelopeGain);
		}
		const float filterEnvelopeGain = StepVoiceEnvelope(
			voiceHandle,
			NV_PAVS_VOICE_CFG_ENV1, NV_PAVS_VOICE_CFG_ENVF,
			NV_PAVS_VOICE_CFG_MISC, NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE,
			NV_PAVS_VOICE_PAR_NEXT, NV_PAVS_VOICE_PAR_NEXT_EFLVL,
			NV_PAVS_VOICE_CUR_ECNT_EFCOUNT, NV_PAVS_VOICE_PAR_STATE_EFCUR);
		if (captureVoiceDiagnostics) {
			diagnostics->maxFilterEnvelopeGain = std::max(diagnostics->maxFilterEnvelopeGain, filterEnvelopeGain);
		}
		uint32_t lfoALevel = APU_LFO_LEVEL_CENTER;
		uint32_t lfoAReverse = 0;
		uint32_t lfoFLevel = APU_LFO_LEVEL_CENTER;
		uint32_t lfoFReverse = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOALVL, lfoALevel);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOADR, lfoAReverse);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFLVL, lfoFLevel);
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFDR, lfoFReverse);
		uint32_t lfoState = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, lfoState);
		bool lfoADescending = lfoAReverse != 0;
		bool lfoFDescending = lfoFReverse != 0;
		StepVoiceLFOLevel(
			IsVoiceLFODelayActive(lfoState, NV_PAVS_VOICE_PAR_STATE_LFOA_DELAYMODE, NV_PAVS_VOICE_PAR_STATE_EACUR)
				? 0u
				: lfoADelta,
			lfoALevel, lfoADescending);
		StepVoiceLFOLevel(
			IsVoiceLFODelayActive(lfoState, NV_PAVS_VOICE_PAR_STATE_LFOF_DELAYMODE, NV_PAVS_VOICE_PAR_STATE_EFCUR)
				? 0u
				: lfoFDelta,
			lfoFLevel, lfoFDescending);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOALVL, lfoALevel);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOADR, lfoADescending ? 1u : 0u);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFLVL, lfoFLevel);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_LFO, NV_PAVS_VOICE_PAR_LFO_LFOFDR, lfoFDescending ? 1u : 0u);
		const float lfoAValue = NormalizeVoiceLFOModulationLevel(lfoALevel);
		const float lfoFValue = NormalizeVoiceLFOModulationLevel(lfoFLevel);
		// Keep tremolo centered around unity so positive and negative swings can
		// both attenuate and boost the decoded sample; clamping at 2.0 is a
		// practical host-side limit that preserves the full signed guest range
		// without allowing runaway amplification in the software mix path.
		const float amplitudeLFOModulation = std::clamp(1.0f + lfoAValue * lfoAmplitudeAmount, 0.0f, 2.0f);
		const float cutoffLFOOctaves = lfoAValue * lfoAmplitudeCutoffAmount;
		const double modulatedPitchStep = pitchStep * std::exp2(
			static_cast<double>(lfoAValue * lfoAmplitudePitchAmount + lfoFValue * lfoPitchAmount));
		const float totalGain = envelopeGain * amplitudeLFOModulation;
		if (captureVoiceDiagnostics) {
			diagnostics->maxAmplitudeLFOModulation = std::max(
				diagnostics->maxAmplitudeLFOModulation, amplitudeLFOModulation);
			diagnostics->maxPitchStep = std::max(diagnostics->maxPitchStep, modulatedPitchStep);
		}

		uint32_t activeState = 0;
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, activeState) ||
			(activeState & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0) {
			break;
		}

		float currentLeft = 0.0f;
		float currentRight = 0.0f;
		const uint32_t frameOffset = multipass ? static_cast<uint32_t>(frame) : currentOffset;
		if (!decodeFrame(baseAddress, frameOffset, currentLeft, currentRight)) {
			// stopVoice clears the cached playback state and guest active bits; break out of
			// the per-frame loop immediately, then let the final voiceStopped return below
			// skip the remaining state writes for this voice.
			failVoice();
			break;
		}
		if (multipass) {
			applyLowPass(currentLeft, currentRight, filterEnvelopeGain, cutoffLFOOctaves);
			if (capture3DHandoff) {
				storeCaptured3DSample(frame, currentLeft, currentRight, totalGain);
			}
			if (hrtfEnabled) {
				ProcessHRTFSample(voiceHandle, currentLeft, currentRight);
			}
			mixSamples(currentLeft, currentRight, totalGain, frame);
			if (captureVoiceDiagnostics) {
				++diagnostics->framesRendered;
			}
			continue;
		}

		float nextLeft = currentLeft;
		float nextRight = currentRight;
		auto previewSSLData = sslData;
		uint32_t previewBaseAddress = baseAddress;
		uint32_t previewEndOffset = endOffset;
		uint32_t previewOffset = currentOffset + offsetStep;
		if (advancePlaybackPosition(previewSSLData, previewBaseAddress, previewEndOffset, previewOffset, false)) {
			if (!decodeFrame(previewBaseAddress, previewOffset, nextLeft, nextRight)) {
				// If the look-ahead sample cannot be decoded, keep the current sample so
				// interpolation continues briefly. Repeated look-ahead failures mean the next
				// sample is persistently unreadable, so stop the voice instead of looping on
				// the same broken preview forever.
				++playbackState.previewDecodeFailures;
				if (playbackState.previewDecodeFailures >= APU_MAX_CONSECUTIVE_PREVIEW_DECODE_FAILURES) {
					failVoice();
					break;
				}
				nextLeft = currentLeft;
				nextRight = currentRight;
			} else {
				playbackState.previewDecodeFailures = 0;
			}
		}

		const float interpolation = static_cast<float>(playbackState.fraction);
		float sampleLeft = currentLeft + (nextLeft - currentLeft) * interpolation;
		float sampleRight = currentRight + (nextRight - currentRight) * interpolation;
		applyLowPass(sampleLeft, sampleRight, filterEnvelopeGain, cutoffLFOOctaves);
		if (capture3DHandoff) {
			storeCaptured3DSample(frame, sampleLeft, sampleRight, totalGain);
		}
		if (hrtfEnabled) {
			ProcessHRTFSample(voiceHandle, sampleLeft, sampleRight);
		}

		mixSamples(sampleLeft, sampleRight, totalGain, frame);
		if (captureVoiceDiagnostics) {
			++diagnostics->framesRendered;
		}

		const double nextPlaybackPosition = playbackState.fraction + modulatedPitchStep;
		const uint32_t wholeFrames = static_cast<uint32_t>(nextPlaybackPosition);
		playbackState.fraction = nextPlaybackPosition - static_cast<double>(wholeFrames);
		for (uint32_t step = 0; step < wholeFrames; ++step) {
			currentOffset += offsetStep;
			if (captureVoiceDiagnostics) {
				diagnostics->offsetAdvance += offsetStep;
			}
			if (!advancePlaybackPosition(sslData, baseAddress, endOffset, currentOffset, true)) {
				frame = frameCount;
				break;
			}
		}
	}
	if (shouldSkipStateWrites) {
		return;
	}

	if (submit3DToAC97) {
		++m_ChunkSubmittedHostSpatialVoiceCount;
		g_AC97->Submit3DVoiceFrames(voiceHandle, hrtfEntryIndex, guestHRTFPan, stereo,
			m_VPHRTFSubmix, hrtfSubmixVolumes, m_VPHRTFHeadroom,
			m_VP3DVoiceCaptureScratch.data(), frameCount);
	}

	if (multipass) {
		return;
	}

	m_VPSSLData[voiceHandle] = sslData;
	playbackState.offset = currentOffset;
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_CBO, currentOffset);
	WriteNotifierValue(voiceHandle, MCPX_HW_NOTIFIER_VOICE_POSITION, currentOffset);
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
	SetRegister32(APU_VP_BASE + APU_VP_FREE, GetVPFifoFreeSlots());
}

uint32_t APUDevice::GetVPFifoFreeSlots() const
{
	return m_VPFifoLevel >= APU_VP_FIFO_CAPACITY
		? 0u
		: (APU_VP_FIFO_CAPACITY - m_VPFifoLevel);
}

void APUDevice::RefreshInterruptStatus()
{
	uint32_t status = GetRegister32(NV_PAPU_ISTS) & ~NV_PAPU_ISTS_GINTSTS;
	if ((GetRegister32(NV_PAPU_FECTL) & NV_PAPU_FECTL_FEMETHMODE) == NV_PAPU_FECTL_FEMETHMODE_TRAPPED) {
		status |= NV_PAPU_ISTS_FETINTSTS;
	}

	if ((GetRegister32(NV_PAPU_IEN) & NV_PAPU_ISTS_GINTSTS) != 0 &&
		((status & ~NV_PAPU_ISTS_GINTSTS) & GetRegister32(NV_PAPU_IEN)) != 0) {
		status |= NV_PAPU_ISTS_GINTSTS;
	}

	SetRegister32(NV_PAPU_ISTS, status);
	HalSystemInterrupts[APU_IRQ].Assert((status & NV_PAPU_ISTS_GINTSTS) != 0);
}
