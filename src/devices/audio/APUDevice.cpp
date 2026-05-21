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
#include "common/AddressRanges.h"
#include "common/audio/XADPCM.h"
#include "core/kernel/support/Emu.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#define LOG_PREFIX CXBXR_MODULE::MCPX

namespace {

constexpr uint32_t APU_VP_BASE = 0x20000;
constexpr uint32_t APU_VP_SIZE = 0x10000;
constexpr uint32_t APU_VP_FREE = 0x10;
constexpr uint32_t APU_VP_FIFO_CAPACITY = 0x80;
constexpr uint32_t APU_VP_STATUS_EMPTY = 0x80;
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
constexpr uint32_t NV_PAVS_VOICE_CFG_ENVF = 0x00000014;
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
constexpr uint32_t NV_PAVS_VOICE_CUR_PSL_START_BA = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_ECNT_EACOUNT = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_CUR_ECNT_EFCOUNT = 0xFFFF0000;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET_CBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_OFFSET_EALVL = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT_EBO = 0x00FFFFFF;
constexpr uint32_t NV_PAVS_VOICE_PAR_NEXT_EFLVL = 0xFF000000;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH = 0xFFFF0000;
constexpr uint32_t NV_PAVS_VOICE_TAR_LFO_ENV_EA_RELEASERATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE = 0x00000FFF;
constexpr uint32_t NV_PAVS_VOICE_CFG_MISC_FMODE = 0x00030000;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCA = 0x00000074;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCB = 0x00000078;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCA_FC0 = 0x0000FFFF;
constexpr uint32_t NV_PAVS_VOICE_TAR_FCA_FC1 = 0xFFFF0000;

constexpr uint32_t APU_VOICE_LIST_INHERIT = 0;
constexpr uint32_t APU_SGE_PAGE_SIZE = 0x1000;
constexpr size_t APU_AUDIO_CHUNK_FRAMES = 256;
constexpr float APU_VOLUME_DECIBEL_DIVISOR = 64.0f * -20.0f;
constexpr uint32_t MCPX_HW_NOTIFIER_BASE_OFFSET = 16;
constexpr uint32_t MCPX_HW_NOTIFIER_COUNT = 16;
constexpr uint32_t MCPX_HW_NOTIFIER_SSLA_DONE = 0;
constexpr uint32_t MCPX_HW_NOTIFIER_SSLB_DONE = 1;
constexpr uint8_t NV1BA0_NOTIFICATION_STATUS_DONE_SUCCESS = 0xFF;
constexpr uint8_t APU_NOTIFY_ENV_STATE_ACTIVE = 1;
constexpr double APU_PITCH_STEP_EXPONENT = 4096.0;
constexpr size_t APU_XADPCM_PCM_SAMPLES_PER_BLOCK = XBOX_ADPCM_DSTSIZE / sizeof(int16_t);
constexpr size_t APU_XADPCM_MAX_CHANNELS = 2;
constexpr size_t APU_XADPCM_MAX_SOURCE_BLOCK_BYTES = XBOX_ADPCM_SRCSIZE * APU_XADPCM_MAX_CHANNELS;
constexpr size_t APU_XADPCM_MAX_DECODED_SAMPLES = APU_XADPCM_PCM_SAMPLES_PER_BLOCK * APU_XADPCM_MAX_CHANNELS;
constexpr size_t APU_MIXBIN_COUNT = 32;
constexpr uint32_t APU_MAX_3D_VOICES = static_cast<uint32_t>(APUDevice::MAX_HRTF_VOICES);
constexpr size_t APU_HRTF_SUBMIX_COUNT = 4;
static_assert(APU_HRTF_SUBMIX_COUNT <= 8, "HRTF submix handoff expects no more than eight voice volumes");
constexpr size_t APU_HRTF_ENTRY_COUNT = 128;
constexpr size_t APU_HRTF_COEFFICIENT_COUNT = APUDevice::HRTF_FILTER_TAPS;
constexpr uint32_t APU_INVALID_HRTF_ENTRY_INDEX = 0xFFFF;
constexpr float APU_HRTF_ITD_SCALE = 512.0f;
constexpr float APU_HRTF_PARAM_SMOOTH_ALPHA = 0.01f;
constexpr float APU_HRTF_NORMALIZATION_EPSILON = 0.000001f;
constexpr float APU_HRTF_MAX_DELAY_SAMPLES_FLOAT = static_cast<float>(APUDevice::HRTF_FILTER_DELAY_SAMPLES);
// Scale normalized floating-point samples to signed 16-bit PCM amplitude.
constexpr float APU_SAMPLE_SCALE_FACTOR = 32767.0f;
constexpr size_t APU_DIAGNOSTIC_MAX_VOICES_TO_LOG = 4;
constexpr size_t APU_DIAGNOSTIC_MAX_ACTIVE_VOICES_TO_LOG = 4;

uint32_t GetFEMethodTargetVoiceOrDefault(uint32_t addr, uint32_t value, uint32_t currentVoiceValue)
{
	switch (addr) {
	case NV1BA0_PIO_SET_CURRENT_VOICE:
		return value & APU_VP_VOICE_MAX_HANDLE;
	case NV1BA0_PIO_VOICE_ON:
	case NV1BA0_PIO_VOICE_OFF:
	case NV1BA0_PIO_VOICE_RELEASE:
	case NV1BA0_PIO_GET_VOICE_POSITION:
	case NV1BA0_PIO_VOICE_PAUSE:
	case NV1BA0_PIO_SET_CURRENT_HRTF_ENTRY:
		return value & APU_VP_VOICE_MAX_HANDLE;
	case NV1BA0_PIO_VOICE_LOCK:
	case NV1BA0_PIO_SET_CONTEXT_DMA_NOTIFY:
	case NV1BA0_PIO_SET_CURRENT_SSL_CONTEXT_DMA:
	case NV1BA0_PIO_SET_CURRENT_SSL:
	case NV1BA0_PIO_SET_HRTF_SUBMIXES:
	case NV1BA0_PIO_SET_HRTF_HEADROOM:
	case NV1BA0_PIO_SET_VOICE_CFG_VBIN:
	case NV1BA0_PIO_SET_VOICE_CFG_FMT:
	case NV1BA0_PIO_SET_VOICE_CFG_ENV0:
	case NV1BA0_PIO_SET_VOICE_CFG_ENVA:
	case NV1BA0_PIO_SET_VOICE_CFG_ENV1:
	case NV1BA0_PIO_SET_VOICE_CFG_ENVF:
	case NV1BA0_PIO_SET_VOICE_CFG_MISC:
	case NV1BA0_PIO_SET_VOICE_TAR_HRTF:
	case NV1BA0_PIO_SET_VOICE_SSL_A:
	case NV1BA0_PIO_SET_VOICE_SSL_B:
	case NV1BA0_PIO_SET_VOICE_TAR_VOLA:
	case NV1BA0_PIO_SET_VOICE_TAR_VOLB:
	case NV1BA0_PIO_SET_VOICE_TAR_VOLC:
	case NV1BA0_PIO_SET_VOICE_LFO_ENV:
	case NV1BA0_PIO_SET_VOICE_TAR_FCA:
	case NV1BA0_PIO_SET_VOICE_TAR_FCB:
	case NV1BA0_PIO_SET_VOICE_TAR_PITCH:
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_BASE:
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_LBO:
	case NV1BA0_PIO_SET_VOICE_BUF_CBO:
	case NV1BA0_PIO_SET_VOICE_CFG_BUF_EBO:
		return currentVoiceValue & APU_VP_VOICE_MAX_HANDLE;
	default:
		return APU_VP_VOICE_MAX_HANDLE;
	}
}
// Match xemu's VP filter bounds: hardware-style cutoff is clamped to 2^-8..1.0.
constexpr float APU_FILTER_MIN_FREQUENCY = 0.003906f;
// Match xemu's minimum stable SVF resonance derived from the MCPX FC1 range.
constexpr float APU_FILTER_MIN_Q = 0.079407f;
// FC1 is a 16-bit fixed-point resonance value normalized against 0x8000.
constexpr float APU_FILTER_Q_NORMALIZER = 32768.0f;

uint32_t AbsoluteMixMagnitude(int32_t value)
{
	const int64_t signedSample = static_cast<int64_t>(value);
	const uint64_t magnitude = signedSample < 0
		? static_cast<uint64_t>(-signedSample)
		: static_cast<uint64_t>(signedSample);
	return static_cast<uint32_t>(magnitude);
}

uint32_t PeakAbsoluteMixAmplitude(const int32_t* samples, size_t sampleCount)
{
	uint32_t peak = 0;
	for (size_t i = 0; i < sampleCount; ++i) {
		peak = std::max(peak, AbsoluteMixMagnitude(samples[i]));
	}
	return peak;
}

uint32_t PeakAbsoluteMixBinAmplitude(const int32_t* mixBins, size_t frameCount, size_t slot)
{
	return PeakAbsoluteMixAmplitude(mixBins + slot * frameCount, frameCount);
}

float ClampUnitSample(float value)
{
	return std::clamp(value, -1.0f, 1.0f);
}

float RunLowPassFilter(float& high, float& band, float& low, float cutoff, float resonance, float input)
{
	// State-variable low-pass filter adapted to the lightweight MCPX VP path.
	// The small bias and cubic damping terms are the same stabilizers used in xemu's SVF implementation.
	const float normalizedInput = std::sqrt(resonance / 2.0f + 0.01f) * input;
	band -= band * band * band * 0.001f;
	high = normalizedInput - low - resonance * band;
	band += cutoff * high;
	low += cutoff * band;
	return low;
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

uint32_t ReadRegisterFragment(uint32_t value, uint32_t byteOffset, unsigned size)
{
	const uint32_t shift = byteOffset * 8;
	if (size >= sizeof(uint32_t)) {
		return value;
	}

	const uint32_t mask = (1u << (size * 8)) - 1;
	return (value >> shift) & mask;
}

uint32_t Ctz32(uint32_t value)
{
	uint32_t shift = 0;
	while (((value >> shift) & 1u) == 0u && shift < 32) {
		++shift;
	}
	return shift;
}

bool IsGuestRangeAccessible(uint32_t guestAddress, uint32_t size)
{
	return size > 0 && guestAddress <= PHYSICAL_MAP_SIZE && size <= (PHYSICAL_MAP_SIZE - guestAddress);
}

float AttenuateVoiceVolume(uint32_t volume)
{
	const uint32_t clamped = volume & 0x0FFF;
	return clamped == 0x0FFF ? 0.0f : std::pow(10.0f, static_cast<float>(clamped) / APU_VOLUME_DECIBEL_DIVISOR);
}

static int32_t ApplyHeadroomToMixSample(int32_t sample, uint8_t headroom)
{
	if (headroom == 0) {
		return sample;
	}

	const int32_t rounding = 1 << (headroom - 1);
	if (sample >= 0) {
		return (sample + rounding) >> headroom;
	}

	return -(((-sample) + rounding) >> headroom);
}

float ConvertUnsigned8(uint8_t value)
{
	return (static_cast<float>(value) - 128.0f) / 128.0f;
}

float ConvertSigned16(int16_t value)
{
	return static_cast<float>(value) / 32768.0f;
}

float ConvertSigned24(uint32_t value)
{
	const int32_t extended = (static_cast<int32_t>(value << 8)) >> 8;
	return static_cast<float>(extended) / 8388608.0f;
}

float ConvertSigned32(int32_t value)
{
	return static_cast<float>(value) / 2147483648.0f;
}

double DecodePitchStep(uint32_t pitch)
{
	const int16_t signedPitch = static_cast<int16_t>(pitch & 0xFFFF);
	return std::exp2(static_cast<double>(signedPitch) / APU_PITCH_STEP_EXPONENT);
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
	float maxEnvelopeGain = 0.0f;
	uint32_t decodedPeak = 0;
	uint32_t mixedPeak = 0;
	bool visited = false;
	bool active = false;
	bool mixed = false;
	bool decodedNonZero = false;
	bool stereoContribution = false;
	bool nonStereoContribution = false;
};

extern AC97Device* g_AC97;

// Basic VP playback and guest-visible buffer plumbing exist here, but full
// GP/EP DSP execution and threaded audio scheduling are still incomplete.

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
	m_VPInputSgeHandle = 0;
	m_VPOutputSgeHandle = 0;
	m_VPNotifyContextDMA = 0;
	m_VPCurrentSSLContextDMA = 0;
	m_VPSSLBasePage = 0;
	m_VPCurrentHRTFEntry = 0;
	m_VPLastVoicePositionHandle = 0;
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
	m_VPOutBufferCursor.fill(0);
	m_VPSSLData.fill(APUDevice::SSLData{});
	m_VPPlaybackState.fill(APUDevice::PlaybackState{});
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
		return ReadRegisterFragment(GetAPUTime(), addr - NV_PAPU_XGSCNT, size);
	}

	return ReadRegister(addr, size);
}

void APUDevice::MMIOWrite(int barIndex, uint32_t addr, uint32_t value, unsigned size)
{
	(void)barIndex;
	SynchronizeAudio();

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

	if (addr >= NV_PAPU_FEMEMDATA && addr < NV_PAPU_FEMEMDATA + sizeof(uint32_t)) {
		WriteRegister(addr, value, size);
		WriteGuestWord(GetRegister32(NV_PAPU_FEMEMADDR), GetRegister32(NV_PAPU_FEMEMDATA));
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
}


uint32_t APUDevice::GPRead(uint32_t addr, unsigned size)
{
	if (addr >= NV_PAPU_GPXMEM && addr < NV_PAPU_GPXMEM + m_GPXMem.size()) {
		return ReadMemoryWindow(m_GPXMem.data(), m_GPXMem.size(), addr - NV_PAPU_GPXMEM, size);
	}
	if (addr >= NV_PAPU_GPMIXBUF && addr < NV_PAPU_GPMIXBUF + m_GPMixBuf.size()) {
		return ReadMemoryWindow(m_GPMixBuf.data(), m_GPMixBuf.size(), addr - NV_PAPU_GPMIXBUF, size);
	}
	if (addr >= NV_PAPU_GPYMEM && addr < NV_PAPU_GPYMEM + m_GPYMem.size()) {
		return ReadMemoryWindow(m_GPYMem.data(), m_GPYMem.size(), addr - NV_PAPU_GPYMEM, size);
	}
	if (addr >= NV_PAPU_GPPMEM && addr < NV_PAPU_GPPMEM + m_GPPMem.size()) {
		return ReadMemoryWindow(m_GPPMem.data(), m_GPPMem.size(), addr - NV_PAPU_GPPMEM, size);
	}
	return ReadRegister(APU_GP_BASE + addr, size);
}

void APUDevice::GPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	if (addr >= NV_PAPU_GPXMEM && addr < NV_PAPU_GPXMEM + m_GPXMem.size()) {
		WriteMemoryWindow(m_GPXMem.data(), m_GPXMem.size(), addr - NV_PAPU_GPXMEM, value, size);
		return;
	}
	if (addr >= NV_PAPU_GPMIXBUF && addr < NV_PAPU_GPMIXBUF + m_GPMixBuf.size()) {
		WriteMemoryWindow(m_GPMixBuf.data(), m_GPMixBuf.size(), addr - NV_PAPU_GPMIXBUF, value, size);
		return;
	}
	if (addr >= NV_PAPU_GPYMEM && addr < NV_PAPU_GPYMEM + m_GPYMem.size()) {
		WriteMemoryWindow(m_GPYMem.data(), m_GPYMem.size(), addr - NV_PAPU_GPYMEM, value, size);
		return;
	}
	if (addr >= NV_PAPU_GPPMEM && addr < NV_PAPU_GPPMEM + m_GPPMem.size()) {
		WriteMemoryWindow(m_GPPMem.data(), m_GPPMem.size(), addr - NV_PAPU_GPPMEM, value, size);
		return;
	}
	WriteRegister(APU_GP_BASE + addr, value, size);
	if (addr == NV_PAPU_GPRST && size == sizeof(uint32_t) &&
		(value & NV_PAPU_GPRST_GPRST) == 0) {
		m_GPXMem.fill(0);
		m_GPMixBuf.fill(0);
		m_GPYMem.fill(0);
		m_GPPMem.fill(0);
	}
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
	ConsumeVPMethod(addr, value, size);
	if (m_VPFifoLevel < APU_VP_FIFO_CAPACITY) {
		++m_VPFifoLevel;
	}
	RefreshVPStatus();
}


uint32_t APUDevice::EPRead(uint32_t addr, unsigned size)
{
	if (addr >= NV_PAPU_EPXMEM && addr < NV_PAPU_EPXMEM + m_EPXMem.size()) {
		return ReadMemoryWindow(m_EPXMem.data(), m_EPXMem.size(), addr - NV_PAPU_EPXMEM, size);
	}
	if (addr >= NV_PAPU_EPYMEM && addr < NV_PAPU_EPYMEM + m_EPYMem.size()) {
		return ReadMemoryWindow(m_EPYMem.data(), m_EPYMem.size(), addr - NV_PAPU_EPYMEM, size);
	}
	if (addr >= NV_PAPU_EPPMEM && addr < NV_PAPU_EPPMEM + m_EPPMem.size()) {
		return ReadMemoryWindow(m_EPPMem.data(), m_EPPMem.size(), addr - NV_PAPU_EPPMEM, size);
	}
	return ReadRegister(APU_EP_BASE + addr, size);
}

void APUDevice::EPWrite(uint32_t addr, uint32_t value, unsigned size)
{
	if (addr >= NV_PAPU_EPXMEM && addr < NV_PAPU_EPXMEM + m_EPXMem.size()) {
		WriteMemoryWindow(m_EPXMem.data(), m_EPXMem.size(), addr - NV_PAPU_EPXMEM, value, size);
		return;
	}
	if (addr >= NV_PAPU_EPYMEM && addr < NV_PAPU_EPYMEM + m_EPYMem.size()) {
		WriteMemoryWindow(m_EPYMem.data(), m_EPYMem.size(), addr - NV_PAPU_EPYMEM, value, size);
		return;
	}
	if (addr >= NV_PAPU_EPPMEM && addr < NV_PAPU_EPPMEM + m_EPPMem.size()) {
		WriteMemoryWindow(m_EPPMem.data(), m_EPPMem.size(), addr - NV_PAPU_EPPMEM, value, size);
		return;
	}
	WriteRegister(APU_EP_BASE + addr, value, size);
	if (addr == NV_PAPU_EPRST && size == sizeof(uint32_t) &&
		(value & NV_PAPU_GPRST_GPRST) == 0) {
		m_EPXMem.fill(0);
		m_EPYMem.fill(0);
		m_EPPMem.fill(0);
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
			NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 1);
		WriteVoiceMask(selectedHandle, NV_PAVS_VOICE_PAR_STATE,
			NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 1);
		m_VPSSLData[selectedHandle].ssl_index = 0;
		m_VPSSLData[selectedHandle].ssl_seg = 0;
		m_VPPlaybackState[selectedHandle] = PlaybackState{};
		ClearHRTFFilterState(selectedHandle);
		InitializeVoiceEnvelopes(selectedHandle, value);
		m_LoggedEmptyVoiceTableDiagnostics = false;
		return;
	}
	case NV1BA0_PIO_VOICE_OFF: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_OFF_HANDLE;
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 0);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
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
	case NV1BA0_PIO_VOICE_LOCK:
		SetVoiceLocked(currentVoice(), (value & 1u) != 0);
		return;
	case NV1BA0_PIO_VOICE_RELEASE: {
		const uint32_t voiceHandle = value & NV1BA0_PIO_VOICE_RELEASE_HANDLE;
		if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
			return;
		}
		BeginVoiceRelease(voiceHandle);
		return;
	}
	case NV1BA0_PIO_GET_VOICE_POSITION:
		m_VPLastVoicePositionHandle = value & NV1BA0_PIO_GET_VOICE_POSITION_HANDLE;
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
		}
		return;
	case NV1BA0_PIO_SET_VOICE_SSL_B:
		if (currentVoice() < m_VPSSLData.size()) {
			m_VPSSLData[currentVoice()].base[1] = (value & NV1BA0_PIO_SET_VOICE_SSL_A_BASE) >> Ctz32(NV1BA0_PIO_SET_VOICE_SSL_A_BASE);
			m_VPSSLData[currentVoice()].count[1] = static_cast<uint8_t>((value & NV1BA0_PIO_SET_VOICE_SSL_A_COUNT) >> Ctz32(NV1BA0_PIO_SET_VOICE_SSL_A_COUNT));
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
		if ((GetRegister32(NV_PAPU_FETFORCE1) & NV_PAPU_FETFORCE1_SE2FE_IDLE_VOICE) != 0) {
			uint32_t fectl = GetRegister32(NV_PAPU_FECTL);
			fectl &= ~(NV_PAPU_FECTL_FEMETHMODE | NV_PAPU_FECTL_FETRAPREASON);
			fectl |= NV_PAPU_FECTL_FEMETHMODE_TRAPPED | NV_PAPU_FECTL_FETRAPREASON_REQUESTED;
			SetRegister32(NV_PAPU_FECTL, fectl);
			RefreshInterruptStatus();
		}
		return;
	default:
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
			const uint32_t sslTableBase = GetRegister32(NV_PAPU_VPSSLADDR);
			if (sslTableBase != 0) {
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
			return;
		}
		if (addr >= NV1BA0_PIO_SET_OUTBUF_LEN && addr < NV1BA0_PIO_SET_OUTBUF_LEN + 0x20 && ((addr - NV1BA0_PIO_SET_OUTBUF_LEN) % 8) == 0) {
			const size_t slot = (addr - NV1BA0_PIO_SET_OUTBUF_LEN) / 8;
			WriteRegister(APU_VP_BASE + addr, value & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE, sizeof(uint32_t));
			if (slot < m_VPOutBufferCursor.size()) {
				m_VPOutBufferCursor[slot] = 0;
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
	if (dest == nullptr || !IsGuestRangeAccessible(guestAddress, static_cast<uint32_t>(size))) {
		return false;
	}

	std::memcpy(dest, reinterpret_cast<const void*>(static_cast<uintptr_t>(CONTIGUOUS_MEMORY_BASE + guestAddress)), size);
	return true;
}

bool APUDevice::WriteGuestWord(uint32_t guestAddress, uint32_t value)
{
	return WriteGuestBytes(guestAddress, &value, sizeof(value));
}

bool APUDevice::WriteGuestBytes(uint32_t guestAddress, const void* src, size_t size)
{
	if (src == nullptr || !IsGuestRangeAccessible(guestAddress, static_cast<uint32_t>(size))) {
		return false;
	}

	std::memcpy(reinterpret_cast<void*>(static_cast<uintptr_t>(CONTIGUOUS_MEMORY_BASE + guestAddress)), src, size);
	return true;
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

bool APUDevice::ReadVoiceMask(uint32_t voiceHandle, uint32_t offset, uint32_t mask, uint32_t& value) const
{
	if (voiceHandle >= APU_VP_VOICE_MAX_HANDLE) {
		return false;
	}

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		if (!m_LoggedVoiceTableReadFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU ReadVoiceMask blocked voiceTableBase=0x00000000 handle=%u offset=0x%08x mask=0x%08x",
				voiceHandle,
				offset,
				mask);
			m_LoggedVoiceTableReadFailure = true;
		}
		return false;
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

	const uint32_t shift = mask == 0xFFFFFFFF ? 0 : Ctz32(mask);
	value = mask == 0xFFFFFFFF ? current : ((current & mask) >> shift);
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

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		if (!m_LoggedVoiceTableWriteFailure) {
			EmuLog(LOG_LEVEL::WARNING,
				"APU WriteVoiceMask blocked voiceTableBase=0x00000000 handle=%u offset=0x%08x mask=0x%08x value=0x%08x",
				voiceHandle,
				offset,
				mask,
				value);
			m_LoggedVoiceTableWriteFailure = true;
		}
		return false;
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

void APUDevice::WriteNotifierStatus(uint32_t voiceHandle, uint32_t notifier, uint8_t status)
{
	const uint32_t notifierBase = GetRegister32(NV_PAPU_FENADDR);
	if (notifierBase != 0) {
		const uint32_t offset = 16 * (MCPX_HW_NOTIFIER_BASE_OFFSET + voiceHandle * MCPX_HW_NOTIFIER_COUNT + notifier);
		WriteGuestBytes(notifierBase + offset + 14, &APU_NOTIFY_ENV_STATE_ACTIVE, sizeof(APU_NOTIFY_ENV_STATE_ACTIVE));
		WriteGuestBytes(notifierBase + offset + 15, &status, sizeof(status));
	}

	SetRegister32(NV_PAPU_ISTS, GetRegister32(NV_PAPU_ISTS) | NV_PAPU_ISTS_FEVINTSTS | NV_PAPU_ISTS_FENINTSTS);
	RefreshInterruptStatus();
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

bool APUDevice::ResolveVoiceAddress(uint32_t linearAddress, uint32_t& guestAddress) const
{
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
		uint32_t count = 0;
		uint32_t level = 0xFF;

		switch (startState) {
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_OFF:
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_DELAY:
			ReadVoiceMask(voiceHandle, reg0, NV_PAVS_VOICE_CFG_ENV0_EA_DELAYTIME, count);
			count *= 16;
			level = 0;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK:
			level = 0;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD:
			ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_HOLDTIME, count);
			count *= 16;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY:
			ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_DECAYRATE, count);
			count *= 16;
			break;
		case NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN:
			ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_SUSTAINLEVEL, level);
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
}

float APUDevice::StepVoiceEnvelope(uint32_t voiceHandle, uint32_t reg0, uint32_t regA,
	uint32_t rrReg, uint32_t rrMask, uint32_t levelRegister, uint32_t levelMask,
	uint32_t countMask, uint32_t stateMask)
{
	uint32_t currentState = 0;
	if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask, currentState)) {
		return 1.0f;
	}

	const bool amplitudeEnvelope = countMask == NV_PAVS_VOICE_CUR_ECNT_EACOUNT;
	const auto stopVoice = [&]() {
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 0);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
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
		} else {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count - 1);
		}
		return 0.0f;
	}
	case NV_PAVS_VOICE_PAR_STATE_EFCUR_ATTACK: {
		uint32_t count = 0;
		uint32_t attackRate = 0;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_ECNT, countMask, count);
		ReadVoiceMask(voiceHandle, reg0, NV_PAVS_VOICE_CFG_ENV0_EA_ATTACKRATE, attackRate);

		const uint32_t attackSpan = attackRate * 16;
		uint32_t level = 0xFF;
		if (attackRate != 0 && attackSpan != 0) {
			level = std::min<uint32_t>(0xFF, static_cast<uint32_t>((count * 0xFFu) / attackSpan));
		}
		WriteVoiceMask(voiceHandle, levelRegister, levelMask, level);

		if (attackRate == 0 || count >= attackSpan) {
			uint32_t holdTime = 0;
			ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_HOLDTIME, holdTime);
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_HOLD);
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
			ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_DECAYRATE, decayRate);
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_DECAY);
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
		ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_DECAYRATE, decayRate);
		ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_SUSTAINLEVEL, sustainLevel);

		if (decayRate == 0 || count == 0) {
			WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, stateMask,
				NV_PAVS_VOICE_PAR_STATE_EFCUR_SUSTAIN);
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
		ReadVoiceMask(voiceHandle, regA, NV_PAVS_VOICE_CFG_ENVA_EA_SUSTAINLEVEL, sustainLevel);
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
	const uint32_t now = GetAPUTime();
	uint32_t remaining = now - m_LastAudioUpdate;
	while (remaining > 0) {
		const size_t chunk = std::min<size_t>(remaining, APU_AUDIO_CHUNK_FRAMES);
		RenderBasicAudioChunk(chunk);
		m_LastAudioUpdate += static_cast<uint32_t>(chunk);
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
			EmuLog(LOG_LEVEL::WARNING,
				"APU render blocked because NV_PAPU_VPVADDR is still zero; FE/VP voice state is not wired to guest memory yet");
			m_LoggedMissingVoiceTableDuringRender = true;
		}
	} else {
		m_LoggedMissingVoiceTableDuringRender = false;
	}

	if (g_AC97 != nullptr) {
		g_AC97->Begin3DVoiceFrameBatch();
	}

	std::vector<int32_t> mixBins(frameCount * APU_MIXBIN_COUNT, 0);
	const size_t visited2D = RenderBasicVoiceList(NV_PAPU_TVL2D, mixBins.data(), frameCount);
	const size_t visited3D = RenderBasicVoiceList(NV_PAPU_TVL3D, mixBins.data(), frameCount);
	const size_t visitedMP = RenderBasicVoiceList(NV_PAPU_TVLMP, mixBins.data(), frameCount);
	const bool hasVoiceActivity = (visited2D + visited3D + visitedMP) != 0;
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		if (!hasVoiceActivity) {
			if (!m_LoggedEmptyVoiceTableDiagnostics) {
				LogVoiceTableDiagnostics();
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
			outBufferPeak[slot] = PeakAbsoluteMixBinAmplitude(mixBins.data(), frameCount, slot);
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
	for (size_t frame = 0; frame < frameCount; ++frame) {
		output[frame * 2] = ClampToInt16(mixBins[frame]);
		output[frame * 2 + 1] = ClampToInt16(mixBins[frameCount + frame]);
	}

	uint32_t stereoPeak = 0;
	if constexpr (audio_diagnostics::kEnableDiagnosticLogging) {
		stereoPeak = audio_diagnostics::PeakAbsoluteSampleAmplitude(output.data(), output.size());
		if (hasVoiceActivity || stereoPeak != 0) {
			EmuLog(LOG_LEVEL::INFO,
				"APU stereo mix peak before SubmitPCMFrames=%u frames=%zu",
				static_cast<unsigned>(stereoPeak),
				frameCount);
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
		const uint32_t outBufferLength = outBufferLengthRegister & NV1BA0_PIO_SET_OUTBUF_LEN_VALUE;
		if (outBufferBase == 0 || outBufferLength < sizeof(int16_t)) {
			m_VPOutBufferCursor[slot] = 0;
			continue;
		}

		for (size_t frame = 0; frame < frameCount; ++frame) {
			output[frame] = ClampToInt16(mixBins[slot * frameCount + frame]);
		}

		WriteGuestCircularBuffer(outBufferBase, outBufferLength, m_VPOutBufferCursor[slot],
			output.data(), output.size() * sizeof(output[0]));
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
		return diagnostics.active &&
			((diagnostics.decodedNonZero && !diagnostics.stereoContribution) ||
			 !diagnostics.decodedNonZero ||
			 (diagnostics.framesRendered != 0 && diagnostics.offsetAdvance == 0 && diagnostics.pitchStep > 0.0));
	};
	for (size_t visited = 0; visited < 1024 && voiceHandle < APU_VP_VOICE_MAX_HANDLE; ++visited) {
		uint32_t nextHandle = APU_VP_VOICE_MAX_HANDLE;
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_PITCH_LINK,
			NV_PAVS_VOICE_TAR_PITCH_LINK_NEXT_VOICE_HANDLE, nextHandle);
		BasicVoiceDiagnosticSummary diagnostics;
		RenderBasicVoice(voiceHandle, mixBins, frameCount, &diagnostics);
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
					"APU voice diag handle=%u offsets=%u+%u frames=%zu pitchStep=%.6f decodedPeak=%u mixedPeak=%u envMax=%.3f stereo=%d other=%d bins=[%u,%u,%u,%u,%u,%u,%u,%u] volumes=[%u,%u,%u,%u,%u,%u,%u,%u] headroom=[%u,%u,%u,%u,%u,%u,%u,%u]",
					diagnostics.voiceHandle,
					diagnostics.startOffset,
					diagnostics.offsetAdvance,
					diagnostics.framesRendered,
					diagnostics.pitchStep,
					diagnostics.decodedPeak,
					diagnostics.mixedPeak,
					diagnostics.maxEnvelopeGain,
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

	const uint32_t voiceTableBase = GetRegister32(NV_PAPU_VPVADDR);
	if (voiceTableBase == 0) {
		EmuLog(LOG_LEVEL::INFO,
			"APU voice table diagnostics voiceTableBase=0x00000000 active=0 paused=0 new=0 handles=[none]");
		LogRecentFEMethodDiagnostics();
		return;
	}

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
			"APU voice table diagnostics voiceTableBase=0x%08x active=%zu paused=%zu new=%zu handles=[none]",
			voiceTableBase,
			activeVoiceCount,
			pausedVoiceCount,
			newVoiceCount);
		LogRecentFEMethodDiagnostics();
		return;
	}

	EmuLog(LOG_LEVEL::INFO,
		"APU voice table diagnostics voiceTableBase=0x%08x active=%zu paused=%zu new=%zu handles=[%u,%u,%u,%u]",
		voiceTableBase,
		activeVoiceCount,
		pausedVoiceCount,
		newVoiceCount,
		activeHandles[0],
		loggedActiveHandles > 1 ? activeHandles[1] : APU_VP_VOICE_MAX_HANDLE,
		loggedActiveHandles > 2 ? activeHandles[2] : APU_VP_VOICE_MAX_HANDLE,
		loggedActiveHandles > 3 ? activeHandles[3] : APU_VP_VOICE_MAX_HANDLE);
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
		(state & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0 ||
		(state & NV_PAVS_VOICE_PAR_STATE_PAUSED) != 0) {
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
	} else if (samplesPerBlock != 1) {
		return;
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
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_PSL_START, NV_PAVS_VOICE_CUR_PSL_START_BA, baseAddress);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_CBO, currentOffset);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_NEXT, NV_PAVS_VOICE_PAR_NEXT_EBO, endOffset);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CUR_PSH_SAMPLE, NV_PAVS_VOICE_CUR_PSH_SAMPLE_LBO, loopOffset);
	ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_TAR_PITCH_LINK, NV_PAVS_VOICE_TAR_PITCH_LINK_PITCH, pitch);
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
	const uint32_t bytesPerFrame = adpcm ? 0u : containerSize * channels;
	const uint32_t bytesPerBlock = containerSize * channels;
	auto& playbackState = m_VPPlaybackState[voiceHandle];
	if (!playbackState.valid || playbackState.offset != currentOffset) {
		playbackState.offset = currentOffset;
		playbackState.fraction = 0.0;
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
	auto stopVoice = [&]() {
		playbackState = PlaybackState{};
		m_VPLowPassState[voiceHandle] = {};
		ClearHRTFFilterState(voiceHandle);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE, 0);
		WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, NV_PAVS_VOICE_PAR_STATE_NEW_VOICE, 0);
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
			if (bins[binIndex] >= APU_MIXBIN_COUNT || volumes[binIndex] == 0) {
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
			const float headroomDivisor = static_cast<float>(1u << headroom);
			const float gain = AttenuateVoiceVolume(volumes[binIndex]) * envelopeGain /
				headroomDivisor;
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
	auto applyLowPass = [&](float& sampleLeft, float& sampleRight) {
		if (!lowPassEnabled) {
			return;
		}
		auto& filterState = m_VPLowPassState[voiceHandle];
		sampleLeft = ClampUnitSample(RunLowPassFilter(filterState[0].high, filterState[0].band, filterState[0].low,
			lowPassCutoff[0], lowPassResonance[0], sampleLeft));
		sampleRight = ClampUnitSample(RunLowPassFilter(filterState[1].high, filterState[1].band, filterState[1].low,
			lowPassCutoff[1], lowPassResonance[1], sampleRight));
	};
	uint32_t hrtfEntryIndex = APU_INVALID_HRTF_ENTRY_INDEX;
	const bool hrtfEnabled = voiceHandle < APU_MAX_3D_VOICES &&
		ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_CFG_HRTF_TARGET,
			NV_PAVS_VOICE_CFG_HRTF_TARGET_HANDLE, hrtfEntryIndex) &&
		hrtfEntryIndex < m_VPHRTFEntries.size();
	const bool capture3DForOpenAL = hrtfEnabled && g_AC97 != nullptr;
	if (capture3DForOpenAL) {
		m_VP3DVoiceCaptureScratch.assign(frameCount * 2, 0);
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
	auto readSampleBytes = [&](uint32_t sampleAddress, void* dest, size_t size) {
		return streaming ? ReadGuestBytes(sampleAddress, dest, size) : ReadVoiceBufferBytes(sampleAddress, dest, size);
	};
	auto loadStreamingSegment = [&](SSLData& voiceSSLData, uint32_t& segmentBaseAddress, uint32_t& segmentEndOffset, uint32_t& segmentCurrentOffset, bool commit) -> bool {
		for (size_t attempts = 0; attempts < 4; ++attempts) {
			if (voiceSSLData.ssl_index > 1) {
				voiceSSLData.ssl_index = 0;
			}
			const uint32_t sslIndex = voiceSSLData.ssl_index;
			if (voiceSSLData.count[sslIndex] == 0) {
				if (commit) {
					stopVoice();
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

			const uint32_t sslTableBase = GetRegister32(NV_PAPU_VPSSLADDR);
			if (sslTableBase == 0) {
				if (!m_LoggedStreamingSSLFailure) {
					EmuLog(LOG_LEVEL::WARNING,
						"APU streaming voice %u needs NV_PAPU_VPSSLADDR but the SSL table base is still zero",
						voiceHandle);
					m_LoggedStreamingSSLFailure = true;
				}
				return false;
			}
			m_LoggedStreamingSSLFailure = false;

			const uint32_t segmentPage = voiceSSLData.base[sslIndex] + static_cast<uint32_t>(voiceSSLData.ssl_seg);
			uint32_t segmentOffset = 0;
			uint32_t segmentLength = 0;
			if (!ReadGuestWord(sslTableBase + segmentPage * 8, segmentOffset) ||
				!ReadGuestWord(sslTableBase + segmentPage * 8 + 4, segmentLength)) {
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
				return false;
			}

			segmentBaseAddress = segmentOffset;
			segmentEndOffset = segmentSamples - 1;
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

		const uint32_t linearAddress = segmentBaseAddress + segmentCurrentOffset * bytesPerFrame;
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
		(void)StepVoiceEnvelope(
			voiceHandle,
			NV_PAVS_VOICE_CFG_ENV1, NV_PAVS_VOICE_CFG_ENVF,
			NV_PAVS_VOICE_CFG_MISC, NV_PAVS_VOICE_CFG_MISC_EF_RELEASERATE,
			NV_PAVS_VOICE_PAR_NEXT, NV_PAVS_VOICE_PAR_NEXT_EFLVL,
			NV_PAVS_VOICE_CUR_ECNT_EFCOUNT, NV_PAVS_VOICE_PAR_STATE_EFCUR);

		uint32_t activeState = 0;
		if (!ReadVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_STATE, 0xFFFFFFFF, activeState) ||
			(activeState & NV_PAVS_VOICE_PAR_STATE_ACTIVE_VOICE) == 0) {
			break;
		}

		float currentLeft = 0.0f;
		float currentRight = 0.0f;
		const uint32_t frameOffset = multipass ? static_cast<uint32_t>(frame) : currentOffset;
		if (!decodeFrame(baseAddress, frameOffset, currentLeft, currentRight)) {
			return;
		}
		if (multipass) {
			applyLowPass(currentLeft, currentRight);
			if (capture3DForOpenAL) {
				storeCaptured3DSample(frame, currentLeft, currentRight, envelopeGain);
			}
			if (hrtfEnabled) {
				ProcessHRTFSample(voiceHandle, currentLeft, currentRight);
			}
			mixSamples(currentLeft, currentRight, envelopeGain, frame);
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
		uint32_t previewOffset = currentOffset + 1;
		if (advancePlaybackPosition(previewSSLData, previewBaseAddress, previewEndOffset, previewOffset, false)) {
			if (!decodeFrame(previewBaseAddress, previewOffset, nextLeft, nextRight)) {
				return;
			}
		}

		const float interpolation = static_cast<float>(playbackState.fraction);
		float sampleLeft = currentLeft + (nextLeft - currentLeft) * interpolation;
		float sampleRight = currentRight + (nextRight - currentRight) * interpolation;
		applyLowPass(sampleLeft, sampleRight);
		if (capture3DForOpenAL) {
			storeCaptured3DSample(frame, sampleLeft, sampleRight, envelopeGain);
		}
		if (hrtfEnabled) {
			ProcessHRTFSample(voiceHandle, sampleLeft, sampleRight);
		}

		mixSamples(sampleLeft, sampleRight, envelopeGain, frame);
		if (captureVoiceDiagnostics) {
			++diagnostics->framesRendered;
		}

		const double nextPlaybackPosition = playbackState.fraction + pitchStep;
		const uint32_t wholeFrames = static_cast<uint32_t>(nextPlaybackPosition);
		playbackState.fraction = nextPlaybackPosition - static_cast<double>(wholeFrames);
		for (uint32_t step = 0; step < wholeFrames; ++step) {
			++currentOffset;
			if (captureVoiceDiagnostics) {
				++diagnostics->offsetAdvance;
			}
			if (!advancePlaybackPosition(sslData, baseAddress, endOffset, currentOffset, true)) {
				frame = frameCount;
				break;
			}
		}
	}

	if (capture3DForOpenAL) {
		g_AC97->Submit3DVoiceFrames(voiceHandle, hrtfEntryIndex, stereo,
			m_VPHRTFSubmix, hrtfSubmixVolumes, m_VPHRTFHeadroom,
			m_VP3DVoiceCaptureScratch.data(), frameCount);
	}

	if (multipass) {
		return;
	}

	m_VPSSLData[voiceHandle] = sslData;
	playbackState.offset = currentOffset;
	WriteVoiceMask(voiceHandle, NV_PAVS_VOICE_PAR_OFFSET, NV_PAVS_VOICE_PAR_OFFSET_CBO, currentOffset);
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
}
