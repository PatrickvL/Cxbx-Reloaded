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
// *  This file is heavily based on code from XQEMU
// *  https://github.com/xqemu/xqemu/tree/master/hw/xbox/nv2a/nv2a_pcrtc.c
// *  Copyright (c) 2012 espes
// *  Copyright (c) 2015 Jannik Vogel
// *  Copyright (c) 2018 Matt Borgerson
// *
// *  Contributions for Cxbx-Reloaded
// *  Copyright (c) 2017-2018 Luke Usher <luke.usher@outlook.com>
// *  Copyright (c) 2018 Patrick van Logchem <pvanlogchem@gmail.com>
// *
// *  All rights reserved
// *
// ******************************************************************

// Read Vertical Display End (visible scanlines) from VGA CRT registers.
// Returns the number of visible lines (e.g. 480 for NTSC, 576 for PAL).
// This is the same value GetFrameHeight() computes, but kept local to PCRTC.

// Needed for protecting VBlank enable when GPU ISR is connected
#include "core\kernel\exports\EmuKrnl.h"

static unsigned int pcrtc_get_visible_lines(NV2AState *d)
{
	unsigned int vde = ((unsigned int)d->prmcio.cr[NV_CIO_CR_VDE_INDEX])
		| (((unsigned int)d->prmcio.cr[NV_CIO_CR_OVL_INDEX] & 0x02) >> 1 << 8)
		| (((unsigned int)d->prmcio.cr[NV_CIO_CR_OVL_INDEX] & 0x40) >> 6 << 9)
		| (((unsigned int)d->prmcio.cr[NV_CIO_CRE_LSR_INDEX] & 0x02) >> 1 << 10)
		| (((unsigned int)d->prmcio.cr[NV_CIO_CRE_EBR_INDEX] & 0x04) >> 2 << 11);
	return vde + 1; // VDE is end value (0-based), so +1 for count
}

// Read Vertical Display Total from VGA CRT registers.
// Returns the total number of scanlines per frame including blanking.
static unsigned int pcrtc_get_total_lines(NV2AState *d)
{
	unsigned int vdt = ((unsigned int)d->prmcio.cr[NV_CIO_CR_VDT_INDEX])
		| (((unsigned int)d->prmcio.cr[NV_CIO_CR_OVL_INDEX] & 0x01) << 8)
		| (((unsigned int)d->prmcio.cr[NV_CIO_CR_OVL_INDEX] & 0x20) >> 5 << 9)
		| (((unsigned int)d->prmcio.cr[NV_CIO_CRE_LSR_INDEX] & 0x01) << 10)
		| (((unsigned int)d->prmcio.cr[NV_CIO_CRE_EBR_INDEX] & 0x01) << 11);
	return vdt + 2; // VDT register value is total-2
}

// Read Horizontal Display Total from VGA CRT registers (in character clocks).
// Returns the total number of character clocks per scanline including blanking.
static unsigned int pcrtc_get_htotal_chars(NV2AState *d)
{
	unsigned int hdt = ((unsigned int)d->prmcio.cr[NV_CIO_CR_HDT_INDEX])
		| (((unsigned int)d->prmcio.cr[NV_CIO_CRE_HEB__INDEX] & 0x01) << 8);
	return hdt + 5; // HDT register value is total-5
}

// Compute refresh rate (Hz) from VPLL pixel clock and CRT timing registers.
// refresh = pixel_clock / (htotal_pixels * vtotal_lines)
// where pixel_clock = (NV2A_CRYSTAL_FREQ * N) / (M * 2^P)
// and htotal_pixels = htotal_chars * 8 (each character clock = 8 pixels)
static unsigned int pcrtc_get_refresh_rate(NV2AState *d, unsigned int totalLines)
{
	uint32_t coeff = d->pramdac.video_clock_coeff;
	unsigned int m = coeff & NV_PRAMDAC_VPLL_COEFF_MDIV;
	unsigned int n = (coeff & NV_PRAMDAC_VPLL_COEFF_NDIV) >> 8;
	unsigned int p = (coeff & NV_PRAMDAC_VPLL_COEFF_PDIV) >> 16;
	if (m == 0 || totalLines == 0) return 60; // avoid division by zero

	// Pixel clock in Hz
	uint64_t pixelClock = ((uint64_t)NV2A_CRYSTAL_FREQ * n) / ((1 << p) * m);

	unsigned int htotalChars = pcrtc_get_htotal_chars(d);
	unsigned int htotalPixels = htotalChars * 8;
	if (htotalPixels == 0) return 60;

	unsigned int refreshRate = (unsigned int)(pixelClock / ((uint64_t)htotalPixels * totalLines));
	// Clamp to sane range (avoid garbage from uninitialized registers)
	if (refreshRate < 24 || refreshRate > 120) return 60;
	return refreshRate;
}

DEVICE_READ32(PCRTC)
{
	DEVICE_READ32_SWITCH() {

	case NV_PCRTC_INTR_0:
		result = d->pcrtc.pending_interrupts;
		break;
	case NV_PCRTC_INTR_EN_0:
		result = d->pcrtc.enabled_interrupts;
		break;
	case NV_PCRTC_START:
		result = d->pcrtc.start;
		break;
	case NV_PCRTC_RASTER: {
		// Test case: Alter Echo, FieldRender
		// Return scanline position relative to the last VBlank interrupt.
		// This ensures games polling NV_PCRTC_RASTER see timing consistent
		// with VBlank interrupt delivery (same time source).
		unsigned int visibleLines = pcrtc_get_visible_lines(d);
		unsigned int totalLines = pcrtc_get_total_lines(d);
		// Guard against uninitialized registers (early boot before AvSetDisplayMode)
		if (visibleLines == 0 || totalLines == 0) {
			visibleLines = 480;
			totalLines = 525;
		}
		unsigned int refreshRate = pcrtc_get_refresh_rate(d, totalLines);
		LARGE_INTEGER now;
		QueryPerformanceCounter(&now);
		// Frame period in QPC ticks (HostQPCFrequency from Timer.h, set once in timer_init)
		LONGLONG frameTicks = HostQPCFrequency / refreshRate;
		// Compute position relative to last VBlank (instead of free-running QPC modulo)
		int64_t lastVBlankQPC = d->vblank_last_qpc.load(std::memory_order_acquire);
		LONGLONG posInFrame;
		if (lastVBlankQPC > 0) {
			posInFrame = now.QuadPart - lastVBlankQPC;
			if (posInFrame < 0) posInFrame = 0;
			if (posInFrame >= frameTicks) posInFrame = frameTicks - 1;
		} else {
			// Fallback before first VBlank fires (early boot)
			posInFrame = now.QuadPart % frameTicks;
		}
		unsigned int scanline = (unsigned int)(posInFrame * totalLines / frameTicks);
		result = scanline & NV_PCRTC_RASTER_POSITION;
		// Bit 16: VERT_BLANK - active when scanline is in the blanking interval
		if (scanline >= visibleLines) {
			result |= NV_PCRTC_RASTER_VERT_BLANK;
		}
		// Bit 20: FIELD - toggles each VBlank for interlaced modes (0=EVEN, 1=ODD)
		if (d->pcrtc.vblank_count & 1) {
			result |= NV_PCRTC_RASTER_FIELD;
		}
	} break;
	default: 
		result = 0;
		//DEVICE_READ32_REG(pcrtc); // Was : DEBUG_READ32_UNHANDLED(PCRTC);
		break;
	}

	DEVICE_READ32_END(PCRTC);
}

DEVICE_WRITE32(PCRTC)
{
	switch (addr) {

	case NV_PCRTC_INTR_0:
		d->pcrtc.pending_interrupts &= ~value;
		update_irq(d);
		break;
	case NV_PCRTC_INTR_EN_0:
		d->pcrtc.enabled_interrupts = value;
		// On real Xbox, the kernel miniport owns NV_PCRTC_INTR_EN_0 and ensures
		// VBlank stays enabled after ISR connection. Since we emulate the miniport
		// at the API level (KeConnectInterrupt), the miniport's MMIO write never
		// executes. Preserve VBlank enable when the GPU ISR is connected to prevent
		// the D3D runtime's init sequence from accidentally disabling it.
		if (EmuInterruptList[3] && EmuInterruptList[3]->Connected) {
			d->pcrtc.enabled_interrupts |= NV_PCRTC_INTR_0_VBLANK;
		}
		update_irq(d);
		break;
	case NV_PCRTC_START:
        value &= 0x07FFFFFF;
        // assert(val < memory_region_size(d->vram));
		d->pcrtc.start = value;

        NV2A_DPRINTF("PCRTC_START - %x %x %x %x\n",
                d->vram_ptr[value+64], d->vram_ptr[value+64+1],
                d->vram_ptr[value+64+2], d->vram_ptr[value+64+3]);
		break;

	default: 
		DEVICE_WRITE32_REG(pcrtc); // Was : DEBUG_WRITE32_UNHANDLED(PCRTC);
		break;
	}

	DEVICE_WRITE32_END(PCRTC);
}
