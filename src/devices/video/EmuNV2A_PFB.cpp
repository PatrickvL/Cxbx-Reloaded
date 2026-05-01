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
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a_pfb.c
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

DEVICE_READ32(PFB)
{
	DEVICE_READ32_SWITCH() {
	case NV_PFB_CFG0:
        /* 3-4 memory partitions. The debug bios checks this. */
		result = 3; // = NV_PFB_CFG0_PART_4
		break;
	case NV_PFB_CSTATUS:
		result = d->vram_size;
		break;
	case NV_PFB_WBC:
		result = 0; // = !NV_PFB_WBC_FLUSH /* Flush not pending. */
		break;
	default:
		DEVICE_READ32_REG(pfb);
		break;
	}

	DEVICE_READ32_END(PFB);
}

// TODO: Remove disabled warning once case are add to PFB switch.
#pragma warning(push)
#pragma warning(disable: 4065)
DEVICE_WRITE32(PFB)
{
	switch (addr) {
	default:
		DEVICE_WRITE32_REG(pfb);

		// Mirror tile region registers to PGRAPH (real hardware does this automatically).
		// NV_PFB_TILE[0..7] at 0x240..0x2BF map to NV_PGRAPH_TTILE[0..7] at 0x900..0x97F.
		if (addr >= 0x240 && addr < 0x2C0) {
			unsigned int pgraph_addr = NV_PGRAPH_TTILE(0) + (addr - 0x240);
			d->pgraph.regs[RI(pgraph_addr)] = value;
		}
		// NV_PFB_ZCOMP[0..7] at 0x300..0x31F map to NV_PGRAPH_ZCOMP[0..7] at 0x980..0x99F.
		// NV_PFB_ZCOMP_OFFSET (0x324) maps to NV_PGRAPH_ZCOMP_OFFSET (0x9A0).
		else if (addr >= 0x300 && addr <= 0x324) {
			unsigned int pgraph_addr = NV_PGRAPH_ZCOMP(0) + (addr - 0x300);
			d->pgraph.regs[RI(pgraph_addr)] = value;
		}
		break;
	}

	DEVICE_WRITE32_END(PFB);
}
#pragma warning(pop)
