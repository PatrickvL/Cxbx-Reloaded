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
// *  https://github.com/xqemu/xqemu/blob/master/hw/xbox/nv2a/nv2a_prmvio.c
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

// TODO: Remove disabled warning once case are add to PRMVIO switch.
#pragma warning(push)
#pragma warning(disable: 4065)

DEVICE_READ32(PRMVIO)
{
	// VGA sequencer and graphics controller indexed I/O
	DEVICE_READ32_SWITCH() {
	case VGA_SEQ_I:
		result = d->prmvio.seq_index;
		break;
	case VGA_SEQ_D:
		result = d->prmvio.seq[d->prmvio.seq_index];
		break;
	case VGA_GFX_I:
		result = d->prmvio.gfx_index;
		break;
	case VGA_GFX_D:
		result = d->prmvio.gfx[d->prmvio.gfx_index];
		break;
	case VGA_MIS_R:
		result = d->prmvio.misc_output;
		break;
	default:
		DEBUG_READ32_UNHANDLED(PRMVIO);
		break;
	}

	DEVICE_READ32_END(PRMVIO);
}
#pragma warning(pop)

// TODO: Remove disabled warning once case are add to PRMVIO switch.
#pragma warning(push)
#pragma warning(disable: 4065)
DEVICE_WRITE32(PRMVIO)
{
	// VGA sequencer and graphics controller indexed I/O
	switch (addr) {
	case VGA_SEQ_I:
		d->prmvio.seq_index = value & 0xFF;
		break;
	case VGA_SEQ_D:
		d->prmvio.seq[d->prmvio.seq_index] = value & 0xFF;
		break;
	case VGA_GFX_I:
		d->prmvio.gfx_index = value & 0xFF;
		break;
	case VGA_GFX_D:
		d->prmvio.gfx[d->prmvio.gfx_index] = value & 0xFF;
		break;
	case VGA_MIS_W:
		d->prmvio.misc_output = value & 0xFF;
		break;
	default:
		DEBUG_WRITE32_UNHANDLED(PRMVIO);
		break;
	}

	DEVICE_WRITE32_END(PRMVIO);
}
#pragma warning(pop)
