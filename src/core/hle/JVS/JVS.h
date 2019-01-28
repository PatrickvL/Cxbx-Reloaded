// ******************************************************************
// *
// *    .,-:::::    .,::      .::::::::.    .,::      .:
// *  ,;;;'````'    `;;;,  .,;;  ;;;'';;'   `;;;,  .,;;
// *  [[[             '[[,,[['   [[[__[[\.    '[[,,[['
// *  $$$              Y$$$P     $$""""Y$$     Y$$$P
// *  `88bo,__,o,    oP"``"Yo,  _88o,,od8P   oP"``"Yo,
// *    "YUMMMMMP",m"       "Mm,""YUMMMP" ,m"       "Mm,
// *
// *   src->core->HLE->JVS->JVS.h
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
// *  (c) 2019 Luke Usher <luke.usher@outlook.com>
// *
// *  All rights reserved
// *
// ******************************************************************
#ifndef JVS_H
#define JVS_H

DWORD WINAPI EMUPATCH(JvsBACKUP_Read)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsBACKUP_Write)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsEEPROM_Read)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsEEPROM_Write)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsFirmwareDownload)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsFirmwareUpload)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsNodeReceivePacket)
(
	DWORD a1,
	DWORD a2,
	DWORD a3
);

DWORD WINAPI EMUPATCH(JvsNodeSendPacket)
(
	DWORD a1,
	DWORD a2,
	DWORD a3
);

DWORD WINAPI EMUPATCH(JvsRTC_Read)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsScFirmwareDownload)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsScFirmwareUpload)
(
	DWORD a1,
	DWORD a2,
	DWORD a3,
	DWORD a4
);

DWORD WINAPI EMUPATCH(JvsScReceiveMidi)
(
	DWORD a1,
	DWORD a2,
	DWORD a3
);

DWORD WINAPI EMUPATCH(JvsScSendMidi)
(
	DWORD a1,
	DWORD a2,
	DWORD a3
);

DWORD WINAPI EMUPATCH(JvsScReceiveRs323c)
(
	DWORD a1,
	DWORD a2,
	DWORD a3
);


DWORD WINAPI EMUPATCH(JvsScSendRs323c)
(
	DWORD a1,
	DWORD a2,
	DWORD a3
);

#endif