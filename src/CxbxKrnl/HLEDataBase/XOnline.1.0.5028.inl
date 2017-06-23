// ******************************************************************
// *
// *    .,-:::::    .,::      .::::::::.    .,::      .:
// *  ,;;;'````'    `;;;,  .,;;  ;;;'';;'   `;;;,  .,;;
// *  [[[             '[[,,[['   [[[__[[\.    '[[,,[['
// *  $$$              Y$$$P     $$""""Y$$     Y$$$P
// *  `88bo,__,o,    oP"``"Yo,  _88o,,od8P   oP"``"Yo,
// *    "YUMMMMMP",m"       "Mm,""YUMMMP" ,m"       "Mm,
// *
// *   Cxbx->Win32->CxbxKrnl->XOnline.1.0.5028.inl
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
// *  (c) 2002-2003 Aaron Robinson <caustik@caustik.com>
// *
// *  All rights reserved
// *
// ******************************************************************

// ******************************************************************
// * XOnline_5028
// ******************************************************************
OOVPATable XOnline_5028[] ={

	REGISTER_OOVPAS(XNetStartup, 4361), // PATCH 
	REGISTER_OOVPAS(WSAStartup, 4361), // PATCH 
	REGISTER_OOVPAS(XnInit, 4627), // XREF
	REGISTER_OOVPAS(socket, 4361), // PATCH 
	REGISTER_OOVPAS(bind, 4361), // PATCH 
	REGISTER_OOVPAS(listen, 4361), // PATCH 
	REGISTER_OOVPAS(ioctlsocket, 4361), // PATCH 
	REGISTER_OOVPAS(XNetGetEthernetLinkStatus, 4627), // PATCH 
	REGISTER_OOVPAS(XoUpdateLaunchNewImageInternal, 4627), // XREF
};

// ******************************************************************
// * XOnline_5028_SIZE
// ******************************************************************
uint32 XOnline_5028_SIZE = sizeof(XOnline_5028);
