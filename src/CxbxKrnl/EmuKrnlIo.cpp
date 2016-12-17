// ******************************************************************
// *
// *    .,-:::::    .,::      .::::::::.    .,::      .:
// *  ,;;;'````'    `;;;,  .,;;  ;;;'';;'   `;;;,  .,;;
// *  [[[             '[[,,[['   [[[__[[\.    '[[,,[['
// *  $$$              Y$$$P     $$""""Y$$     Y$$$P
// *  `88bo,__,o,    oP"``"Yo,  _88o,,od8P   oP"``"Yo,
// *    "YUMMMMMP",m"       "Mm,""YUMMMP" ,m"       "Mm,
// *
// *   Cxbx->Win32->CxbxKrnl->EmuKrnlIo.cpp
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
// *  (c) 2016 Patrick van Logchem <pvanlogchem@gmail.com>
// *
// *  All rights reserved
// *
// ******************************************************************
#define _CXBXKRNL_INTERNAL
#define _XBOXKRNL_DEFEXTRN_

// prevent name collisions
namespace xboxkrnl
{
#include <xboxkrnl/xboxkrnl.h> // For IoCompletionObjectType, etc.
};

#include "Logging.h" // For LOG_FUNC()
#include "EmuKrnlLogging.h"
#include "CxbxKrnl.h" // For CxbxKrnlCleanup
#include "Emu.h" // For EmuWarning()
#include "EmuFile.h" // For CxbxCreateSymbolicLink(), etc.

// ******************************************************************
// * 0x0040 - IoCompletionObjectType
// ******************************************************************
XBSYSAPI EXPORTNUM(64) xboxkrnl::OBJECT_TYPE xboxkrnl::IoCompletionObjectType =
{
	/*
	ExAllocatePoolWithTag,
	ExFreePool,
	NULL,
	IopDeleteIoCompletion,
	NULL,
	*/
	NULL, // &ObpDefaultObject,
	'pmoC' // = first four characters of "Completion" in reverse
};

// ******************************************************************
// * 0x0042 - IoCreateFile()
// ******************************************************************
XBSYSAPI EXPORTNUM(66) xboxkrnl::NTSTATUS NTAPI xboxkrnl::IoCreateFile
(
	OUT PHANDLE             FileHandle,
	IN  ACCESS_MASK         DesiredAccess,
	IN  POBJECT_ATTRIBUTES  ObjectAttributes,
	OUT PIO_STATUS_BLOCK    IoStatusBlock,
	IN  PLARGE_INTEGER      AllocationSize OPTIONAL,
	IN  ULONG               FileAttributes,
	IN  ULONG               ShareAccess,
	IN  ULONG               Disposition,
	IN  ULONG               CreateOptions,
	IN  ULONG               Options
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG_OUT(FileHandle)
		LOG_FUNC_ARG(DesiredAccess)
		LOG_FUNC_ARG(ObjectAttributes)
		LOG_FUNC_ARG_OUT(IoStatusBlock)
		LOG_FUNC_ARG(AllocationSize)
		LOG_FUNC_ARG(FileAttributes)
		LOG_FUNC_ARG(ShareAccess)
		LOG_FUNC_ARG(Disposition)
		LOG_FUNC_ARG(CreateOptions)
		LOG_FUNC_ARG(Options)
		LOG_FUNC_END;

	NativeObjectAttributes nativeObjectAttributes;

	NTSTATUS ret = CxbxObjectAttributesToNT(ObjectAttributes, nativeObjectAttributes, "IoCreateFile"); /*var*/

	// redirect to NtCreateFile
	ret = NtDll::NtCreateFile(
		FileHandle, 
		DesiredAccess | GENERIC_READ, 
		nativeObjectAttributes.NtObjAttrPtr, 
		NtDll::PIO_STATUS_BLOCK(IoStatusBlock), 
		NtDll::PLARGE_INTEGER(AllocationSize), 
		FileAttributes, 
		ShareAccess, 
		Disposition, 
		CreateOptions, 
		NULL, 
		0);

	if (FAILED(ret))
	{
		DbgPrintf("EmuKrnl (0x%X): IoCreateFile Failed! (0x%.08X)\n", GetCurrentThreadId(), ret);

		// TODO: refactor this ugly mess

		// It seems that some titles will attempt to write files to the utility drive
		// without creating the necessary parent folders.  I'm not sure if this is
		// intended behavior, but whenever a file creation attempt fails because its
		// parent folder does not exist (and if the file was supposed to be created),
		// the emulator will create the necessary parent folders first and try the
		// operation again.  It would be nice to confirm this behavior on a real Xbox.

		// If there was an attempt to create a file, check if the parent folders exist
		if ((0xC000003A == ret) // TODO: constantify to STATUS_OBJECT_PATH_NOT_FOUND
			&& (FILE_CREATE == Disposition || FILE_OPEN_IF == Disposition || FILE_OVERWRITE_IF == Disposition)
			&& (CreateOptions & FILE_NON_DIRECTORY_FILE))
		{
			// Convert the path to a wstring
			std::wstring path(nativeObjectAttributes.NtObjAttrPtr->ObjectName->Buffer, nativeObjectAttributes.NtObjAttrPtr->ObjectName->Length / sizeof(NtDll::WCHAR));
			
			// Check if it has parent folders
			size_t indexOfSlash = path.rfind('\\');
			if (indexOfSlash >= 0)
			{
				DbgPrintf("EmuKrnl (0x%X): IoCreateFile - attempting to create parent folders...\n", GetCurrentThreadId());
				
				// It does; create them as necessary
				ret = CxbxCreateParentFolders(nativeObjectAttributes.NtObjAttrPtr->RootDirectory, path.substr(0, indexOfSlash));
				if (FAILED(ret))
				{
					DbgPrintf("EmuKrnl (0x%X): IoCreateFile - could not create parent folders! (0x%.08X)\n", GetCurrentThreadId(), ret);

					// Make sure we return the original error code
					ret = 0xC000003A; // TODO: constantify to STATUS_OBJECT_PATH_NOT_FOUND
				}
				else
				{
					DbgPrintf("EmuKrnl (0x%X): IoCreateFile - reattempting operation...\n", GetCurrentThreadId());

					// Try creating the file again
					ret = NtDll::NtCreateFile(
						FileHandle,
						DesiredAccess | GENERIC_READ,
						nativeObjectAttributes.NtObjAttrPtr,
						NtDll::PIO_STATUS_BLOCK(IoStatusBlock),
						NtDll::PLARGE_INTEGER(AllocationSize),
						FileAttributes,
						ShareAccess,
						Disposition,
						CreateOptions,
						NULL,
						0);

					if (FAILED(ret))
					{
						DbgPrintf("EmuKrnl (0x%X): IoCreateFile Failed! (0x%.08X)\n", GetCurrentThreadId(), ret);
					}
					else
					{
						DbgPrintf("EmuKrnl (0x%X): IoCreateFile = 0x%.08X\n", GetCurrentThreadId(), *FileHandle);
					}
				}
			}
			else
			{
				// No parent directories; log the error
				DbgPrintf("EmuKrnl (0x%X): IoCreateFile Failed! (0x%.08X)\n", GetCurrentThreadId(), ret);
			}
		}
	}
	else
	{
		DbgPrintf("EmuKrnl (0x%X): IoCreateFile = 0x%.08X\n", GetCurrentThreadId(), *FileHandle);
	}

	RETURN(ret);
}

// ******************************************************************
// * 0x0043 - IoCreateSymbolicLink()
// ******************************************************************
XBSYSAPI EXPORTNUM(67) xboxkrnl::NTSTATUS NTAPI xboxkrnl::IoCreateSymbolicLink
(
	IN PSTRING SymbolicLinkName,
	IN PSTRING DeviceName
)
{
	LOG_FUNC_BEGIN
		LOG_FUNC_ARG(SymbolicLinkName)
		LOG_FUNC_ARG(DeviceName)
		LOG_FUNC_END;

	NTSTATUS ret = CxbxCreateSymbolicLink(std::string(SymbolicLinkName->Buffer, SymbolicLinkName->Length), std::string(DeviceName->Buffer, DeviceName->Length));

	RETURN(ret);
}

// ******************************************************************
// * 0x0045 - IoDeleteSymbolicLink()
// ******************************************************************
XBSYSAPI EXPORTNUM(69) xboxkrnl::NTSTATUS NTAPI xboxkrnl::IoDeleteSymbolicLink
(
	IN PSTRING SymbolicLinkName
)
{
	LOG_FUNC_ONE_ARG(SymbolicLinkName);

	EmuNtSymbolicLinkObject* symbolicLink = FindNtSymbolicLinkObjectByName(std::string(SymbolicLinkName->Buffer, SymbolicLinkName->Length));

	NTSTATUS ret = STATUS_OBJECT_NAME_NOT_FOUND;

	if ((symbolicLink != NULL))
		ret = symbolicLink->NtClose();

	RETURN(ret);
}

// ******************************************************************
// * 0x0046 - IoDeviceObjectType
// ******************************************************************
XBSYSAPI EXPORTNUM(70) xboxkrnl::OBJECT_TYPE xboxkrnl::IoDeviceObjectType =
{
	/*
	ExAllocatePoolWithTag,
	ExFreePool,
	NULL,
	NULL,
	IoParseDevice,
	*/
	NULL, // &ObpDefaultObject,
	'iveD' // = first four characters of "Device" in reverse
};

// ******************************************************************
// * 0x0047 - IoFileObjectType
// ******************************************************************
XBSYSAPI EXPORTNUM(71) xboxkrnl::OBJECT_TYPE xboxkrnl::IoFileObjectType =
{
	/*
	ExAllocatePoolWithTag,
	ExFreePool,
	IopCloseFile,
	IopDeleteFile,
	IopParseFile,
	*/
	NULL, // (PVOID)FIELD_OFFSET(FILE_OBJECT, Event.Header),
	'eliF' // = "File" in reverse
};

// ******************************************************************
// * 0x005A - IoDismountVolume()
// ******************************************************************
XBSYSAPI EXPORTNUM(90) xboxkrnl::NTSTATUS NTAPI xboxkrnl::IoDismountVolume
(
	IN PDEVICE_OBJECT DeviceObject
)
{
	LOG_FUNC_ONE_ARG(DeviceObject);

	NTSTATUS ret = STATUS_SUCCESS;

	LOG_UNIMPLEMENTED();

	RETURN(ret);
}

// ******************************************************************
// * 0x005B - IoDismountVolumeByName()
// ******************************************************************
XBSYSAPI EXPORTNUM(91) xboxkrnl::NTSTATUS NTAPI xboxkrnl::IoDismountVolumeByName
(
	IN PSTRING VolumeName
)
{
	LOG_FUNC_ONE_ARG(VolumeName);

	NTSTATUS ret = STATUS_SUCCESS;

	// TODO: Anything?
	LOG_UNIMPLEMENTED();

	RETURN(ret);
}

