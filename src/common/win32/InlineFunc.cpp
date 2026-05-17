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
// *  (c) 2018 RadWolfie
// *
// *  All rights reserved
// *
// ******************************************************************

#if defined(_WIN32) || defined(WIN32)

#include <windows.h>
#include <optional>
#include "common/util/cliConfig.hpp"
#include "Util.h"

// Source: https://stackoverflow.com/questions/8046097/how-to-check-if-a-process-has-the-administrative-rights
bool CxbxrIsElevated() {
	bool fRet = false;
	HANDLE hToken = NULL;
	if (OpenProcessToken(GetCurrentProcess(), TOKEN_QUERY, &hToken)) {
		TOKEN_ELEVATION Elevation;
		DWORD cbSize = sizeof(TOKEN_ELEVATION);
		if (GetTokenInformation(hToken, TokenElevation, &Elevation, sizeof(Elevation), &cbSize)) {
			fRet = Elevation.TokenIsElevated > 0;
		}
	}
	if (hToken) {
		CloseHandle(hToken);
	}
	return fRet;
}

std::optional<std::string> CxbxrExec(bool useDebugger, void** hProcess, bool requestHandleProcess, bool isReboot) {

	STARTUPINFO startupInfo = { 0 };
	PROCESS_INFORMATION processInfo = { 0 };
	std::string szProcArgsBuffer;
	if (!cli_config::GenCMD(szProcArgsBuffer)) {
		return std::make_optional<std::string>("Failed to retrieve the command line arguments to launch the new emulation process");
	}

	// TODO: Set a configuration variable for this. For now it will be within the same folder as Cxbx.exe
	if (useDebugger) {
		szProcArgsBuffer = "cxbxr-debugger.exe " + szProcArgsBuffer;
	}

	/* NOTE: CreateProcess's 2nd parameter (lpCommandLine) is char*, not const char*. Plus it has ability to change the input buffer data.
		Source: https://msdn.microsoft.com/en-us/library/ms682425.aspx

		Using ShellExecute has proper implement. Unfortunately, we need created process handle for Debugger monitor.
		Plus ShellExecute is high level whilst CreateProcess is low level. We want to use official low level functions as possible to reduce
		cpu load cycles to get the task done.

		Without the DETACHED_PROCESS flag, the default behavior would be for the new process to inherit the console of the parent process,
		which is wrong since we want the emulation process to have its own console allocated with AllocConsole instead.
	*/
	if (CreateProcess(nullptr, const_cast<LPSTR>(szProcArgsBuffer.c_str()), nullptr, nullptr, false, DETACHED_PROCESS, nullptr, nullptr, &startupInfo, &processInfo) == 0) {
		return std::make_optional<std::string>("Failed to create the new emulation process. CreateProcess failed because: " + WinError2Str());
	}

	// Place the child process in a Job Object so that it is automatically
	// terminated when cxbx.exe exits (for any reason: graceful close, crash,
	// killed via Task Manager).  This is necessary because the render window
	// uses WS_POPUP (owned) instead of WS_CHILD, so Windows no longer
	// auto-destroys it when the parent window/process goes away.
	// During reboot the child inherits the GUI's job via normal process
	// inheritance, so we must not create a second job here (its
	// KILL_ON_JOB_CLOSE would kill the child when this process exits).
	if (!isReboot) {
		static HANDLE s_hJob = NULL;
		if (!s_hJob) {
			s_hJob = CreateJobObject(NULL, NULL);
			if (s_hJob) {
				JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = {};
				jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
				SetInformationJobObject(s_hJob, JobObjectExtendedLimitInformation, &jeli, sizeof(jeli));
			}
		}
		if (s_hJob) {
			AssignProcessToJobObject(s_hJob, processInfo.hProcess);
		}
	}

	// Allow the child process to call SetForegroundWindow so it can claim
	// foreground status after creating its render window (WS_POPUP owned by us).
	AllowSetForegroundWindow(processInfo.dwProcessId);

	CloseHandle(processInfo.hThread);

	if (requestHandleProcess) {
		*hProcess = processInfo.hProcess;
	}
	else {
		CloseHandle(processInfo.hProcess);
	}

	return std::nullopt;
}

#endif
