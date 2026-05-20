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
// *  (c) 2017 RadWolfie
// *
// *  All rights reserved
// *
// ******************************************************************

#define LOG_PREFIX CXBXR_MODULE::GUI

#include "common\Settings.hpp" // for g_Settings
#include "common/Logging.h"

#include <AL/alc.h>

#include <cstring>
#include <string>
#include <vector>

#include "DlgAudioConfig.h"
#include "resource/ResCxbx.h"

/*! windows dialog procedure */
static INT_PTR CALLBACK DlgAudioConfigProc(HWND hWndDlg, UINT uMsg, WPARAM wParam, LPARAM lParam);

/*! audio configuration */
static Settings::s_audio g_XBAudio;
/*! changes flag */
static BOOL g_bHasChanges = FALSE;

namespace {

constexpr const char* AUDIO_DEFAULT_DEVICE_LABEL = "<Default OpenAL device>";
constexpr size_t AUDIO_DEVICE_NAME_LIMIT = 4096;
constexpr size_t AUDIO_DEVICE_COUNT_LIMIT = 256;

std::vector<std::string> GetAvailableOpenALDevices()
{
    const ALCchar* deviceList = nullptr;
    if (alcIsExtensionPresent(nullptr, "ALC_ENUMERATE_ALL_EXT") == ALC_TRUE) {
        deviceList = alcGetString(nullptr, ALC_ALL_DEVICES_SPECIFIER);
    } else if (alcIsExtensionPresent(nullptr, "ALC_ENUMERATION_EXT") == ALC_TRUE) {
        deviceList = alcGetString(nullptr, ALC_DEVICE_SPECIFIER);
    }

    std::vector<std::string> devices;
    if (deviceList == nullptr) {
        return devices;
    }

    for (const ALCchar* current = deviceList; *current != '\0' && devices.size() < AUDIO_DEVICE_COUNT_LIMIT;) {
        const size_t nameLength = strnlen(current, AUDIO_DEVICE_NAME_LIMIT);
        if (nameLength == AUDIO_DEVICE_NAME_LIMIT) {
            EmuLog(LOG_LEVEL::WARNING, "OpenAL device enumeration stopped after an unexpectedly long device name");
            break;
        }
        devices.emplace_back(current, nameLength);
        current += nameLength + 1;
    }
    return devices;
}

void PopulateAudioDeviceList(HWND hWndDlg)
{
    HWND deviceCombo = GetDlgItem(hWndDlg, IDC_AC_AUDIO_ADAPTER);
    const std::string configuredDevice = g_Settings->GetAudioOutputDevice();
    const std::vector<std::string> devices = GetAvailableOpenALDevices();

    SendMessage(deviceCombo, CB_RESETCONTENT, 0, 0);
    SendMessage(deviceCombo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(AUDIO_DEFAULT_DEVICE_LABEL));

    int selectedIndex = 0;
    bool configuredDeviceFound = configuredDevice.empty();
    for (const std::string& device : devices) {
        const int index = static_cast<int>(SendMessage(deviceCombo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(device.c_str())));
        if (!configuredDevice.empty() && device == configuredDevice) {
            selectedIndex = index;
            configuredDeviceFound = true;
        }
    }

    if (!configuredDevice.empty() && !configuredDeviceFound) {
        selectedIndex = static_cast<int>(SendMessage(deviceCombo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(configuredDevice.c_str())));
    }

    SendMessage(deviceCombo, CB_SETCURSEL, selectedIndex, 0);
}

std::string GetSelectedAudioDevice(HWND hWndDlg)
{
    HWND deviceCombo = GetDlgItem(hWndDlg, IDC_AC_AUDIO_ADAPTER);
    const LRESULT selectedIndex = SendMessage(deviceCombo, CB_GETCURSEL, 0, 0);
    if (selectedIndex <= 0) {
        return "";
    }

    const LRESULT textLength = SendMessage(deviceCombo, CB_GETLBTEXTLEN, selectedIndex, 0);
    if (textLength < 0) {
        return "";
    }

    std::vector<char> selectedDevice(static_cast<size_t>(textLength) + 1, '\0');
    SendMessage(deviceCombo, CB_GETLBTEXT, selectedIndex, reinterpret_cast<LPARAM>(selectedDevice.data()));
    return std::string(selectedDevice.data());
}

} // namespace

void ShowAudioConfig(HWND hwnd)
{
    /*! reset changes flag */
    g_bHasChanges = FALSE;

    /*! retrieve audio configuration */
    g_XBAudio = g_Settings->m_audio;

    /*! show dialog box */
    DialogBox(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDD_AUDIO_CFG), hwnd, DlgAudioConfigProc);
}

INT_PTR CALLBACK DlgAudioConfigProc(HWND hWndDlg, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch(uMsg)
    {
        case WM_INITDIALOG:
        {
            /*! set window icon */
            SetClassLong(hWndDlg, GCL_HICON, (LONG)LoadIcon(GetModuleHandle(nullptr), MAKEINTRESOURCE(IDI_CXBX)));

            /*! check appropriate options */
            SendMessage(GetDlgItem(hWndDlg, IDC_AC_PCM), BM_SETCHECK, (WPARAM)g_XBAudio.codec_pcm, 0);
            SendMessage(GetDlgItem(hWndDlg, IDC_AC_XADPCM), BM_SETCHECK, (WPARAM)g_XBAudio.codec_xadpcm, 0);
            SendMessage(GetDlgItem(hWndDlg, IDC_AC_UNKNOWN_CODEC), BM_SETCHECK, (WPARAM)g_XBAudio.codec_unknown, 0);
            SendMessage(GetDlgItem(hWndDlg, IDC_AC_MUTE_WHEN_UNFOCUS), BM_SETCHECK, (WPARAM)g_XBAudio.mute_on_unfocus, 0);
            PopulateAudioDeviceList(hWndDlg);
        }
        break;

        case WM_CLOSE:
        {
            /*! if changes have been made, check if the user wants to save them */
            if(g_bHasChanges)
            {
                PopupReturn ret = PopupQuestion(hWndDlg, "Do you wish to apply your changes?");

                switch(ret)
                {
                    case PopupReturn::Yes:
                        PostMessage(hWndDlg, WM_COMMAND, IDC_AC_ACCEPT, 0);
                        break;
                    case PopupReturn::No:
                        PostMessage(hWndDlg, WM_COMMAND, IDC_AC_CANCEL, 0);
                        break;
                }
                break;
            }

            PostMessage(hWndDlg, WM_COMMAND, IDC_AC_CANCEL, 0);
        }
        break;

        case WM_COMMAND:
        {
            switch(LOWORD(wParam))
            {
                case IDC_AC_CANCEL:
                    EndDialog(hWndDlg, wParam);
                    break;

                case IDC_AC_ACCEPT:
                {
                    /*! save PCM/XADPCM/UnknownCodec options */
                    LRESULT lRet = SendMessage(GetDlgItem(hWndDlg, IDC_AC_PCM), BM_GETCHECK, 0, 0);
                    g_XBAudio.codec_pcm = (lRet == BST_CHECKED);

                    lRet = SendMessage(GetDlgItem(hWndDlg, IDC_AC_XADPCM), BM_GETCHECK, 0, 0);
                    g_XBAudio.codec_xadpcm = (lRet == BST_CHECKED);

                    lRet = SendMessage(GetDlgItem(hWndDlg, IDC_AC_UNKNOWN_CODEC), BM_GETCHECK, 0, 0);
                    g_XBAudio.codec_unknown = (lRet == BST_CHECKED);

                    lRet = SendMessage(GetDlgItem(hWndDlg, IDC_AC_MUTE_WHEN_UNFOCUS), BM_GETCHECK, 0, 0);
                    g_XBAudio.mute_on_unfocus = (lRet == BST_CHECKED);

                    /*! save audio configuration */
                    g_Settings->m_audio = g_XBAudio;
                    g_Settings->SetAudioOutputDevice(GetSelectedAudioDevice(hWndDlg));

                    EndDialog(hWndDlg, wParam);
                }
                break;
            }
        }
        break;
    }
    return FALSE;
}
