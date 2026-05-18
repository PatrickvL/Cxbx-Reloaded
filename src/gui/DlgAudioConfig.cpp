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

#include "DlgAudioConfig.h"
#include "resource/ResCxbx.h"

/*! windows dialog procedure */
static INT_PTR CALLBACK DlgAudioConfigProc(HWND hWndDlg, UINT uMsg, WPARAM wParam, LPARAM lParam);

/*! audio configuration */
static Settings::s_audio g_XBAudio;
/*! changes flag */
static BOOL g_bHasChanges = FALSE;

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
            SendMessage(GetDlgItem(hWndDlg, IDC_AC_MUTE_WHEN_UNFOCUS), BM_SETCHECK, (WPARAM)(g_XBAudio.mute_on_unfocus == true), 0);
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

                    EndDialog(hWndDlg, wParam);
                }
                break;
            }
        }
        break;
    }
    return FALSE;
}
