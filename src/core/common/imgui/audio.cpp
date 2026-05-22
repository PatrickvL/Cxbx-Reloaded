// Copyright 2021 Cxbx-Reloaded Project
// Licensed under GPLv2+
// Refer to the COPYING file included.

#define LOG_PREFIX CXBXR_MODULE::GUI

#include <thread>

#include "audio.hpp"

#include "EmuShared.h"

bool ImGuiAudio::Initialize()
{
	g_EmuShared->GetImGuiAudioWindows(&m_windows);
	return true;
}

void ImGuiAudio::Shutdown()
{
	g_EmuShared->SetImGuiAudioWindows(&m_windows);
}

void ImGuiAudio::DrawMenu()
{
	if (ImGui::BeginMenu("Audio")) {
		ImGui::MenuItem("Debug General Cache Stats", NULL, &m_windows.cache_stats_general);
		ImGui::MenuItem("SoundBuffer Visualization", NULL, &m_windows.cache_visualization);
		ImGui::EndMenu();
	}
}

void ImGuiAudio::DrawWidgets(bool is_focus, ImGuiWindowFlags input_handler)
{
	(void)is_focus;
	(void)input_handler;
}
