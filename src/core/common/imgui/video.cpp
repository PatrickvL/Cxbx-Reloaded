// Copyright 2021 Cxbx-Reloaded Project
// Licensed under GPLv2+
// Refer to the COPYING file included.

#define LOG_PREFIX CXBXR_MODULE::GUI

#include <thread>

#include "video.hpp"
#include "ui.hpp"

#include "EmuShared.h"

#include "core/kernel/init/CxbxKrnl.h"

const ImColor ImGuiVideo::m_laser_col[4] = {
		ImColor(ImVec4(1.0f, 0.0f, 0.0f, 1.0f)), // player1: red
		ImColor(ImVec4(0.0f, 1.0f, 0.0f, 1.0f)), // player2: green
		ImColor(ImVec4(0.0f, 0.0f, 1.0f, 1.0f)), // player3: blue
		ImColor(ImVec4(1.0f, 1.0f, 0.0f, 1.0f))  // player4: yellow
};

bool ImGuiVideo::Initialize()
{
	g_EmuShared->GetImGuiVideoWindows(&m_windows);
	return true;
}

void ImGuiVideo::Shutdown()
{
	g_EmuShared->SetImGuiVideoWindows(&m_windows);
}

void ImGuiVideo::DrawMenu()
{
//	if (ImGui::BeginMenu("Video")) {
//		ImGui::EndMenu();
//	}
}

void ImGuiVideo::DrawWidgets(bool is_focus, ImGuiWindowFlags input_handler)
{
	// Render the lightgun laser
	for (int port = PORT_1; port < XBOX_NUM_PORTS; ++port) {
		if (g_devs[port].type == XBOX_INPUT_DEVICE::LIGHTGUN && g_devs[port].info.ligthgun.laser) {
			ImGui::Begin("Laser", nullptr, ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoInputs | ImGuiWindowFlags_NoDecoration);
			ImGui::GetForegroundDrawList()->AddCircleFilled(g_InputDeviceManager.CalcLaserPos(port), 5, m_laser_col[port], 0);
			ImGui::End();
		}
	}
}
