cmake_minimum_required (VERSION 3.12)

project(Cxbx-Reloaded)

# Allow devs to disable regeneration option.
# Suppress extra stuff from generated solution
#set(CMAKE_SUPPRESS_REGENERATION true)

# Output all binary files into one folder
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
 set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${Cxbx-Reloaded_BINARY_DIR}/bin")
endif()

set(SUBHOOK_TESTS OFF)
set(SUBHOOK_INSTALL OFF)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/import/subhook" EXCLUDE_FROM_ALL)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/import/XbSymbolDatabase")

# Not require since only include a header file
#add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/import/simpleini")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/import/SDL2" EXCLUDE_FROM_ALL)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/import/mio" EXCLUDE_FROM_ALL)

# Optional: detect Vulkan SDK for future backend work
find_package(Vulkan QUIET)
if(Vulkan_FOUND)
    message(STATUS "Vulkan SDK found: ${Vulkan_LIBRARY}")
endif()

# Cxbx-Reloaded projects
set(CXBXR_ROOT_DIR ${CMAKE_CURRENT_LIST_DIR})

add_custom_target(misc-batch
  ${CMAKE_COMMAND} -DTargetRunTimeDir=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>
  -P ${CXBXR_ROOT_DIR}/projects/misc/batch.cmake
  WORKING_DIRECTORY ${CXBXR_ROOT_DIR}
)

# Custom CMake projects since some import libraries doesn't have it.

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/projects/libtom")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/projects/imgui")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/projects/libusb")

set(nv2a_vsh_cpu_UNIT_TEST OFF)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/import/nv2a_vsh_cpu" EXCLUDE_FROM_ALL)

# Split the files into group for which project is likely
# going to be used for both header and source files.
# Then move only specific project files into their
# cmakelist file in finalize stage.
# Any of files being shared between two projects need to
# leave in parent cmakelist (here) for better sync across
# child projects.

# Common (GUI and Emulator)
file (GLOB CXBXR_HEADER_COMMON
 "${CXBXR_ROOT_DIR}/src/common/AddressRanges.h"
 "${CXBXR_ROOT_DIR}/src/common/crypto/EmuDes.h"
 "${CXBXR_ROOT_DIR}/src/common/crypto/EmuRsa.h"
 "${CXBXR_ROOT_DIR}/src/common/crypto/EmuSha.h"
 "${CXBXR_ROOT_DIR}/src/common/crypto/LibRc4.h"
 "${CXBXR_ROOT_DIR}/src/common/CxbxDebugger.h"
 "${CXBXR_ROOT_DIR}/src/common/cxbxr.hpp"
 "${CXBXR_ROOT_DIR}/src/common/EmuEEPROM.h"
 "${CXBXR_ROOT_DIR}/src/common/Error.h"
 "${CXBXR_ROOT_DIR}/src/common/FilePaths.hpp"
 "${CXBXR_ROOT_DIR}/src/common/input/DInputKeyboardCodes.h"
 "${CXBXR_ROOT_DIR}/src/common/input/DInputKeyboardMouse.h"
 "${CXBXR_ROOT_DIR}/src/common/input/layout_xbox_device.h"
 "${CXBXR_ROOT_DIR}/src/common/input/LibusbDevice.h"
 "${CXBXR_ROOT_DIR}/src/common/input/InputDevice.h"
 "${CXBXR_ROOT_DIR}/src/common/input/InputManager.h"
 "${CXBXR_ROOT_DIR}/src/common/input/SdlJoystick.h"
 "${CXBXR_ROOT_DIR}/src/common/input/XInputPad.h"
 "${CXBXR_ROOT_DIR}/src/common/input/RawDevice.h"
 "${CXBXR_ROOT_DIR}/src/common/IPCHybrid.hpp"
 "${CXBXR_ROOT_DIR}/src/common/Logging.h"
 "${CXBXR_ROOT_DIR}/src/common/Settings.hpp"
 "${CXBXR_ROOT_DIR}/src/common/util/cliConfig.hpp"
 "${CXBXR_ROOT_DIR}/src/common/util/cliConverter.hpp"
 "${CXBXR_ROOT_DIR}/src/common/util/CPUID.h"
 "${CXBXR_ROOT_DIR}/src/common/util/CxbxUtil.h"
 "${CXBXR_ROOT_DIR}/src/common/util/std_extend.hpp"
 "${CXBXR_ROOT_DIR}/src/common/util/strConverter.hpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/AlignPosfix1.h"
 "${CXBXR_ROOT_DIR}/src/common/win32/AlignPrefix1.h"
 "${CXBXR_ROOT_DIR}/src/common/win32/EmuShared.h"
 "${CXBXR_ROOT_DIR}/src/common/win32/Mutex.h"
 "${CXBXR_ROOT_DIR}/src/common/win32/Util.h"
 "${CXBXR_ROOT_DIR}/src/common/win32/WineEnv.h"
 "${CXBXR_ROOT_DIR}/src/common/xbdm/CxbxXbdm.h"
 "${CXBXR_ROOT_DIR}/src/common/xbe/Xbe.h"
 "${CXBXR_ROOT_DIR}/src/common/xbe/XbePrinter.h"
 "${CXBXR_ROOT_DIR}/src/common/xbox/Logging.hpp"
 "${CXBXR_ROOT_DIR}/src/common/xbox/Types.hpp"
 "${CXBXR_ROOT_DIR}/src/common/xbox_types.h"
 "${CXBXR_ROOT_DIR}/src/common/xdvdfs-tools/buffered_io.h"
 "${CXBXR_ROOT_DIR}/src/common/xdvdfs-tools/xdvdfs.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbConvert.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbD3D8Types.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbD3D8Types_D3D11.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbD3D8Types_Enums.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbD3D8Types_Resources.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/XAPI/Xapi.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlLogging.h"
 "${CXBXR_ROOT_DIR}/src/Cxbx.h"
 "${CXBXR_ROOT_DIR}/src/CxbxVersion.h"
 "${CXBXR_ROOT_DIR}/src/version.h"
)

# GUI v1
file (GLOB CXBXR_HEADER_GUIv1
 "${CXBXR_ROOT_DIR}/src/common/input/Button.h"
 "${CXBXR_ROOT_DIR}/src/common/input/EmuDevice.h"
 "${CXBXR_ROOT_DIR}/src/common/input/InputWindow.h"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgDukeControllerConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgLibusbControllerConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgLightgunConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgSBControllerConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgAbout.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgAudioConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgInputConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgEepromConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgLoggingConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgNetworkConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/DlgVideoConfig.h"
 "${CXBXR_ROOT_DIR}/src/gui/Wnd.h"
 "${CXBXR_ROOT_DIR}/src/gui/WndMain.h"
)

# Emulator (module)
set(IMGUI_GFX_HEADER "${CXBXR_ROOT_DIR}/import/imgui/backends/imgui_impl_dx11.h")
set(IMGUI_GFX_SOURCE "${CXBXR_ROOT_DIR}/import/imgui/backends/imgui_impl_dx11.cpp")
file (GLOB CXBXR_HEADER_EMU_IMPORT
 "${IMGUI_GFX_HEADER}"
 "${CXBXR_ROOT_DIR}/import/imgui/backends/imgui_impl_win32.h"
)
file (GLOB CXBXR_HEADER_EMU
 "${CXBXR_ROOT_DIR}/src/common/audio/converter.hpp"
 "${CXBXR_ROOT_DIR}/src/common/audio/XADPCM.h"
 "${CXBXR_ROOT_DIR}/src/common/Timer.h"
 "${CXBXR_ROOT_DIR}/src/common/win32/Threads.h"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/audio.hpp"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/ui.hpp"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/settings.h"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/video.hpp"
 "${CXBXR_ROOT_DIR}/src/core/common/video/RenderBase.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxPixelShaderInput.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxPixelShaderFunctions.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxScreenspaceTransform.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxVertexShaderCommon.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxVertexFetch.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxVertexFetchLayout.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxVertexOutputFooter.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/RenderGlobals.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/NV2A_PGRAPH_Helpers.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/EmuD3D8_common.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxNV2AMathHelpers.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxNV2APixelShaderConstants.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxNV2AVshOps.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxPGRAPHRegs.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxPSAuxFromPGRAPH.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxRegisterCombinerInterpreter.hlsl"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxRegisterCombinerInterpreterState.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxVertexShaderInterpreter.hlsl"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxVertexShaderInterpreterState.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxFixedFunctionVertexShader.hlsl"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/CxbxFixedFunctionVertexShaderState.hlsli"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_Internal.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/PatchDraw.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders/Shader.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/WalkIndexBuffer.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/FixedFunctionState.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/ResourceTracker.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbD3D8Logging.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbNV2A.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbPixelShader.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbPushBuffer.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbVertexBuffer.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbVertexShader.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/common/windows/WFXformat.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSound.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundGlobal.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/ApuPlayCursor.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundInline.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundLogging.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DSStream_PacketManager.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/XbDSoundFuncs.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/XbDSoundLogging.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/XbDSoundTypes.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/common/XbInternalStruct.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/Intercept.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/Patches.hpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/XACTENG/XactEng.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/XGRAPHIC/XGraphic.h"
 "${CXBXR_ROOT_DIR}/src/core/hle/XONLINE/XOnline.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/dbg.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/ex.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/hal.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/io.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/kernel.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/mm.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/nt.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/ob.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/ps.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/rtl.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/strings.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/types.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/common/xbox.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlAvModes.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlIo.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlKe.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlKi.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlNt.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlPs.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/init/CxbxKrnl.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/init/KrnlPatches.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/memory-manager/PhysicalMemory.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/memory-manager/PoolManager.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/memory-manager/VMManager.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/Emu.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuCdRom.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuDisk.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuFile.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuFS.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuMediaBoard.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuMu.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuNtDll.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/NativeHandle.h"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/PatchRdtsc.hpp"
 "${CXBXR_ROOT_DIR}/src/devices/ADM1032Device.h"
 "${CXBXR_ROOT_DIR}/src/devices/Chihiro/JvsIo.h"
 "${CXBXR_ROOT_DIR}/src/devices/Chihiro/MediaBoard.h"
 "${CXBXR_ROOT_DIR}/src/devices/EEPROMDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/audio/APUDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/network/NVNetDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/LED.h"
 "${CXBXR_ROOT_DIR}/src/devices/MCPXDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/PCIBus.h"
 "${CXBXR_ROOT_DIR}/src/devices/PCIDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/SMBus.h"
 "${CXBXR_ROOT_DIR}/src/devices/SMCDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/SMDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/usb/Hub.h"
 "${CXBXR_ROOT_DIR}/src/devices/usb/OHCI.h"
 "${CXBXR_ROOT_DIR}/src/devices/usb/UsbCommon.h"
 "${CXBXR_ROOT_DIR}/src/devices/usb/USBDevice.h"
 "${CXBXR_ROOT_DIR}/src/devices/usb/XidGamepad.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/nv2a.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/nv2a_debug.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/nv2a_int.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/nv2a_regs.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/qemu-thread.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/queue.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/swizzle.h"
 "${CXBXR_ROOT_DIR}/src/devices/video/vga.h"
 "${CXBXR_ROOT_DIR}/src/devices/x86/EmuX86.h"
 "${CXBXR_ROOT_DIR}/src/devices/Xbox.h"
)

# filter hlsl files into its own list
# excluding hlsli file(s) (which are only needed as #include dependencies)
set(CXBXR_HEADER_HLSL ${CXBXR_HEADER_EMU})
list(FILTER CXBXR_HEADER_HLSL INCLUDE REGEX ".*\\.hlsl$")

# Common (GUI and Emulator)
file (GLOB CXBXR_SOURCE_COMMON
 "${CXBXR_ROOT_DIR}/src/common/AddressRanges.cpp"
 "${CXBXR_ROOT_DIR}/src/common/crypto/EmuDes.cpp"
 "${CXBXR_ROOT_DIR}/src/common/crypto/EmuRsa.cpp"
 "${CXBXR_ROOT_DIR}/src/common/crypto/EmuSha.cpp"
 "${CXBXR_ROOT_DIR}/src/common/crypto/LibRc4.cpp"
 "${CXBXR_ROOT_DIR}/src/common/CxbxDebugger.cpp"
 "${CXBXR_ROOT_DIR}/src/common/cxbxr.cpp"
 "${CXBXR_ROOT_DIR}/src/common/EmuEEPROM.cpp"
 "${CXBXR_ROOT_DIR}/src/common/Error.cpp"
 "${CXBXR_ROOT_DIR}/src/common/FilePaths.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/DInputKeyboardMouse.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/InputDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/InputManager.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/LibusbDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/SdlJoystick.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/XInputPad.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/RawDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/common/Logging.cpp"
 "${CXBXR_ROOT_DIR}/src/common/Settings.cpp"
 "${CXBXR_ROOT_DIR}/src/common/util/cliConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/common/util/cliConverter.cpp"
 "${CXBXR_ROOT_DIR}/src/common/util/CxbxUtil.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/EmuShared.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/InlineFunc.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/IPCWindows.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/Mutex.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/Util.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/WineEnv.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xbdm/CxbxXbdm.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xbe/Xbe.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xbe/XbePrinter.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xbox/Logging.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xbox/Types.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xdvdfs-tools/buffered_io.cpp"
 "${CXBXR_ROOT_DIR}/src/common/xdvdfs-tools/xdvdfs.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbConvert.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbConvert_ColorFormats.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/XAPI/Xapi.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlLogging.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlXbox.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlXc.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlXe.cpp"
 "${CXBXR_ROOT_DIR}/src/CxbxVersion.cpp"
 "${CXBXR_ROOT_DIR}/src/HighPerformanceGraphicsEnabler.c"
)

# GUI v1
file (GLOB CXBXR_SOURCE_GUIv1
 "${CXBXR_ROOT_DIR}/src/common/input/Button.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/EmuDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/common/input/InputWindow.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgDukeControllerConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgLibusbControllerConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgLightgunConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/input/DlgSBControllerConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgAbout.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgAudioConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgInputConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgEepromConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgLoggingConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgNetworkConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/DlgVideoConfig.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/WinMain.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/Wnd.cpp"
 "${CXBXR_ROOT_DIR}/src/gui/WndMain.cpp"
)

# Emulator (module)
# TODO: Merge CXBXR_KRNL_CPP into CXBXR_SOURCE_EMU once cxbx project is removed.
file (GLOB CXBXR_KRNL_CPP
 "${CXBXR_ROOT_DIR}/src/core/kernel/init/CxbxKrnl.cpp"
)
file (GLOB CXBXR_SOURCE_EMU_IMPORT
 "${IMGUI_GFX_SOURCE}"
 "${CXBXR_ROOT_DIR}/import/imgui/backends/imgui_impl_win32.cpp"
)
file (GLOB CXBXR_SOURCE_EMU
 "${CXBXR_KRNL_CPP}"
 "${CXBXR_ROOT_DIR}/src/common/Timer.cpp"
 "${CXBXR_ROOT_DIR}/src/common/VerifyAddressRanges.cpp"
 "${CXBXR_ROOT_DIR}/src/common/win32/Threads.cpp"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/audio.cpp"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/ui.cpp"
 "${CXBXR_ROOT_DIR}/src/core/common/imgui/video.cpp"
 "${CXBXR_ROOT_DIR}/src/core/common/video/RenderBase.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_State.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_Compute.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_Draw.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_VertexFetch.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_PageTracker.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Shading/VertexShaderCache.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Shading/PixelShaderCache.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Shading/ShaderDiskCache.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Shading/Shader.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/RenderGlobals.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/NV2A_PGRAPH_Helpers.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostImGui.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostRender.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostResource.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostResourceCreate.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostResourceUpload.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostSync.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/HostWindow.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/IndexBufferConvert.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/PatchDraw.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/WalkIndexBuffer.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/FixedFunctionState.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/ResourceTracker.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbD3D8Logging.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/XbPushBuffer.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_VertexShader.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_PixelShader.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSound.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSound3DCalculator.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundBuffer.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundGlobal.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundLogging.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DirectSoundStream.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/DSStream_PacketManager.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/DirectSound/XFileMediaObject.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/XbDSoundLogging.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/common/XbInternalDSVoice.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/DSOUND/common/XbInternalStruct.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/Intercept.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/JVS/JVS.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/Patches.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/XACTENG/XactEng.cpp"
 "${CXBXR_ROOT_DIR}/src/core/hle/XGRAPHIC/XGraphic.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnl.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlAv.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlDbg.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlEx.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlEx.hpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlFs.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlHal.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlIo.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlKd.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlKe.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlKi.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlMm.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlNt.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlOb.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlPhy.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlPs.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/EmuKrnlRtl.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/exports/KernelThunk.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/memory-manager/PhysicalMemory.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/memory-manager/PoolManager.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/memory-manager/VMManager.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/Emu.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuCdRom.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuDisk.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuFile.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuFS.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuMediaBoard.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuMu.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/EmuNtDll.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/NativeHandle.cpp"
 "${CXBXR_ROOT_DIR}/src/core/kernel/support/PatchRdtsc.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/ADM1032Device.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/Chihiro/JvsIo.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/Chihiro/MediaBoard.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/EEPROMDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/audio/APUDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/network/NVNetDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/MCPXDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/PCIBus.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/PCIDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/SMBus.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/SMCDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/SMDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/usb/Hub.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/usb/OHCI.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/usb/USBDevice.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/usb/XidGamepad.cpp"

 #NOTE: These files are direct includes inside nv2a.cpp file. If we comment them out, we will have compile errors.
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_DEBUG.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PBUS.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PCOUNTER.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PCRTC.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PFB.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PFIFO.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PGRAPH.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PMC.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRAMDAC.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRAMIN.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRMA.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRMCIO.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRMDIO.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRMFB.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PRMVIO.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PSTRAPS.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PTIMER.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PTV.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PVIDEO.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_PVPE.cpp"
 #"${CXBXR_ROOT_DIR}/src/devices/video/EmuNV2A_USER.cpp"

 "${CXBXR_ROOT_DIR}/src/devices/video/nv2a.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/video/nv2a_debug.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/video/qemu-thread-win32.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/video/swizzle.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/x86/EmuX86.cpp"
 "${CXBXR_ROOT_DIR}/src/devices/Xbox.cpp"
 # Temporary usage for need ReserveAddressRanges func with cxbx.exe's emulation.
 "${CXBXR_ROOT_DIR}/src/common/ReserveAddressRanges.cpp"
)


# Debugger build disabled — capstone submodule requires CMake <4.0
#option(BUILD_CXBXR_DEBUGGER "Build cxbxr-debugger tool (with cheat table support)")
set(BUILD_CXBXR_DEBUGGER OFF)

# Set shader output directory early so cxbxr-emu can find it during configure
set(CXBXR_SHADER_GEN_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/shaders" CACHE INTERNAL "")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/projects/cxbx")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/projects/cxbxr-ldr")

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/projects/cxbxr-emu")

set(cxbxr_INSTALL_files "COPYING" "README.md")

# Cxbx-Reloaded project with third-party libraries
set_target_properties(cxbx cxbxr-ldr cxbxr-emu misc-batch SDL2 subhook libXbSymbolDatabase libtommath libtomcrypt imgui libusb
 PROPERTIES FOLDER Cxbx-Reloaded
)

set_target_properties(nv2a_vsh_emulator nv2a_vsh_disassembler nv2a_vsh_cpu
 PROPERTIES FOLDER Cxbx-Reloaded/nv2a_vsh
)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
 # Configure startup project
 set_property(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}" PROPERTY VS_STARTUP_PROJECT cxbx)
endif()

# Check if generator is Visual Studio then enable Cxbxr-Debugger project.
# Since C# is currently supported with Visual Studio for now.
if(${CMAKE_GENERATOR} MATCHES "Visual Studio ([^9]|[9][0-9])" AND BUILD_CXBXR_DEBUGGER)
 # Issues with compile (the same with develop branch) and
 # for some reason did not put the files into virtual folder?
 # Might need to put the list in the source folder for workaround fix.
 add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/src/CxbxDebugger")

 set_target_properties(Tests_cs_x86
  PROPERTIES EXCLUDE_FROM_ALL TRUE
 )

 # Suppress warnings in external capstone code
 if ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  target_compile_options(capstone-shared PRIVATE
   /wd4018 # signed/unsigned mismatch
   /wd4244 # conversion, possible loss of data
   /wd4267 # conversion from size_t, possible loss of data
  )
 endif()

 # Cxbx-Debugger project with third-party libraries
 set_target_properties(cxbxr-debugger cs_x86 Tests_cs_x86 capstone-shared
  PROPERTIES FOLDER Cxbx-Reloaded/debugger
 )
endif()

install(FILES ${cxbxr_INSTALL_files}
  DESTINATION bin
)

# Copy HLSL files to the output directory, which are loaded at runtime
set(CXBXR_HLSL_FILES ${CXBXR_HEADER_EMU})
list(FILTER CXBXR_HLSL_FILES INCLUDE REGEX ".*/src/core/hle/D3D8/Rendering/Shaders/[^/]+\\.hlsli?")
add_custom_command(
    TARGET misc-batch POST_BUILD
    COMMENT "Copying any updated hlsl files..."
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/hlsl
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CXBXR_HLSL_FILES} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/hlsl"
    # These files can be edited.
    # Create backup copies for convenience of restoring original shader behaviour.
    COMMAND ${CMAKE_COMMAND} -E make_directory  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/hlsl/backup
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CXBXR_HLSL_FILES} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/hlsl/backup"
)
install(DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/hlsl DESTINATION bin)

# Build-time HLSL compilation — compile shaders with fxc then embed the
# compiled bytecode directly into cxbxr-emu.dll via generated C headers.
# Dependency chain: .hlsl/.hlsli -> fxc -> .cso -> gen_cso_header -> .h -> .dll
# Editing ANY .hlsl or .hlsli triggers a full recompile+relink automatically.
set(CXBXR_HLSL_DIR "${CXBXR_ROOT_DIR}/src/core/hle/D3D8/Rendering/Shaders")
set(CSO_HEADER_SCRIPT "${CXBXR_HLSL_DIR}/gen_cso_header.cmake")
set(SHADER_OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/shaders")
file(MAKE_DIRECTORY "${SHADER_OUT_DIR}")

# All .hlsli headers that the HLSL files #include (for dependency tracking)
set(CXBXR_HLSLI_DEPS
    "${CXBXR_HLSL_DIR}/CxbxRegisterCombinerInterpreterState.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxVertexShaderInterpreterState.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxFixedFunctionVertexShaderState.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxPixelShaderFunctions.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxPixelShaderInput.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxVertexShaderCommon.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxVertexFetch.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxVertexFetchLayout.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxVertexOutputFooter.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxScreenspaceTransform.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxNV2AMathHelpers.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxNV2APixelShaderConstants.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxNV2AVshOps.hlsli"
    "${CXBXR_HLSL_DIR}/CxbxPGRAPHRegs.hlsli"
)

# Map shader files to their fxc profile and array names
set(SHADER_SOURCES
    CxbxRegisterCombinerInterpreter.hlsl
    CxbxVertexShaderInterpreter.hlsl
    CxbxFixedFunctionVertexShader.hlsl
    CxbxBlitVS.hlsl
    CxbxBlitPS.hlsl
    CxbxPointSpriteGS.hlsl
    CxbxThickLineGS.hlsl
    CxbxUnswizzleCS.hlsl
    CxbxUnswizzleBGRA_CS.hlsl
    CxbxIndexConvertCS.hlsl
    CxbxPaletteExpandCS.hlsl
    CxbxFormatConvertCS.hlsl
    CxbxVertexConvertCS.hlsl
)
set(SHADER_PROFILES
    ps_5_0
    vs_5_0
    vs_5_0
    vs_5_0
    ps_5_0
    gs_5_0
    gs_5_0
    cs_5_0
    cs_5_0
    cs_5_0
    cs_5_0
    cs_5_0
    cs_5_0
)
set(SHADER_ARRAY_NAMES
    cso_CxbxRCInterpreterPS
    cso_CxbxVSInterpreterVS
    cso_CxbxFixedFunctionVS
    cso_CxbxBlitVS
    cso_CxbxBlitPS
    cso_CxbxPointSpriteGS
    cso_CxbxThickLineGS
    cso_CxbxUnswizzleCS
    cso_CxbxUnswizzleBGRA_CS
    cso_CxbxIndexConvertCS
    cso_CxbxPaletteExpandCS
    cso_CxbxFormatConvertCS
    cso_CxbxVertexConvertCS
)

set(CSO_HEADERS)
list(LENGTH SHADER_SOURCES SHADER_COUNT)
math(EXPR SHADER_LAST "${SHADER_COUNT} - 1")
foreach(I RANGE 0 ${SHADER_LAST})
    list(GET SHADER_SOURCES ${I} SHADER_FILE)
    list(GET SHADER_PROFILES ${I} PROFILE)
    list(GET SHADER_ARRAY_NAMES ${I} ARRAY_NAME)
    set(HLSL_PATH "${CXBXR_HLSL_DIR}/${SHADER_FILE}")
    set(CSO_PATH "${SHADER_OUT_DIR}/${ARRAY_NAME}.cso")
    set(HDR_PATH "${SHADER_OUT_DIR}/${ARRAY_NAME}.h")

    # Step 1: HLSL -> CSO (fxc)
    add_custom_command(
        OUTPUT "${CSO_PATH}"
        COMMAND fxc /nologo /T ${PROFILE} /E main /O1
                /I "${CXBXR_HLSL_DIR}"
                /Fo "${CSO_PATH}"
                "${HLSL_PATH}"
        DEPENDS "${HLSL_PATH}" ${CXBXR_HLSLI_DEPS}
        COMMENT "HLSL->CSO: ${SHADER_FILE} (${PROFILE})"
        VERBATIM
    )

    # Step 2: CSO -> C header
    add_custom_command(
        OUTPUT "${HDR_PATH}"
        COMMAND ${CMAKE_COMMAND}
                -DCSO_FILE=${CSO_PATH}
                -DHEADER_FILE=${HDR_PATH}
                -DARRAY_NAME=${ARRAY_NAME}
                -P "${CSO_HEADER_SCRIPT}"
        DEPENDS "${CSO_PATH}" "${CSO_HEADER_SCRIPT}"
        COMMENT "CSO->header: ${ARRAY_NAME}.h"
        VERBATIM
    )

    list(APPEND CSO_HEADERS "${HDR_PATH}")
endforeach()

# ---------------------------------------------------------------------------
# Specialized shader variants — compile the same HLSL source with different
# /D defines to produce optimized variants for common configurations.
# ---------------------------------------------------------------------------

# Helper: compile a variant with extra /D defines
# Args: HLSL_FILE PROFILE ARRAY_NAME DEFINES (semicolon-separated /Dkey=val)
function(compile_shader_variant HLSL_FILE PROFILE ARRAY_NAME)
    set(EXTRA_DEFINES ${ARGN})
    set(HLSL_PATH "${CXBXR_HLSL_DIR}/${HLSL_FILE}")
    set(CSO_PATH "${SHADER_OUT_DIR}/${ARRAY_NAME}.cso")
    set(HDR_PATH "${SHADER_OUT_DIR}/${ARRAY_NAME}.h")

    add_custom_command(
        OUTPUT "${CSO_PATH}"
        COMMAND fxc /nologo /T ${PROFILE} /E main /O1
                /I "${CXBXR_HLSL_DIR}"
                ${EXTRA_DEFINES}
                /Fo "${CSO_PATH}"
                "${HLSL_PATH}"
        DEPENDS "${HLSL_PATH}" ${CXBXR_HLSLI_DEPS}
        COMMENT "HLSL->CSO variant: ${ARRAY_NAME} (${PROFILE})"
        VERBATIM
    )
    add_custom_command(
        OUTPUT "${HDR_PATH}"
        COMMAND ${CMAKE_COMMAND}
                -DCSO_FILE=${CSO_PATH}
                -DHEADER_FILE=${HDR_PATH}
                -DARRAY_NAME=${ARRAY_NAME}
                -P "${CSO_HEADER_SCRIPT}"
        DEPENDS "${CSO_PATH}" "${CSO_HEADER_SCRIPT}"
        COMMENT "CSO->header: ${ARRAY_NAME}.h"
        VERBATIM
    )
    # Propagate to parent scope
    set(CSO_HEADERS ${CSO_HEADERS} "${HDR_PATH}" PARENT_SCOPE)
endfunction()

# RC Pixel Shader — single base variant (JIT replaces specialized variants)
# VS Interpreter — single base variant (JIT replaces NoCtx optimization)

add_custom_target(compile-shaders ALL DEPENDS ${CSO_HEADERS})
set_target_properties(compile-shaders PROPERTIES FOLDER "Cxbx-Reloaded/tools")

install(PROGRAMS $<TARGET_FILE_DIR:cxbx>/${CMAKE_SHARED_LIBRARY_PREFIX}subhook${CMAKE_SHARED_LIBRARY_SUFFIX}
  DESTINATION bin
)

install(PROGRAMS $<TARGET_FILE_DIR:cxbx>/${CMAKE_SHARED_LIBRARY_PREFIX}SDL2${CMAKE_SHARED_LIBRARY_SUFFIX}
  DESTINATION bin
)
