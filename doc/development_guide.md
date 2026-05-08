# Development Guide

## Building

Both targets must match — always build both:
```powershell
cmake --build . --target cxbxr-emu --config Release -- /m
cmake --build . --target cxbx --config Release -- /m
```
If only cxbxr-emu.dll is rebuilt: "Mismatch detected between EmuShared and cxbx.exe/cxbxr-emu.dll"

Build directory: `build/` (relative to repo root)

## Launching XBEs

**GUI with arg:**
```powershell
Start-Process .\cxbx.exe -ArgumentList '"path\to\game.xbe"'
```

**CLI:**
```powershell
cd build\bin\Release
cmd /c '"$((Resolve-Path .\cxbxr-ldr.exe).Path)" /load "path\to\game.xbe"'
```
- MUST use full path to exe (cliConverter.cpp argv[0] handling)
- MUST use `/load` argument

## Test Samples

The test suite uses compiled XBE samples. Set the environment variable `CXBX_XBE_SAMPLES` to point
to your local samples directory. Expected layout:
- `%CXBX_XBE_SAMPLES%/Dolphin/Dolphin.xbe`
- `%CXBX_XBE_SAMPLES%/Meshes/Meshes.xbe`
- etc.

Do NOT use `default.xbe` — use the named XBE (e.g., `Dolphin.xbe`).

## Shader Cache

Clear before testing shader changes:
```powershell
Remove-Item "$env:APPDATA\Cxbx-Reloaded\ShaderCache" -Recurse -Force
```

## Profiler

- Header: `src/core/hle/D3D8/Rendering/Backend/Backend_D3D11_Profiler.h`
- Output: `%TEMP%\CxbxProfiler.log` + printf to stdout + OutputDebugStringA
- Triggers once per second after first Present (via `CxbxProfilerFrameTick()`)
- MMIO count: bare `InterlockedIncrement` in VEH — no QPC (too hot)

## VEH Architecture

- Single VEH: `lleException` in Emu.cpp — registered with priority first
  1. Checks PageTracker faults (AV in contiguous 0x80000000 or tiled 0xF0000000 region)
  2. Falls through to LLE MMIO decode (EmuX86_DecodeException)
- `EmuException`: SetUnhandledExceptionFilter — last resort popup
- `IsXboxCodeAddress(addr)`: checks `addr >= XBE_IMAGE_BASE && addr <= XBE_MAX_VA`

---

## D3D11 Texture Upload Pitfalls

### UpdateSubresource vs Map/Unmap
- `UpdateSubresource()` is **INVALID for D3D11_USAGE_DYNAMIC** — silently dropped
- DYNAMIC textures MUST use `Map(D3D11_MAP_WRITE_DISCARD)/Unmap`
- DEFAULT textures can use either `UpdateSubresource` or staging copy

### CS Unswizzle UAV Compatibility
- Formats like B4G4R4A4_UNORM do NOT support typed UAV access
- `CheckFormatSupport()` with `D3D11_FORMAT_SUPPORT_TYPED_UNORDERED_ACCESS_VIEW` returns false
- Fallback: texture created as DYNAMIC, CPU `EmuUnswizzleBox` + Map/Unmap

### Texture Upload Path (HostResourceUpload.cpp)
1. CS unswizzle: requires DEFAULT+UAV texture, dispatches compute shader
2. CS fallback: CPU EmuUnswizzleBox + Map/Unmap (DYNAMIC) or UpdateSubresource (DEFAULT)
3. Normal path: Map/Unmap for DYNAMIC, staging buffer for DEFAULT (multi-mip)

See [nv2a_resource_management.md](nv2a_resource_management.md) §1.1 for the complete texture lifecycle.

### Texture Creation (HostResourceCreate.cpp)
- Single-mip swizzled textures: check `bCanUAV` via CheckFormatSupport
  - UAV supported → DEFAULT + BIND_SHADER_RESOURCE | BIND_UNORDERED_ACCESS
  - UAV not supported → DYNAMIC + BIND_SHADER_RESOURCE (CPU_ACCESS_WRITE)
