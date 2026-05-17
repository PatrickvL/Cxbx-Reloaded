# capture_emulator.ps1 - Capture screenshots from the Cxbx-Reloaded emulator window
# Usage: .\capture_emulator.ps1 -XbePath "path\to\game.xbe" -OutDir "path\to\output" [-DelayMs 5000] [-Frame2DelayMs 800]
#
# Launches the emulator, waits for rendering, captures two screenshots:
#   1) After initial delay (default 5s after window appears) - "early" frame
#   2) After additional delay (default 800ms later) - "late" frame
# This dual-capture helps detect black-screen flicker issues.

param(
    [Parameter(Mandatory=$true)]
    [string]$XbePath,

    [string]$OutDir,

    [string]$EmulatorPath,

    # Milliseconds to wait after window appears before first capture
    [int]$DelayMs = 8000,

    # Milliseconds between first and second capture
    [int]$Frame2DelayMs = 1500,

    # If true, clear shader cache before launch (cold start)
    [switch]$ClearCache,

    # If true, kill any existing emulator instances first
    [switch]$KillExisting,

    # If true, stop the emulator after capturing
    [switch]$StopAfter,

    # Override the output file base name (defaults to XBE filename without extension)
    [string]$Name
)

$RepoRoot = (Resolve-Path "$PSScriptRoot\..").Path
if (-not $OutDir) { $OutDir = "$RepoRoot\test_screenshots" }
if (-not $EmulatorPath) { $EmulatorPath = "$RepoRoot\build\bin\Release\cxbx.exe" }

Add-Type @"
using System;
using System.Runtime.InteropServices;
using System.Drawing;
using System.Drawing.Imaging;

public class WindowCapture {
    [DllImport("user32.dll")]
    public static extern IntPtr FindWindow(string lpClassName, string lpWindowName);

    [DllImport("user32.dll")]
    public static extern bool GetWindowRect(IntPtr hWnd, out RECT lpRect);

    [DllImport("user32.dll")]
    public static extern bool GetClientRect(IntPtr hWnd, out RECT lpRect);

    [DllImport("user32.dll")]
    public static extern IntPtr GetDC(IntPtr hWnd);

    [DllImport("user32.dll")]
    public static extern int ReleaseDC(IntPtr hWnd, IntPtr hDC);

    [DllImport("gdi32.dll")]
    public static extern bool BitBlt(IntPtr hdcDest, int xDest, int yDest, int wDest, int hDest,
        IntPtr hdcSrc, int xSrc, int ySrc, int rop);

    [DllImport("user32.dll")]
    public static extern bool ClientToScreen(IntPtr hWnd, ref POINT point);

    [DllImport("user32.dll")]
    public static extern bool SetForegroundWindow(IntPtr hWnd);

    [DllImport("user32.dll")]
    public static extern bool ShowWindow(IntPtr hWnd, int nCmdShow);

    [DllImport("user32.dll")]
    public static extern bool BringWindowToTop(IntPtr hWnd);

    [DllImport("user32.dll")]
    public static extern IntPtr GetForegroundWindow();

    [DllImport("user32.dll")]
    public static extern uint GetWindowThreadProcessId(IntPtr hWnd, out uint processId);

    [DllImport("user32.dll")]
    public static extern bool AttachThreadInput(uint idAttach, uint idAttachTo, bool fAttach);

    [DllImport("kernel32.dll")]
    public static extern uint GetCurrentThreadId();

    [DllImport("user32.dll")]
    public static extern bool IsWindowVisible(IntPtr hWnd);

    [DllImport("user32.dll")]
    public static extern bool PrintWindow(IntPtr hWnd, IntPtr hdcBlt, uint nFlags);

    [DllImport("user32.dll")]
    public static extern int GetWindowTextLength(IntPtr hWnd);

    [DllImport("user32.dll", CharSet = CharSet.Auto)]
    public static extern int GetWindowText(IntPtr hWnd, System.Text.StringBuilder lpString, int nMaxCount);

    [DllImport("user32.dll")]
    public static extern bool EnumWindows(EnumWindowsProc lpEnumFunc, IntPtr lParam);

    public delegate bool EnumWindowsProc(IntPtr hWnd, IntPtr lParam);

    [StructLayout(LayoutKind.Sequential)]
    public struct RECT {
        public int Left, Top, Right, Bottom;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct POINT {
        public int X, Y;
    }

    public const int SRCCOPY = 0x00CC0020;
    public const int SW_RESTORE = 9;
    public const uint PW_RENDERFULLCONTENT = 0x00000002;

    // Force a window to the foreground using the AttachThreadInput trick
    public static void ForceForeground(IntPtr hWnd) {
        IntPtr foreWnd = GetForegroundWindow();
        uint foreThread, curThread;
        uint pid;
        foreThread = GetWindowThreadProcessId(foreWnd, out pid);
        curThread = GetCurrentThreadId();

        if (foreThread != curThread) {
            AttachThreadInput(curThread, foreThread, true);
            ShowWindow(hWnd, SW_RESTORE);
            BringWindowToTop(hWnd);
            SetForegroundWindow(hWnd);
            AttachThreadInput(curThread, foreThread, false);
        } else {
            ShowWindow(hWnd, SW_RESTORE);
            BringWindowToTop(hWnd);
            SetForegroundWindow(hWnd);
        }
    }

    // Find the emulator render window by process ID
    public static IntPtr FindWindowByPid(uint pid) {
        IntPtr found = IntPtr.Zero;
        EnumWindows(delegate(IntPtr hWnd, IntPtr lParam) {
            uint wpid;
            GetWindowThreadProcessId(hWnd, out wpid);
            if (wpid == pid && IsWindowVisible(hWnd)) {
                int len = GetWindowTextLength(hWnd);
                if (len > 0) {
                    var sb = new System.Text.StringBuilder(len + 1);
                    GetWindowText(hWnd, sb, sb.Capacity);
                    string title = sb.ToString();
                    // The emulator window title typically contains "Cxbx" or "CxbxReloaded"
                    if (title.IndexOf("Cxbx", StringComparison.OrdinalIgnoreCase) >= 0) {
                        found = hWnd;
                        return false; // stop enumerating
                    }
                }
            }
            return true;
        }, IntPtr.Zero);
        return found;
    }

    // Capture a window's client area to a bitmap file
    public static bool CaptureWindow(IntPtr hWnd, string outputPath) {
        if (hWnd == IntPtr.Zero) return false;

        RECT windowRect;
        GetWindowRect(hWnd, out windowRect);
        int width = windowRect.Right - windowRect.Left;
        int height = windowRect.Bottom - windowRect.Top;

        if (width <= 0 || height <= 0) return false;

        // Always use BitBlt from screen — PrintWindow doesn't capture D3D11 swapchain content
        using (Bitmap bmp = new Bitmap(width, height, PixelFormat.Format32bppArgb)) {
            using (Graphics g = Graphics.FromImage(bmp)) {
                IntPtr hdcScreen = GetDC(IntPtr.Zero);
                IntPtr hdcBmp = g.GetHdc();
                BitBlt(hdcBmp, 0, 0, width, height, hdcScreen, windowRect.Left, windowRect.Top, SRCCOPY);
                g.ReleaseHdc(hdcBmp);
                ReleaseDC(IntPtr.Zero, hdcScreen);
            }
            bmp.Save(outputPath, ImageFormat.Png);
        }
        return true;
    }
}
"@ -ReferencedAssemblies System.Drawing

# Ensure output directory exists
if (!(Test-Path $OutDir)) {
    New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
}

# Kill existing instances if requested
if ($KillExisting) {
    Get-Process -Name cxbx, cxbxr-ldr, cxbxr-emu -ErrorAction SilentlyContinue | Stop-Process -Force
    Start-Sleep -Milliseconds 500
}

# Clear shader cache if requested
if ($ClearCache) {
    $cachePath = Join-Path (Split-Path $EmulatorPath) "ShaderCache"
    Remove-Item "$cachePath\*" -Recurse -Force -ErrorAction SilentlyContinue
    Write-Host "Shader cache cleared"
}

# Deploy latest HLSL
$hlslSrc = "$RepoRoot\src\core\hle\D3D8\Rendering\Shaders"
$hlslDst = Split-Path $EmulatorPath
$hlslDstDir = Join-Path $hlslDst "hlsl"
if (Test-Path $hlslSrc) {
    if (!(Test-Path $hlslDstDir)) { New-Item -ItemType Directory -Path $hlslDstDir -Force | Out-Null }
    Get-ChildItem "$hlslSrc\*" -Include "*.hlsl","*.hlsli" | Copy-Item -Destination "$hlslDstDir\" -Force
    Write-Host "HLSL deployed"
}

# Extract sample name for filenames
if ($Name) {
    $sampleName = $Name
} else {
    $sampleName = [System.IO.Path]::GetFileNameWithoutExtension($XbePath)
}

# Launch emulator
Write-Host "Launching $sampleName..."
$proc = Start-Process -FilePath $EmulatorPath -ArgumentList "`"$XbePath`"" -PassThru

# Wait for the emulator render window to appear
Write-Host "Waiting for emulator window..."
$hWnd = [IntPtr]::Zero
$maxWait = 60  # seconds
$waited = 0
while ($hWnd -eq [IntPtr]::Zero -and $waited -lt $maxWait) {
    Start-Sleep -Milliseconds 500
    $waited += 0.5

    # Search all emulator-related processes for a visible "Cxbx" window
    foreach ($procName in @("cxbx", "cxbxr-ldr", "cxbxr-emu")) {
        $procs = Get-Process -Name $procName -ErrorAction SilentlyContinue
        foreach ($p in $procs) {
            $hWnd = [WindowCapture]::FindWindowByPid([uint32]$p.Id)
            if ($hWnd -ne [IntPtr]::Zero) { break }
        }
        if ($hWnd -ne [IntPtr]::Zero) { break }
    }
}

if ($hWnd -eq [IntPtr]::Zero) {
    Write-Error "Could not find emulator window after ${maxWait}s"
    exit 1
}

Write-Host "Found emulator window (handle=$hWnd) after ${waited}s"

# Wait for rendering to stabilize
Write-Host "Waiting ${DelayMs}ms for rendering to stabilize..."
Start-Sleep -Milliseconds $DelayMs

# Capture first screenshot (early)
$ts = Get-Date -Format "yyyyMMdd_HHmmss"
$earlyPath = Join-Path $OutDir "${sampleName}_early_${ts}.png"
[WindowCapture]::ForceForeground($hWnd)
Start-Sleep -Milliseconds 500

$ok = [WindowCapture]::CaptureWindow($hWnd, $earlyPath)
if ($ok) {
    Write-Host "Early capture saved: $earlyPath"
} else {
    Write-Error "Failed to capture early screenshot"
}

# Wait and capture second screenshot (late)
Write-Host "Waiting ${Frame2DelayMs}ms for late capture..."
Start-Sleep -Milliseconds $Frame2DelayMs

$latePath = Join-Path $OutDir "${sampleName}_late_${ts}.png"
$ok = [WindowCapture]::CaptureWindow($hWnd, $latePath)
if ($ok) {
    Write-Host "Late capture saved: $latePath"
} else {
    Write-Error "Failed to capture late screenshot"
}

# Stop emulator if requested
if ($StopAfter) {
    Get-Process -Name cxbx, cxbxr-ldr, cxbxr-emu -ErrorAction SilentlyContinue | Stop-Process -Force
    Write-Host "Emulator stopped"
}

# Output paths for easy consumption
Write-Output "EARLY=$earlyPath"
Write-Output "LATE=$latePath"
