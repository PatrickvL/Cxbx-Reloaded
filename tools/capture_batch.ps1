# capture_batch.ps1 - Capture screenshots for multiple XDK samples
param(
    [string[]]$Samples = @("Vertices", "Textures", "BumpEarth", "Dolphin"),
    [string]$SamplesDir = "C:\Users\patrick.vanlogchem\Downloads\Xbox\XDK_Samples\Compiled",
    [int]$DelayMs = 12000,
    [int]$Frame2DelayMs = 2000,
    [switch]$ClearCache
)

$scriptPath = "c:\Workspaces\Mine\Cxbx-Reloaded\tools\capture_emulator.ps1"

foreach ($sample in $Samples) {
    $xbePath = Join-Path $SamplesDir "$sample\$sample.xbe"
    if (!(Test-Path $xbePath)) {
        # Try default.xbe
        $xbePath = Join-Path $SamplesDir "$sample\default.xbe"
        if (!(Test-Path $xbePath)) {
            Write-Warning "Skipping $sample - no XBE found"
            continue
        }
    }
    Write-Host "`n=== Capturing $sample ===" -ForegroundColor Cyan
    $args = @(
        "-ExecutionPolicy", "Bypass",
        "-File", $scriptPath,
        "-XbePath", $xbePath,
        "-KillExisting",
        "-StopAfter",
        "-DelayMs", $DelayMs,
        "-Frame2DelayMs", $Frame2DelayMs
    )
    if ($ClearCache) { $args += "-ClearCache" }
    & powershell @args
    Start-Sleep -Seconds 2  # cooldown between launches
}

Write-Host "`n=== All captures complete ===" -ForegroundColor Green
Get-ChildItem "c:\Workspaces\Mine\Cxbx-Reloaded\build-x86\bin\Release\screenshots\*.png" | 
    Sort-Object LastWriteTime | 
    Format-Table Name, LastWriteTime -AutoSize
