# Install the latest searchSNP release on Windows.
#
#   irm https://raw.githubusercontent.com/padpadpadpad/searchSNP/main/install.ps1 | iex
#
# Environment variable:
#   INSTALL_DIR   where to put searchSNP.exe (default: %LOCALAPPDATA%\searchSNP\bin)
$ErrorActionPreference = 'Stop'

$repo = 'padpadpadpad/searchSNP'
$bin  = 'searchSNP'

$arch = [System.Runtime.InteropServices.RuntimeInformation]::OSArchitecture
if ($arch -ne 'X64') {
    throw "Unsupported architecture: $arch. Only x64 builds are published (ARM64 Windows can run the x64 build under emulation)."
}

$target = 'x86_64-pc-windows-msvc'
$url = "https://github.com/$repo/releases/latest/download/$bin-$target.zip"

if ($env:INSTALL_DIR) {
    $installDir = $env:INSTALL_DIR
} else {
    $installDir = Join-Path $env:LOCALAPPDATA "$bin\bin"
}

$tmpZip = Join-Path ([System.IO.Path]::GetTempPath()) "$bin-$target.zip"

Write-Host "Downloading $url"
Invoke-WebRequest -Uri $url -OutFile $tmpZip -UseBasicParsing

New-Item -ItemType Directory -Force $installDir | Out-Null
Expand-Archive -Path $tmpZip -DestinationPath $installDir -Force
Remove-Item $tmpZip

$exe = Join-Path $installDir "$bin.exe"
Write-Host "Installed $exe"

$userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
if (-not $userPath) { $userPath = '' }
if (($userPath -split ';') -notcontains $installDir) {
    $newPath = if ($userPath) { "$userPath;$installDir" } else { $installDir }
    [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')
    $env:Path = "$env:Path;$installDir"
    Write-Host "Added $installDir to your user PATH. Open a new terminal for it to take effect."
}

& $exe --version
