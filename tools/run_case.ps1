# Usage: powershell -File tools\run_case.ps1 -Log run02 -Code work\OpenXY\Code -Data work\pt10x10 -Stem Pt10x10 [-TimeoutSec 900]
# Runs run_openxy_headless under matlab -batch with stdout+stderr captured to logs\<Log>.log;
# kills MATLAB if it exceeds the timeout. Exit code of MATLAB is appended to the log.
param([string]$Log, [string]$Code, [string]$Data, [string]$Stem, [int]$TimeoutSec = 900)
$r = (Get-Location).Path
$out = "$r\logs\$Log.stdout"; $err = "$r\logs\$Log.stderr"
$cmd = "try, addpath('$r\tools'); run_openxy_headless('$r\$Code','$r\$Data','$Stem','$Log.ang'); catch ME, fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); exit(3); end; exit(0)"
$p = Start-Process -FilePath "C:\Program Files\MATLAB\R2024a\bin\matlab.exe" -ArgumentList "-batch `"$cmd`"" `
     -PassThru -WindowStyle Hidden -RedirectStandardOutput $out -RedirectStandardError $err
$timedOut = $false
if (-not $p.WaitForExit($TimeoutSec * 1000)) {
    $timedOut = $true
    Get-Process MATLAB -ErrorAction SilentlyContinue | Stop-Process -Force
}
$p.WaitForExit()
$lines = @("=== run_case $Log  $(Get-Date -Format s)  code=$Code data=$Data stem=$Stem ===", "--- stdout ---") + (Get-Content $out -ErrorAction SilentlyContinue) + "--- stderr ---" + (Get-Content $err -ErrorAction SilentlyContinue)
if ($timedOut) { $lines += "*** TIMEOUT after $TimeoutSec s - MATLAB killed ***" } else { $lines += "*** MATLAB exit code $($p.ExitCode) ***" }
$lines | Set-Content "$r\logs\$Log.log" -Encoding utf8
Remove-Item $out, $err -ErrorAction SilentlyContinue
$lines
