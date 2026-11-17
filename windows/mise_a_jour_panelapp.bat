@echo off
REM ===================================================================
REM  Telecharge PanelApp Australie (tous les panels, confiance de chaque
REM  gene) dans %WORK%\panelapp, puis il suffit de refaire le rapport.
REM  1) par Python ; 2) a defaut, par PowerShell (proxy / certificats
REM  de l'hopital mieux geres par Windows).
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
python -m splicefast panelapp --update -w "%WORK%"
if not errorlevel 1 goto ok
echo.
echo  Echec en Python : nouvel essai via PowerShell...
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0panelapp_download.ps1" -OutDir "%WORK%\panelapp\json"
if errorlevel 1 (
  echo  PanelApp injoignable. L'analyse continue sans annotation de panels.
  if not "%~1"=="nopause" pause
  exit /b 0
)
python -m splicefast panelapp --from-json "%WORK%\panelapp\json" -w "%WORK%"
:ok
echo  PanelApp a jour. Refaire le rapport : refaire_rapport.bat
if not "%~1"=="nopause" pause
