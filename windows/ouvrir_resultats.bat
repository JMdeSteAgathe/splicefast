@echo off
setlocal
cd /d "%~dp0"
call config.bat
start "" "%WORK%\results\SpliceFast_viewer.html"
