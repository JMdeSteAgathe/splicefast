@echo off
REM ===================================================================
REM  Regenere UNIQUEMENT le rapport (quelques minutes, aucun test refait).
REM  Usage : refaire_rapport.bat [padj=0.05] [delta=0.1]
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
set PADJ=%~1
if "%PADJ%"=="" set PADJ=0.05
set DELTA=%~2
if "%DELTA%"=="" set DELTA=0.1
python -m splicefast report -w "%WORK%" --padj %PADJ% --delta %DELTA%
if errorlevel 1 ( pause & exit /b 1 )
start "" "%WORK%\results\SpliceFast_viewer.html"
