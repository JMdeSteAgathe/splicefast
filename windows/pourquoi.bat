@echo off
REM ===================================================================
REM  Pourquoi un evenement attendu n'est-il pas dans le rapport ?
REM  Usage : pourquoi.bat ECH02-PAX_S61 NIPBL
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
if "%~2"=="" ( echo Usage : pourquoi.bat ECHANTILLON GENE & pause & exit /b 1 )
python -m splicefast why -w "%WORK%" --sample %1 --gene %2
pause
