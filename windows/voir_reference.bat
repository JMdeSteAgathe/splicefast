@echo off
REM Liste les echantillons : reference (norme) ou teste seulement, avec leurs reads.
setlocal
cd /d "%~dp0"
call config.bat
python -m splicefast reference -w "%WORK%" --list
pause
