@echo off
REM ===================================================================
REM  Le disque tient-il la charge ? Lecture continue pendant 20 minutes.
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
echo  Avant : Options d'alimentation ^> Disque dur ^> Arreter apres : Jamais
echo          Parametres USB ^> Suspension selective USB : Desactive
python -m splicefast bench --run "%RUN_A%" "%RUN_B%" -w "%WORK%" -t %THREADS% --stress 20
pause
