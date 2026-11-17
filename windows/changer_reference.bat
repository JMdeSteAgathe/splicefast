@echo off
REM ===================================================================
REM  Change la reference (REF_MATCH / MIN_READS de config.bat) SANS rien
REM  recompter, puis refait tests et rapport (quelques dizaines de minutes).
REM  La liste reference / hors reference est affichee avant.
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
python -m splicefast reference -w "%WORK%" --reference-match %REF_MATCH% --min-reads %MIN_READS% --list
if errorlevel 1 ( pause & exit /b 1 )
echo.
echo  Appuyez sur une touche pour refaire l'analyse avec cette reference (Ctrl-C pour annuler).
pause >nul
python -m splicefast analyze -w "%WORK%" && python -m splicefast report -w "%WORK%"
if errorlevel 1 ( pause & exit /b 1 )
start "" "%WORK%\results\SpliceFast_viewer.html"
pause
