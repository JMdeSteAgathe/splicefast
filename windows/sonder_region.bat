@echo off
REM ===================================================================
REM  Relit une region DIRECTEMENT dans le BAM d'un echantillon et la compare
REM  a ce que SpliceFast a compte (evenement vu dans IGV mais absent).
REM  Usage : sonder_region.bat ECH02-PAX_S61 chr5:37025000-37028000
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
if "%~2"=="" ( echo Usage : sonder_region.bat ECHANTILLON chr5:37025000-37028000 & pause & exit /b 1 )
python -m splicefast probe -w "%WORK%" --sample %1 --region %2
echo.
echo  Si BAM, comptage et matrice concordent : pourquoi.bat %1 GENE
pause
