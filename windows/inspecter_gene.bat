@echo off
REM ===================================================================
REM  Couverture intron/exon d'UN gene, tous echantillons, sans seuil.
REM  Usage : inspecter_gene.bat SMC1A [echantillon]
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat
if "%~1"=="" ( echo Usage : inspecter_gene.bat SMC1A [echantillon] & pause & exit /b 1 )
if "%~2"=="" (
  python -m splicefast gene -w "%WORK%" --gene %1
) else (
  python -m splicefast gene -w "%WORK%" --gene %1 --sample %2
)
pause
