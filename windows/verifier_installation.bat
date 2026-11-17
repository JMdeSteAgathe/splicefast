@echo off
REM Verifie Python, numpy, scipy, le lecteur C++ sfscan.exe, puis l'auto-test.
setlocal
cd /d "%~dp0"
python -c "import sys; print('Python', sys.version)" || (echo Python introuvable & pause & exit /b 1)
python -c "import numpy, scipy; print('numpy', numpy.__version__, '/ scipy', scipy.__version__)" || (
  echo.
  echo  numpy ou scipy manquant. Installation sans droits admin :
  echo      python -m pip install --user numpy scipy
  pause & exit /b 1
)
python -m splicefast --version
"%~dp0splicefast\bin\sfscan.exe" --version || echo  ATTENTION : sfscan.exe ne se lance pas (antivirus ?) : lecteur Python, beaucoup plus lent
echo.
echo  === Auto-test sur donnees simulees ===
python -m splicefast selftest -w "%TEMP%\splicefast_selftest"
echo.
echo  Vous devez lire "Auto-test : 4/4 evenements retrouves"
echo  et "lecteur C++ vs Python sur 4 echantillons : IDENTIQUES".
pause
