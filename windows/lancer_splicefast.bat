@echo off
REM ===================================================================
REM  SpliceFast 2.5 - analyse complete des runs de config.bat
REM  Comptage + couverture des introns en UNE lecture des BAM (lecteur
REM  C++ sfscan.exe fourni), puis matrices, tests, PanelApp et rapport.
REM  Relancer ce fichier apres un incident : le travail fait est garde.
REM ===================================================================
setlocal
cd /d "%~dp0"
call config.bat

if not exist "%GTF%" (
  echo.
  echo  ERREUR : annotation introuvable : %GTF%
  echo  Editez config.bat, ligne "set GTF=..."
  pause & exit /b 1
)
for %%R in ("%RUN_A%" "%RUN_B%") do if not exist "%%~R" (
  echo  ERREUR : dossier de run introuvable : %%~R  ^(editez config.bat^)
  pause & exit /b 1
)

echo.
echo  === Estimation de la duree ===
python -m splicefast bench --run "%RUN_A%" "%RUN_B%" -w "%WORK%" -t %THREADS%
echo.
echo  Appuyez sur une touche pour lancer l'analyse complete (Ctrl-C pour annuler).
pause >nul

REM  PanelApp : tous les panels (une fois ; mise_a_jour_panelapp.bat pour rafraichir)
if not exist "%WORK%\panelapp\panelapp_genes.tsv.gz" (
  echo  === Telechargement PanelApp ===
  call "%~dp0mise_a_jour_panelapp.bat" nopause
)

python -m splicefast all ^
  --run "%RUN_A%" "%RUN_B%" ^
  --gtf "%GTF%" ^
  -w "%WORK%" ^
  -t %THREADS% -j %JOBS% ^
  --retention ^
  --reference-match %REF_MATCH% --min-reads %MIN_READS% ^
  --io-retries 30 --io-wait 5

if errorlevel 1 (
  echo.
  echo  L'analyse s'est arretee sur une erreur. Relancez ce fichier :
  echo  les etapes deja terminees ne seront pas refaites.
  echo  Journal detaille : %WORK%\logs\splicefast.log
  pause & exit /b 1
)
echo.
echo  Termine. Visualiseur : %WORK%\results\SpliceFast_viewer.html
echo  (dans la page : "Open results folder" puis choisir %WORK%\results)
start "" "%WORK%\results\SpliceFast_viewer.html"
pause
