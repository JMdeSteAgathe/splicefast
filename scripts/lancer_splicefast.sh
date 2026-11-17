#!/usr/bin/env bash
# Analyse complete : comptage + couverture des introns en UNE lecture des BAM,
# puis matrices, tests et rapport. Relancer la meme commande apres un incident :
# le travail deja fait est garde.
source "$(dirname "$0")/_commun.sh"
[ -f "$GTF" ] || { echo "ERREUR : annotation introuvable : $GTF (editer scripts/config.sh)" >&2; exit 1; }
[ -x splicefast/bin/sfscan ] || echo "ATTENTION : sfscan non compile (make) -> lecteur Python, bien plus lent"

echo "=== Estimation de la duree ==="
sf bench --run "${RUNS[@]}" -w "$WORK" -t "$THREADS" || true
read -r -p "Entree pour lancer l'analyse complete (Ctrl-C pour annuler) " _

mapfile -d '' P < <(panels)
sf all --run "${RUNS[@]}" --gtf "$GTF" -w "$WORK" \
       -t "$THREADS" -j "$JOBS" --retention "${P[@]}" \
       --io-retries 30 --io-wait 5 || {
  echo; echo "Arret sur erreur. Journal : $WORK/logs/splicefast.log"
  echo "Relancez ce script : les etapes terminees ne seront pas refaites."; exit 1; }
ouvrir "$WORK/results/SpliceFast_viewer.html"
