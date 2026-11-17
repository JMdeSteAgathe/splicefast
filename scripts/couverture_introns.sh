#!/usr/bin/env bash
# Module retention d'intron par la couverture, sur TOUS les introns.
# Inutile si l'analyse a ete lancee avec --retention (lancer_splicefast.sh le
# fait) : les mesures sont deja en cache, seul l'assemblage est refait.
source "$(dirname "$0")/_commun.sh"
sf retention -w "$WORK" -t "$THREADS" -j "$JOBS"
mapfile -d '' P < <(panels)
sf report -w "$WORK" "${P[@]}"
ouvrir "$WORK/results/SpliceFast_viewer.html"
