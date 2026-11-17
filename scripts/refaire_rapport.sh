#!/usr/bin/env bash
# Regenere UNIQUEMENT le rapport (aucun test recalcule) : pour changer seuils/panels.
# Usage : scripts/refaire_rapport.sh [padj=0.05] [delta=0.1]
source "$(dirname "$0")/_commun.sh"
mapfile -d '' P < <(panels)
sf report -w "$WORK" --padj "${1:-0.05}" --delta "${2:-0.1}" "${P[@]}"
ouvrir "$WORK/results/SpliceFast_viewer.html"
