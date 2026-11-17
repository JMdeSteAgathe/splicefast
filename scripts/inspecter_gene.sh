#!/usr/bin/env bash
# Couverture intron/exon d'UN gene pour tous les echantillons, sans seuil.
# Usage : scripts/inspecter_gene.sh SMC1A [identifiant_echantillon]
source "$(dirname "$0")/_commun.sh"
[ $# -ge 1 ] || { echo "Usage : $0 SMC1A [echantillon]" >&2; exit 1; }
if [ $# -ge 2 ]; then sf gene -w "$WORK" --gene "$1" --sample "$2"
else sf gene -w "$WORK" --gene "$1"; fi
