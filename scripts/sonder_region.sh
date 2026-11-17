#!/usr/bin/env bash
# Relit une region dans le BAM et la compare au comptage stocke.
# Usage : scripts/sonder_region.sh ECH02-PAX_S61 chr5:37025000-37028000
source "$(dirname "$0")/_commun.sh"
sf probe -w "$WORK" --sample "$1" --region "$2"
