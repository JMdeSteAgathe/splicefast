#!/usr/bin/env bash
# Fabrique les deux archives de distribution depuis le depot :
#   dist/splicefast-<version>-linux.zip    (sources, scripts .sh, Makefile, tests)
#   dist/splicefast-<version>-windows.zip  (.bat a la racine + paquet splicefast\ avec sfscan.exe)
# Usage : scripts/make_dist.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
V=$(sed -n 's/^__version__ = "\(.*\)"/\1/p' splicefast/__init__.py)
TAG=${V//./_}
OUT="$ROOT/dist"; TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"

# fichiers suivis par git uniquement (pas de __pycache__, third_party, resultats...)
files() { git ls-files "$@"; }

# --- Linux : un dossier splicefast/ comme l'archive d'origine
L="$TMP/linux/splicefast"; mkdir -p "$L"
files | grep -v -E '^(windows/|dist/|README\.md$|splicefast/bin/sfscan\.exe$)' | while read -r f; do
  mkdir -p "$L/$(dirname "$f")"; cp -p "$f" "$L/$f"; done
# fichiers locaux non versionnes (vrais identifiants) : gardes dans l'archive privee
for f in scripts/*.local.tsv; do [ -f "$f" ] && cp -p "$f" "$L/$f"; done
(cd "$TMP/linux" && rm -f "$OUT/splicefast-$TAG-linux.zip" && zip -qr "$OUT/splicefast-$TAG-linux.zip" splicefast)

# --- Windows : .bat et docs a la racine, paquet Python dans splicefast\
W="$TMP/windows"; mkdir -p "$W"
cp -p windows/* "$W/"
cp -p scripts/panelapp_download.ps1 README_FR.md CHANGELOG_2.6.md CHANGELOG_2.7.md "$W/"
files splicefast | grep -v -E '^splicefast/bin/sfscan$' | while read -r f; do
  mkdir -p "$W/$(dirname "$f")"; cp -p "$f" "$W/$f"; done
[ -f "$W/splicefast/bin/sfscan.exe" ] || echo "ATTENTION : splicefast/bin/sfscan.exe absent (scripts/build_sfscan_windows.sh)" >&2
(cd "$W" && rm -f "$OUT/splicefast-$TAG-windows.zip" && zip -qr "$OUT/splicefast-$TAG-windows.zip" .)

ls -l "$OUT"/splicefast-"$TAG"-*.zip
