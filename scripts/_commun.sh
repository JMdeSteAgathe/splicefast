# shellcheck shell=bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROG="$(dirname "$HERE")"
# shellcheck source=config.sh
source "$HERE/config.sh"
cd "$PROG"
sf() { python3 -m splicefast "$@"; }
panels() { local p out=(); for p in "${PANELS[@]}"; do [ -f "$p" ] && out+=("$p"); done
           [ ${#out[@]} -gt 0 ] && printf '%s\0' --panel "${out[@]}"; return 0; }
ouvrir() {  # ouvre le rapport si une session graphique est la
  echo "Rapport : $1"
  if [ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ] && command -v xdg-open >/dev/null; then
    xdg-open "$1" >/dev/null 2>&1 &
  fi
}
