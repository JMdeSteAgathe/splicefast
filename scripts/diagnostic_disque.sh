#!/usr/bin/env bash
# Le disque tient-il la charge ? Lecture continue des BAM pendant N minutes
# (20 par defaut), chaque coupure est signalee.   Usage : $0 [minutes]
source "$(dirname "$0")/_commun.sh"
echo "Pendant le test, dans un autre terminal :  sudo dmesg -w | grep -iE 'usb|reset|I/O error'"
sf bench --run "${RUNS[@]}" -w "$WORK" -t 2 --engine python --stress "${1:-20}"
