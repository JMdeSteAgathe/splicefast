#!/usr/bin/env bash
# Telecharge PanelApp Australie (tous les panels) puis refait le rapport.
source "$(dirname "$0")/_commun.sh"
sf panelapp --update -w "$WORK" "$@" && sf report -w "$WORK"
