# SpliceFast

Détection d'événements d'épissage aberrants à partir de RNA-seq (BAM), dans
l'esprit de FRASER / FRASER 2.0 (psi5, psi3, theta, Intron Jaccard Index).
L'outil vise de gros runs (70 à 80 échantillons, plusieurs centaines de Go de
BAM) sur une machine modeste, sans droits administrateur, sous Linux ou Windows.

Documentation complète : **[README_FR.md](README_FR.md)**. Installation sous
Windows : [windows/INSTALL_WINDOWS.md](windows/INSTALL_WINDOWS.md). Nouveautés :
`CHANGELOG_2.*.md`.

## Contenu du dépôt

| Dossier / fichier | Rôle |
|---|---|
| `splicefast/` | le programme (paquet Python), dont le visualiseur `viewer.html` |
| `splicefast/bin/sfscan.exe` | lecteur BAM compilé pour Windows (aucune DLL requise) |
| `src/sfscan.cpp`, `Makefile`, `install_htslib.sh` | lecteur BAM C++ : `make` sous Linux |
| `scripts/*.sh` | lanceurs Linux (`config.sh` à éditer) |
| `windows/*.bat` | lanceurs Windows (`config.bat` à éditer) |
| `tests/` | tests de non-régression : `python3 tests/test_recurrence.py`, etc. |
| `scripts/make_dist.sh` | fabrique les archives `dist/splicefast-<version>-linux.zip` et `-windows.zip` |

## Démarrage rapide (Linux)

```bash
make                                   # ou ./install_htslib.sh && make, sans droits root
python3 -m splicefast selftest -w /tmp/sf_selftest
```

Le binaire Linux `splicefast/bin/sfscan` est compilé sur place et n'est pas versionné.

## Distribuer une version

```bash
scripts/make_dist.sh                   # -> dist/splicefast-X_Y_Z-{linux,windows}.zip
```

L'archive Windows s'installe comme les précédentes : les `.bat` à la racine du
dossier, le paquet dans `splicefast\`.
