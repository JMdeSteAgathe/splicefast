#!/usr/bin/env bash
# Verifie Python/numpy/scipy et le lecteur C++, puis lance l'auto-test.
source "$(dirname "$0")/_commun.sh"
python3 -c "import sys; print('Python', sys.version.split()[0])"
python3 -c "import numpy, scipy; print('numpy', numpy.__version__, '/ scipy', scipy.__version__)" || {
  echo "numpy ou scipy manquant :  python3 -m pip install --user numpy scipy" >&2; exit 1; }
if [ -x splicefast/bin/sfscan ]; then splicefast/bin/sfscan --version
else echo "sfscan non compile : lancer 'make' (ou ./install_htslib.sh sans droits root)"; fi
sf selftest -w "${TMPDIR:-/tmp}/splicefast_selftest"
echo 'Attendu : "4/4 evenements retrouves" et "lecteur C++ vs Python ... IDENTIQUES".'
