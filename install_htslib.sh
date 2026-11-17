#!/usr/bin/env bash
# =============================================================================
#  Compile htslib (+ libdeflate si cmake est la) dans ./third_party, SANS
#  droits root, puis compile sfscan.
#
#  A utiliser si ni `apt install libhts-dev` ni `conda install htslib` ne sont
#  possibles. Il faut : gcc/g++, make, curl ou wget, et zlib (en-tetes
#  zlib1g-dev ; sinon `conda install zlib`).
#
#    ./install_htslib.sh            # versions par defaut
#    HTSLIB_VERSION=1.21 ./install_htslib.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$PWD
TP=$ROOT/third_party
HTSLIB_VERSION=${HTSLIB_VERSION:-1.21}
LIBDEFLATE_VERSION=${LIBDEFLATE_VERSION:-1.22}
JOBS=$(nproc 2>/dev/null || echo 2)
mkdir -p "$TP" && cd "$TP"

fetch() {  # url fichier
  if [ -s "$2" ]; then return; fi
  if command -v curl >/dev/null; then curl -fL --retry 3 -o "$2.part" "$1"
  else wget -O "$2.part" "$1"; fi
  mv "$2.part" "$2"
}

echo "== verification de zlib"
if ! printf '#include <zlib.h>\nint main(){return 0;}\n' | ${CC:-gcc} -x c - -o /dev/null 2>/dev/null; then
  echo "ERREUR : en-tetes zlib absents (paquet zlib1g-dev, ou conda install zlib)" >&2
  exit 1
fi

CPPFLAGS_EXTRA=""; LDFLAGS_EXTRA=""; DEFLATE_OPT="--without-libdeflate"
if command -v cmake >/dev/null; then
  echo "== libdeflate $LIBDEFLATE_VERSION (decompression ~2x plus rapide)"
  fetch "https://github.com/ebiggers/libdeflate/archive/refs/tags/v$LIBDEFLATE_VERSION.tar.gz" libdeflate.tar.gz
  rm -rf libdeflate-src && mkdir libdeflate-src
  tar xzf libdeflate.tar.gz -C libdeflate-src --strip-components=1
  cmake -S libdeflate-src -B libdeflate-src/build -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$TP/libdeflate" -DCMAKE_INSTALL_LIBDIR=lib \
        -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON >/dev/null
  cmake --build libdeflate-src/build -j "$JOBS" >/dev/null
  cmake --install libdeflate-src/build >/dev/null
  CPPFLAGS_EXTRA="-I$TP/libdeflate/include"; LDFLAGS_EXTRA="-L$TP/libdeflate/lib"
  DEFLATE_OPT="--with-libdeflate"
else
  echo "   (cmake absent : libdeflate sautee -- 'pip install --user cmake' pour l'avoir)"
fi

echo "== htslib $HTSLIB_VERSION"
fetch "https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2" htslib.tar.bz2
rm -rf htslib && mkdir htslib
tar xjf htslib.tar.bz2 -C htslib --strip-components=1
cd htslib
# bz2 / lzma / curl ne servent qu'a certains CRAM ou aux fichiers distants :
# actives seulement s'ils sont presents sur la machine.
has() { printf "#include <$1>\nint main(){return 0;}\n" | ${CC:-gcc} -x c - -o /dev/null 2>/dev/null; }
OPTS="$DEFLATE_OPT"
has bzlib.h || OPTS="$OPTS --disable-bz2"
has lzma.h  || OPTS="$OPTS --disable-lzma"
has curl/curl.h || OPTS="$OPTS --disable-libcurl"
CPPFLAGS="$CPPFLAGS_EXTRA" LDFLAGS="$LDFLAGS_EXTRA" ./configure $OPTS >/dev/null
make -j "$JOBS" lib-static htslib_static.mk >/dev/null
cd "$ROOT"

echo "== sfscan"
make HTSLIB="$TP/htslib"
"$ROOT/splicefast/bin/sfscan" --version
echo
echo "Termine. Verification complete : python3 -m splicefast selftest -w /tmp/splicefast_selftest"
