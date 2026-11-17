#!/usr/bin/env bash
# Compile splicefast/bin/sfscan.exe pour Windows depuis Linux (mingw-w64),
# binaire statique : aucune DLL a installer sur le PC Windows.
#   sudo apt install g++-mingw-w64-x86-64-posix libz-mingw-w64-dev
#   scripts/build_sfscan_windows.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
B="$ROOT/third_party/win"; mkdir -p "$B/stub"
cp "$ROOT/src/mingw_regex_stub.h" "$B/stub/regex.h"   # htslib : regex non utilise par sfscan
cd "$B"
[ -d htslib-1.21 ] || { curl -fL -o htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2; tar xjf htslib.tar.bz2; }
cd htslib-1.21
./configure --host=x86_64-w64-mingw32 CC=x86_64-w64-mingw32-gcc-posix \
  --disable-bz2 --disable-lzma --disable-libcurl --disable-gcs --disable-s3 \
  --disable-plugins --without-libdeflate CFLAGS=-O2 >/dev/null
make -j"$(nproc)" lib-static CPPFLAGS="-I$B/stub" >/dev/null
x86_64-w64-mingw32-g++-posix -O2 -std=c++17 -I. "$ROOT/src/sfscan.cpp" \
  -o "$ROOT/splicefast/bin/sfscan.exe" libhts.a -static -static-libgcc \
  -static-libstdc++ -lz -lws2_32 -lpthread
x86_64-w64-mingw32-strip "$ROOT/splicefast/bin/sfscan.exe"
echo "OK : $ROOT/splicefast/bin/sfscan.exe"
