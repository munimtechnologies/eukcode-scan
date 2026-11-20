#!/bin/bash
# Fetch and build tRNAscan-SE 2.0.12 (Chan & Lowe, GPL-3.0) and Infernal 1.1.5 (Eddy lab, BSD)
# into ./trna-tools next to eukcode-scan. Neither is bundled with this repository.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
set -euo pipefail
cd "$(dirname "$0")"
PREFIX="$(pwd)/trna-tools"
BUILD=/tmp/trna-build
mkdir -p "$PREFIX/bin" "$BUILD"
LOG="$(pwd)/trna-build.log"
: > "$LOG"

# Infernal: use the precompiled Linux/Intel binaries; build from the bundled source if they do not run.
cd "$BUILD"
curl -fsSL -o infernal.tar.gz http://eddylab.org/infernal/infernal-1.1.5-linux-intel-gcc.tar.gz
tar xzf infernal.tar.gz
IDIR="$(find "$BUILD" -maxdepth 1 -type d -name 'infernal-1.1.5*' | head -n 1 || true)"
CMS="$(find "$IDIR" -maxdepth 3 -type f -name cmsearch | head -n 1 || true)"
if [ -n "$CMS" ] && "$CMS" -h > /dev/null 2>&1; then
  echo "using precompiled Infernal binaries from $(dirname "$CMS")"
  find "$(dirname "$CMS")" -maxdepth 1 -type f -exec cp {} "$PREFIX/bin/" \;
else
  echo "precompiled Infernal did not run; building from source"
  cd "$IDIR"
  ./configure --prefix="$PREFIX" >> "$LOG" 2>&1
  make -j2 >> "$LOG" 2>&1
  make install >> "$LOG" 2>&1
  cd "$BUILD"
fi
"$PREFIX/bin/cmsearch" -h > "$BUILD/cmsearch-h.txt"
head -n 2 "$BUILD/cmsearch-h.txt"

# tRNAscan-SE 2.0.12, installed with the same prefix so its config finds Infernal in $PREFIX/bin.
cd "$BUILD"
curl -fsSL -o trnascan.tar.gz http://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz
tar xzf trnascan.tar.gz
TDIR="$(find "$BUILD" -maxdepth 1 -type d -iname 'trnascan*' | head -n 1 || true)"
cd "$TDIR"
./configure --prefix="$PREFIX" >> "$LOG" 2>&1
make >> "$LOG" 2>&1
make install >> "$LOG" 2>&1
if [ ! -x "$PREFIX/bin/tRNAscan-SE" ]; then
  echo "tRNAscan-SE was not installed"; tail -n 60 "$LOG"; exit 1
fi
grep -i infernal "$PREFIX/bin/tRNAscan-SE.conf" | head -n 5 || true
PATH="$PREFIX/bin:$PATH" "$PREFIX/bin/tRNAscan-SE" -h 2>&1 | head -n 5 || true
rm -rf "$BUILD"
echo "tRNAscan-SE build OK"
