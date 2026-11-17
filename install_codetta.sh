#!/bin/bash
# Fetch and build Codetta 2.0 (Shulgina & Eddy, GPL-3.0) next to eukcode-scan. Not bundled.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
set -euo pipefail
cd "$(dirname "$0")"
if [ ! -d codetta ]; then
  git clone --depth 1 https://github.com/kshulgina/codetta.git codetta
fi
cd codetta
git rev-parse HEAD
bash setup.sh > ../codetta-setup.log 2>&1 || true
tail -n 5 ../codetta-setup.log
# setup.sh exits 0 even when its requirements check (pkg_resources) or a build step fails,
# so if the binaries are missing, build the bundled HMMER the same way it would have.
if [ ! -x hmmer-3.3.2/bin/hmmscan ]; then
  echo "setup.sh did not build HMMER; building it directly"
  tar xf hmmer.tar.gz
  cd hmmer-3.3.2
  ./configure --prefix="$(pwd)" > ../../hmmer-build.log 2>&1
  make -j4 >> ../../hmmer-build.log 2>&1
  make install >> ../../hmmer-build.log 2>&1
  (cd easel && make install >> ../../../hmmer-build.log 2>&1)
  cd ..
fi
for b in hmmscan hmmpress esl-sfetch; do
  if [ ! -x "hmmer-3.3.2/bin/$b" ]; then
    echo "missing hmmer-3.3.2/bin/$b"; tail -n 60 ../codetta-setup.log; exit 1
  fi
done
echo "Codetta build OK"
