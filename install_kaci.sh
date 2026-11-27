#!/usr/bin/env bash
# eukcode-scan addendum A8: build lib-kaci (MIT, Artem V. Melnykov) at a pinned commit.
# Copyright 2026 MV2 of Munim, Inc. Apache-2.0. Code written by AI agents working in MV2's name.
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
COMMIT=cc8601f2121ddb2546b69c5d494652c38294c02a
cd "$HERE"
if [ ! -d lib-kaci ]; then
  git clone --quiet https://github.com/artmeln/lib-kaci.git lib-kaci
fi
git -C lib-kaci checkout --quiet "$COMMIT"
# The python image (buildpack-deps) ships gcc, make and zlib headers; install only if missing.
if ! command -v make >/dev/null 2>&1 || [ ! -f /usr/include/zlib.h ]; then
  apt-get update > /tmp/apt.log 2>&1
  apt-get install -y build-essential zlib1g-dev >> /tmp/apt.log 2>&1
fi
(cd lib-kaci/src && make libkaci > /tmp/kaci_make.log 2>&1) || { tail -n 40 /tmp/kaci_make.log; exit 1; }
test -f lib-kaci/src/libkaci.so
echo "lib-kaci $COMMIT built: $(ls -l lib-kaci/src/libkaci.so)"
