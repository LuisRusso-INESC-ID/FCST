#!/usr/bin/env bash
# Fetch and build SDSL-lite and CST-CN into deps/ (used by `make bench_cst`).
set -euo pipefail
cd "$(dirname "$0")/.."
DEPS="$PWD/deps"
mkdir -p "$DEPS"
cd "$DEPS"

if [ ! -d sdsl-lite ]; then
    git clone https://github.com/simongog/sdsl-lite.git
    (cd sdsl-lite && git checkout c32874c && git submodule update --init external/libdivsufsort)
fi
if [ ! -f sdsl/lib/libsdsl.a ]; then
    mkdir -p sdsl-lite/build
    (cd sdsl-lite/build &&
     cmake .. -DCMAKE_INSTALL_PREFIX="$DEPS/sdsl" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-Wno-error -Wno-template-body" &&
     make -j"$(nproc)" sdsl divsufsort divsufsort64 && make install)
fi
if [ ! -d CST-CN ]; then
    git clone https://github.com/rcanovas/CST-CN.git
    (cd CST-CN && git checkout b8273c9)
fi
echo "dependencies ready in $DEPS"
