#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/env.sh"
variant="${1:-2.0.0}"
case "$variant" in 2.0.0|2.2.3|2.0.0-label-patch) ;; *) echo 'Choose 2.0.0, 2.2.3, or 2.0.0-label-patch' >&2; exit 2;; esac
build="$REPO_DIR/build/$variant"
cmake -S "$REPO_DIR/vendor/evtgen-$variant" -B "$build/evtgen" \
  -DCMAKE_BUILD_TYPE=Release -DEVTGEN_HEPMC3=OFF \
  -DEVTGEN_PYTHIA=ON -DEVTGEN_PHOTOS=ON -DEVTGEN_TAUOLA=ON \
  -DHEPMC2_ROOT_DIR="$HEPMC_PREFIX" -DPYTHIA8_ROOT_DIR="$PYTHIA_PREFIX" \
  -DPHOTOSPP_ROOT_DIR="$PHOTOS_PREFIX" -DTAUOLAPP_ROOT_DIR="$TAUOLA_PREFIX"
cmake --build "$build/evtgen" -j "${JOBS:-2}"
flags=()
if [[ "$variant" == 2.0.0-label-patch ]]; then flags+=(-DLABEL_PATCH); fi
"${CXX:-g++}" -std=c++11 -O2 "${flags[@]}" "$REPO_DIR/src/generate.cc" \
  -I"$REPO_DIR/tests/cp_study" -I"$REPO_DIR/vendor/evtgen-$variant" \
  -I"$PYTHIA_PREFIX/include" -I"$HEPMC_PREFIX/include" \
  -L"$build/evtgen/lib" -L"$PYTHIA_PREFIX/lib" -L"$HEPMC_PREFIX/lib" \
  -Wl,-rpath,"$build/evtgen/lib" -Wl,-rpath,"$PYTHIA_PREFIX/lib" -Wl,-rpath,"$HEPMC_PREFIX/lib" \
  -lEvtGenExternal -lEvtGen -lpythia8 -lHepMC -ldl -o "$build/generate"
