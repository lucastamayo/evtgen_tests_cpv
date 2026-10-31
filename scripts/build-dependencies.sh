#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/env.sh"
work="$REPO_DIR/build/dependencies"
mkdir -p "$work"
for prefix in "$HEPMC_PREFIX" "$PYTHIA_PREFIX" "$PHOTOS_PREFIX" "$TAUOLA_PREFIX"; do
  if [[ -e "$prefix" ]]; then echo "Prefix exists; choose an empty DEPS_PREFIX: $prefix" >&2; exit 2; fi
done
command -v gfortran >/dev/null
for archive in hepmc2.06.10.tgz pythia8306.tgz PHOTOS.3.64.tar.gz TAUOLA.1.1.8-LHC.tar.gz; do
  tar -xf "$REPO_DIR/archives/$archive" -C "$work"
done
cmake -S "$work/HepMC-2.06.10" -B "$work/hepmc-build" \
  -DCMAKE_INSTALL_PREFIX="$HEPMC_PREFIX" -Dmomentum=GEV -Dlength=MM
cmake --build "$work/hepmc-build" -j "${JOBS:-2}"
cmake --install "$work/hepmc-build"
(cd "$work/pythia8306"; ./configure --enable-shared --prefix="$PYTHIA_PREFIX"; make -j "${JOBS:-2}"; make install)
(cd "$work/PHOTOS"; ./configure --without-hepmc3 --with-hepmc="$HEPMC_PREFIX" --prefix="$PHOTOS_PREFIX"; make -j "${JOBS:-2}"; make install)
(cd "$work/TAUOLA"; ./configure --without-hepmc3 --with-hepmc="$HEPMC_PREFIX" --prefix="$TAUOLA_PREFIX"; make -j "${JOBS:-2}"; make install)
