#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

python3 cp_study/hepmc_cp_asymmetry.py \
  --sample SVS=eventsSVS.hepmc \
  --sample SSD_CP=events_SSDCP2.hepmc \
  --outdir cp_study/results_existing \
  "$@"
