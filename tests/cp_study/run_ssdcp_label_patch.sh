#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

events="${1:-50000}"
outdir="cp_study/results_ssdcp_label_patch"
evtgen_lib="/home/darthlucas/programs/test/cp_study/evtgen_patch/build/lib"

mkdir -p "${outdir}"

LD_LIBRARY_PATH="${evtgen_lib}:${LD_LIBRARY_PATH:-}" \
  cp_study/generator_ssdcp_label_patch "${events}" \
  > "${outdir}/generation_${events}.log" 2>&1

python3 cp_study/hepmc_cp_asymmetry.py \
  --sample SVS=eventsSVS.hepmc \
  --sample SSD_CP_LABEL_PATCH="${outdir}/events_SSDCP_label_patch.hepmc" \
  --outdir "${outdir}/analysis_${events}" \
  --max-events "${events}"

python3 cp_study/ssdcp_truth_asymmetry.py \
  "${outdir}/ssdcp_truth_tags.csv" \
  --outdir "${outdir}/analysis_${events}"
