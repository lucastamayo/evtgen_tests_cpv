#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/env.sh"
variant="${1:-2.0.0}"; events="${2:-1000}"; seed="${3:-12345}"
case "$variant" in 2.0.0|2.2.3|2.0.0-label-patch) ;; *) exit 2;; esac
out="$REPO_DIR/runs/${variant}_${events}_${seed}"
if [[ -e "$out" ]]; then echo "Run already exists: $out" >&2; exit 2; fi
mkdir -p "$out"
export LD_LIBRARY_PATH="$REPO_DIR/build/$variant/evtgen/lib:$LD_LIBRARY_PATH"
for model in SVS SSD_CP; do
  (cd "$out"; "$REPO_DIR/build/$variant/generate" "$model" "$events" "$out/$model.hepmc" "$out/${model}_truth.csv" "$seed" 1 "$REPO_DIR/tests") > "$out/$model.log" 2>&1
done
python3 "$REPO_DIR/tests/cp_study/hepmc_cp_asymmetry.py" \
  --sample "SVS=$out/SVS.hepmc" --sample "SSD_CP=$out/SSD_CP.hepmc" --outdir "$out/analysis"
if [[ "$variant" == 2.0.0-label-patch ]]; then
  python3 "$REPO_DIR/tests/cp_study/ssdcp_truth_asymmetry.py" "$out/SSD_CP_truth.csv" --outdir "$out/truth_analysis"
fi
printf 'Results: %s\n' "$out"
