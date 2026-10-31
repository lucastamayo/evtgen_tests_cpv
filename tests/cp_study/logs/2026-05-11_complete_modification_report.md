# 2026-05-11 - Complete Modification Report

This document summarizes the changes made during the `Bd -> J/psi K_S` CP/mixing investigation. The work was kept under `programs/test/cp_study` whenever possible, and the original installed EvtGen at `/home/darthlucas/local/evtgen200` was not replaced.

## High-Level Goal

The study started with two questions:

1. Compare the existing `SVS` and `SSD_CP` HepMC samples.
2. Understand why `SSD_CP` produced significant same-flavour `B -> B -> X` chains, classified as `weird`, while `SVS` did not.

The final state of the study has two parallel tracks:

1. Analysis of the original HepMC samples.
2. A separate experimental patched generation path for `SSD_CP`, with clean exported B labels and a separate truth-tag CSV.

## Main Conclusions

For the original samples:

```text
SVS:
  unmixed = 35210
  mixed   = 8152
  weird   = 0

SSD_CP:
  unmixed = 35050
  mixed   = 6888
  weird   = 1235
```

The `weird` chains in the original `SSD_CP` sample come from the way `EvtSSDCP` uses `EvtCPUtil::OtherIncoherentB()`. That function can flip B flavour IDs inside the EvtGen tree. When the Pythia8 EvtGen interface writes the tree to the Pythia/HepMC record, this can appear as same-flavour `B -> B -> X` chains.

For the patched `SSD_CP` sample:

```text
SSD_CP_LABEL_PATCH:
  signals = 43426
  initial_B0_fraction = 0.5018
  unmixed = 35085
  mixed   = 8341
  weird   = 0
```

The patched generation removes the `weird` contamination from the exported HepMC B labels. The internal `SSD_CP` truth tag is exported separately to:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

The truth-tag asymmetry closes against the expected reference:

```text
RMS vs 0.700 sin(delta_m t) = 0.0322
```

## Existing Analysis Code Changes

### `cp_study/hepmc_cp_asymmetry.py`

Purpose: parse accepted HepMC2 events and make CP/mixing plots from the final event record.

Main changes:

- Added `mplhep` CMS styling for all plots.
- Added a direct HepMC2 text parser, independent of `pyhepmc_ng`.
- Added `SignalEvent` records with:
  - event number
  - proper time
  - root B flavour
  - final/signal B flavour
  - root barcode
  - signal barcode
  - chain classification
- Added B-chain classification:

```text
unmixed: root B is also the signal B
mixed:   root B and signal B have opposite flavour in a B -> B -> X chain
weird:   root B and signal B have the same flavour in a B -> B -> X chain
```

- Added output CSVs:
  - `SVS_signal_events.csv`
  - `SSD_CP_signal_events.csv`
  - `asymmetry_bins.csv`
  - `mixing_counts.csv`
  - `mixing_closure.csv`
  - `summary.json`
- Added plots:
  - `cp_tag_asymmetry.png`
  - `initial_final_flavor_asymmetry.png`
  - `record_mixing_asymmetry.png`
  - `mixing_components.png`
  - `proper_time_distribution.png`
- Added `record_mixing` observable:

```text
mixed-like = mixed + weird
(N(unmixed) - N(mixed-like)) / sum
```

- Added closure metrics against `cos(delta_m t)`:

```text
mixing_closure.csv
```

Important interpretation:

- `initial_final_flavor_asymmetry.png` works as a mixing diagnostic for `SVS`.
- It does not work as a clean mixing observable for original `SSD_CP`.
- `record_mixing_asymmetry.png` fixes the integrated `mixed + weird` counting, but original `SSD_CP` still does not close against `cos(delta_m t)` bin-by-bin.

### `mixing_from_hepmc.py`

Purpose: standalone mixed/unmixed/weird time histogram from HepMC.

Main changes:

- Added/kept CMS-style plotting with `mplhep` when available.
- Added `--max` / `--max-events` support.
- Added a fallback to `cp_study.hepmc_cp_asymmetry.read_hepmc` when `pyhepmc_ng` is unavailable.
- Avoided double counting by selecting root B candidates.
- Plots separated time distributions for:
  - sum
  - mixed
  - unmixed
  - weird

This script was used to cross-check that `SSD_CP` had a nonzero `weird` component while `SVS` did not.

### `cp_study/fwlite_cp_asymmetry.py`

Purpose: FWLite-based CP asymmetry plotting.

Main changes:

- Updated plot style to CMS/mplhep where possible.
- Made FWLite imports optional so the script can be imported/checked in environments where CMSSW/FWLite is not available.

No FWLite production output was overwritten as part of the later EvtGen-label patch work.

### `cp_study/run_existing_samples.sh`

Purpose: convenience wrapper for analyzing existing HepMC samples.

Current behavior:

```bash
python3 cp_study/hepmc_cp_asymmetry.py \
  --sample SVS=eventsSVS.hepmc \
  --sample SSD_CP=events_SSDCP2.hepmc \
  --outdir cp_study/results_existing \
  "$@"
```

This allows commands such as:

```bash
bash cp_study/run_existing_samples.sh --max-events 1000000 --outdir cp_study/results_1M_mixing
```

## Documentation Files Added or Updated

### `cp_study/README.md`

Updated with:

- How to run existing-sample analysis.
- Meaning of each plot.
- Why `initial_final_flavor_asymmetry.png` is not reliable for original `SSD_CP`.
- Meaning of `record_mixing_asymmetry.png`.
- Description of the experimental `SSD_CP` label patch.
- How to run the patched generation.
- Meaning of the new truth-tag CSV.

### `cp_study/logs/2026-05-11_mixing_chain_log.md`

Records the mixing-chain study and summary of mixed/unmixed/weird counting.

### `cp_study/logs/2026-05-11_weird_chain_investigation.md`

Records the source-level investigation of `weird` chains:

- `SVS` model does not call `EvtCPUtil`.
- `SSD_CP` calls `EvtCPUtil::OtherIncoherentB()`.
- `OtherIncoherentB()` can flip flavour IDs.
- The Pythia8 EvtGen interface exports the mutated tree to HepMC.

Later extended with:

- The label-cleaning patch.
- The truth-tag export.
- Numerical before/after comparisons.

### `cp_study/logs/2026-05-11_complete_modification_report.md`

This file.

## Generated Analysis Directories

### `cp_study/results_1M_mixing`

Generated from the original existing samples with:

```bash
bash cp_study/run_existing_samples.sh --max-events 1000000 --outdir cp_study/results_1M_mixing
```

Important files:

- `summary.json`
- `mixing_counts.csv`
- `mixing_closure.csv`
- `cp_tag_asymmetry.png`
- `initial_final_flavor_asymmetry.png`
- `record_mixing_asymmetry.png`
- `mixing_components.png`
- `proper_time_distribution.png`
- `SVS_signal_events.csv`
- `SSD_CP_signal_events.csv`

Important numerical result:

```text
SVS:
  signals = 43362
  mixed_chain_fraction = 0.1879987
  weird = 0

SSD_CP:
  signals = 43173
  mixed_chain_fraction = 0.1881500
  mixed = 6888
  weird = 1235
```

Closure against `cos(delta_m t)`:

```text
SVS initial_final:       RMS = 0.0263
SVS record_mixing:       RMS = 0.0263
SSD_CP initial_final:    RMS = 1.1306
SSD_CP record_mixing:    RMS = 1.0938
```

Interpretation:

- `SVS` closes against the expected mixing curve.
- Original `SSD_CP` does not, even after counting `weird` as mixed-like.

### `cp_study/results_smoke_record_mixing`

Small validation output used while developing `record_mixing` and `mixing_closure`.

### `cp_study/results_mixing0_test`

Temporary test generated with:

```text
cp_study/generator_ssdcp_mixing0_test.cc
```

Important caveat:

- The EvtGen constructor argument `0` is not "no mixing".
- In this code path, `0` maps to `EvtCPUtil::Coherent`.

Result:

```text
SSD_CP_MIX0:
  signals = 8553
  unmixed = 6930
  mixed   = 1390
  weird   = 233
```

Conclusion:

- This was not a clean no-mixing test.
- It showed that simply changing the constructor integer did not remove the same-flavour chains.

## Experimental SSD_CP Label Patch

This was done separately from the installed EvtGen.

### Local EvtGen Copy

Created:

```text
cp_study/evtgen_patch/R02-00-00
```

This is a copy of:

```text
/home/darthlucas/programs/EvtGen/R02-00-00
```

Built under:

```text
cp_study/evtgen_patch/build
```

Important built libraries:

```text
cp_study/evtgen_patch/build/lib/libEvtGen.so.2
cp_study/evtgen_patch/build/lib/libEvtGenExternal.so.2
```

The original installed libraries under `/home/darthlucas/local/evtgen200` were not replaced.

### `cp_study/evtgen_patch/R02-00-00/src/EvtGenBase/EvtCPUtil.cpp`

Purpose of patch:

1. Do not re-sample the parent lifetime when the B already has a generic incoherent-mixing parent.
2. Choose the internal `other_b` flavour tag for `SSD_CP` without rewriting the B IDs in the exported decay tree.

Original problematic behavior:

```text
if B was already mixed:
  p->getParent()->setLifetime()
  p->getParent()->setId(chargeConj(...))
  p->setId(chargeConj(...))
```

Problem:

- Re-sampling lifetime breaks the correlation between the mixing decision and the exported time.
- Rewriting `p` and `p->getParent()` contaminates the HepMC B labels.

Patched behavior:

```text
if B was already mixed:
  use p->getParent()->getLifetime()

choose effectiveTagSource with random flip
otherb = chargeConj(effectiveTagSource)
do not modify p->getId()
do not modify p->getParent()->getId()
```

Effect:

- The internal `SSD_CP` amplitude still receives `other_b`.
- The exported HepMC B labels remain physically clean.

### `cp_study/evtgen_patch/R02-00-00/src/EvtGenModels/EvtSSDCP.cpp`

Purpose of patch:

- Export the accepted `SSD_CP` truth tag through `EvtParticle` attributes.

Added attributes:

```text
ssdcp_truth
ssdcp_parent_pdg
ssdcp_other_b_pdg
ssdcp_flip
ssdcp_t_mm
```

These are set after:

```cpp
EvtCPUtil::getInstance()->OtherIncoherentB(p, t, other_b, 0.5);
```

Important detail:

- Attributes are stored on the accepted EvtGen particle tree.
- This avoids relying on the old `ssdcp_debug.txt`, which can include accept/reject trials and is not a clean accepted-event ntuple.

### `cp_study/EvtGenCleanLabels.h`

Purpose:

- Local copy of the Pythia8 EvtGen interface.
- Used only by the patched generator.

Changes:

- Added a public struct:

```cpp
struct SSDCPTruth {
  int pyIndex, parentPdg, otherBPdg, flip;
  double tMm;
};
```

- Added:

```cpp
vector<SSDCPTruth> ssdcpTruth;
```

- Cleared `ssdcpTruth` at the start of each `EvtGenDecays::decay()` call.
- Copied `ssdcp_*` attributes from the accepted EvtGen particle to `ssdcpTruth`.
- Synchronized B labels in the Pythia record with the EvtGen mother ID for neutral B mesons.

This header is not installed globally. It is included only by:

```text
cp_study/generator_ssdcp_label_patch.cc
```

### `cp_study/generator_ssdcp_label_patch.cc`

Purpose:

- Generate a new `SSD_CP` sample with:
  - clean exported HepMC B labels
  - a separate accepted truth-tag CSV

Outputs:

```text
cp_study/results_ssdcp_label_patch/events_SSDCP_label_patch.hepmc
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

CLI:

```bash
cp_study/generator_ssdcp_label_patch 50000
```

The argument is the number of accepted Pythia events to attempt/write.

### `cp_study/run_ssdcp_label_patch.sh`

Purpose:

- Reproduce the patched generation and analysis in one command.

Usage:

```bash
bash cp_study/run_ssdcp_label_patch.sh 50000
```

It does three things:

1. Runs `generator_ssdcp_label_patch` using the patched local EvtGen library via `LD_LIBRARY_PATH`.
2. Runs `hepmc_cp_asymmetry.py` on:
   - existing `SVS`
   - patched `SSD_CP`
3. Runs `ssdcp_truth_asymmetry.py` on the truth-tag CSV.

### `cp_study/ssdcp_truth_asymmetry.py`

Purpose:

- Analyze the accepted `SSD_CP` truth tag exported to CSV.

Input:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

Outputs:

```text
ssdcp_truth_tag_asymmetry.png
ssdcp_truth_asymmetry_bins.csv
```

Observable:

```text
signal_initial = -other_b
(N(truth initial B0) - N(truth initial anti-B0)) / sum
```

Reference curve:

```text
0.700 sin(delta_m t)
```

For the 50k patched sample:

```text
truth rows = 43426
events with truth = 33942
RMS vs 0.700 sin(delta_m t) = 0.0322
```

## Patched Output Directory

### `cp_study/results_ssdcp_label_patch`

Main files:

```text
events_SSDCP_label_patch.hepmc
ssdcp_truth_tags.csv
generation_5000.log
generation_50000.log
analysis_5000/
analysis_50000/
```

### `cp_study/results_ssdcp_label_patch/analysis_50000`

Important files:

```text
summary.json
mixing_counts.csv
mixing_closure.csv
cp_tag_asymmetry.png
initial_final_flavor_asymmetry.png
record_mixing_asymmetry.png
mixing_components.png
proper_time_distribution.png
ssdcp_truth_tag_asymmetry.png
ssdcp_truth_asymmetry_bins.csv
```

Important result:

```text
SVS:
  signals = 43362
  initial_B0_fraction = 0.4985
  unmixed = 35210
  mixed   = 8152
  weird   = 0

SSD_CP_LABEL_PATCH:
  signals = 43426
  initial_B0_fraction = 0.5018
  unmixed = 35085
  mixed   = 8341
  weird   = 0
```

Interpretation:

- The exported HepMC labels are clean.
- Same-flavour `B -> B -> X` chains are removed.
- The internal `SSD_CP` CP tag is now in `ssdcp_truth_tags.csv`, not in the exported B PDG labels.

## Important Interpretation Changes

### Before Patch

Original `SSD_CP` used the exported HepMC B labels for two different roles:

1. Physical B-chain/mixing labels.
2. Internal CP-tag information from `EvtCPUtil`.

Those roles got mixed together and produced:

```text
B -> B -> X
```

same-flavour chains.

### After Patch

The two concepts are separated:

1. HepMC B labels describe the exported B chain.
2. `ssdcp_truth_tags.csv` stores the internal `SSD_CP` truth tag.

That is why:

```text
weird = 0
```

and:

```text
ssdcp_truth_tag_asymmetry.png
```

is now the correct place to inspect the `SSD_CP` truth-tag CP asymmetry.

## Commands Used to Reproduce Main Outputs

Original-sample analysis:

```bash
cd /home/darthlucas/programs/test
bash cp_study/run_existing_samples.sh --max-events 1000000 --outdir cp_study/results_1M_mixing
```

Patched `SSD_CP` generation and analysis:

```bash
cd /home/darthlucas/programs/test
bash cp_study/run_ssdcp_label_patch.sh 50000
```

Standalone truth-tag analysis:

```bash
python3 cp_study/ssdcp_truth_asymmetry.py \
  cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv \
  --outdir cp_study/results_ssdcp_label_patch/analysis_50000
```

## Build Notes

The patched EvtGen was configured with CMake under:

```text
cp_study/evtgen_patch/build
```

It links against the local installations:

```text
/home/darthlucas/local/pythia8306
/home/darthlucas/local/hepmc2.06.10
/home/darthlucas/local/photos364
/home/darthlucas/local/tauola118
```

The patched generator must be run with the patched library first in `LD_LIBRARY_PATH`:

```bash
LD_LIBRARY_PATH=/home/darthlucas/programs/test/cp_study/evtgen_patch/build/lib:$LD_LIBRARY_PATH
```

`run_ssdcp_label_patch.sh` handles this automatically.

## Side Effects and Temporary Files

### `ssdcp_debug.txt`

The original/local `EvtSSDCP.cpp` contains:

```cpp
static std::ofstream fout("ssdcp_debug.txt");
```

So running `SSD_CP` generation can create or overwrite:

```text
/home/darthlucas/programs/test/ssdcp_debug.txt
```

This file is not used as the final truth source because it can contain accept/reject trial information. The accepted-event truth source is now:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

### Python Cache Files

Running `py_compile` and Python scripts may create:

```text
cp_study/__pycache__/
```

These are incidental Python cache files.

### Temporary Mixing-0 Test

The following files were created for the constructor-argument test:

```text
cp_study/generator_ssdcp_mixing0_test.cc
cp_study/generator_ssdcp_mixing0_test
cp_study/results_mixing0_test/
```

They are kept for provenance. The result should be interpreted with the caveat that constructor argument `0` is `Coherent`, not "mixing off".

## Files Added by This Work

Primary source/scripts:

```text
cp_study/EvtGenCleanLabels.h
cp_study/generator_ssdcp_label_patch.cc
cp_study/run_ssdcp_label_patch.sh
cp_study/ssdcp_truth_asymmetry.py
cp_study/generator_ssdcp_mixing0_test.cc
```

Primary logs/docs:

```text
cp_study/logs/2026-05-11_mixing_chain_log.md
cp_study/logs/2026-05-11_weird_chain_investigation.md
cp_study/logs/2026-05-11_complete_modification_report.md
```

Primary generated output directories:

```text
cp_study/results_1M_mixing/
cp_study/results_smoke_record_mixing/
cp_study/results_mixing0_test/
cp_study/results_ssdcp_label_patch/
```

Local EvtGen copy/build:

```text
cp_study/evtgen_patch/R02-00-00/
cp_study/evtgen_patch/build/
```

## Files Modified by This Work

```text
cp_study/README.md
cp_study/hepmc_cp_asymmetry.py
cp_study/fwlite_cp_asymmetry.py
mixing_from_hepmc.py
```

Inside the local copied EvtGen patch:

```text
cp_study/evtgen_patch/R02-00-00/src/EvtGenBase/EvtCPUtil.cpp
cp_study/evtgen_patch/R02-00-00/src/EvtGenModels/EvtSSDCP.cpp
```

No installed EvtGen file under `/home/darthlucas/local/evtgen200` was modified.

## Recommended Usage Going Forward

For original samples:

- Use `cp_tag_asymmetry.png` for the original `SSD_CP` CP behavior.
- Treat `initial_final_flavor_asymmetry.png` as unreliable for original `SSD_CP`.
- Use `mixing_components.png`, `mixing_counts.csv`, and `mixing_closure.csv` as diagnostics.

For patched samples:

- Use HepMC labels for clean B-chain inspection.
- Use `ssdcp_truth_tags.csv` and `ssdcp_truth_tag_asymmetry.png` for the internal `SSD_CP` CP truth tag.
- Do not infer the `SSD_CP` internal CP tag from exported B PDG labels after the patch.

The patched setup is the cleaner design because it separates:

```text
physical exported event record labels
```

from:

```text
internal CP truth tag used by SSD_CP
```
