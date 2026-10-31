# 2026-05-06 CP Study Log

## Scope

Investigate Bd -> J/psi K_S asymmetries for the existing SVS and SSD_CP HepMC
samples in `/home/darthlucas/programs/test`.

## Files Created

- `cp_study/hepmc_cp_asymmetry.py`
  - Reason: parse accepted HepMC2 events, reconstruct the signal B chain, compute
    proper time in ps, and plot CP-tag and initial/final flavour asymmetries.
- `cp_study/run_existing_samples.sh`
  - Reason: one-command runner for the existing `eventsSVS.hepmc` and
    `events_SSDCP2.hepmc` samples.
- `cp_study/README.md`
  - Reason: document the observables and why debug output written inside
    `EvtDecayAmp::decay()` is not an accepted-event ntuple.
- `cp_study/results_smoke/*`
  - Reason: smoke-test output using `--max-events 1000`.
- `cp_study/results_existing/*`
  - Reason: full accepted-event analysis for 50k generated events per sample.

No existing source, decay, generator, or installed EvtGen file was modified in
this pass, so no backup was needed. If a source edit becomes necessary later,
backups should be made before touching the file.

## Commands Run

```bash
cd /home/darthlucas/programs/test
bash cp_study/run_existing_samples.sh --max-events 1000 --outdir cp_study/results_smoke
bash cp_study/run_existing_samples.sh
python3 -m py_compile cp_study/hepmc_cp_asymmetry.py
```

## Results

Full samples:

- SVS: 43,362 accepted Bd -> J/psi K_S candidates from 50,000 events.
- SSD_CP: 43,173 accepted Bd -> J/psi K_S candidates from 50,000 events.

Key observation:

- `cp_tag_asymmetry.png`: SSD_CP is not flat when using the accepted-event CP
  observable `(N(initial B0) - N(initial anti-B0)) / sum`; SVS is consistent
  with zero, as expected for this observable.
- `initial_final_flavor_asymmetry.png`: SVS follows the expected
  same/opposite flavour mixing curve. SSD_CP is approximately flat here. This
  same/opposite observable is not the CP asymmetry for a CP eigenstate.

Interpretation:

The flat SSD_CP curve seen in same/opposite flavour plots is an observable
choice issue, not evidence by itself that SSD_CP failed to generate the CP
asymmetry. Also, `B -> B -> X` chains are used in the HepMC record to represent
neutral-B flavour changes; the parser walks back to the root B before computing
the initial flavour and lifetime.

## Generated Outputs

- `cp_study/results_existing/cp_tag_asymmetry.png`
- `cp_study/results_existing/initial_final_flavor_asymmetry.png`
- `cp_study/results_existing/proper_time_distribution.png`
- `cp_study/results_existing/asymmetry_bins.csv`
- `cp_study/results_existing/SVS_signal_events.csv`
- `cp_study/results_existing/SSD_CP_signal_events.csv`
- `cp_study/results_existing/summary.json`

## EvtGen Source Status

No SSD_CP model file was edited. The current accepted HepMC analysis does not
trigger the requested "SSD_CP is flat" condition for the CP-tag observable.
