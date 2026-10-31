# Bd -> J/psi K_S CP Study

This directory keeps the asymmetry study separate from the existing generators.

## Existing Samples

Run the accepted-event HepMC analysis with:

```bash
cd /home/darthlucas/programs/test
bash cp_study/run_existing_samples.sh
```

For a quick smoke test:

```bash
bash cp_study/run_existing_samples.sh --max-events 1000 --outdir cp_study/results_smoke
```

For the largest available accepted-event pass in the current files:

```bash
bash cp_study/run_existing_samples.sh --max-events 1000000 --outdir cp_study/results_1M_mixing
```

## Observables

`cp_tag_asymmetry.png` uses the root B flavour in the accepted HepMC event:

```text
(N(initial B0) - N(initial anti-B0)) / (N(initial B0) + N(initial anti-B0))
```

This is the CP observable to inspect for SSD_CP in `Bd -> J/psi K_S`.

`initial_final_flavor_asymmetry.png` uses:

```text
(N(initial flavour == final flavour) - N(initial flavour != final flavour)) / sum
```

This is useful as a mixing cross-check. It is not the CP asymmetry for a CP
eigenstate. For `SSD_CP`, this plot is only a flavour-label diagnostic because
some oscillation-like chains are written as same-flavour `B -> B -> X`.

`record_mixing_asymmetry.png` uses the reconstructed B-chain classification:

```text
(N(unmixed) - N(mixed-like)) / sum
mixed-like = mixed + weird
```

This is the preferred record-level diagnostic when same-flavour `B -> B -> X`
chains are present. In the current `SSD_CP` sample it still does not close
against `cos(delta_m t)`, which means the exported B-chain labels should not be
used as a flavour-specific mixing observable for that model.

`mixing_closure.csv` stores the RMS and maximum absolute residual of both
mixing diagnostics relative to `cos(delta_m t)`.

`mixing_counts.csv` and `mixing_components.png` classify the accepted signal
chain as:

```text
unmixed: root B is the signal B
mixed:   root B and signal B have opposite flavour in a B -> B -> X chain
weird:   root B and signal B have the same flavour in a B -> B -> X chain
```

The `summary.json` also stores the expected mixed count from
`P_mix(t) = 0.5 * (1 - cos(delta_m * t))`, using the accepted candidate proper
times.

`../mixing_from_hepmc.py` can also make the standalone mixed/unmixed/weird time
histogram. It uses `pyhepmc_ng` when available, and otherwise falls back to this
study's HepMC2 text parser.

The parser reconstructs `B -> B -> X` chains by walking from the signal B back
to the root B. The proper time is computed from the root B flight time:

```text
c tau = delta(t_vertex) * m_B / E_B
t[ps] = c tau[mm] / c[mm/ps]
```

Debug ntuples written inside an EvtGen decay model should not be interpreted as
accepted events, because `EvtDecayAmp::makeDecay()` calls `decay()` during
accept/reject trials.

## Experimental SSD_CP Label Patch

The directory `evtgen_patch/` contains a local copy of EvtGen used only for an
SSD_CP label-cleaning test. It does not replace `/home/darthlucas/local/evtgen200`.

The patch changes `EvtCPUtil::OtherIncoherentB()` so that the internal SSD_CP
flavour-tag flip is used to choose `other_b`, but does not rewrite the B IDs in
the exported decay tree. It also keeps the lifetime already sampled by the
generic incoherent-mixing step instead of re-sampling it inside SSD_CP.

The generator `generator_ssdcp_label_patch.cc` also uses a local copy of the
Pythia8 EvtGen interface, `EvtGenCleanLabels.h`, for this isolated test.

Run the patched sample and analysis with:

```bash
bash cp_study/run_ssdcp_label_patch.sh 50000
```

Outputs are written under `cp_study/results_ssdcp_label_patch/`.

In the 50k-event patched test, the exported B-chain labels are cleaned:

```text
original SSD_CP: weird=1235, initial_B0_fraction=0.6163
patched SSD_CP:  weird=0,    initial_B0_fraction=0.5018
```

The `initial_final`/`record_mixing` curves still should not be interpreted as
the CP observable for SSD_CP; after the label patch they are clean record-level
diagnostics, while the internal SSD_CP flavour tag is no longer encoded in the
exported B PDG labels.

The generator now writes the accepted SSD_CP truth tag separately:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

Columns:

```text
hepmc_event, py_index, parent_pdg, other_b_pdg, flip, t_mm, t_ps
```

`ssdcp_truth_asymmetry.py` reads this CSV and writes:

```text
ssdcp_truth_tag_asymmetry.png
ssdcp_truth_asymmetry_bins.csv
```

For the 50k patched sample, the truth-tag asymmetry closes against the expected
`0.700 sin(delta_m t)` curve with RMS residual about `0.032`.
