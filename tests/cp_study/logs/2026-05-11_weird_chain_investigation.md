# 2026-05-11 - SSD_CP weird-chain investigation

## Question

Why does the SVS sample have no `weird` chains, while the SSD_CP sample has a non-negligible number of `B -> B -> X` chains?

## Main result

The `weird` class is almost certainly an event-record/flavour-bookkeeping effect from the `SSD_CP` model plus the Pythia8 EvtGen interface, not a new physical population.

In the 50k-event samples:

| model | signals | unmixed | mixed | weird | mixed + weird |
|---|---:|---:|---:|---:|---:|
| SVS | 43362 | 35210 | 8152 | 0 | 8152 |
| SSD_CP | 43173 | 35050 | 6888 | 1235 | 8123 |

For SSD_CP, `mixed + weird = 8123`, close to the SVS mixed count (`8152`) and to the delta-m expectation in the analysis (`7923.1`). This is why `mixed` alone undercounts the oscillation-like population in SSD_CP.

## Source-level evidence

- `Bd_JpsiKs_SVS.dec` uses the plain `SVS` model:
  - `Decay B0`
  - `1.000 J/psi K_S0 SVS ;`
- `Bd_JpsiKs_SSDCP.dec` uses `SSD_CP`:
  - `Decay B0`
  - `1.000 K_S0 J/psi SSD_CP 0.502e12 0.0 1.0 -0.775 1.0 0. 1.0 0.;`
- `EvtSVS::decay()` only initializes phase space and vertices. It does not call `EvtCPUtil` and does not flip the B flavour ID.
- `EvtSSDCP::decay()` calls:
  - `EvtCPUtil::getInstance()->OtherIncoherentB(p, t, other_b, 0.5);`
- `EvtCPUtil::OtherIncoherentB()` explicitly enables flavour flipping. If the particle is already a mixed B node, it can flip both the parent and the child IDs:
  - `p->getParent()->setId(chargeConj(parent_id))`
  - `p->setId(chargeConj(p_id))`
- Generic EvtGen incoherent mixing inserts a one-daughter oscillation node only when the daughter is the charge conjugate of the mother.
- The Pythia8 EvtGen interface flags oscillation only when an EvtGen particle has exactly one daughter and the daughter ID differs from the mother ID. It appends the EvtGen daughter to the original Pythia particle; it does not retroactively rewrite the original Pythia mother ID in the HepMC record.

## Concrete HepMC pattern

SVS mixed example:

```text
B0 -> anti-B0 -> J/psi K_S0
```

SSD_CP weird example:

```text
anti-B0 -> anti-B0 -> K_S0 J/psi
```

The second case is counted as `weird` by the parser because there is a B-to-B one-daughter chain, but the initial and final B IDs are the same. This appears when the SSD_CP flavour flip changes the EvtGen-side IDs around a generic oscillation node.

## Cross-check

A temporary SSD_CP generator was built with the constructor mixing argument changed from `1` to `0`:

- `cp_study/generator_ssdcp_mixing0_test.cc`
- `cp_study/results_mixing0_test/events_SSDCP_mixing0.hepmc`

Result for 10k events:

| model | signals | unmixed | mixed | weird |
|---|---:|---:|---:|---:|
| SSD_CP_MIX0 | 8553 | 6930 | 1390 | 233 |

The weird fraction stayed similar (`2.72%` vs `2.86%`). Important caveat: in this EvtGen code path, `0` maps to `EvtCPUtil::Coherent`, not to "no mixing", so this is not a clean no-mixing test.

## Conclusion

SVS has no `weird` chains because it leaves the generic EvtGen mixing record untouched: mixed events look like opposite-flavour `B -> anti-B` chains, and unmixed events are direct.

SSD_CP has `weird` chains because `EvtSSDCP` uses `OtherIncoherentB`, which can flip flavour IDs for CP tagging. In combination with the Pythia8 EvtGen interface, this can leave the exported HepMC record with same-flavour `B -> B -> X` chains. For SSD_CP, use `mixed + weird` as the oscillation-like chain count; using `mixed` alone removes a real part of the record-level mixing population.

## Follow-up: binned mixing closure

Adding `weird` back into the mixed-like count fixes the integrated count, but it does not make the SSD_CP B-chain labels close against `cos(delta_m t)` bin-by-bin.

The updated analysis now writes:

- `record_mixing_asymmetry.png`
- `mixing_closure.csv`

For the current 50k-event samples:

| model | observable | RMS vs cos(delta_m t) | max abs residual |
|---|---|---:|---:|
| SVS | initial_final | 0.0263 | 0.0719 |
| SVS | record_mixing | 0.0263 | 0.0719 |
| SSD_CP | initial_final | 1.1306 | 1.7197 |
| SSD_CP | record_mixing | 1.0938 | 1.6601 |

So the SSD_CP HepMC B-chain labels are not a reliable flavour-specific mixing observable as a function of proper time. The CP-tag asymmetry remains the relevant SSD_CP observable.

## Experimental label-cleaning patch

A separate EvtGen copy was made under:

- `cp_study/evtgen_patch/R02-00-00`

The installed EvtGen under `/home/darthlucas/local/evtgen200` was not replaced.

Two experimental changes were tested:

- In `EvtCPUtil::OtherIncoherentB()`, do not re-sample the parent lifetime when the B already has a generic incoherent-mixing parent.
- Use the SSD_CP random flavour flip only to choose the internal `other_b` tag, without rewriting the B IDs in the exported decay tree.

A local Pythia8 EvtGen interface copy was also added:

- `cp_study/EvtGenCleanLabels.h`

The patched generator is:

- `cp_study/generator_ssdcp_label_patch.cc`

Run script:

- `cp_study/run_ssdcp_label_patch.sh`

50k-event patched result:

| sample | signals | initial B0 fraction | unmixed | mixed | weird |
|---|---:|---:|---:|---:|---:|
| SSD_CP original | 43173 | 0.6163 | 35050 | 6888 | 1235 |
| SSD_CP label patch | 43426 | 0.5018 | 35085 | 8341 | 0 |

This confirms that the exported B-chain labels can be cleaned. However, after this patch the internal SSD_CP CP tag is no longer encoded in the exported B PDG labels. Therefore `cp_tag_asymmetry.png` from this patched sample should not be used as the SSD_CP CP asymmetry unless a separate truth tag is exported.

## Accepted truth-tag export

The patch now stores the accepted SSD_CP truth tag as EvtParticle attributes and copies them out through the local Pythia8 EvtGen interface. The generator writes:

- `cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv`

Columns:

| column | meaning |
|---|---|
| `hepmc_event` | zero-based HepMC event number written by the generator |
| `py_index` | Pythia event-record index for the decayed B |
| `parent_pdg` | clean exported B PDG label after the label patch |
| `other_b_pdg` | SSD_CP internal flavour tag used in the amplitude |
| `flip` | non-CP final-state daughter flip flag |
| `t_mm` | SSD_CP internal time in mm |
| `t_ps` | SSD_CP internal time in ps |

The script `cp_study/ssdcp_truth_asymmetry.py` plots the truth-tag asymmetry.
For the 50k patched sample:

- truth rows: `43426`
- events with truth: `33942`
- RMS residual vs `0.700 sin(delta_m t)`: `0.0322`

The output plot is:

- `cp_study/results_ssdcp_label_patch/analysis_50000/ssdcp_truth_tag_asymmetry.png`
