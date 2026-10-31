# CP-Violation Investigation with SSD_CP in EvtGen/Pythia

This document summarizes the investigation of the `SSD_CP` model behavior for
the channel:

```text
Bd -> J/psi K_S
```

in a standalone `EvtGen 2.0.0 + Pythia 8.306` setup, with HepMC output. The
goal is to document the motivation, observed symptoms, current hypothesis,
tested modifications, results, and remaining open questions.

## Motivation

The study started because the sample generated with `SSD_CP` showed strange B
chains in the HepMC record, of the form:

```text
B -> B -> X
```

or, more concretely:

```text
anti-B0 -> anti-B0 -> J/psi K_S
```

This pattern did not appear in the `SVS` control sample, even though both
samples were produced with the same Pythia configuration and the same
EvtGen/Pythia interface.

The main question was whether these chains were:

1. an expected representation of mixing in EvtGen;
2. a bookkeeping artifact between EvtGen, Pythia and HepMC;
3. an indication that the `SSD_CP` model was not being used correctly in a pp
   collision context.

## Setup

The main samples were:

```text
SVS:
  generator: programs/test/generator3.cc
  decay file: Bd_JpsiKs_SVS.dec
  output: eventsSVS.hepmc

Original SSD_CP:
  generator: programs/test/generator2.cc
  decay file: Bd_JpsiKs_SSDCP.dec
  output: events_SSDCP2.hepmc

SSD_CP label patch:
  generator: cp_study/generator_ssdcp_label_patch.cc
  output: cp_study/results_ssdcp_label_patch/events_SSDCP_label_patch.hepmc
  truth tag: cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

The common Pythia configuration was:

```text
HardQCD:gg2bbbar = on
HardQCD:qqbar2bbbar = on
PartonLevel:MPI = off
PhaseSpace:pTHatMin = 8.0
ParticleDecays:mixB = off
511:mayDecay = off
-511:mayDecay = off
```

The B mixing used in this study came from the EvtGen/Pythia interface:

```cpp
EvtGenDecays(..., 1)
```

where the last argument enables the incoherent-mixing treatment in the EvtGen
interface path.

## Initial Observation

When walking through the HepMC tree of the accepted `B0/anti-B0 -> J/psi K_S`
candidates, the B chains were classified into three categories:

```text
unmixed: root B and signal B are the same particle
mixed:   there is a B -> B -> X chain with opposite flavours
weird:   there is a B -> B -> X chain with the same flavour
```

In the original 50k-event sample:

```text
SVS:
  signals = 43362
  unmixed = 35210
  mixed   = 8152
  weird   = 0

Original SSD_CP:
  signals = 43173
  unmixed = 35050
  mixed   = 6888
  weird   = 1235
```

The important point is that, in the original `SSD_CP` sample:

```text
mixed + weird = 8123
```

which was very close to the number of mixed events in the `SVS` sample:

```text
SVS mixed = 8152
```

This indicated that the `weird` category was probably not a new physical
population. It looked like part of the oscillation-like population whose label
had been contaminated.

## Current Hypothesis

The current hypothesis is that the same-flavour `B -> B -> X` chains arise from
an overlap between two mechanisms:

1. the generic incoherent B-mixing step in the EvtGen/Pythia interface;
2. the flavour-flip logic in `EvtCPUtil::OtherIncoherentB()`, used by
   CP-violation decay models such as `SSD_CP`.

In `SSD_CP`, the code calls:

```cpp
EvtCPUtil::getInstance()->OtherIncoherentB(p, t, other_b, 0.5);
```

This call selects the effective tag used by the CP model to compute the
amplitude. In the original code, however, this decision could also modify the
particle IDs in the EvtGen decay tree itself.

The problematic original lines were:

```cpp
p->getParent()->setId(EvtPDL::chargeConj(p->getParent()->getId()));
p->setId(EvtPDL::chargeConj(p->getId()));
```

In other words, an internal tag/flavour decision could be written directly into
the B-chain labels. When that chain already contained a generic mixing node, the
HepMC output could end up with an ambiguous or inconsistent genealogy, appearing
as a same-flavour `B -> B -> X` chain.

## Experimental Modification

A local copy of EvtGen was created to test an isolated modification:

```text
cp_study/evtgen_patch/R02-00-00
```

This copy does not replace the globally installed EvtGen.

### EvtCPUtil::OtherIncoherentB()

Modified file:

```text
cp_study/evtgen_patch/R02-00-00/src/EvtGenBase/EvtCPUtil.cpp
```

The main modification was to replace the flavour flip applied directly to the
tree IDs by a flip applied only to a local variable:

```cpp
EvtId effectiveTagSource = isMixed ? p->getParent()->getId() : p->getId();

if (isFlipped) {
    effectiveTagSource = EvtPDL::chargeConj(effectiveTagSource);
}

otherb = EvtPDL::chargeConj(effectiveTagSource);
```

Thus:

```text
before: the flip changed p->getId() and/or p->getParent()->getId()
after:  the flip changes only effectiveTagSource
```

The parent lifetime was also no longer re-sampled for already mixed candidates.
Before:

```cpp
p->getParent()->setLifetime();
t = p->getParent()->getLifetime();
```

After:

```cpp
t = p->getParent()->getLifetime();
```

The idea is that the generic mixing step had already assigned that lifetime, so
re-sampling it inside `OtherIncoherentB()` could decorrelate the mixing decision
from the exported time.

### EvtSSDCP.cpp

Modified file:

```text
cp_study/evtgen_patch/R02-00-00/src/EvtGenModels/EvtSSDCP.cpp
```

After the call to `OtherIncoherentB()`, debug attributes were added to the
`EvtParticle`:

```cpp
p->setAttribute("ssdcp_truth", 1);
p->setAttribute("ssdcp_parent_pdg", EvtPDL::getStdHep(p->getId()));
p->setAttribute("ssdcp_other_b_pdg", EvtPDL::getStdHep(other_b));
p->setAttribute("ssdcp_flip", flip ? 1 : 0);
p->setAttributeDouble("ssdcp_t_mm", t);
```

These attributes make it possible to export the accepted internal `SSD_CP` tag
without encoding it in the PDG labels of the B chain.

### Local EvtGen/Pythia Interface

A local copy of the interface was added:

```text
cp_study/EvtGenCleanLabels.h
```

It copies the `ssdcp_*` attributes out of the EvtGen tree and allows the
generator to write:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

with the columns:

```text
hepmc_event,py_index,parent_pdg,other_b_pdg,flip,t_mm,t_ps
```

## Result of the Modification

In the patched 50k-event sample:

```text
SSD_CP_LABEL_PATCH:
  signals = 43426
  unmixed = 35085
  mixed   = 8341
  weird   = 0
```

Integrated comparison:

```text
original SSD_CP: weird = 1235
patched SSD_CP:  weird = 0
```

This shows that the modification removed the visible symptom of same-flavour
`B -> B -> X` chains in the HepMC record.

The patched sample also recovered an initial `B0` fraction close to 50%:

```text
SSD_CP_LABEL_PATCH initial_B0_fraction = 0.5018
```

while the fraction inferred from the original HepMC labels was biased:

```text
Original SSD_CP initial_B0_fraction ~= 0.616
```

## Truth-Level CP Asymmetry

After the patch, the internal `SSD_CP` tag is no longer encoded in the HepMC PDG
labels. Therefore, the truth-level CP asymmetry was computed using:

```text
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
```

The convention used was:

```text
truth_signal_initial_pdg = -other_b_pdg
```

In time bins, the asymmetry is computed as:

```text
A_CP(t) =
  (N(truth initial B0, t) - N(truth initial anti-B0, t))
  /
  (N(truth initial B0, t) + N(truth initial anti-B0, t))
```

The corresponding script is:

```text
cp_study/ssdcp_truth_asymmetry.py
```

Main output:

```text
cp_study/results_ssdcp_label_patch/analysis_50000/ssdcp_truth_tag_asymmetry.png
```

For the patched 50k-event sample:

```text
truth rows = 43426
RMS residual vs 0.700 sin(delta_m t) ~= 0.032
```

This indicates that the exported internal tag reproduces the expected
truth-level CP curve:

```text
A_CP(t) ~= 0.700 sin(delta_m t)
```

with:

```text
delta_m = 0.502 ps^-1
```

## Why record_mixing_asymmetry Is Not the CP Curve

The plot:

```text
record_mixing_asymmetry.png
```

measures:

```text
A_record(t) = (N_unmixed - N_mixed) / sum
```

or, in the original diagnostic:

```text
mixed_like = mixed + weird
```

This observable tests the consistency of the HepMC chain as a flavour-specific
mixing measurement. For a pure mixing sample, one expects:

```text
A_mix(t) ~= cos(delta_m t)
```

In `SVS`, this closure works well:

```text
SVS record_mixing RMS vs cos = 0.026
```

In `SSD_CP_LABEL_PATCH`, even with no `weird` chains, it does not close:

```text
SSD_CP_LABEL_PATCH record_mixing RMS vs cos = 0.561
```

Current interpretation:

```text
the patched HepMC record is clean as a particle record, but the mixed/unmixed
classification of the signal B in SSD_CP should not be interpreted as a
physical time-dependent mixing observable.
```

The conceptual reason is that `Bd -> J/psi K_S` is a common final state / CP
eigenstate:

```text
B0      -> J/psi K_S
anti-B0 -> J/psi K_S
```

The `SSD_CP` model computes a coherent amplitude involving mixed and unmixed
components. Therefore, the relevant information for the CP asymmetry is the
flavour tag used by the model, not the `mixed/unmixed` classification inferred
afterwards from the HepMC tree of the signal B.

## Open Questions

The points that remain open are:

1. Confirm with the EvtGen developers whether the interpretation of the overlap
   between incoherent mixing and the flavour flip in
   `EvtCPUtil::OtherIncoherentB()` is correct.

2. Understand whether, for `SSD_CP` in a standalone pp setup, the
   `mixed/unmixed` asymmetry extracted directly from the HepMC record should or
   should not reproduce a `cos(delta_m t)` curve.

3. Decide the correct way to export the internal `SSD_CP` tag in a real
   production workflow. The current CSV is useful for debugging, but it is not a
   final production solution.

4. Evaluate whether the final solution should use:

```text
HepMC3 attributes
CMSSW EDM truth product
or another metadata structure in the EvtGen/Pythia interface
```

5. Check whether the modification in `OtherIncoherentB()` is physically
   acceptable in general, or whether it should remain restricted to a specific
   test/debug mode.

6. Clarify the EvtGen documentation regarding CP models in incoherent
   production. The documentation mentions that CP models may modify the flavour
   of the chain in order to produce the correct fractions, but it does not
   clearly describe the consequence for HepMC export in a pp context.

## Main Files

Documentation:

```text
cp_study/SSD_CP_modifications_explained.md
cp_study/logs/2026-05-11_weird_chain_investigation.md
cp_study/logs/2026-05-11_complete_modification_report.md
```

Generators:

```text
programs/test/generator2.cc
programs/test/generator3.cc
cp_study/generator_ssdcp_label_patch.cc
cp_study/generator_ssdcp_mixing0_test.cc
```

Experimental patch:

```text
cp_study/evtgen_patch/R02-00-00/src/EvtGenBase/EvtCPUtil.cpp
cp_study/evtgen_patch/R02-00-00/src/EvtGenModels/EvtSSDCP.cpp
cp_study/EvtGenCleanLabels.h
```

Analysis:

```text
cp_study/hepmc_cp_asymmetry.py
cp_study/ssdcp_truth_asymmetry.py
mixing_from_hepmc.py
```

Main results:

```text
cp_study/results_ssdcp_label_patch/events_SSDCP_label_patch.hepmc
cp_study/results_ssdcp_label_patch/ssdcp_truth_tags.csv
cp_study/results_ssdcp_label_patch/analysis_50000/
```

## One-Sentence Summary

The experimental patch separates two pieces of information that were entangled
in the original sample: the HepMC labels of the B chain and the internal CP tag
used by `SSD_CP`. With this separation, the same-flavour `B -> B -> X` chains
disappear from the HepMC record, and the truth-level CP asymmetry can be
extracted from the exported internal tag; it remains open whether the
`mixed/unmixed` asymmetry of the `SSD_CP` HepMC record has a physical
interpretation as a pure mixing curve.
