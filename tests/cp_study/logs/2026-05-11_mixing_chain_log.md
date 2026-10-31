# 2026-05-11 Mixing Chain Follow-Up

## Scope

Add the `mixing_from_hepmc.py` style mixed/unmixed/weird chain classification
to the accepted HepMC CP study for the existing `eventsSVS.hepmc` and
`events_SSDCP2.hepmc` samples.

## Implementation

- `cp_study/hepmc_cp_asymmetry.py`
  - Added `mixing_class` per accepted signal candidate.
  - Added `mixing_counts.csv`.
  - Added `mixing_components.png`.
  - Added expected mixed count/fraction from
    `P_mix(t) = 0.5 * (1 - cos(delta_m * t))`.
- `mixing_from_hepmc.py`
  - Added `--max-events` / `--max`.
  - Restricted counting to root B candidates to avoid double-counting the same
    `B -> B -> X` chain.
  - Added a fallback to the study's HepMC2 text parser when `pyhepmc_ng` is not
    installed.

## Commands Run

```bash
cd /home/darthlucas/programs/test
python3 -m py_compile cp_study/hepmc_cp_asymmetry.py mixing_from_hepmc.py
python3 mixing_from_hepmc.py eventsSVS.hepmc -o cp_study/results_smoke_mixing/SVS_mixing_from_hepmc.png --max 1000 --debug-events 0
python3 mixing_from_hepmc.py events_SSDCP2.hepmc -o cp_study/results_smoke_mixing/SSD_CP_mixing_from_hepmc.png --max 1000 --debug-events 0
python3 cp_study/hepmc_cp_asymmetry.py --sample SVS=eventsSVS.hepmc --sample SSD_CP=events_SSDCP2.hepmc --max-events 1000 --outdir cp_study/results_smoke_mixing
python3 cp_study/hepmc_cp_asymmetry.py --sample SVS=eventsSVS.hepmc --sample SSD_CP=events_SSDCP2.hepmc --max-events 1000000 --outdir cp_study/results_1M_mixing
python3 mixing_from_hepmc.py eventsSVS.hepmc -o cp_study/results_1M_mixing/SVS_mixing_from_hepmc.png --max 1000000 --debug-events 0
python3 mixing_from_hepmc.py events_SSDCP2.hepmc -o cp_study/results_1M_mixing/SSD_CP_mixing_from_hepmc.png --max 1000000 --debug-events 0
```

The environment does not have `pyhepmc_ng` installed, so the standalone
`mixing_from_hepmc.py` commands used the text-parser fallback.

## Results

The `--max-events 1000000` pass read all events available in the current files:
50,000 events per sample.

```text
SVS:
  signals: 43,362
  unmixed: 35,210
  mixed:    8,152
  weird:        0

SSD_CP:
  signals: 43,173
  unmixed: 35,050
  mixed:    6,888
  weird:    1,235
```

The same-flavour `B -> B -> X` chains are therefore present in SSD_CP at
`1235 / 43173 = 2.86%` of accepted signal candidates, while SVS has none in
this sample.
