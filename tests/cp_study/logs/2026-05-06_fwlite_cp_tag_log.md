# 2026-05-06 FWLite CP-Tag Script Log

## Scope

Add a FWLite ROOT version of the accepted-record CP-tag reconstruction used in
the HepMC study.

## Files Created

- `cp_study/fwlite_cp_asymmetry.py`
  - Reason: apply the same concept to `std::vector<reco::GenParticle>` in a
    FWLite ROOT file.

No existing file was modified. Therefore no backup was needed for this step.

## What The Script Does

- Selects terminal `B0`/`anti-B0` candidates with descendants `J/psi` and
  `K_S0`.
- Walks upward through `B0`/`anti-B0` mothers to reconstruct the root B flavour.
- Computes proper time from the root-B production vertex to the terminal-B
  decay vertex.
- Produces:
  - `*_cp_tag.png`: initial `B0` vs initial `anti-B0` asymmetry.
  - `*_mixing.png`: same initial/final flavour vs opposite initial/final
    flavour asymmetry.
  - `*_candidates.csv`: selected candidate truth rows.
  - `*_asymmetry_bins.csv`: binned asymmetry values.

## Important Distinction

The FWLite script still reconstructs the `OtherB`-equivalent tag from the final
accepted truth record. It does not read the internal `other_b` variable from
EvtGen.

## Update: Missing `endVertex()`

After testing on lxplus, `reco::GenParticle` did not expose `endVertex()`.

Backup made before editing:

- `cp_study/fwlite_cp_asymmetry.py.bak_20260506_112638`

Change:

- Added `decay_vertex(p)`, which uses `p.endVertex()` when available and falls
  back to the production vertex of the first daughter.
- Updated `particle_key()` and `proper_time_ps()` to use this helper.

Reason:

- Standard FWLite `reco::GenParticle` often stores a particle's production
  vertex as `vertex()`, but not the parent's decay vertex as `endVertex()`.
  For a decayed B, the daughters' `vertex()` is the B decay vertex.
