#!/usr/bin/env python3
"""FWLite truth-level Bd -> J/psi K_S asymmetries.

This applies the same accepted-record idea used in the HepMC study to a CMS
FWLite ROOT file: reconstruct the terminal signal B, walk back to the root B,
then compare either initial B0 vs initial anti-B0 (CP tag) or same vs opposite
initial/final flavour (mixing cross-check).
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

try:
    from DataFormats.FWLite import Events, Handle
except ImportError:
    Events = None
    Handle = None


plt.style.use(hep.style.CMS)

MB = 5.27963
C_CM_PER_PS = 0.0299792458
B0 = 511
B0BAR = -511
JPSI = 443
KS0 = 310


def daughters(p):
    if p is None:
        return []
    return [p.daughter(i) for i in range(p.numberOfDaughters()) if p.daughter(i) is not None]


def mothers(p):
    if p is None:
        return []
    if hasattr(p, "numberOfMothers") and p.numberOfMothers() > 0:
        return [p.mother(i) for i in range(p.numberOfMothers()) if p.mother(i) is not None]
    m = p.mother()
    return [] if m is None else [m]


def has_desc(p, target_abs, depth=0, maxd=30):
    if p is None or depth > maxd:
        return False
    for d in daughters(p):
        if abs(d.pdgId()) == target_abs:
            return True
        if has_desc(d, target_abs, depth + 1, maxd):
            return True
    return False


def has_signal_bd_daughter(p):
    for d in daughters(p):
        if abs(d.pdgId()) == B0 and has_desc(d, JPSI) and has_desc(d, KS0):
            return True
    return False


def is_terminal_signal_bd(p):
    return (
        p is not None
        and abs(p.pdgId()) == B0
        and has_desc(p, JPSI)
        and has_desc(p, KS0)
        and not has_signal_bd_daughter(p)
    )


def root_bd(b):
    """Walk through Bd/Bdbar mothers to the first B in the signal chain."""
    cur = b
    seen = set()
    for _ in range(80):
        if cur is None:
            return None
        key = particle_key(cur)
        if key in seen:
            return cur
        seen.add(key)

        bd_mothers = [m for m in mothers(cur) if abs(m.pdgId()) == B0]
        if not bd_mothers:
            return cur

        chosen = None
        for m in bd_mothers:
            if has_desc(m, JPSI) and has_desc(m, KS0):
                chosen = m
                break
        cur = chosen if chosen is not None else bd_mothers[0]
    return cur


def particle_key(p):
    """Approximate stable key for de-duplicating GenParticle handles."""
    v = p.vertex()
    ev = decay_vertex(p)
    return (
        int(p.pdgId()),
        round(float(p.px()), 10),
        round(float(p.py()), 10),
        round(float(p.pz()), 10),
        vertex_key(v),
        vertex_key(ev),
    )


def vertex_key(v):
    if v is None:
        return None
    return (round(float(v.x()), 10), round(float(v.y()), 10), round(float(v.z()), 10))


def decay_vertex(p):
    """Return the decay vertex for FWLite GenParticle variants.

    Some ROOT files expose HepMC-like `endVertex()`, while standard
    `reco::GenParticle` often only exposes the production vertex of each
    daughter. For a decayed parent, the daughters' production vertex is the
    parent's decay vertex.
    """
    if p is None:
        return None

    end_vertex = getattr(p, "endVertex", None)
    if callable(end_vertex):
        try:
            v = end_vertex()
            if v is not None:
                return v
        except AttributeError:
            pass

    for d in daughters(p):
        v = d.vertex()
        if v is not None:
            return v

    return None


def proper_time_ps(b_root, b_term):
    """Proper time from root-B production vertex to terminal-B decay vertex."""
    if b_root is None or b_term is None:
        return None

    v0 = b_root.vertex()
    v1 = decay_vertex(b_term)
    if v0 is None or v1 is None:
        return None

    p = b_term.p()
    if p is None or p <= 0.0:
        return None

    dx = v1.x() - v0.x()
    dy = v1.y() - v0.y()
    dz = v1.z() - v0.z()
    projected_l = (dx * b_term.px() + dy * b_term.py() + dz * b_term.pz()) / float(p)
    return (MB / float(p)) * abs(projected_l) / C_CM_PER_PS


def asymmetry_and_err(n_plus, n_minus):
    total = n_plus + n_minus
    if total <= 0:
        return math.nan, math.nan
    asym = (n_plus - n_minus) / total
    err = math.sqrt(max(0.0, 1.0 - asym * asym) / total)
    return asym, err


def output_path(prefix, suffix):
    prefix = Path(prefix)
    if prefix.suffix:
        return prefix.with_name(f"{prefix.stem}_{suffix}{prefix.suffix}")
    return prefix.with_name(f"{prefix.name}_{suffix}.png")


def add_cms_label(ax):
    hep.cms.label(data=False, com=13.6, ax=ax)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Make FWLite Bd -> J/psi K_S CP-tag and mixing asymmetry plots."
    )
    parser.add_argument("input", help="Input FWLite ROOT file")
    parser.add_argument("-o", "--output-prefix", default="fwlite_asymmetry")
    parser.add_argument("--nbins", type=int, default=25)
    parser.add_argument("--tmin", type=float, default=0.0)
    parser.add_argument("--tmax", type=float, default=10.0)
    parser.add_argument("--dm", type=float, default=None, help="Optional delta_m overlay in ps^-1")
    parser.add_argument("--label", default=("genParticles", "", "GEN"), nargs=3)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--debug-events", type=int, default=0)
    return parser.parse_args()


def collect_candidates(fname, label, max_events=None, debug_events=0):
    if Events is None or Handle is None:
        raise RuntimeError("DataFormats.FWLite is not available. Run this script inside a CMSSW/FWLite environment.")

    events = Events(fname)
    handle = Handle("std::vector<reco::GenParticle>")

    candidates = []
    seen_per_event = set()

    for iev, ev in enumerate(events):
        if max_events is not None and iev >= max_events:
            break

        ev.getByLabel(tuple(label), handle)
        parts = handle.product()
        seen_per_event.clear()

        for p in parts:
            if not is_terminal_signal_bd(p):
                continue

            term = p
            root = root_bd(term)
            if root is None:
                continue

            key = particle_key(term)
            if key in seen_per_event:
                continue
            seen_per_event.add(key)

            t = proper_time_ps(root, term)
            if t is None:
                continue

            initial = int(root.pdgId())
            final = int(term.pdgId())
            if abs(initial) != B0 or abs(final) != B0:
                continue

            candidates.append(
                {
                    "event": iev,
                    "t_ps": t,
                    "initial_pdg": initial,
                    "final_pdg": final,
                    "mixed_chain": int(initial != final),
                    "other_b_equiv": -initial,
                }
            )

            if iev < debug_events:
                print(
                    f"event={iev} t={t:.4f} ps "
                    f"initial={initial} final={final} other_b_equiv={-initial}"
                )

    return candidates


def make_rows(candidates, bins, observable):
    rows = []
    t = np.array([c["t_ps"] for c in candidates])
    initial = np.array([c["initial_pdg"] for c in candidates], dtype=int)
    final = np.array([c["final_pdg"] for c in candidates], dtype=int)

    for low, high in zip(bins[:-1], bins[1:]):
        mask = (t >= low) & (t < high)
        if observable == "cp_tag":
            plus = int(np.sum(mask & (initial == B0)))
            minus = int(np.sum(mask & (initial == B0BAR)))
            plus_name = "initial_B0"
            minus_name = "initial_antiB0"
        elif observable == "mixing":
            plus = int(np.sum(mask & (initial == final)))
            minus = int(np.sum(mask & (initial != final)))
            plus_name = "same_initial_final"
            minus_name = "opposite_initial_final"
        else:
            raise ValueError(observable)

        asym, err = asymmetry_and_err(plus, minus)
        rows.append(
            {
                "observable": observable,
                "bin_low_ps": low,
                "bin_high_ps": high,
                "bin_center_ps": 0.5 * (low + high),
                "count": plus + minus,
                plus_name: plus,
                minus_name: minus,
                "asymmetry": asym,
                "error": err,
            }
        )
    return rows


def save_csv(path, rows):
    if not rows:
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows, outpng, observable, tmin, tmax, dm):
    centers = np.array([r["bin_center_ps"] for r in rows])
    asym = np.array([r["asymmetry"] for r in rows])
    err = np.array([r["error"] for r in rows])
    valid = np.isfinite(asym)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.errorbar(centers[valid], asym[valid], yerr=err[valid], fmt="o", markersize=4, capsize=0)

    if dm is not None:
        tt = np.linspace(tmin, tmax, 600)
        if observable == "cp_tag":
            ax.plot(tt, 0.700 * np.sin(dm * tt), "-", linewidth=1, label="0.700 sin(delta_m t)")
        else:
            ax.plot(tt, np.cos(dm * tt), "-", linewidth=1, label="cos(delta_m t)")
        ax.legend(frameon=False, loc="best")

    ax.axhline(0.0, linewidth=1)
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("t (ps)")
    if observable == "cp_tag":
        ax.set_ylabel("(N(initial B0) - N(initial anti-B0)) / sum")
    else:
        ax.set_ylabel("(N(same) - N(opposite)) / sum")

    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    add_cms_label(ax)
    fig.tight_layout()
    fig.savefig(outpng, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    bins = np.linspace(args.tmin, args.tmax, args.nbins + 1)

    candidates = collect_candidates(
        args.input,
        args.label,
        max_events=args.max_events,
        debug_events=args.debug_events,
    )

    prefix = Path(args.output_prefix)
    save_csv(output_path(prefix, "candidates").with_suffix(".csv"), candidates)

    all_rows = []
    for observable in ("cp_tag", "mixing"):
        rows = make_rows(candidates, bins, observable)
        all_rows.extend(rows)
        plot_rows(
            rows,
            output_path(prefix, observable),
            observable,
            args.tmin,
            args.tmax,
            args.dm,
        )

    save_csv(output_path(prefix, "asymmetry_bins").with_suffix(".csv"), all_rows)

    n_initial_b0 = sum(c["initial_pdg"] == B0 for c in candidates)
    n_initial_b0bar = sum(c["initial_pdg"] == B0BAR for c in candidates)
    n_same = sum(c["initial_pdg"] == c["final_pdg"] for c in candidates)
    n_opposite = sum(c["initial_pdg"] != c["final_pdg"] for c in candidates)

    print(f"Input: {args.input}")
    print(f"Selected terminal Bd->J/psi K_S candidates: {len(candidates)}")
    print(f"Initial B0      : {n_initial_b0}")
    print(f"Initial anti-B0 : {n_initial_b0bar}")
    print(f"Same init/final : {n_same}")
    print(f"Opposite        : {n_opposite}")
    print(f"Saved: {output_path(prefix, 'cp_tag')}")
    print(f"Saved: {output_path(prefix, 'mixing')}")
    print(f"Saved: {output_path(prefix, 'asymmetry_bins').with_suffix('.csv')}")


if __name__ == "__main__":
    main()
