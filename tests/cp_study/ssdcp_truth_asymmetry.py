#!/usr/bin/env python3
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


plt.style.use(hep.style.CMS)

B0 = 511
B0BAR = -511


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot SSD_CP truth-tag asymmetry.")
    parser.add_argument("truth_csv", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bins", type=int, default=16)
    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=8.0)
    parser.add_argument("--delta-m", type=float, default=0.502)
    parser.add_argument("--reference-sin-coeff", type=float, default=0.700)
    return parser.parse_args()


def read_truth(path: Path) -> list[dict]:
    with path.open() as handle:
        return [
            {
                "hepmc_event": int(row["hepmc_event"]),
                "py_index": int(row["py_index"]),
                "parent_pdg": int(row["parent_pdg"]),
                "other_b_pdg": int(row["other_b_pdg"]),
                "flip": int(row["flip"]),
                "t_mm": float(row["t_mm"]),
                "t_ps": float(row["t_ps"]),
            }
            for row in csv.DictReader(handle)
        ]


def truth_rows(events: list[dict], bins: np.ndarray) -> list[dict]:
    rows = []
    t = np.array([event["t_ps"] for event in events])
    other_b = np.array([event["other_b_pdg"] for event in events])

    # In EvtSSDCP, other_b == anti-B0 means the signal-side B was B0 at t=0.
    signal_initial = -other_b
    for low, high in zip(bins[:-1], bins[1:]):
        mask = (t >= low) & (t < high)
        initial_b0 = int(np.sum(mask & (signal_initial == B0)))
        initial_antib0 = int(np.sum(mask & (signal_initial == B0BAR)))
        total = initial_b0 + initial_antib0
        asym = float("nan")
        err = float("nan")
        if total > 0:
            asym = (initial_b0 - initial_antib0) / total
            err = math.sqrt(max(0.0, 1.0 - asym * asym) / total)
        rows.append(
            {
                "bin_low_ps": low,
                "bin_high_ps": high,
                "bin_center_ps": 0.5 * (low + high),
                "count": total,
                "truth_initial_B0": initial_b0,
                "truth_initial_antiB0": initial_antib0,
                "asymmetry": asym,
                "error": err,
            }
        )
    return rows


def write_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows: list[dict], outpath: Path, delta_m: float, reference_sin_coeff: float) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 8.0))
    centers = np.array([row["bin_center_ps"] for row in rows])
    asym = np.array([row["asymmetry"] for row in rows])
    err = np.array([row["error"] for row in rows])
    valid = np.isfinite(asym)
    ax.errorbar(centers[valid], asym[valid], yerr=err[valid], fmt="o", capsize=2, label="SSD_CP truth tag")

    t_ref = np.linspace(rows[0]["bin_low_ps"], rows[-1]["bin_high_ps"], 500)
    ax.plot(
        t_ref,
        reference_sin_coeff * np.sin(delta_m * t_ref),
        "k--",
        linewidth=1.2,
        label=f"{reference_sin_coeff:.3f} sin(delta_m t)",
    )
    ax.axhline(0.0, color="0.4", linewidth=0.8)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("Asymmetry")
    ax.set_ylim(-1.1, 1.1)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    events = read_truth(args.truth_csv)
    bins = np.linspace(args.t_min, args.t_max, args.bins + 1)
    rows = truth_rows(events, bins)
    write_rows(args.outdir / "ssdcp_truth_asymmetry_bins.csv", rows)
    plot_rows(rows, args.outdir / "ssdcp_truth_tag_asymmetry.png", args.delta_m, args.reference_sin_coeff)
    print(f"Wrote SSD_CP truth-tag results to {args.outdir}")
    print(f"truth_rows={len(events)} events_with_truth={len({event['hepmc_event'] for event in events})}")


if __name__ == "__main__":
    main()
