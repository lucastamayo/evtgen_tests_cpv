#!/usr/bin/env python3
"""Plot Bd -> J/psi K_S asymmetries from accepted HepMC2 events.

The SSD_CP model chooses the flavour tag and proper time inside the EvtGen
accept/reject machinery. Debug lines written from model::decay therefore include
rejected trials; this script reads the final HepMC event record instead.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np


plt.style.use(hep.style.CMS)

C_MM_PER_PS = 0.299792458
B0 = 511
B0BAR = -511
JPSI = 443
KS0 = 310


def add_cms_label(ax) -> None:
    return None


@dataclass
class Particle:
    barcode: int
    pid: int
    px: float
    py: float
    pz: float
    e: float
    mass: float
    status: int
    prod_vertex: int
    end_vertex: int


@dataclass
class SignalEvent:
    model: str
    event: int
    t_ps: float
    ctau_mm: float
    initial_pdg: int
    final_pdg: int
    root_barcode: int
    signal_barcode: int
    root_end_vertex: int
    signal_end_vertex: int
    mixed_chain: bool
    mixing_class: str


def parse_sample_arg(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("sample must be LABEL=PATH")
    label, path = value.split("=", 1)
    if not label:
        raise argparse.ArgumentTypeError("sample label is empty")
    return label, Path(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract accepted Bd -> J/psi K_S candidates from HepMC2 and plot asymmetries."
    )
    parser.add_argument(
        "--sample",
        action="append",
        type=parse_sample_arg,
        required=True,
        help="Sample as LABEL=PATH. Repeat for multiple models.",
    )
    parser.add_argument("--outdir", type=Path, default=Path("cp_study/results"))
    parser.add_argument("--bins", type=int, default=16)
    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=8.0)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--delta-m", type=float, default=0.502, help="Bd delta_m in ps^-1.")
    parser.add_argument(
        "--reference-sin-coeff",
        type=float,
        default=0.700,
        help="Reference S coefficient for the SSD_CP CP asymmetry overlay.",
    )
    return parser.parse_args()


def finish_event(
    model: str,
    event_number: int,
    particles: dict[int, Particle],
    vertices: dict[int, dict],
) -> list[SignalEvent]:
    if not particles:
        return []

    for particle in particles.values():
        if particle.end_vertex in vertices:
            vertices[particle.end_vertex]["in"].append(particle.barcode)

    signals: list[SignalEvent] = []
    for particle in particles.values():
        if abs(particle.pid) != B0 or particle.end_vertex not in vertices:
            continue

        daughter_pids = [
            particles[barcode].pid
            for barcode in vertices[particle.end_vertex]["out"]
            if barcode in particles
        ]
        if JPSI not in daughter_pids or KS0 not in daughter_pids:
            continue

        root = find_root_b(particle, particles, vertices)
        prod_pos = vertices.get(root.prod_vertex, {}).get("pos", (0.0, 0.0, 0.0, 0.0))
        decay_pos = vertices.get(root.end_vertex, {}).get("pos")
        if decay_pos is None:
            decay_pos = vertices.get(particle.end_vertex, {}).get("pos")
        if decay_pos is None:
            continue

        ctau_mm = proper_ctau_mm(root, prod_pos, decay_pos)
        if not math.isfinite(ctau_mm):
            continue

        signals.append(
            SignalEvent(
                model=model,
                event=event_number,
                t_ps=ctau_mm / C_MM_PER_PS,
                ctau_mm=ctau_mm,
                initial_pdg=root.pid,
                final_pdg=particle.pid,
                root_barcode=root.barcode,
                signal_barcode=particle.barcode,
                root_end_vertex=root.end_vertex,
                signal_end_vertex=particle.end_vertex,
                mixed_chain=(root.barcode != particle.barcode or root.pid != particle.pid),
                mixing_class=classify_mixing_chain(root, particle),
            )
        )

    return signals


def classify_mixing_chain(root: Particle, signal: Particle) -> str:
    if root.barcode == signal.barcode:
        return "unmixed"
    if root.pid == -signal.pid:
        return "mixed"
    return "weird"


def find_root_b(
    particle: Particle, particles: dict[int, Particle], vertices: dict[int, dict]
) -> Particle:
    root = particle
    seen = {root.barcode}
    while True:
        prod_vertex = vertices.get(root.prod_vertex)
        if prod_vertex is None:
            return root
        parent_bs = [
            particles[barcode]
            for barcode in prod_vertex["in"]
            if barcode in particles and abs(particles[barcode].pid) == B0
        ]
        if len(parent_bs) != 1 or parent_bs[0].barcode in seen:
            return root
        root = parent_bs[0]
        seen.add(root.barcode)


def proper_ctau_mm(
    particle: Particle, prod_pos: tuple[float, float, float, float], decay_pos: tuple[float, float, float, float]
) -> float:
    dt = decay_pos[3] - prod_pos[3]
    if particle.e > 0.0:
        return dt * particle.mass / particle.e

    dx = decay_pos[0] - prod_pos[0]
    dy = decay_pos[1] - prod_pos[1]
    dz = decay_pos[2] - prod_pos[2]
    p2 = particle.px * particle.px + particle.py * particle.py + particle.pz * particle.pz
    if p2 <= 0.0:
        return float("nan")
    return (dx * particle.px + dy * particle.py + dz * particle.pz) * particle.mass / p2


def read_hepmc(model: str, path: Path, max_events: int | None) -> tuple[list[SignalEvent], dict]:
    events: list[SignalEvent] = []
    particles: dict[int, Particle] = {}
    vertices: dict[int, dict] = {}
    current_vertex: int | None = None
    event_number = -1
    skipped_short_particles = 0

    def flush_current_event() -> bool:
        nonlocal particles, vertices
        events.extend(finish_event(model, event_number, particles, vertices))
        particles = {}
        vertices = {}
        return max_events is not None and event_number + 1 >= max_events

    with path.open() as handle:
        for line in handle:
            if line.startswith("E "):
                if event_number >= 0 and flush_current_event():
                    break
                event_number += 1
                current_vertex = None
                continue

            if line.startswith("V "):
                fields = line.split()
                if len(fields) < 7:
                    continue
                current_vertex = int(fields[1])
                vertices[current_vertex] = {
                    "pos": tuple(float(value) for value in fields[3:7]),
                    "in": [],
                    "out": [],
                }
                continue

            if line.startswith("P ") and current_vertex is not None:
                fields = line.split()
                if len(fields) < 12:
                    skipped_short_particles += 1
                    continue
                barcode = int(fields[1])
                particle = Particle(
                    barcode=barcode,
                    pid=int(fields[2]),
                    px=float(fields[3]),
                    py=float(fields[4]),
                    pz=float(fields[5]),
                    e=float(fields[6]),
                    mass=float(fields[7]),
                    status=int(fields[8]),
                    prod_vertex=current_vertex,
                    end_vertex=int(fields[11]),
                )
                particles[barcode] = particle
                vertices[current_vertex]["out"].append(barcode)

    if event_number >= 0:
        events.extend(finish_event(model, event_number, particles, vertices))

    summary = {
        "model": model,
        "path": str(path),
        "events_seen": event_number + 1,
        "signals": len(events),
        "skipped_short_particles": skipped_short_particles,
    }
    return events, summary


def asymmetry_rows(
    model: str,
    events: list[SignalEvent],
    bins: np.ndarray,
    observable: str,
) -> list[dict]:
    rows = []
    t = np.array([event.t_ps for event in events])
    initial = np.array([event.initial_pdg for event in events])
    final = np.array([event.final_pdg for event in events])
    mixing_class = np.array([event.mixing_class for event in events])

    for low, high in zip(bins[:-1], bins[1:]):
        mask = (t >= low) & (t < high)
        if observable == "cp_tag":
            plus = int(np.sum(mask & (initial == B0)))
            minus = int(np.sum(mask & (initial == B0BAR)))
            plus_label = "initial_B0"
            minus_label = "initial_antiB0"
        elif observable == "initial_final":
            plus = int(np.sum(mask & (initial == final)))
            minus = int(np.sum(mask & (initial != final)))
            plus_label = "same_initial_final"
            minus_label = "opposite_initial_final"
        elif observable == "record_mixing":
            mixed_like = (mixing_class == "mixed") | (mixing_class == "weird")
            plus = int(np.sum(mask & (mixing_class == "unmixed")))
            minus = int(np.sum(mask & mixed_like))
            plus_label = "unmixed_chain"
            minus_label = "mixed_like_chain"
        else:
            raise ValueError(f"unknown observable: {observable}")

        total = plus + minus
        asym = float("nan")
        err = float("nan")
        if total > 0:
            asym = (plus - minus) / total
            err = math.sqrt(max(0.0, 1.0 - asym * asym) / total)

        rows.append(
            {
                "model": model,
                "observable": observable,
                "bin_low_ps": low,
                "bin_high_ps": high,
                "bin_center_ps": 0.5 * (low + high),
                "count": total,
                plus_label: plus,
                minus_label: minus,
                "asymmetry": asym,
                "error": err,
            }
        )

    return rows


def write_signal_csv(path: Path, events: list[SignalEvent]) -> None:
    if not events:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(events[0]).keys()))
        writer.writeheader()
        for event in events:
            writer.writerow(asdict(event))


def write_rows_csv(path: Path, rows: list[dict]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_mixing_counts(path: Path, events_by_model: dict[str, list[SignalEvent]]) -> None:
    rows = []
    for model, events in events_by_model.items():
        total = len(events)
        for category in ("unmixed", "mixed", "weird"):
            count = sum(event.mixing_class == category for event in events)
            rows.append(
                {
                    "model": model,
                    "mixing_class": category,
                    "count": count,
                    "fraction": count / total if total else float("nan"),
                }
            )
    write_rows_csv(path, rows)


def mixing_closure_rows(
    rows_by_observable: dict[str, dict[str, list[dict]]],
    delta_m: float,
) -> list[dict]:
    rows = []
    for observable in ("initial_final", "record_mixing"):
        for model, binned_rows in rows_by_observable[observable].items():
            residuals = []
            pulls = []
            for row in binned_rows:
                asym = row["asymmetry"]
                err = row["error"]
                if not math.isfinite(asym):
                    continue
                expected = math.cos(delta_m * row["bin_center_ps"])
                residual = asym - expected
                residuals.append(residual)
                if math.isfinite(err) and err > 0:
                    pulls.append(residual / err)

            rows.append(
                {
                    "model": model,
                    "observable": observable,
                    "bins": len(residuals),
                    "rms_vs_cos": math.sqrt(np.mean(np.square(residuals))) if residuals else float("nan"),
                    "max_abs_vs_cos": max((abs(value) for value in residuals), default=float("nan")),
                    "pull_rms_vs_cos": math.sqrt(np.mean(np.square(pulls))) if pulls else float("nan"),
                }
            )
    return rows


def plot_asymmetry(
    rows_by_model: dict[str, list[dict]],
    observable: str,
    outpath: Path,
    delta_m: float,
    reference_sin_coeff: float,
) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 8.0))
    for model, rows in rows_by_model.items():
        centers = np.array([row["bin_center_ps"] for row in rows])
        asym = np.array([row["asymmetry"] for row in rows])
        err = np.array([row["error"] for row in rows])
        valid = np.isfinite(asym)
        ax.errorbar(centers[valid], asym[valid], yerr=err[valid], fmt="o", capsize=2, label=model)

    t_min = min(row["bin_low_ps"] for rows in rows_by_model.values() for row in rows)
    t_max = max(row["bin_high_ps"] for rows in rows_by_model.values() for row in rows)
    t_ref = np.linspace(t_min, t_max, 500)
    if observable == "cp_tag":
        ax.plot(
            t_ref,
            reference_sin_coeff * np.sin(delta_m * t_ref),
            "k--",
            linewidth=1.2,
            label=f"{reference_sin_coeff:.3f} sin(delta_m t)",
        )
        ylabel = "(N(initial B0) - N(initial anti-B0)) / sum"
    elif observable == "initial_final":
        ax.plot(t_ref, np.cos(delta_m * t_ref), "k--", linewidth=1.2, label="cos(delta_m t)")
        ylabel = "(N(same) - N(opposite)) / sum"
    elif observable == "record_mixing":
        ax.plot(t_ref, np.cos(delta_m * t_ref), "k--", linewidth=1.2, label="cos(delta_m t)")
        ylabel = "(N(unmixed) - N(mixed)) / sum"
    else:
        raise ValueError(f"unknown observable: {observable}")

    ax.axhline(0.0, color="0.4", linewidth=0.8)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel(ylabel)
    ax.set_ylim(-1.1, 1.1)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    add_cms_label(ax)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def plot_lifetime(events_by_model: dict[str, list[SignalEvent]], outpath: Path, t_min: float, t_max: float, bins: int) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 8.0))
    for model, events in events_by_model.items():
        t = [event.t_ps for event in events if t_min <= event.t_ps <= t_max]
        ax.hist(t, bins=bins, range=(t_min, t_max), histtype="step", linewidth=1.4, label=model)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("accepted signal candidates")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    add_cms_label(ax)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def plot_mixing_components(
    model: str,
    events: list[SignalEvent],
    outpath: Path,
    t_min: float,
    t_max: float,
    bins: int,
) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 8.0))
    edges = np.linspace(t_min, t_max, bins + 1)
    category_labels = {
        "unmixed": "Unmixed",
        "mixed": "Mixed",
        "weird": "Weird (B->B->X)",
    }

    all_values = [event.t_ps for event in events if t_min <= event.t_ps <= t_max]
    counts, _ = np.histogram(all_values, bins=edges)
    ax.hist(edges[:-1], bins=edges, weights=counts, histtype="step", linewidth=1.4, label="Sum")

    for category in ("unmixed", "mixed", "weird"):
        values = [
            event.t_ps
            for event in events
            if event.mixing_class == category and t_min <= event.t_ps <= t_max
        ]
        counts, _ = np.histogram(values, bins=edges)
        ax.hist(
            edges[:-1],
            bins=edges,
            weights=counts,
            histtype="step",
            linewidth=1.4,
            label=category_labels[category],
        )

    ax.set_yscale("log")
    ax.set_xlim(t_min, t_max)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("Events")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    events_by_model: dict[str, list[SignalEvent]] = {}
    summaries = []
    for model, path in args.sample:
        events, summary = read_hepmc(model, path, args.max_events)
        events_by_model[model] = events
        if events:
            initial = np.array([event.initial_pdg for event in events])
            final = np.array([event.final_pdg for event in events])
            times = np.array([event.t_ps for event in events])
            expected_mixed = np.sum(0.5 * (1.0 - np.cos(args.delta_m * times)))
            summary.update(
                {
                    "mean_t_ps": float(np.mean(times)),
                    "initial_B0_fraction": float(np.mean(initial == B0)),
                    "mixed_chain_fraction": float(np.mean([event.mixed_chain for event in events])),
                    "record_mixed_like_fraction": float(np.mean([event.mixing_class != "unmixed" for event in events])),
                    "initial_final_opposite_fraction": float(np.mean(initial != final)),
                    "expected_mixed_count_from_delta_m": float(expected_mixed),
                    "expected_mixed_fraction_from_delta_m": float(expected_mixed / len(events)),
                    "mixing_counts": {
                        category: int(sum(event.mixing_class == category for event in events))
                        for category in ("unmixed", "mixed", "weird")
                    },
                }
            )
        summaries.append(summary)
        write_signal_csv(args.outdir / f"{model}_signal_events.csv", events)

    bins = np.linspace(args.t_min, args.t_max, args.bins + 1)
    all_rows: list[dict] = []
    rows_by_observable: dict[str, dict[str, list[dict]]] = {
        "cp_tag": {},
        "initial_final": {},
        "record_mixing": {},
    }
    for model, events in events_by_model.items():
        for observable in rows_by_observable:
            rows = asymmetry_rows(model, events, bins, observable)
            rows_by_observable[observable][model] = rows
            all_rows.extend(rows)

    write_rows_csv(args.outdir / "asymmetry_bins.csv", all_rows)
    write_mixing_counts(args.outdir / "mixing_counts.csv", events_by_model)
    closure_rows = mixing_closure_rows(rows_by_observable, args.delta_m)
    write_rows_csv(args.outdir / "mixing_closure.csv", closure_rows)
    summary_by_model = {summary["model"]: summary for summary in summaries}
    for row in closure_rows:
        summary = summary_by_model[row["model"]]
        prefix = f"{row['observable']}_closure"
        summary[f"{prefix}_rms_vs_cos"] = row["rms_vs_cos"]
        summary[f"{prefix}_max_abs_vs_cos"] = row["max_abs_vs_cos"]

    plot_asymmetry(
        rows_by_observable["cp_tag"],
        "cp_tag",
        args.outdir / "cp_tag_asymmetry.png",
        args.delta_m,
        args.reference_sin_coeff,
    )
    plot_asymmetry(
        rows_by_observable["initial_final"],
        "initial_final",
        args.outdir / "initial_final_flavor_asymmetry.png",
        args.delta_m,
        args.reference_sin_coeff,
    )
    plot_asymmetry(
        rows_by_observable["record_mixing"],
        "record_mixing",
        args.outdir / "record_mixing_asymmetry.png",
        args.delta_m,
        args.reference_sin_coeff,
    )
    plot_lifetime(events_by_model, args.outdir / "proper_time_distribution.png", args.t_min, args.t_max, args.bins * 4)
    for model, events in events_by_model.items():
        plot_mixing_components(
            model,
            events,
            args.outdir / f"{model}_mixing_components.png",
            args.t_min,
            args.t_max,
            args.bins * 4,
        )

    with (args.outdir / "summary.json").open("w") as handle:
        json.dump(summaries, handle, indent=2)

    print(f"Wrote results to {args.outdir}")
    for summary in summaries:
        print(
            f"{summary['model']}: signals={summary['signals']} "
            f"events_seen={summary['events_seen']} "
            f"initial_B0_fraction={summary.get('initial_B0_fraction', float('nan')):.4f} "
            f"mixed_chain_fraction={summary.get('mixed_chain_fraction', float('nan')):.4f} "
            f"mixing_counts={summary.get('mixing_counts', {})}"
        )
    for row in closure_rows:
        print(
            f"{row['model']} {row['observable']}: "
            f"rms_vs_cos={row['rms_vs_cos']:.4f} "
            f"max_abs_vs_cos={row['max_abs_vs_cos']:.4f}"
        )


if __name__ == "__main__":
    main()
