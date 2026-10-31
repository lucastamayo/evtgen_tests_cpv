#!/usr/bin/env python3
import math
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

try:
    import pyhepmc_ng as pyhepmc
except ImportError:
    pyhepmc = None

try:
    import mplhep as hep

    plt.style.use(hep.style.CMS)
except ImportError:
    hep = None

MB = 5.27963
C_MM_PER_PS = 0.299792458

# cortes numéricos
MIN_L = 1e-6
MIN_P = 1e-6
MIN_T = 1e-6


# ----------------------------
# Helpers
# ----------------------------
def _val(x):
    return x() if callable(x) else x

def fourvec_x(v):
    return _val(getattr(v, "x"))

def fourvec_y(v):
    return _val(getattr(v, "y"))

def fourvec_z(v):
    return _val(getattr(v, "z"))

def particle_pid(p):
    return _val(getattr(p, "pid"))

def particle_momentum(p):
    return p.momentum

def particle_px(p):
    return fourvec_x(particle_momentum(p))

def particle_py(p):
    return fourvec_y(particle_momentum(p))

def particle_pz(p):
    return fourvec_z(particle_momentum(p))

def particle_p(p):
    px, py, pz = particle_px(p), particle_py(p), particle_pz(p)
    return math.sqrt(px * px + py * py + pz * pz)

def prod_vertex(p):
    return p.production_vertex

def end_vertex(p):
    return p.end_vertex

def vertex_pos(v):
    return v.position

def vertex_x(v):
    return fourvec_x(vertex_pos(v))

def vertex_y(v):
    return fourvec_y(vertex_pos(v))

def vertex_z(v):
    return fourvec_z(vertex_pos(v))

def children(p):
    ev = end_vertex(p)
    if ev and hasattr(ev, "particles_out"):
        return list(ev.particles_out)
    return []

def parents(p):
    pv = prod_vertex(p)
    if pv and hasattr(pv, "particles_in"):
        return list(pv.particles_in)
    return []



def has_desc(p, target_abs, depth=0, maxd=20):
    if depth > maxd:
        return False

    for d in children(p):
        if abs(particle_pid(d)) == target_abs:
            return True
        if has_desc(d, target_abs, depth + 1, maxd):
            return True

    return False


def terminal_bd(b):
    """
    Pega o último B0/B0bar da cadeia que realmente decai.
    """
    cur = b
    for _ in range(100):
        bd_daus = [d for d in children(cur) if abs(particle_pid(d)) == 511]

        if not bd_daus:
            return cur

        chosen = None
        for d in bd_daus:
            if end_vertex(d) is not None and len(children(d)) > 0:
                chosen = d
                break

        cur = chosen if chosen else bd_daus[0]

    return cur


def proper_time_ps(b_prod, b_decay):
    v0 = prod_vertex(b_prod)
    v1 = end_vertex(b_decay)

    if v0 is None or v1 is None:
        return None

    daus = children(b_decay)
    if len(daus) == 0:
        return None

    p = particle_p(b_decay)
    if p < MIN_P:
        return None

    dx = vertex_x(v1) - vertex_x(v0)
    dy = vertex_y(v1) - vertex_y(v0)
    dz = vertex_z(v1) - vertex_z(v0)

    L_vec = math.sqrt(dx * dx + dy * dy + dz * dz)
    if L_vec < MIN_L:
        return None

    px = particle_px(b_decay)
    py = particle_py(b_decay)
    pz = particle_pz(b_decay)

    L = (dx * px + dy * py + dz * pz) / p

    t = (MB / p) * abs(L) / C_MM_PER_PS

    if t < MIN_T:
        return None

    return t


def classify_candidate(p_term):
    moms = [m for m in parents(p_term) if abs(particle_pid(m)) == 511]

    if len(moms) == 0:
        return "unmixed"

    m = moms[0]
    if particle_pid(m) == -particle_pid(p_term):
        return "mixed"

    return "weird"


def is_root_bd(p):
    return len([m for m in parents(p) if abs(particle_pid(m)) == 511]) == 0


# ----------------------------
# Debug
# ----------------------------
def debug_candidate(p, p_term, t):
    print("\n=== DEBUG ===")
    print(f"Initial B: {particle_pid(p)}")
    print(f"Final   B: {particle_pid(p_term)}")

    v0 = prod_vertex(p)
    v1 = end_vertex(p_term)

    if v0 is not None:
        print(f"v_prod = ({vertex_x(v0):.4e}, {vertex_y(v0):.4e}, {vertex_z(v0):.4e})")
    if v1 is not None:
        print(f"v_dec  = ({vertex_x(v1):.4e}, {vertex_y(v1):.4e}, {vertex_z(v1):.4e})")

    print(f"t = {t:.6f} ps")


# ----------------------------
# Main
# ----------------------------
def main(fname, outpng, nbins, tmin, tmax, max_events, debug_events):
    bins = np.linspace(tmin, tmax, nbins + 1)

    t_mixed = []
    t_unmixed = []
    t_weird = []
    events_seen = 0

    if pyhepmc is None:
        from cp_study.hepmc_cp_asymmetry import read_hepmc

        events, summary = read_hepmc("sample", Path(fname), max_events)
        events_seen = summary["events_seen"]
        for event in events:
            if event.mixing_class == "mixed":
                t_mixed.append(event.t_ps)
            elif event.mixing_class == "unmixed":
                t_unmixed.append(event.t_ps)
            else:
                t_weird.append(event.t_ps)
    else:
        with pyhepmc.open(fname) as f:
            for ievt, event in enumerate(f):
                if max_events is not None and ievt >= max_events:
                    break
                events_seen += 1
                seen_terminals = set()
                for p in event.particles:

                    if abs(particle_pid(p)) != 511:
                        continue

                    if not is_root_bd(p):
                        continue

                    if not (has_desc(p, 443) and has_desc(p, 310)):
                        continue

                    p_term = terminal_bd(p)
                    terminal_key = id(p_term)
                    if terminal_key in seen_terminals:
                        continue
                    seen_terminals.add(terminal_key)

                    if len(children(p_term)) == 0:
                        continue

                    t = proper_time_ps(p, p_term)
                    if t is None:
                        continue

                    cls = classify_candidate(p_term)

                    if ievt < debug_events:
                        debug_candidate(p, p_term, t)

                    if cls == "mixed":
                        t_mixed.append(t)
                    elif cls == "unmixed":
                        t_unmixed.append(t)
                    else:
                        t_weird.append(t)

    hm, _ = np.histogram(t_mixed, bins=bins)
    hu, _ = np.histogram(t_unmixed, bins=bins)
    hw, _ = np.histogram(t_weird, bins=bins)
    hs = hm + hu + hw

    fig, ax = plt.subplots(figsize=(8, 8))

    ax.hist(bins[:-1], bins=bins, weights=hs, histtype="step", linewidth=1.8, label="Sum")
    ax.hist(bins[:-1], bins=bins, weights=hm, histtype="step", linewidth=1.8, label="Mixed")
    ax.hist(bins[:-1], bins=bins, weights=hu, histtype="step", linewidth=1.8, label="Unmixed")
    ax.hist(bins[:-1], bins=bins, weights=hw, histtype="step", linewidth=1.8, label="Weird (B->B->X)")

    ax.set_yscale("log")
    ax.set_xlim(tmin, tmax)
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("Events")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()

    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(outpng, dpi=150)
    # plt.show()

    print("\n================ SUMMARY ================")
    print(f"Input   : {fname}")
    print(f"Events  : {events_seen}")
    print(f"Mixed   : {len(t_mixed)}")
    print(f"Unmixed : {len(t_unmixed)}")
    print(f"Weird   : {len(t_weird)}")
    print(f"Saved   : {outpng}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot proper-time distributions for mixed/unmixed B0 -> J/psi Ks candidates from HepMC."
    )
    parser.add_argument("input", help="Input HepMC file")
    parser.add_argument("-o", "--output", default="mixing_time.png", help="Output image name")
    parser.add_argument("--nbins", type=int, default=50, help="Number of bins")
    parser.add_argument("--tmin", type=float, default=0.0, help="Minimum t (ps)")
    parser.add_argument("--tmax", type=float, default=10.0, help="Maximum t (ps)")
    parser.add_argument("--max-events", "--max", dest="max_events", type=int, default=None, help="Maximum number of events to read")
    parser.add_argument("--debug-events", type=int, default=3, help="Number of initial events to debug-print")

    args = parser.parse_args()

    main(
        fname=args.input,
        outpng=args.output,
        nbins=args.nbins,
        tmin=args.tmin,
        tmax=args.tmax,
        max_events=args.max_events,
        debug_events=args.debug_events,
    )
