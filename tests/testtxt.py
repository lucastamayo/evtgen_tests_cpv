#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

# =========================
# INPUT
# =========================
if len(sys.argv) < 2:
    print("Usage: python analyze_ssdcp_correct.py ssdcp_debug.txt")
    sys.exit(1)

fname = sys.argv[1]

# =========================
# CONFIG
# =========================
NBINS = 10      # menos bins = mais estabilidade
TMIN = 0
TMAX = 5

DELTA_M = 0.5065  # ps^-1 (Bd)

# =========================
# LOAD DATA
# =========================
print(f"Reading file: {fname}")

data = np.loadtxt(fname)

t     = data[:,0]
pid   = data[:,1].astype(int)
other = data[:,2].astype(int)
flip  = data[:,3].astype(int)

print("\n===== BASIC INFO =====")
print("Total events:", len(t))

# =========================
# SANITY CHECKS
# =========================
print("\n===== SANITY CHECKS =====")

print("\nPID:")
print("B0       :", np.mean(pid == 511))
print("anti-B0  :", np.mean(pid == -511))

print("\nother_b:")
print("B0       :", np.mean(other == 511))
print("anti-B0  :", np.mean(other == -511))

print("\nflip:")
print("flip=1   :", np.mean(flip == 1))

# =========================
# HIST t
# =========================
plt.figure(figsize=(6,5))
plt.hist(t, bins=50)
plt.xlabel("t [ps]")
plt.ylabel("Entries")
plt.title("Proper time distribution")
plt.tight_layout()
plt.savefig("t_distribution.png")
print("Saved: t_distribution.png")

# =========================
# ASYMMETRY (CORRETA)
# =========================

bins = np.linspace(TMIN, TMAX, NBINS+1)
centers = 0.5 * (bins[:-1] + bins[1:])

A = np.zeros_like(centers)
err = np.zeros_like(centers)

for i in range(len(centers)):

    mask = (t >= bins[i]) & (t < bins[i+1])

    # SAME vs OPPOSITE flavor
    same = np.sum((pid == other) & mask)
    opp  = np.sum((pid != other) & mask)

    N = same + opp

    if N > 0:
        A[i] = (same - opp) / N
        err[i] = np.sqrt((1 - A[i]**2) / N)
    else:
        A[i] = np.nan
        err[i] = np.nan

# =========================
# TEORIA
# =========================
t_th = np.linspace(TMIN, TMAX, 500)
A_th = np.cos(DELTA_M * t_th)

# =========================
# PLOT
# =========================
plt.figure(figsize=(6,5))

mask = np.isfinite(A)

plt.errorbar(centers[mask], A[mask], yerr=err[mask], fmt='o', label="EvtGen")

plt.plot(t_th, A_th, '-', label=r'$\cos(\Delta m t)$')

plt.axhline(0)
plt.xlabel("t [ps]")
plt.ylabel("A(t)")
plt.title("Asymmetry from SSD_CP (correct observable)")
plt.legend()
plt.ylim(-1.1, 1.1)

plt.tight_layout()
plt.savefig("asymmetry_correct.png")
print("Saved: asymmetry_correct.png")

print("\n===== DONE =====")