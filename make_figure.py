#!/usr/bin/env python3
"""Figure for the note on (9.2.5).  Output: fig/hex_925.png, fig/hex_925.pdf
Left: primitive points of Z[omega]; the half-plane H = {mu >= 1}, the sector a, b >= 1, and the six associates of w = 2 + omega.
Right: 3^{s/2} R(s) from Theorem 2 (fast evaluation); dots: closed forms R(6k); dashed: s = 2, 4 (mod 6), where Theorem 1 says
nothing about R(s) (there R(s) is taken from the lattice sum)."""
import math
import os

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42  # TrueType (arXiv flags Type 3 fonts)
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt  # noqa: E402
import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402

os.makedirs("fig", exist_ok=True)
mp.mp.dps = 30
q0 = mp.e ** (-mp.pi * mp.sqrt(3))
SQ3 = math.sqrt(3)


def Lq(s):
    tot, n = mp.mpf(0), 1
    while True:
        term = mp.fsum(mp.mpf(d) ** (s - 1) for d in range(1, n + 1) if n % d == 0) * q0 ** n
        tot += (-1) ** (n + 1) * term
        if term < mp.mpf(10) ** -32:
            return tot
        n += 1


def R_fast(s):
    s = mp.mpf(s)
    den = 2 * mp.cos(mp.pi * (s + 1) / 6) * mp.cos(mp.pi * (s - 1) / 6)
    return (-Lq(s) * (2 * mp.pi) ** s / (2 * mp.gamma(s) * mp.zeta(s)) - mp.cos(mp.pi * s / 6)) / den


mm = np.arange(1, 1501)
A, B = np.meshgrid(mm, mm, indexing="ij")
g = np.gcd(A, B) == 1
A, B = A[g].astype(float), B[g].astype(float)
QN, PSI = A * A + A * B + B * B, np.arctan((B - A) / (SQ3 * (A + B)))


def R_lattice(s):
    return float(np.sum(np.cos(s * PSI) / QN ** (s / 2)))


svals = np.linspace(3.0, 24.0, 1400)
Rv = []
for s in svals:
    near_degenerate = min(abs(s - d) for d in (4, 8, 10, 14, 16, 20, 22)) < 0.08
    Rv.append(R_lattice(s) if near_degenerate else float(R_fast(s)))
Rv = np.array(Rv)
G = mp.gamma(mp.mpf(1) / 3)
E6 = 27 * G ** 18 / (512 * mp.pi ** 12)
C = [mp.mpf(1), mp.mpf(250) / 691, mp.mpf(5500) / 43867, mp.mpf(10285000) / 236364091]
R6k = [float((-1) ** k * (C[k - 1] * E6 ** k - 3) / 3) for k in range(1, 5)]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1, 1.25]})

# left panel
w = complex(0.5, SQ3 / 2)
pts = [(a, b) for a in range(-6, 7) for b in range(-6, 7) if (a, b) != (0, 0) and math.gcd(a, b) == 1 and abs(a + b * w) <= 5.6]
z = np.array([a + b * w for a, b in pts])
ax1.fill([0, 6, 6, 5.6 / SQ3], [0, 0, 5.6, 5.6], color="#1f77b4", alpha=0.12, lw=0)  # open sector 0 < arg < pi/3
t = np.linspace(-2 * math.pi / 3, math.pi / 3, 2)
ax1.plot(6.2 * np.cos(t), 6.2 * np.sin(t), color="#555", lw=0.8, ls="--")
ax1.scatter(z.real, z.imag, s=7, color="#999", zorder=2)
orb = [(2, 1)]
for _ in range(5):
    a, b = orb[-1]
    orb.append((-b, a + b))
for j, (a, b) in enumerate(orb):
    p = a + b * w
    inH = a >= 1
    ax1.scatter([p.real], [p.imag], s=46, zorder=3, color="#d62728" if inH else "white", edgecolor="#d62728", lw=1.2)
    lab = ["$w$", r"$\omega w$", r"$\omega^2 w$", r"$-w$", r"$\omega^{-2}w$", r"$\omega^{-1}w$"][j]
    ax1.annotate(lab, (p.real, p.imag), xytext=(6, 4), textcoords="offset points", fontsize=9, color="#d62728")
ax1.text(3.3, 1.2, r"$a,b\geq 1$", color="#1f77b4", fontsize=10)
ax1.text(-4.3, -5.3, r"$\mu=0$", color="#555", fontsize=9)
ax1.text(2.2, -4.6, r"$H=\{\mu\geq1\}$", color="#555", fontsize=9)
ax1.set_aspect("equal")
ax1.set_xlim(-6, 6)
ax1.set_ylim(-5.6, 5.6)
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_title(r"Lemma 1: the orbit of $w=2+\omega$ meets $H$ three times (filled)", fontsize=10)

# right panel
ax2.plot(svals, 3 ** (svals / 2) * Rv, color="#1f77b4", lw=1.3, label=r"$3^{s/2}R(s)$ (Theorem 2)")
ks = np.array([6, 12, 18, 24])
ax2.scatter(ks, 3 ** (ks / 2) * np.array(R6k), color="#d62728", zorder=3, s=30,
            label=r"$R(6k)=(-1)^k(c_kE_6(\rho)^k-3)/3$")
for d in (4, 8, 10, 14, 16, 20, 22):
    ax2.axvline(d, color="#999", lw=0.6, ls=":")
ax2.axhline(1, color="#555", lw=0.6)
ax2.set_xlabel("s")
ax2.set_title(r"dotted: $s\equiv2,4\ (\mathrm{mod}\ 6)$, where the coefficient of $R(s)$ vanishes", fontsize=10)
ax2.legend(fontsize=8.5, frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig("fig/hex_925.png", dpi=200)
fig.savefig("fig/hex_925.pdf")
fig.savefig("paper/hex_925.pdf")  # flat copy for the arXiv source

# consistency of the plotted data: lattice sum vs Theorem 2 away from degenerate s, and dots on the curve
far = (5.0, 7.0, 11.5, 17.3)
worst = max(abs(R_lattice(s) - float(R_fast(s))) / abs(float(R_fast(s))) for s in far)
worst_dot = max(abs(float(R_fast(6 * k)) - R6k[k - 1]) / abs(R6k[k - 1]) for k in range(1, 5))
print(f"lattice (N=1500) vs Theorem 2 at s = 5, 7, 11.5, 17.3: max rel. diff {worst:.1e}")
print(f"closed-form dots vs Theorem 2 at s = 6, 12, 18, 24: max rel. diff {worst_dot:.1e}")
print(f"3^(s/2) R(s) range on [3, 24]: {(3 ** (svals / 2) * Rv).min():.4f} .. {(3 ** (svals / 2) * Rv).max():.4f}")
