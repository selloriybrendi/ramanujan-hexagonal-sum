"""LN IV (9.2.5): Ramanujan yozgan ayniyatdagi aniqlanmagan K, λ (Andrews–Berndt: "unable either to identify them").
Gipoteza (Eyzenshteyn panjarasining 6-karrali simmetriyasidan): K = μ−ν, λ = μ+ν, μ,ν ≥ 1, (μ,ν)=1.

Ramanujan (9.2.5):  σ_{s-1}(1)e^{-π√3} − σ_{s-1}(2)e^{-2π√3} + ...
   = −2Γ(s)/((2π)^s ζ(s)) · { cos(πs/6) + 2cos(π(s+1)/6)cos(π(s−1)/6) Σ cos(s·atan(K/(λ√3)))/(μ²+μν+ν²)^{s/2} }
Nazorat: (a) Andrews–Berndt isbotlagan (9.3.3)-dan keyingi shakl; (b) K,λ ning boshqa "tabiiy" tanlovlari.
"""
import math
import numpy as np
import mpmath as mp

mp.mp.dps = 30


def chap(s, nmax=60):
    # σ_{s-1}(n) e^{-nπ√3}, ishora (+,−,+,...) — e^{-π√3}≈0.0043, 60 had yetarli
    jami = mp.mpf(0)
    for n in range(1, nmax + 1):
        sig = sum(mp.mpf(d) ** (s - 1) for d in range(1, n + 1) if n % d == 0)
        jami += (-1) ** (n + 1) * sig * mp.e ** (-n * mp.pi * mp.sqrt(3))
    return jami


def panjara(N):
    m = np.arange(1, N + 1)
    M, V = np.meshgrid(m, m, indexing="ij")
    g = np.gcd(M, V) == 1
    return M[g].astype(float), V[g].astype(float)


def ong(s, K, L, M, V):
    Q = M * M + M * V + V * V
    S = np.sum(np.cos(s * np.arctan(K / (L * math.sqrt(3)))) / Q ** (s / 2))
    c = math.cos
    qavs = c(math.pi * s / 6) + 2 * c(math.pi * (s + 1) / 6) * c(math.pi * (s - 1) / 6) * S
    return -2 * math.gamma(s) * float(mp.zeta(s)) / (2 * math.pi) ** s * qavs   # ζ(s) SURATDA (OCR 'Γ(s) (2π)sζ(s)' ni noto'g'ri o'qigan edim)


def andrews_berndt(s, N):
    # Σ_{μ≥1, ν∈Z, (μ,ν)=1} cos(s·atan((μ+2ν)/(√3μ)))/(μ²+μν+ν²)^{s/2}; ular LHS'ining ishorasi (−1)^n → bizniki −(...)
    mu = np.arange(1, N + 1)
    nu = np.arange(-N, N + 1)
    M, V = np.meshgrid(mu, nu, indexing="ij")
    g = np.gcd(M, np.abs(V)) == 1
    M, V = M[g].astype(float), V[g].astype(float)
    Q = M * M + M * V + V * V
    S = np.sum(np.cos(s * np.arctan((M + 2 * V) / (math.sqrt(3) * M))) / Q ** (s / 2))
    return -math.gamma(s) * float(mp.zeta(s)) / (2 * math.pi) ** s * S


N = 2500
M, V = panjara(N)
print(f"panjara N={N}: {M.size} ta o'zaro tub juftlik")
tanlov = {
    "K=μ−ν, λ=μ+ν (GIPOTEZA)": (M - V, M + V),
    "K=μ+2ν, λ=μ   (A–B burchagi)": (M + 2 * V, M),
    "K=μ−ν, λ=μ    (nazorat)": (M - V, M),
    "K=ν,   λ=μ    (nazorat)": (V, M),
}
for s in (4.0, 6.0, 7.5, 10.0):
    L = float(chap(s))
    ab = andrews_berndt(s, 1200)
    print(f"\ns = {s}:  CHAP = {L:.15e}")
    print(f"   Andrews–Berndt (isbotlangan) shakl: {ab:.15e}  nisbiy farq {abs(ab - L) / abs(L):.2e}")
    for nom, (K, Lm) in tanlov.items():
        R = ong(s, K, Lm, M, V)
        print(f"   {nom:32s}: {R:.15e}  nisbiy farq {abs(R - L) / abs(L):.2e}")
