"""#43 dan qurol: olti burchakli (Eyzenshteyn) panjara yig'indisi
    R(s) = Σ_{μ,ν≥1,(μ,ν)=1} cos(s·tan⁻¹((μ−ν)/((μ+ν)√3))) / (μ²+μν+ν²)^{s/2}
ni to'liq (9.2.5) orqali q-qatordan tez hisoblash:
    R(s) = [ −L(s)·(2π)^s / (2Γ(s)ζ(s)) − cos(πs/6) ] / (2cos(π(s+1)/6)cos(π(s−1)/6)),
    L(s) = Σ (−1)^{n+1} σ_{s−1}(n) e^{−nπ√3}   (e^{−π√3} ≈ 0.0043 — har had ~2.4 xona beradi).
Maxraj 0 bo'lgan s ≡ 2,4 (mod 6) da ayniyat R haqida ma'lumot bermaydi.
Nazorat: to'g'ridan-to'g'ri 2D yig'indi (N gacha) — s ga yaqin 2 da juda sekin (xato ~N^{2−s}).
"""
import math
import time
import numpy as np
import mpmath as mp


def L(s, dps=40):
    mp.mp.dps = dps
    s = mp.mpf(s)
    q = mp.e ** (-mp.pi * mp.sqrt(3))
    jami, n = mp.mpf(0), 1
    while True:
        sig = mp.fsum(mp.mpf(d) ** (s - 1) for d in range(1, n + 1) if n % d == 0)
        had = sig * q ** n
        jami += (-1) ** (n + 1) * had
        if had < mp.mpf(10) ** (-dps - 5):
            return jami, n
        n += 1


def R_tez(s, dps=40):
    mp.mp.dps = dps
    Ls, n = L(s, dps)
    s = mp.mpf(s)
    maxraj = 2 * mp.cos(mp.pi * (s + 1) / 6) * mp.cos(mp.pi * (s - 1) / 6)
    if abs(maxraj) < mp.mpf(10) ** (-dps // 2):
        return None, n
    return (-Ls * (2 * mp.pi) ** s / (2 * mp.gamma(s) * mp.zeta(s)) - mp.cos(mp.pi * s / 6)) / maxraj, n


def R_brute(s, N):
    m = np.arange(1, N + 1)
    M, V = np.meshgrid(m, m, indexing="ij")
    g = np.gcd(M, V) == 1
    M, V = M[g].astype(float), V[g].astype(float)
    return float(np.sum(np.cos(s * np.arctan((M - V) / ((M + V) * math.sqrt(3)))) / (M * M + M * V + V * V) ** (s / 2)))


print("=== 1) Tezlik va aniqlik: q-qator (9.2.5) vs to'g'ridan-to'g'ri 2D yig'indi")
for s in (2.5, 3.0, 3.3, 5.5, 7.5):
    t0 = time.time()
    r, n = R_tez(s)
    t1 = time.time()
    qator = []
    for N in (500, 2000, 4000):
        t2 = time.time()
        b = R_brute(s, N)
        t3 = time.time()
        qator.append(f"N={N}: farq {abs(b - float(r)) / abs(float(r)):.1e} ({t3 - t2:.1f}s)")
    print(f" s={s}: R = {mp.nstr(r, 30)}  [q-qator {n} had, {1000 * (t1 - t0):.0f} ms]")
    print("        brute: " + " · ".join(qator))

print("\n=== 2) Yopiq formulalar: s=6 va s=12 da L(s) = E_s(ρ) orqali (ρ = e^{2πi/3}, E₄(ρ)=0)")
mp.mp.dps = 40
G = mp.gamma(mp.mpf(1) / 3)
tau = mp.mpc(-0.5, mp.sqrt(3) / 2)             # ρ
q = mp.exp(2j * mp.pi * tau)                     # = −e^{−π√3}
E6 = 1 - 504 * mp.nsum(lambda n: n ** 5 * q ** n / (1 - q ** n), [1, mp.inf])
print(" E6(ρ) =", mp.nstr(E6.real, 30), "| E6(ρ) - 27Γ(1/3)^18/(512π^12) =", mp.nstr(E6.real - 27 * G ** 18 / (512 * mp.pi ** 12), 3))
for s in (6, 12):
    r, _ = R_tez(s)
    closed = 1 - 9 * G ** 18 / (512 * mp.pi ** 12) if s == 6 else 30375 * G ** 36 / (90570752 * mp.pi ** 24) - 1
    print(f" R({s}) = {mp.nstr(r, 30)} · closed form = {mp.nstr(closed, 30)} · difference = {mp.nstr(abs(r - closed), 3)}")
