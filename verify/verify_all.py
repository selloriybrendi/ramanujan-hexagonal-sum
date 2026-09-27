#!/usr/bin/env python3
"""verify_all.py — one-command check of the note on (9.2.5) of Ramanujan's Lost Notebook (Andrews–Berndt, Part IV, Sect. 9.2).
Requires numpy, mpmath, sympy.   Run: python3 verify_all.py   (exit code 0 = all checks pass)

Claim: for Re s > 2,  sum (-1)^{n+1} sigma_{s-1}(n) e^{-n pi sqrt3}
   = -2 Gamma(s) zeta(s)/(2 pi)^s * { cos(pi s/6) + 2 cos(pi(s+1)/6) cos(pi(s-1)/6) R(s) },
   R(s) = sum_{mu,nu>=1, gcd=1} cos(s * atan((mu-nu)/((mu+nu) sqrt3))) / (mu^2+mu nu+nu^2)^{s/2},  i.e. K = mu-nu, lambda = mu+nu.
"""
import math
import sys

import mpmath as mp
import numpy as np
import sympy as sp

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  —  {detail}" if detail else ""))


# ---------- 1. symbolic identities used in the proof ----------
a, b, s, psi = sp.symbols("a b s psi", positive=True)
w = sp.expand(a + b * sp.exp(sp.I * sp.pi / 3))
theta = sp.atan2(sp.im(w), sp.re(w))
e1 = sp.simplify(sp.expand_trig(sp.tan(theta - sp.pi / 6)) - (b - a) / (sp.sqrt(3) * (a + b)))
L3 = sp.cos(s * (psi + sp.pi / 3)) + sp.cos(s * psi) + sp.cos(s * (psi - sp.pi / 3))
e2 = sp.simplify(sp.expand_trig(sp.expand(L3 - sp.cos(s * psi) * (1 + 2 * sp.cos(sp.pi * s / 3)))))
e3 = sp.simplify(sp.expand_trig(1 + 2 * sp.cos(sp.pi * s / 3) - 4 * sp.cos(sp.pi * (s + 1) / 6) * sp.cos(sp.pi * (s - 1) / 6)))
e4 = sp.nsimplify(sp.expand_complex(sp.expand(w * sp.conjugate(w)) - (a * a + a * b + b * b)))
check("1 symbolic: tan(arg w - pi/6), three rotations, product formula, |w|^2 — all differences simplify to 0",
      e1 == 0 and e2 == 0 and e3 == 0 and e4 == 0)


# ---------- 2. orbit lemma ----------
def associates(x, y):  # w = x + y omega, omega*w = -y + (x+y) omega
    out = []
    for _ in range(6):
        out.append((x, y))
        x, y = -y, x + y
    return out


bad, total, unit_k = 0, 0, set()
for x in range(-120, 121):
    for y in range(-120, 121):
        if (x, y) == (0, 0) or math.gcd(x, y) != 1:
            continue
        orb = associates(x, y)
        k = sum(1 for (m, _) in orb if m >= 1)
        if x * x + x * y + y * y == 1:
            unit_k.add(k)
        elif k != 3 or any(m == 0 for (m, _) in orb):
            bad += 1
        total += 1
rep_bad = sum(1 for x in range(-40, 41) for y in range(-40, 41)
              if (x, y) != (0, 0) and math.gcd(x, y) == 1 and x * x + x * y + y * y > 1
              and sum(1 for (m, n) in associates(x, y) if m >= 1 and n >= 1) != 1)
check("2 orbit lemma: every non-unit primitive orbit has exactly 3 elements with mu >= 1 and exactly 1 with mu,nu >= 1; unit orbit has 2",
      bad == 0 and rep_bad == 0 and unit_k == {2}, f"{total} primitive w checked")

# ---------- 3. the identity with K = mu - nu, lambda = mu + nu, against controls ----------
mp.mp.dps = 30
q0 = mp.e ** (-mp.pi * mp.sqrt(3))


def left(sv, nmax=60):
    return sum((-1) ** (n + 1) * sum(mp.mpf(d) ** (sv - 1) for d in range(1, n + 1) if n % d == 0) * q0 ** n for n in range(1, nmax + 1))


N = 2500
m_ = np.arange(1, N + 1)
M, V = np.meshgrid(m_, m_, indexing="ij")
g = np.gcd(M, V) == 1
M, V = M[g].astype(float), V[g].astype(float)
Q = M * M + M * V + V * V


def right(sv, K, Lm):
    S = np.sum(np.cos(sv * np.arctan(K / (Lm * math.sqrt(3)))) / Q ** (sv / 2))
    c = math.cos
    return -2 * math.gamma(sv) * float(mp.zeta(sv)) / (2 * math.pi) ** sv * (c(math.pi * sv / 6) + 2 * c(math.pi * (sv + 1) / 6) * c(math.pi * (sv - 1) / 6) * S)


svals = (4.7, 5.5, 6.0, 7.5, 8.7, 9.25, 13.1)
rel_main, rel_ctrl, rel_s6 = [], [], None
for sv in svals:
    Lv = float(left(sv))
    rel_main.append(abs(right(sv, M - V, M + V) - Lv) / abs(Lv))
    for name, K, Lm in (("mu+2nu,mu", M + 2 * V, M), ("mu-nu,mu", M - V, M), ("nu,mu", V, M)):
        r = abs(right(sv, K, Lm) - Lv) / abs(Lv)
        if sv == 6.0 and name == "mu+2nu,mu":
            rel_s6 = r  # its angles are pi/3 - psi_w, and cos(6(pi/3 - psi)) = cos(6 psi)
        else:
            rel_ctrl.append(r)
check("3 identity holds with K = mu-nu, lambda = mu+nu at s = 4.7, 5.5, 6, 7.5, 8.7, 9.25, 13.1 (N = 2500)",
      max(rel_main) < 1e-11, f"{M.size} coprime pairs, rel. diff {min(rel_main):.1e} .. {max(rel_main):.1e}")
check("4 other natural choices of K, lambda fail (relative miss >= 5e-6)", min(rel_ctrl) >= 5e-6,
      f"miss {min(rel_ctrl):.1e} .. {max(rel_ctrl):.1e}")
refl = max(abs(math.atan((m + 2 * v) / (m * math.sqrt(3))) + math.atan((m - v) / ((m + v) * math.sqrt(3))) - math.pi / 3)
           for m in range(1, 200) for v in range(1, 200))
check("4b at s = 6 the choice K = mu+2nu, lambda = mu also agrees: its angles are pi/3 - psi_w and cos(6(pi/3 - psi)) = cos(6 psi)",
      rel_s6 < 1e-11 and refl < 1e-12, f"rel. diff {rel_s6:.1e}, max |angle - (pi/3 - psi_w)| {refl:.1e}")

# ---------- 5. fast evaluation of R(s) and closed forms ----------
mp.mp.dps = 40
q0 = mp.e ** (-mp.pi * mp.sqrt(3))


def Lq(sv):
    sv = mp.mpf(sv)
    tot, n = mp.mpf(0), 1
    while True:
        term = mp.fsum(mp.mpf(d) ** (sv - 1) for d in range(1, n + 1) if n % d == 0) * q0 ** n
        tot += (-1) ** (n + 1) * term
        if term < mp.mpf(10) ** -45:
            return tot
        n += 1


def R_fast(sv):
    sv = mp.mpf(sv)
    den = 2 * mp.cos(mp.pi * (sv + 1) / 6) * mp.cos(mp.pi * (sv - 1) / 6)
    return (-Lq(sv) * (2 * mp.pi) ** sv / (2 * mp.gamma(sv) * mp.zeta(sv)) - mp.cos(mp.pi * sv / 6)) / den


def R_brute(sv, n_):
    mm = np.arange(1, n_ + 1)
    A_, B_ = np.meshgrid(mm, mm, indexing="ij")
    gg = np.gcd(A_, B_) == 1
    A_, B_ = A_[gg].astype(float), B_[gg].astype(float)
    return float(np.sum(np.cos(sv * np.arctan((A_ - B_) / ((A_ + B_) * math.sqrt(3)))) / (A_ * A_ + A_ * B_ + B_ * B_) ** (sv / 2)))


rb = R_brute(5.5, 2000)
rf = float(R_fast(5.5))
check("5 fast q-series R(5.5) agrees with the direct lattice sum (N = 2000)", abs(rb - rf) / abs(rf) < 1e-5, f"rel. diff {abs(rb - rf) / abs(rf):.1e}")

G = mp.gamma(mp.mpf(1) / 3)
tau = mp.mpc(-0.5, mp.sqrt(3) / 2)
qn = mp.exp(2j * mp.pi * tau)


def E_num(k):
    return (1 - 2 * k / mp.bernoulli(k) * mp.nsum(lambda n: n ** (k - 1) * qn ** n / (1 - qn ** n), [1, mp.inf])).real


E6r = E_num(6)
check("6 E_6(rho) = 27 Gamma(1/3)^18 / (512 pi^12) (Chowla–Selberg)", abs(E6r - 27 * G ** 18 / (512 * mp.pi ** 12)) < mp.mpf(10) ** -35)
d6 = abs(R_fast(6) - (1 - 9 * G ** 18 / (512 * mp.pi ** 12)))
d12 = abs(R_fast(12) - (30375 * G ** 36 / (90570752 * mp.pi ** 24) - 1))
check("7 R(6) = 1 - 9 Gamma(1/3)^18/(512 pi^12) and R(12) = 30375 Gamma(1/3)^36/(90570752 pi^24) - 1",
      d6 < mp.mpf(10) ** -30 and d12 < mp.mpf(10) ** -30, f"{mp.nstr(d6, 3)}, {mp.nstr(d12, 3)}")


# ---------- 8. c_k and R(6k) = (-1)^k (c_k E_6(rho)^k - 3)/3 ----------
def eis(k, Mq):
    Bk = sp.bernoulli(k)
    return [sp.Integer(1)] + [sp.Rational(-2 * k, 1) / Bk * sp.divisor_sigma(n, k - 1) for n in range(1, Mq)]


def mul(x, y, Mq):
    return [sum(x[i] * y[n - i] for i in range(n + 1)) for n in range(Mq)]


def c_exact(k):
    wt = 6 * k
    basis = [(i, (wt - 4 * i) // 6) for i in range(wt // 4 + 1) if (wt - 4 * i) % 6 == 0]
    Mq = len(basis)
    E4, E6 = eis(4, Mq), eis(6, Mq)
    cols = []
    for i, j in basis:
        v = [sp.Integer(1)] + [sp.Integer(0)] * (Mq - 1)
        for _ in range(i):
            v = mul(v, E4, Mq)
        for _ in range(j):
            v = mul(v, E6, Mq)
        cols.append(v)
    coef = sp.Matrix(cols).T.solve(sp.Matrix(eis(wt, Mq)))
    return sp.nsimplify(coef[[t for t, (i, _) in enumerate(basis) if i == 0][0]])


ck = [c_exact(k) for k in range(1, 13)]
check("8 c_1..c_5 = 1, 250/691, 5500/43867, 10285000/236364091, 26021050000/1723168255201 (as in the note)",
      ck[:5] == [sp.Integer(1), sp.Rational(250, 691), sp.Rational(5500, 43867), sp.Rational(10285000, 236364091),
                 sp.Rational(26021050000, 1723168255201)])
worst_c = max(abs(mp.mpf(c.p) / c.q - E_num(6 * k) / E6r ** k) / abs(E_num(6 * k) / E6r ** k) for k, c in zip(range(1, 13), ck))
check("9 c_k = E_6k(rho)/E_6(rho)^k for k = 1..12", worst_c < mp.mpf(10) ** -35, f"max rel. diff {mp.nstr(worst_c, 3)}")
worst_R = max(abs(mp.mpf((-1) ** k) * (mp.mpf(c.p) / c.q * E6r ** k - 3) / 3 - R_fast(6 * k)) for k, c in zip(range(1, 6), ck))
check("10 R(6k) = (-1)^k (c_k E_6(rho)^k - 3)/3 agrees with the q-series for k = 1..5", worst_R < mp.mpf(10) ** -30, f"max diff {mp.nstr(worst_R, 3)}")
bn = [abs(int(sp.fraction(sp.bernoulli(6 * k))[0])) for k in range(1, 13)]
check("11 the denominator of c_k divides the numerator of B_6k (k = 1..12)", all(bb % c.q == 0 for bb, c in zip(bn, ck)))

TABLE = ("0.039486299736351458", "0.0013602406593698314", "5.0757704420300926e-5", "1.88165499540916937e-6", "6.96920718649738464e-8")
tab_ok = True
for k, txt in enumerate(TABLE, 1):
    mant, _, ex = txt.partition("e")
    digits = len(mant.replace(".", "").lstrip("0"))
    val = R_fast(6 * k)
    tab_ok &= mp.nstr(val, digits + 3, min_fixed=-mp.inf, max_fixed=mp.inf).startswith(mant) if not ex else \
        mp.nstr(val / mp.mpf(10) ** int(ex), digits + 3).startswith(mant)
check("12 the digits of R(6k), k = 1..5, printed in Table 1 of the note are correct (truncated, not rounded)", tab_ok)

e36, e60 = E_num(36), E_num(60)
check("13 remark (consistency): E_36(rho) = 3.00000000774..., E_60(rho) = 3.00000000000001... (truncated)",
      mp.nstr(e36, 20).startswith("3.00000000774") and mp.nstr(e60, 20).startswith("3.00000000000001"),
      f"{mp.nstr(e36, 15)}, {mp.nstr(e60, 17)}")

print(f"\n{sum(RESULTS)}/{len(RESULTS)} checks passed")
sys.exit(0 if all(RESULTS) else 1)
