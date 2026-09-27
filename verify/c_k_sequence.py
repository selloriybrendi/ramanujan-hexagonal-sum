#!/usr/bin/env python3
"""c_k = E_{6k}(rho) / E_6(rho)^k, the coefficient of E_6^k when E_{6k} is written as a polynomial in E_4, E_6
(rho = e^{2 pi i/3}, E_4(rho) = 0).  Computes c_1..c_K exactly (SymPy) and checks each against the numerical value
E_{6k}(rho)/E_6(rho)^k (mpmath, 50 digits).  Usage: python3 c_k_sequence.py [K]   (default K = 12)"""
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 50


def eis(k, M):  # q-expansion of E_k to O(q^M), normalized E_k = 1 - (2k/B_k) sum sigma_{k-1}(n) q^n
    Bk = sp.bernoulli(k)
    return [sp.Integer(1)] + [sp.Rational(-2 * k, 1) / Bk * sp.divisor_sigma(n, k - 1) for n in range(1, M)]


def mul(x, y, M):
    return [sum(x[i] * y[n - i] for i in range(n + 1)) for n in range(M)]


def c_exact(k):
    w = 6 * k
    basis = [(i, (w - 4 * i) // 6) for i in range(w // 4 + 1) if (w - 4 * i) % 6 == 0]
    M = len(basis)
    E4, E6 = eis(4, M), eis(6, M)
    cols = []
    for i, j in basis:
        v = [sp.Integer(1)] + [sp.Integer(0)] * (M - 1)
        for _ in range(i):
            v = mul(v, E4, M)
        for _ in range(j):
            v = mul(v, E6, M)
        cols.append(v)
    coef = sp.Matrix(cols).T.solve(sp.Matrix(eis(w, M)))
    return sp.nsimplify(coef[[t for t, (i, _) in enumerate(basis) if i == 0][0]])


tau = mp.mpc(-0.5, mp.sqrt(3) / 2)
qn = mp.exp(2j * mp.pi * tau)


def E_num(k):
    return (1 - 2 * k / mp.bernoulli(k) * mp.nsum(lambda n: n ** (k - 1) * qn ** n / (1 - qn ** n), [1, mp.inf])).real


if __name__ == "__main__":
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    E6r = E_num(6)
    ck = [c_exact(k) for k in range(1, K + 1)]
    worst = max(abs(mp.mpf(c.p) / c.q - E_num(6 * k) / E6r ** k) / abs(E_num(6 * k) / E6r ** k) for k, c in zip(range(1, K + 1), ck))
    nums = [int(c.p) for c in ck]
    dens = [int(c.q) for c in ck]
    print("c_k numerators  :", ",".join(map(str, nums)))
    print("c_k denominators:", ",".join(map(str, dens)))
    print(f"max relative difference c_k vs E_6k(rho)/E_6(rho)^k (k <= {K}): {mp.nstr(worst, 3)}")
    # relation of the denominators with numerators of Bernoulli numbers B_{6k}
    bn = [abs(int(sp.fraction(sp.bernoulli(6 * k))[0])) for k in range(1, K + 1)]
    print("denominator divides |numerator(B_6k)|:", all(b % dd == 0 for b, dd in zip(bn, dens)),
          "| quotients:", [b // dd for b, dd in zip(bn, dens)])
