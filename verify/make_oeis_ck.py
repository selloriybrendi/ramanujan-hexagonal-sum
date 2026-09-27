#!/usr/bin/env python3
"""Writes OEIS drafts for the numerators and denominators of c(n) = E_{6n}(rho)/E_6(rho)^n (oeis/*.txt).
Terms are computed exactly (c_k_sequence.c_exact) and the %S line is cut at the OEIS limit of about 260 characters."""
import os

import sympy as sp

from c_k_sequence import c_exact

LIMIT = 260
os.makedirs("../oeis", exist_ok=True)
ck = [c_exact(k) for k in range(1, 21)]

# assertions behind the comments of the drafts
E4, E6 = sp.symbols("E4 E6")
assert ck[1] == sp.Rational(250, 691)  # E_12 = (441 E_4^3 + 250 E_6^2)/691
for k, c in enumerate(ck, 1):
    assert abs(sp.fraction(sp.bernoulli(6 * k))[0]) % c.q == 0, k


def terms(seq):
    out, n = [], 0
    for t in seq:
        s = ",".join(out + [str(t)])
        if len(s) > LIMIT:
            break
        out.append(str(t))
        n += 1
    return ",".join(out), n


PROG = """%o A000000 (Python)
%o A000000 from sympy import bernoulli, divisor_sigma, Matrix, Rational
%o A000000 def eis(k, M): return [1] + [Rational(-2*k)/bernoulli(k)*divisor_sigma(n, k-1) for n in range(1, M)]
%o A000000 def mul(x, y): return [sum(x[i]*y[n-i] for i in range(n+1)) for n in range(len(x))]
%o A000000 def c(k):  # coefficient of E6^k in E_{6k} as a polynomial in E4, E6
%o A000000     B = [(i, (6*k-4*i)//6) for i in range(6*k//4+1) if (6*k-4*i) % 6 == 0]; M = len(B)
%o A000000     cols = []
%o A000000     for i, j in B:
%o A000000         v = [1] + [0]*(M-1)
%o A000000         for _ in range(i): v = mul(v, eis(4, M))
%o A000000         for _ in range(j): v = mul(v, eis(6, M))
%o A000000         cols.append(v)
%o A000000     return Matrix(cols).T.solve(Matrix(eis(6*k, M)))[[t for t, b in enumerate(B) if b[0] == 0][0]]
%o A000000 print([c(k).{part} for k in range(1, 13)])"""

COMMON = """%C A000000 c(n) is the coefficient of E_6^n when E_{{6n}} is written as a polynomial in E_4 and E_6; every other monomial contains E_4, and E_4(rho) = 0.
%C A000000 E_6(rho) = 27*Gamma(1/3)^18/(512*Pi^12) (Chowla-Selberg). As n -> infinity, E_{{6n}}(rho) -> 3, since the six units of Z[rho] contribute 1 each (halved).
%C A000000 c(n) gives closed forms for the sector lattice sum R(s) = Sum_{{mu,nu >= 1, gcd(mu,nu)=1}} cos(s*arctan((mu-nu)/((mu+nu)*sqrt(3))))/(mu^2+mu*nu+nu^2)^(s/2), which supplies the constants K = mu-nu and lambda = mu+nu in Ramanujan's identity (9.2.5) of the Lost Notebook: R(6n) = (-1)^n*(c(n)*E_6(rho)^n - 3)/3.
%C A000000 The denominator of c(n) divides the numerator of B_{{6n}} (verified for n <= 20); the quotient is 1 for n <= 4 and first differs from 1 at n = 5 (quotient 5).
%D A000000 G. E. Andrews and B. C. Berndt, Ramanujan's Lost Notebook, Part IV, Springer, 2013, Sections 9.2-9.3.
%H A000000 O. U. Kenjaev, <a href="https://doi.org/10.5281/zenodo.22987480">The constants K and lambda in Ramanujan's identity (9.2.5) of the Lost Notebook</a>, Zenodo, 2026.
%e A000000 c(2) = 250/691 because E_12 = (441*E_4^3 + 250*E_6^2)/691.
{prog}
%Y A000000 Cf. A004009 (E_4), A013973 (E_6), A029828 (691*E_12), A000367 (numerators of B_2n), {other}.
%K A000000 nonn,frac
%O A000000 1,2
%A A000000 _Otakhon U. Kenjaev_, Sep 28 2026
"""

quot = [abs(sp.fraction(sp.bernoulli(6 * k))[0]) // c.q for k, c in enumerate(ck, 1)]
assert quot[:4] == [1, 1, 1, 1] and quot[4] == 5, quot

for part, name, other in (("p", "Numerators", "denominators: A000000"), ("q", "Denominators", "numerators: A000000")):
    S, n = terms(int(getattr(c, part)) for c in ck)
    txt = (f"%I A000000\n%S A000000 {S}\n"
           f"%N A000000 {name} of c(n) = E_{{6n}}(rho)/E_6(rho)^n, where E_k is the normalized Eisenstein series and rho = exp(2*Pi*i/3).\n"
           + COMMON.format(prog=PROG.replace("{part}", part), other=other))
    fn = f"../oeis/c_k-{name.lower()}.txt"
    open(fn, "w").write(txt)
    print(f"{fn}: {n} terms, %S length {len(S)}")
