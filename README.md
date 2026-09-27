# The constants K and λ in Ramanujan's identity (9.2.5)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22987480.svg)](https://doi.org/10.5281/zenodo.22987480)

Live calculator: https://selloriybrendi.github.io/ramanujan-hexagonal-sum/

On pages 270–271 of the Lost Notebook, Ramanujan wrote an identity for the divisor sums
σ<sub>s−1</sub>(n)e<sup>−nπ√3</sup> whose right-hand side contains two undefined constants K and λ.
Andrews and Berndt (*Ramanujan's Lost Notebook, Part IV*, Springer 2013, §9.2–9.3, pp. 215–216) write:
"we are unable either to identify them or to obtain an identity of the form given by Ramanujan".

**Result.** For every s > 2,

    σ_{s−1}(1)e^{−π√3} − σ_{s−1}(2)e^{−2π√3} + σ_{s−1}(3)e^{−3π√3} − ···
      = −2 Γ(s)ζ(s)/(2π)^s · { cos(πs/6) + 2cos(π(s+1)/6)cos(π(s−1)/6) · R(s) },

    R(s) = Σ_{μ,ν≥1, gcd(μ,ν)=1} cos( s·tan⁻¹( (μ−ν) / ((μ+ν)√3) ) ) / (μ²+μν+ν²)^{s/2},

that is, **K = μ − ν and λ = μ + ν** — the hexagonal counterpart of tan⁻¹((μ−ν)/(μ+ν)) in the
companion identity (9.2.4).

**Proof sketch.** Andrews and Berndt prove (p. 217) that the left side equals
Γ(s)ζ(s)/(2π)^s · Σ_{μ≥1, ν∈ℤ, gcd=1} cos(s·arg z)/|z|^s with z = e^{iπ/6}(μ+νω), ω = e^{iπ/3}.
The half-plane μ ≥ 1 is the sector −2π/3 < arg(μ+νω) < π/3. The six units ±1, ±ω, ±ω² permute primitive
points; for each primitive w = a+bω with a, b ≥ 1 exactly w, ω⁻¹w, ω⁻²w lie in the half-plane. With
ψ = arg w − π/6 their terms add to cos(s(ψ+π/3)) + cos(sψ) + cos(s(ψ−π/3)) = 4cos(π(s+1)/6)cos(π(s−1)/6)cos(sψ),
where |w|² = a²+ab+b² and tan ψ = (b−a)/(√3(a+b)). The orbit of w = 1 gives 2cos(πs/6). ∎
(`verify/sympy_check.py` checks the three identities symbolically.)

**Numerical check** (`verify/numeric_check_K_lambda.py`): at s = 4.7, 5.5, 6, 7.5, 8.7, 9.25, 13.1 the two
sides agree to relative 1e-12 … 1e-16 (3.8 million coprime pairs; left side to 30 digits); other natural
choices of K, λ miss by 5e-6 … 2e-1. At s ≡ 2, 4 (mod 6) the coefficient of R(s) vanishes.

**What it gives.** R(s) is a slowly convergent lattice sum (error ≈ N^{2−s} after N² terms); the identity
evaluates it from about twenty terms of a q-series with q = e^{−π√3} ≈ 0.0043. With E₄(ρ) = 0 and
E₆(ρ) = 27Γ(1/3)¹⁸/(512π¹²) (Chowla–Selberg) it yields closed forms, e.g.

    R(6)  = 1 − 9Γ(1/3)¹⁸ / (512π¹²)
    R(12) = 30375Γ(1/3)³⁶ / (90570752π²⁴) − 1

(`verify/fast_R_and_closed_forms.py`, agreement to 40 digits).

## Cite as

Kenjaev, O. U. (2026). *The constants K and λ in Ramanujan's identity (9.2.5) of the Lost Notebook, with a hexagonal lattice sum calculator* (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.22987480

(10.5281/zenodo.22987480 always resolves to the latest version; v1.0.1 itself is 10.5281/zenodo.22987481.)

Author: Otakhon U. Kenjaev · ORCID https://orcid.org/0009-0009-3566-9285 · Telegram https://t.me/sheki · channel https://t.me/nizomliy
