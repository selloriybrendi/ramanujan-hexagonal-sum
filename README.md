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
sides agree to relative 8e-16 … 3.5e-12 (3.8 million coprime pairs; left side to 30 digits); other natural
choices of K, λ miss by 5e-6 … 1.5e-1 (at s = 6 the choice K = μ+2ν, λ = μ also agrees: its angles are π/3 − ψ,
and cos(6(π/3 − ψ)) = cos(6ψ)). At s ≡ 2, 4 (mod 6) the coefficient of R(s) vanishes.

**What it gives.** R(s) is a slowly convergent lattice sum (error ≈ N^{2−s} after N² terms); the identity
evaluates it from about twenty terms of a q-series with q = e^{−π√3} ≈ 0.0043. With E₄(ρ) = 0 and
E₆(ρ) = 27Γ(1/3)¹⁸/(512π¹²) (Chowla–Selberg) it yields closed forms, e.g.

    R(6)  = 1 − 9Γ(1/3)¹⁸ / (512π¹²)
    R(12) = 30375Γ(1/3)³⁶ / (90570752π²⁴) − 1

(`verify/fast_R_and_closed_forms.py`, agreement to 40 digits).

## Version 1.1 (2026-09-27)

- `paper/Kenjaev_Ramanujan_9_2_5.pdf` — 4-page paper: Lemma 1 (orbit partition, full proof), Lemma 2 (sign),
  Theorem 1 (K, λ), Theorem 2 (fast evaluation), Corollary 1 (closed forms), numerical evidence.
- **Closed forms for all k:** R(6k) = (−1)^k (c_k E₆(ρ)^k − 3)/3, E₆(ρ) = 27Γ(1/3)¹⁸/(512π¹²), where c_k is the coefficient of
  E₆^k in E_{6k} written in E₄, E₆: c₁ = 1, c₂ = 250/691, c₃ = 5500/43867, c₄ = 10285000/236364091,
  c₅ = 26021050000/1723168255201 (checked to 40 digits for k ≤ 5).
- `verify/deepening_checks.py` (comments in Uzbek): exhaustive check of Lemma 1 (35088 primitive w, 0 violations),
  40 random s in (2.1, 20) with tail-corrected lattice sums (40/40 within the error estimate), complex s, R(6k).
- Wording: numerical agreement is empirical, not a rigorous bound (relative 8e-16 … 3.5e-12 at the fixed s; 5.4e-17 … 5.8e-11
  at the three complex s).

## Version 1.3 (2026-09-28)

- `verify/verify_all.py` — one command, 14 pass/fail checks, exit code 0 iff all pass (≈ 2 s): the symbolic identities,
  Lemma 1, Theorem 1 at seven s with controls, Theorem 2 against the lattice sum, Chowla–Selberg, R(6), R(12), c₁…c₁₂,
  R(6k) for k ≤ 5, every digit of Table 1, and the two values in Remark 2.
- `verify/mutation_test.py` — 17 deliberate alterations of stated formulas/constants; each one makes its check fail (17/17).
- `make_figure.py` → `fig/hex_925.pdf` (Figure 1 of the paper): the orbit of 2+ω and the curve 3^{s/2}R(s) with the points R(6k).
- `verify/c_k_sequence.py` — c₁…c₁₂ exactly; `verify/make_oeis_ck.py` → `oeis/` drafts for the numerators and denominators of c_k
  (neither is in the OEIS as of September 2026).
- Paper (6 pages): Figure 1, Remark 3 on c_k, reproducibility paragraph.
- Corrections found by the new checks: Table 1 now shows truncated digits (k = 1, 3, 5 had been rounded in the last place), and
  the explanation of the coincidence at s = 6 is a reflection ψ ↦ π/3 − ψ, not a shift by π/3.
- v1.3.1: Zenodo description (`.zenodo.json`) brought in line with the paper (it still quoted the v1.1 wording
  "1e-12 to 1e-20" and "4-page paper").

## Cite as

Kenjaev, O. U. (2026). *The constants K and λ in Ramanujan's identity (9.2.5) of the Lost Notebook, with a hexagonal lattice sum calculator* (v1.3.1). Zenodo. https://doi.org/10.5281/zenodo.22987480

(10.5281/zenodo.22987480 always resolves to the latest version.)

Author: Otakhon U. Kenjaev · ORCID https://orcid.org/0009-0009-3566-9285 · Telegram https://t.me/sheki · channel https://t.me/nizomliy

---
**Author:** Otakhon U. Kenjaev (also written *Otaxon Kenjayev* / *Отахон Кенжаев*) · [otakhonkenjaev.com](https://otakhonkenjaev.com/) · ORCID [0009-0009-3566-9285](https://orcid.org/0009-0009-3566-9285)
