"""Symbolic check of the three identities used in the proof (all differences must simplify to 0)."""
import sympy as sp

a, b, s, psi = sp.symbols("a b s psi", positive=True)
w = sp.expand(a + b * sp.exp(sp.I * sp.pi / 3))
theta = sp.atan2(sp.im(w), sp.re(w))
print("1) tan(arg w - pi/6) - (b-a)/(sqrt3 (a+b)) =",
      sp.simplify(sp.expand_trig(sp.tan(theta - sp.pi / 6)) - (b - a) / (sp.sqrt(3) * (a + b))))
L = sp.cos(s * (psi + sp.pi / 3)) + sp.cos(s * psi) + sp.cos(s * (psi - sp.pi / 3))
print("2) three rotations - cos(s psi)(1+2cos(pi s/3)) =",
      sp.simplify(sp.expand_trig(sp.expand(L - sp.cos(s * psi) * (1 + 2 * sp.cos(sp.pi * s / 3))))))
print("3) 1+2cos(pi s/3) - 4cos(pi(s+1)/6)cos(pi(s-1)/6) =",
      sp.simplify(sp.expand_trig(1 + 2 * sp.cos(sp.pi * s / 3)
                                 - 4 * sp.cos(sp.pi * (s + 1) / 6) * sp.cos(sp.pi * (s - 1) / 6))))
print("4) |w|^2 - (a^2+ab+b^2) =", sp.nsimplify(sp.expand_complex(sp.expand(w * sp.conjugate(w)) - (a * a + a * b + b * b))))
