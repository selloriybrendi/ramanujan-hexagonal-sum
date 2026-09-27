"""#43 chuqurlashtirish (retsenziya 5 savoli + yangi xulosalar).

Belgilar: w = a + b*omega, omega = e^{i pi/3}; |w|^2 = a^2+ab+b^2; A–B (LN IV b.217) yig'indisi mu>=1, nu in Z, gcd=1,
z = e^{i pi/6}(mu + nu*omega). R(s) = SUM_{a,b>=1, gcd=1} cos(s psi)/|w|^s, tan psi = (b-a)/(sqrt3 (a+b)).
Ayniyat (9.2.5):  L(s) := SUM (-1)^{n+1} sigma_{s-1}(n) e^{-n pi sqrt3}
                 = -2 Gamma(s) zeta(s)/(2pi)^s * { cos(pi s/6) + 2 cos(pi(s+1)/6) cos(pi(s-1)/6) R(s) }.

A) LEMMA (orbita bo'linishi): har primitiv w uchun 6 ta assotsiatdan AYNAN 3 tasi mu>=1 yarim tekislikda,
   faqat birlik orbitasi uchun 2 ta (1 va omega^{-1}). Chegara nurlarida (mu=0) faqat +-omega yotadi.
B) RANDOM TEKSHIRUV: s ~ U(2.1, 20), 40 ta; panjara aylana |w|<=Rr bilan kesiladi + analitik dum tuzatmasi
   T(R) = 12/(sqrt3 pi^2) * 2 sin(s pi/6)/s * R^{2-s}/(s-2); xato bahosi = |natija(Rr) - natija(Rr/2)|.
C) KOMPLEKS s: ayniyat Re s > 2 da analitik — s = 5+3i, 7.3-2i, 3.5+1i.
D) YANGI XULOSA: R(6k) = (-1)^k (E_{6k}(rho) - 3)/3, E_{6k}(rho) = c_k E_6(rho)^k (E_4(rho)=0), c_k ratsional.
E) Sodda hol: E_{6k}(rho) -> 3 (6 birlik), demak R(6k) -> 0 — mos kelishi kerak.
"""
import math
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
rng = np.random.default_rng(20260927)

# ---------- A) LEMMA: to'liq sanab tekshirish ----------
def assotsiatlar(a, b):
    # w = a + b*omega; omega*w = -b + (a+b)*omega  (omega^2 = omega - 1)
    out, x, y = [], a, b
    for _ in range(6):
        out.append((x, y))
        x, y = -y, x + y
    return out

N = 120
yomon, jami, birlik = 0, 0, []
for a in range(-N, N + 1):
    for b in range(-N, N + 1):
        if (a, b) == (0, 0) or math.gcd(a, b) != 1:
            continue
        orb = assotsiatlar(a, b)
        k = sum(1 for (m, n) in orb if m >= 1)
        chegara = [(m, n) for (m, n) in orb if m == 0]
        if a * a + a * b + b * b == 1:
            birlik.append((k, sorted(chegara)))
        elif k != 3 or chegara:
            yomon += 1
        jami += 1
print(f"A) lemma: |a|,|b|<={N} dagi {jami} primitiv w; birlik bo'lmaganlarda 'aynan 3' buzilishi = {yomon}")
print(f"   birlik orbitasi: mu>=1 dagilar soni va chegara nurlari = {birlik[0]}")
# a,b>=1 vakil har orbitada AYNAN bitta ekanligi
vakil_yomon = 0
for a in range(-40, 41):
    for b in range(-40, 41):
        if (a, b) != (0, 0) and math.gcd(a, b) == 1 and a * a + a * b + b * b > 1:
            if sum(1 for (m, n) in assotsiatlar(a, b) if m >= 1 and n >= 1) != 1:
                vakil_yomon += 1
print(f"   har birlik bo'lmagan orbitada a,b>=1 vakil aynan 1 ta: buzilish = {vakil_yomon}")

# ---------- yordamchilar ----------
def L_q(s):
    s = mp.mpc(s)
    q0 = mp.e ** (-mp.pi * mp.sqrt(3))
    tot = mp.mpf(0)
    for n in range(1, 80):
        sig = mp.fsum(mp.mpf(d) ** (s - 1) for d in range(1, n + 1) if n % d == 0)
        tot += (-1) ** (n + 1) * sig * q0 ** n
    return tot

def olti(s):
    s = mp.mpc(s)
    pre = -2 * mp.gamma(s) * mp.zeta(s) / (2 * mp.pi) ** s
    return pre, mp.cos(mp.pi * s / 6), 2 * mp.cos(mp.pi * (s + 1) / 6) * mp.cos(mp.pi * (s - 1) / 6)

Rr = 2400
a_ = np.arange(1, Rr + 1)
A, B = np.meshgrid(a_, a_, indexing="ij")
Q = (A * A + A * B + B * B).astype(np.float64)
msk = (np.gcd(A, B) == 1) & (Q <= Rr * Rr)
A1, B1, Q1 = A[msk].astype(float), B[msk].astype(float), Q[msk]
PSI = np.arctan((B1 - A1) / (math.sqrt(3) * (A1 + B1)))
Rad = np.sqrt(Q1)

def R_panjara(s, Rmax, K=None, Lm=None):
    m = Rad <= Rmax
    psi = PSI[m] if K is None else np.arctan(K[m] / (math.sqrt(3) * Lm[m]))
    val = np.sum(np.cos(s * psi) * Q1[m] ** (-s / 2))
    tail = 12 / (math.sqrt(3) * math.pi ** 2) * 2 * math.sin(s * math.pi / 6) / s * Rmax ** (2 - s) / (s - 2)
    return val + tail

# ---------- B) random s ----------
print("\nB) random s in (2.1, 20), 40 ta: nisbiy farq(K=mu-nu, lam=mu+nu) va xato bahosi; nazorat (K=mu-nu, lam=mu)")
ss = np.sort(rng.uniform(2.1, 20, 40))
rows = []
for s in ss:
    Lq = L_q(s)
    pre, c1, c2 = olti(s)
    if abs(c2) < 1e-6:
        continue
    r1, r2 = R_panjara(s, Rr), R_panjara(s, Rr / 2)
    rhs1, rhs2 = pre * (c1 + c2 * r1), pre * (c1 + c2 * r2)
    fr = float(abs(rhs1 - Lq) / abs(Lq))
    est = float(abs(rhs1 - rhs2) / abs(Lq))
    rn = np.sum(np.cos(s * np.arctan((A1 - B1) / (math.sqrt(3) * A1))) * Q1 ** (-s / 2))
    fn = float(abs(pre * (c1 + c2 * rn) - Lq) / abs(Lq))
    rows.append((s, fr, est, fn))
ok = sum(1 for s, fr, est, fn in rows if fr <= 10 * est + 1e-13)
for s, fr, est, fn in rows[::4]:
    print(f"   s={s:7.4f}: farq={fr:.2e}  xato-bahosi={est:.2e}  nazorat-farq={fn:.2e}")
print(f"   jami {len(rows)} ta s: farq <= 10*xato-bahosi (+1e-13) bo'lganlar = {ok}/{len(rows)}; "
      f"nazorat eng kichik farqi = {min(r[3] for r in rows):.2e}")

# ---------- C) kompleks s ----------
print("\nC) kompleks s:")
for s in (mp.mpc(5, 3), mp.mpc(7.3, -2), mp.mpc(3.5, 1)):
    Lq = L_q(s)
    pre, c1, c2 = olti(s)
    sc = complex(s)
    r1 = np.sum(np.cos(sc * PSI) * Q1 ** (-sc / 2)) + 12 / (math.sqrt(3) * math.pi ** 2) * 2 * np.sin(sc * math.pi / 6) / sc * Rr ** (2 - sc) / (sc - 2)
    rhs = pre * (c1 + c2 * mp.mpc(r1))
    print(f"   s={complex(s)}: |RHS-L|/|L| = {float(abs(rhs - Lq) / abs(Lq)):.2e}")

# ---------- D) R(6k) yopiq formula ----------
print("\nD) R(6k) = (-1)^k (E_6k(rho) - 3)/3,  E_6k(rho) = c_k E_6(rho)^k:")
qs = sp.symbols('q')
def eis(k, M):
    Bk = sp.bernoulli(k)
    return [sp.Integer(1)] + [sp.Rational(-2 * k, 1) / Bk * sum(d ** (k - 1) for d in sp.divisors(n)) for n in range(1, M)]
def mul(x, y, M):
    return [sum(x[i] * y[n - i] for i in range(n + 1)) for n in range(M)]
tau = mp.mpc(-0.5, mp.sqrt(3) / 2)
qn = mp.exp(2j * mp.pi * tau)
def E_num(k):
    return 1 - 2 * k / mp.bernoulli(k) * mp.nsum(lambda n: n ** (k - 1) * qn ** n / (1 - qn ** n), [1, mp.inf])
E6r = E_num(6).real
G = mp.gamma(mp.mpf(1) / 3)
print(f"   E6(rho) - 27 Gamma(1/3)^18/(512 pi^12) = {mp.nstr(E6r - 27 * G ** 18 / (512 * mp.pi ** 12), 5)}")
for k in range(1, 6):
    w = 6 * k
    M = w // 12 + 3
    # M_w bazasi: E4^i E6^j, 4i+6j=w
    basis = [(i, (w - 4 * i) // 6) for i in range(w // 4 + 1) if (w - 4 * i) % 6 == 0]
    E4, E6 = eis(4, M + 4), eis(6, M + 4)
    vecs = []
    for i, j in basis:
        v = [sp.Integer(1)] + [0] * (M + 3)
        for _ in range(i):
            v = mul(v, E4, M + 4)
        for _ in range(j):
            v = mul(v, E6, M + 4)
        vecs.append(v[:len(basis)])
    target = eis(w, M + 4)[:len(basis)]
    coef = sp.Matrix(vecs).T.solve(sp.Matrix(target))
    ck_sym = sp.nsimplify(coef[[b for b in range(len(basis)) if basis[b][0] == 0][0]])
    ck = mp.mpf(ck_sym.p) / ck_sym.q
    Ew_num = E_num(w).real
    R_closed = mp.mpf((-1) ** k) * (ck * E6r ** k - 3) / 3
    s = mp.mpf(w)
    pre, c1, c2 = olti(s)
    R_from_L = ((L_q(s) / pre) - c1) / c2
    print(f"   k={k}: c_k={ck_sym}  |c_k E6^k - E_{w}(rho)| = {mp.nstr(abs(ck * E6r ** k - Ew_num), 3)}  "
          f"R({w}) = {mp.nstr(R_closed, 18)}  |yopiq - q-qator| = {mp.nstr(abs(R_closed - R_from_L), 3)}  "
          f"|yopiq - panjara| = {abs(float(R_closed) - R_panjara(w, Rr)):.1e}")

# ---------- E) sodda hol ----------
print("\nE) E_6k(rho) -> 3 (6 ta birlik, har biri u^{-6k} = 1):")
for k in (1, 3, 6, 10):
    print(f"   k={k}: E_{6*k}(rho) = {mp.nstr(E_num(6 * k).real, 15)}")
