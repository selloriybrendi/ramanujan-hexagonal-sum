#!/usr/bin/env python3
"""Mutation test for verify_all.py: each mutant alters one stated formula or constant; the named check must then fail.
Run: python3 mutation_test.py   (exit code 0 = every mutant is detected)"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, "verify_all.py")).read()
MUTANTS = [  # (description, original text, altered text, check that must fail)
    ("tan psi sign", "(b - a) / (sp.sqrt(3) * (a + b))", "(a - b) / (sp.sqrt(3) * (a + b))", "1"),
    ("product formula", "4 * sp.cos(sp.pi * (s + 1) / 6)", "4 * sp.cos(sp.pi * (s + 2) / 6)", "1"),
    ("orbit rotation", "x, y = -y, x + y", "x, y = -y, x - y", "2"),
    ("K = mu+nu", "right(sv, M - V, M + V)", "right(sv, M + V, M + V)", "3"),
    ("lambda = mu", "right(sv, M - V, M + V)", "right(sv, M - V, M)", "3"),
    ("unit term cos -> sin", "(c(math.pi * sv / 6) + 2", "(math.sin(math.pi * sv / 6) + 2", "3"),
    ("control replaced by the main choice", '("mu-nu,mu", M - V, M)', '("mu-nu,mu", M - V, M + V)', "4"),
    ("reflection -> shift", "+ math.atan((m - v) / ((m + v)", "- math.atan((m - v) / ((m + v)", "4b"),
    ("E6(rho): 27 -> 28", "27 * G ** 18 / (512", "28 * G ** 18 / (512", "6"),
    ("R(6): 9 -> 8", "1 - 9 * G ** 18", "1 - 8 * G ** 18", "7"),
    ("R(12): 30375 -> 30376", "30375 * G ** 36", "30376 * G ** 36", "7"),
    ("c_2: 250 -> 251", "sp.Rational(250, 691)", "sp.Rational(251, 691)", "8"),
    ("R(6k) sign", "mp.mpf((-1) ** k) * (", "mp.mpf((-1) ** (k + 1)) * (", "10"),
    ("R(6k): -3 -> -2", "E6r ** k - 3) / 3", "E6r ** k - 2) / 3", "10"),
    ("B_6k -> B_6k+2", "sp.bernoulli(6 * k))[0]", "sp.bernoulli(6 * k + 2))[0]", "11"),
    ("table digit", '"1.88165499540916937e-6"', '"1.88165499540916938e-6"', "12"),
    ("E_60 digit", 'startswith("3.00000000000001")', 'startswith("3.00000000000002")', "13"),
]
killed = 0
with tempfile.TemporaryDirectory() as tmp:
    for i, (name, old, new, must_fail) in enumerate(MUTANTS):
        assert old in SRC, name
        path = os.path.join(tmp, f"mutant_{i}.py")
        open(path, "w").write(SRC.replace(old, new, 1))
        r = subprocess.run([sys.executable, path], capture_output=True, text=True, timeout=600)
        failed = [ln.split("] ")[1].split(" ")[0] for ln in r.stdout.splitlines() if ln.startswith("[FAIL]")]
        ok = r.returncode != 0 and must_fail in failed
        killed += ok
        print(f"[{'KILLED' if ok else 'SURVIVED'}] {name:36s} failed checks: {failed}")
print(f"\n{killed}/{len(MUTANTS)} mutants detected")
sys.exit(0 if killed == len(MUTANTS) else 1)
