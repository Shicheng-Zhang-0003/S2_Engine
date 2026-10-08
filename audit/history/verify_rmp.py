#!/usr/bin/env python3
"""HISTORICAL (pre-v9R4), superseded; DO NOT RUN. Kept for provenance only.
See v9R4/verify_scripts.sh and audit/run_audit.sh for current verifiers.
Running this file today would re-mutate a fixed tree.
"""
"""
Verify most probable radius calculation against analytical values.
"""

import numpy as np
from scipy.optimize import minimize_scalar

BOHR_TO_ANGSTROM = 0.529177210903

def factorial(n):
    if n < 0:
        return 0
    result = 1
    for i in range(2, n+1):
        result *= i
    return result

def laguerre(p, q, x):
    if p == 0:
        return 1.0
    if p == 1:
        return 1.0 + q - x
    L_prev2 = 1.0
    L_prev1 = 1.0 + q - x
    for k in range(2, p+1):
        L_curr = ((2*k - 1 + q - x) * L_prev1 - (k - 1 + q) * L_prev2) / k
        L_prev2, L_prev1 = L_prev1, L_curr
    return L_prev1

def radial_probability(n, l, Z_eff, r_ang):
    if r_ang < 0 or n < 1 or l < 0 or l >= n or Z_eff <= 0:
        return 0.0
    a0 = BOHR_TO_ANGSTROM
    scale = 2.0 * Z_eff / (n * a0)
    rho = scale * r_ang
    num = factorial(n - l - 1)
    den = 2.0 * n * factorial(n + l)
    N2 = scale**3 * num / den
    N = -np.sqrt(N2)
    lag = laguerre(n - l - 1, 2*l + 1, rho)
    R = N * np.exp(-rho / 2.0) * (rho**l) * lag
    return r_ang**2 * R**2

def analytic_r_mp(n, l, Z_eff):
    """Analytical most probable radius for hydrogen-like orbitals"""
    a0 = BOHR_TO_ANGSTROM
    if l == n - 1:
        # For l = n-1 (circular orbits), r_mp = n^2 * a0 / Z_eff
        return n * n * a0 / Z_eff
    else:
        # For other cases, need to solve d/dr (r^2 R^2) = 0
        # This is more complex; use numerical for now
        return None

def numeric_r_mp(n, l, Z_eff):
    """Find most probable radius numerically using golden-section search"""
    a = 0.0001
    b = 30.0 * n * n / Z_eff
    
    # Partition into subintervals to find global max
    best_r = a
    best_p = -1.0
    NSUB = 64
    
    for s in range(NSUB):
        lo = a + (b - a) * s / NSUB
        hi = a + (b - a) * (s + 1) / NSUB
        phi = 0.6180339887
        c = hi - phi * (hi - lo)
        d = lo + phi * (hi - lo)
        for _ in range(60):
            if abs(hi - lo) < 1e-10:
                break
            if radial_probability(n, l, Z_eff, c) < radial_probability(n, l, Z_eff, d):
                lo = c
            else:
                hi = d
            c = hi - phi * (hi - lo)
            d = lo + phi * (hi - lo)
        r = (lo + hi) / 2.0
        p = radial_probability(n, l, Z_eff, r)
        if p > best_p:
            best_p = p
            best_r = r
    return best_r

# Analytical values for comparison
# H 1s: r_mp = a0 = 0.529177 Å
# H 2s: r_mp = 5.24 a0 = 2.77 Å (has two peaks, global at ~5.24 a0)
# H 2p: r_mp = 4 a0 = 2.12 Å
# H 3s: r_mp ≈ 13.07 a0 = 6.92 Å
# H 3p: r_mp = 12 a0 = 6.35 Å
# H 3d: r_mp = 9 a0 = 4.76 Å

test_cases = [
    (1, 0, 1.0, 0.529177),   # H 1s
    (2, 0, 1.0, 2.77),       # H 2s (global max)
    (2, 1, 1.0, 2.1167),     # H 2p = 4*a0
    (3, 0, 1.0, 6.92),       # H 3s (global max ~13.07*a0)
    (3, 1, 1.0, 6.35),       # H 3p = 12*a0
    (3, 2, 1.0, 4.76),       # H 3d = 9*a0
]

print("Most probable radius verification:")
print("=" * 70)
print(f"  {'n':>2} {'l':>2} {'Z_eff':>6} {'Numeric (Å)':>12} {'Analytic (Å)':>12} {'Diff (Å)':>10} {'Status'}")
print("-" * 70)

for n, l, Z_eff, expected in test_cases:
    numeric = numeric_r_mp(n, l, Z_eff)
    diff = abs(numeric - expected)
    status = "✓" if diff < 0.01 else "✗" if diff < 0.1 else "✗✗"
    print(f"  {n:>2} {l:>2} {Z_eff:>6.2f} {numeric:>12.4f} {expected:>12.4f} {diff:>10.4f} {status}")

print("\nNote: For n=2, l=0 (2s), there are two peaks (at ~0.76 and ~5.24 a0).")
print("The global maximum is at ~5.24 a0 = 2.77 Å.")