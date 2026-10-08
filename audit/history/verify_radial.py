#!/usr/bin/env python3
"""HISTORICAL (pre-v9R4), superseded; DO NOT RUN. Kept for provenance only.
See v9R4/verify_scripts.sh and audit/run_audit.sh for current verifiers.
Running this file today would re-mutate a fixed tree.
"""
"""
Verify quantum radial wavefunction normalization.
"""

import numpy as np
from scipy.integrate import quad

# Constants
BOHR_TO_ANGSTROM = 0.529177210903

def factorial(n):
    if n < 0:
        return 0
    result = 1
    for i in range(2, n+1):
        result *= i
    return result

def laguerre(p, q, x):
    """Associated Laguerre polynomial L_p^q(x)"""
    if p == 0:
        return 1.0
    if p == 1:
        return 1.0 + q - x
    
    L_prev2 = 1.0
    L_prev1 = 1.0 + q - x
    L_curr = 0.0
    
    for k in range(2, p+1):
        L_curr = ((2*k - 1 + q - x) * L_prev1 - (k - 1 + q) * L_prev2) / k
        L_prev2 = L_prev1
        L_prev1 = L_curr
    return L_curr

def radial_wavefunction(n, l, Z_eff, r_ang):
    """Hydrogen-like radial wavefunction R_nl(r)"""
    if r_ang < 0.0 or n < 1 or l < 0 or l >= n or Z_eff <= 0.0:
        return 0.0
    
    a0_ang = BOHR_TO_ANGSTROM
    scale = 2.0 * Z_eff / (n * a0_ang)
    rho = scale * r_ang
    
    # Normalization: N^2 = scale^3 * (n-l-1)! / (2n (n+l)!)
    num = factorial(n - l - 1)
    den = 2.0 * n * factorial(n + l)
    N2 = scale**3 * num / den
    N = -np.sqrt(N2)  # Griffiths sign convention
    
    lag = laguerre(n - l - 1, 2*l + 1, rho)
    radial = N * np.exp(-rho / 2.0) * (rho**l) * lag
    
    return radial

def radial_probability(n, l, Z_eff, r_ang):
    R = radial_wavefunction(n, l, Z_eff, r_ang)
    return r_ang**2 * R**2

def test_normalization():
    print("Testing radial wavefunction normalization...")
    print("=" * 60)
    
    test_cases = [
        (1, 0, 1.0),   # H 1s
        (2, 0, 1.0),   # H 2s
        (2, 1, 1.0),   # H 2p
        (3, 0, 1.0),   # H 3s
        (3, 1, 1.0),   # H 3p
        (3, 2, 1.0),   # H 3d
        (1, 0, 2.0),   # He+ 1s
        (2, 0, 2.0),   # He+ 2s
        (2, 1, 2.0),   # He+ 2p
    ]
    
    for n, l, Z_eff in test_cases:
        # Integrate r^2 |R|^2 dr from 0 to infinity
        # In atomic units: ∫_0^∞ r^2 |R|^2 dr = 1
        # In Angstrom: same integral, just different units
        integrand = lambda r: radial_probability(n, l, Z_eff, r)
        result, error = quad(integrand, 0, 100*n**2/Z_eff, limit=1000, epsabs=1e-12, epsrel=1e-12)
        
        # Also test in atomic units (Bohr radii)
        a0 = BOHR_TO_ANGSTROM
        integrand_au = lambda r_au: (r_au*a0)**2 * radial_wavefunction(n, l, Z_eff, r_au*a0)**2 * a0
        result_au, error_au = quad(integrand_au, 0, 100*n**2/Z_eff, limit=1000, epsabs=1e-12, epsrel=1e-12)
        
        status = "✓" if abs(result - 1.0) < 1e-10 else "✗"
        print(f"  n={n}, l={l}, Z_eff={Z_eff}: ∫P(r)dr = {result:.10f} (err={error:.2e}) {status}")
        if abs(result - 1.0) > 1e-6:
            print(f"    WARNING: Normalization error = {abs(result - 1.0):.2e}")

if __name__ == "__main__":
    test_normalization()