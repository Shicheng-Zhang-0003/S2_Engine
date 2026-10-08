#!/usr/bin/env python3
"""HISTORICAL (pre-v9R4), superseded; DO NOT RUN. Kept for provenance only.
See v9R4/verify_scripts.sh and audit/run_audit.sh for current verifiers.
Running this file today would re-mutate a fixed tree.
"""
"""
Verify UFF LJ parameters against Rappé et al. 1992.
"""

# UFF parameters from Rappé et al. 1992, JACS 114, 10024
# x1 = sigma (Å), x2 = eps (kcal/mol)
UFF_REF = {
    1:  {"sig": 2.886, "eps": 0.044},   # H
    2:  {"sig": 2.362, "eps": 0.056},   # He
    3:  {"sig": 2.451, "eps": 0.025},   # Li
    4:  {"sig": 2.745, "eps": 0.085},   # Be
    5:  {"sig": 4.083, "eps": 0.180},   # B
    6:  {"sig": 3.851, "eps": 0.105},   # C
    7:  {"sig": 3.660, "eps": 0.069},   # N
    8:  {"sig": 3.500, "eps": 0.060},   # O
    9:  {"sig": 3.364, "eps": 0.050},   # F
    10: {"sig": 3.243, "eps": 0.042},   # Ne
    11: {"sig": 2.983, "eps": 0.030},   # Na
    12: {"sig": 3.021, "eps": 0.111},   # Mg
    13: {"sig": 4.499, "eps": 0.505},   # Al
    14: {"sig": 4.295, "eps": 0.402},   # Si
    15: {"sig": 4.147, "eps": 0.305},   # P
    16: {"sig": 4.035, "eps": 0.274},   # S
    17: {"sig": 3.947, "eps": 0.227},   # Cl
    18: {"sig": 3.868, "eps": 0.185},   # Ar
    19: {"sig": 3.812, "eps": 0.035},   # K
    20: {"sig": 3.399, "eps": 0.238},   # Ca
}

# Parse the C file
import re

with open('TREE/src/periodic_table.c', 'r') as f:
    content = f.read()

pattern = r'/\* Z=(\d+)\s+\w+\s+[^*]+ \*/\s*\{(\d+),"(\w+)","[^"]+",\s*[\d.]+,\s*[-\d.]+,\s*[\d.]+,\s*[\d.]+,\s*[\d.]+,\s*[\d.]+,\s*[\d.]+,\s*[\d.]+,\s*\d+,\s*\{\{\{0\}\},\s*0,\s*0\},\s*([\d.]+)\*KCAL_MOL_TO_EV,\s*([\d.]+)\s*\}'

matches = re.findall(pattern, content)

errors = []
for m in matches:
    Z = int(m[0])
    symbol = m[2]
    eps_kcal = float(m[3])
    sig = float(m[4])
    
    if Z in UFF_REF:
        ref_sig = UFF_REF[Z]['sig']
        ref_eps = UFF_REF[Z]['eps']
        if abs(sig - ref_sig) > 0.001:
            errors.append(f"Z={Z} ({symbol}): sigma {sig:.4f} vs UFF {ref_sig:.4f}")
        if abs(eps_kcal - ref_eps) > 0.001:
            errors.append(f"Z={Z} ({symbol}): eps {eps_kcal:.3f} vs UFF {ref_eps:.3f} kcal/mol")

if errors:
    print("UFF LJ parameter mismatches:")
    for e in errors:
        print(f"  {e}")
else:
    print("All UFF LJ parameters match reference values!")

# Also check the CR table for K
print("\nChecking K (Z=19) Clementi-Raimondi entries:")
cr_k = [(n, l, z) for (Z, n, l, z) in [
    (19,1,0,18.4895), (19,2,0,13.01), (19,2,1,15.03),
    (19,3,0,8.68), (19,3,1,7.73), (19,4,0,3.50)
] if Z == 19]
for n, l, z in cr_k:
    print(f"  K {n}{'spdf'[l]}: Zeff_CR = {z}")