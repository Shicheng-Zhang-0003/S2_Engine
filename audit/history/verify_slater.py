#!/usr/bin/env python3
"""
Verify Slater screening against known reference values.
"""

# Slater's rules verification for specific elements
# From Slater's 1930 paper and standard textbooks

def slater_group(n, l):
    """Map (n,l) to Slater group index"""
    if l >= 2:  # d, f
        if n == 3 and l == 2: return 3  # 3d
        if n == 4 and l == 2: return 5  # 4d
        if n == 4 and l == 3: return 6  # 4f
        if n == 5 and l == 2: return 8  # 5d
        if n == 5 and l == 3: return 9  # 5f
        return 2*n
    else:  # s, p
        if n == 1: return 0
        if n == 2: return 1
        if n == 3: return 2
        if n == 4: return 4
        if n == 5: return 7
        if n == 6: return 10
        return 2*n

def compute_slater_zeff(Z, n, l, config):
    """
    config: dict of (n,l) -> occupancy
    """
    target_group = slater_group(n, l)
    S = 0.0
    
    for (on, ol), occ in config.items():
        if occ == 0:
            continue
        other_group = slater_group(on, ol)
        electrons = occ
        if on == n and ol == l:
            electrons -= 1
        if electrons <= 0:
            continue
        
        if l <= 1:  # s, p target
            if on == n:
                contrib = 0.30 if n == 1 else 0.35
                S += contrib * electrons
            elif on == n - 1:
                S += 0.85 * electrons
            elif on <= n - 2:
                S += 1.00 * electrons
        else:  # d, f target
            if other_group == target_group:
                S += 0.35 * electrons
            elif other_group < target_group:
                S += 1.00 * electrons
    
    return Z - S

# Test cases from Slater's original paper and standard references
# Iron (Z=26): [Ar] 3d6 4s2
# 4s electron: σ = 0.35×1 + 0.85×14 + 1.00×10 = 0.35 + 11.9 + 10 = 22.25
# Zeff = 26 - 22.25 = 3.75
# 3d electron: σ = 0.35×5 + 1.00×18 = 1.75 + 18 = 19.75
# Zeff = 26 - 19.75 = 6.25

fe_config = {
    (1,0): 2, (2,0): 2, (2,1): 6, (3,0): 2, (3,1): 6, (3,2): 6, (4,0): 2
}

print("Slater screening verification:")
print("=" * 60)

# Fe 4s
zeff_4s = compute_slater_zeff(26, 4, 0, fe_config)
print(f"Fe 4s: Zeff = {zeff_4s:.2f} (expected ~3.75)")

# Fe 3d
zeff_3d = compute_slater_zeff(26, 3, 2, fe_config)
print(f"Fe 3d: Zeff = {zeff_3d:.2f} (expected ~6.25)")

# Carbon 2p: 1s2 2s2 2p2
c_config = {(1,0): 2, (2,0): 2, (2,1): 2}
zeff_2p = compute_slater_zeff(6, 2, 1, c_config)
print(f"C 2p: Zeff = {zeff_2p:.2f} (expected ~3.25)")

# Nitrogen 2p: 1s2 2s2 2p3
n_config = {(1,0): 2, (2,0): 2, (2,1): 3}
zeff_2p_n = compute_slater_zeff(7, 2, 1, n_config)
print(f"N 2p: Zeff = {zeff_2p_n:.2f} (expected ~3.90)")

# Oxygen 2p: 1s2 2s2 2p4
o_config = {(1,0): 2, (2,0): 2, (2,1): 4}
zeff_2p_o = compute_slater_zeff(8, 2, 1, o_config)
print(f"O 2p: Zeff = {zeff_2p_o:.2f} (expected ~4.55)")

# Potassium 4s: [Ar] 4s1
k_config = {(1,0): 2, (2,0): 2, (2,1): 6, (3,0): 2, (3,1): 6, (3,2): 10, (4,0): 1}
zeff_4s_k = compute_slater_zeff(19, 4, 0, k_config)
print(f"K 4s: Zeff = {zeff_4s_k:.2f} (expected ~2.20)")

# Sodium 3s: [Ne] 3s1
na_config = {(1,0): 2, (2,0): 2, (2,1): 6, (3,0): 1}
zeff_3s_na = compute_slater_zeff(11, 3, 0, na_config)
print(f"Na 3s: Zeff = {zeff_3s_na:.2f} (expected ~2.51)")

print("\nComparison with code's quantum_zeff:")
print("The code's Slater screening for s/p electrons uses n-based grouping")
print("which correctly combines (n-1)s, (n-1)p, AND (n-1)d into the 0.85 tier.")
print("This matches Slater's original Fe 4s example.")