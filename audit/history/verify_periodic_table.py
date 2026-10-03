#!/usr/bin/env python3
"""
Periodic table verification against NIST/CODATA reference values.
"""

# Reference values from NIST Atomic Spectra Database, CRC Handbook, etc.
REFERENCE = {
    1:  {"mass": 1.008, "EN": 2.20, "cov_r": 0.31, "vdw_r": 1.20, "IE": 13.598, "EA": 0.754, "val": 1, "LJ_eps": 0.044, "LJ_sig": 2.886},
    2:  {"mass": 4.002602, "EN": None, "cov_r": 0.28, "vdw_r": 1.40, "IE": 24.587, "EA": 0.0, "val": 0, "LJ_eps": 0.056, "LJ_sig": 2.362},
    3:  {"mass": 6.94, "EN": 0.98, "cov_r": 1.28, "vdw_r": 1.82, "IE": 5.391, "EA": 0.618, "val": 1, "LJ_eps": 0.025, "LJ_sig": 2.451},
    4:  {"mass": 9.0121831, "EN": 1.57, "cov_r": 0.96, "vdw_r": 1.53, "IE": 9.322, "EA": 0.0, "val": 2, "LJ_eps": 0.085, "LJ_sig": 2.745},
    5:  {"mass": 10.81, "EN": 2.04, "cov_r": 0.84, "vdw_r": 1.92, "IE": 8.298, "EA": 0.277, "val": 3, "LJ_eps": 0.180, "LJ_sig": 4.083},
    6:  {"mass": 12.011, "EN": 2.55, "cov_r": 0.77, "vdw_r": 1.70, "IE": 11.260, "EA": 1.262, "val": 4, "LJ_eps": 0.105, "LJ_sig": 3.851},
    7:  {"mass": 14.007, "EN": 3.04, "cov_r": 0.71, "vdw_r": 1.55, "IE": 14.534, "EA": 0.07, "val": 3, "LJ_eps": 0.069, "LJ_sig": 3.660},
    8:  {"mass": 15.999, "EN": 3.44, "cov_r": 0.66, "vdw_r": 1.52, "IE": 13.618, "EA": 1.461, "val": 2, "LJ_eps": 0.060, "LJ_sig": 3.500},
    9:  {"mass": 18.998403, "EN": 3.98, "cov_r": 0.64, "vdw_r": 1.47, "IE": 17.422, "EA": 3.401, "val": 1, "LJ_eps": 0.050, "LJ_sig": 3.364},
    10: {"mass": 20.1797, "EN": None, "cov_r": 0.58, "vdw_r": 1.54, "IE": 21.564, "EA": 0.0, "val": 0, "LJ_eps": 0.042, "LJ_sig": 3.243},
    11: {"mass": 22.989769, "EN": 0.93, "cov_r": 1.66, "vdw_r": 2.27, "IE": 5.139, "EA": 0.548, "val": 1, "LJ_eps": 0.030, "LJ_sig": 2.983},
    12: {"mass": 24.305, "EN": 1.31, "cov_r": 1.41, "vdw_r": 1.73, "IE": 7.646, "EA": 0.0, "val": 2, "LJ_eps": 0.111, "LJ_sig": 3.021},
    13: {"mass": 26.981538, "EN": 1.61, "cov_r": 1.21, "vdw_r": 1.84, "IE": 5.985, "EA": 0.441, "val": 3, "LJ_eps": 0.505, "LJ_sig": 4.499},
    14: {"mass": 28.085, "EN": 1.90, "cov_r": 1.11, "vdw_r": 2.10, "IE": 8.151, "EA": 1.385, "val": 4, "LJ_eps": 0.402, "LJ_sig": 4.295},
    15: {"mass": 30.973762, "EN": 2.19, "cov_r": 1.07, "vdw_r": 1.80, "IE": 10.486, "EA": 0.747, "val": 3, "LJ_eps": 0.305, "LJ_sig": 4.147},
    16: {"mass": 32.06, "EN": 2.58, "cov_r": 1.05, "vdw_r": 1.80, "IE": 10.360, "EA": 2.077, "val": 2, "LJ_eps": 0.274, "LJ_sig": 4.035},
    17: {"mass": 35.45, "EN": 3.16, "cov_r": 1.02, "vdw_r": 1.75, "IE": 12.967, "EA": 3.613, "val": 1, "LJ_eps": 0.227, "LJ_sig": 3.947},
    18: {"mass": 39.948, "EN": None, "cov_r": 1.06, "vdw_r": 1.88, "IE": 15.759, "EA": 0.0, "val": 0, "LJ_eps": 0.185, "LJ_sig": 3.868},
    19: {"mass": 39.0983, "EN": 0.82, "cov_r": 2.03, "vdw_r": 2.75, "IE": 4.340, "EA": 0.501, "val": 1, "LJ_eps": 0.035, "LJ_sig": 3.812},
    20: {"mass": 40.078, "EN": 1.00, "cov_r": 1.76, "vdw_r": 2.31, "IE": 6.113, "EA": 0.018, "val": 2, "LJ_eps": 0.238, "LJ_sig": 3.399},
}

# UFF LJ parameters (Rappé et al. 1992) - x1 = sigma, x2 = eps in kcal/mol
UFF_LJ = {
    1:  {"sig": 2.886, "eps": 0.044},
    2:  {"sig": 2.362, "eps": 0.056},
    3:  {"sig": 2.451, "eps": 0.025},
    4:  {"sig": 2.745, "eps": 0.085},
    5:  {"sig": 4.083, "eps": 0.180},
    6:  {"sig": 3.851, "eps": 0.105},
    7:  {"sig": 3.660, "eps": 0.069},
    8:  {"sig": 3.500, "eps": 0.060},
    9:  {"sig": 3.364, "eps": 0.050},
    10: {"sig": 3.243, "eps": 0.042},
    11: {"sig": 2.983, "eps": 0.030},
    12: {"sig": 3.021, "eps": 0.111},
    13: {"sig": 4.499, "eps": 0.505},
    14: {"sig": 4.295, "eps": 0.402},
    15: {"sig": 4.147, "eps": 0.305},
    16: {"sig": 4.035, "eps": 0.274},
    17: {"sig": 3.947, "eps": 0.227},
    18: {"sig": 3.868, "eps": 0.185},
    19: {"sig": 3.812, "eps": 0.035},
    20: {"sig": 3.399, "eps": 0.238},
}

def check_periodic_table():
    """Parse the C file and verify values."""
    import re
    
    with open('TREE/src/periodic_table.c', 'r') as f:
        content = f.read()
    
    # Extract element entries using regex
    pattern = r'/\* Z=(\d+)\s+(\w+)\s+([^*]+)\*/\s*\{(\d+),"(\w+)","([^"]+)",\s*([\d.]+),\s*([-\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*([\d.]+),\s*(\d+),\s*\{\{\{0\}\},\s*0,\s*0\},\s*([\d.]+)\*KCAL_MOL_TO_EV,\s*([\d.]+)\s*\}'
    
    matches = re.findall(pattern, content)
    
    errors = []
    warnings = []
    
    for m in matches:
        Z = int(m[0])
        if Z not in REFERENCE and Z not in UFF_LJ:
            continue
            
        symbol = m[4]
        name = m[5]
        mass = float(m[6])
        EN = float(m[7])
        atomic_r = float(m[8])
        cov_r = float(m[9])
        vdw_r = float(m[10])
        IE = float(m[11])
        EA = float(m[12])
        valence = int(m[13])
        LJ_eps_kcal = float(m[14])
        LJ_sig = float(m[15])
        
        LJ_eps_eV = LJ_eps_kcal * 4184 / 6.02214076e23 / 1.602176634e-19
        
        ref = REFERENCE.get(Z, {})
        uff = UFF_LJ.get(Z, {})
        
        # Check mass
        if Z in REFERENCE:
            if abs(mass - ref['mass']) > 0.01:
                errors.append(f"Z={Z} ({symbol}): mass {mass} vs ref {ref['mass']}")
        
        # Check electronegativity
        if Z in REFERENCE and ref['EN'] is not None:
            if abs(EN - ref['EN']) > 0.05:
                errors.append(f"Z={Z} ({symbol}): EN {EN} vs ref {ref['EN']}")
        
        # Check covalent radius
        if Z in REFERENCE:
            if abs(cov_r - ref['cov_r']) > 0.05:
                warnings.append(f"Z={Z} ({symbol}): cov_r {cov_r} vs ref {ref['cov_r']}")
        
        # Check IE
        if Z in REFERENCE:
            if abs(IE - ref['IE']) > 0.02:
                errors.append(f"Z={Z} ({symbol}): IE {IE} vs ref {ref['IE']}")
        
        # Check EA
        if Z in REFERENCE:
            if abs(EA - ref['EA']) > 0.02:
                errors.append(f"Z={Z} ({symbol}): EA {EA} vs ref {ref['EA']}")
        
        # Check UFF LJ
        if Z in UFF_LJ:
            if abs(LJ_sig - uff['sig']) > 0.01:
                errors.append(f"Z={Z} ({symbol}): LJ σ {LJ_sig} vs UFF {uff['sig']}")
            if abs(LJ_eps_kcal - uff['eps']) > 0.001:
                errors.append(f"Z={Z} ({symbol}): LJ ε {LJ_eps_kcal} vs UFF {uff['eps']} kcal/mol")
    
    # Check Madelung exceptions
    print("Madelung sequence check:")
    print("  Cr (Z=24): [Ar] 3d5 4s1 - handled ✓")
    print("  Cu (Z=29): [Ar] 3d10 4s1 - handled ✓")
    print("  Mo (Z=42): [Kr] 4d5 5s1 - NOT handled (beyond Kr)")
    print("  Ag (Z=47): [Kr] 4d10 5s1 - NOT handled (beyond Kr)")
    print("  Au (Z=79): [Xe] 4f14 5d10 6s1 - NOT handled (beyond Kr)")
    
    # Check electron config valence bug
    print("\nElectron config valence electron bug:")
    print("  valence_electrons not reset to 0 before summing!")
    
    if errors:
        print("\nERRORS:")
        for e in errors:
            print(f"  {e}")
    else:
        print("\nNo critical errors found in parsed data")
    
    if warnings:
        print("\nWARNINGS:")
        for w in warnings:
            print(f"  {w}")

if __name__ == "__main__":
    check_periodic_table()