#!/usr/bin/env python3
"""HISTORICAL (pre-v9R4), superseded; DO NOT RUN. Kept for provenance only.
See v9R4/verify_scripts.sh and audit/run_audit.sh for current verifiers.
Running this file today would re-mutate a fixed tree.
"""
"""
Add static assertions to types.h
"""

with open('TREE/include/types.h', 'r') as f:
    content = f.read()

old = '''/* ═══════════════════════════════════════════════════════════════════════════
 * Return codes
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    SIM_OK               =  0,
    SIM_ERR_ALLOC        = -1,
    SIM_ERR_BADATOM      = -2,
    SIM_ERR_BADBOND      = -3,
    SIM_ERR_OVERFLOW     = -4,
    SIM_ERR_BADPARAM     = -5
} SimError;

#endif /* TYPES_H */'''

new = '''/* ════════════════════════════════════════════════════════════════════════════
 * Compile-time assertions
 * ══════════════════════════════════════════════════════════════════════════ */

STATIC_ASSERT(sizeof(Vec3) == 3 * sizeof(double), "Vec3 must be exactly 3 doubles");
STATIC_ASSERT(sizeof(QuantumNumbers) == 3 * sizeof(int) + sizeof(double), "QuantumNumbers layout");
STATIC_ASSERT(sizeof(Orbital) >= sizeof(QuantumNumbers) + sizeof(double) + sizeof(int), "Orbital layout");
STATIC_ASSERT(sizeof(ElectronConfig) == MAX_SHELLS * 4 * sizeof(int) + 2 * sizeof(int), "ElectronConfig layout");
STATIC_ASSERT(sizeof(Element) >= 8 * sizeof(double) + 4 * sizeof(int) + 36 + sizeof(ElectronConfig), "Element layout");
STATIC_ASSERT(sizeof(Atom) >= 3 * sizeof(Vec3) + 6 * sizeof(double) + 4 * sizeof(int) + sizeof(ElectronConfig) + MAX_ORBITALS * sizeof(Orbital) + MAX_BONDS_PER_ATOM * 2 * sizeof(int), "Atom layout");
STATIC_ASSERT(sizeof(Bond) == 4 * sizeof(int) + 2 * sizeof(double), "Bond layout");
STATIC_ASSERT(sizeof(Angle) == 3 * sizeof(int) + 2 * sizeof(double), "Angle layout");
STATIC_ASSERT(sizeof(Dihedral) == 4 * sizeof(int) + 3 * sizeof(double), "Dihedral layout");
STATIC_ASSERT(sizeof(SimBox) == sizeof(Vec3) + 3 * sizeof(int), "SimBox layout");
STATIC_ASSERT(sizeof(Thermostat) == sizeof(int) + 4 * sizeof(double), "Thermostat layout");
STATIC_ASSERT(sizeof(Simulation) >= 7 * sizeof(void*) + 10 * sizeof(int) + 10 * sizeof(double) + sizeof(SimBox) + sizeof(Thermostat) + sizeof(uint64_t), "Simulation layout");

/* Verify enum values fit in int */
STATIC_ASSERT(SIM_OK == 0, "SIM_OK must be 0");
STATIC_ASSERT(SIM_ERR_ALLOC < 0, "Error codes must be negative");
STATIC_ASSERT(THERMOSTAT_NONE == 0, "THERMOSTAT_NONE must be 0");
STATIC_ASSERT(QN_L_S == 0 && QN_L_P == 1 && QN_L_D == 2 && QN_L_F == 3, "Quantum number l values");

/* Verify constants */
STATIC_ASSERT(MAX_ELECTRONS == 128, "MAX_ELECTRONS");
STATIC_ASSERT(MAX_ATOMS == 100000, "MAX_ATOMS");
STATIC_ASSERT(MAX_BONDS == 200000, "MAX_BONDS");
STATIC_ASSERT(MAX_BONDS_PER_ATOM == 8, "MAX_BONDS_PER_ATOM");
STATIC_ASSERT(MAX_SHELLS == 7, "MAX_SHELLS");
STATIC_ASSERT(MAX_ELEMENTS == 118, "MAX_ELEMENTS");
STATIC_ASSERT(MAX_ORBITALS == 32, "MAX_ORBITALS");

/* ════════════════════════════════════════════════════════════════════════════
 * Return codes
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    SIM_OK               =  0,
    SIM_ERR_ALLOC        = -1,
    SIM_ERR_BADATOM      = -2,
    SIM_ERR_BADBOND      = -3,
    SIM_ERR_OVERFLOW     = -4,
    SIM_ERR_BADPARAM     = -5
} SimError;

#endif /* TYPES_H */'''

content = content.replace(old, new)

with open('TREE/include/types.h', 'w') as f:
    f.write(content)

print("Static assertions added to types.h!")