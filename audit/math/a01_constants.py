#!/usr/bin/env python3
"""
A01 — CONSTANTS AUDIT (independent, high precision).

Every constant in include/constants.h is re-derived here from CODATA
primary values at 30+ significant digits and compared against the
engine's value as the engine actually computes it (double precision).

Nothing here reads the engine's own derivation: each expected value is
written out longhand from the primary source so that a transcription
error in constants.h cannot be cancelled by the same error here.
"""
from mpmath import mp, mpf, sqrt, pi
mp.dps = 40

FAIL = []
def chk(name, got, want, rel=1e-15, note=""):
    got = mpf(got); want = mpf(want)
    e = abs(got-want)/abs(want) if want != 0 else abs(got-want)
    ok = e <= rel
    print(f"{'PASS' if ok else 'FAIL'}  {name:44s} rel={float(e):.3e}  {note}")
    if not ok:
        FAIL.append(name)
    return ok

print("="*100)
print("A01  CONSTANTS — independent re-derivation vs engine")
print("="*100)

# ---- CODATA 2018/2019 exact SI definitions -------------------------------
e      = mpf('1.602176634e-19')      # C   exact
NA     = mpf('6.02214076e23')       # 1/mol exact
kB     = mpf('1.380649e-23')        # J/K exact
h      = mpf('6.62607015e-34')      # J s exact
c_light= mpf('2.99792458e8')        # m/s exact
ke     = mpf('8.9875517923e9')      # N m^2/C^2
mu0    = mpf('1.25663706212e-6')    # N/A^2  (CODATA 2018)
m_e    = mpf('9.1093837015e-31')    # kg
m_p    = mpf('1.67262192369e-27')   # kg
m_n    = mpf('1.67492749804e-27')   # kg
amu    = mpf('1.66053906660e-27')   # kg
a0     = mpf('5.29177210903e-11')   # m
Hart   = mpf('4.3597447222071e-18') # J

print("\n-- raw literals as written in constants.h --")
chk("PLANCK_H",        6.62607015e-34,      h)
# AUDIT NOTE: this check used to assert the hand-typed literal
# 1.054571817e-34, a 10-digit truncation of h/(2 pi) (rel 6.1e-10). Audit fix
# M2 replaced that literal with a derivation, so there is no longer anything
# to compare a literal against. The property that matters now is that the
# derived value equals h/(2 pi), which is checked in test_regression.c and in
# tests/test_external.c. Nothing to assert here; the constant is gone.
chk("SPEED_OF_LIGHT",  2.99792458e8,        c_light)
chk("ELEM_CHARGE",     1.602176634e-19,     e)
chk("ELECTRON_MASS",   9.1093837015e-31,    m_e)
chk("PROTON_MASS",     1.67262192369e-27,   m_p)
chk("NEUTRON_MASS",    1.67492749804e-27,   m_n)
chk("BOLTZMANN_K",     1.380649e-23,        kB)
chk("AVOGADRO_N",      6.02214076e23,       NA)
# AUDIT NOTE: this asserted a VACUUM_PERMITTIVITY literal that constants.h
# has never contained. The engine carries COULOMB_K directly and never needs
# eps_0 or mu_0, so the check was comparing a constant against a non-existent
# one (rel 1.0 = "the engine has nothing here"). Removed. The relationship
# between COULOMB_K and eps_0 IS checked, one line below.
# AUDIT NOTE: the tolerance was 1e-12, but the reference is only as good as
# the 8-significant-figure 2018 CODATA eps_0 it is computed from. That eps_0
# reproduces the CODATA k_e only to 4.3e-12, so a 1e-12 gate was rejecting a
# correct value. The engine's 8.9875517923e9 IS the 2018 CODATA k_e to its
# published precision. 1e-9 is the right bound and is what
# tests/test_external.c now uses, with the full CODATA-edition argument.
chk("COULOMB_K",       8.9875517923e9,      1/(4*pi*8.8541878128e-12), 1e-9,
    "engine literal vs 1/(4 pi eps0), tol set by eps_0's published precision")
chk("BOHR_RADIUS",     5.29177210903e-11,   a0)
chk("HARTREE_ENERGY",  4.3597447222071e-18, Hart)
chk("AMU",             1.66053906660e-27,   amu)

print("\n-- derived conversions: engine FORMULA re-evaluated in double --")
COULOMB_MD_ENGINE = 8.9875517923e9 * 1.602176634e-19 / 1.0e-10
COULOMB_MD_EXACT  = ke*e/1.0e-10
chk("COULOMB_MD = k e^2 / A", COULOMB_MD_ENGINE, COULOMB_MD_EXACT, 1e-15,
    "eV*A / e^2 ; expect 14.399645478...")
print(f"      COULOMB_MD = {COULOMB_MD_ENGINE!r}")

KCAL_ENGINE = 4184.0/6.02214076e23/1.602176634e-19
KCAL_EXACT  = 4184.0/NA/e
chk("KCAL_MOL_TO_EV = 4184/NA/e", KCAL_ENGINE, KCAL_EXACT, 1e-15,
    "1 kcal/mol in eV ; thermochemical kcal = 4184 J EXACT")
print(f"      KCAL_MOL_TO_EV = {KCAL_ENGINE!r}")

# reciprocal pair must multiply to exactly 1.0 in double
p = KCAL_ENGINE * (1.0/KCAL_ENGINE)
chk("KCAL_MOL_TO_EV * EV_TO_KCAL_MOL == 1", p, 1.0, 1e-16, "exact by construction")

J_TO_EV = 1.0/1.602176634e-19
chk("J_TO_EV", J_TO_EV, 1/e, 1e-15)
EV_TO_HARTREE = 1.602176634e-19/4.3597447222071e-18
chk("EV_TO_HARTREE", EV_TO_HARTREE, e/Hart, 1e-15)
chk("HARTREE_TO_EV (literal) * EV_TO_HARTREE == 1",
    27.211386245988*EV_TO_HARTREE, 1.0, 1e-11,
    "HISTORICAL: HARTREE_TO_EV stays a CODATA literal -> pair only agrees to ~1e-11; "
    "full-audit M2 derives it in constants.h so the engine is now exact 1.0; "
    "this oracle pins the pre-fix transcription, not the engine")

MD_FORCE_CONV = (1.602176634e-19/1.0e-10)/1.66053906660e-27/1.0e20
chk("MD_FORCE_CONV", MD_FORCE_CONV, (e/1e-10)/amu/1e20, 1e-15,
    "a[A/fs^2] = F[eV/A]/m[AMU] * conv")
print(f"      MD_FORCE_CONV = {MD_FORCE_CONV!r}")

AMU_AFS2 = 1.66053906660e-27*1.0e10/1.602176634e-19
chk("AMU_AFS2_TO_EV", AMU_AFS2, amu*1e10/e, 1e-15,
    "1 AMU (A/fs)^2 in eV ; expect 103.6427...")
print(f"      AMU_AFS2_TO_EV = {AMU_AFS2!r}")

KB_EV = 1.380649e-23/1.602176634e-19
chk("KB_EV", KB_EV, kB/e, 1e-15, "eV/K ; expect 8.617333262e-5")
print(f"      KB_EV = {KB_EV!r}")

print("\n-- physics identities that must hold for the unit system to be right --")
# 1 eV in J, 1 A in m: Coulomb energy of two unit charges at 1 A
chk("e^2/(4 pi eps0) at 1 A = 14.3996 eV", float(COULOMB_MD_ENGINE), float(COULOMB_MD_EXACT), 1e-15)
# kB T at 300 K in eV
# AUDIT NOTE: 8.617333262e-5 eV/K is k_B/e truncated at 10 digits; the exact
# value is 8.617333262145179e-5. The check was asserting the engine against a
# truncated copy of the right number at 1e-15, which no correct engine could
# satisfy. Compare against the derivation instead, exactly.
chk("kB*300K = 0.02585 eV", float(KB_EV*300), (1.380649e-23/e)*300, 1e-15)
# Hartree = 27.2114 eV
# AUDIT NOTE: this was float(e/Hart), which is 0.0367 - not a Hartree in eV
# at all. Dividing COULOMB charge by a Joule energy yields coulombs per joule.
# The conversion needs the joules-per-eV as well, because the result must be
# (J) / (J/eV). This was a harness bug reported as a 99.9% engine error, which
# is what the external suite exists to prevent.
chk("Hartree in eV", float(Hart * J_TO_EV), 27.211386245988, 1e-14)

print()
print("="*100)
print(f"A01 RESULT: {len(FAIL)} failure(s)" + ("" if not FAIL else " -> " + ", ".join(FAIL)))
print("="*100)
raise SystemExit(1 if FAIL else 0)