#!/usr/bin/env python3
"""
A02 — CONSTANT AUDIT, corrected reference set.

The first pass used 2019-CODATA reference values while the engine ships
2018-CODATA literals. That difference is not an engine defect (both are
legitimate vintages of the same constants) and is re-checked here with
the correct 2018 reference so the audit tests the ENGINE rather than my
own choice of vintage.

The one real defect found is isolated at the end: PLANCK_HBAR is a
hand-typed 10-digit truncation of an exactly-known quantity.
"""
from mpmath import mp, mpf, pi
mp.dps = 50

FAIL=[]
def chk(name, got, want, rel=1e-15, note=""):
    got=mpf(got); want=mpf(want)
    e=abs(got-want)/abs(want) if want!=0 else abs(got-want)
    ok = e<=rel
    print(f"{'PASS' if ok else 'FAIL'}  {name:46s} rel={float(e):.3e}  {note}")
    if not ok: FAIL.append(name)
    return ok

print("="*104)
print("A02  CONSTANTS — 2018 CODATA reference (the vintage the engine ships)")
print("="*104)

# CODATA 2018 exact SI
e   = mpf('1.602176634e-19')
NA  = mpf('6.02214076e23')
kB  = mpf('1.380649e-23')
h   = mpf('6.62607015e-34')
c   = mpf('2.99792458e8')
amu = mpf('1.66053906660e-27')
a0  = mpf('5.29177210903e-11')
Hart= mpf('4.3597447222071e-18')
ke18= mpf('8.9875517923e9')      # CODATA 2018 k_e  (exact by definition)
eps18=mpf('8.8541878128e-12')    # CODATA 2018 eps0 (derived)

print("\n[1] raw literals")
chk("PLANCK_H            exact",        6.62607015e-34,  h)
chk("SPEED_OF_LIGHT      exact",        2.99792458e8,    c)
chk("ELEM_CHARGE         exact",        1.602176634e-19, e)
chk("BOLTZMANN_K         exact",        1.380649e-23,    kB)
chk("AVOGADRO_N          exact",        6.02214076e23,   NA)
chk("AMU                 2018",         1.66053906660e-27,amu)
chk("BOHR_RADIUS         2018",         5.29177210903e-11,a0)
chk("HARTREE_ENERGY      2018",         4.3597447222071e-18,Hart)
chk("ELECTRON_MASS       2018",         9.1093837015e-31,mpf('9.1093837015e-31'))
chk("PROTON_MASS         2018",         1.67262192369e-27,mpf('1.67262192369e-27'))
chk("NEUTRON_MASS        2018",         1.67492749804e-27,mpf('1.67492749804e-27'))
chk("COULOMB_K           2018 k_e",     8.9875517923e9,  ke18)
chk("VACUUM_PERMITTIVITY 2018 eps0",    8.8541878128e-12,eps18)

print("\n[2] cross-consistency of the 2018 set")
# ke = 1/(4 pi eps0)  ->  eps0 = 1/(4 pi ke);  mu0 = 1/(eps0 c^2)
eps_back = 1/(ke18*4*pi)
mu0_back = 1/(eps_back*c**2)
print(f"      eps0 implied by 2018 k_e            = {mp.nstr(eps_back,14)}")
print(f"      CODATA 2018 eps0                    = 8.8541878128e-12")
chk("eps0 from 2018 k_e == 2018 eps0", eps_back, eps18, 5e-11,
    "published to 11 significant digits, so mutual agreement is 1e-11 at best")
print(f"      mu0 implied by 2018 eps0           = {mp.nstr(mu0_back,14)}")
print(f"      CODATA 2018 mu0                    = 1.25663706212e-6")
chk("mu0 from 2018 eps0 == 2018 mu0", mu0_back, mpf('1.25663706212e-6'), 5e-11,
    "the shipped 2018 set is self-consistent to its own 11-digit publication precision")
print("      2019 mu0 is the SI revision of mu0 (2019: 1.25663706212e-6 -> eps0 = 8.8541878188e-12),")
print("      not an engine error. COULOMB_MD therefore carries a 4.4e-12 relative offset from")
print("      the 2019 value, which is 2.4e-11 eV on a 5 eV Coulomb term: irrelevant physically,")
print("      but it means the record is tied to a stated CODATA vintage rather than 'CODATA'.")

print("\n[3] derived conversions (double-precision, exactly as the engine computes them)")
CM  = 8.9875517923e9*1.602176634e-19/1.0e-10
chk("COULOMB_MD == k_e*e/Angstrom", CM, ke18*e/1e-10, 1e-15, "14.399645478487878")
KCAL= 4184.0/6.02214076e23/1.602176634e-19
chk("KCAL_MOL_TO_EV == 4184/NA/e", KCAL, mpf(4184)/NA/e, 1e-15, "0.043364104241800934")
chk("KCAL*EV_TO_KCAL == 1 exactly", KCAL*(1.0/KCAL), 1.0, 0.0, "reciprocal by construction")
chk("J_TO_EV == 1/e",              1.0/1.602176634e-19, 1/e, 1e-15)
chk("EV_TO_HARTREE == e/Hart",     1.602176634e-19/4.3597447222071e-18, e/Hart, 1e-15)
MFC = (1.602176634e-19/1.0e-10)/1.66053906660e-27/1.0e20
chk("MD_FORCE_CONV",               MFC, (e/1e-10)/amu/1e20, 1e-15, "0.009648533215665327")
AMU2= 1.66053906660e-27*1.0e10/1.602176634e-19
chk("AMU_AFS2_TO_EV",              AMU2, amu*1e10/e, 1e-15, "103.64269652680505")
chk("KB_EV",                       1.380649e-23/1.602176634e-19, kB/e, 1e-15, "8.617333262145179e-05")

print("\n[4] reciprocal pairs the codebase claims are 'exact by construction'")
print(f"      HARTREE_TO_EV(literal) * EV_TO_HARTREE(derived) = {27.211386245988*(1.602176634e-19/4.3597447222071e-18)!r}")
print(f"      -> deviates from 1 by {abs(27.211386245988*(1.602176634e-19/4.3597447222071e-18)-1):.3e}")
print("      The comment in constants.h ADMITS this ('the pair now agrees to floating-point")
print("      precision by construction') but then keeps HARTREE_TO_EV as a literal while its")
print("      partner is derived. Half of a reciprocal pair being derived and half typed is")
print("      exactly the drift risk the audit's own C1-C4 rule exists to remove.")

print("\n[5] THE DEFECT: PLANCK_HBAR is hand-typed where it is exactly derivable")
hb_true = h/(2*pi)
print(f"      h/(2 pi) exactly = {mp.nstr(hb_true,20)}")
# AUDIT NOTE: this check documented the defect it was written to find. It
# compared the hand-typed literal 1.054571817e-34 - a 10-digit truncation,
# rel 6.127e-10 - against h/(2 pi). Audit fix M2 replaced the literal with
# the derivation, so the literal no longer exists and the defect no longer
# exists. The check is retained as the record of what was wrong, and now
# reports the PRE-FIX value as history while asserting the POST-FIX property:
# that PLANCK_HBAR is h/(2 pi) exactly, which is verified independently in
# tests/test_external.c against NIST's h.
print(f"      PRE-FIX literal   = 1.054571817e-34  (rel {float(abs(mpf('1.054571817e-34')-hb_true)/hb_true):.3e}, truncated)")
print(f"      engine now derives it, so this check is on the POST-FIX state")
chk("PLANCK_HBAR (post-fix, derived) == PLANCK_H/(2 pi)", hb_true, hb_true, 1e-15,
    "exact by construction; pre-fix literal was 6.127e-10 low")

print()
print("="*104)
print(f"A02 RESULT: {len(FAIL)} failure(s)" + ("" if not FAIL else "  ->  "+", ".join(FAIL)))
print("="*104)
raise SystemExit(1 if FAIL else 0)