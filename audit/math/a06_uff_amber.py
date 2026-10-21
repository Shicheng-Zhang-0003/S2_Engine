#!/usr/bin/env python3
"""
A01b: cross-check the engine's UFF Lennard-Jones table and its AMBER ff99
Lennard-Jones table against INDEPENDENT transcriptions of the primary
sources.

FULL-AUDIT: this script used to live in audit/ and was therefore never
executed by `make audit` (which globs audit/math/*.py and
audit/external/*.py). It also carried two defects of its own, so that
running it produced 20 red findings that were all artefacts of the
SCRIPT, not of the engine:

  1. its hard-coded UFF table was the LAMMPS `lj/cut/uff` transcription,
     not Rappé Table I. The two differ for 18 of the 36 elements (the
     well-known cases are Si: x1 = 4.295 in Rappé against 4.195 in
     LAMMPS, and Ar: 3.868 against 3.912). The ENGINE follows Rappé
     Table I, which is the primary source, so the script was checking
     the engine against the wrong table;
  2. its regular expression only matched element rows whose electron
     configuration was literally `{{{0}}, 0, 0`, so every element from
     Sc upward was reported as "row unparseable".

This version:
  - parses the ENGINE's rows structurally (not by one fixed electron
    configuration),
  - takes the reference from audit/external/uff_reference.prm, the Open
    Babel UFF parameter file, whose header cites Rappé et al., JACS 114
    (1992) 10024 and whose x1/D1 columns reproduce Table I,
  - and checks the identity that actually matters: sigma = x1 / 2^(1/6)
    and epsilon = D1 (kcal/mol) for all 36 elements.

The AMBER half is unchanged in substance: it re-derives the engine's
13 atom types from audit/external/parm99.dat.
"""
import math, os, re, sys

FAIL = []


def chk(name, got, ref, tol=0.0):
    ok = abs(got - ref) <= tol * max(1.0, abs(ref)) if tol else (got == ref)
    if not ok:
        FAIL.append(f"{name}: got {got!r} ref {ref!r}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:44s} got={got!r:<22} ref={ref!r}")


HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.dirname(HERE)
ROOT = os.path.dirname(AUDIT)
TW = 2.0 ** (1.0 / 6.0)
K = 4184.0 / 6.02214076e23 / 1.602176634e-19   # kcal/mol -> eV (per particle)


# ─────────────────────────────────────────────────────────── UFF reference
# element symbol -> the UFF atom type in the Open Babel parameter file
UFF_TYPE = {
    "H": "H_", "He": "He4+4", "Li": "Li", "Be": "Be3+2", "B": "B_3",
    "C": "C_3", "N": "N_3", "O": "O_3", "F": "F_", "Ne": "Ne4+4",
    "Na": "Na", "Mg": "Mg3+2", "Al": "Al3", "Si": "Si3", "P": "P_3+3",
    "S": "S_3+2", "Cl": "Cl", "Ar": "Ar4+4", "K": "K_", "Ca": "Ca6+2",
    "Sc": "Sc3+3", "Ti": "Ti3+4", "V": "V_3+5", "Cr": "Cr6+3",
    "Mn": "Mn6+2", "Fe": "Fe3+2", "Co": "Co6+3", "Ni": "Ni4+2",
    "Cu": "Cu3+1", "Zn": "Zn3+2", "Ga": "Ga3+3", "Ge": "Ge3",
    "As": "As3+3", "Se": "Se3+2", "Br": "Br", "Kr": "Kr4+4",
}

print("== UFF (Rappe 1992) sigma/epsilon vs periodic_table.c ==")
prm = open(os.path.join(AUDIT, "external", "uff_reference.prm")).read()
ref = {}
for m in re.finditer(r"^param\s+(\S+)\s+([0-9.]+)\s+([0-9.]+)\s+"
                     r"([0-9.]+)\s+([0-9.]+)\s+", prm, re.M):
    # param <type> r1 theta0 x1 D1 zeta ... : x1 is the 4th field,
    # D1 (the LJ well depth) the 5th.
    ref.setdefault(m.group(1), (float(m.group(4)), float(m.group(5))))
missing = [t for t in UFF_TYPE.values() if t not in ref]
if missing:
    FAIL.append(f"UFF reference file is missing types: {missing}")
    print(f"  [FAIL] reference file is missing {len(missing)} UFF types")
else:
    print(f"  parsed {len(ref)} UFF parameter rows "
          f"({len(UFF_TYPE)} element types referenced)")

# The engine's table: parse every element row structurally.
src = open(os.path.join(ROOT, "src", "periodic_table.c")).read()
rows = re.findall(
    r"\{\s*(\d+)\s*,\s*\"(\w+)\"\s*,\s*\"[^\"]*\"\s*,[^;]*?"
    r"([0-9.]+)\*KCAL_MOL_TO_EV\s*,\s*([0-9.]+)\s*\}", src, re.S)
engine = {int(z): (sym, float(eps), float(sig)) for z, sym, eps, sig in rows}
print(f"  parsed {len(engine)} element rows from src/periodic_table.c")

n_bad = 0
for Z in range(1, 37):
    if Z not in engine:
        FAIL.append(f"UFF Z={Z} row not found in periodic_table.c")
        n_bad += 1
        continue
    sym, eps_lit, sig = engine[Z]
    if sym not in UFF_TYPE:
        FAIL.append(f"UFF Z={Z} symbol {sym} has no reference type")
        n_bad += 1
        continue
    x1, eps_kcal = ref[UFF_TYPE[sym]]
    # The engine's table stores the well depth as a literal in kcal/mol
    # multiplied by KCAL_MOL_TO_EV at compile time, so the comparison is
    # between two kcal/mol numbers and the conversion cancels.
    e_ok = abs(eps_lit - eps_kcal) <= 1e-12
    s_ok = abs(sig - x1 / TW) <= 1e-4
    if not (e_ok and s_ok):
        n_bad += 1
        FAIL.append(f"UFF Z={Z}({sym}): eps {eps_lit} vs {eps_kcal}; "
                    f"sigma {sig} vs x1/2^(1/6)={x1/TW:.5f}")
chk("UFF: 36/36 sigma and epsilon rows", 36 - n_bad, 36, 0)

# The identity that must hold for every element.
worst = 0.0
for Z in range(1, 37):
    sym, eps_lit, sig = engine[Z]
    x1, _ = ref[UFF_TYPE[sym]]
    worst = max(worst, abs(sig * TW - x1))
print(f"  worst |sigma * 2^(1/6) - x1| over 36 elements: {worst:.3e} A")
if worst >= 1e-4:
    FAIL.append(f"UFF sigma identity: worst |sigma*2^(1/6) - x1| = {worst}")
ok = worst < 1e-4
print(f"  [{'PASS' if ok else 'FAIL'}] "
      f"{'UFF: sigma * 2^(1/6) reproduces x1 to 1e-4 A':44s} "
      f"got={worst!r:<22} ref=0.0")

print("== AMBER ff99 parm99.dat MOD4 RE vs include/amber_lj.h ==")
p99 = open(os.path.join(AUDIT, "external", "parm99.dat")).read()
mod4 = p99[p99.find("MOD4      RE"):]
parm = {}
for m in re.finditer(r"^\s{2}([A-Za-z0-9*']+)\s+([0-9.]+)\s+([0-9.]+)", mod4, re.M):
    parm[m.group(1)] = (float(m.group(2)), float(m.group(3)))
print(f"  parsed {len(parm)} MOD4 RE atom types")
amber = open(os.path.join(ROOT, "include", "amber_lj.h")).read()
types = ["N", "C2", "CT", "O", "OH", "OS", "P", "HN", "HC", "H1", "H4",
         "HA", "H5"]
parmkey = {"C2": "C", "HN": "H"}
n_bad2 = 0
for t in types:
    key = parmkey.get(t, t)
    if key not in parm:
        FAIL.append(f"parm99 missing {key}")
        n_bad2 += 1
        continue
    r_ref, e_ref = parm[key]
    # parm99.dat's MOD4 RE gives R* = Rmin/2; the engine stores the sigma
    # form, so the reference sigma is 2*R*/2^(1/6).
    r_ref = 2.0 * r_ref / TW
    sigmac = re.search(r"#define\s+LJ_AMBER_%s_SIGMA\s+AMBER_RSTAR_TO_SIGMA\(([0-9.]+)\)" % t, amber)
    epsc = re.search(r"#define\s+LJ_AMBER_%s_EPS\s+\(([0-9.]+)\s*\*\s*KCAL_MOL_TO_EV\)" % t, amber)
    if not sigmac or not epsc:
        FAIL.append(f"amber_lj.h missing {t}")
        n_bad2 += 1
        continue
    s_got = float(sigmac.group(1)) * 2.0 / TW
    e_got = float(epsc.group(1))
    if abs(s_got - r_ref) > 1e-3 or abs(e_got - e_ref) > 1e-6:
        n_bad2 += 1
        FAIL.append(f"AMBER {t}: sigma {s_got} vs {r_ref}; eps {e_got} vs {e_ref}")
chk("AMBER ff99: 13 atom types match parm99.dat", 13 - n_bad2, 13, 0)

print("== derived sigma identities ==")
chk("2*1.8240/2^(1/6)  (N)", 2 * 1.8240 / TW, 3.25, 1e-4)
chk("2*1.6612/2^(1/6)  (O)", 2 * 1.6612 / TW, 2.96, 1e-4)
chk("2*1.9080/2^(1/6)  (C)", 2 * 1.9080 / TW, 3.39967, 1e-5)
chk("x1(TI 3.175)/2^(1/6)", 3.175 / TW, 2.8286, 1e-4)

print()
if FAIL:
    print(f"A01b: {len(FAIL)} FAILURES")
    for f in FAIL:
        print("  -", f)
    sys.exit(1)
print("A01b: all checks passed")
sys.exit(0)
