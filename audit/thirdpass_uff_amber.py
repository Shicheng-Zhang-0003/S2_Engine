#!/usr/bin/env python3
"""A01b: table verification, fixed parsing against real file formats."""
import math, re, sys
FAIL = []
def chk(name, got, ref, tol):
    ok = abs(got - ref) <= tol * max(1.0, abs(ref))
    if not ok: FAIL.append(f"{name}: got {got!r} ref {ref!r}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:44s} got={got!r:<22} ref={ref!r}")

ROOT = "/home/magi-01/Desktop/work/projects/4179-Magi"
TW = 2.0**(1.0/6.0)

print("== UFF (Rappe 1992) sigma/epsilon vs periodic_table.c ==")
UFF = {  # Z: x1(Rmin A), eps kcal/mol, published sigma
 1:(2.886,.044,2.571), 2:(2.362,.056,2.104), 3:(2.451,.025,2.184), 4:(2.745,.085,2.445),
 5:(4.083,.180,3.637), 6:(3.851,.105,3.431), 7:(3.660,.069,3.260), 8:(3.500,.060,3.118),
 9:(3.364,.050,2.997),10:(3.243,.042,2.889),11:(2.983,.030,2.657),12:(3.021,.111,2.691),
13:(4.499,.505,4.008),14:(4.195,.402,3.737),15:(4.147,.305,3.694),16:(4.035,.274,3.594),
17:(3.947,.227,3.516),18:(3.912,.185,3.485),19:(3.812,.035,3.396),20:(3.399,.238,3.028),
21:(3.295,.056,2.935),22:(3.175,.056,2.828),23:(3.144,.056,2.802),24:(3.023,.056,2.693),
25:(2.961,.056,2.638),26:(2.912,.056,2.594),27:(2.872,.056,2.559),28:(2.834,.056,2.525),
29:(3.495,.024,3.113),30:(2.763,.055,2.462),31:(4.383,.238,3.904),32:(4.280,.379,3.813),
33:(4.244,.309,3.781),34:(4.205,.291,3.746),35:(4.189,.251,3.732),36:(4.176,.220,3.720),
}
src = open(f"{ROOT}/src/periodic_table.c").read()
nbad = 0
for Z,(x1, eps, sigref) in UFF.items():
    m = re.search(r"\{\s*%d\s*,\s*\"(\w+)\",\s*\"[^\"]+\",\s*([0-9.]+),\s*-?[0-9.]+,\s*([0-9.]+),\s*([0-9.]+),\s*([0-9.]+),\s*[0-9.]+,\s*-?[0-9.]+,\s*\d,\s*\n\s*\{\{\{0\}\}, 0, 0\},\s*([0-9.]+)\*KCAL_MOL_TO_EV,\s*([0-9.]+)\s*\}" % Z, src)
    if not m:
        FAIL.append(f"UFF Z={Z} row unparseable"); nbad += 1; continue
    sym, mass, arad, cov, vdw = m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4)), float(m.group(5))
    eps_got, sig_got = float(m.group(6)), float(m.group(7))
    e_ok = abs(eps_got - eps) < 1e-12
    s_ok = abs(sig_got - x1/TW) < 1e-4 and abs(sig_got - sigref) < 5e-3
    if not (e_ok and s_ok):
        nbad += 1
        FAIL.append(f"UFF Z={Z}({sym}): eps {eps_got} vs {eps}; sigma {sig_got} vs x1/2^(1/6)={x1/TW:.5f} published={sigref}")
chk("UFF: 36/36 sigma and epsilon rows", 36-nbad, 36, 0)

print("== AMBER ff99 parm99.dat MOD4 RE vs include/amber_lj.h ==")
p99 = open(f"{ROOT}/audit/external/parm99.dat").read()
mod4 = p99[p99.find("MOD4      RE"):]
parm = {}
for m in re.finditer(r"^\s{2}([A-Za-z0-9*']+)\s+([0-9.]+)\s+([0-9.]+)", mod4, re.M):
    parm[m.group(1)] = (float(m.group(2)), float(m.group(3)))
print(f"  parsed {len(parm)} MOD4 RE atom types")
amber = open(f"{ROOT}/include/amber_lj.h").read()
def amber_val(name):
    m = re.search(r"#define\s+LJ_AMBER_%s_(EPS|SIGMA)\s+([^\n]+)" % name, amber)
    return None
types = ["N","C2","CT","O","OH","OS","P","HN","HC","H1","H4","HA","H5"]
parmkey = {"C2":"C","HN":"H"}
nbad2 = 0
for t in types:
    key = parmkey.get(t,t)
    if key not in parm: FAIL.append(f"parm99 missing {key}"); nbad2 += 1; continue
    r_ref, e_ref = parm[key]
    sigmac = re.search(r"#define\s+LJ_AMBER_%s_SIGMA\s+AMBER_RSTAR_TO_SIGMA\(([0-9.]+)\)" % t, amber)
    epsc   = re.search(r"#define\s+LJ_AMBER_%s_EPS\s+\(([0-9.]+)\s*\*\s*KCAL_MOL_TO_EV\)" % t, amber)
    if not sigmac or not epsc:
        FAIL.append(f"amber_lj.h {t} unparseable"); nbad2 += 1; continue
    rstar, eps = float(sigmac.group(1)), float(epsc.group(1))
    ok = abs(rstar - r_ref) < 1e-9 and abs(eps - e_ref) < 1e-9
    if not ok:
        nbad2 += 1
        FAIL.append(f"AMBER {t}: code Rstar={rstar} eps={eps} vs parm99 Rstar={r_ref} eps={e_ref}")
chk("AMBER ff99: 13 atom types match parm99.dat", 13-nbad2, 13, 0)

print("== derived sigma identities ==")
chk("2*1.8240/2^(1/6)", 2*1.8240/TW, 3.25000, 1e-5)   # code comment says 3.24979 -> doc bug
chk("2*1.6612/2^(1/6)", 2*1.6612/TW, 2.96000, 1e-5)
chk("2*1.9080/2^(1/6)", 2*1.9080/TW, 3.39967, 1e-5)
chk("x1(TI 3.175)/2^(1/6)", 3.175/TW, 2.82860, 1e-5)

print()
if FAIL:
    print(f"A01b: {len(FAIL)} FAILURES"); [print("  -",f) for f in FAIL]; sys.exit(1)
print("A01b: all table checks pass")
