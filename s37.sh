#!/usr/bin/env bash
set -euo pipefail
# s37: KcsA §1 — add the dehydration penalty (the missing thermodynamic leg).
#
# The vacuum model computes only the filter-binding leg. Real selectivity
# also depends on the cost of stripping the ion's hydration shell, which a
# vacuum model structurally cannot represent. This adds that leg as a
# constant correction using literature hydration free energies.
#
# Cycle:  Ion(aq) --dG_dehyd--> Ion(gas) --dG_filter--> Ion(filter)
#   dG_binding(ion) = dG_dehyd(ion) + dG_filter(ion)
#   ddG_select = [dG_filter(K)-dG_filter(Na)] + [dG_dehyd(K)-dG_dehyd(Na)]
#                 └─ vacuum computes this ─┘    └─ THIS SCRIPT ADDS THIS ─┘
#
# Source: Marcus, Y. "Thermodynamics of solvation of ions. Part 5."
#   J. Chem. Soc. Faraday Trans. 87, 2995 (1991).  1 eV = 96.485 kJ/mol.
#   K+:  dG_hyd=-295 kJ/mol=-3.057 eV  ->  dG_dehyd=+3.057 eV
#   Na+: dG_hyd=-365 kJ/mol=-3.783 eV  ->  dG_dehyd=+3.783 eV
#   Na+ pays 0.726 eV MORE to dehydrate — this is the selectivity term.
#
# RECORD-MOVING: Demo 12 gains a corrected-selectivity section.
echo "=== s37: KcsA dehydration penalty correction ==="

python3 - <<'PYEOF'
import sys

with open('src/main.c') as f:
    src = f.read()

# ── 1. Constants, inserted after KCSA_RING_Z_SEP ───────────────────────
anchor = '#define KCSA_RING_Z_SEP 3.084'
if anchor not in src:
    print('FAIL: KCSA_RING_Z_SEP anchor not found'); sys.exit(1)

constants = anchor + '''

/* ══ KcsA dehydration penalty (s37) ═══════════════════════════════════
* Hydration free energies -> dehydration cost. The thermodynamic leg the
* vacuum model cannot represent. Source: Marcus, Y. "Thermodynamics of
* solvation of ions. Part 5." J. Chem. Soc. Faraday Trans. 87, 2995
* (1991). Conversion: 1 eV = 96.485 kJ/mol.
*   K+:  dG_hyd = -295 kJ/mol = -3.057 eV  ->  dG_dehyd = +3.057 eV
*   Na+: dG_hyd = -365 kJ/mol = -3.783 eV  ->  dG_dehyd = +3.783 eV
* Na+ pays 0.726 eV MORE to dehydrate - this is the selectivity term. */
#define KCSA_DEHYD_K_EV   3.057
#define KCSA_DEHYD_NA_EV  3.783
/* Experimental K+/Na+ selectivity ~1000:1 for KcsA. At 300 K the free
* energy is -kT*ln(1000) = -0.179 eV (K+ favored). Used to validate the
* corrected magnitude. k_B from CODATA 2018. */
#define KCSA_KB_EV        8.617333262e-5
#define KCSA_T_KELVIN     300.0
#define KCSA_EXPT_RATIO   1000.0'''

src = src.replace(anchor, constants, 1)
print('  [OK] Added dehydration constants.')

# ── 2. Correction block at end of demo_kcsa_filter() ──────────────────
fn = src.find('static void demo_kcsa_filter(void) {')
if fn == -1:
    print('FAIL: demo_kcsa_filter not found'); sys.exit(1)
brace = src.find('{', fn)
depth, pos = 0, brace
while pos < len(src):
    if   src[pos] == '{': depth += 1
    elif src[pos] == '}':
        depth -= 1
        if depth == 0: break
    pos += 1
# pos = closing brace of demo_kcsa_filter

block = '''
    /* ══ DEHYDRATION-CORRECTED SELECTIVITY (s37) ══════════════════════
    * The vacuum tests above compute only the filter-binding leg. The
    * dehydration cost is a property of the ION (site-independent), so
    * the same correction applies to every site. Applied here to the
    * antiprism result, the most complete cage model. */
    {
        double k_filt  = kcsa_antiprism_energy(19, 1.0, "K+ ", 2.70, 2.83, NULL);
        double na_filt = kcsa_antiprism_energy(11, 1.0, "Na+", 2.70, 2.83, NULL);
        double vac_ddG  = k_filt - na_filt;                       /* + = Na+ favored */
        double dehyd    = KCSA_DEHYD_K_EV - KCSA_DEHYD_NA_EV;     /* -0.726 eV */
        double corr_ddG = vac_ddG + dehyd;                        /* - = K+ favored */
        double expt_ddG = -KCSA_KB_EV * KCSA_T_KELVIN * log(KCSA_EXPT_RATIO);

        printf("\\n");
        printf("-- Dehydration-corrected selectivity (s37) --\\n");
        printf("Filter binding (antiprism): K+ = %.6f eV  Na+ = %.6f eV\\n",
               k_filt, na_filt);
        printf("Vacuum selectivity dG(K)-dG(Na) = %+.4f eV (%s)\\n",
               vac_ddG, vac_ddG > 0 ? "Na+ favored, wrong direction"
                                    : "K+ favored");
        printf("Dehydration penalty: K+ = +%.3f eV  Na+ = +%.3f eV (Marcus 1991)\\n",
               KCSA_DEHYD_K_EV, KCSA_DEHYD_NA_EV);
        printf("Corrected selectivity = %+.4f eV (%s)\\n",
               corr_ddG, corr_ddG < 0 ? "K+ favored, CORRECT direction"
                                      : "Na+ favored, still wrong");
        printf("Experimental (1000:1 at 300 K) = %+.4f eV\\n", expt_ddG);
        printf("Deviation from experiment: %.4f eV\\n", corr_ddG - expt_ddG);
    }
'''

src = src[:pos] + block + src[pos:]
print('  [OK] Added dehydration correction block.')

with open('src/main.c', 'w') as f:
    f.write(src)
PYEOF

echo "[2/3] Building..."
make clean >/dev/null 2>&1
make > /tmp/s37_build.log 2>&1 || { echo "BUILD FAILED:"; cat /tmp/s37_build.log; exit 1; }
echo "  Build clean."

echo "[3/3] Corrected output:"
./carbonsim 2>/dev/null | grep -A 7 "Dehydration-corrected"

echo "=== s37 complete (record moved; re-anchor via s40) ==="
