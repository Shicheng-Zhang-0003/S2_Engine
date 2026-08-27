#!/usr/bin/env bash
# verify_scripts.sh — read-only. Does NOT modify anything. Run from v9R4/.
set -uo pipefail
ok(){ printf '  %-52s %s\n' "$1" "$2"; }

echo "=== s11b: KcsA ring z-separation wiring ==="
ok "KCSA_RING_Z_SEP in main.c        (want >=2)" "$(grep -Fc 'KCSA_RING_Z_SEP' src/main.c)"
ok "half_sep in main.c               (want >=3)" "$(grep -Fc 'half_sep' src/main.c)"
ok "6-arg call with &min_oo          (want 1)"   "$(grep -Fc '&min_oo' src/main.c)"
ok "'sourced z-separation 3.084' text(want 1)"  "$(grep -Fc 'sourced z-separation 3.084 A now' src/main.c)"

echo; echo "=== s13: readme SHA + Demo-12 bullet sync ==="
ok "CURRENT_BASELINE_SHA.txt exists?" "$( [ -f CURRENT_BASELINE_SHA.txt ] && echo YES || echo 'NO -> s13 gate cannot run' )"
ok "old SHA 875a2c0c… in readme      (want 0)"  "$(grep -Fc '875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935' readme.md)"
ok "'antiprism block is now a real result' (>=1)" "$(grep -Fc 'antiprism block is now a real result' readme.md)"

echo; echo "=== s14: three readme antiprism refs ==="
ok "fix1 'Until then…cited' gone     (want 0)"  "$(grep -Fc 'Until then its numbers must not be cited' readme.md)"
ok "fix2 '*(done)*' present          (want >=1)" "$(grep -Fc "Complete Demo 12's antiprism block** *(done)*" readme.md)"
ok "fix3 'symmetry-validated in s10/s10b' (>=1)" "$(grep -Fc 'symmetry-validated in s10/s10b' readme.md)"

echo; echo "=== s15: embedded-block intro paragraph ==="
ok "old 'labels…placeholder' gone    (want 0)"  "$(grep -Fc 'labels its antiprism block a placeholder' readme.md)"
ok "new 'block now a real result'    (want >=1)" "$(grep -Fc 'its antiprism block now a real result' readme.md)"

echo; echo "=== s16: demo_dna_duplex printf \\n escapes ==="
ok "demo_dna_duplex present          (want 1)"  "$(grep -Fc 'static void demo_dna_duplex' src/main.c)"

echo; echo "=== s17: proper WC base-pair geometry ==="
ok "g_N1_to_O6_norm placement math   (want >=1)" "$(grep -Fc 'g_N1_to_O6_norm' src/main.c)"
ok "thymine at vec3(2.9, 0.0, 3.4)   (want 1)"  "$(grep -Fc 'vec3(2.9, 0.0, 3.4)' src/main.c)"

echo; echo "=== s18: string.h + unused-var casts ==="
ok "#include <string.h>              (want 1)"  "$(grep -Fc '#include <string.h>' src/main.c)"
ok "(void)a_to_t;                    (want 1)"  "$(grep -Fc '(void)a_to_t;' src/main.c)"
ok "'c_C4 used implicitly' comment   (want 1)"  "$(grep -Fc 'c_C4 used implicitly' src/main.c)"

echo; echo "=== s19: DNA placement (expected STALE / no-op) ==="
ok "old C 'temp, far' placement      (want 0)"  "$(grep -Fc 'vec3(15.0, 0.0, 0.0)); /* temp, far */' src/main.c)"
ok "old T vec3(15.0, 0.0, 3.4)       (want 0)"  "$(grep -Fc 'vec3(15.0, 0.0, 3.4)' src/main.c)"
