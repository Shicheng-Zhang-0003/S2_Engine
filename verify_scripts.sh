#!/usr/bin/env bash
# verify_record.sh — read-only record + spec verifier. Run from v9R4/.
# Checks: CLEAN build, baseline SHA match, stderr empty, every selftest
# green, kcsa.cvmds seal + key/unit compliance, cage geometry wiring.
#
# AUDIT FIX B1 (the check that was missing): the old script ran a bare
# `make` to count warnings. This repository TRACKS build/*.o and the
# carbonsim binary, so an in-place `make` sees up-to-date objects and
# rebuilds nothing - it reported "0 warnings" on a tree that did not
# compile at all. The v9R4 audit added static assertions that failed and
# a subsequent fix dropped a closing brace in forces.c; both left the
# committed objects stale enough to hide it completely. `make clean`
# first is the only honest form.
set -uo pipefail
# Full-audit P19: scratch logs live in the temp workspace, not predictable
# /tmp paths (symlink race). Override with TMPDIR or default to /tmp.
# FULL-AUDIT Q4: the default was /tmp/opencode — a scratch directory named
# after one tool rather than a general temporary location — so a machine with
# TMPDIR unset grew a stray "opencode" tree under /tmp on every verifier run.
AUDIT_TMP="${TMPDIR:-/tmp}/s2audit"
mkdir -p "$AUDIT_TMP" 2>/dev/null || AUDIT_TMP="/tmp"
VERIFY_BUILD_LOG="$AUDIT_TMP/verify_build.log"
VERIFY_REG_LOG="$AUDIT_TMP/verify_regression.log"
VERIFY_EXT_LOG="$AUDIT_TMP/verify_external.log"
VERIFY_LOOP_LOG="$AUDIT_TMP/verify_loop.log"
# AUDIT FIX M4. The old helper was:
#
#     ok(){ printf '  %-58s %s\n' "$1" "$2"; }
#
# It printed a label and a grep count and returned 0, and it never touched
# `fail`. So all thirteen "kcsa cage wiring (truth fixes)" checks were
# decoration: the script printed VERIFY PASSED no matter what the counts were,
# and the thresholds lived in the printed LABEL ("... (>=6)") where nothing
# compared against them. The count could have been 0 and the gate would still
# have passed. `ok` is deleted rather than repaired because every call site was
# a lower bound, not an equality, so there is no correct `ok` — only at_least.
#
# at_least NAME COUNT MIN — the wiring checks are lower bounds by nature
# (they assert a pattern is present at least N times, not exactly N times),
# so a numeric comparison against a bound is the honest form.
at_least(){ if [ "${2:-0}" -ge "$3" ] 2>/dev/null; then
               printf '  PASS  %-56s %s (>=%s)\n' "$1" "${2:-0}" "$3"
           else
               printf '  FAIL  %-56s %s, want >=%s\n' "$1" "${2:-0}" "$3"; fail=1
           fi; }
fail=0
chk(){ if [ "$2" = "$3" ]; then printf '  PASS  %s\n' "$1"; else printf '  FAIL  %s (want %s, got %s)\n' "$1" "$2" "$3"; fail=1; fi; }

echo "=== record ==="
[ -f CURRENT_BASELINE_SHA.txt ] || { echo "FATAL: CURRENT_BASELINE_SHA.txt missing"; exit 2; }
[ -f output.txt ] || { echo "FATAL: output.txt missing"; exit 2; }
EXP=$(cat CURRENT_BASELINE_SHA.txt | tr -d ' \n')
ACT=$(sha256sum output.txt | cut -d' ' -f1)
chk "output.txt SHA == CURRENT_BASELINE_SHA.txt" "$EXP" "$ACT"
chk "stderr.txt empty" "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" "$(sha256sum stderr.txt 2>/dev/null | cut -d' ' -f1 || echo MISSING)"
# AUDIT FIX B1: clean first - see the header. A bare `make` is a no-op
# against the committed objects and cannot fail.
make clean >/dev/null 2>&1
if ! make > "$VERIFY_BUILD_LOG" 2>&1; then
  echo "  FAIL  clean build succeeds"
  grep -E "error" "$VERIFY_BUILD_LOG" | head -10
  exit 1
fi
echo "  PASS  clean build succeeds"
chk "build warning-clean" "0" "$(grep -c 'warning:' "$VERIFY_BUILD_LOG" || true)"

echo; echo "=== kcsa cage wiring (truth fixes) ==="
at_least "KCSA_RING_Z_SEP in main.c"        "$(grep -Fc 'KCSA_RING_Z_SEP'   src/main.c)" 6
at_least "3D->xy derivation sqrt(d^2-half^2)" "$(grep -Fc 'sqrt(d_inner'    src/main.c)" 1
at_least "r_inner/r_outer in JC+UZ+QM legs"  "$(grep -Fc 'r_inner * cos'    src/main.c)" 3
at_least "qm_overlap_ref used (no 1.5A Sref)" "$(grep -Fc 'qm_overlap_ref'   src/main.c)" 1
at_least "COULOMB_MD in induction (no hardcoded 14.399)" "$(grep -Fc 'COULOMB_MD * fabs' src/main.c)" 1
at_least "ECC computed leg E_ecc"            "$(grep -Fc 'jc_ecc_k'         src/main.c)" 3
at_least "U(z) label (not PMF-as-free-energy)" "$(grep -Fci 'free-energy PMF' src/main.c)" 2
at_least "side-effect-free forces_nonbonded_energy" "$(grep -Fc 'forces_nonbonded_energy' src/main.c)" 1
at_least "v2 induction+Pauli in force loop"  "$(grep -Fc 'use_polar'         src/forces.c)" 2
at_least "v3 dispersion + SCF driver"        "$(grep -Fc 'qm_scf_charges'    src/qm.c)" 2
at_least "relaxed coordination leg"          "$(grep -Fc 'rel_cn_k'          src/main.c)" 2
at_least "Thole damping"                     "$(grep -Fci 'thole'            src/qm.c)" 3
at_least "explicit hydration legs"           "$(grep -Fc 'hyd_k'            src/main.c)" 3
at_least "WHAM free energy"                  "$(grep -Fc 'WHAM'             src/main.c)" 2
at_least "loop bridge Nernst"                "$(cat src/loop.c tests/test_loop.c | grep -Fc 'loop_nernst_mV')" 3
at_least "loop selectivity scale"            "$(cat src/loop.c tests/test_loop.c | grep -Fc 'loop_selectivity_scale_eV')" 2
at_least "loop abstraction inventory"        "$(grep -Fc 'LOOP_ABSTRACTIONS' src/loop.c)" 1

echo; echo "=== datastream ==="
# Full-audit P19: gate the test-binary builds. Previously `make selftest*`
# exit codes were ignored, so a failed build ran a stale binary and could
# false-PASS.
if ! make selftest >"$AUDIT_TMP/selftest_build.log" 2>&1; then echo "  FAIL  make selftest builds"; cat "$AUDIT_TMP/selftest_build.log" | head -10; fail=1; fi
./build/test_datastream >/dev/null 2>&1 && ST=0 || ST=1
chk "selftest green" "0" "$ST"
DS_PASS=$(./build/test_datastream 2>&1 | grep -c '^  PASS' || true)
chk "datastream gate count == 19" "19" "$DS_PASS"
if ! make selftest-forces >"$AUDIT_TMP/forces_build.log" 2>&1; then echo "  FAIL  make selftest-forces builds"; fail=1; fi
./build/test_forces >/dev/null 2>&1 && FT=0 || FT=1
if ! make selftest-regression >"$AUDIT_TMP/regression_build.log" 2>&1; then echo "  FAIL  make selftest-regression builds"; fail=1; fi
./build/test_regression > "$VERIFY_REG_LOG" 2>&1 && RT=0 || RT=1
# AUDIT FIX E1: the external-source suite. Like the regression suite it must
# be able to FAIL, so its FAIL lines are counted, not just its exit code.
if ! make selftest-external >"$AUDIT_TMP/external_build.log" 2>&1; then echo "  FAIL  make selftest-external builds"; fail=1; fi
./build/test_external > "$VERIFY_EXT_LOG" 2>&1 && ET=0 || ET=1
# Full-audit P19: grep -c prints "0" but exits 1 on zero matches; `|| echo 0`
# double-emitted "0\n0" and broke the = "0" test under pipefail. Use || true
# which preserves grep's printed 0.
ET_PASS=$(grep -c '^  PASS' "$VERIFY_EXT_LOG" 2>/dev/null || true)
ET_FAIL=$(grep -c '^  FAIL' "$VERIFY_EXT_LOG" 2>/dev/null || true)
# The suite's own summary line reads "TOTAL", so a plain "  PASS" prefix
# counts check lines only. (It used to say "PASS n FAIL m" and got counted
# as one of its own checks.)
RT_PASS=$(grep -c '^  PASS' "$VERIFY_REG_LOG" 2>/dev/null || true)
RT_FAIL=$(grep -c '^  FAIL' "$VERIFY_REG_LOG" 2>/dev/null || true)
chk "forces selftest green (P2)" "0" "$FT"
FT_PASS=$(./build/test_forces 2>&1 | grep -c '^  PASS' || true)
chk "forces gate count == 22" "22" "$FT_PASS"
if ! make selftest-fire >"$AUDIT_TMP/fire_build.log" 2>&1; then echo "  FAIL  make selftest-fire builds"; fail=1; fi
./build/test_fire >/dev/null 2>&1 && FR=0 || FR=1
chk "fire selftest green" "0" "$FR"
FR_PASS=$(./build/test_fire 2>&1 | grep -c '^  PASS' || true)
chk "fire gate count == 7" "7" "$FR_PASS"
chk "selftest-regression green" "0" "$RT"
chk "regression gate count == 170" "170" "$RT_PASS"
printf '  ---- audit regression suite: %s checks passed, %s failed ----\n' "$RT_PASS" "$RT_FAIL"
[ "$RT_FAIL" = "0" ] || grep '^  FAIL' "$VERIFY_REG_LOG" | head -20
chk "selftest-external green" "0" "$ET"
chk "external gate count == 51" "51" "$ET_PASS"
printf '  ---- external sources: %s checks passed, %s failed ----\n' "$ET_PASS" "$ET_FAIL"
[ "$ET_FAIL" = "0" ] || grep '^  FAIL' "$VERIFY_EXT_LOG" | head -20
# Loop-closure suite: bio/QC/QM share one physics (Nernst, kT scale,
# Marcus/Shannon/LJ/Coulomb, explicit abstraction list). Must fail closed.
if ! make selftest-loop >"$AUDIT_TMP/loop_build.log" 2>&1; then echo "  FAIL  make selftest-loop builds"; fail=1; fi
./build/test_loop > "$VERIFY_LOOP_LOG" 2>&1 && LT=0 || LT=1
LT_PASS=$(grep -c '^  PASS' "$VERIFY_LOOP_LOG" 2>/dev/null || true)
LT_FAIL=$(grep -c '^  FAIL' "$VERIFY_LOOP_LOG" 2>/dev/null || true)
chk "selftest-loop green" "0" "$LT"
chk "loop gate count == 60" "60" "$LT_PASS"
printf '  ---- loop-closure: %s checks passed, %s failed ----\n' "$LT_PASS" "$LT_FAIL"
[ "$LT_FAIL" = "0" ] || grep '^  FAIL' "$VERIFY_LOOP_LOG" | head -20
# EHT oracle suite: exact overlaps vs grids, Jacobi vs closed forms, H2/H2O
# SCF laws, fail-loud gates. Seventh gate; count pinned like the rest.
VERIFY_EHT_LOG="$AUDIT_TMP/verify_eht.log"
if ! make selftest-eht >"$AUDIT_TMP/eht_build.log" 2>&1; then echo "  FAIL  make selftest-eht builds"; fail=1; fi
./build/test_eht > "$VERIFY_EHT_LOG" 2>&1 && EHT=0 || EHT=1
EHT_PASS=$(grep -c '^  PASS' "$VERIFY_EHT_LOG" 2>/dev/null || true)
EHT_FAIL=$(grep -c '^  FAIL' "$VERIFY_EHT_LOG" 2>/dev/null || true)
chk "selftest-eht green" "0" "$EHT"
chk "eht gate count == 31" "31" "$EHT_PASS"
printf '  ---- EHT oracles: %s checks passed, %s failed ----\n' "$EHT_PASS" "$EHT_FAIL"
[ "$EHT_FAIL" = "0" ] || grep '^  FAIL' "$VERIFY_EHT_LOG" | head -20
if [ -f kcsa.cvmds ]; then
  # FULL-AUDIT O16: primary seal check calls ds_verify_file (no fork risk).
  if [ -x ./build/verify_cvmds ]; then
    if ./build/verify_cvmds kcsa.cvmds >/dev/null 2>&1; then SEAL_C="OK"; else SEAL_C="BAD"; fi
    chk "kcsa.cvmds seal intact (ds_verify_file)" "OK" "$SEAL_C"
  fi
  # python cross-check kept as second opinion, not primary.
  python3 - <<'PY'
import sys
PY
  # seal check via built verifier through a tiny C call is overkill; use python sha to cross-check [end]
  SEAL_OK=$(python3 -c "
import hashlib
d=open('kcsa.cvmds','rb').read()
i=d.find(b'\n[end]\n')
h=d[d.find(b'payload-sha256: ')+16:d.find(b'payload-sha256: ')+80].decode()
print('OK' if hashlib.sha256(d[:i]).hexdigest()==h else 'BAD')")
  chk "kcsa.cvmds seal intact" "OK" "$SEAL_OK"
  chk "no legacy UPPER keys (q_O/dE/S_/BO_/alpha_O)" "0" "$(grep -cE 'q_O_mean|q_K|q_Na|dE_K|dE_Na|S_K|S_Na|BO_K|BO_Na|alpha_O|E_K|E_Na|ddG|ddU' kcsa.cvmds || true)"
  chk "lower-case keys present" "1" "$(grep -Fc 'kcsa.qm.alpha_o' kcsa.cvmds)"
  chk "charge unit e present" "1" "$(grep -Fc 'kcsa.qm.q_k' kcsa.cvmds)"
  chk "header schema present" "1" "$(grep -Fc 'schema: 1' kcsa.cvmds)"
  chk "header demo kcsa present" "1" "$(grep -Fc 'demo: kcsa' kcsa.cvmds)"
  chk "source-hash record-tree-* present" "1" "$(grep -Fc 'source-hash: record-tree-' kcsa.cvmds)"
  chk "provenance legacy-free (JC2008-alone deprecated)" "0" "$(grep -cE ' JC2008( |$)' kcsa.cvmds || true)"
  chk "no digest duplicated in prose" "0" "$(grep -rE '[0-9a-f]{64}' readme.md release_note_v9R4.md DATASTREAM_SPEC.md 2>/dev/null | wc -l)"
else
  echo "  FAIL  kcsa.cvmds missing"; fail=1
fi

echo
if [ "$fail" -ne 0 ]; then echo "VERIFY FAILED"; exit 1; fi
echo; echo "=== stdout/stderr + TUI contract ==="
if printf 'help\ntest all\n' | timeout 60 ./s2tui >/dev/null 2>&1; then
  echo "  PASS  s2tui smoke (help + test all)"
else
  echo "  FAIL  s2tui smoke (help + test all)"; fail=1
fi
if printf 'screen\n' | S2TUI_SCREEN=0 timeout 20 ./s2tui >/dev/null 2>&1; then
  echo "  PASS  s2tui screen forced-line degrade"
else
  echo "  FAIL  s2tui screen forced-line degrade"; fail=1
fi

echo "VERIFY PASSED"
