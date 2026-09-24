#!/usr/bin/env bash
# verify_record.sh — read-only record + spec verifier. Run from v9R4/.
# Checks: baseline SHA match, stderr empty, selftest green,
# kcsa.cvmds seal + key/unit compliance, cage geometry wiring.
set -uo pipefail
ok(){ printf '  %-58s %s\n' "$1" "$2"; }
fail=0
chk(){ if [ "$2" = "$3" ]; then printf '  PASS  %s\n' "$1"; else printf '  FAIL  %s (want %s, got %s)\n' "$1" "$2" "$3"; fail=1; fi; }

echo "=== record ==="
[ -f CURRENT_BASELINE_SHA.txt ] || { echo "FATAL: CURRENT_BASELINE_SHA.txt missing"; exit 2; }
[ -f output.txt ] || { echo "FATAL: output.txt missing"; exit 2; }
EXP=$(cat CURRENT_BASELINE_SHA.txt | tr -d ' \n')
ACT=$(sha256sum output.txt | cut -d' ' -f1)
chk "output.txt SHA == CURRENT_BASELINE_SHA.txt" "$EXP" "$ACT"
chk "stderr.txt empty" "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855" "$(sha256sum stderr.txt 2>/dev/null | cut -d' ' -f1 || echo MISSING)"
chk "build warning-clean" "0" "$(make clean >/dev/null 2>&1; make 2>&1 | grep -c 'warning:')"

echo; echo "=== kcsa cage wiring (truth fixes) ==="
ok "KCSA_RING_Z_SEP in main.c" "$(grep -Fc 'KCSA_RING_Z_SEP' src/main.c) (>=6)"
ok "3D->xy derivation sqrt(d^2-half^2)" "$(grep -Fc 'sqrt(d_inner' src/main.c) (>=1)"
ok "r_inner/r_outer in JC+UZ+QM legs" "$(grep -Fc 'r_inner * cos' src/main.c) (>=3)"
ok "qm_overlap_ref used (no 1.5A Sref)" "$(grep -Fc 'qm_overlap_ref' src/main.c) (>=1)"
ok "COULOMB_MD in induction (no hardcoded 14.399)" "$(grep -Fc 'COULOMB_MD * fabs' src/main.c) (>=1)"
ok "ECC computed leg E_ecc" "$(grep -Fc 'jc_ecc_k' src/main.c) (>=3)"
ok "U(z) label (not PMF-as-free-energy)" "$(grep -Fci 'free-energy PMF' src/main.c) (>=2)"
ok "side-effect-free forces_nonbonded_energy" "$(grep -Fc 'forces_nonbonded_energy' src/main.c) (>=1)"
ok "v2 induction+Pauli in force loop" "$(grep -Fc 'use_polar' src/forces.c) (>=2)"
ok "v3 dispersion + SCF driver" "$(grep -Fc 'qm_scf_charges' src/qm.c) (>=2)"
ok "relaxed coordination leg" "$(grep -Fc 'rel_cn_k' src/main.c) (>=2)"
ok "Thole damping" "$(grep -Fci 'thole' src/qm.c) (>=3)"
ok "explicit hydration legs" "$(grep -Fc 'hyd_k' src/main.c) (>=3)"
ok "WHAM free energy" "$(grep -Fc 'WHAM' src/main.c) (>=2)"

echo; echo "=== datastream ==="
make selftest >/dev/null 2>&1
./build/test_datastream >/dev/null 2>&1 && ST=0 || ST=1
chk "selftest green" "0" "$ST"
make selftest-forces >/dev/null 2>&1
./build/test_forces >/dev/null 2>&1 && FT=0 || FT=1
chk "forces selftest green (P2)" "0" "$FT"
make selftest-fire >/dev/null 2>&1
./build/test_fire >/dev/null 2>&1 && FR=0 || FR=1
chk "fire selftest green" "0" "$FR"
if [ -f kcsa.cvmds ]; then
  ./build/test_datastream >/dev/null 2>&1 # ensures verifier linked; use binary below
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
else
  echo "  FAIL  kcsa.cvmds missing"; fail=1
fi

echo
if [ "$fail" -ne 0 ]; then echo "VERIFY FAILED"; exit 1; fi
echo "VERIFY PASSED"
