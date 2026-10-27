#!/usr/bin/env bash
# s01_verify_record.sh
#
# Record-verification harness - the foundation every later script
# delegates to. Rebuilds the tree in two configurations (the makefile's
# own normal CFLAGS, and an AddressSanitizer+UBSan build), runs the
# binary in each, and asserts that the fresh output reproduces the
# recorded record byte-for-byte and that both program stderrs are
# empty.
#
# This script does NOT modify any source. It only rebuilds, runs, and
# compares. Safe to run at any time.

set -euo pipefail

# This script lives INSIDE the tree it verifies (moved in from
# biological/ during the v9R4 audit), so HERE is the tree itself. That is
# what makes it safe to run from anywhere: it can no longer resolve to a
# sibling directory the way it did while it sat one level up.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TREE="$HERE"
SCRATCH="$HERE/audit/scratch/s01"
# Record SHA is the single source of truth in the tree's
# CURRENT_BASELINE_SHA.txt.
# Read live so this harness cannot fork from the tree it verifies.
EXPECTED_SHA="$(cat "$TREE/CURRENT_BASELINE_SHA.txt" | tr -d ' \n')"
EMPTY_SHA="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
# AUDIT FIX M7: this script used to live one directory up, at
# biological/s01_verify_record.sh, with TREE="$HERE/v9R4" pointing into
# the sibling tree. Running a copy from anywhere therefore still operated
# on the real tree, which is how an audit pass once rebuilt it
# unannounced. It now lives in the tree it verifies.
#
# AUDIT FIX M8: output.asan.txt used to be tracked and this harness asserted
# the fresh ASan stdout against it. That made a *derived byproduct* a second
# source of truth for the record: every legitimate record change required a
# second artifact to be refreshed in lockstep, and a forgotten refresh failed
# the gate while telling you nothing about the engine. It was also the only
# check that could fail purely because of git hygiene.
#
# The ASan build is still verified here, against the single source of truth:
# it must reproduce the recorded SHA byte-for-byte and leave stderr empty.
# Byte-equality with the normal build is then implied, and is asserted
# explicitly below because it is the property the archive was standing in for.
# Full-audit O5: dead ARCHIVED var removed (M8 leftover, never used).

[ -d "$TREE" ] || { echo "FATAL: $TREE not found"; exit 2; }
command -v sha256sum >/dev/null || { echo "FATAL: sha256sum not found"; exit 2; }

rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"

# FULL-AUDIT: ds_open() stamps wallclock unless SOURCE_DATE_EPOCH is set.
# Leave it live and the normal and sanitizer sidecars can never be compared
# meaningfully; the record is supposed to be reproducible. Default to the
# LAST SOURCE COMMIT's time (not HEAD's): HEAD moves when the record is
# committed, while the source commit is stable.
if [ -z "${SOURCE_DATE_EPOCH:-}" ]; then
  SDE="$(git -C "$TREE" log -1 --format=%ct -- :/src :/include :/tests :/makefile 2>/dev/null || true)"
  [ -n "$SDE" ] || SDE="$(date +%s)"
  SOURCE_DATE_EPOCH="$SDE"
  export SOURCE_DATE_EPOCH
fi

# FULL-AUDIT: the sanitizer build was missing several things the record
# build has, and each omission weakened the gate:
#   -ffp-contract=off         the record build pins FP contraction; the
#                              sanitizer build did not, so it was the ONE
#                              build in the tree that could fuse FMAs
#                              and change the low-order digits.
#   -fno-sanitize-recover=undefined
#                              UBSan diagnostics are RECOVERABLE by
#                              default: the process continues after
#                              printing one. The gate only caught them
#                              because it also asserts stderr is empty -
#                              a happy accident, not a design.
#   -fno-omit-frame-pointer   a report that does fire then has a usable
#                              stack trace.
#   -Werror                    warnings are gated in the record build
#                              but not in the sanitizer build.
# The ASAN_OPTIONS / UBSAN_OPTIONS below are EXPORTED rather than
# inherited: an ambient ASAN_OPTIONS=log_path=/tmp/x would move the
# report out of stderr and the "stderr is empty" check would then pass
# with the error sitting in a file.
ASAN_FLAGS="-g -O0 -fno-omit-frame-pointer -Wall -Wextra -Werror -std=c11 -Iinclude -ffp-contract=off -fsanitize=address,undefined -fno-sanitize-recover=undefined"
export ASAN_OPTIONS="detect_leaks=1:halt_on_error=1:abort_on_error=1:detect_stack_use_after_return=1"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1"

echo "[1/6] clean + normal build (makefile's own CFLAGS)"
make -C "$TREE" clean >/dev/null
make -C "$TREE" >"$SCRATCH/build_normal.log" 2>&1

echo "[2/6] run normal build"
# FULL-AUDIT: the engine writes kcsa.cvmds relative to the CURRENT
# directory, so running it from the caller's directory overwrote the
# TRACKED kcsa.cvmds with a wallclock-stamped one. The runs happen in
# the scratch directory; the tracked sidecar is restored afterwards.
( cd "$SCRATCH" && "$TREE/carbonsim" ) >"$SCRATCH/output.txt" 2>"$SCRATCH/stderr_normal.txt"
# Preserve the normal build's sidecar before the sanitizer run overwrites it.
cp "$SCRATCH/kcsa.cvmds" "$SCRATCH/kcsa.normal.cvmds"

echo "[3/6] clean + ASan build"
make -C "$TREE" clean >/dev/null
# FULL-AUDIT: overriding CFLAGS on the command line replaced the WHOLE
# variable, so the datastream build-identity defines (-DDS_COMPILER_ID /
# -DDS_CFLAGS_ID) were dropped and the sealed kcsa.cvmds produced by the
# sanitizer build recorded "build-compiler: unknown". The identity is
# now passed explicitly, so the sanitizer artifact carries usable
# provenance and the two builds can be compared.
make -C "$TREE" CFLAGS="$ASAN_FLAGS" \
     DS_COMPILER_ID="gcc-asan $(gcc -dumpversion)" \
     DS_CFLAGS_ID="$ASAN_FLAGS" \
     >"$SCRATCH/build_asan.log" 2>&1

echo "[4/6] run ASan build"
( cd "$SCRATCH" && "$TREE/carbonsim" ) >"$SCRATCH/output.asan.txt" 2>"$SCRATCH/stderr_asan.txt"

echo "[5/6] restore tracked build (tree dirtied by clean)"
# Full-audit O6: the final `make clean` deleted TRACKED build/*.o + carbonsim,
# leaving the tree git-dirty by design. Rebuild normal so the working tree is
# left in its committed state (clean build, zero warnings) rather than dirty.
# FULL-AUDIT s01-hardening: clean-first, not plain `make`. A plain make trusts
# ambient .o files, which are ASan-tainted at this point; flag changes are
# invisible to mtime comparison, so the relink fails with __asan undefined
# references (observed 2026-10-08: restore linked ASan .o without -fsanitize).
make -C "$TREE" clean >/dev/null
make -C "$TREE" >"$SCRATCH/build_restore.log" 2>&1 || { echo "FATAL: restore build failed"; exit 1; }

echo "[6/6] assertions"
norm_sha=$(sha256sum "$SCRATCH/output.txt"        | cut -d' ' -f1)
asan_sha=$(sha256sum "$SCRATCH/output.asan.txt"   | cut -d' ' -f1)
norm_err=$(sha256sum "$SCRATCH/stderr_normal.txt" | cut -d' ' -f1)
asan_err=$(sha256sum "$SCRATCH/stderr_asan.txt"   | cut -d' ' -f1)

fail=0
check() {  # $1 label, $2 expected, $3 actual
  if [ "$2" = "$3" ]; then
    printf '  PASS  %s\n' "$1"
  else
    printf '  FAIL  %s\n        expected: %s\n        actual:   %s\n' "$1" "$2" "$3"
    fail=1
  fi
}
check "normal output SHA == recorded"              "$EXPECTED_SHA" "$norm_sha"
check "ASan output SHA == recorded"                "$EXPECTED_SHA" "$asan_sha"
check "ASan output == normal output (byte-identical)" "$norm_sha"  "$asan_sha"
check "ASan stderr empty"                          "$EMPTY_SHA"     "$asan_err"
check "normal stderr empty (no fallback warnings)" "$EMPTY_SHA"     "$norm_err"

# FULL-AUDIT: the sealed sidecar is half the record. Only its seal was ever
# checked here; its CLAIMS could have drifted away from the freshly computed
# values while stdout happened to match. Compare the [claims] sections of all
# three copies: the normal build, the sanitizer build, and the tracked file.
claims_section() { sed -n '/^\[claims\]/,/^\[end\]/p' "$1" 2>/dev/null | grep -v '^\['; }
claims_section "$SCRATCH/kcsa.normal.cvmds" > "$SCRATCH/claims.normal"
claims_section "$SCRATCH/kcsa.cvmds"        > "$SCRATCH/claims.asan"
claims_section "$TREE/kcsa.cvmds"           > "$SCRATCH/claims.tracked"
norm_claims_sha=$(sha256sum "$SCRATCH/claims.normal"  | cut -d' ' -f1)
asan_claims_sha=$(sha256sum "$SCRATCH/claims.asan"    | cut -d' ' -f1)
track_claims_sha=$(sha256sum "$SCRATCH/claims.tracked" | cut -d' ' -f1)
norm_claims_n=$(grep -vc '^#\|^$' "$SCRATCH/claims.normal" || true)
check "fresh sidecar has claim rows" "yes" \
      "$([ "${norm_claims_n:-0}" -gt 0 ] && echo yes || echo no)"
check "ASan sidecar claims == normal sidecar claims" "$norm_claims_sha" "$asan_claims_sha"
check "tracked sidecar claims == fresh normal claims" "$norm_claims_sha" "$track_claims_sha"

if [ -x "$TREE/build/verify_cvmds" ]; then
  "$TREE/build/verify_cvmds" "$SCRATCH/kcsa.normal.cvmds" >/dev/null 2>&1 && NS=0 || NS=1
  "$TREE/build/verify_cvmds" "$SCRATCH/kcsa.cvmds"        >/dev/null 2>&1 && AS=0 || AS=1
  "$TREE/build/verify_cvmds" "$TREE/kcsa.cvmds"           >/dev/null 2>&1 && TS=0 || TS=1
  check "fresh normal sidecar seal verifies" "0" "$NS"
  check "fresh ASan sidecar seal verifies"   "0" "$AS"
  check "tracked sidecar seal verifies"      "0" "$TS"
else
  check "build/verify_cvmds present" "yes" "no"
fi

if [ "$fail" -ne 0 ]; then
  echo
  echo "Record NOT intact. Inspect $SCRATCH for details."
  exit 1
fi
echo
echo "Record intact. Fresh evidence kept in $SCRATCH/"
