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
# Record SHA is the single source of truth in v9R4/CURRENT_BASELINE_SHA.txt.
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
ARCHIVED="$TREE/output.asan.txt"

[ -d "$TREE" ] || { echo "FATAL: $TREE not found"; exit 2; }
command -v sha256sum >/dev/null || { echo "FATAL: sha256sum not found"; exit 2; }

rm -rf "$SCRATCH"; mkdir -p "$SCRATCH"

ASAN_FLAGS="-g -O0 -Wall -Wextra -std=c11 -Iinclude -fsanitize=address,undefined"

echo "[1/6] clean + normal build (makefile's own CFLAGS)"
make -C "$TREE" clean >/dev/null
make -C "$TREE" >"$SCRATCH/build_normal.log" 2>&1

echo "[2/6] run normal build"
"$TREE/carbonsim" >"$SCRATCH/output.txt" 2>"$SCRATCH/stderr_normal.txt"

echo "[3/6] clean + ASan build"
make -C "$TREE" clean >/dev/null
make -C "$TREE" CFLAGS="$ASAN_FLAGS" >"$SCRATCH/build_asan.log" 2>&1

echo "[4/6] run ASan build"
"$TREE/carbonsim" >"$SCRATCH/output.asan.txt" 2>"$SCRATCH/stderr_asan.txt"

echo "[5/6] leave tree source-clean"
make -C "$TREE" clean >/dev/null

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

if [ "$fail" -ne 0 ]; then
  echo
  echo "Record NOT intact. Inspect $SCRATCH for details."
  exit 1
fi
echo
echo "Record intact. Fresh evidence kept in $SCRATCH/"
