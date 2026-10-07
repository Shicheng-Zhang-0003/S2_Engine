#!/usr/bin/env bash
#
# audit/run_audit.sh — build and run everything in audit/.
#
# This directory used to live outside the engine tree, so its scripts carried
# absolute paths and no single command built all of it. Everything is now
# resolved from this script's own location, so it works from any directory.
#
#   ./audit/run_audit.sh             every oracle and harness
#   ./audit/run_audit.sh --revert    ...and link the suite against the pre-fix
#                                    engine, to prove the regression tests
#                                    still go red on the defects they claim
#                                    to catch
#
# Nothing here is part of `make test`; this is the audit's own working
# material, kept so the findings stay reproducible.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TREE="$(cd "$HERE/.." && pwd)"
BUILD="$HERE/build"
REVERT=0
[ "${1:-}" = "--revert" ] && REVERT=1

cd "$TREE"
mkdir -p "$BUILD"

# Engine objects, minus main.o. Rebuilt if needed; `make` is the honest way
# to ask, and a stale object is how this tree once hid a build failure.
# Excludes tui.o/tui_view.o (separate binary with its own main) — same
# filter as makefile ENG_OBJS — otherwise every harness link fails with
# multiple definition of `main`. tui_screen.o is likewise TUI-only (and
# would drag -lncursesw into harness links); harnesses test the engine.
OBJS=()
for o in "$TREE"/build/*.o; do
    b="$(basename "$o")"
    case "$b" in main.o|tui.o|tui_view.o|tui_screen.o|test_*|*.d) continue ;; esac
    OBJS+=("$o")
done
if [ "${#OBJS[@]}" -eq 0 ]; then
    echo "FATAL: no engine objects in $TREE/build — run 'make' first" >&2
    exit 2
fi

pass=0; fail=0
note() { printf '\n=== %s ===\n' "$1"; }
verdict() {
    # $1 = label, $2 = exit status
    if [ "$2" -eq 0 ]; then pass=$((pass+1)); printf '  PASS  %s\n' "$1"
    else fail=$((fail+1)); printf '  FAIL  %s (exit %s)\n' "$1" "$2"; fi
}

# ── Python oracles ─────────────────────────────────────────────────────────
# Two kinds, and conflating them would be a lie:
#
#   GATES      test the CURRENT engine. Must pass.
#   FORENSIC   characterise a defect that has since been fixed. Their failures
#              ARE the output: a05_induction.py contains a pure-Python
#              reimplementation of the PRE-FIX qm_solve_dipoles - it does not
#              link the engine at all - so it still reports a large dipole
#              residual, which is the bug it was written to demonstrate.
#              Running it against a fixed engine and calling the result a
#              failure would be exactly backwards. It is run for its record,
#              not gated on.
FORENSIC="a05_induction.py"
for f in "$HERE"/math/*.py; do
    [ -e "$f" ] || continue
    b="$(basename "$f")"
    case " $FORENSIC " in *" $b "*) forensic=1; label=" (forensic, not gated)" ;; *)
                              forensic=0; label="" ;; esac
    note "math/$b$label"
    rc=0
    python3 "$f" >"$BUILD/$b.log" 2>&1 || rc=$?
    if [ "$forensic" -eq 1 ]; then
        nfail=$(grep -c '^FAIL' "$BUILD/$b.log")
        printf '  INFO  %s reported %s red checks against the PRE-FIX solver (expected)\n' "$b" "$nfail"
    else
        verdict "$b" "$rc"
        [ "$rc" -eq 0 ] || sed -n '1,15p' "$BUILD/$b.log"
    fi
done

# ── External oracles ───────────────────────────────────────────────────────
for f in "$HERE"/external/*.py; do
    [ -e "$f" ] || continue
    note "external/$(basename "$f")"
    rc=0
    python3 "$f" >"$BUILD/$(basename "$f").log" 2>&1 || rc=$?
    verdict "$(basename "$f")" "$rc"
    [ "$rc" -eq 0 ] || sed -n '1,15p' "$BUILD/$(basename "$f").log"
done

# ── C harnesses ────────────────────────────────────────────────────────────
# Each harness includes ../../include/*.h relative to its own directory,
# which is why it sits two levels under the tree root.
note "harness/ (linked against the engine objects)"
for c in "$HERE"/harness/*.c; do
    [ -e "$c" ] || continue
    b="$(basename "$c" .c)"
    if ! gcc -O2 -std=c11 -I"$TREE/include" -o "$BUILD/$b" "$c" "${OBJS[@]}" -lm \
         >"$BUILD/$b.build.log" 2>&1; then
        verdict "$b (build)" 1
        sed -n '1,8p' "$BUILD/$b.build.log"
        continue
    fi
    rc=0
    "$BUILD/$b" >"$BUILD/$b.log" 2>&1 || rc=$?
    verdict "$b" "$rc"
    [ "$rc" -eq 0 ] || sed -n '1,20p' "$BUILD/$b.log"
done

# ── The revert test ────────────────────────────────────────────────────────
# orig/*.c are byte-identical to their commits and keep their original
# ../../-free `#include "../include/..."` lines, which resolved from src/.
# Build them in a scratch directory at the same depth so those includes still
# mean what they meant, without editing the evidence.
if [ "$REVERT" -eq 1 ]; then
    note "revert: the new suite against the PRE-FIX engine"
    STAGE="$BUILD/revert"
    rm -rf "$STAGE"; mkdir -p "$STAGE"
    cp "$TREE"/*.h "$TREE/include"/*.h "$STAGE/" 2>/dev/null
    mkdir -p "$STAGE/include" "$STAGE/src"
    cp "$TREE/include"/*.h "$STAGE/include/"
    cp "$HERE/orig/qm_orig.c" "$STAGE/src/qm.c"
    cp "$HERE/orig/integrator_orig.c" "$STAGE/src/integrator.c"
    cp "$HERE/orig/constants_orig.h" "$STAGE/include/constants.h"

    ROUT=()
    for o in "${OBJS[@]}"; do
        b="$(basename "$o")"
        case "$b" in qm.o|integrator.o) continue ;; esac
        ROUT+=("$o")
    done
    if gcc -O2 -std=c11 -I"$STAGE/include" -c -o "$STAGE/src/qm.o" "$STAGE/src/qm.c" \
          >"$STAGE/qm.build.log" 2>&1 &&
       gcc -O2 -std=c11 -I"$STAGE/include" -c -o "$STAGE/src/integrator.o" "$STAGE/src/integrator.c" \
          >>"$STAGE/qm.build.log" 2>&1 &&
       gcc -O2 -std=c11 -I"$TREE/include" -o "$STAGE/test_regression_orig" \
          "$TREE/tests/test_regression.c" "$STAGE/src/qm.o" "$STAGE/src/integrator.o" \
          "${ROUT[@]}" -lm >>"$STAGE/qm.build.log" 2>&1; then
        "$STAGE/test_regression_orig" >"$STAGE/regression.log" 2>&1
        rc=$?
        np=$(grep -c '^  PASS' "$STAGE/regression.log")
        nf=$(grep -c '^  FAIL' "$STAGE/regression.log")
        printf '  suite: %s passed, %s failed against the PRE-FIX engine\n' "$np" "$nf"
        grep '^  FAIL' "$STAGE/regression.log" | sed 's/^/    /'
        # A revert run that fails everything is broken; one that passes
        # everything means the suite no longer tests the defects.
        if [ "$rc" -ne 0 ] && [ "$nf" -ge 1 ]; then
            printf '  PASS  the pre-fix engine is still caught by the suite (%s red checks)\n' "$nf"
            pass=$((pass+1))
        else
            printf '  FAIL  the suite did not catch the pre-fix engine\n'
            fail=$((fail+1))
        fi
    else
        verdict "revert build" 1
        sed -n '1,12p' "$STAGE/qm.build.log"
    fi
fi

note "summary"
printf '  %s passed, %s failed   (logs in %s)\n' "$pass" "$fail" "$BUILD"
[ "$fail" -eq 0 ] || exit 1
exit 0