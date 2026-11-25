#!/bin/sh
# Build and run the engine-path KcsA probe.
#
# This does NOT modify the sealed engine. It copies src/qm_fluctuation.c
# into a temp directory, makes the timestep / friction / dielectric
# compile-time overrides, and links the copy against the engine's own
# objects. The temp directory is removed on exit.
#
# Usage:
#   side_fixes/probes/run_kcsa_engine_probe.sh          # 200 fs matrix
#   side_fixes/probes/run_kcsa_engine_probe.sh --long   # + 300/500/900 fs
#
# The numbers are diagnostic only. The engine's UNRESOLVED KcsA verdict
# stands; see side_fixes/KCSA_ENGINE_PLUGIN.md.

set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
cd "$ROOT"

if [ ! -f build/sim.o ]; then
    make -s
fi

TMP=$(mktemp -d "${TMPDIR:-/tmp}/kcsa-engine-probe.XXXXXX")
cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT INT TERM

python3 - "$ROOT/src/qm_fluctuation.c" "$TMP/qm_fluctuation_plug.c" <<'PY'
import sys

src, dst = sys.argv[1], sys.argv[2]
text = open(src, encoding="utf-8").read()

subs = [
    ("#define QM_KCSA_GAMMA (1.0 / 1000.0)",
     "#ifndef QM_KCSA_GAMMA\n"
     "#define QM_KCSA_GAMMA (1.0 / 1000.0)\n"
     "#endif"),
    ("#define QM_KCSA_DT 1.0",
     "#ifndef QM_KCSA_DT\n"
     "#define QM_KCSA_DT 1.0\n"
     "#endif\n\n"
     "/* TEMPORARY PROBE: engine default 1.0; side fix uses 4.0. */\n"
     "#ifndef QM_KCSA_EPS\n"
     "#define QM_KCSA_EPS 1.0\n"
     "#endif"),
    ("    sim->use_bonds = 1; sim->use_angles = 1; sim->use_dihedrals = 1;",
     "    sim->use_bonds = 1; sim->use_angles = 1; sim->use_dihedrals = 1;\n"
     "    sim->dielectric = QM_KCSA_EPS;   /* TEMPORARY PROBE */"),
]

for old, new in subs:
    if old not in text:
        sys.exit("patch target not found: " + old)
    text = text.replace(old, new, 1)

open(dst, "w", encoding="utf-8").write(text)
PY

grep -q 'sim->dielectric = QM_KCSA_EPS;' "$TMP/qm_fluctuation_plug.c" || {
    echo "probe: patch did not apply" >&2
    exit 1
}

OBJS=$(find build -maxdepth 1 -name '*.o' \
    ! -name 'main.o' ! -name 'qm_fluctuation.o' \
    ! -name 'tui.o' ! -name 'tui_view.o' ! -name 'tui_screen.o' \
    ! -name 'test_*.o' | sort)
if [ -z "$OBJS" ]; then
    echo "probe: no engine objects in build/" >&2
    exit 1
fi

build() {
    name=$1
    shift
    # shellcheck disable=SC2086
    cc -O2 -g -Wall -Wextra -Werror -ffp-contract=off -std=c11 \
       -Iinclude -I. "$@" \
       "$HERE/kcsa_engine_probe.c" "$TMP/qm_fluctuation_plug.c" \
       $OBJS -lm -o "$TMP/$name"
}

run() {
    label=$1
    bin=$2
    tether=$3
    ne=$4
    ns=$5
    printf '\n=== %s (tether=%s eV/A^2, %s equil + %s sample steps) ===\n' \
        "$label" "$tether" "$ne" "$ns"
    "$TMP/$bin" "$tether" "$ne" "$ns"
}

build baseline
build eps4    -DQM_KCSA_EPS=4.0
build dt05    -DQM_KCSA_DT=0.5
build gam02   -DQM_KCSA_GAMMA=0.02
build eps_dt  -DQM_KCSA_EPS=4.0 -DQM_KCSA_DT=0.5
build eps_gam -DQM_KCSA_EPS=4.0 -DQM_KCSA_GAMMA=0.02
build dt_gam  -DQM_KCSA_DT=0.5 -DQM_KCSA_GAMMA=0.02
build full    -DQM_KCSA_EPS=4.0 -DQM_KCSA_DT=0.5 -DQM_KCSA_GAMMA=0.02

run "engine baseline: eps=1, dt=1.0, gamma=0.001" baseline 8.0 200 200
run "eps=4 only" eps4 8.0 200 200
run "dt=0.5 only" dt05 8.0 200 200
run "gamma=0.02 only" gam02 8.0 200 200
run "eps=4 + dt=0.5" eps_dt 8.0 200 200
run "eps=4 + gamma=0.02" eps_gam 8.0 200 200
run "dt=0.5 + gamma=0.02" dt_gam 8.0 200 200
run "full side-fix set: eps=4, dt=0.5, gamma=0.02" full 8.0 200 200

if [ "${1:-}" = "--long" ]; then
    run "full set, 300 fs" full 8.0 200 400
    run "full set, 500 fs" full 8.0 200 800
    run "full set, 900 fs" full 8.0 200 1600
    run "eps+gamma, 300 fs" eps_gam 8.0 200 400
    run "eps+gamma, 500 fs" eps_gam 8.0 200 800
fi
