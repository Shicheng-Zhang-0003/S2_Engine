#!/usr/bin/env python3
"""
s41b.py — HISTORICAL script (s41 era), superseded 2026-09-24.

Original purpose: fix the datastream.h comment stray '*/' and gate on
the s41 behavior-neutral record (datastream not yet wired).
DO NOT RUN as a verifier today: Gate-2's "not yet wired" premise is false
(s42 wired kcsa.cvmds in Demo 12), Gate-5's s01 EXPECTED_SHA is now read
live from CURRENT_BASELINE_SHA.txt, and the record has intentionally moved
(truth-fix 4356c970...). Use v9R4/verify_scripts.sh instead.
Kept for provenance only.
"""
import hashlib, pathlib, subprocess, sys

# This file used to live at biological/s41b.py and reach into a sibling
# tree via TREE = ROOT / "v9R4". It now lives inside the tree itself, so
# TREE is just ROOT. Kept historical and unrunnable-by-design; see the
# module docstring.
ROOT = pathlib.Path(__file__).resolve().parent
TREE = ROOT

BROKEN = "ds_add_atom*/step*/claim*"
FIXED  = "ds_add_atom* -> ds_add_step* -> ds_add_claim*"

def die(msg):
    print(f"FATAL: {msg}"); sys.exit(2)

def read(p):  return pathlib.Path(p).read_text()
def write(p, s): pathlib.Path(p).write_text(s)

def run(cmd, **kw):
    return subprocess.run(cmd, cwd=str(TREE), capture_output=True, text=True, **kw)

hdr = TREE / "include" / "datastream.h"
if not hdr.is_file(): die("include/datastream.h missing — run s41 first")
mk = read(TREE / "makefile")
if "selftest:" not in mk: die("makefile has no selftest target — run s41 first")

# ── PRE / DO: patch the comment idempotently ─────────────────────────
src = read(hdr)
if BROKEN in src:
    write(hdr, src.replace(BROKEN, FIXED, 1))
    print(f"patched include/datastream.h: removed stray '*/' from comment")
elif FIXED in src:
    print("comment already fixed — nothing to patch")
else:
    die("neither broken nor fixed pattern found in datastream.h — inspect manually")

# ── Gate 1: normal build, zero warnings ─────────────────────────────
r = run(["make", "clean"]); r = run(["make"])
log = r.stdout + r.stderr
if r.returncode != 0: die(f"make failed:\n{log}")
warns = [l for l in log.splitlines() if "warning:" in l]
if warns: die("normal build produced warnings:\n" + "\n".join(warns))
print("GATE 1  normal build: zero warnings")

# ── Gate 2: demo output still matches the record (behavior-neutral) ──
expected_sha = read(TREE / "CURRENT_BASELINE_SHA.txt").strip()
r = subprocess.run(["./carbonsim"], cwd=str(TREE), capture_output=True)
if r.returncode != 0: die(f"carbonsim exited {r.returncode}")
if r.stderr: die("carbonsim stderr not empty:\n" + r.stderr.decode(errors="replace"))
actual_sha = hashlib.sha256(r.stdout).hexdigest()
if actual_sha != expected_sha:
    die(f"record BROKEN: demo stdout {actual_sha} != record {expected_sha}")
print("GATE 2  demo stdout byte-identical to record (datastream not yet wired in)")

# ── Gate 3: selftest builds ─────────────────────────────────────────
r = run(["make", "selftest"])
log = r.stdout + r.stderr
if r.returncode != 0: die(f"make selftest failed:\n{log}")
warns = [l for l in log.splitlines() if "warning:" in l]
if warns: die("selftest build produced warnings:\n" + "\n".join(warns))
print("GATE 3  selftest build: zero warnings")

# ── Gate 4: selftest passes ─────────────────────────────────────────
r = run(["./build/test_datastream"])
if r.returncode != 0: die(f"test_datastream FAILED:\n{r.stdout}\n{r.stderr}")
if "all passed" not in r.stdout:
    die(f"test_datastream did not report all-pass:\n{r.stdout}")
print("GATE 4  test_datastream: all checks passed")

# ── Gate 5: authoritative record harness ────────────────────────────
r = subprocess.run(["bash", str(ROOT / "s01_verify_record.sh")],
                   cwd=str(ROOT), capture_output=True, text=True)
print(r.stdout)
if r.returncode != 0 or "Record intact" not in r.stdout:
    die("s01_verify_record.sh did NOT report 'Record intact'")
print("GATE 5  s01 harness: Record intact")

print("\nPASS s41b: datastream writer builds clean, selftest green, record intact")
