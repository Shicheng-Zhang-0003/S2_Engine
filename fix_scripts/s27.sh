#!/usr/bin/env bash
set -euo pipefail
# s27: fix A-T pair stacking direction.
# The helix axis is along Y (ring normals point ±Y), but the A-T pair
# was placed at z=rise instead of y=rise. This put the A-T pair in the
# same XZ slab as G-C, causing 134 inter-base clashes.
#
# Fix: place A-T pair at y=rise, and place T far away along X before
# aligning (same approach as G-C pair uses for cytosine).
echo "=== s27: fix A-T stacking direction (z=rise -> y=rise) ==="

python3 - <<'PYEOF'
import sys
with open('src/main.c') as f:
    src = f.read()

# Fix the A-T pair placement in demo_dna_duplex.
# Change: vec3(0.0, 0.0, rise) -> vec3(0.0, rise, 0.0) for adenine
# Change: vec3(2.9, 0.0, rise) -> vec3(15.0, rise, 0.0) for thymine
#   (place T far away along X, at the same Y as A, then translate)

old_a = 'int a = sim_place_adenine(sim, vec3(0.0, 0.0, rise));'
new_a = 'int a = sim_place_adenine(sim, vec3(0.0, rise, 0.0));'

old_t = 'int t = sim_place_thymine(sim, vec3(2.9, 0.0, rise));'
new_t = 'int t = sim_place_thymine(sim, vec3(15.0, rise, 0.0));'

if old_a in src:
    src = src.replace(old_a, new_a, 1)
    print("  [OK] Adenine placement: z=rise -> y=rise")
else:
    print("  WARN: adenine placement pattern not found (may already be fixed)")

if old_t in src:
    src = src.replace(old_t, new_t, 1)
    print("  [OK] Thymine placement: (2.9,0,rise) -> (15,rise,0)")
else:
    print("  WARN: thymine placement pattern not found (may already be fixed)")

with open('src/main.c', 'w') as f:
    f.write(src)
PYEOF

echo ""
echo "Building..."
make clean >/dev/null
make 2>&1 | tail -3

echo ""
echo "Running DNA demo..."
./carbonsim --dna

echo ""
echo "Running diagnostic to verify..."
./carbonsim --diag 2>/dev/null | head -60

echo "=== s27 complete ==="
