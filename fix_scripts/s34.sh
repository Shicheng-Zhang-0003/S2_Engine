#!/usr/bin/env bash
set -euo pipefail
# s34: incorporate Demo 17 into the main carbonsim executable.
# Currently demo_dna_duplex() is only reachable via --dna flag.
# This moves it into the main demo sequence (after Demo 12) and
# removes the --dna flag, so it runs as part of the standard flow.
# RECORD-MOVING: main executable output changes.
echo "=== s34: incorporate Demo 17 into main executable ==="

python3 - <<'PYEOF'
import sys
with open('src/main.c') as f:
    src = f.read()

# 1. Remove the --dna flag check from main()
old_dna_flag = '''    if (argc > 1 && strcmp(argv[1], "--dna") == 0) {
        demo_dna_duplex();
        return 0;
    }
'''
if old_dna_flag in src:
    src = src.replace(old_dna_flag, '', 1)
    print('  [OK] Removed --dna flag from main().')
else:
    print('  WARN: --dna flag pattern not found (may already be removed).')

# 2. Add demo_dna_duplex() call after Demo 12 in the main flow.
# Find the last demo call (demo_kcsa_filter) and add after it.
anchor = 'demo_kcsa_filter();'
if anchor in src:
    insertion = anchor + '\n    demo_dna_duplex();'
    src = src.replace(anchor, insertion, 1)
    print('  [OK] Added demo_dna_duplex() after Demo 12 in main flow.')
else:
    print('  FAIL: could not find demo_kcsa_filter() call in main().')
    sys.exit(1)

with open('src/main.c', 'w') as f:
    f.write(src)
PYEOF

echo ""
echo "Building..."
make clean >/dev/null
make 2>&1 | tail -3

echo ""
echo "Running main executable (first 50 lines + last 30 lines)..."
./carbonsim 2>/dev/null | head -50
echo "  ..."
./carbonsim 2>/dev/null | tail -30

echo ""
echo "=== s34 complete ==="
echo "Demo 17 is now part of the main executable."
echo "Run './carbonsim' to see all 13 demos."
echo "The --dna flag is no longer needed."
