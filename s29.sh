#!/usr/bin/env bash
set -euo pipefail
echo "=== s29: remove stale --diag dispatch (function deleted by s28) ==="

python3 - <<'PYEOF'
with open('src/main.c') as f:
    src = f.read()

old = '''    if (argc > 1 && strcmp(argv[1], "--diag") == 0) {
        demo_dna_duplex_diag();
        return 0;
    }
    '''

if old in src:
    src = src.replace(old, '', 1)
    print("  [OK] Removed --diag dispatch.")
else:
    # Try alternate whitespace
    old2 = '''if (argc > 1 && strcmp(argv[1], "--diag") == 0) {
        demo_dna_duplex_diag();
        return 0;
    }'''
    if old2 in src:
        src = src.replace(old2, '', 1)
        print("  [OK] Removed --diag dispatch (alt whitespace).")
    else:
        print("  WARN: --diag dispatch not found (may already be removed).")

with open('src/main.c', 'w') as f:
    f.write(src)
PYEOF

echo "Rebuilding..."
make clean >/dev/null
make 2>&1 | tail -3
echo "Running DNA demo..."
./carbonsim --dna
echo "=== s29 complete ==="
