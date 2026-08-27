#!/usr/bin/env bash
set -euo pipefail
# s26_diag: print the diagnostic values needed to fix demo_dna_duplex
# Inserts diagnostics, builds, runs, then reverts main.c.

echo "=== s26_diag: extract placement diagnostics ==="

cp src/main.c src/main.c.bak

python3 - <<'PYEOF'
import sys

with open('src/main.c') as f:
    src = f.read()

# Insert diagnostics right after the last nb_transform_rigid for thymine,
# before the "Placed:" printf.
anchor = 'printf("  Placed: G-C pair at z=0, A-T pair at z=3.4\\n");'
assert anchor in src, "anchor not found"

diag = '''
/* ══ DIAGNOSTIC BLOCK (s26_diag) ══ */
{
    printf("\\n══ DIAGNOSTICS ══\\n");

    /* 1. Ring normals of all four bases */
    int g_r[3] = {g+0, g+1, g+6};
    int c_r[3] = {c+0, c+1, c+2};
    int a_r[3] = {a+0, a+1, a+6};
    int t_r[3] = {t+0, t+1, t+3};
    Vec3 nG = nb_ring_normal(sim, g_r);
    Vec3 nC = nb_ring_normal(sim, c_r);
    Vec3 nA = nb_ring_normal(sim, a_r);
    Vec3 nT = nb_ring_normal(sim, t_r);
    printf("  Ring normal G: (%+.6f, %+.6f, %+.6f)\\n", nG.x, nG.y, nG.z);
    printf("  Ring normal C: (%+.6f, %+.6f, %+.6f)\\n", nC.x, nC.y, nC.z);
    printf("  Ring normal A: (%+.6f, %+.6f, %+.6f)\\n", nA.x, nA.y, nA.z);
    printf("  Ring normal T: (%+.6f, %+.6f, %+.6f)\\n", nT.x, nT.y, nT.z);

    /* 2. Positions of the clashing atoms (6 and 29) */
    printf("  Atom 6  (G:N1): (%+.6f, %+.6f, %+.6f)\\n",
           sim->atoms[6].position.x, sim->atoms[6].position.y, sim->atoms[6].position.z);
    printf("  Atom 29 (A:N9): (%+.6f, %+.6f, %+.6f)\\n",
           sim->atoms[29].position.x, sim->atoms[29].position.y, sim->atoms[29].position.z);
    printf("  Distance 6-29: %.6f A\\n",
           vec3_dist(sim->atoms[6].position, sim->atoms[29].position));

    /* 3. Centroid of each base */
    Vec3 cG = vec3_zero(), cC = vec3_zero(), cA = vec3_zero(), cT = vec3_zero();
    for (int i = g; i < g+16; i++) vec3_iadd(&cG, sim->atoms[i].position);
    for (int i = c; i < c+13; i++) vec3_iadd(&cC, sim->atoms[i].position);
    for (int i = a; i < a+15; i++) vec3_iadd(&cA, sim->atoms[i].position);
    for (int i = t; i < t+15; i++) vec3_iadd(&cT, sim->atoms[i].position);
    cG = vec3_scale(cG, 1.0/16.0);
    cC = vec3_scale(cC, 1.0/13.0);
    cA = vec3_scale(cA, 1.0/15.0);
    cT = vec3_scale(cT, 1.0/15.0);
    printf("  Centroid G: (%+.6f, %+.6f, %+.6f)\\n", cG.x, cG.y, cG.z);
    printf("  Centroid C: (%+.6f, %+.6f, %+.6f)\\n", cC.x, cC.y, cC.z);
    printf("  Centroid A: (%+.6f, %+.6f, %+.6f)\\n", cA.x, cA.y, cA.z);
    printf("  Centroid T: (%+.6f, %+.6f, %+.6f)\\n", cT.x, cT.y, cT.z);

    /* 4. Also print the 10 closest inter-pair atom pairs */
    printf("\\n  10 closest inter-pair contacts:\\n");
    double dists[59*59]; int pi[59*59], pj[59*59]; int np = 0;
    for (int i = 0; i < sim->num_atoms; i++)
        for (int j = i+1; j < sim->num_atoms; j++) {
            dists[np] = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
            pi[np] = i; pj[np] = j; np++;
        }
    /* simple selection sort for top 10 */
    for (int k = 0; k < 10 && k < np; k++) {
        int best = k;
        for (int m = k+1; m < np; m++)
            if (dists[m] < dists[best]) best = m;
        double tmp = dists[k]; dists[k] = dists[best]; dists[best] = tmp;
        int ti = pi[k]; pi[k] = pi[best]; pi[best] = ti;
        int tj = pj[k]; pj[k] = pj[best]; pj[best] = tj;
    }
    for (int k = 0; k < 10 && k < np; k++)
        printf("    atoms %2d-%2d: %.4f A\\n", pi[k], pj[k], dists[k]);

    printf("══ END DIAGNOSTICS ══\\n\\n");
}
/* ══ END DIAGNOSTIC BLOCK ══ */

'''

src = src.replace(anchor, diag + anchor, 1)
with open('src/main.c', 'w') as f:
    f.write(src)
print("  [OK] Diagnostics inserted.")
PYEOF

echo "[2/4] Building..."
make clean >/dev/null
make 2>&1 | grep -E "error" && { echo "BUILD FAILED"; mv src/main.c.bak src/main.c; exit 1; }
echo "  Build OK."

echo "[3/4] Running --dna..."
./carbonsim --dna

echo "[4/4] Reverting main.c..."
mv src/main.c.bak src/main.c
echo "  Reverted."
echo "=== s26_diag complete. Paste the DIAGNOSTICS block output above. ==="
