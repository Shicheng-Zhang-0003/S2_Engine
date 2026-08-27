#!/usr/bin/env bash
set -euo pipefail
# s26_diag: diagnostic dump for DNA duplex geometry
# Adds a --diag mode that prints ring normals, centroids, key atom
# positions, and all close-contact pairs BEFORE any minimization/MD.
# Run: bash s26_diag.sh

echo "=== s26_diag: adding diagnostic mode to demo_dna_duplex ==="

python3 - <<'PYEOF'
import sys

path = 'src/main.c'
src = open(path).read()

# Find the end of demo_dna_duplex to insert the diagnostic function after it
marker = 'int main(int argc, char *argv[]) {'
pos = src.find(marker)
if pos == -1:
    print('FAIL: main() not found'); sys.exit(1)

diag_func = '''
/* ════════════════════════════════════════════════════════════════════════════
* DIAGNOSTIC: DNA duplex geometry dump — no MD, no minimization.
* Prints ring normals, centroids, key atom positions, and all close pairs.
* ════════════════════════════════════════════════════════════════════════════ */
static void print_v3(const char *label, Vec3 v) {
    printf("  %-30s (%8.4f, %8.4f, %8.4f)\\n", label, v.x, v.y, v.z);
}

static void demo_dna_duplex_diag(void) {
    banner("DIAG: DNA duplex geometry dump (no MD)");
    Simulation *sim = sim_create(64, 64);
    sim->dielectric = 4.0;
    double rise = 3.4;

    /* ── G-C pair (same placement as production code) ─────────────────── */
    int g = sim_place_guanine(sim, vec3_zero());
    Vec3 g_N1 = sim->atoms[g+6].position;
    Vec3 g_O6 = sim->atoms[g+5].position;
    Vec3 g_N1_to_O6_norm = vec3_normalize(vec3_sub(g_O6, g_N1));
    Vec3 c_origin = vec3_add(g_N1, vec3_scale(g_N1_to_O6_norm, 2.9));
    int c = sim_place_cytosine(sim, c_origin);

    int g_ring[3] = {g+0, g+1, g+6};
    int c_ring[3] = {c+0, c+1, c+2};
    Vec3 n_G = nb_ring_normal(sim, g_ring);
    Vec3 n_C = nb_ring_normal(sim, c_ring);
    double cosang = vec3_dot(n_G, n_C);
    Vec3 axis = vec3_cross(n_C, n_G);
    double angle;
    if (vec3_norm(axis) < 1.0e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(n_C, vec3(1,0,0));
        if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
    } else {
        angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    }
    Vec3 c_pivot = sim->atoms[c+0].position;
    nb_transform_rigid(sim, c, 13, c_pivot, axis, angle, vec3_zero());

    g_N1 = sim->atoms[g+6].position;
    g_O6 = sim->atoms[g+5].position;
    Vec3 c_N3 = sim->atoms[c+0].position;
    Vec3 c_N4 = sim->atoms[c+5].position;
    double az = nb_signed_inplane_angle(vec3_sub(c_N4, c_N3),
                                        vec3_sub(g_O6, g_N1), n_G);
    nb_transform_rigid(sim, c, 13, c_pivot, n_G, az, vec3_zero());

    g_N1 = sim->atoms[g+6].position;
    c_N3 = sim->atoms[c+0].position;
    Vec3 to_c = vec3_sub(c_N3, g_N1);
    Vec3 target_c = vec3_add(g_N1, vec3_scale(vec3_normalize(to_c), 2.9));
    nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0,
                       vec3_sub(target_c, c_N3));

    /* ── A-T pair (same placement as production code) ─────────────────── */
    int a = sim_place_adenine(sim, vec3(0.0, 0.0, rise));
    int t = sim_place_thymine(sim, vec3(15.0, 0.0, rise));

    int a_ring[3] = {a+0, a+1, a+6};
    int t_ring[3] = {t+0, t+1, t+3};
    Vec3 n_A = nb_ring_normal(sim, a_ring);
    Vec3 n_T = nb_ring_normal(sim, t_ring);
    cosang = vec3_dot(n_A, n_T);
    axis = vec3_cross(n_T, n_A);
    if (vec3_norm(axis) < 1.0e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(n_T, vec3(1,0,0));
        if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
    } else {
        angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    }
    Vec3 t_pivot = sim->atoms[t+3].position;
    nb_transform_rigid(sim, t, 15, t_pivot, axis, angle, vec3_zero());

    Vec3 a_N1 = sim->atoms[a+6].position;
    Vec3 a_N6 = sim->atoms[a+5].position;
    Vec3 t_N3 = sim->atoms[t+3].position;
    Vec3 t_O4 = sim->atoms[t+5].position;
    double az2 = nb_signed_inplane_angle(vec3_sub(t_O4, t_N3),
                                         vec3_sub(a_N6, a_N1), n_A);
    nb_transform_rigid(sim, t, 15, t_pivot, n_A, az2, vec3_zero());

    a_N1 = sim->atoms[a+6].position;
    t_N3 = sim->atoms[t+3].position;
    Vec3 to_t = vec3_sub(t_N3, a_N1);
    Vec3 target_t = vec3_add(a_N1, vec3_scale(vec3_normalize(to_t), 2.9));
    nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0,
                       vec3_sub(target_t, t_N3));

    printf("  Total atoms: %d\\n\\n", sim->num_atoms);

    /* ── 1. RING NORMALS ─────────────────────────────────────────────── */
    printf("=== RING NORMALS (unit vectors) ===\\n");
    /* Recompute from current positions */
    n_G = nb_ring_normal(sim, g_ring);
    n_C = nb_ring_normal(sim, c_ring);
    n_A = nb_ring_normal(sim, a_ring);
    n_T = nb_ring_normal(sim, t_ring);
    print_v3("G ring normal:", n_G);
    print_v3("C ring normal:", n_C);
    print_v3("A ring normal:", n_A);
    print_v3("T ring normal:", n_T);

    /* ── 2. BASE CENTROIDS ───────────────────────────────────────────── */
    printf("\\n=== BASE CENTROIDS ===\\n");
    /* Guanine: 16 atoms, indices g..g+15 */
    Vec3 cg = vec3_zero();
    for (int i = g; i < g+16; i++) cg = vec3_add(cg, sim->atoms[i].position);
    cg = vec3_scale(cg, 1.0/16.0);
    print_v3("G centroid:", cg);

    /* Cytosine: 13 atoms, indices c..c+12 */
    Vec3 cc = vec3_zero();
    for (int i = c; i < c+13; i++) cc = vec3_add(cc, sim->atoms[i].position);
    cc = vec3_scale(cc, 1.0/13.0);
    print_v3("C centroid:", cc);

    /* Adenine: 15 atoms, indices a..a+14 */
    Vec3 ca = vec3_zero();
    for (int i = a; i < a+15; i++) ca = vec3_add(ca, sim->atoms[i].position);
    ca = vec3_scale(ca, 1.0/15.0);
    print_v3("A centroid:", ca);

    /* Thymine: 15 atoms, indices t..t+14 */
    Vec3 ct = vec3_zero();
    for (int i = t; i < t+15; i++) ct = vec3_add(ct, sim->atoms[i].position);
    ct = vec3_scale(ct, 1.0/15.0);
    print_v3("T centroid:", ct);

    /* ── 3. KEY ATOM POSITIONS ───────────────────────────────────────── */
    printf("\\n=== KEY ATOM POSITIONS ===\\n");
    /* Guanine: N9=0,C8=1,N7=2,C5=3,C6=4,O6=5,N1=6,C2=7,N2=8,N3=9,C4=10,
     *          HN9=11,H8=12,HN1=13,HN21=14,HN22=15 */
    print_v3("G:N1",  sim->atoms[g+6].position);
    print_v3("G:O6",  sim->atoms[g+5].position);
    print_v3("G:N2",  sim->atoms[g+8].position);
    print_v3("G:HN1", sim->atoms[g+13].position);

    /* Cytosine: N3=0,C4=1,N1=2,C2=3,O2=4,N4=5,C5=6,C6=7,
     *           HN41=8,HN42=9,H5=10,H6=11,HN1=12 */
    print_v3("C:N3",  sim->atoms[c+0].position);
    print_v3("C:N4",  sim->atoms[c+5].position);
    print_v3("C:O2",  sim->atoms[c+4].position);

    /* Adenine: N9=0,C8=1,N7=2,C5=3,C6=4,N6=5,N1=6,C2=7,N3=8,C4=9,
     *          HN9=10,H8=11,HN61=12,HN62=13,H2=14 */
    print_v3("A:N1",  sim->atoms[a+6].position);
    print_v3("A:N6",  sim->atoms[a+5].position);

    /* Thymine: N1=0,C2=1,O2=2,N3=3,C4=4,O4=5,C5=6,C6=7,
     *          HN1=8,HN3=9,H6=10,CM=11,HM1=12,HM2=13,HM3=14 */
    print_v3("T:N3",  sim->atoms[t+3].position);
    print_v3("T:O4",  sim->atoms[t+5].position);
    print_v3("T:HN3", sim->atoms[t+9].position);

    /* ── 4. INTER-BASE H-BOND DISTANCES ──────────────────────────────── */
    printf("\\n=== INTER-BASE H-BOND DISTANCES ===\\n");
    printf("  G:N1 to C:N3  = %.4f A\\n",
           vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
    printf("  G:O6 to C:N4  = %.4f A\\n",
           vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));
    printf("  G:N2 to C:O2  = %.4f A\\n",
           vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
    printf("  A:N1 to T:N3  = %.4f A\\n",
           vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position));
    printf("  A:N6 to T:O4  = %.4f A\\n",
           vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position));

    /* ── 5. ALL PAIR DISTANCES < 2.5 A (potential clashes) ───────────── */
    printf("\\n=== ALL INTER-BASE PAIR DISTANCES < 2.5 A ===\\n");
    int clash_count = 0;
    for (int i = 0; i < sim->num_atoms - 1; i++) {
        for (int j = i + 1; j < sim->num_atoms; j++) {
            /* Skip intra-base pairs (same molecule) */
            int bi = -1, bj = -1;
            if (i >= g && i < g+16) bi = 0;
            else if (i >= c && i < c+13) bi = 1;
            else if (i >= a && i < a+15) bi = 2;
            else if (i >= t && i < t+15) bi = 3;
            if (j >= g && j < g+16) bj = 0;
            else if (j >= c && j < c+13) bj = 1;
            else if (j >= a && j < a+15) bj = 2;
            else if (j >= t && j < t+15) bj = 3;
            if (bi == bj) continue; /* same base, skip */
            double d = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
            if (d < 2.5) {
                printf("  atom %2d (%s, base %d) to atom %2d (%s, base %d) = %.4f A\\n",
                       i, sim->atoms[i].element->symbol, bi,
                       j, sim->atoms[j].element->symbol, bj, d);
                clash_count++;
            }
        }
    }
    printf("  Total clashes (< 2.5 A): %d\\n", clash_count);

    /* ── 6. CLOSEST 20 INTER-BASE PAIRS (sorted) ─────────────────────── */
    printf("\\n=== CLOSEST 20 INTER-BASE PAIRS ===\\n");
    typedef struct { int i, j; double d; } Pair;
    Pair pairs[4096];
    int np = 0;
    for (int i = 0; i < sim->num_atoms - 1; i++) {
        for (int j = i + 1; j < sim->num_atoms; j++) {
            int bi = -1, bj = -1;
            if (i >= g && i < g+16) bi = 0;
            else if (i >= c && i < c+13) bi = 1;
            else if (i >= a && i < a+15) bi = 2;
            else if (i >= t && i < t+15) bi = 3;
            if (j >= g && j < g+16) bj = 0;
            else if (j >= c && j < c+13) bj = 1;
            else if (j >= a && j < a+15) bj = 2;
            else if (j >= t && j < t+15) bj = 3;
            if (bi == bj) continue;
            double d = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
            if (np < 4096) { pairs[np].i = i; pairs[np].j = j; pairs[np].d = d; np++; }
        }
    }
    /* Simple insertion sort for top 20 */
    for (int i = 0; i < np - 1; i++) {
        int min_idx = i;
        for (int j = i + 1; j < np; j++)
            if (pairs[j].d < pairs[min_idx].d) min_idx = j;
        if (min_idx != i) { Pair tmp = pairs[i]; pairs[i] = pairs[min_idx]; pairs[min_idx] = tmp; }
        if (i >= 19) break;
    }
    for (int i = 0; i < 20 && i < np; i++) {
        printf("  atom %2d (%s) to atom %2d (%s) = %.4f A\\n",
               pairs[i].i, sim->atoms[pairs[i].i].element->symbol,
               pairs[i].j, sim->atoms[pairs[i].j].element->symbol,
               pairs[i].d);
    }

    sim_destroy(sim);
    printf("\\n=== DIAG COMPLETE ===\\n");
}

'''

src = src[:pos] + diag_func + src[pos:]

# Add --diag to main()
old_main_check = 'if (argc > 1 && strcmp(argv[1], "--dna") == 0) {'
new_main_check = '''if (argc > 1 && strcmp(argv[1], "--diag") == 0) {
        demo_dna_duplex_diag();
        return 0;
    }
    if (argc > 1 && strcmp(argv[1], "--dna") == 0) {'''
src = src.replace(old_main_check, new_main_check, 1)

open(path, 'w').write(src)
print('  [OK] Diagnostic function added.')
PYEOF

echo ""
echo "Building..."
make clean >/dev/null
make 2>&1 | tail -3

echo ""
echo "Running diagnostic..."
./carbonsim --diag

echo ""
echo "=== s26_diag complete ==="
