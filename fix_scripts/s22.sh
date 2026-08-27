#!/usr/bin/env bash
set -euo pipefail
# s22: fix DNA duplex initial geometry — add B-DNA helical twist
#
# s21 correctly placed both base pairs with proper Watson-Crick geometry
# (placement check passed: G-C N1...N3=2.950, A-T N1...N3=2.900), and
# the final result after minimization+MD was correct (all H-bonds ~2.9 Å,
# both pairs PAIRED). But the initial PE was 2.6e11 eV — a catastrophic
# steric clash.
#
# Root cause: the A-T pair was placed with the same in-plane orientation
# as the G-C pair. Adenine sits directly above guanine at z=3.4 with no
# rotational offset, so ring atoms from the two pairs land nearly on top
# of each other. B-DNA has ~36° helical twist per base pair; applying
# that twist to the A-T pair prevents the in-plane overlap.
#
# RECORD-MOVING: Demo 17 output changes.

echo "=== s22: fix DNA duplex helical twist ==="

echo "[1/5] Editing src/main.c..."
python3 - <<'PYEOF'
import sys

path = 'src/main.c'
src = open(path).read()

start = src.find('static void demo_dna_duplex(void) {')
if start == -1:
    print('ERROR: demo_dna_duplex not found'); sys.exit(1)

end = src.find('\nstatic void ', start + 1)
if end == -1:
    end = src.find('\nint main(', start + 1)
if end == -1:
    print('ERROR: could not find end of demo_dna_duplex'); sys.exit(1)

new_func = r'''static void demo_dna_duplex(void) {
banner("DEMO 17: DNA duplex - minimal G-C and A-T base pair stack");
Simulation *sim = sim_create(64, 64);
sim->dielectric = 4.0;

/* ── G-C pair (z = 0) ─────────────────────────────────────────────── */
int g = sim_place_guanine(sim, vec3_zero());

Vec3 g_N1 = sim->atoms[g+6].position;
Vec3 g_O6 = sim->atoms[g+5].position;
Vec3 g_N1_to_O6_norm = vec3_normalize(vec3_sub(g_O6, g_N1));
Vec3 c_origin = vec3_add(g_N1, vec3_scale(g_N1_to_O6_norm, 2.9));
int c = sim_place_cytosine(sim, c_origin);

/* Align ring normals */
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

/* Azimuthal rotation: align C's N3->N4 with G's N1->O6 */
g_N1 = sim->atoms[g+6].position;
g_O6 = sim->atoms[g+5].position;
Vec3 c_N3 = sim->atoms[c+0].position;
Vec3 c_N4 = sim->atoms[c+5].position;
double az_angle = nb_signed_inplane_angle(vec3_sub(c_N4, c_N3),
                                          vec3_sub(g_O6, g_N1), n_G);
nb_transform_rigid(sim, c, 13, c_pivot, n_G, az_angle, vec3_zero());

/* Translate C so N3 lands at 2.9 A from G's N1 */
g_N1 = sim->atoms[g+6].position;
c_N3 = sim->atoms[c+0].position;
Vec3 to_c = vec3_sub(c_N3, g_N1);
Vec3 target_c = vec3_add(g_N1, vec3_scale(vec3_normalize(to_c), 2.9));
nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0,
                   vec3_sub(target_c, c_N3));

/* ── A-T pair (z = 3.4, B-DNA rise) ──────────────────────────────── */
double rise = 3.4;
int a = sim_place_adenine(sim, vec3(0.0, 0.0, rise));
int t = sim_place_thymine(sim, vec3(2.9, 0.0, rise));

/* Align ring normals */
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

/* Azimuthal rotation: align T's N3->O4 with A's N1->N6 */
Vec3 a_N1 = sim->atoms[a+6].position;
Vec3 a_N6 = sim->atoms[a+5].position;
Vec3 t_N3 = sim->atoms[t+3].position;
Vec3 t_O4 = sim->atoms[t+5].position;
double az_angle2 = nb_signed_inplane_angle(vec3_sub(t_O4, t_N3),
                                           vec3_sub(a_N6, a_N1), n_A);
nb_transform_rigid(sim, t, 15, t_pivot, n_A, az_angle2, vec3_zero());

/* Translate T so N3 lands at 2.9 A from A's N1 */
a_N1 = sim->atoms[a+6].position;
t_N3 = sim->atoms[t+3].position;
Vec3 to_t = vec3_sub(t_N3, a_N1);
Vec3 target_t = vec3_add(a_N1, vec3_scale(vec3_normalize(to_t), 2.9));
nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0,
                   vec3_sub(target_t, t_N3));

/* ── Helical twist: rotate A-T pair ~36 deg around z-axis ────────── */
/* B-DNA has ~36 deg twist per base pair. Without this the A-T pair
* sits directly above G-C with the same in-plane orientation and
* ring atoms from the two pairs nearly coincide (initial PE ~2.6e11
* eV in s21). The twist is a rigid rotation of both bases together
* around the helix axis, so intra-pair WC geometry is preserved. */
double twist = 36.0 * 3.14159265358979323846 / 180.0;
Vec3 twist_axis  = vec3(0.0, 0.0, 1.0);
Vec3 twist_pivot = vec3(0.0, 0.0, rise);
nb_transform_rigid(sim, a, 15, twist_pivot, twist_axis, twist, vec3_zero());
nb_transform_rigid(sim, t, 15, twist_pivot, twist_axis, twist, vec3_zero());

printf("  Placed: G-C at z=0, A-T at z=%.1f (36-deg helical twist)\n",
       rise);
printf("  Total atoms: %d\n", sim->num_atoms);

/* Validate initial geometry */
forces_calculate(sim);
double initial_pe = sim->potential_energy;
printf("  Initial PE: %.6f eV\n", initial_pe);

/* Steric-clash diagnostic */
double min_dist = 1.0e9; int ci = -1, cj = -1;
for (int i = 0; i < sim->num_atoms - 1; i++)
    for (int j = i + 1; j < sim->num_atoms; j++) {
        double d = vec3_dist(sim->atoms[i].position,
                             sim->atoms[j].position);
        if (d < min_dist) { min_dist = d; ci = i; cj = j; }
    }
printf("  Closest atom pair: %d-%d at %.3f A\n", ci, cj, min_dist);
if (min_dist < 1.0)
    printf("  WARNING: steric clash detected (< 1.0 A)\n");

/* Minimize if needed */
if (initial_pe > 100.0) {
    double min_pe = integrator_minimize(sim, 5000, 0.001, 0.01);
    printf("  After minimization: PE = %.6f eV\n", min_pe);
}

/* MD */
sim->dt = 0.5;
sim->thermostat.type               = THERMOSTAT_BERENDSEN;
sim->thermostat.target_temperature = 50.0;
sim->thermostat.tau                = 50.0;
integrator_maxwell_boltzmann(sim, 50.0, 42UL);
forces_calculate(sim);
sim->kinetic_energy = integrator_kinetic_energy(sim);
sim->total_energy   = sim->kinetic_energy + sim->potential_energy;
sim->temperature    = integrator_temperature(sim);

int N_steps = 800;
for (int step = 0; step < N_steps; step++) integrator_step(sim);
printf("  After %d steps: PE=%.6f eV  T=%.2f K\n",
       N_steps, sim->potential_energy, sim->temperature);

/* H-bond report */
printf("\nG-C H-bonds:");
printf(" N1...N3=%.3f", vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
printf("  O6...N4=%.3f", vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));
printf("  N2...O2=%.3f", vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
printf("\nA-T H-bonds:");
printf(" N1...N3=%.3f", vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position));
printf("  N6...O4=%.3f", vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position));

/* Planarity */
int gr[] = {g+0, g+1, g+6};
int cr[] = {c+0, c+1, c+2};
int ar[] = {a+0, a+1, a+6};
int tr[] = {t+0, t+1, t+3};
printf("\nPlanarity: G=%.4f  C=%.4f  A=%.4f  T=%.4f\n",
       nb_planarity_deviation(sim, gr, 3),
       nb_planarity_deviation(sim, cr, 3),
       nb_planarity_deviation(sim, ar, 3),
       nb_planarity_deviation(sim, tr, 3));

/* Verdict */
double gc1 = vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position);
double gc2 = vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position);
double gc3 = vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position);
double at1 = vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position);
double at2 = vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position);
int gc_ok = (gc1 < 3.3) && (gc2 < 3.3) && (gc3 < 3.3);
int at_ok = (at1 < 3.3) && (at2 < 3.3);
printf("\nVERDICT: G-C %s, A-T %s\n",
       gc_ok ? "PAIRED" : "NOT PAIRED",
       at_ok ? "PAIRED" : "NOT PAIRED");
if (gc_ok && at_ok)
    printf("--> A real, H-bonded two-base-pair DNA stack.\n");
else
    printf("--> Bases drifted apart or sheared; NOT a true duplex.\n");

sim_destroy(sim);
}
'''

src = src[:start] + new_func + src[end:]
open(path, 'w').write(src)
print("  [OK] demo_dna_duplex replaced with helical-twist version.")
PYEOF

echo "[2/5] Building..."
make clean >/dev/null
if ! make 2>&1; then
    echo "ERROR: Build failed"; exit 1
fi
echo "  Build clean."

echo "[3/5] Running DNA demo..."
./carbonsim --dna

echo "[4/5] Checking record..."
if [ -f CURRENT_BASELINE_SHA.txt ]; then
    BASELINE=$(cat CURRENT_BASELINE_SHA.txt)
    ./carbonsim > /tmp/s22_record.txt 2>/dev/null
    NOW=$(sha256sum /tmp/s22_record.txt | awk '{print $1}')
    if [ "$NOW" = "$BASELINE" ]; then
        echo "  [OK] Record unchanged."
    else
        echo "  Record moved (expected: Demo 17 output changed)."
        echo "  Old: $BASELINE"
        echo "  New: $NOW"
    fi
else
    echo "  (no CURRENT_BASELINE_SHA.txt — skipping record check)"
fi

echo "[5/5] Done."
echo "=== s22 complete ==="
