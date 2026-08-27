#!/usr/bin/env bash
set -euo pipefail
# s21: rebuild demo_dna_duplex using Demo 7's VALIDATED Watson-Crick
# placement (the only pairing logic in this tree proven to form a real
# H-bonded pair), and add a HARD verification gate so an unpaired
# result can no longer pass silently.
echo "=== s21: DNA duplex with validated WC placement + verification gate ==="

python3 - <<'PYEOF'
import sys
path = 'src/main.c'
src = open(path).read()

start = src.find('static void demo_dna_duplex(void) {')
if start == -1:
    print('FAIL: demo_dna_duplex not found'); sys.exit(1)
end = src.find('\nint main(', start + 1)
if end == -1:
    print('FAIL: could not find end of demo_dna_duplex'); sys.exit(1)

new_func = r'''static void demo_dna_duplex(void) {
banner("DEMO 17: DNA duplex - minimal G-C and A-T base pair stack");
Simulation *sim = sim_create(64, 64);
sim->dielectric = 4.0;

/* ── Helper: place a validated Watson-Crick pair ─────────────────────
* Ports the EXACT placement strategy from demo_basepairing (Demo 7),
* which is the only pairing code in this tree verified to produce a
* real H-bonded pair (G-C N1...N3 = 2.878 A there). The duplex demo
* previously used an incorrect "place along N1->O6" heuristic that
* left the WC edges mis-facing and the bases ~4.5 A apart.        */

/* --- G-C pair (3 H-bonds) --- */
int g = sim_place_guanine(sim, vec3_zero());
int c = sim_place_cytosine(sim, vec3(15.0, 0.0, 0.0));
{
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
    } else angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    Vec3 c_pivot = sim->atoms[c+0].position;
    nb_transform_rigid(sim, c, 13, c_pivot, axis, angle, vec3_zero());
    /* azimuthal: align C's N3->N4 with G's N1->O6 */
    Vec3 g_N1 = sim->atoms[g+6].position, g_O6 = sim->atoms[g+5].position;
    Vec3 c_N3 = sim->atoms[c+0].position, c_N4 = sim->atoms[c+5].position;
    Vec3 v_G = vec3_sub(g_O6, g_N1);
    Vec3 v_C = vec3_sub(c_N4, c_N3);
    double az = nb_signed_inplane_angle(v_C, v_G, n_G);
    nb_transform_rigid(sim, c, 13, c_pivot, n_G, az, vec3_zero());
    /* translate: C:N3 to 2.95 A along G's N1->HN1 donor direction */
    Vec3 g_HN1 = sim->atoms[g+13].position;
    Vec3 donor_dir = vec3_normalize(vec3_sub(g_HN1, g_N1));
    Vec3 target = vec3_add(g_N1, vec3_scale(donor_dir, 2.95));
    Vec3 trans = vec3_sub(target, sim->atoms[c+0].position);
    nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0, trans);
}

/* --- A-T pair (2 H-bonds), stacked at z = 3.4 A (B-DNA rise) --- */
double rise = 3.4;
int a = sim_place_adenine(sim, vec3(0.0, 0.0, rise));
int t = sim_place_thymine(sim, vec3(15.0, 0.0, rise));
{
    int a_ring[3] = {a+0, a+1, a+6};
    int t_ring[3] = {t+0, t+1, t+3};
    Vec3 n_A = nb_ring_normal(sim, a_ring);
    Vec3 n_T = nb_ring_normal(sim, t_ring);
    double cosang = vec3_dot(n_A, n_T);
    Vec3 axis = vec3_cross(n_T, n_A);
    double angle;
    if (vec3_norm(axis) < 1.0e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(n_T, vec3(1,0,0));
        if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
    } else angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    Vec3 t_pivot = sim->atoms[t+3].position;
    nb_transform_rigid(sim, t, 15, t_pivot, axis, angle, vec3_zero());
    /* azimuthal: align T's N3->O4 with A's N1->N6 */
    Vec3 a_N1 = sim->atoms[a+6].position, a_N6 = sim->atoms[a+5].position;
    Vec3 t_N3 = sim->atoms[t+3].position, t_O4 = sim->atoms[t+5].position;
    Vec3 v_A = vec3_sub(a_N6, a_N1);
    Vec3 v_T = vec3_sub(t_O4, t_N3);
    double az = nb_signed_inplane_angle(v_T, v_A, n_A);
    nb_transform_rigid(sim, t, 15, t_pivot, n_A, az, vec3_zero());
    /* translate: T:N3 so its H points at A:N1 at 2.90 A */
    Vec3 t_HN3 = sim->atoms[t+9].position;
    Vec3 donor_dir = vec3_normalize(vec3_sub(t_HN3, t_N3));
    Vec3 target_N3 = vec3_sub(a_N1, vec3_scale(donor_dir, 2.90));
    Vec3 trans = vec3_sub(target_N3, sim->atoms[t+3].position);
    nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0, trans);
}

printf("  Placed: G-C pair at z=0, A-T pair at z=%.1f\n", rise);
printf("  Total atoms: %d\n", sim->num_atoms);

/* ── HARD VERIFICATION GATE (before any dynamics) ────────────────────
* If the placement didn't put the primary WC heavy-atom contacts in
* the real H-bond window, STOP and say so. A duplex demo that prints
* a pretty table for unpaired bases is a silent physics failure.  */
{
    double gc = vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position);
    double at = vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position);
    printf("  Placement check: G-C N1...N3=%.3f A   A-T N1...N3=%.3f A\n", gc, at);
    if (!(gc > 2.4 && gc < 3.4) || !(at > 2.4 && at < 3.4)) {
        printf("  FAIL: initial placement is NOT Watson-Crick paired.\n");
        printf("        Expected primary contacts in [2.4, 3.4] A.\n");
        printf("        Got G-C=%.3f, A-T=%.3f. Aborting before MD.\n", gc, at);
        sim_destroy(sim);
        return;
    }
    printf("  [OK] Placement is WC-paired. Proceeding.\n");
}

/* Relax residual clashes, then short MD */
forces_calculate(sim);
printf("  Initial PE: %.6f eV\n", sim->potential_energy);
double min_pe = integrator_minimize(sim, 5000, 0.001, 0.01);
printf("  After minimization: PE = %.6f eV\n", min_pe);

sim->dt = 0.5;
sim->thermostat.type               = THERMOSTAT_BERENDSEN;
sim->thermostat.target_temperature = 50.0;
sim->thermostat.tau                = 50.0;
integrator_maxwell_boltzmann(sim, 50.0, 42UL);
forces_calculate(sim);
sim->kinetic_energy = integrator_kinetic_energy(sim);
sim->temperature    = integrator_temperature(sim);

int N_steps = 800;
for (int step = 0; step < N_steps; step++) integrator_step(sim);
printf("  After %d steps: PE=%.6f eV  T=%.2f K\n",
       N_steps, sim->potential_energy, sim->temperature);

/* Report H-bond distances and GATE them too */
double gc1 = vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position);
double gc2 = vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position);
double gc3 = vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position);
double at1 = vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position);
double at2 = vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position);
printf("\nG-C H-bonds: N1...N3=%.3f  O6...N4=%.3f  N2...O2=%.3f\n", gc1, gc2, gc3);
printf("A-T H-bonds: N1...N3=%.3f  N6...O4=%.3f\n", at1, at2);

int gc_ok = (gc1 < 3.3) && (gc2 < 3.3) && (gc3 < 3.3);
int at_ok = (at1 < 3.3) && (at2 < 3.3);
printf("\nVERDICT: G-C %s, A-T %s\n",
       gc_ok ? "PAIRED (all contacts < 3.3 A)" : "NOT PAIRED",
       at_ok ? "PAIRED (all contacts < 3.3 A)" : "NOT PAIRED");
if (gc_ok && at_ok)
    printf("--> A real, H-bonded two-base-pair DNA stack.\n");
else
    printf("--> Bases drifted apart or sheared; NOT a true duplex.\n");

int g_ring2[] = {g+0, g+1, g+6};
int c_ring2[] = {c+0, c+1, c+2};
int a_ring2[] = {a+0, a+1, a+6};
int t_ring2[] = {t+0, t+1, t+3};
printf("Planarity: G=%.4f  C=%.4f  A=%.4f  T=%.4f\n",
       nb_planarity_deviation(sim, g_ring2, 3),
       nb_planarity_deviation(sim, c_ring2, 3),
       nb_planarity_deviation(sim, a_ring2, 3),
       nb_planarity_deviation(sim, t_ring2, 3));
sim_destroy(sim);
}
'''

src = src[:start] + new_func + src[end:]
open(path, 'w').write(src)
print("  [OK] demo_dna_duplex replaced with validated-placement + gated version.")
PYEOF

echo "[2/3] Building..."
make clean >/dev/null
if ! make 2>&1 | tee /tmp/s21_build.log | grep -E "error|warning"; then
    echo "  Build clean."
fi
if grep -qE "error" /tmp/s21_build.log; then
    echo "  BUILD FAILED:"; cat /tmp/s21_build.log; exit 1
fi

echo "[3/3] Running DNA demo..."
./carbonsim --dna
echo "=== s21 complete ==="
