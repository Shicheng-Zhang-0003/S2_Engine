#!/usr/bin/env bash
set -euo pipefail
# s20: replace demo_dna_duplex with a corrected version.
# Fixes: initial geometry clashes, A-T placement, MD parameters, output.

echo "=== s20: fix demo_dna_duplex ==="

python3 - <<'PYEOF'
import sys, re

with open('src/main.c') as f:
    src = f.read()

start = src.find('static void demo_dna_duplex(void) {')
if start == -1:
    print("ERROR: demo_dna_duplex not found")
    sys.exit(1)

end = src.find('\nstatic void demo_', start + 1)
if end == -1:
    end = src.find('\nint main(', start + 1)
if end == -1:
    print("ERROR: could not find end of demo_dna_duplex")
    sys.exit(1)

new_func = '''static void demo_dna_duplex(void) {
    banner("DEMO 17: DNA duplex - minimal G-C and A-T base pair stack");

    Simulation *sim = sim_create(64, 64);
    sim->dielectric = 4.0;

    /* ── G-C pair ─────────────────────────────────────────────────────── */
    int g = sim_place_guanine(sim, vec3_zero());

    /* Place cytosine so its N3 (index c+0) faces guanine's N1 (index g+6).
     * Initial placement: put C's origin along the G:N1->G:O6 direction,
     * offset by ~2.9 A from G's N1. Ring alignment + azimuthal rotation
     * below will orient the H-bond donors/acceptors correctly. */
    Vec3 g_N1  = sim->atoms[g+6].position;
    Vec3 g_O6  = sim->atoms[g+5].position;
    Vec3 dir   = vec3_normalize(vec3_sub(g_O6, g_N1));
    Vec3 c_pos = vec3_add(g_N1, vec3_scale(dir, 2.9));
    int c = sim_place_cytosine(sim, c_pos);

    /* Align ring normals */
    int g_ring[3] = {g+0, g+1, g+6};
    int c_ring[3] = {c+0, c+1, c+2};
    Vec3 n_G = nb_ring_normal(sim, g_ring);
    Vec3 n_C = nb_ring_normal(sim, c_ring);
    double cosang = vec3_dot(n_G, n_C);
    Vec3 axis = vec3_cross(n_C, n_G);
    double angle;
    if (vec3_norm(axis) < 1e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(n_C, vec3(1,0,0));
        if (vec3_norm(axis) < 1e-8) axis = vec3(0,1,0);
    } else {
        angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    }
    Vec3 c_pivot = sim->atoms[c+0].position;
    nb_transform_rigid(sim, c, 13, c_pivot, axis, angle, vec3_zero());

    /* Azimuthal rotation: align C:N3->N4 with G:N1->O6 */
    g_N1 = sim->atoms[g+6].position;
    g_O6 = sim->atoms[g+5].position;
    Vec3 c_N3 = sim->atoms[c+0].position;
    Vec3 c_N4 = sim->atoms[c+5].position;
    double az = nb_signed_inplane_angle(
        vec3_sub(c_N4, c_N3), vec3_sub(g_O6, g_N1), n_G);
    nb_transform_rigid(sim, c, 13, c_pivot, n_G, az, vec3_zero());

    /* Translate so C:N3 lands at 2.9 A from G:N1 */
    g_N1 = sim->atoms[g+6].position;
    c_N3 = sim->atoms[c+0].position;
    Vec3 to_c = vec3_sub(c_N3, g_N1);
    Vec3 target_c = vec3_add(g_N1, vec3_scale(vec3_normalize(to_c), 2.9));
    nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0,
                       vec3_sub(target_c, c_N3));

    /* ── A-T pair (stacked at z = 3.4 A, B-DNA rise) ──────────────────── */
    int a = sim_place_adenine(sim, vec3(0.0, 0.0, 3.4));
    int t = sim_place_thymine(sim, vec3(2.9, 0.0, 3.4));

    /* Align ring normals */
    int a_ring[3] = {a+0, a+1, a+6};
    int t_ring[3] = {t+0, t+1, t+3};
    Vec3 n_A = nb_ring_normal(sim, a_ring);
    Vec3 n_T = nb_ring_normal(sim, t_ring);
    cosang = vec3_dot(n_A, n_T);
    axis = vec3_cross(n_T, n_A);
    if (vec3_norm(axis) < 1e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(n_T, vec3(1,0,0));
        if (vec3_norm(axis) < 1e-8) axis = vec3(0,1,0);
    } else {
        angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    }
    Vec3 t_pivot = sim->atoms[t+3].position;
    nb_transform_rigid(sim, t, 15, t_pivot, axis, angle, vec3_zero());

    /* Azimuthal rotation: align T:N3->O4 with A:N1->N6 */
    Vec3 a_N1 = sim->atoms[a+6].position;
    Vec3 a_N6 = sim->atoms[a+5].position;
    Vec3 t_N3 = sim->atoms[t+3].position;
    Vec3 t_O4 = sim->atoms[t+5].position;
    double az2 = nb_signed_inplane_angle(
        vec3_sub(t_O4, t_N3), vec3_sub(a_N6, a_N1), n_A);
    nb_transform_rigid(sim, t, 15, t_pivot, n_A, az2, vec3_zero());

    /* Translate so T:N3 lands at 2.9 A from A:N1 */
    a_N1 = sim->atoms[a+6].position;
    t_N3 = sim->atoms[t+3].position;
    Vec3 to_t = vec3_sub(t_N3, a_N1);
    Vec3 target_t = vec3_add(a_N1, vec3_scale(vec3_normalize(to_t), 2.9));
    nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0,
                       vec3_sub(target_t, t_N3));

    printf("  Placed: G-C at z=0, A-T at z=3.4\\n");
    printf("  Total atoms: %d\\n", sim->num_atoms);

    /* Relax initial clashes */
    forces_calculate(sim);
    double pe0 = sim->potential_energy;
    printf("  Initial PE: %.4f eV\\n", pe0);
    if (pe0 > 100.0) {
        double pe_min = integrator_minimize(sim, 5000, 0.001, 0.01);
        printf("  After minimization: PE = %.4f eV\\n", pe_min);
    }

    /* MD at 50 K */
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
    printf("  After %d steps: PE=%.4f eV  T=%.1f K\\n",
           N_steps, sim->potential_energy, sim->temperature);

    /* H-bond distances */
    printf("\\nG-C H-bonds:");
    printf(" N1...N3=%.3f", vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
    printf("  O6...N4=%.3f", vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));
    printf("  N2...O2=%.3f", vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
    printf("\\nA-T H-bonds:");
    printf(" N1...N3=%.3f", vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position));
    printf("  N6...O4=%.3f", vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position));

    /* Planarity */
    int gr[] = {g+0, g+1, g+6};
    int cr[] = {c+0, c+1, c+2};
    int ar[] = {a+0, a+1, a+6};
    int tr[] = {t+0, t+1, t+3};
    printf("\\nPlanarity: G=%.4f  C=%.4f  A=%.4f  T=%.4f\\n",
           nb_planarity_deviation(sim, gr, 3),
           nb_planarity_deviation(sim, cr, 3),
           nb_planarity_deviation(sim, ar, 3),
           nb_planarity_deviation(sim, tr, 3));

    sim_destroy(sim);
}
'''

src = src[:start] + new_func + src[end:]
with open('src/main.c', 'w') as f:
    f.write(src)
print("  [OK] demo_dna_duplex replaced.")
PYEOF

echo "[2/3] Building..."
make clean >/dev/null
make 2>&1 | head -20
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "ERROR: build failed"
    exit 1
fi

echo "[3/3] Running DNA demo..."
./carbonsim --dna

echo "=== s20 complete ==="
