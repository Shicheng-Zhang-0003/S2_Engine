#!/usr/bin/env bash
set -euo pipefail
# s30: fix the primary bug — pyrimidines placed on wrong side of purines.
# C:N3 must be on the OPPOSITE side of G:N1 from G:O6.
# T:N3 must be on the OPPOSITE side of A:N1 from A:N6.
# Also: place C and T far away first, stack A-T along Y, coplanarity fix.
echo "=== s30: fix wrong-side placement (the primary bug) ==="

python3 - <<'PYEOF'
import sys

with open('src/main.c') as f:
    src = f.read()

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
double rise = 3.4;

/* ══════════════════════════════════════════════════════════════════
* G-C PAIR
* C:N3 must end up on the OPPOSITE side of G:N1 from G:O6.
* ══════════════════════════════════════════════════════════════════ */
int g = sim_place_guanine(sim, vec3_zero());
int c = sim_place_cytosine(sim, vec3(15.0, 0.0, 0.0));  /* far away */

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

/* Coplanarity: translate C in Y to match G's ring plane */
{
    Vec3 gc = vec3_zero(), cc = vec3_zero();
    for (int i = 0; i < 3; i++) {
        gc = vec3_add(gc, sim->atoms[g_ring[i]].position);
        cc = vec3_add(cc, sim->atoms[c_ring[i]].position);
    }
    gc = vec3_scale(gc, 1.0/3.0);
    cc = vec3_scale(cc, 1.0/3.0);
    nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0,
                       vec3(0.0, gc.y - cc.y, 0.0));
}

/* Azimuthal: align C:N3->C:N4 with G:N1->G:O6 */
Vec3 g_N1 = sim->atoms[g+6].position;
Vec3 g_O6 = sim->atoms[g+5].position;
Vec3 c_N3 = sim->atoms[c+0].position;
Vec3 c_N4 = sim->atoms[c+5].position;
double az = nb_signed_inplane_angle(vec3_sub(c_N4, c_N3),
                                    vec3_sub(g_O6, g_N1), n_G);
nb_transform_rigid(sim, c, 13, c_pivot, n_G, az, vec3_zero());

/* Translate C:N3 to 2.9 A from G:N1, OPPOSITE side from G:O6 */
g_N1 = sim->atoms[g+6].position;
g_O6 = sim->atoms[g+5].position;
c_N3 = sim->atoms[c+0].position;
Vec3 dir_gc = vec3_normalize(vec3_sub(g_N1, g_O6));   /* away from O6 */
Vec3 tgt_gc = vec3_add(g_N1, vec3_scale(dir_gc, 2.9));
nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0,
                   vec3_sub(tgt_gc, c_N3));

/* ══════════════════════════════════════════════════════════════════
* A-T PAIR  (stacked along Y = ring-normal direction)
* T:N3 must end up on the OPPOSITE side of A:N1 from A:N6.
* ══════════════════════════════════════════════════════════════════ */
int a = sim_place_adenine(sim, vec3(0.0, rise, 0.0));
int t = sim_place_thymine(sim, vec3(15.0, rise, 0.0));  /* far away */

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

/* Coplanarity for A-T */
{
    Vec3 ac = vec3_zero(), tc = vec3_zero();
    for (int i = 0; i < 3; i++) {
        ac = vec3_add(ac, sim->atoms[a_ring[i]].position);
        tc = vec3_add(tc, sim->atoms[t_ring[i]].position);
    }
    ac = vec3_scale(ac, 1.0/3.0);
    tc = vec3_scale(tc, 1.0/3.0);
    nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0,
                       vec3(0.0, ac.y - tc.y, 0.0));
}

/* Azimuthal: align T:N3->T:O4 with A:N1->A:N6 */
Vec3 a_N1 = sim->atoms[a+6].position;
Vec3 a_N6 = sim->atoms[a+5].position;
Vec3 t_N3 = sim->atoms[t+3].position;
Vec3 t_O4 = sim->atoms[t+5].position;
double az2 = nb_signed_inplane_angle(vec3_sub(t_O4, t_N3),
                                     vec3_sub(a_N6, a_N1), n_A);
nb_transform_rigid(sim, t, 15, t_pivot, n_A, az2, vec3_zero());

/* Translate T:N3 to 2.9 A from A:N1, OPPOSITE side from A:N6 */
a_N1 = sim->atoms[a+6].position;
a_N6 = sim->atoms[a+5].position;
t_N3 = sim->atoms[t+3].position;
Vec3 dir_at = vec3_normalize(vec3_sub(a_N1, a_N6));   /* away from N6 */
Vec3 tgt_at = vec3_add(a_N1, vec3_scale(dir_at, 2.9));
nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0,
                   vec3_sub(tgt_at, t_N3));

/* ── Report ─────────────────────────────────────────────────────── */
printf("  Placed: G-C at y=0, A-T at y=%.1f\n", rise);
printf("  Total atoms: %d\n", sim->num_atoms);

forces_calculate(sim);
double initial_pe = sim->potential_energy;
printf("  Initial PE: %.6f eV\n", initial_pe);

double min_dist = 1.0e9; int ci = -1, cj = -1;
for (int i = 0; i < sim->num_atoms - 1; i++)
    for (int j = i + 1; j < sim->num_atoms; j++) {
        double d = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
        if (d < min_dist) { min_dist = d; ci = i; cj = j; }
    }
printf("  Closest pair: %d-%d at %.3f A\n", ci, cj, min_dist);
if (min_dist < 1.0) printf("  WARNING: steric clash (< 1.0 A)\n");

if (initial_pe > 100.0) {
    double mp = integrator_minimize(sim, 5000, 0.001, 0.01);
    printf("  After minimization: PE = %.6f eV\n", mp);
}

printf("  Post-min H-bonds:\n");
printf("    G-C: N1...N3=%.3f  O6...N4=%.3f  N2...O2=%.3f\n",
       vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position),
       vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position),
       vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
printf("    A-T: N1...N3=%.3f  N6...O4=%.3f\n",
       vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position),
       vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position));

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

printf("\nG-C H-bonds:");
printf(" N1...N3=%.3f", vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
printf("  O6...N4=%.3f", vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));
printf("  N2...O2=%.3f", vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
printf("\nA-T H-bonds:");
printf(" N1...N3=%.3f", vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position));
printf("  N6...O4=%.3f", vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position));

int gr[]={g+0,g+1,g+6}, cr[]={c+0,c+1,c+2};
int ar[]={a+0,a+1,a+6}, tr[]={t+0,t+1,t+3};
printf("\nPlanarity: G=%.4f  C=%.4f  A=%.4f  T=%.4f\n",
       nb_planarity_deviation(sim,gr,3), nb_planarity_deviation(sim,cr,3),
       nb_planarity_deviation(sim,ar,3), nb_planarity_deviation(sim,tr,3));

double gc1=vec3_dist(sim->atoms[g+6].position,sim->atoms[c+0].position);
double gc2=vec3_dist(sim->atoms[g+5].position,sim->atoms[c+5].position);
double gc3=vec3_dist(sim->atoms[g+8].position,sim->atoms[c+4].position);
double at1=vec3_dist(sim->atoms[a+6].position,sim->atoms[t+3].position);
double at2=vec3_dist(sim->atoms[a+5].position,sim->atoms[t+5].position);
int gc_ok=(gc1<3.3)&&(gc2<3.3)&&(gc3<3.3);
int at_ok=(at1<3.3)&&(at2<3.3);
printf("\nVERDICT: G-C %s, A-T %s\n",
       gc_ok?"PAIRED":"NOT PAIRED", at_ok?"PAIRED":"NOT PAIRED");
if (gc_ok && at_ok)
    printf("--> A real, H-bonded two-base-pair DNA stack.\n");
else
    printf("--> Bases drifted apart or sheared; NOT a true duplex.\n");

sim_destroy(sim);
}
'''

src = src[:start] + new_func + src[end:]
with open('src/main.c', 'w') as f:
    f.write(src)
print('  [OK] demo_dna_duplex replaced with correct-side placement.')
PYEOF

echo ""
echo "Building..."
make clean >/dev/null
make 2>&1 | tail -3

echo ""
echo "Running DNA demo..."
./carbonsim --dna

echo "=== s30 complete ==="
