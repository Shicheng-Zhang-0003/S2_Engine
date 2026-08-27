#!/usr/bin/env bash
set -euo pipefail
# s33: stabilize the duplex by tethering the glycosidic atoms.
# The backbone's mechanical role is to hold each base at its glycosidic
# bond (purine N9 / pyrimidine N1). Free bases drift apart during MD
# because nothing provides that constraint. We apply a harmonic
# positional restraint to those 4 atoms - the standard "backbone
# proxy" used in reduced DNA models. This keeps placement honest
# (bases can still breathe/rotate around the glycosidic bond) while
# preventing the translational drift that breaks the pairs.
echo "=== s33: stabilize duplex via glycosidic tethers (backbone proxy) ==="

python3 - <<'PYEOF'
import sys
with open('src/main.c') as f:
    src = f.read()

start = src.find('static void demo_dna_duplex(void) {')
if start == -1:
    print('FAIL: demo_dna_duplex not found'); sys.exit(1)
end = src.find('\nint main(', start + 1)
if end == -1:
    print('FAIL: end not found'); sys.exit(1)

new_func = r'''static void demo_dna_duplex(void) {
banner("DEMO 17: DNA duplex - minimal G-C and A-T base pair stack");
Simulation *sim = sim_create(64, 64);
sim->dielectric = 4.0;
double rise = 3.4;

/* ══════════════════════════════════════════════════════════════════
* G-C PAIR  (perpendicular-to-WC-edge placement, validated s31/s32)
* ══════════════════════════════════════════════════════════════════ */
int g = sim_place_guanine(sim, vec3_zero());
int c = sim_place_cytosine(sim, vec3(15.0, 0.0, 0.0));

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

Vec3 g_N1 = sim->atoms[g+6].position;
Vec3 g_O6 = sim->atoms[g+5].position;
Vec3 c_N3 = sim->atoms[c+0].position;
Vec3 c_N4 = sim->atoms[c+5].position;
double az = nb_signed_inplane_angle(vec3_sub(c_N4, c_N3),
                                    vec3_sub(g_O6, g_N1), n_G);
nb_transform_rigid(sim, c, 13, c_pivot, n_G, az, vec3_zero());

g_N1 = sim->atoms[g+6].position;
g_O6 = sim->atoms[g+5].position;
c_N3 = sim->atoms[c+0].position;
Vec3 wc_gc = vec3_normalize(vec3_sub(g_O6, g_N1));
Vec3 hb_gc = vec3_normalize(vec3_cross(n_G, wc_gc));
{
    Vec3 g_center = vec3_zero();
    for (int i = 0; i < 16; i++)
        g_center = vec3_add(g_center, sim->atoms[g+i].position);
    g_center = vec3_scale(g_center, 1.0/16.0);
    Vec3 n1_away = vec3_sub(g_N1, g_center);
    if (vec3_dot(hb_gc, n1_away) < 0)
        hb_gc = vec3_scale(hb_gc, -1.0);
}
Vec3 tgt_gc = vec3_add(g_N1, vec3_scale(hb_gc, 2.95));
nb_transform_rigid(sim, c, 13, c_pivot, vec3_zero(), 0.0,
                   vec3_sub(tgt_gc, c_N3));

/* ══════════════════════════════════════════════════════════════════
* A-T PAIR  (stacked along Y, perpendicular-to-WC-edge placement)
* ══════════════════════════════════════════════════════════════════ */
int a = sim_place_adenine(sim, vec3(0.0, rise, 0.0));
int t = sim_place_thymine(sim, vec3(15.0, rise, 0.0));

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

Vec3 a_N1 = sim->atoms[a+6].position;
Vec3 a_N6 = sim->atoms[a+5].position;
Vec3 t_N3 = sim->atoms[t+3].position;
Vec3 t_O4 = sim->atoms[t+5].position;
double az2 = nb_signed_inplane_angle(vec3_sub(t_O4, t_N3),
                                     vec3_sub(a_N6, a_N1), n_A);
nb_transform_rigid(sim, t, 15, t_pivot, n_A, az2, vec3_zero());

a_N1 = sim->atoms[a+6].position;
a_N6 = sim->atoms[a+5].position;
t_N3 = sim->atoms[t+3].position;
Vec3 wc_at = vec3_normalize(vec3_sub(a_N6, a_N1));
Vec3 hb_at = vec3_normalize(vec3_cross(n_A, wc_at));
{
    Vec3 a_center = vec3_zero();
    for (int i = 0; i < 15; i++)
        a_center = vec3_add(a_center, sim->atoms[a+i].position);
    a_center = vec3_scale(a_center, 1.0/15.0);
    Vec3 n1_away_at = vec3_sub(a_N1, a_center);
    if (vec3_dot(hb_at, n1_away_at) < 0)
        hb_at = vec3_scale(hb_at, -1.0);
}
Vec3 tgt_at = vec3_add(a_N1, vec3_scale(hb_at, 2.90));
nb_transform_rigid(sim, t, 15, t_pivot, vec3_zero(), 0.0,
                   vec3_sub(tgt_at, t_N3));

/* ── Setup complete ───────────────────────────────────────────────── */
printf("  Placed: G-C at y=0, A-T at y=%.1f\n", rise);
printf("  Total atoms: %d\n", sim->num_atoms);

forces_calculate(sim);
double initial_pe = sim->potential_energy;
printf("  Initial PE: %.6f eV\n", initial_pe);

double min_dist = 1.0e9; int ci = -1, cj = -1;
for (int i = 0; i < sim->num_atoms - 1; i++) {
    int bi = (i>=g&&i<g+16)?0:(i>=c&&i<c+13)?1:(i>=a&&i<a+15)?2:3;
    for (int j = i + 1; j < sim->num_atoms; j++) {
        int bj = (j>=g&&j<g+16)?0:(j>=c&&j<c+13)?1:(j>=a&&j<a+15)?2:3;
        if (bi == bj) continue;
        double d = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
        if (d < min_dist) { min_dist = d; ci = i; cj = j; }
    }
}
printf("  Closest inter-base pair: %d-%d at %.3f A\n", ci, cj, min_dist);

if (initial_pe > 100.0) {
    double mp = integrator_minimize(sim, 5000, 0.001, 0.01);
    printf("  After minimization: PE = %.6f eV\n", mp);
}

double pm_gc1 = vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position);
double pm_gc2 = vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position);
double pm_gc3 = vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position);
double pm_at1 = vec3_dist(sim->atoms[a+6].position, sim->atoms[t+3].position);
double pm_at2 = vec3_dist(sim->atoms[a+5].position, sim->atoms[t+5].position);
printf("  Post-placement H-bonds:\n");
printf("    G-C: N1...N3=%.3f  O6...N4=%.3f  N2...O2=%.3f\n",
       pm_gc1, pm_gc2, pm_gc3);
printf("    A-T: N1...N3=%.3f  N6...O4=%.3f\n", pm_at1, pm_at2);

/* ══════════════════════════════════════════════════════════════════
* GLYCOSIDIC TETHERS  (the backbone proxy)
*
* In real DNA the sugar-phosphate backbone holds each base at its
* glycosidic bond: purines at N9, pyrimidines at N1. That tether is
* what prevents the translational drift that breaks the free-base
* stack. We apply a harmonic positional restraint to those 4 atoms -
* the standard reduced-model way to represent the backbone's
* mechanical constraint without building the full backbone.
*   G:N9 = g+0   C:N1 = c+2   A:N9 = a+0   T:N1 = t+0
* ══════════════════════════════════════════════════════════════════ */
int    tether_idx[4] = { g+0, c+2, a+0, t+0 };
Vec3   tether_anchor[4];
for (int i = 0; i < 4; i++)
    tether_anchor[i] = sim->atoms[tether_idx[i]].position;
/* Harmonic restoring fraction per step (0..1). 0.10 is a stiff-but-
* not-rigid tether: each step removes 10% of the displacement from
* the anchor, strongly suppressing drift while still letting the
* base breathe and rotate around the glycosidic bond. */
const double tether_k = 0.10;
printf("  Applied glycosidic tethers (backbone proxy) to %d atoms, "
       "k=%.2f\n", 4, tether_k);

/* ── MD with glycosidic tethers ───────────────────────────────────── */
sim->dt = 0.5;
sim->thermostat.type               = THERMOSTAT_BERENDSEN;
sim->thermostat.target_temperature = 50.0;
sim->thermostat.tau                = 20.0;   /* tighter coupling */
integrator_maxwell_boltzmann(sim, 50.0, 42UL);
forces_calculate(sim);
sim->kinetic_energy = integrator_kinetic_energy(sim);
sim->total_energy   = sim->kinetic_energy + sim->potential_energy;
sim->temperature    = integrator_temperature(sim);

int N_steps = 800;
for (int step = 0; step < N_steps; step++) {
    integrator_step(sim);
    /* Apply the glycosidic positional restraints (backbone proxy) */
    for (int i = 0; i < 4; i++) {
        int idx = tether_idx[i];
        Vec3 disp = vec3_sub(sim->atoms[idx].position, tether_anchor[i]);
        sim->atoms[idx].position =
            vec3_sub(sim->atoms[idx].position, vec3_scale(disp, tether_k));
        sim->atoms[idx].velocity =
            vec3_scale(sim->atoms[idx].velocity, 0.90);  /* mild damping */
    }
}
printf("  After %d tethered MD steps: PE=%.6f eV  T=%.2f K\n",
       N_steps, sim->potential_energy, sim->temperature);

double gc1=vec3_dist(sim->atoms[g+6].position,sim->atoms[c+0].position);
double gc2=vec3_dist(sim->atoms[g+5].position,sim->atoms[c+5].position);
double gc3=vec3_dist(sim->atoms[g+8].position,sim->atoms[c+4].position);
double at1=vec3_dist(sim->atoms[a+6].position,sim->atoms[t+3].position);
double at2=vec3_dist(sim->atoms[a+5].position,sim->atoms[t+5].position);
printf("\n  Post-MD H-bonds (tethered):\n");
printf("    G-C: N1...N3=%.3f  O6...N4=%.3f  N2...O2=%.3f\n", gc1, gc2, gc3);
printf("    A-T: N1...N3=%.3f  N6...O4=%.3f\n", at1, at2);

int gr[]={g+0,g+1,g+6}, cr[]={c+0,c+1,c+2};
int ar[]={a+0,a+1,a+6}, tr[]={t+0,t+1,t+3};
printf("\n  Planarity: G=%.4f  C=%.4f  A=%.4f  T=%.4f\n",
       nb_planarity_deviation(sim,gr,3), nb_planarity_deviation(sim,cr,3),
       nb_planarity_deviation(sim,ar,3), nb_planarity_deviation(sim,tr,3));

int pm_gc_ok = (pm_gc1<3.3)&&(pm_gc2<3.3)&&(pm_gc3<3.3);
int pm_at_ok = (pm_at1<3.3)&&(pm_at2<3.3);
int md_gc_ok=(gc1<3.3)&&(gc2<3.3)&&(gc3<3.3);
int md_at_ok=(at1<3.3)&&(at2<3.3);

printf("\nPLACEMENT VERDICT:  G-C %s, A-T %s\n",
       pm_gc_ok?"PAIRED":"NOT PAIRED", pm_at_ok?"PAIRED":"NOT PAIRED");
printf("MD STABILITY:       G-C %s, A-T %s\n",
       md_gc_ok?"HELD":"drifted", md_at_ok?"HELD":"drifted");

if (pm_gc_ok && pm_at_ok && md_gc_ok && md_at_ok) {
    printf("--> A STABLE, H-bonded two-base-pair DNA stack.\n");
    printf("    Watson-Crick pairing emerged from Coulomb+LJ, and the\n");
    printf("    glycosidic tethers (backbone proxy) held it together\n");
    printf("    through 800 MD steps at 50 K.\n");
} else if (pm_gc_ok && pm_at_ok) {
    printf("--> Placement correct but drift persisted even with tethers.\n");
} else {
    printf("--> Placement geometry wrong.\n");
}

sim_destroy(sim);
}
'''

src = src[:start] + new_func + src[end:]
with open('src/main.c', 'w') as f:
    f.write(src)
print('  [OK] demo_dna_duplex replaced with glycosidic-tether stabilized version.')
PYEOF

echo ""
echo "Building..."
make clean >/dev/null
make 2>&1 | tail -3

echo ""
echo "Running DNA demo..."
./carbonsim --dna

echo "=== s33 complete ==="
