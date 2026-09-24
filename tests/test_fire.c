#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/sim.h"
#include "../include/forces.h"
#include "../include/integrator.h"
#include "../include/periodic_table.h"

/*
 * test_fire.c — FIRE minimizer validation.
 *
 * Contract: on a fixed 7-atom Lennard-Jones cluster (argon-like,
 * deterministic displaced-cubic start), FIRE and steepest descent must
 * reach the same minimum to 5e-3 eV, both downhill from the start, with
 * finite energies and intact topology. Velocities and dt must be
 * restored by FIRE (MD state untouched).
 */

static int failures = 0;
static void check(int cond, const char *name) {
    if (cond) printf("  PASS  %s\n", name);
    else { printf("  FAIL  %s\n", name); failures++; }
}

static Simulation *make_cluster(void) {
    Simulation *sim = sim_create(16, 16);
    if (!sim) return NULL;
    /* Deterministic perturbed cube + center (argon LJ via UFF table). */
    static const double xyz[7][3] = {
        {0.0, 0.0, 0.0},
        {3.4, 0.1, -0.1}, {-0.1, 3.5, 0.1}, {0.1, -0.1, 3.3},
        {3.5, 3.4, 0.2}, {3.3, 0.2, 3.6}, {0.0, 3.6, 3.4}
    };
    for (int i = 0; i < 7; i++)
        sim_add_atom(sim, 18, vec3(xyz[i][0], xyz[i][1], xyz[i][2]), 0.0);
    sim->cutoff = 100.0;
    sim->dielectric = 1.0;
    return sim;
}

int main(void) {
    printf("fire selftest (vs steepest descent oracle)\n");
    Simulation *a = make_cluster();
    Simulation *b = make_cluster();
    check(a && b, "cluster allocation");
    if (!a || !b) return 1;
    forces_calculate(a);
    double e0 = a->potential_energy;
    /* Give both sims identical nonzero velocities to check restore. */
    for (int i = 0; i < 7; i++) {
        a->atoms[i].velocity = vec3(0.01 * (i + 1), -0.005 * i, 0.002 * i);
        b->atoms[i].velocity = a->atoms[i].velocity;
    }
    Vec3 v0 = a->atoms[3].velocity;
    double dt0 = 0.5;
    a->dt = dt0; b->dt = dt0;

    double es = integrator_minimize(a, 2000, 0.01, 1e-4);
    double ef = integrator_fire(b, 2000, 0.5, 1e-4);
    printf("  steepest: %.6f eV   FIRE: %.6f eV   (start %.6f eV)\n", es, ef, e0);
    check(isfinite(es) && isfinite(ef), "both finite");
    check(es < e0 && ef < e0, "both downhill");
    check(fabs(es - ef) < 5e-3, "minima agree to 5e-3 eV");
    check(fabs(b->atoms[3].velocity.x - v0.x) < 1e-12 &&
          fabs(b->atoms[3].velocity.y - v0.y) < 1e-12 &&
          fabs(b->atoms[3].velocity.z - v0.z) < 1e-12,
          "FIRE restores velocities");
    check(fabs(b->dt - dt0) < 1e-12, "FIRE restores dt");
    /* Topology intact: no atom teleported (all within 20 A of origin). */
    {
        int ok = 1;
        for (int i = 0; i < 7; i++)
            if (vec3_norm(b->atoms[i].position) > 20.0) ok = 0;
        check(ok, "FIRE cluster intact (no tunneling)");
    }
    sim_destroy(a);
    sim_destroy(b);
    if (failures) { printf("fire selftest: %d FAILURE(S)\n", failures); return 1; }
    printf("fire selftest: all passed\n");
    return 0;
}
