#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/forces.h"
#include "../include/types.h"

/*
 * test_forces.c — analytic dihedral gradient validation (audit P2).
 *
 * Contract: forces_dihedral() (analytic chain rule) must agree with
 * forces_dihedral_fd() (finite-difference oracle) to 1e-6 eV/A on
 * generic, helical, and near-planar geometries, with net-zero total
 * force (translation invariance) to 1e-9, and identical energies.
 * Collinear plane-normals must return finite energy with zero forces
 * (no NaN) under the analytic path.
 */

static int failures = 0;
static void check(int cond, const char *name) {
    if (cond) printf("  PASS  %s\n", name);
    else { printf("  FAIL  %s\n", name); failures++; }
}

static void run_case(const char *name, Vec3 p[4], double k, int n, double delta) {
    Atom aa[4], af[4];
    memset(aa, 0, sizeof aa);
    memset(af, 0, sizeof af);
    for (int i = 0; i < 4; i++) { aa[i].position = p[i]; af[i].position = p[i]; }
    Dihedral dh = {0, 1, 2, 3, k, n, delta};
    double Ea = forces_dihedral(aa, 4, &dh);
    double Ef = forces_dihedral_fd(af, 4, &dh);
    double maxd = 0.0;
    Vec3 sum = vec3_zero();
    for (int i = 0; i < 4; i++) {
        double dx = aa[i].force.x - af[i].force.x;
        double dy = aa[i].force.y - af[i].force.y;
        double dz = aa[i].force.z - af[i].force.z;
        double d = sqrt(dx*dx + dy*dy + dz*dz);
        if (d > maxd) maxd = d;
        vec3_iadd(&sum, aa[i].force);
    }
    double net = vec3_norm(sum);
    char buf[160];
    snprintf(buf, sizeof buf, "%s: analytic-vs-FD max|dF|=%.2e (<1e-6)", name, maxd);
    check(maxd < 1e-6 && isfinite(maxd), buf);
    snprintf(buf, sizeof buf, "%s: net force %.2e (<1e-9)", name, net);
    check(net < 1e-9 && isfinite(net), buf);
    snprintf(buf, sizeof buf, "%s: energies %.6f vs %.6f", name, Ea, Ef);
    check(fabs(Ea - Ef) < 1e-12 && isfinite(Ea), buf);
}

int main(void) {
    printf("forces selftest (dihedral analytic vs FD oracle)\n");
    {
        Vec3 p[4] = {vec3(0,0,0), vec3(1,0,0), vec3(1,1,0), vec3(1,1,1)};
        run_case("orthogonal", p, 0.1, 1, 0.0);
        run_case("orthogonal-n2", p, 0.2, 2, 1.5707963267948966);
        run_case("orthogonal-n3", p, 0.15, 3, 3.141592653589793);
    }
    {
        Vec3 p[4] = {vec3(0,0,0), vec3(1.45,0,0), vec3(2.0,1.2,0.3), vec3(1.2,2.0,-0.4)};
        run_case("helical", p, 0.1, 1, 0.5);
        run_case("helical-n2", p, 0.3, 2, 0.0);
    }
    {
        Vec3 p[4] = {vec3(0,0,0), vec3(1.4,0,0), vec3(2.8,0.05,0), vec3(4.2,0.02,0.5)};
        run_case("near-planar", p, 0.1, 1, 0.0);
    }
    {
        Vec3 p[4] = {vec3(-1.2,0.3,0.1), vec3(0,0,0), vec3(1.3,0.2,-0.1), vec3(2.5,-0.4,0.3)};
        run_case("generic", p, 0.25, 1, 2.0);
    }
    /* Collinear guard: energy finite, forces zero, no NaN. */
    {
        Atom ac[4];
        memset(ac, 0, sizeof ac);
        ac[0].position = vec3(0,0,0); ac[1].position = vec3(1,0,0);
        ac[2].position = vec3(2,0,0); ac[3].position = vec3(3,0,0);
        Dihedral dh = {0, 1, 2, 3, 0.1, 1, 0.0};
        double E = forces_dihedral(ac, 4, &dh);
        double f = 0.0;
        for (int i = 0; i < 4; i++) f += vec3_norm(ac[i].force);
        check(isfinite(E) && f == 0.0, "collinear: finite energy, zero forces");
    }
    if (failures) { printf("forces selftest: %d FAILURE(S)\n", failures); return 1; }
    printf("forces selftest: all passed\n");
    return 0;
}
