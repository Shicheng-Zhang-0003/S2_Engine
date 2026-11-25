/*
 * TEMPORARY DIAGNOSTIC, NOT PART OF THE ENGINE BUILD OR THE RECORD.
 *
 * Drives the engine's own KcsA sampler (qm_kcsa_binding_series) after the
 * side-fix parameters have been compiled into a throwaway copy of
 * src/qm_fluctuation.c. The copy is made by run_kcsa_engine_probe.sh in a
 * temp directory; no sealed file is modified.
 *
 * Protocol: real 1K4C TVGYG tetramer + side-chain hydrogens, S2 site,
 * 300 K, K+ vs Na+, 200 equilibration steps + N sampled steps, Langevin.
 * dt, gamma and epsilon are compile-time overrides of the copied source:
 *
 *   QM_KCSA_DT     default 1.0 fs   (side fix uses 0.5)
 *   QM_KCSA_GAMMA  default 0.001    (side fix uses 0.02)
 *   QM_KCSA_EPS    default 1.0      (side fix uses 4.0)
 *
 * The tether stiffness is a runtime argument.
 *
 * Sign convention follows include/qm_fluctuation.h:
 *   ddG = 0.5*beta*(var_K - var_Na), NEGATIVE = K+ favoured.
 * The side-fix module originally returned the opposite sign; the probe
 * prints both so the difference cannot be missed.
 */
#include "qm_fluctuation.h"
#include "kcsa_filter.h"
#include "constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void show(const char *tag, int Z, Vec3 site, double tether,
                 int n_equil, int n_sample, unsigned long seed,
                 QmKcsaSeries *o) {
    int rc = qm_kcsa_binding_series(Z, site, 300.0, n_equil, n_sample, seed,
                                    tether, NULL, NULL, 0, o);
    printf("%s rc=%2d ok=%d held=%d n=%3d  CN=%6.3f  <ion-O>=%.3f  <O-O>=%.3f  "
           "<Ebind>=%9.4f  var=%9.6f  sem=%8.6f  semvar=%9.6f  "
           "<Efilt>=%9.4f  varF=%9.6f  ooLJ=%9.6f\n",
           tag, rc, o->ok, o->site_maintained, o->n_samples,
           o->cn_mean, o->mean_iono_A, o->mean_oo_A,
           o->mean_bind_eV, o->var_bind_eV2, o->sem_bind_eV,
           o->sem_var_bind_eV2, o->mean_filter_eV, o->var_filter_eV2,
           o->oo_lj_eV);
}

int main(int argc, char **argv) {
    double tether = (argc > 1) ? atof(argv[1]) : 8.0;
    int n_equil   = (argc > 2) ? atoi(argv[2]) : 200;
    int n_sample  = (argc > 3) ? atoi(argv[3]) : 200;

    if (!(tether > 0.0) || n_equil < 0 || n_sample < 2) {
        fprintf(stderr, "usage: %s [tether_k] [n_equil] [n_sample]\n", argv[0]);
        return 2;
    }

    Vec3 sites[4];
    if (kcsa_ion_sites(4, sites) < 2) {
        fprintf(stderr, "no ion sites\n");
        return 1;
    }

    QmKcsaSeries K, Na;
    show("K+", 19, sites[1], tether, n_equil, n_sample, 424242UL, &K);
    show("Na+", 11, sites[1], tether, n_equil, n_sample, 848484UL, &Na);

    double kBT = BOLTZMANN_K * 300.0 * J_TO_EV;
    double beta = 1.0 / kBT;
    double d_engine = 0.5 * beta * (K.var_bind_eV2 - Na.var_bind_eV2);
    double d_sidefix = 0.5 * beta * (Na.var_bind_eV2 - K.var_bind_eV2);
    double err = 0.5 * beta * sqrt(K.sem_var_bind_eV2 * K.sem_var_bind_eV2 +
                                   Na.sem_var_bind_eV2 * Na.sem_var_bind_eV2);
    printf("ddG engine sign   = %+.4f +/- %.4f eV = %+.2f kcal/mol "
           "(0.5*beta*(var_K-var_Na); NEGATIVE = K+ favoured)\n",
           d_engine, err, d_engine / KCAL_MOL_TO_EV);
    printf("ddG side-fix sign = %+.4f +/- %.4f eV = %+.2f kcal/mol "
           "(0.5*beta*(var_Na-var_K); OPPOSITE of engine convention)\n",
           d_sidefix, err, d_sidefix / KCAL_MOL_TO_EV);
    printf("literature target = 5-6 kcal/mol = %.3f-%.3f eV (Dixit-Asthagiri)\n",
           5.0 * KCAL_MOL_TO_EV, 6.0 * KCAL_MOL_TO_EV);
    printf("beta = %.6f 1/eV, kBT = %.6f eV\n", beta, kBT);
    return 0;
}
