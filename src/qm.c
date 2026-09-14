#include <math.h>
#include <string.h>
#include "../include/qm.h"
#include "../include/quantum.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Real spherical harmonics ──────────────────────────────────────── */

double qm_Y_real(int l, int ml, double theta, double phi) {
    double ct = cos(theta), st = sin(theta);
    const double s_norm = 0.28209479177387814;   /* sqrt(1/4pi) */
    const double p_norm = 0.4886025119029199;    /* sqrt(3/4pi) */
    const double d_norm_a = 0.31539156525252005; /* sqrt(5/16pi), z2 + x2-y2/xy set */
    const double d_norm_b = 1.0925484305920792;  /* sqrt(15/4pi) for xz/yz */
    const double d_norm_c = 0.5462742152960396;  /* sqrt(15/16pi) for x2-y2/xy */
    if (l == 0) { (void)ml; (void)theta; (void)phi; return s_norm; }
    if (l == 1) {
        if (ml == 0) return p_norm * ct;                    /* pz */
        if (ml == 1) return p_norm * st * cos(phi);         /* px */
        if (ml == -1) return p_norm * st * sin(phi);        /* py */
        return 0.0;
    }
    if (l == 2) {
        double st2 = st * st;
        if (ml == 0) return d_norm_a * (3.0 * ct * ct - 1.0);          /* dz2 */
        if (ml == 1) return d_norm_b * st * ct * cos(phi);             /* dxz */
        if (ml == -1) return d_norm_b * st * ct * sin(phi);            /* dyz */
        if (ml == 2) return d_norm_c * st2 * cos(2.0 * phi);           /* dx2-y2 */
        if (ml == -2) return d_norm_c * st2 * sin(2.0 * phi);          /* dxy */
        return 0.0;
    }
    return 0.0; /* f and above: documented cutoff in v1 */
}

double qm_psi(int n, int l, int ml, double Zeff, Vec3 pos) {
    double r = vec3_norm(pos);
    double R = quantum_radial_wavefunction(n, l, Zeff, r);
    double theta = 0.0, phi = 0.0;
    if (r > 1e-12) {
        theta = acos(pos.z / r);
        phi = atan2(pos.y, pos.x);
    }
    return R * qm_Y_real(l, ml, theta, phi);
}

/* ── Hybridization ─────────────────────────────────────────────────── */

QmHybrid qm_hybridization(const Atom *atom) {
    QmHybrid h;
    memset(&h, 0, sizeof h);
    int s_count = atom->electron_config.config[1][0]; /* 2s census */
    int p_count = atom->electron_config.config[1][1]; /* 2p census */
    int n_bonds = atom->num_bonds;
    if (atom->Z == 1) {
        snprintf(h.label, sizeof h.label, "1s");
        h.n_lobes = (n_bonds > 0) ? 1 : 0;
        if (h.n_lobes) h.lobes[0] = vec3(0, 0, 1);
        h.n_lone_pairs = 0;
        return h;
    }
    /* Free atom (no bonds): hybridization is a molecular concept; report
     * the atomic valence shell instead of forcing sp labels. */
    if (n_bonds == 0) {
        int has_s = atom->electron_config.config[1][0] > 0;
        int has_p = atom->electron_config.config[1][1] > 0;
        if (has_s && has_p) snprintf(h.label, sizeof h.label, "atomic-sp");
        else if (has_p) snprintf(h.label, sizeof h.label, "atomic-p");
        else snprintf(h.label, sizeof h.label, "atomic-s");
        h.n_lobes = 0;
        int v = atom->electron_config.valence_electrons;
        h.n_lone_pairs = (v > 0) ? v / 2 : 0;
        return h;
    }
    /* Second-period logic from s/p census + coordination. Carbon with
     * 4 sigma partners -> sp3; 3 -> sp2; 2 -> sp; lone pairs fill the
     * remainder of the 4 tetrahedral slots. Generic for N/O as well. */
    int steric = n_bonds;
    /* Count lone pairs from valence octet: 8 - (bonding e + nonbonding) */
    int v = atom->electron_config.valence_electrons;
    int lone = 0;
    if (v > 0) {
        int bonding_e = n_bonds; /* ~1 e contributed per sigma bond */
        int rem = v - bonding_e;
        lone = (rem > 0) ? rem / 2 : 0;
        steric = n_bonds + lone;
    }
    if (steric >= 4 || (s_count > 0 && p_count >= 3)) {
        snprintf(h.label, sizeof h.label, "sp3");
        h.n_lobes = 4;
        h.lobes[0] = vec3_normalize(vec3(1, 1, 1));
        h.lobes[1] = vec3_normalize(vec3(1, -1, -1));
        h.lobes[2] = vec3_normalize(vec3(-1, 1, -1));
        h.lobes[3] = vec3_normalize(vec3(-1, -1, 1));
    } else if (steric == 3 || p_count == 2) {
        snprintf(h.label, sizeof h.label, "sp2");
        h.n_lobes = 3;
        for (int i = 0; i < 3; i++) {
            double a = i * 2.0 * M_PI / 3.0;
            h.lobes[i] = vec3(cos(a), sin(a), 0);
        }
    } else if (steric == 2) {
        snprintf(h.label, sizeof h.label, "sp");
        h.n_lobes = 2;
        h.lobes[0] = vec3(1, 0, 0);
        h.lobes[1] = vec3(-1, 0, 0);
    } else if (p_count > 0) {
        snprintf(h.label, sizeof h.label, "p");
        h.n_lobes = (n_bonds > 0) ? n_bonds : 1;
        for (int i = 0; i < h.n_lobes && i < 4; i++)
            h.lobes[i] = vec3(0, 0, 1);
    } else {
        snprintf(h.label, sizeof h.label, "s");
        h.n_lobes = 0;
    }
    h.n_lone_pairs = lone;
    (void)s_count;
    return h;
}

/* ── chi/J and alpha ───────────────────────────────────────────────── */

void qm_chi_J(const Element *el, double *chi_out, double *J_out) {
    double ie = el->ionization_energy;
    double ea = el->electron_affinity; /* table uses 0 for negative */
    if (ea < 0) ea = 0;
    if (chi_out) *chi_out = 0.5 * (ie + ea);
    if (J_out) *J_out = 0.5 * (ie - ea);
}

double qm_alpha(const Atom *atom) {
    ElectronConfig cfg = atom->electron_config;
    /* Valence (n,l) with highest (energy, n, l) — same rule as Demo 1. */
    int bn = 1, bl = 0, found = 0;
    double be = -1e300;
    static const int MN[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
    static const int ML[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
    for (int i = 0; i < 19; i++) {
        int n = MN[i], l = ML[i];
        if (cfg.config[n-1][l] == 0) continue;
        double e = quantum_orbital_energy(atom->Z, n, l, &cfg);
        if (!found || e > be || (e == be && (n > bn || (n == bn && l > bl)))) {
            found = 1; be = e; bn = n; bl = l;
        }
    }
    if (!found) return 0.0;
    double zeff = quantum_zeff(atom->Z, bn, bl, &cfg);
    double rmp = quantum_most_probable_radius(bn, bl, zeff);
    return rmp * rmp * rmp; /* sphere-volume order of magnitude */
}

/* ── Overlap / bond order / Pauli ──────────────────────────────────── */

double qm_overlap(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b) {
    if (R_ang < 1e-9) return 1.0;
    /* Effective decay from each atom's valence shell. */
    ElectronConfig ca = a->electron_config, cb = b->electron_config;
    int na = 1, la = 0, nb = 1, lb = 0;
    double ea = -1e300, eb = -1e300;
    int fa = 0, fb = 0;
    static const int MN[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
    static const int ML[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
    for (int i = 0; i < 19; i++) {
        int n = MN[i], l = ML[i];
        if (ca.config[n-1][l]) {
            double e = quantum_orbital_energy(a->Z, n, l, &ca);
            if (!fa || e > ea) { fa = 1; ea = e; na = n; la = l; }
        }
        if (cb.config[n-1][l]) {
            double e = quantum_orbital_energy(b->Z, n, l, &cb);
            if (!fb || e > eb) { fb = 1; eb = e; nb = n; lb = l; }
        }
    }
    double za = quantum_zeff(a->Z, na, la, &ca);
    double zb = quantum_zeff(b->Z, nb, lb, &cb);
    double a0 = BOHR_TO_ANGSTROM;
    double zeta = 0.5 * (za / ((double)na * a0) + zb / ((double)nb * a0));
    double radial = exp(-zeta * R_ang);
    /* Angular factor: lobe alignment along the bond axis. */
    Vec3 u = vec3_normalize(dir_a_to_b);
    double th = (vec3_norm(u) < 1e-12) ? 0.0 : acos(u.z);
    double ph = atan2(u.y, u.x);
    /* Use sigma-type (ml=0) projection as the orientation-averaged
     * directional measure; s shells are isotropic by construction. */
    double ya = (la == 0) ? 1.0 : fabs(qm_Y_real(la, 0, th, ph)) / 0.4886025119029199;
    double yb = (lb == 0) ? 1.0 : fabs(qm_Y_real(lb, 0, M_PI - th, ph)) / 0.4886025119029199;
    if (ya > 1.0) ya = 1.0;
    if (yb > 1.0) yb = 1.0;
    return radial * ya * yb;
}

double qm_bond_order(double S, double S_ref) {
    if (S_ref <= 1e-12) return 0.0;
    double bo = 3.0 * S / S_ref;
    if (bo < 0) bo = 0;
    if (bo > 3) bo = 3;
    return bo;
}

double qm_pauli(double S, double J_a, double J_b) {
    double A = 0.5 * (J_a + J_b);
    if (A < 0) A = 0;
    return A * S * S;
}

/* ── QEq solver ────────────────────────────────────────────────────── */

int qm_qeq(const Simulation *sim, double total_q, double dielectric,
           double *out_q) {
    int n = sim->num_atoms;
    if (n < 1 || n > 128 || !out_q) return -1;
    if (dielectric < 1e-9) dielectric = 1.0;
    /* Augmented system (n+1)x(n+1): [A 1; 1^T 0] [q; -mu] = [-chi; total]. */
    static double M[129*129], rhs[129], sol[129];
    int N = n + 1;
    double chi[128], J[128];
    for (int i = 0; i < n; i++) {
        qm_chi_J(sim->atoms[i].element, &chi[i], &J[i]);
    }
    for (int i = 0; i < N * N; i++) M[i] = 0.0;
    for (int i = 0; i < n; i++) {
        M[i * N + i] = J[i];
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            double r = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
            if (r < 0.2) r = 0.2; /* regularize coincidence */
            M[i * N + j] = COULOMB_MD / (r * dielectric);
        }
        M[i * N + n] = 1.0;
        M[n * N + i] = 1.0;
        rhs[i] = -chi[i];
    }
    M[n * N + n] = 0.0;
    rhs[n] = total_q;
    /* Gaussian elimination with partial pivot. */
    for (int c = 0; c < N; c++) {
        int piv = c;
        double best = fabs(M[c * N + c]);
        for (int r = c + 1; r < N; r++) {
            double v = fabs(M[r * N + c]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-9) return -1;
        if (piv != c) {
            for (int k = c; k < N; k++) {
                double t = M[c * N + k]; M[c * N + k] = M[piv * N + k]; M[piv * N + k] = t;
            }
            double t = rhs[c]; rhs[c] = rhs[piv]; rhs[piv] = t;
        }
        double d = M[c * N + c];
        for (int r = c + 1; r < N; r++) {
            double f = M[r * N + c] / d;
            if (f == 0.0) continue;
            for (int k = c; k < N; k++) M[r * N + k] -= f * M[c * N + k];
            rhs[r] -= f * rhs[c];
        }
    }
    for (int r = N - 1; r >= 0; r--) {
        double acc = rhs[r];
        for (int k = r + 1; k < N; k++) acc -= M[r * N + k] * sol[k];
        if (fabs(M[r * N + r]) < 1e-12) return -1;
        sol[r] = acc / M[r * N + r];
    }
    for (int i = 0; i < n; i++) out_q[i] = sol[i];
    return 0;
}

int qm_qeq_pinned(const Simulation *sim, double total_q, double dielectric,
                  int pinned_idx, double pinned_q, double *out_q) {
    int n = sim->num_atoms;
    if (n < 2 || n > 128 || !out_q) return -1;
    if (pinned_idx < 0 || pinned_idx >= n) return -1;
    if (dielectric < 1e-9) dielectric = 1.0;
    /* Free-atom index map. */
    int idx[128], m = 0;
    for (int i = 0; i < n; i++) if (i != pinned_idx) idx[m++] = i;
    double target = total_q - pinned_q;
    /* Pinned atom contributes a fixed background potential:
     * rhs_i = -chi_i - C_i,pinned * pinned_q. */
    int N = m + 1;
    static double M[129*129], rhs[129], sol[129];
    double chi[128], J[128];
    for (int i = 0; i < n; i++) qm_chi_J(sim->atoms[i].element, &chi[i], &J[i]);
    for (int i = 0; i < N * N; i++) M[i] = 0.0;
    for (int a = 0; a < m; a++) {
        int i = idx[a];
        M[a * N + a] = J[i];
        for (int b = 0; b < m; b++) {
            int j = idx[b];
            if (i == j) continue;
            double r = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
            if (r < 0.2) r = 0.2;
            M[a * N + b] = COULOMB_MD / (r * dielectric);
        }
        M[a * N + m] = 1.0;
        M[m * N + a] = 1.0;
        double rp = vec3_dist(sim->atoms[i].position, sim->atoms[pinned_idx].position);
        if (rp < 0.2) rp = 0.2;
        rhs[a] = -chi[i] - COULOMB_MD / (rp * dielectric) * pinned_q;
    }
    M[m * N + m] = 0.0;
    rhs[m] = target;
    for (int c = 0; c < N; c++) {
        int piv = c;
        double best = fabs(M[c * N + c]);
        for (int r = c + 1; r < N; r++) {
            double v = fabs(M[r * N + c]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-9) return -1;
        if (piv != c) {
            for (int k = c; k < N; k++) {
                double t = M[c * N + k]; M[c * N + k] = M[piv * N + k]; M[piv * N + k] = t;
            }
            double t = rhs[c]; rhs[c] = rhs[piv]; rhs[piv] = t;
        }
        double d = M[c * N + c];
        for (int r = c + 1; r < N; r++) {
            double f = M[r * N + c] / d;
            if (f == 0.0) continue;
            for (int k = c; k < N; k++) M[r * N + k] -= f * M[c * N + k];
            rhs[r] -= f * rhs[c];
        }
    }
    for (int r = N - 1; r >= 0; r--) {
        double acc = rhs[r];
        for (int k = r + 1; k < N; k++) acc -= M[r * N + k] * sol[k];
        if (fabs(M[r * N + r]) < 1e-12) return -1;
        sol[r] = acc / M[r * N + r];
    }
    for (int a = 0; a < m; a++) out_q[idx[a]] = sol[a];
    out_q[pinned_idx] = pinned_q;
    return 0;
}

int qm_refresh_charges(Simulation *sim, double total_q, double dielectric) {
    double q[128];
    if (sim->num_atoms > 128) return -1;
    if (qm_qeq(sim, total_q, dielectric, q) != 0) return -1;
    for (int i = 0; i < sim->num_atoms; i++)
        sim->atoms[i].partial_charge = q[i];
    return 0;
}
