#include <math.h>
#include <string.h>
#include "../include/qm.h"
#include "../include/quantum.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Forward: SCF-consistent spatial Zeff (defined with screening block). */
static double qm_spatial_zeff(const Atom *atom, int n, int l,
                              const ElectronConfig *cfg);

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
    if (l == 3) {
        /* Real cubic f harmonics (Condon-Shortley, normalized):
         * ml=0: z3 ~ Y30 = sqrt(7/16pi)(5cos^3-3cos)
         * ml=±1: xz2, yz2 ~ sqrt(21/32pi) sin(2t)... use standard set:
         * f_z3, f_xz2, f_yz2, f_zx2y2 (ml=±2 pair), f_xyz, f_xxxyyy. */
        double st2 = st * st, ct2 = ct * ct;
        const double f0 = 0.3731763325901154;   /* sqrt(7/16pi), z3 */
        const double f1 = 0.4570457990904658;   /* sqrt(21/32pi), xz2/yz2 */
        const double f2 = 1.445305721320277;    /* (1/4)*sqrt(105/pi), verified by direct norm integral */
        const double f3 = 0.5900435899266435;   /* sqrt(35/32pi), xxx/yyy */
        if (ml == 0) return f0 * (5.0 * ct2 * ct - 3.0 * ct);            /* fz3 */
        if (ml == 1) return f1 * (5.0 * ct2 - 1.0) * st * cos(phi);      /* fxz2 */
        if (ml == -1) return f1 * (5.0 * ct2 - 1.0) * st * sin(phi);     /* fyz2 */
        if (ml == 2) return f2 * st2 * ct * cos(2.0 * phi);              /* fzx2-y2 */
        if (ml == -2) return f2 * st2 * ct * sin(2.0 * phi);             /* fzyz-like */
        if (ml == 3) return f3 * st2 * st * cos(3.0 * phi);              /* fxxx-yyy */
        if (ml == -3) return f3 * st2 * st * sin(3.0 * phi);             /* fxyy-xxx */
        return 0.0;
    }
    return 0.0; /* g and above: no chemistry in H..Kr demos */
}

double qm_psi(int n, int l, int ml, double Zeff, Vec3 pos) {
    if (n < 1 || l < 0 || l >= n) return 0.0;
    if (!(Zeff > 0.0) || !isfinite(Zeff)) return 0.0;
    double r = vec3_norm(pos);
    if (!isfinite(r)) return 0.0;
    double R = quantum_radial_wavefunction(n, l, Zeff, r);
    double theta = 0.0, phi = 0.0;
    if (r > 1e-12) {
        double cz = pos.z / r;
        if (cz > 1.0) cz = 1.0;
        else if (cz < -1.0) cz = -1.0;
        theta = acos(cz);
        phi = atan2(pos.y, pos.x);
    }
    return R * qm_Y_real(l, ml, theta, phi);
}

/* ── Hybridization ─────────────────────────────────────────────────── */

QmHybrid qm_hybridization(const Atom *atom) {
    QmHybrid h;
    memset(&h, 0, sizeof h);
    if (!atom || !atom->element) {
        snprintf(h.label, sizeof h.label, "none");
        return h;
    }
    /* Valence-shell s/p census: use the outermost occupied shell, not
     * hardcoded n=2, so Si/P/S/Cl (n=3) are classified, not silently
     * dropped to the steric fallback. */
    int vshell = 0;
    for (int s = 0; s < MAX_SHELLS; s++) {
        int c = 0;
        for (int sub = 0; sub < 4; sub++) c += atom->electron_config.config[s][sub];
        if (c > 0) vshell = s;
    }
    int s_count = atom->electron_config.config[vshell][0];
    int p_count = atom->electron_config.config[vshell][1];
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
        int has_s = s_count > 0;
        int has_p = p_count > 0;
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
    if (!el) return;
    double ie = el->ionization_energy;
    double ea = el->electron_affinity; /* table uses 0 for negative/unknown */
    /* NOTE: clamping EA<0 to 0 biases J=(IE-EA)/2 low for species whose
     * true EA is negative (noble gases, Be/Mg-like s2). Kept for table
     * consistency; do not cite J for those elements beyond order of
     * magnitude. */
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
    if (!atom) return 0.0;
    /* SCF-consistent spatial Zeff (CR where tabulated) for the radius;
     * energies elsewhere stay Slater (documented split). */
    double zeff = qm_spatial_zeff(atom, bn, bl, &cfg);
    /* Charge-responsive radius: live partial charge shifts screening.
     * Neutral atoms (q=0, e.g. fresh diagnostics) unaffected. */
    if (isfinite(atom->partial_charge) && atom->partial_charge != 0.0) {
        double g = qm_gamma_atom(atom, bn, bl);
        zeff += g * atom->partial_charge;
        if (!(zeff > 0.5)) zeff = 0.5;
    }
    double rmp = quantum_most_probable_radius(bn, bl, zeff);
    return rmp * rmp * rmp; /* sphere-volume order of magnitude */
}

/* ── Overlap / bond order / Pauli ──────────────────────────────────── */

/* Charge-response slope dZeff/dq from integer Slater configs:
 * gamma = (Zeff(N+1) - Zeff(N-1))/2 with N = cfg total electrons
 * (q = Z - N so dZeff/dq = -dZeff/dN... note Zeff(N+1) < Zeff(N-1)
 * since added electrons screen: gamma is typically negative*... check:
 * more electrons → more screening → LOWER Zeff, so Zeff(N+1)-Zeff(N-1)
 * < 0, gamma < 0, and Zeff(q) = Zeff(0) + gamma*q gives anions
 * (q<0) HIGHER... wait that is backwards. Redo: anion (extra electron)
 * should be MORE diffuse = LOWER Zeff. q=-1 → Zeff + gamma*(-1) must
 * decrease → gamma > 0. But finite difference gives negative?!
 * Resolution: N+1 config ADDS screening (+0.35ish) so Zeff drops;
 * dZeff/dN intra-shell ≈ -0.35; dq = -dN → dZeff/dq = +0.35 > 0. The
 * formula below computes (Zeff(N-1) - Zeff(N+1))/2 = +0.35. Correct. */
static double qm_screening_gamma(int Z, int n, int l, const ElectronConfig *cfg) {
    if (!cfg) return 0.0;
    int N = cfg->total_electrons;
    ElectronConfig cm, cp;
    int Nm = (N - 1 >= 0) ? N - 1 : 0;
    int Np = (N + 1 <= MAX_ELECTRONS) ? N + 1 : MAX_ELECTRONS;
    if (Np <= Nm) return 0.0;
    pt_electron_config_n(Z, Nm, &cm);
    pt_electron_config_n(Z, Np, &cp);
    double zm = quantum_zeff_raw(Z, n, l, &cm);
    double zp = quantum_zeff_raw(Z, n, l, &cp);
    double g = 0.5 * (zm - zp);
    if (!isfinite(g)) return 0.0;
    return g;
}

/* Public wrapper: screening slope for an atom's (n,l) shell. */
double qm_gamma_atom(const Atom *atom, int n, int l) {
    if (!atom) return 0.0;
    return qm_screening_gamma(atom->Z, n, l, &atom->electron_config);
}

/* Spatial Zeff for ranges/sizes: Clementi-Raimondi SCF exponent where
 * tabulated (neutral ground state; valid for ions on inner shells
 * since higher shells don't screen — e.g. K-3p 7.73 for K+, Na-2p
 * 6.802 for Na+), Slater fallback otherwise. NOT for energies (CR
 * exponents fit shapes, not IEs — documented in quantum.h). */
static double qm_spatial_zeff(const Atom *atom, int n, int l,
                              const ElectronConfig *cfg) {
    if (!atom) return quantum_zeff(1, n, l, cfg);
    double zcr = quantum_zeff_cr(atom->Z, n, l);
    if (zcr > 0.0 && isfinite(zcr)) return zcr;
    if (!cfg) return 1.0;
    return quantum_zeff(atom->Z, n, l, cfg);
}

/* Shell-max |Y_lm| over ml and direction (for lobe-max overlap):
 * searched once on a deterministic grid (no hand constants). */
static double qm_shell_ymax(int l) {
    static double cache[4] = {-1,-1,-1,-1};
    if (l < 0 || l > 3) return 1.0;
    if (cache[l] > 0) return cache[l];
    double best = 0.0;
    for (int it = 0; it < 40; it++) {
        double th = (it + 0.5) * 3.14159265358979323846 / 40.0;
        for (int ip = 0; ip < 80; ip++) {
            double ph = (ip + 0.5) * 2.0 * 3.14159265358979323846 / 80.0;
            for (int m = -l; m <= l; m++) {
                double v = fabs(qm_Y_real(l, m, th, ph));
                if (v > best) best = v;
            }
        }
    }
    if (!(best > 0)) best = 1.0;
    cache[l] = best;
    return best;
}

/* Max-ml projection of shell l onto direction (th,ph): the lobe most
 * aligned with the bond carries sigma; returns share in [0,1]. */
static double qm_lobe_max(int l, double th, double ph) {
    if (l == 0) return 1.0;
    double ym = qm_shell_ymax(l);
    double best = 0.0;
    for (int m = -l; m <= l; m++) {
        double v = fabs(qm_Y_real(l, m, th, ph)) / ym;
        if (v > best) best = v;
    }
    if (best > 1.0) best = 1.0;
    if (!isfinite(best)) return 0.0;
    return best;
}

/* Full overlap engine: charge-responsive Zeff + sigma/pi channels. */
static double qm_overlap_full(const Atom *a, const Atom *b, double R_ang,
                              Vec3 dir_a_to_b, double qa, double qb, int want_pi) {
    if (!a || !b) return 0.0;
    if (R_ang < 1e-9) return want_pi ? 0.0 : 1.0;
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
            /* Tie-break toward higher (n,l): Slater-degenerate s/p picks
             * the chemistry-relevant p shell (same rule as Demo 1/qm_alpha). */
            if (!fa || e > ea || (e == ea && (n > na || (n == na && l > la)))) { fa = 1; ea = e; na = n; la = l; }
        }
        if (cb.config[n-1][l]) {
            double e = quantum_orbital_energy(b->Z, n, l, &cb);
            if (!fb || e > eb || (e == eb && (n > nb || (n == nb && l > lb)))) { fb = 1; eb = e; nb = n; lb = l; }
        }
    }
    double za = qm_spatial_zeff(a, na, la, &ca);
    double zb = qm_spatial_zeff(b, nb, lb, &cb);
    /* Charge response: anions diffuse (lower Zeff), cations contract.
     * gamma from integer-config finite differences (no tables). qa=qb=0
     * reproduces the legacy neutral result exactly. */
    if (isfinite(qa) && qa != 0.0) {
        double g = qm_screening_gamma(a->Z, na, la, &ca);
        za += g * qa;
        if (!(za > 0.5)) za = 0.5;
    }
    if (isfinite(qb) && qb != 0.0) {
        double g = qm_screening_gamma(b->Z, nb, lb, &cb);
        zb += g * qb;
        if (!(zb > 0.5)) zb = 0.5;
    }
    double a0 = BOHR_TO_ANGSTROM;
    double zeta = 0.5 * (za / ((double)na * a0) + zb / ((double)nb * a0));
    if (!isfinite(zeta) || zeta <= 0) return 0.0;
    double radial = exp(-zeta * R_ang);
    /* Angular factor: lobe alignment along the bond axis.
     * Normalize each l by its own sigma-type (ml=0) maximum so s/p/d
     * directionality is comparable: p max = p_norm, d(z2) max = 2*d_norm_a
     * (=0.6308). Dividing d by p_norm saturated at 1.0 and destroyed
     * d-directionality; per-l maxima preserve it. s shells isotropic. */
    Vec3 u = vec3_normalize(dir_a_to_b);
    double un = vec3_norm(u);
    double th, ph;
    if (un < 1e-12) { th = 0.0; ph = 0.0; }
    else {
        double cz = u.z;
        if (cz > 1.0) cz = 1.0;
        else if (cz < -1.0) cz = -1.0;
        th = acos(cz);
        ph = atan2(u.y, u.x);
    }
    double ya, yb;
    /* Lobe-max projection: the ml lobe most aligned with the bond axis
     * carries sigma (a px lobe along x gives full sigma, not zero).
     * The old ml=0-only projection scored correctly-aligned chemistry
     * as zero whenever the lab z differed from the bond (e.g. O-O
     * along x with pz lobes: sigma ~1e-36). Sum-rule pi = remainder. */
    ya = (la <= 3) ? qm_lobe_max(la, th, ph) : 0.0;
    {
        double thb = M_PI - th;
        yb = (lb <= 3) ? qm_lobe_max(lb, thb, ph) : 0.0;
    }
    if (!want_pi) return radial * ya * yb;
    /* Pi channel via the l-shell sum rule: normalized sigma shares
     * ya^2+..., the remainder bounds pi. Per-atom pi share:
     * pi_a = sqrt(max(0, C_l - ya^2)) with C_l = (2l+1)/4pi normalized
     * to the sigma-maximum convention... simpler documented bound:
     * pi share = sqrt(max(0,1-ya^2)) for p (ya<=1 by clamp), times the
     * same for b. s shells (ya=1) give pi=0 correctly. */
    if (la != 1 || lb != 1) return 0.0;
    double pa = sqrt(ya >= 1.0 ? 0.0 : 1.0 - ya * ya);
    double pb = sqrt(yb >= 1.0 ? 0.0 : 1.0 - yb * yb);
    return radial * pa * pb;
}

double qm_overlap(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b) {
    return qm_overlap_full(a, b, R_ang, dir_a_to_b, 0.0, 0.0, 0);
}

double qm_overlap_q(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b,
                    double qa, double qb) {
    return qm_overlap_full(a, b, R_ang, dir_a_to_b, qa, qb, 0);
}

double qm_overlap_pi(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b) {
    return qm_overlap_full(a, b, R_ang, dir_a_to_b, 0.0, 0.0, 1);
}

double qm_overlap_ref(const Atom *a, const Atom *b, Vec3 dir_a_to_b) {
    if (!a || !b || !a->element || !b->element) return 0.0;
    double rref = a->element->covalent_radius + b->element->covalent_radius;
    if (!(rref > 0.0) || !isfinite(rref)) return 0.0;
    return qm_overlap(a, b, rref, dir_a_to_b);
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
    if (!sim) return -1;
    int n = sim->num_atoms;
    if (n < 1 || n > 128 || !out_q) return -1;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    if (!isfinite(total_q)) return -1;
    /* Augmented system (n+1)x(n+1): [A 1; 1^T 0] [q; -mu] = [-chi; total].
     * Thread-local so concurrent calls do not race on shared statics. */
    static _Thread_local double M[129*129], rhs[129], sol[129];
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
    if (!sim) return -1;
    int n = sim->num_atoms;
    if (n < 2 || n > 128 || !out_q) return -1;
    if (pinned_idx < 0 || pinned_idx >= n) return -1;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    if (!isfinite(total_q) || !isfinite(pinned_q)) return -1;
    /* Free-atom index map. */
    int idx[128], m = 0;
    for (int i = 0; i < n; i++) if (i != pinned_idx) idx[m++] = i;
    double target = total_q - pinned_q;
    /* Pinned atom contributes a fixed background potential:
     * rhs_i = -chi_i - C_i,pinned * pinned_q. */
    int N = m + 1;
    static _Thread_local double M[129*129], rhs[129], sol[129];
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
    if (!sim) return -1;
    if (sim->num_atoms > 128) return -1;
    if (!isfinite(total_q)) return -1;
    if (qm_qeq(sim, total_q, dielectric, q) != 0) return -1;
    for (int i = 0; i < sim->num_atoms; i++) {
        if (!isfinite(q[i])) return -1;
        sim->atoms[i].partial_charge = q[i];
    }
    return 0;
}

/* ── v2: polarizability ──────────────────────────────────────────── */

double qm_polarizability(const Atom *atom) {
    if (!atom || !atom->element) return 0.0;
    /* Closed-shell ions use ion data (consistent with qm_alpha_IE):
     * neutral-atom fallback would give K+ the 4s r_mp^3 volume and Na+
     * a 3s volume — qualitatively wrong (the valence shell is GONE).
     * K+: 0.83 A^3, Na+: 0.18 A^3 (Pauling crystal values). This is THE
     * term that lets off-center K+ out-polarize Na+ (4.6x): at symmetric
     * sites the ion field cancels and it sleeps, which is why every
     * symmetric test showed identical induction. Asymmetry wakes it. */
    if (atom->formal_charge != 0) {
        if (atom->Z == 19) return 0.83;
        if (atom->Z == 11) return 0.18;
    }
    /* Applequist isotropic volumes (A^3), standard set: H 0.42, C 1.35,
     * N 1.10, O 0.84. Carbonyl-O uses 0.84 (matches KcsA legacy leg).
     * Other elements: r_mp^3 sphere fallback (order-of-magnitude). */
    switch (atom->Z) {
        case 1: return 0.42;
        case 6: return 1.35;
        case 7: return 1.10;
        case 8: return 0.84;
        default: break;
    }
    return qm_alpha(atom);
}

/* Thole damping width (A) for induction: f(r) = 1 - exp(-(r/a)^3).
 * At r >> a, f→1 (exact point field); at r→0, f~(r/a)^3 kills the
 * 1/r^2 singularity (E_d ~ r, U finite). a=2.0 keeps first-shell
 * ion-O (2.35-2.75 A, f=0.80-0.93) near-exact while preventing the
 * polarization catastrophe that collapsed the undamped hydration
 * minimization (E_pol=-28000 eV). Standard Thole form, documented. */
#define QM_THOLE_A 2.0
static double qm_thole_f(double r) {
    double u = r / QM_THOLE_A;
    return 1.0 - exp(-u * u * u);
}
static double qm_thole_df(double r) {
    double u = r / QM_THOLE_A;
    return 3.0 * u * u * exp(-u * u * u) / QM_THOLE_A;
}

/* Field at each atom from all other point charges (V/A), Thole-damped. */
static void qm_fields(const Simulation *sim, double dielectric, Vec3 *Eout) {
    int N = sim->num_atoms;
    for (int i = 0; i < N; i++) Eout[i] = vec3_zero();
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    for (int i = 0; i < N; i++) {
        Vec3 E = vec3_zero();
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double qj = sim->atoms[j].partial_charge;
            if (!isfinite(qj) || fabs(qj) < 1e-12) continue;
            Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
            double r2 = vec3_norm2(d);
            if (r2 < 1e-8 || !isfinite(r2)) continue;
            if (r2 > sim->cutoff * sim->cutoff) continue;
            double r = sqrt(r2);
            double f = COULOMB_MD * qj / (r2 * r * dielectric) * qm_thole_f(r);
            vec3_iadd(&E, vec3_scale(d, f));
        }
        Eout[i] = E;
    }
}

double qm_induction_energy(const Simulation *sim, double dielectric) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    int N = sim->num_atoms;
    if (N > 256) return 0.0;
    Vec3 E[256];
    qm_fields(sim, dielectric, E);
    const double C = 0.069446; /* eV per (V^2·A), see KcsA derivation */
    double U = 0.0;
    for (int i = 0; i < N; i++) {
        double a = qm_polarizability(&sim->atoms[i]);
        if (!(a > 0.0) || !isfinite(a)) continue;
        double e2 = vec3_norm2(E[i]);
        if (!isfinite(e2)) continue;
        U += -0.5 * C * a * e2;
    }
    return U;
}

/* Shared induction force kernel: accumulates F_k = Σ_i V_i·dE0_i/dR_k
 * scaled by s_i, with Thole-damped dE0/dR (f and df terms).
 * First-order: V=E0, s=C*alpha. Coupled (Hellmann-Feynman): V=mu, s=1. */
static void qm_induction_apply(const Simulation *sim, double dielectric,
                               const Vec3 *V, const double *s) {
    int N = sim->num_atoms;
    for (int k = 0; k < N; k++) {
        Vec3 Fk = vec3_zero();
        for (int i = 0; i < N; i++) {
            if (s[i] == 0.0) continue;
            Vec3 Vi = V[i];
            if (i == k) {
                for (int j = 0; j < N; j++) {
                    if (j == i) continue;
                    double qj = sim->atoms[j].partial_charge;
                    if (!isfinite(qj) || fabs(qj) < 1e-12) continue;
                    Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
                    double r2 = vec3_norm2(d);
                    if (r2 < 1e-8 || !isfinite(r2)) continue;
                    double r = sqrt(r2);
                    if (r > sim->cutoff) continue;
                    double pre = COULOMB_MD * qj / dielectric;
                    double fth = qm_thole_f(r), dfth = qm_thole_df(r);
                    double r3 = r2 * r, r5 = r3 * r2;
                    double dVi = vec3_dot(d, Vi);
                    Vec3 TE = vec3_sub(vec3_scale(Vi, fth / r3),
                                       vec3_scale(d, 3.0 * fth * dVi / r5));
                    vec3_iadd(&TE, vec3_scale(d, dfth * dVi / (r2 * r2)));
                    vec3_iadd(&Fk, vec3_scale(TE, s[i] * pre));
                }
            } else {
                double qk = sim->atoms[k].partial_charge;
                if (!isfinite(qk) || fabs(qk) < 1e-12) continue;
                Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[k].position);
                double r2 = vec3_norm2(d);
                if (r2 < 1e-8 || !isfinite(r2)) continue;
                double r = sqrt(r2);
                if (r > sim->cutoff) continue;
                double pre = COULOMB_MD * qk / dielectric;
                double fth = qm_thole_f(r), dfth = qm_thole_df(r);
                double r3 = r2 * r, r5 = r3 * r2;
                double dVi = vec3_dot(d, Vi);
                Vec3 TE = vec3_sub(vec3_scale(Vi, fth / r3),
                                   vec3_scale(d, 3.0 * fth * dVi / r5));
                vec3_iadd(&TE, vec3_scale(d, dfth * dVi / (r2 * r2)));
                vec3_isub(&Fk, vec3_scale(TE, s[i] * pre));
            }
        }
        if (isfinite(Fk.x) && isfinite(Fk.y) && isfinite(Fk.z))
            vec3_iadd(&sim->atoms[k].force, Fk);
    }
}

double qm_induction_forces(Simulation *sim, double dielectric) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    int N = sim->num_atoms;
    if (N > 256) return 0.0;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    Vec3 E[256];
    qm_fields(sim, dielectric, E);
    const double C = 0.069446;
    double alpha[256];
    for (int i = 0; i < N; i++) {
        alpha[i] = qm_polarizability(&sim->atoms[i]);
        if (!(alpha[i] > 0.0) || !isfinite(alpha[i])) alpha[i] = 0.0;
    }
    double U = 0.0;
    for (int i = 0; i < N; i++) U += -0.5 * C * alpha[i] * vec3_norm2(E[i]);
    if (!isfinite(U)) U = 0.0;
    Vec3 V[256];
    double s[256];
    for (int i = 0; i < N; i++) { V[i] = E[i]; s[i] = C * alpha[i]; }
    qm_induction_apply(sim, dielectric, V, s);
    return U;
}

/* ── v4: self-consistent dipoles (dipole-dipole coupling) ─────────────── */

/* Solve mu = alpha*E(mu) with Thole-damped dipole-dipole coupling:
 * E_i = E0_i + Σ_j f(r_ij)*k*[3(d·mu_j)d/r^5 - mu_j/r^3]/D.
 * mu in e·A (mu = alpha*E/C_MD). Plain iteration + 0.5 mixing, tol
 * 1e-6 relative, max 100. Returns iters, or -1 on non-convergence
 * (mu_out holds last iterate; caller degrades gracefully). */
int qm_solve_dipoles(const Simulation *sim, double dielectric, Vec3 *mu_out) {
    if (!sim || !sim->atoms || !mu_out) return -1;
    int N = sim->num_atoms;
    if (N < 1 || N > 256) return -1;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    Vec3 E0[256];
    qm_fields(sim, dielectric, E0);
    double alpha[256];
    for (int i = 0; i < N; i++) {
        alpha[i] = qm_polarizability(&sim->atoms[i]);
        if (!(alpha[i] > 0.0) || !isfinite(alpha[i])) alpha[i] = 0.0;
        mu_out[i] = vec3_scale(E0[i], alpha[i] / COULOMB_MD);
    }
    Vec3 new_mu[256];
    for (int it = 0; it < 100; it++) {
        double maxd = 0.0, maxm = 0.0;
        for (int i = 0; i < N; i++) {
            Vec3 E = E0[i];
            for (int j = 0; j < N; j++) {
                if (j == i || alpha[j] == 0.0) continue;
                Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
                double r2 = vec3_norm2(d);
                if (r2 < 1e-8 || !isfinite(r2)) continue;
                if (r2 > sim->cutoff * sim->cutoff) continue;
                double r = sqrt(r2);
                double fth = qm_thole_f(r);
                double r3 = r2 * r, r5 = r3 * r2;
                double dmu = vec3_dot(d, mu_out[j]);
                Vec3 Tmu = vec3_sub(vec3_scale(d, 3.0 * dmu / r5),
                                    vec3_scale(mu_out[j], 1.0 / r3));
                vec3_iadd(&E, vec3_scale(Tmu, COULOMB_MD * fth / dielectric));
            }
            Vec3 target = vec3_scale(E, alpha[i] / COULOMB_MD);
            /* 0.5 mixing for stability. */
            new_mu[i] = vec3_add(vec3_scale(mu_out[i], 0.5), vec3_scale(target, 0.5));
            double dm = vec3_norm(vec3_sub(new_mu[i], mu_out[i]));
            double mm = vec3_norm(new_mu[i]);
            if (dm > maxd) maxd = dm;
            if (mm > maxm) maxm = mm;
        }
        for (int i = 0; i < N; i++) mu_out[i] = new_mu[i];
        if (!isfinite(maxd) || !isfinite(maxm)) return -1;
        if (maxd < 1e-6 * (maxm > 1e-9 ? maxm : 1e-9)) return it + 1;
    }
    return -1;
}

/* SCF induction energy + Hellmann-Feynman analytic forces.
 * U = -1/2 Σ mu·E0 (variational in mu: dmu/dR terms cancel at
 * convergence, so forces need only dE0/dR — the same TE kernel with
 * contraction vectors mu and unit scalars). Falls back to first-order
 * weights if the dipole solver fails to converge (flagged in print). */
double qm_induction_scf_forces(Simulation *sim, double dielectric, int *iters_out) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    int N = sim->num_atoms;
    if (N > 256) return 0.0;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    Vec3 E0[256], mu[256];
    qm_fields(sim, dielectric, E0);
    int it = qm_solve_dipoles(sim, dielectric, mu);
    if (iters_out) *iters_out = it;
    double s[256];
    Vec3 V[256];
    if (it < 0) {
        /* Graceful degradation: first-order weights. */
        const double C = 0.069446;
        for (int i = 0; i < N; i++) {
            double a = qm_polarizability(&sim->atoms[i]);
            if (!(a > 0.0) || !isfinite(a)) a = 0.0;
            V[i] = E0[i]; s[i] = C * a;
        }
    } else {
        for (int i = 0; i < N; i++) { V[i] = mu[i]; s[i] = 1.0; }
    }
    double U = 0.0;
    if (it < 0) {
        const double C = 0.069446;
        for (int i = 0; i < N; i++) {
            double a = qm_polarizability(&sim->atoms[i]);
            if (!(a > 0.0)) continue;
            U += -0.5 * C * a * vec3_norm2(E0[i]);
        }
    } else {
        for (int i = 0; i < N; i++) U += -0.5 * vec3_dot(mu[i], E0[i]);
    }
    if (!isfinite(U)) U = 0.0;
    qm_induction_apply(sim, dielectric, V, s);
    return U;
}

/* ── v2: Pauli energy + FD forces ────────────────────────────────── */
static int qm_pair_excluded(const Simulation *sim, int i, int j) {
    const Atom *ai = &sim->atoms[i];
    for (int p = 0; p < ai->num_bonds; p++)
        if (ai->bond_partners[p] == j) return 1;
    for (int p = 0; p < ai->num_bonds; p++) {
        int k = ai->bond_partners[p];
        if (k < 0 || k >= sim->num_atoms) continue;
        const Atom *ak = &sim->atoms[k];
        for (int q = 0; q < ak->num_bonds; q++)
            if (ak->bond_partners[q] == j) return 1;
    }
    for (int a = 0; a < sim->num_angles; a++) {
        const Angle *ang = &sim->angles[a];
        if ((ang->atom_a == i && ang->atom_c == j) ||
            (ang->atom_a == j && ang->atom_c == i)) return 1;
    }
    return 0;
}

static double qm_pair_pauli(const Simulation *sim, int i, int j) {
    Vec3 d = vec3_sub(sim->atoms[j].position, sim->atoms[i].position);
    double r2 = vec3_norm2(d);
    if (!(r2 > 1e-18) || !isfinite(r2)) return 0.0;
    if (r2 > sim->cutoff * sim->cutoff) return 0.0;
    double r = sqrt(r2);
    Vec3 u = vec3_normalize(d);
    if (vec3_norm(u) < 1e-12) u = vec3(0, 0, 1);
    /* Live partial charges feed the charge-responsive screening, so the
     * Pauli wall breathes with SCF charge flow (neutral fallback when
     * charges are exactly zero, e.g. fresh diagnostics). */
    double S = qm_overlap_q(&sim->atoms[i], &sim->atoms[j], r, u,
                            sim->atoms[i].partial_charge,
                            sim->atoms[j].partial_charge);
    if (!isfinite(S) || S <= 0.0) return 0.0;
    double chi_a = 0.0, J_a = 0.0, chi_b = 0.0, J_b = 0.0;
    qm_chi_J(sim->atoms[i].element, &chi_a, &J_a);
    qm_chi_J(sim->atoms[j].element, &chi_b, &J_b);
    if (!isfinite(J_a) || !isfinite(J_b)) return 0.0;
    return qm_pauli(S, J_a, J_b);
}

double qm_pauli_energy(const Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 2) return 0.0;
    double E = 0.0;
    for (int i = 0; i < sim->num_atoms - 1; i++)
        for (int j = i + 1; j < sim->num_atoms; j++) {
            if (qm_pair_excluded(sim, i, j)) continue;
            E += qm_pair_pauli(sim, i, j);
        }
    return isfinite(E) ? E : 0.0;
}

double qm_pauli_forces(Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 2) return 0.0;
    const double h = 1.0e-5;
    int N = sim->num_atoms;
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            if (qm_pair_excluded(sim, i, j)) continue;
            /* Zero-contribution skip: pairs with |E| < 1e-12 eV (overlap
             * below ~1e-6, e.g. all filter-range contacts) carry forces
             * below 2*zeta*E ~ 1e-11 eV/A (E = A*S^2 >= 0, dE/dx bounded
             * by 2*zeta*E) — unrepresentable next to 1e-9 test floors.
             * Skipping their 12 FD evals is exact to 1e-12 eV (typical
             * cage systems skip nearly all pairs). */
            if (fabs(qm_pair_pauli(sim, i, j)) < 1e-12) continue;
            int idx[2] = {i, j};
            for (int a = 0; a < 2; a++) {
                Atom *at = &sim->atoms[idx[a]];
                double *cc[3] = {&at->position.x, &at->position.y, &at->position.z};
                for (int c = 0; c < 3; c++) {
                    double o = *cc[c];
                    *cc[c] = o + h;
                    double Ep = qm_pair_pauli(sim, i, j);
                    *cc[c] = o - h;
                    double Em = qm_pair_pauli(sim, i, j);
                    *cc[c] = o;
                    if (!isfinite(Ep) || !isfinite(Em)) continue;
                    double F = -(Ep - Em) / (2.0 * h);
                    if (!isfinite(F)) continue;
                    /* Full -dE/dx belongs to the perturbed atom. */
                    if (c == 0) at->force.x += F;
                    else if (c == 1) at->force.y += F;
                    else at->force.z += F;
                }
            }
        }
    }
    /* Positions restored exactly; return exact total. */
    return qm_pauli_energy(sim);
}

/* ── v3: ion-aware alpha/IE, Slater-Kirkwood C6, TT damping, dispersion ─ */

void qm_alpha_IE(const Atom *atom, double *alpha_out, double *IE_out) {
    double a = 0.0, ie = 0.0;
    if (!atom || !atom->element) { if (alpha_out) *alpha_out = 0; if (IE_out) *IE_out = 0; return; }
    /* Closed-shell ions: neutral-atom table values are qualitatively
     * wrong (first IE ejects the already-missing electron). Use ion
     * data: K+ α=0.83 A^3 (Pauling), IE2=31.63 eV; Na+ α=0.18,
     * IE2=47.29 eV (NIST). Detected as formal_charge≠0 on K/Na. */
    if (atom->formal_charge != 0) {
        if (atom->Z == 19) { a = 0.83; ie = 31.63; }
        else if (atom->Z == 11) { a = 0.18; ie = 47.29; }
    }
    if (ie <= 0.0) {
        a = qm_polarizability(atom);
        ie = atom->element->ionization_energy;
    }
    if (alpha_out) *alpha_out = a;
    if (IE_out) *IE_out = ie;
}

double qm_c6(const Atom *a, const Atom *b) {
    double aa, ia, ab, ib;
    qm_alpha_IE(a, &aa, &ia);
    qm_alpha_IE(b, &ab, &ib);
    if (!(aa > 0) || !(ab > 0) || !(ia > 0) || !(ib > 0)) return 0.0;
    if (!isfinite(aa) || !isfinite(ab) || !isfinite(ia) || !isfinite(ib)) return 0.0;
    return 1.5 * aa * ab * ia * ib / (ia + ib);
}

/* Mean Slater zeta of the pair's valence shells (1/A). */
static double qm_pair_zeta(const Atom *a, const Atom *b) {
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
            /* Tie-break toward higher (n,l): Slater-degenerate s/p picks
             * the chemistry-relevant p shell (same rule as Demo 1/qm_alpha). */
            if (!fa || e > ea || (e == ea && (n > na || (n == na && l > la)))) { fa = 1; ea = e; na = n; la = l; }
        }
        if (cb.config[n-1][l]) {
            double e = quantum_orbital_energy(b->Z, n, l, &cb);
            if (!fb || e > eb || (e == eb && (n > nb || (n == nb && l > lb)))) { fb = 1; eb = e; nb = n; lb = l; }
        }
    }
    double za = qm_spatial_zeff(a, na, la, &ca);
    double zb = qm_spatial_zeff(b, nb, lb, &cb);
    /* Charge-responsive: live partial charges shift screening (same
     * gamma machinery as overlap). Neutral atoms unaffected. */
    if (isfinite(a->partial_charge) && a->partial_charge != 0.0) {
        double g = qm_screening_gamma(a->Z, na, la, &ca);
        za += g * a->partial_charge;
        if (!(za > 0.5)) za = 0.5;
    }
    if (isfinite(b->partial_charge) && b->partial_charge != 0.0) {
        double g = qm_screening_gamma(b->Z, nb, lb, &cb);
        zb += g * b->partial_charge;
        if (!(zb > 0.5)) zb = 0.5;
    }
    double a0 = BOHR_TO_ANGSTROM;
    double zeta = 0.5 * (za / ((double)na * a0) + zb / ((double)nb * a0));
    if (!isfinite(zeta) || zeta <= 0) return 2.0;
    return zeta;
}

double qm_tt_f6(double b, double r, double *dfdr_out) {
    if (!(b > 0) || !(r > 0) || !isfinite(b) || !isfinite(r)) {
        if (dfdr_out) *dfdr_out = 0.0;
        return 0.0;
    }
    double x = b * r;
    /* Σ_{k=0..6} x^k/k! via recurrence */
    double sum = 1.0, term = 1.0;
    for (int k = 1; k <= 6; k++) { term *= x / k; sum += term; }
    double e = exp(-x);
    double f = 1.0 - e * sum;
    if (dfdr_out) {
        /* df/dr = e*b*x^6/720 */
        double x6 = term * 720.0; /* term ended as x^6/6! */
        *dfdr_out = e * b * x6 / 720.0;
    }
    return f;
}

static double qm_pair_disp(const Simulation *sim, int i, int j) {
    Vec3 d = vec3_sub(sim->atoms[j].position, sim->atoms[i].position);
    double r2 = vec3_norm2(d);
    if (!(r2 > 1e-18) || !isfinite(r2) || r2 > sim->cutoff * sim->cutoff) return 0.0;
    double r = sqrt(r2);
    double c6 = qm_c6(&sim->atoms[i], &sim->atoms[j]);
    if (!(c6 > 0) || !isfinite(c6)) return 0.0;
    double b = qm_pair_zeta(&sim->atoms[i], &sim->atoms[j]);
    double f = qm_tt_f6(b, r, NULL);
    return -c6 * f / (r * r * r * r * r * r);
}

double qm_dispersion_energy(const Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 2) return 0.0;
    double E = 0.0;
    for (int i = 0; i < sim->num_atoms - 1; i++)
        for (int j = i + 1; j < sim->num_atoms; j++) {
            if (qm_pair_excluded(sim, i, j)) continue;
            E += qm_pair_disp(sim, i, j);
        }
    return isfinite(E) ? E : 0.0;
}

double qm_dispersion_forces(Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 2) return 0.0;
    int N = sim->num_atoms;
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            if (qm_pair_excluded(sim, i, j)) continue;
            Vec3 d = vec3_sub(sim->atoms[j].position, sim->atoms[i].position);
            double r2 = vec3_norm2(d);
            if (!(r2 > 1e-18) || !isfinite(r2)) continue;
            double r = sqrt(r2);
            if (r > sim->cutoff) continue;
            double c6 = qm_c6(&sim->atoms[i], &sim->atoms[j]);
            if (!(c6 > 0) || !isfinite(c6)) continue;
            double b = qm_pair_zeta(&sim->atoms[i], &sim->atoms[j]);
            double df;
            double f = qm_tt_f6(b, r, &df);
            /* dE/dr = -c6*(df/r^6 - 6f/r^7); F_i = -dE/dr*d/r */
            double r6 = r2 * r2 * r2, r7 = r6 * r;
            double dEdr = -c6 * (df / r6 - 6.0 * f / r7);
            if (!isfinite(dEdr)) continue;
            Vec3 Fi = vec3_scale(d, -dEdr / r);
            if (!isfinite(Fi.x + Fi.y + Fi.z)) continue;
            vec3_iadd(&sim->atoms[i].force, Fi);
            vec3_isub(&sim->atoms[j].force, Fi);
        }
    }
    return qm_dispersion_energy(sim);
}

/* SCF: pinned-QEq with dipole reaction-field chi shift. */
static int qm_scf_run(Simulation *sim, double total_q, double dielectric,
                      const int *pidx, const double *pinq, int npin) {
    if (!sim || !sim->atoms) return -1;
    int n = sim->num_atoms;
    if (n < 2 || n > 128) return -1;
    if (npin < 1 || npin > 4) return -1;
    for (int pp = 0; pp < npin; pp++) {
        if (pidx[pp] < 0 || pidx[pp] >= n) return -1;
        if (!isfinite(pinq[pp])) return -1;
    }
    if (!isfinite(total_q)) return -1;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    /* Snapshot base chi (neutral table) for shift reference. */
    double chi0[128], J0[128];
    for (int i = 0; i < n; i++) {
        chi0[i] = 0.0; J0[i] = 1.0;
        qm_chi_J(sim->atoms[i].element, &chi0[i], &J0[i]);
    }
    const double C = 0.069446;
    double qprev[128];
    for (int i = 0; i < n; i++) qprev[i] = sim->atoms[i].partial_charge;
    int it;
    for (it = 0; it < 5; it++) {
        /* Dipole field from current charges (Thole-damped). */
        Vec3 E[128];
        for (int i = 0; i < n; i++) E[i] = vec3_zero();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                double qj = qprev[j];
                for (int pp = 0; pp < npin; pp++)
                    if (j == pidx[pp]) { qj = pinq[pp]; break; }
                if (!isfinite(qj) || fabs(qj) < 1e-12) continue;
                Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
                double r2 = vec3_norm2(d);
                if (r2 < 1e-8 || !isfinite(r2)) continue;
                double r = sqrt(r2);
                if (r > sim->cutoff) continue;
                double uu = r / QM_THOLE_A;
                double fth = 1.0 - exp(-uu * uu * uu);
                vec3_iadd(&E[i], vec3_scale(d, COULOMB_MD * qj * fth / (r2 * r * dielectric)));
            }
        }
        /* Reaction potential V_i = -C*Σ_j alpha_j E_j·T_ji damped.
         * Simplified to dipole-charge form: shift chi_i by the dipole
         * potential at i from converged dipoles mu_j = alpha_j*E_j:
         * V_i = Σ_j mu_j·d_ji/(r^3)*fth (V), chi shift = e*V (eV). */
        double Vpot[128] = {0};
        for (int i = 0; i < n; i++) {
            double V = 0.0;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                double aj = qm_polarizability(&sim->atoms[j]);
                if (!(aj > 0)) continue;
                Vec3 mu = vec3_scale(E[j], aj); /* e·A (alpha*E with E in V/A, e implicit) */
                Vec3 dji = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
                double r2 = vec3_norm2(dji);
                if (r2 < 1e-8 || !isfinite(r2)) continue;
                double r = sqrt(r2);
                if (r > sim->cutoff) continue;
                double uu = r / QM_THOLE_A;
                double fth = 1.0 - exp(-uu * uu * uu);
                V += vec3_dot(mu, dji) / (r2 * r) * fth;
            }
            Vpot[i] = V; /* V (volts) since mu in e·A over A^2 */
        }
        /* Dipole reaction field feeds back into electronegativity:
         * chi_eff = chi0 + g*Vpot with feedback gain g=0.5. Full gain
         * (1.0) is a positive loop (closer O → bigger field → bigger
         * dipole → more negative O → stronger pull → closer) that runs
         * to the charge clamps inside dynamics (observed -51 eV Na+
         * collapse); g=0.5 keeps sign/order, converges stably. */
        /* Build (m+1) system like qm_qeq_pinned with shifted chi. */
        int idx[128], m = 0;
        for (int i = 0; i < n; i++) {
            int ispin = 0;
            for (int pp = 0; pp < npin; pp++) if (i == pidx[pp]) { ispin = 1; break; }
            if (!ispin) idx[m++] = i;
        }
        double target = total_q;
        for (int pp = 0; pp < npin; pp++) target -= pinq[pp];
        int Nsys = m + 1;
        static _Thread_local double M[129*129], rhs[129], sol[129];
        double Jv[128];
        for (int i = 0; i < n; i++) { double c; qm_chi_J(sim->atoms[i].element, &c, &Jv[i]); }
        for (int i = 0; i < Nsys * Nsys; i++) M[i] = 0.0;
        for (int a = 0; a < m; a++) {
            int i = idx[a];
            M[a * Nsys + a] = Jv[i];
            for (int bb = 0; bb < m; bb++) {
                int j = idx[bb];
                if (i == j) continue;
                double r = vec3_dist(sim->atoms[i].position, sim->atoms[j].position);
                if (r < 0.2) r = 0.2;
                M[a * Nsys + bb] = COULOMB_MD / (r * dielectric);
            }
            M[a * Nsys + m] = 1.0;
            M[m * Nsys + a] = 1.0;
            double bg = 0.0;
            for (int pp = 0; pp < npin; pp++) {
                double rp = vec3_dist(sim->atoms[i].position, sim->atoms[pidx[pp]].position);
                if (rp < 0.2) rp = 0.2;
                bg += COULOMB_MD / (rp * dielectric) * pinq[pp];
            }
            rhs[a] = -(chi0[i] + 0.5 * Vpot[i]) - bg;
        }
        M[m * Nsys + m] = 0.0;
        rhs[m] = target;
        int ok = 1;
        for (int c = 0; c < Nsys; c++) {
            int piv = c;
            double best = fabs(M[c * Nsys + c]);
            for (int rr = c + 1; rr < Nsys; rr++) {
                double v = fabs(M[rr * Nsys + c]);
                if (v > best) { best = v; piv = rr; }
            }
            if (best < 1e-9) { ok = 0; break; }
            if (piv != c) {
                for (int k = c; k < Nsys; k++) {
                    double t2 = M[c * Nsys + k]; M[c * Nsys + k] = M[piv * Nsys + k]; M[piv * Nsys + k] = t2;
                }
                double t2 = rhs[c]; rhs[c] = rhs[piv]; rhs[piv] = t2;
            }
            double dd = M[c * Nsys + c];
            for (int rr = c + 1; rr < Nsys; rr++) {
                double f2 = M[rr * Nsys + c] / dd;
                if (f2 == 0.0) continue;
                for (int k = c; k < Nsys; k++) M[rr * Nsys + k] -= f2 * M[c * Nsys + k];
                rhs[rr] -= f2 * rhs[c];
            }
        }
        if (!ok) return -1;
        for (int rr = Nsys - 1; rr >= 0; rr--) {
            double acc = rhs[rr];
            for (int k = rr + 1; k < Nsys; k++) acc -= M[rr * Nsys + k] * sol[k];
            if (fabs(M[rr * Nsys + rr]) < 1e-12) return -1;
            sol[rr] = acc / M[rr * Nsys + rr];
        }
        double maxd = 0.0;
        for (int a = 0; a < m; a++) {
            double raw = sol[a];
            /* Under-relaxation + clamp: dipole feedback can overshoot
             * (charge sloshing → collapse in dynamics). Mix 70/30 and
             * cap |q|≤2 e (beyond any physical carbonyl charge). */
            if (raw > 2.0) raw = 2.0;
            else if (raw < -2.0) raw = -2.0;
            double mixed = 0.3 * raw + 0.7 * qprev[idx[a]];
            double dch = fabs(mixed - qprev[idx[a]]);
            if (dch > maxd) maxd = dch;
            qprev[idx[a]] = mixed;
        }
        for (int pp = 0; pp < npin; pp++) qprev[pidx[pp]] = pinq[pp];
        if (maxd < 1e-4) {
            for (int i = 0; i < n; i++) {
                if (!isfinite(qprev[i])) return -1;
                sim->atoms[i].partial_charge = qprev[i];
            }
            (void)C;
            return it + 1;
        }
    }
    for (int i = 0; i < n; i++) {
        if (!isfinite(qprev[i])) return -1;
        sim->atoms[i].partial_charge = qprev[i];
    }
    return 5;
}

int qm_scf_charges(Simulation *sim, double total_q, double dielectric,
                   int pinned_idx, double pinned_q) {
    return qm_scf_run(sim, total_q, dielectric, &pinned_idx, &pinned_q, 1);
}

int qm_scf_charges_2pin(Simulation *sim, double total_q, double dielectric,
                        int pinned0, double q0, int pinned1, double q1) {
    if (pinned0 == pinned1) return -1;
    int pidx[2] = {pinned0, pinned1};
    double pinq[2] = {q0, q1};
    return qm_scf_run(sim, total_q, dielectric, pidx, pinq, 2);
}
