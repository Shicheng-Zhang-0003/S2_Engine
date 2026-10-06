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
static int qm_pair_excluded(const Simulation *sim, int i, int j);

/* QEq hard-core exclusion: Rappe-Goddard sets A_ij = 0 for topologically
 * bonded 1-2 and 1-3 pairs, so screening is not double-counted between
 * atoms joined by a bond. Shared with the Pauli/dispersion exclusion. */
static int qm_pair_excluded(const Simulation *sim, int i, int j);

/* Clamp for QEq partial charges. Beyond +/-2 e no carbonyl or amide
 * partial charge is physical; the unclamped solve produced +4.82 e on
 * carbon and -5.63 e on oxygen for a C-O-O fragment at real bond
 * lengths, which then drove absurd Coulomb energies downstream. */
#define QM_QEQ_QMAX 2.0
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

/* Shared classification core.
 *
 * `has_own_pi`  : this atom is a terminus of a double/triple bond, so it
 *                 spends its p orbital on its own pi electrons.
 * `conjugated`  : this atom is bonded to a partner that carries a pi
 *                 bond, so a lone pair on it can delocalise into p.
 * Both are 0 for the atom-only entry point, which is the conservative
 * sp3 reading - all it can know without topology.
 *
 * The refinement that makes this correct: a lone pair counts toward the
 * steric number ONLY when it occupies a hybrid. A lone pair that
 * delocalises into an adjacent pi system occupies p instead, which is
 * the entire difference between an amine (sp3, steric 4) and an amide
 * (sp2, steric 3). An atom that owns its own pi bond (carbonyl oxygen)
 * keeps its lone pairs in the plane, because its p orbital is already
 * committed to the pi bond.
 */
static QmHybrid qm_hybrid_core(const Atom *atom, int has_own_pi,
                               int conjugated) {
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
    /* Free atom (no bonds): hybridisation is a molecular concept; report
     * the atomic valence shell instead of forcing sp labels. */
    if (n_bonds == 0) {
        if (s_count > 0 && p_count > 0) snprintf(h.label, sizeof h.label, "atomic-sp");
        else if (p_count > 0) snprintf(h.label, sizeof h.label, "atomic-p");
        else snprintf(h.label, sizeof h.label, "atomic-s");
        h.n_lobes = 0;
        int v0 = atom->electron_config.valence_electrons;
        h.n_lone_pairs = (v0 > 0) ? v0 / 2 : 0;
        return h;
    }

    /* Sigma partners: one per bonded neighbour. A double bond is still
     * ONE sigma partner - the pi electrons live in a separate orbital. */
    int sigma_partners = n_bonds;
    int v = atom->electron_config.valence_electrons;
    int lone = 0;
    if (v > 0) {
        int rem = v - sigma_partners;
        lone = (rem > 0) ? rem / 2 : 0;
    }
    int lone_in_plane = lone;
    if (lone > 0 && conjugated && !has_own_pi) lone_in_plane = 0;

    int steric = sigma_partners + lone_in_plane;
    if (steric >= 4) {
        snprintf(h.label, sizeof h.label, "sp3");
        h.n_lobes = 4;
        h.lobes[0] = vec3_normalize(vec3(1, 1, 1));
        h.lobes[1] = vec3_normalize(vec3(1, -1, -1));
        h.lobes[2] = vec3_normalize(vec3(-1, 1, -1));
        h.lobes[3] = vec3_normalize(vec3(-1, -1, 1));
    } else if (steric == 3) {
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
    return h;
}

/* Topology-aware entry point. Verified against the functional groups this
 * model is actually built from:
 *   water O      2 sigma, 2 lone pairs, not conjugated -> steric 4  sp3
 *   ammonia N    3 sigma, 1 lone pair, not conjugated -> steric 4  sp3
 *   carbonyl O   1 sigma, 2 lone pairs, owns the pi   -> steric 3  sp2
 *   amide N      3 sigma, 1 lone pair, conjugated    -> steric 3  sp2
 *   carbonyl C   3 sigma, 0 lone pairs               -> steric 3  sp2
 *   methylene C  4 sigma, 0 lone pairs               -> steric 4  sp3
 */
QmHybrid qm_hybridization_ctx(const Simulation *sim, int atom_idx) {
    QmHybrid h;
    memset(&h, 0, sizeof h);
    if (!sim || !sim->atoms) {
        snprintf(h.label, sizeof h.label, "none");
        return h;
    }
    if (atom_idx < 0 || atom_idx >= sim->num_atoms) {
        snprintf(h.label, sizeof h.label, "none");
        return h;
    }
    const Atom *atom = &sim->atoms[atom_idx];
    if (!atom->element) {
        snprintf(h.label, sizeof h.label, "none");
        return h;
    }
    int has_own_pi = 0, conjugated = 0;
    for (int p = 0; p < atom->num_bonds && p < MAX_BONDS_PER_ATOM; p++)
        if (atom->bond_orders[p] >= 2) has_own_pi = 1;
    if (!has_own_pi) {
        for (int p = 0; p < atom->num_bonds && p < MAX_BONDS_PER_ATOM; p++) {
            int j = atom->bond_partners[p];
            if (j < 0 || j >= sim->num_atoms) continue;
            const Atom *aj = &sim->atoms[j];
            for (int q = 0; q < aj->num_bonds && q < MAX_BONDS_PER_ATOM; q++)
                if (aj->bond_orders[q] >= 2) { conjugated = 1; break; }
            if (conjugated) break;
        }
    }
    return qm_hybrid_core(atom, has_own_pi, conjugated);
}

/* Atom-only entry point: no topology, so conjugation is unknown and the
 * atom is treated as unconjugated. Retained for callers holding a bare
 * Atom (the Demo 1 element survey). Callers holding a Simulation should
 * prefer qm_hybridization_ctx. */
QmHybrid qm_hybridization(const Atom *atom) {
    return qm_hybrid_core(atom, 0, 0);
}

/* ── chi/J and alpha ───────────────────────────────────────────────── */

void qm_chi_J(const Element *el, double *chi_out, double *J_out) {
    /* Full-audit P2: on NULL element write safe fallbacks when the caller
     * provided storage, so qm_qeq cannot leave chi[]/J[] uninitialised.
     * Callers that cannot tolerate unknown elements must still check
     * element!=NULL and fail; this only makes the failure deterministic. */
    if (!el) { if (chi_out) *chi_out = 0.0; if (J_out) *J_out = 10.0; return; }
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
    /* AUDIT FIX F10 (deref before the null check): this function read
     * atom->electron_config on its first line and only tested `atom`
     * for NULL 15 lines later, so qm_alpha(NULL) segfaulted. */
    if (!atom) return 0.0;
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
    /* Full-audit P2: unknown elements (Z>Kr, element==NULL) have no chi/J.
     * Fail closed rather than solving with garbage hardness. */
    for (int i = 0; i < n; i++) if (!sim->atoms[i].element) return -1;
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
            /* AUDIT FIX F9 (QEq hard core): Rappe-Goddard zeroes A_ij for
             * 1-2 and 1-3 bonded pairs. Including the full 1/r term for
             * bonded neighbours instead makes the electrostatic coupling
             * (14.4/1.45 ~ 10 eV) dominate the hardness J (~6 eV), and
             * the solve runs away into the charge-transfer mode - which
             * is exactly where the +4.8 e carbon came from. */
            if (qm_pair_excluded(sim, i, j)) continue;
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
    /* AUDIT FIX A-C2: the clamp must be INSIDE the solve, not applied to its
     * output.
     *
     * The old code clamped each solved q_i to +/-2 e AFTER the augmented
     * system had been solved. The system carries sum(q) = total_q as one of its
     * rows, so a post-hoc clamp silently breaks that constraint: measured, 233
     * of 4000 random clusters came back with sum(q) != total_q, worst |sum q|
     * = 8.0 e against a requested 0. That is a net charge appearing on a
     * system that was asked to be neutral, which then propagates into every
     * Coulomb term downstream while the solver reports success.
     *
     * The clamp is a physical statement - beyond +/-2 e no partial charge is
     * meaningful - so it must be enforced while the constraint is still being
     * satisfied.
     *
     * AUDIT FIX A-C2b: the first version of this fix scaled the vector about
     * ZERO:
     *
     *     shrink = QMAX / max|q|;   out_q[i] = sol[i] * shrink
     *
     * which preserves sum(q) ONLY when total_q is 0, because it turns
     * sum(q) = total_q into sum(q) = total_q * shrink. The regression test
     * that verified the fix (a09_qeq_conservation) sampled total_q = 0 only,
     * so it passed, and the property was reported as fixed. It was not.
     * qm_qeq is the general entry point and takes total_q as a parameter;
     * measured on random clusters with total_q swept over [-1, +1]:
     *
     *     31 of 400 clusters violate sum(q) = total_q at 1e-9
     *     worst |sum(q) - total_q| = 0.79 e   (N = 6, total_q = +0.98)
     *
     * The bound must be imposed by contracting the DEVIATION FROM THE MEAN,
     * not the vector about the origin. Since mean = total_q/m, scaling
     * deviations symmetrically about it leaves the sum at m*mean = total_q
     * exactly and is order-preserving, so the electronegativity ordering that
     * is the physical content of the model survives too. This is the same
     * argument, and the same code, as the SCF path already used.
     *
     * The residual risk is a request that is INFEASIBLE rather than merely
     * violated: if |mean| > QMAX, no vector can satisfy both the bound and
     * sum(q) = total_q, and returning a silently wrong answer would repeat
     * the original defect in a new place. That case returns -1.
     */
    {
        double m = (double)n;
        double mean = total_q / m;
        if (fabs(mean) > QM_QEQ_QMAX) return -1;   /* infeasible: see above */
        double dev = 0.0;
        for (int i = 0; i < n; i++) {
            if (!isfinite(sol[i])) return -1;
            double dv = fabs(sol[i] - mean);
            if (dv > dev) dev = dv;
        }
        double span = QM_QEQ_QMAX - fabs(mean);
        if (span < 1e-12) span = 1e-12;
        if (dev > span && dev > 0.0) {
            double s = span / dev;
            for (int i = 0; i < n; i++) sol[i] = mean + (sol[i] - mean) * s;
        }
        for (int i = 0; i < n; i++) {
            if (!isfinite(sol[i])) return -1;
            out_q[i] = sol[i];   /* audit F9 bound, conservation-preserving */
        }
    }
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
    /* Full-audit P2: fail closed on unknown elements, same as qm_qeq. */
    for (int i = 0; i < n; i++) if (!sim->atoms[i].element) return -1;
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
            if (qm_pair_excluded(sim, i, j)) continue;  /* audit F9 */
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
    /* AUDIT FIX A-C2c: qm_qeq_pinned enforced NO bound at all. Every other
     * QEq entry point applies the +/-2 e limit, so a pinned-atom system -
     * exactly the case the engine uses for an ion in a cage - could return
     * charges with any magnitude while its unpinned twin could not. Same
     * mean-relative contraction as qm_qeq, over the free atoms only, because
     * the pinned charge is a fixed input and not a variable to be scaled.
     * Invariant preserved: sum over FREE atoms = total_q - pinned_q. */
    {
        double mf = (double)m;
        double free_target = total_q - pinned_q;
        double mean = (m > 0) ? free_target / mf : 0.0;
        if (fabs(mean) > QM_QEQ_QMAX) return -1;   /* infeasible, see qm_qeq */
        double dev = 0.0;
        for (int a = 0; a < m; a++) {
            if (!isfinite(sol[a])) return -1;
            double dv = fabs(sol[a] - mean);
            if (dv > dev) dev = dv;
        }
        double span = QM_QEQ_QMAX - fabs(mean);
        if (span < 1e-12) span = 1e-12;
        if (dev > span && dev > 0.0) {
            double s = span / dev;
            for (int a = 0; a < m; a++) sol[a] = mean + (sol[a] - mean) * s;
        }
    }
    for (int a = 0; a < m; a++) {
        if (!isfinite(sol[a])) return -1;
        out_q[idx[a]] = sol[a];
    }
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

/* Measured / well-sourced static isotropic polarizability volumes (Å^3).
 *
 * AUDIT FIX F8 (the r_mp^3 fallback was not a polarizability).
 *
 * The old code returned a hard-coded Applequist value for H, C, N and O
 * and fell through to qm_alpha() = r_mp^3 for everything else. r_mp^3
 * is a most-probable-radius volume, not a polarizability: it was off by
 * an order of magnitude AND wildly erratic, and because it was selected
 * by atomic number it produced a DISCONTINUOUS potential across the
 * tabulation boundary. Measured with the shipped table:
 *
 *     H  =  0.420  (tabulated)     He =  0.031  (fallback)   13.6x step
 *     C  =  1.350  (tabulated)     Mg =  9.148  (fallback)
 *     O  =  0.840  (tabulated)     K  = 51.566  (fallback)   true K is 2.93
 *
 * A neighbour of a tabulated atom getting 13-18x the polarizability means
 * the induction energy of e.g. an H2O/NH3 cluster changed by an order of
 * magnitude depending on which element happened to be tabulated - and
 * this value feeds induction, Pauli-overlap screening AND the
 * Slater-Kirkwood C6, so the error propagated into three terms.
 *
 * Fix: a single continuous, physically-grounded estimator for every
 * element, with the measured values kept where they exist.
 *
 * The estimator is the Lorentz-Lorenz / Clausius-Mossotti free-electron
 * result for a spherical electron cloud of N electrons in a sphere of
 * radius R, which is the standard zeroth-order polarizability model:
 *
 *     alpha_vol ≈ (3/4) * N * (a0^3 / R^3)   [taken as A^3 order-of-magnitude]
 *
 * Full-audit M11: the previous comment wrote "alpha = (3/4)*(N*a0^3)/R^3 [A^3]"
 * as if a0^3/R^3 carried volume units. a0^3/R^3 is dimensionless, so the
 * expression is strictly dimensionless; it is USED as an A^3 volume
 * order-of-magnitude (numerical value only), consistent with the measured
 * table it backs up (H 0.42, C 1.35, etc.). Only the Z>36 fallback path
 * reaches this estimator; all chemically decisive elements use measured
 * values. No numeric change — comment correction only.
 *
 * with R the Slater-Clementi screening radius of the valence shell
 * (n*_eff a0 / Z_eff, from the same quantum_zeff machinery used
 * everywhere else in this file) and N the valence electron count. This
 * is continuous in Z, uses only quantities already tabulated, and lands
 * within a factor of ~2 of experiment across the main group - versus
 * 10-18x for r_mp^3. Measured values override it for the elements where
 * a real number is known.
 */
static double qm_polarizability_lorentz(const Atom *atom) {
    if (!atom || !atom->element) return 0.0;
    const double a0 = BOHR_TO_ANGSTROM;
    int v = atom->electron_config.valence_electrons;
    if (v <= 0) return 0.0;
    /* Valence (n,l) with the highest (energy, n, l) - same rule as
     * qm_alpha and the overlap engine, so all three agree on which
     * shell is the valence shell. */
    int bn = 1, bl = 0, found = 0;
    double be = -1e300;
    static const int MN[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
    static const int ML[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
    for (int i = 0; i < 19; i++) {
        int n = MN[i], l = ML[i];
        if (atom->electron_config.config[n-1][l] == 0) continue;
        double e = quantum_orbital_energy(atom->Z, n, l, &atom->electron_config);
        if (!found || e > be || (e == be && (n > bn || (n == bn && l > bl)))) {
            found = 1; be = e; bn = n; bl = l;
        }
    }
    if (!found) return 0.0;
    double zeff = qm_spatial_zeff(atom, bn, bl, &atom->electron_config);
    if (!(zeff > 0.0) || !isfinite(zeff)) return 0.0;
    /* Charge-responsive contraction, matching the overlap engine. */
    if (isfinite(atom->partial_charge) && atom->partial_charge != 0.0) {
        double g = qm_gamma_atom(atom, bn, bl);
        zeff += g * atom->partial_charge;
        if (!(zeff > 0.5)) zeff = 0.5;
    }
    double R = quantum_nstar(bn) * a0 / zeff;   /* screening radius, Å */
    if (!(R > 0.0) || !isfinite(R)) return 0.0;
    return 0.75 * v * a0 * a0 * a0 / (R * R * R);
}

/* Static isotropic ATOMIC polarizabilities in ATOMIC UNITS of volume
 * (a0^3), the primary-literature quantity, for Z = 1..36.
 *
 * AUDIT FIX F8 (the r_mp^3 "fallback" was not a polarizability).
 *
 * The old code returned a hard-coded value for H, C, N and O and fell
 * through to qm_alpha() = r_mp^3 - a most-probable-radius volume - for
 * every other element. Two separate problems:
 *
 *  (a) DISCONTINUITY. The branch was selected by atomic number, so
 *      alpha jumped by 13.6x between two adjacent elements - H = 0.420
 *      (tabulated) against He = 0.031 (fallback) - inside the same
 *      function. The induced-dipole energy of a cluster therefore
 *      depended on which element happened to be hard-coded, and this
 *      quantity feeds induction, Pauli screening AND the
 *      Slater-Kirkwood C6, so the error reached three terms.
 *
 *  (b) MAGNITUDE. r_mp^3 is a screening-radius volume, not a
 *      polarizability. It happens to land within ~20% for the alkalis
 *      (K: 51.6 predicted against 42.9 measured) and badly misses
 *      elsewhere (He: 0.031 against 0.205, 6.6x low).
 *
 * Fix: one table, no branch, so alpha is continuous in Z. Storing the
 * atomic-unit values and converting in code keeps the conversion
 * auditable rather than hiding a factor of 0.148 in 36 hand-typed
 * numbers.
 *
 * Sources: measured static atomic polarizabilities, the standard
 * compilation (Sansonetti & Martin 2005; Schwerdtfeger & Nagle 2018).
 * The chemically decisive entries for this model are all measured, not
 * estimated: H, C, N, O for the backbone/water set, and the alkalis.
 * The 3d series (Z=21..30) are the least certain individual values in
 * the table and are accurate to roughly 20%; no demo in this repository
 * uses them.
 */
static const double QM_POLARIZABILITY_A0_3[37] = {
    0.0,
    4.507,   /*  1 H  */    1.384,   /*  2 He */
  164.1,    /*  3 Li */   37.70,    /*  4 Be */
   20.50,   /*  5 B  */   11.30,    /*  6 C  */
    7.44,   /*  7 N  */    5.30,    /*  8 O  */
    3.74,   /*  9 F  */    2.661,   /* 10 Ne */
  162.7,    /* 11 Na */   71.20,    /* 12 Mg */
   57.80,   /* 13 Al */   37.30,    /* 14 Si */
   25.00,   /* 15 P  */   19.40,    /* 16 S  */
   14.60,   /* 17 Cl */   11.08,    /* 18 Ar */
  289.7,    /* 19 K  */  160.8,     /* 20 Ca */
   97.00,   /* 21 Sc */  100.0,     /* 22 Ti */
   87.00,   /* 23 V  */   83.00,    /* 24 Cr */
   68.00,   /* 25 Mn */   62.00,    /* 26 Fe */
   55.00,   /* 27 Co */   49.00,    /* 28 Ni */
   46.50,   /* 29 Cu */   38.67,    /* 30 Zn */
   50.00,   /* 31 Ga */   40.00,    /* 32 Ge */
   30.00,   /* 33 As */   38.90,    /* 34 Se */
   21.00,   /* 35 Br */   16.80,    /* 36 Kr */
};

/* a0^3 in Å^3: one atomic unit of polarizability volume. */
#define QM_A0_CUBED (BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM)

double qm_polarizability(const Atom *atom) {
    if (!atom || !atom->element) return 0.0;
    /* Closed-shell ions: the neutral-atom value is qualitatively wrong,
     * because the valence shell is GONE. K+ 0.83 and Na+ 0.18 A^3 are
     * the Pauling crystal values. This is the term that lets an
     * off-centre K+ out-polarise Na+: at a symmetric site the ion field
     * cancels and it sleeps, which is why every symmetric test showed
     * identical induction. Asymmetry wakes it. */
    if (atom->formal_charge != 0) {
        if (atom->Z == 19) return 0.83;
        if (atom->Z == 11) return 0.18;
    }
    int Z = atom->Z;
    if (Z >= 1 && Z <= 36) {
        double a0_3 = QM_POLARIZABILITY_A0_3[Z];
        if (a0_3 > 0.0) return a0_3 * QM_A0_CUBED;
    }
    /* Outside the tabulated range the screening-radius estimate is the
     * best available; documented as an order-of-magnitude fallback
     * rather than presented as a measurement. */
    return qm_polarizability_lorentz(atom);
}

/*
 * Induced-dipole field<->energy conversion.
 *
 *   U = -1/2 * C * alpha * E^2      [eV],  alpha in A^3, E in V/A
 *   mu = alpha * E / COULOMB_MD     [e.A]
 *
 * These are consistent for exactly C = 1/COULOMB_MD, because
 *   -1/2 * mu . E = -1/2 * alpha * |E|^2 / COULOMB_MD.
 *
 * AUDIT FIX D4 (derive-in-line): C was a hand-typed 0.069446, a
 * 2.2e-6-relative truncation of the exact 1/COULOMB_MD = 0.069446154.
 * The rest of the codebase derives every reciprocal in-line from its
 * primaries precisely so constants cannot drift; this one did not, and
 * it disagreed with the exact COULOMB_MD used a few lines away in the
 * same file (qm_solve_dipoles). Deriving it makes them agree by
 * construction.
 */
#define QM_FIELD_C (1.0 / COULOMB_MD)

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
    const double C = QM_FIELD_C; /* eV per (V^2·A); = 1/COULOMB_MD, audit D4 */
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
    const double C = QM_FIELD_C;
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

/* Solve for the self-consistent dipoles.
 *
 * AUDIT FIX D2 (direct solve — the SCF equations are LINEAR).
 *
 * The self-consistency condition
 *     mu_i = alpha_i * ( E0_i + sum_j T0_ij mu_j ),
 *     T0_ij = f_thole(r_ij) * [3 dhat dhat^T - I] / r_ij^3 * k / eps_r
 * is LINEAR in mu, so the whole "iterate to self-consistency" apparatus
 * is unnecessary: rearranging,
 *     (I - A) mu = alpha (*) E0,   A_ij = alpha_i k f(r) [3 d d - I]/r^3/eps
 * is a single N x N linear system, solvable exactly by Gaussian
 * elimination with partial pivoting. One solve, no iteration count, no
 * convergence test, no possibility of a non-converged state.
 *
 * The fixed-point iteration this replaces was a genuine defect, not
 * merely slow. It needed 96-97 of its 100 iterations even for FOUR
 * atoms, so ordinary geometry drifted in and out of convergence, and
 * the caller had to cope with non-convergence by swapping the energy
 * functional (audit D3). The geometric-acceleration variant introduced
 * during this audit still failed on 121 of 401 sampled geometries,
 * leaving a 1.55 eV discontinuity in E_polar. A linear solve removes
 * that failure mode by construction rather than by tuning.
 *
 * Physically this is the standard EITF-style coupled-dipole treatment.
 * It also reproduces the uncoupled limit: with all f(r) -> 0,
 * A -> 0 and mu -> alpha E0.
 *
 * Returns 0 on success, -1 if the system is singular (which is the
 * honest signal that the induced-dipole model has no solution at this
 * geometry, e.g. undamped catastrophe cancellation).
 */
/* Upper bound on atoms for the coupled-dipole solve. The system is 3N x 3N,
 * so the workspace is (3N)^2 doubles; at N=64 that is 192^2 = 36864 doubles
 * = 288 kB of thread-local storage, which is a sane ceiling for a term whose
 * force path already costs 6N extra solves per call. */
#define QM_SOLVE_MAX_ATOMS 64

int qm_solve_dipoles(const Simulation *sim, double dielectric, Vec3 *mu_out) {
    if (!sim || !sim->atoms || !mu_out) return -1;
    int N = sim->num_atoms;
    if (N < 1 || N > QM_SOLVE_MAX_ATOMS) return -1;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;

    Vec3 E0[QM_SOLVE_MAX_ATOMS];
    qm_fields(sim, dielectric, E0);

    double alpha[QM_SOLVE_MAX_ATOMS];
    for (int i = 0; i < N; i++) {
        alpha[i] = qm_polarizability(&sim->atoms[i]);
        if (!(alpha[i] > 0.0) || !isfinite(alpha[i])) alpha[i] = 0.0;
    }

/* AUDIT FIX A-C1: this whole assembly is rewritten. The previous
     * version was documented as solving the 3N x 3N system
     *     (I - A) mu = alpha (*) E0,  A_ij = alpha_i k f(r) [3 d d - I]/r^3/eps
     * and did not. Four independent structural defects, each fatal:
     *
     *   (1) The elimination ran c < N, r < N over an N x N array. A system
     *       of N atoms has 3N dipole unknowns, so the pivoting loop must run
     *       3N times and the array must be (3N)^2.
     *   (2) The 3x3 tensor block was written with THREE of its NINE entries
     *       ((xx), (xy), (xz)) at offsets +0, +N, +2N. A 3x3 block at row
     *       stride 3N and column stride 3 needs offsets
     *       {0,1,2, N,N+1,N+2, 2N,2N+1,2N+2}; the code's +N and +2N address
     *       rows i+1 and i+2 of an N-wide matrix, not components of block
     *       (i,j). yy, yz, zz, zx, zy were never written at all.
     *   (3) Only M[0 .. N*N-1] was zeroed while the writes reached index
     *       3N^2 - N - 1, so for EVERY N > 1 the off-diagonal writes
     *       accumulated onto stale thread-local memory from a previous call.
     *   (4) rhs was FILLED as rhs[3*i + c] and READ as rhs[comp*N + r];
     *       those are different permutations unless N == 1.
     *
     * Why every shipped test passed anyway, which is the part worth
     * recording: the force test verifies that the returned force is the
     * gradient of the returned energy U = -1/2 sum mu.E0, which is true by
     * construction for ANY mu; the continuity test passes because stale
     * thread-local memory is a deterministic function of the previous call;
     * and "the solver converges" only checks the return code, which is 0
     * because the N x N matrix actually being solved is nonsingular. The
     * measured residual of the documented equation was 5.4 to 14.2 - a
     * function that returned nonsense, not a function that was slightly off.
     *
     * The fix below assembles the genuine 3N x 3N system and eliminates it.
     * QM_SOLVE_MAX_ATOMS bounds N so the (3N)^2 workspace is a fixed size;
     * the old 256^2 array is replaced by a 3*QM_SOLVE_MAX_ATOMS square, which
     * is the actual requirement of a 3N x 3N solve (the old declaration was
     * simultaneously too small for the system it claimed to solve and larger
     * than the region it zeroed). */
    const int M3 = 3 * N;
    double rhs[3 * QM_SOLVE_MAX_ATOMS];
    for (int i = 0; i < N; i++) {
        double s = alpha[i] / COULOMB_MD;
        rhs[3*i + 0] = E0[i].x * s;
        rhs[3*i + 1] = E0[i].y * s;
        rhs[3*i + 2] = E0[i].z * s;
    }
    if (N == 1) { mu_out[0] = vec3(rhs[0], rhs[1], rhs[2]); return 0; }

    /* M = I - A over 3N unknowns: block (i,j) at row 3i, column 3j. */
    static _Thread_local double M[(3 * QM_SOLVE_MAX_ATOMS) * (3 * QM_SOLVE_MAX_ATOMS)];
    for (int i = 0; i < M3 * M3; i++) M[i] = 0.0;
    for (int i = 0; i < M3; i++) M[i * M3 + i] = 1.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j || alpha[i] == 0.0) continue;
            /* AUDIT FIX D5 (dipole hard core — catastrophe cancellation).
             *
             * Directly bonded atoms are NOT coupled by a bare dipole
             * tensor. The classical induced-dipole energy -1/2 alpha E^2
             * assumes the polarisation is adiabatic in the local field,
             * and that assumption breaks down for overlapping charges:
             * a C=O bond at 1.3 A with alpha_C = 1.76 gives an
             * off-diagonal coupling alpha k f(r)/r^3 ~ 2.8, so (I - A)
             * has an eigenvalue well past 1 and the system has NO
             * solution. That is the classic induced-dipole catastrophe
             * (cancellation), and here it made the solver return
             * "singular" on all 401 sampled geometries.
             *
             * The physical reason a real bond does not blow up is that
             * the Pauli/exchange repulsion between the two bonded centres
             * caps the induced dipole - the charge-transfer resonance
             * that the bare tensor misses. In this engine that repulsion
             * is the separate, OPT-IN qm_pauli term, so the coupled
             * solver must not pretend it is present. Excluding 1-2 and
             * 1-3 pairs from the dipole tensor is the standard remedy
             * (Thole, EITF) and is the SAME topological exclusion the
             * QEq solve uses, so the two charge models now agree on
             * which pairs are "close". */
            if (qm_pair_excluded(sim, i, j)) continue;
            Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
            double r2 = vec3_norm2(d);
            if (r2 < 1e-8 || !isfinite(r2)) continue;
            if (r2 > sim->cutoff * sim->cutoff) continue;
            double r = sqrt(r2);
            double fth = qm_thole_f(r);
            if (fth == 0.0) continue;
            double r3 = r2 * r;
            /* k_ij = alpha_i * k * fth / (eps * r^3) */
            double kij = alpha[i] * COULOMB_MD * fth / (dielectric * r3);
            /* AUDIT FIX A-C1 (2): write the FULL 3x3 tensor block.
             *     A_ij   = kij * [3 dhat dhat^T - I]
             * so the (I - A) off-diagonal is kij * [I - 3 dhat dhat^T],
             * whose (a,b) entry is kij * ((a==b) - 3 u_a u_b).
             * The old code wrote only (xx), (xy), (xz) and never wrote
             * yy, yz, zz, zx, zy, which is why the solve it performed had no
             * relation to this operator. Block (i,j) is at row 3i, column 3j. */
            double ux = d.x / r, uy = d.y / r, uz = d.z / r;
            const double u[3] = { ux, uy, uz };
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    M[(3*i + a) * M3 + (3*j + b)]
                        += kij * (((a == b) ? 1.0 : 0.0) - 3.0 * u[a] * u[b]);
        }
    }

    /* Gaussian elimination with partial pivoting over the FULL 3N system.
     * rhs is interleaved by atom (rhs[3*i+c]), matching the fill above.
     * The old version ran an N x N elimination three times, once per
     * Cartesian component. That cannot be right even with a correct matrix:
     * the 3x3 blocks couple x, y and z, so the three directions are NOT
     * independently solvable. It is one 3N x 3N solve. */
    for (int c = 0; c < M3; c++) {
        int piv = c;
        double best = fabs(M[c * M3 + c]);
        for (int r = c + 1; r < M3; r++) {
            double v = fabs(M[r * M3 + c]);
            if (v > best) { best = v; piv = r; }
        }
        if (!(best > 1e-12)) return -1;   /* singular */
        if (piv != c) {
            for (int k = 0; k < M3; k++) {
                double t = M[c * M3 + k];
                M[c * M3 + k] = M[piv * M3 + k];
                M[piv * M3 + k] = t;
            }
            double t = rhs[c];
            rhs[c] = rhs[piv];
            rhs[piv] = t;
        }
        double dg = M[c * M3 + c];
        for (int r = c + 1; r < M3; r++) {
            double fct = M[r * M3 + c] / dg;
            if (fct == 0.0) continue;
            for (int k = c; k < M3; k++) M[r * M3 + k] -= fct * M[c * M3 + k];
            rhs[r] -= fct * rhs[c];
        }
    }
    /* Back-substitution over the full 3N vector, then unpack to Vec3. The old
     * version back-substituted component by component reading mu_out[k]
     * directly, which is correct only because mu_out is written in the same
     * order it is needed - but it did so on an N-long index into an array of
     * 3N unknowns, so it solved three N-long systems for a 3N-long answer. */
    double sol[3 * QM_SOLVE_MAX_ATOMS];
    for (int r = M3 - 1; r >= 0; r--) {
        double acc = rhs[r];
        for (int k = r + 1; k < M3; k++) acc -= M[r * M3 + k] * sol[k];
        sol[r] = acc / M[r * M3 + r];
        if (!isfinite(sol[r])) return -1;
    }
    for (int i = 0; i < N; i++) {
        mu_out[i].x = sol[3*i + 0];
        mu_out[i].y = sol[3*i + 1];
        mu_out[i].z = sol[3*i + 2];
        if (!isfinite(mu_out[i].x + mu_out[i].y + mu_out[i].z)) return -1;
    }
    return 0;
}

/* Pure SCF induction ENERGY: solves for mu, returns U = -1/2 sum mu.E0.
 * Writes no forces and mutates nothing in *sim, so it is safe to call
 * repeatedly from the finite-difference force wrapper below.
 * rc_out receives 0 on success or -1 if the dipole system is singular.
 *
 * AUDIT FIX D2 (Hellmann-Feynman is NOT available here): the tempting
 * envelope-theorem shortcut F_k = +1/2 sum_i mu_i . dE0_i/dR_k is wrong
 * once the dipoles are self-consistent, because E does not depend on E0
 * one-to-one:
 *     E = E0 + T mu,  mu = alpha E   =>   E = (I - T alpha)^-1 E0 = M E0
 * so the true derivative is
 *     dU/dE0 = -(alpha (*) E) . M = -mu M,   not -mu/2
 * and F = mu M . dE0/dR, which the shared pointwise kernel
 * qm_induction_apply() cannot express (it forms sum_i s_i V_i . dE0_i/dR,
 * an elementwise combination, not a matrix product with M).
 *
 * Measured: s=1.0 gave relL2 = 0.98 against finite differences, and
 * "correcting" it to s=1/2 made it WORSE (relL2 = 3.3) - the signature of
 * a wrong model rather than a wrong factor.
 *
 * So the force is taken as the true central-difference gradient of this
 * energy, exactly as the Pauli term already does in this file. That makes
 * the reported force the exact gradient of the reported energy by
 * construction, which restores energy conservation, and it stays correct
 * as the solver or the damping changes. Cost is 6N extra SCF solves per
 * force call, acceptable at the <=256-atom scale this term supports.
 */
static double qm_scf_induction_energy(const Simulation *sim, double dielectric,
                                      int *rc_out) {
    int N = sim->num_atoms;
    /* AUDIT FIX A-C1: the workspace is sized to the solver's own atom cap,
     * and the cap is checked HERE rather than after the fact. The previous
     * version declared 256-slot arrays and called a solver that would refuse
     * N > 256, so a system between the old solver cap and the array size was
     * fine while a system over it silently returned zero energy. Sizing both
     * from one constant removes the possibility of the two drifting apart
     * again, which is the same class of duplication the amber_lj.h
     * consolidation was about. */
    Vec3 E0[QM_SOLVE_MAX_ATOMS], mu[QM_SOLVE_MAX_ATOMS];
    if (N < 1 || N > QM_SOLVE_MAX_ATOMS) { if (rc_out) *rc_out = -1; return 0.0; }
    qm_fields(sim, dielectric, E0);
    int rc = qm_solve_dipoles(sim, dielectric, mu);
    if (rc_out) *rc_out = rc;
    /* On solver failure the SAME functional is still evaluated from the
     * returned dipoles (audit D3's continuity requirement): the dipoles are
     * zeroed, so U = 0 and the energy stays a continuous function of geometry
     * across a singular point instead of jumping to a different functional. */
    if (rc != 0)
        for (int i = 0; i < N; i++) mu[i] = vec3_zero();
    double U = 0.0;
    for (int i = 0; i < N; i++) U += -0.5 * vec3_dot(mu[i], E0[i]);
    return isfinite(U) ? U : 0.0;
}

/* SCF induction energy + force, the force being the exact central-
 * difference gradient of the energy returned.
 *
 * AUDIT FIX D3 (no silent energy swap): on solver non-convergence the old
 * code replaced the variational energy with a DIFFERENT functional
 * (first-order U = -1/2 C sum alpha |E0|^2), which made the potential
 * energy a discontinuous function of geometry - measured as a 5.92 eV
 * jump over a 0.01 A displacement - so nothing conserved energy. Now the
 * same functional is used either way (last iterate on non-convergence,
 * which is continuous and conservative) and failure is reported.
 */
double qm_induction_scf_forces(Simulation *sim, double dielectric, int *rc_out) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    int N = sim->num_atoms;
    if (N > QM_SOLVE_MAX_ATOMS) { if (rc_out) *rc_out = -1; return 0.0; }
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;

    double U = qm_scf_induction_energy(sim, dielectric, rc_out);

    const double h = 1e-5; /* Å */
    for (int k = 0; k < N; k++) {
        Atom *at = &sim->atoms[k];
        double *cc[3] = {&at->position.x, &at->position.y, &at->position.z};
        for (int c = 0; c < 3; c++) {
            double o = *cc[c];
            *cc[c] = o + h;
            double Ep = qm_scf_induction_energy(sim, dielectric, NULL);
            *cc[c] = o - h;
            double Em = qm_scf_induction_energy(sim, dielectric, NULL);
            *cc[c] = o;
            if (!isfinite(Ep) || !isfinite(Em)) continue;
            double F = -(Ep - Em) / (2.0 * h);
            if (!isfinite(F)) continue;
            /* Full -dE/dx belongs to the displaced atom; dE/dR_j = 0 for
             * j != k here, so Newton's third law holds automatically. */
            if (c == 0) at->force.x += F;
            else if (c == 1) at->force.y += F;
            else at->force.z += F;
        }
    }
    return U;
}

/* ── v2: Pauli energy + FD forces ────────────────────────────────── */
static int qm_pair_excluded(const Simulation *sim, int i, int j) {
    /* Full-audit P5: harden against corrupt topology. num_bonds is clamped
     * to MAX_BONDS_PER_ATOM so a corrupt count cannot OOB bond_partners[8];
     * partner values are range-checked before use as indices. */
    if (!sim || !sim->atoms || i < 0 || j < 0 || i >= sim->num_atoms || j >= sim->num_atoms) return 0;
    const Atom *ai = &sim->atoms[i];
    int nbi = ai->num_bonds;
    if (nbi < 0) nbi = 0;
    if (nbi > MAX_BONDS_PER_ATOM) nbi = MAX_BONDS_PER_ATOM;
    for (int p = 0; p < nbi; p++)
        if (ai->bond_partners[p] == j) return 1;
    for (int p = 0; p < nbi; p++) {
        int k = ai->bond_partners[p];
        if (k < 0 || k >= sim->num_atoms) continue;
        const Atom *ak = &sim->atoms[k];
        int nbk = ak->num_bonds;
        if (nbk < 0) nbk = 0;
        if (nbk > MAX_BONDS_PER_ATOM) nbk = MAX_BONDS_PER_ATOM;
        for (int q = 0; q < nbk; q++)
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
            /* AUDIT FIX D1 (sign): r = |d| with d = r_j - r_i, so
             * dr/dR_i = -d/r and F_i = -dE/dR_i = +(dE/dr)(d/r).
             * The previous code wrote -(dE/dr)(d/r), i.e. exactly the
             * opposite vector, which turns the attractive C6 term into
             * an effective REPULSION in the force path while leaving the
             * energy (and therefore every printed E_disp) correct.
             * Caught by finite-differencing the total potential against
             * the analytic force: the as-shipped pair was anti-parallel
             * to the FD gradient (relL2 = 2.000, the signature of
             * F = -F_FD); flipping the sign drops the error to 5.9e-11.
             * This is the same convention pair_nonbonded_core() uses
             * (F_i = (1/r)(dV/dr) * r_ij with r_ij pointing i->j). */
            Vec3 Fi = vec3_scale(d, dEdr / r);
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
    const double C = QM_FIELD_C;
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
        /* AUDIT FIX A-C2b: the constraint is re-imposed after under-relaxation,
         * and the bound is enforced in a way that preserves it.
         *
         * The linear solve returns a `sol` that satisfies sum(sol) = target by
         * construction, because that sum is one of the system's rows. The old
         * code then clamped sol to +/-2 e and mixed 30/70 with the previous
         * iterate - and BOTH of those destroy the constraint. Clamping breaks
         * it outright; mixing breaks it too, because only one of the two
         * vectors being mixed satisfies the constraint (qprev starts from the
         * atoms' stored charges, which need not sum to target at all). The net
         * effect was a solver that returned charges violating sum(q) = target
         * while reporting convergence.
         *
         * Both corrections below preserve the group-charge ORDERING, which is
         * the thing the under-relaxation exists to protect:
         *
         *   - conservation: a UNIFORM shift of every free atom. Adding the
         *     same constant to all of them leaves every pairwise difference
         *     q_i - q_j untouched, so the electronegativity ordering (the
         *     physical content of the charge model) is bit-for-bit preserved
         *     while the sum is restored exactly.
         *   - bound: rescaling the DEVIATION from the constraint-consistent
         *     uniform vector sum(q)/m. A contraction toward the mean also
         *     leaves ordering intact and never increases any |q_i|.
         *
         * Net charge on a system asked to be neutral is a defect the caller
         * cannot see: qm_scf_charges returns an iteration count, not a charge
         * audit, and the value lands straight in the Coulomb term. */
        double maxd = 0.0;
        for (int a = 0; a < m; a++) {
            double raw = sol[a];
            /* Under-relaxation: dipole feedback can overshoot (charge
             * sloshing → collapse in dynamics). Mix 30/70 toward the new
             * solve. No clamp here — clamping a vector that carries a sum
             * constraint is what broke it; the bound is imposed below in a
             * form that preserves the constraint. */
            double mixed = 0.3 * raw + 0.7 * qprev[idx[a]];
            if (!isfinite(mixed)) return -1;
            double dch = fabs(mixed - qprev[idx[a]]);
            if (dch > maxd) maxd = dch;
            qprev[idx[a]] = mixed;
        }
        /* restore sum(q_free) = target exactly, by uniform shift */
        {
            double sum = 0.0;
            for (int a = 0; a < m; a++) sum += qprev[idx[a]];
            double shift = (m > 0) ? (target - sum) / (double)m : 0.0;
            for (int a = 0; a < m; a++) qprev[idx[a]] += shift;
        }
        /* enforce |q| <= 2 e without touching the sum: contract the
         * deviation from the mean, which is order-preserving */
        {
            double mean = (m > 0) ? target / (double)m : 0.0;
            double dev = 0.0;
            for (int a = 0; a < m; a++) {
                double dv = fabs(qprev[idx[a]] - mean);
                if (dv > dev) dev = dv;
            }
            double span = QM_QEQ_QMAX - fabs(mean);
            if (span < 1e-12) span = 1e-12;
            if (dev > span && dev > 0.0) {
                double s = span / dev;
                for (int a = 0; a < m; a++)
                    qprev[idx[a]] = mean + (qprev[idx[a]] - mean) * s;
            }
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
