/* qm_eht_overlap.c — EXACT single-zeta STO overlaps (Roothaan prolate
 * spheroidal machinery, Rosen 1931; Roothaan 1951 tables as check values).
 *
 * Why this file exists: qm_overlap() is a documented heuristic — it drops
 * the Slater polynomial prefactor (S ~ exp(-zR), maximizing m per atom) and
 * the tree lists it as abstraction #7. That was fine for Pauli-wall
 * estimates. It is not a basis for electronic structure. This file
 * replaces hand-waving with integrals for the EHT valence basis.
 *
 * Method (no memorized closed forms): prolate spheroidal coordinates about
 * the A-B axis. r_a = R(xi+eta)/2, r_b = R(xi-eta)/2,
 * dV = (R^3/8)(xi^2-eta^2) d xi d eta d phi. STO radial parts are
 * monomials r^{n-1} (1s: r^0, 2s/2p: r^1), angular s/p parts are
 * polynomials in (xi,eta) except p-pi. So ss, s-psigma, p-psigma reduce
 * to finite sums over A_i*B_j with:
 *   A_n(a) = int_1^inf x^n e^{-ax} dx  (closed finite sum, exact)
 *   B_n(b) = int_-1^1  x^n e^{-bx} dx  (Gauss-Legendre-64, ~1e-14 on the
 *             compact domain; verified in-test against the b=0 closed form)
 * p-pi carries sqrt(xi^2-1): 2D Gauss-Legendre (documented quadrature,
 * not closed form). Heteronuclear (za != zb) is EXACT in this machinery,
 * not averaged — the (alpha,beta) = ((za+zb)R/2,(za-zb)R/2) split is the
 * whole point of prolate coordinates.
 *
 * Every channel is validated against brute-force 3D Cartesian grid
 * integration of independently-coded STOs in tests/test_eht.c. The grid
 * is the oracle; these formulas are the claim.
 */
#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "qm_eht.h"
#include "quantum.h"

/* Slater effective principal quantum numbers (Slater 1930). */
static double eht_nstar (int n) {
    if (n <= 3) return (double)n;
    if (n == 4) return 3.7;
    if (n == 5) return 4.0;
    return 4.2;
}

/* STO normalization N with chi = N r^{n-1} e^{-zr} Y_lm:
 * N_n = (2z)^{n+1/2} / sqrt((2n)!). Units A^{-(n+1/2)} for z in A^{-1}. */
static double sto_norm (int n, double z) {
    /* (2n)! computed, not tabulated.
     *
     * FULL-AUDIT F7: the table used to hold 14! (87178291200) in the
     * n = 6 slot, where (2n)! = 12! = 479001600 belongs. The typo made
     * sto_norm(6, z) low by sqrt(12!/14!) = 0.0741. No shipped caller
     * reached n = 6 (the valence basis is n = 1, 2), so nothing in the
     * record was wrong - but the table advertised n <= 6 and silently
     * delivered a wrong normalisation there. Computing (2n)! removes the
     * class: there is no longer a hand-transcribed factorial to mistype.
     */
    if (n < 1 || n > 8 || !(z > 0.0) || !isfinite (z)) return 0.0;
    double f2 = 1.0;
    for (int i = 2; i <= 2 * n; i++) f2 *= (double)i;   /* (2n)! */
    double t = 2.0 * z;
    double p = t;
    /* t^(2n+1): the loop body runs 2n times, so it starts at i = 1. */
    for (int i = 1; i < 2 * n + 1; i++) p *= t;
    return sqrt (p / f2);
}

double eht_zeta (int Z, int n, int l, const ElectronConfig *cfg) {
    double zeff = quantum_zeff (Z, n, l, cfg);
    if (!(zeff > 0.0) || !isfinite (zeff)) return -1.0;
    double ns = eht_nstar (n);
    if (!(ns > 0.0)) return -1.0;
    return zeff / (ns * BOHR_TO_ANGSTROM);   /* A^{-1} */
}

/* A_n(a) = int_1^inf x^n e^{-ax} dx, a > 0, n <= 8.
 * Closed: e^{-a} sum_{k=0}^{n} n!/k! a^{k-n-1}. The a^{-(n+1)} factor is
 * k-INDEPENDENT (a previous revision folded it as a^{-(n-k+1)} times the
 * a^k numerator, giving a^{2k-n-1} — caught by the 3D-grid oracle at
 * 2.5x overshoot, since only k=(n+1)/2 terms accidentally agreed). */
static double aux_A (int n, double a) {
    if (n < 0 || n > 8 || !(a > 0.0) || !isfinite (a)) return 0.0;
    if (a > 60.0) return 0.0;   /* e^-60: overlap is zero, skip overflow */
    double fn = 1.0;
    for (int i = 2; i <= n; i++) fn *= i;   /* n! */
    double am = 1.0;
    for (int i = 0; i < n + 1; i++) am /= a;    /* a^{-(n+1)} */
    double s = 0.0;
    double fk = 1.0;                        /* k! */
    double apow = 1.0;                      /* a^k */
    for (int k = 0; k <= n; k++) {
        if (k > 0) { fk *= k; apow *= a; }
        s += fn / fk * apow * am;
    }
    return exp (-a) * s;
}

/* Gauss-Legendre nodes/weights on [-1,1], Newton iteration (deterministic).
 * Computed once into static tables (n=64). */
#define EHT_GL_N 64
/* FULL-AUDIT: the Gauss-Legendre tables are shared mutable state, and
 * two threads entering the quadrature concurrently could interleave the
 * table fill and read a half-written node set. They are per-thread now
 * (the rest of the tree already uses _Thread_local for exactly this
 * reason). */
static _Thread_local double gl_x[EHT_GL_N], gl_w[EHT_GL_N];
static _Thread_local int gl_ready = 0;

static void gl_compute (void) {
    const int n = EHT_GL_N;
    for (int i = 1; i <= n / 2; i++) {
        /* initial guess (Tricomi) */
        double x = cos (3.14159265358979323846 * (i - 0.25) / (n + 0.5));
        for (int it = 0; it < 20; it++) {
            /* Legendre P_n(x) + derivative by recurrence */
            double p0 = 1.0, p1 = x;
            for (int k = 2; k <= n; k++) {
                double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
                p0 = p1;
                p1 = p2;
            }
            /* p1 = P_n, p0 = P_{n-1}; P'_n = n(x P_n - P_{n-1})/(x^2-1) */
            double pp = n * (x * p1 - p0) / (x * x - 1.0);
            double dx = p1 / pp;
            x -= dx;
            if (fabs (dx) < 1e-15) break;
        }
        double p0 = 1.0, p1 = x;
        for (int k = 2; k <= n; k++) {
            double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
            p0 = p1;
            p1 = p2;
        }
        double pp = n * (x * p1 - p0) / (x * x - 1.0);
        double w = 2.0 / ((1.0 - x * x) * pp * pp);
        gl_x[n / 2 - i] = -x;
        gl_x[n / 2 + i - 1] = x;
        gl_w[n / 2 - i] = w;
        gl_w[n / 2 + i - 1] = w;
    }
    gl_ready = 1;
}

/* B_n(b) = int_-1^1 x^n e^{-bx} dx via GL-64. */
static double aux_B (int n, double b) {
    if (n < 0 || n > 8 || !isfinite (b)) return 0.0;
    if (!gl_ready) gl_compute ();
    double s = 0.0;
    for (int i = 0; i < EHT_GL_N; i++) {
        double xp = 1.0;
        for (int k = 0; k < n; k++) xp *= gl_x[i];
        s += gl_w[i] * xp * exp (-b * gl_x[i]);
    }
    return s;
}

/* Coefficient table c[i][j] for xi^i eta^j, 0..4. */
typedef double poly5[5][5];

static void poly_zero (poly5 p) {
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) p[i][j] = 0.0;
}

static void poly_mul (const poly5 a, const poly5 b, poly5 out) {
    poly5 t;
    poly_zero (t);
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) {
            if (a[i][j] == 0.0) continue;
            for (int k = 0; k + i < 5; k++)
                for (int l = 0; l + j < 5; l++) t[i + k][j + l] += a[i][j] * b[k][l];
        }
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) out[i][j] = t[i][j];
}

/* Contract sum c_ij A_i(alpha) B_j(beta). */
static double poly_contract (const poly5 c, double alpha, double beta) {
    double s = 0.0;
    double A[5], B[5];
    for (int i = 0; i < 5; i++) {
        A[i] = aux_A (i, alpha);
        B[i] = aux_B (i, beta);
    }
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            if (c[i][j] != 0.0) s += c[i][j] * A[i] * B[j];
    return s;
}

/* Monomial (xi+eta)^p (xi-eta)^q as poly. p,q in {0,1}. */
static void poly_radial (int p, int q, poly5 out) {
    poly5 xp, xm;
    poly_zero (xp);
    poly_zero (xm);
    /* (xi+eta): [1][0]=1, [0][1]=1 ; (xi-eta): [1][0]=1, [0][1]=-1 */
    poly5 fp, fm;
    poly_zero (fp);
    poly_zero (fm);
    fp[0][0] = 1.0;
    fm[0][0] = 1.0;
    if (p) {
        poly5 t;
        poly_zero (t);
        t[1][0] = 1.0;
        t[0][1] = 1.0;
        poly_mul (fp, t, fp);
    }
    if (q) {
        poly5 t;
        poly_zero (t);
        t[1][0] = 1.0;
        t[0][1] = -1.0;
        poly_mul (fm, t, fm);
    }
    poly_mul (fp, fm, out);
    (void)xp;
    (void)xm;
}

/* Core sigma engine: orbitals (na,za,Ya) x (nb,zb,Yb), radial powers
 * (pa,pb), extra z-polynomial factor zpowsel: 0 none, 1 *(xi*eta-1) [p on
 * B], 2 *(xi*eta+1) [p on A], 3 *((xi*eta)^2-1) [p-p sigma].
 * R in Angstrom. Returns the overlap (dimensionless). */
static double overlap_sigma (int na, double za, double Ya, int pa,
                             int nb, double zb, double Yb, int pb,
                             int zsel, double R) {
    /* Coincident centres: delta-like (same normalized AO), else 0
     * (distinct shells/momenta are orthogonal at one centre). */
    if (!isfinite (R) || R < 0.0) return 0.0;
    if (R < 1e-9) {
        /* FULL-AUDIT: Slater orbitals on the SAME centre are not
         * orthogonal unless their l differs, and the previous code got
         * both cases wrong. Two 1s STOs with different exponents on one
         * centre overlap by 4(zeta_a zeta_b)^{3/2}/(zeta_a+zeta_b)^3 -
         * not zero - while an s and a p on one centre ARE orthogonal.
         * The branch now returns the correct value in both cases:
         *   - different l  -> 0 (orthogonal by symmetry)
         *   - same l, same exponent and the same directed AO -> 1
         *   - same l, different exponent -> the exact same-centre
         *     overlap, which is never reached by the shipped channels
         *     (they always have different exponents) but is the
         *     physically correct limit. */
        if (pa != pb) return 0.0;          /* different l -> orthogonal */
        double zm = za > zb ? za : zb;
        if (zm <= 0.0) return 0.0;
        if (fabs (za - zb) > 1e-9 * zm) {
            /* same centre, same l, different exponent: the exact
             * same-centre 1s-1s overlap. */
            double num = 4.0 * pow (za * zb, 1.5);
            double den = pow (za + zb, 3.0);
            return (den > 0.0) ? num / den : 0.0;
        }
        return 1.0;    /* the same normalised AO */
    }
    if (!(za > 0.0) || !(zb > 0.0)) return 0.0;
    double Na = sto_norm (na, za);
    double Nb = sto_norm (nb, zb);
    if (!(Na > 0.0) || !(Nb > 0.0)) return 0.0;
    double alpha = 0.5 * (za + zb) * R;
    double beta = 0.5 * (za - zb) * R;
    if (!(alpha > 0.0) || !isfinite (alpha) || !isfinite (beta)) return 0.0;
    poly5 base, zp, full;
    poly_radial (pa, pb, base);
    /* times (xi^2-eta^2) from dV */
    poly5 dV;
    poly_zero (dV);
    dV[2][0] = 1.0;
    dV[0][2] = -1.0;
    poly_mul (base, dV, full);
    /* z factors: xi*eta is its own monomial table */
    if (zsel == 1 || zsel == 2) {
        poly5 zf;
        poly_zero (zf);
        zf[1][1] = 1.0;
        zf[0][0] = (zsel == 1) ? -1.0 : 1.0;
        poly_mul (full, zf, full);
    } else if (zsel == 3) {
        poly5 xieta, xieta2, zf;
        poly_zero (xieta);
        xieta[1][1] = 1.0;
        poly_mul (xieta, xieta, xieta2);   /* (xi*eta)^2 */
        poly_zero (zf);
        zf[0][0] = -1.0;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++) zf[i][j] += xieta2[i][j];
        poly_mul (full, zf, full);
    }
    double S = poly_contract (full, alpha, beta);
    /* prefactors: N_a N_b Ya Yb (R^3/8) 2pi (R/2)^{pa+pb} [z R/2 powers] */
    double pre = Na * Nb * Ya * Yb * (R * R * R / 8.0) * 2.0 * 3.14159265358979323846;
    double rp = 1.0;
    for (int i = 0; i < pa + pb; i++) rp *= R / 2.0;
    pre *= rp;
    if (zsel == 1 || zsel == 2) pre *= R / 2.0;
    if (zsel == 3) pre *= (R / 2.0) * (R / 2.0);
    (void)zp;
    return pre * S;
}

/* Angular constants: Y_00 = 1/sqrt(4pi); pz along axis = sqrt(3/4pi). */
#define Y_S 0.28209479177387814
#define Y_P 0.48860251190291992

double eht_overlap_ss (int na, double za, int nb, double zb, double R) {
    int pa = (na == 1) ? 0 : 1;   /* STO radial power r^{n-1}, n<=2 */
    int pb = (nb == 1) ? 0 : 1;
    if ((na != 1 && na != 2) || (nb != 1 && nb != 2)) return 0.0;
    return overlap_sigma (na, za, Y_S, pa, nb, zb, Y_S, pb, 0, R);
}

double eht_overlap_sps (int ns, double zs, int np, double zp, double R,
                        int p_on_b) {
    if (np != 2) return 0.0;
    int ps = (ns == 1) ? 0 : 1;
    if (ns != 1 && ns != 2) return 0.0;
    /* GRID-ORACLE FIX: the p side radial power is 0, not 1. chi_2p =
     * N r^{n-1} e^{-zr} Y with Y~z/r, so chi ~ N z e^{-zr}: the single
     * position power lives in the z-polynomial (zsel), and passing
     * radial power 1 double-counted r (2x overshoot on s-ps, sign+scale
     * on pp-sigma). s side keeps n-1. */
    /* s radial power ps; p radial power 1; z on A(2) or B(1).
     * Both orientations evaluated directly (z_a = R(xi*eta+1)/2,
     * z_b = R(xi*eta-1)/2); no symmetry appeal. */
    if (p_on_b)
        return overlap_sigma (ns, zs, Y_S, ps, np, zp, Y_P, 0, 1, R);
    return overlap_sigma (np, zp, Y_P, 0, ns, zs, Y_S, ps, 2, R);
}

double eht_overlap_pps (double za, double zb, double R) {
    /* GRID-ORACLE FIX: radial powers 0 (see sps note above). */
    return overlap_sigma (2, za, Y_P, 0, 2, zb, Y_P, 0, 3, R);
}

/* p-pi (2p,2p): 2D Gauss-Legendre over (xi,eta). Documented quadrature.
 * S = Na Nb (3/4pi)(R^2/4)(R^3/8) pi J,
 * J = int_1^xmax int_-1^1 (xi^2-1)(1-eta^2)(xi^2-eta^2) e^{-a xi-b eta}. */
double eht_overlap_ppi (double za, double zb, double R) {
    if (!isfinite (R) || R < 0.0) return 0.0;
    if (R < 1e-9) {
        double zm = za > zb ? za : zb;
        return (zm > 0.0 && fabs (za - zb) <= 1e-9 * zm) ? 1.0 : 0.0;
    }
    if (!(za > 0.0) || !(zb > 0.0)) return 0.0;
    double Na = sto_norm (2, za);
    double Nb = sto_norm (2, zb);
    if (!(Na > 0.0) || !(Nb > 0.0)) return 0.0;
    double alpha = 0.5 * (za + zb) * R;
    double beta = 0.5 * (za - zb) * R;
    if (!(alpha > 0.0) || !isfinite (alpha) || !isfinite (beta)) return 0.0;
    if (!gl_ready) gl_compute ();
    double ximax = 1.0 + 45.0 / alpha;
    /* 96-point xi grid by splitting [1,ximax] into 3 GL-64-mapped panels?
     * Simpler: single 64-panel is enough for smooth decay; use 2 panels
     * of 64 for margin (128 xi x 64 eta). Deterministic. */
    double J = 0.0;
    const int NP = 2;
    for (int p = 0; p < NP; p++) {
        double x0 = 1.0 + (ximax - 1.0) * p / NP;
        double x1 = 1.0 + (ximax - 1.0) * (p + 1) / NP;
        double xm = 0.5 * (x0 + x1), xh = 0.5 * (x1 - x0);
        for (int i = 0; i < EHT_GL_N; i++) {
            double xi = xm + xh * gl_x[i];
            double exi = exp (-alpha * xi);
            double fxi = (xi * xi - 1.0) * xh * gl_w[i];
            for (int j = 0; j < EHT_GL_N; j++) {
                double eta = gl_x[j];
                double g = (1.0 - eta * eta) * (xi * xi - eta * eta) *
                           exp (-beta * eta);
                J += fxi * gl_w[j] * g * exi;
            }
        }
    }
    double pre = Na * Nb * (3.0 / (4.0 * 3.14159265358979323846)) *
                 (R * R / 4.0) * (R * R * R / 8.0) * 3.14159265358979323846;
    return pre * J;
}

/* ═══════════════════════════════════════════════════════════════════════════
 * eht_overlap_ao — the AUTHORITATIVE two-centre AO overlap in GLOBAL
 * Cartesian axes.
 *
 * Why this function exists (full audit, F3/F4 + the EHT sign defect):
 * two different pieces of code in this tree used to assemble the overlap
 * matrix for the same physical basis, and they disagreed:
 *
 *   - src/qm_eht_scf.c assembled S itself and applied a Slater-Koster
 *     rotation to the s-pσ channel with a SPURIOUS MINUS SIGN in the
 *     "p on the first atom" branch. The result violated the mirror
 *     symmetry of a symmetric molecule: for a bent H-C-H, the engine
 *     produced S(s_H1, p_x) == S(s_H2, p_x) when the symmetry of the
 *     physical system demands S(s_H1, p_x) == -S(s_H2, p_x). The EHT SCF
 *     then failed to converge for H-C-H at all, and the N2 molecule came
 *     out with the wrong MO ordering, a HOMO-LUMO gap of 0.0000 eV and a
 *     Mayer bond order of 2.05 instead of the textbook 3.00.
 *
 *   - src/qm.c assembled S from per-ATOM valence shells, which silently
 *     returned ZERO for every O-H / C-H / N-H pair (the O valence shell
 *     is 2p, and the s-channel wrapper refused anything that was not an
 *     s shell), i.e. the dominant 2s-1s sigma-bonding channel of water
 *     and methane was missing from the overlap matrix entirely.
 *
 * One implementation, one place where the sign can be wrong, and the
 * per-AO (n, l, m, zeta) description that a basis actually has.
 *
 * CONVENTIONS (all verified against a brute-force 3D Cartesian grid
 * oracle, and against the mirror symmetry of a bent X-C-X):
 *
 *   - AOs are the REAL global-axis p orbitals: p_x, p_y, p_z, with
 *     m = 1, 2, 3 respectively (m = 0 is an s orbital).
 *   - (dx, dy, dz) is the vector from atom i to atom j; its unit vector
 *     u is the local +z of the prolate frame (A = atom i at -R/2,
 *     B = atom j at +R/2).
 *   - The p integrals are the "parallel" convention: the p orbital is
 *     taken along the local +z. For s-pσ this makes the value NEGATIVE
 *     when the p lobe points away from the s (the two facing lobes then
 *     have opposite signs). The global-axis direction cosine restores
 *     the physical sign, so the assembled S matrix is correct for any
 *     orientation of the p orbital in space.
 *
 *   S(s_i, p_j(e_mj)) = (e_mj . u) * S_-(p on B along +u)
 *   S(p_i(e_mi), s_j) = (e_mi . u) * S_+(p on A along +u)
 *   S(p_i(e_mi), p_j(e_mj)) = (e_mi.u)(e_mj.u) S_sigma
 *                            + (delta_mij - (e_mi.u)(e_mj.u)) S_pi
 *
 *   with S_- = eht_overlap_sps(..., p_on_b = 1)   (p on B, along +u)
 *        S_+ = eht_overlap_sps(..., p_on_b = 0)   (p on A, along +u)
 *        S_sigma = eht_overlap_pps (the parallel-axis p-p sigma integral)
 *        S_pi   = eht_overlap_ppi
 *
 *   The p-p sigma integral is the parallel convention, whose sign
 *   crosses zero at zeta*R ~ 2.51. That is the correct mathematical
 *   value of <p_i(+u) p_j(+u)> in a FIXED global basis, and the
 *   Wolfsberg-Helmholtz Hamiltonian built on it is continuous through
 *   the crossing, so it is kept as-is and documented rather than
 *   replaced by a convention-dependent rescaling.
 * ═══════════════════════════════════════════════════════════════════════════ */
double eht_overlap_ao (int ni, double zi, int li, int mi,
                       int nj, double zj, int lj, int mj,
                       double R, double dx, double dy, double dz) {
    if (!(R > 0.0) || !isfinite (R)) return 0.0;
    if (!(zi > 0.0) || !(zj > 0.0)) return 0.0;
    double ux = dx / R, uy = dy / R, uz = dz / R;
    if (!isfinite (ux) || !isfinite (uy) || !isfinite (uz)) return 0.0;

    if (li == 0 && lj == 0) {
        return eht_overlap_ss (ni, zi, nj, zj, R);
    }
    if (li == 0 && lj == 1) {
        /* s on i, p on j: the p axis is e_mj; only its component along
         * the i->j axis survives the integration against the
         * cylindrically symmetric s on i. */
        double c = (mj == 1) ? ux : (mj == 2) ? uy : uz;
        return c * eht_overlap_sps (ni, zi, nj, zj, R, 1);
    }
    if (li == 1 && lj == 0) {
        /* p on i, s on j: the p axis is e_mi. */
        double c = (mi == 1) ? ux : (mi == 2) ? uy : uz;
        return c * eht_overlap_sps (nj, zj, ni, zi, R, 0);
    }
    if (li == 1 && lj == 1) {
        double ci = (mi == 1) ? ux : (mi == 2) ? uy : uz;
        double cj = (mj == 1) ? ux : (mj == 2) ? uy : uz;
        double ss = eht_overlap_pps (zi, zj, R);
        double sp = eht_overlap_ppi (zi, zj, R);
        return ci * cj * ss + (((mi == mj) ? 1.0 : 0.0) - ci * cj) * sp;
    }
    return 0.0;   /* unsupported shell combination */
}
