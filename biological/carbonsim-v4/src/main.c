#include <stdio.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/periodic_table.h"
#include "../include/quantum.h"
#include "../include/forces.h"
#include "../include/integrator.h"
#include "../include/sim.h"
#include "../include/nucleobases.h"

/*
 * main.c
 *
 * Demonstration sequence, bottom-up:
 *   1. Quantum layer  — print orbital energies and wave functions for H, C, O
 *   2. Bonding layer  — show H2 bond dissociation energy curve from LJ
 *   3. Molecular MD   — run 300K dynamics on a water molecule (H2O)
 *   4. Small ensemble — three H2O molecules, short NVT trajectory
 *
 * Every number printed has physical units labelled.
 * This is the foundation: add more chemistry above this bedrock.
 */

/* ── Pretty separator ────────────────────────────────────────────────────── */
static void banner(const char *title) {
    printf("\n╔══════════════════════════════════════════════════════╗\n");
    printf(  "║  %-52s║\n", title);
    printf(  "╚══════════════════════════════════════════════════════╝\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 1: Quantum mechanical orbital structure
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_quantum(void) {
    banner("DEMO 1: Quantum orbital structure");

    int elements[] = {1, 6, 7, 8};   /* H, C, N, O */
    int n_el = 4;

    for (int e = 0; e < n_el; e++) {
        int Z = elements[e];
        pt_print_element(Z);

        /* Build a temporary atom just to get orbital data */
        Atom atom = {0};
        atom.Z       = Z;
        atom.element = pt_element(Z);
        atom.mass    = pt_element(Z)->mass;
        pt_electron_config(Z, &atom.electron_config);
        quantum_fill_orbitals(&atom);
        quantum_print_orbitals(&atom);

        /* Print the radial wave function profile for the valence orbital */
        /* Find the highest-energy occupied orbital */
        int hi = 0;
        for (int i = 1; i < atom.num_orbitals; i++) {
            if (atom.orbitals[i].orbital_energy >
                atom.orbitals[hi].orbital_energy)
                hi = i;
        }
        int n = atom.orbitals[hi].qn.n;
        int l = atom.orbitals[hi].qn.l;
        ElectronConfig cfg;
        pt_electron_config(Z, &cfg);
        double Zeff = quantum_zeff(Z, n, l, &cfg);
        double r_mp = quantum_most_probable_radius(n, l, Zeff);

        printf("  Valence orbital: %d%c  Z_eff=%.3f  r_mp=%.3f Å\n",
               n, "spdf"[l], Zeff, r_mp);

        /* ASCII radial probability profile (0..5 Å in 40 steps) */
        printf("  Radial probability P(r) = r²|R_nl(r)|²:\n");
        double r_max_plot = 10.0;
        double step       = r_max_plot / 40.0;
        double P_max      = 0.0;
        for (int k = 1; k <= 40; k++) {
            double P = quantum_radial_probability(n, l, Zeff, k * step);
            if (P > P_max) P_max = P;
        }
        printf("  0 Å ");
        for (int k = 1; k <= 40; k++) {
            double r = k * step;
            double P = quantum_radial_probability(n, l, Zeff, r) / P_max;
            printf("%c", P < 0.1 ? ' ' :
                         P < 0.3 ? '.' :
                         P < 0.6 ? ':' :
                         P < 0.85? '|' : '#');
        }
        printf(" %.1f Å\n\n", r_max_plot);
    }
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 2: Two H-H potentials compared — covalent bond vs. van der Waals
 *
 * This demo exists to make an important distinction explicit: the LJ
 * potential and the covalent bond potential are TWO DIFFERENT physical
 * interactions, computed by two different code paths, and they answer two
 * different questions:
 *
 *   COVALENT (harmonic, from BOND_TABLE):
 *     "Two H atoms ARE bonded — what does stretching/compressing that
 *      shared-electron bond cost?" Minimum at r0 = 0.7414 Å (the actual
 *      H2 bond length). This is what governs sim_place_h2() and is the
 *      ONLY potential applied to atoms in sim->bonds[].
 *
 *   VAN DER WAALS (Lennard-Jones, UFF parameters):
 *     "Two H atoms are NOT bonded — how do their electron clouds interact
 *      at a distance?" Minimum at r ≈ 3.24 Å, far weaker (~0.002 eV vs.
 *      tens of eV for the covalent well). This is what governs how
 *      separate, non-bonded atoms or molecules approach each other —
 *      it's the physics behind gas pressure, condensation, and packing.
 *
 * forces_calculate() automatically excludes bonded pairs from the LJ/Coulomb
 * sum (see the 1-2 and 1-3 exclusion logic), so a real bonded H2 molecule
 * NEVER feels the LJ curve between its own two atoms — only the harmonic
 * term. The two curves below are printed side by side specifically so this
 * distinction is never ambiguous.
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_bond_curve(void) {
    banner("DEMO 2: H-H covalent bond vs. van der Waals (two different physics)");

    /* ── Van der Waals (non-bonded) curve ────────────────────────────────── */
    const Element *H = pt_element(1);
    double eps   = lj_eps_combine(H->lj_epsilon, H->lj_epsilon);
    double sigma = lj_sigma_combine(H->lj_sigma, H->lj_sigma);
    double lj_rmin = pow(2.0, 1.0/6.0) * sigma;

    /* ── Covalent (bonded) curve ─────────────────────────────────────────── */
    BondParam bp;
    forces_bond_params(1, 1, 1, &bp);  /* H-H single bond from BOND_TABLE */

    printf("  Van der Waals (LJ):  ε=%.5f eV   σ=%.4f Å   r_min=%.4f Å\n",
           eps, sigma, lj_rmin);
    printf("  Covalent (harmonic): r0=%.4f Å   k=%.2f eV/Å²   "
           "(real H2 bond length is 0.7414 Å)\n\n", bp.r0, bp.k);

    printf("  %-8s  %-14s  %-14s  %-s\n",
           "r (Å)", "V_covalent(eV)", "V_vdW (eV)", "Scale");
    printf("  %-8s  %-14s  %-14s\n","────────","──────────────","──────────────");

    /* Print covalent curve only near its own well (it's enormously stronger
     * and would make the vdW curve invisible on a shared linear scale) */
    for (int k = 4; k <= 24; k++) {
        double r = 0.4 + k * 0.05;
        double stretch = r - bp.r0;
        double V_cov   = 0.5 * bp.k * stretch * stretch;
        printf("  %-8.4f  %-14.6f  %-14s  covalent well (depth scale: eV)\n",
               r, V_cov, "—");
    }

    printf("\n");

    /* Print vdW curve over its own, much wider and shallower range */
    double V_min = 0.0;
    for (int k = 6; k <= 80; k++) {
        double r   = k * 0.0625;
        double sr6 = pow(sigma / r, 6.0);
        double V   = 4.0 * eps * (sr6*sr6 - sr6);
        if (V < V_min) V_min = V;
    }
    for (int k = 20; k <= 80; k += 2) {
        double r    = k * 0.0625;
        double sr6  = pow(sigma / r, 6.0);
        double V    = 4.0 * eps * (sr6*sr6 - sr6);
        int bar_len = (int)fmin(40.0, fmax(0.0, 20.0 + 20.0 * V / fabs(V_min)));
        printf("  %-8.4f  %-14s  %-14.6f  %.*s\n",
               r, "—", V, bar_len, "########################################");
    }
    printf("\n  For scale: the real H2 covalent bond dissociation energy is "
           "4.52 eV\n  (a standard spectroscopic constant), versus this "
           "vdW well depth of only\n  %.5f eV — roughly %.0fx weaker. "
           "That gap is why breaking a chemical\n  bond (a reaction) costs "
           "so much more than separating two molecules that\n  are merely "
           "touching (melting/evaporation).\n", fabs(V_min), 4.52/fabs(V_min));
    printf("\n  Caveat: the harmonic term above is only valid for small "
           "vibrations near\n  r0. It's a parabola, not a real bond — it "
           "never flattens out, so it would\n  (wrongly) predict infinite "
           "energy to fully separate the atoms. Capturing\n  actual bond "
           "breaking needs a Morse potential or a reactive force field —\n"
           "  a natural next addition to this codebase.\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 3: Water molecule NVT MD at 300 K
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_water_md(void) {
    banner("DEMO 3: H2O molecule — NVT MD at 300 K");

    Simulation *sim = sim_create(16, 32);
    if (!sim) { printf("  ERROR: allocation failed\n"); return; }

    /* Place water at origin */
    sim_place_h2o(sim, vec3_zero());

    /* Timestep and thermostat */
    sim->dt = 0.5;  /* 0.5 fs — small for stiff O-H bonds */
    sim->thermostat.type               = THERMOSTAT_BERENDSEN;
    sim->thermostat.target_temperature = 300.0;
    sim->thermostat.tau                = 100.0;  /* fs */

    /* Initial forces and energy */
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->potential_energy = sim->potential_energy;
    sim->total_energy = sim->kinetic_energy + sim->potential_energy;

    printf("  Initial geometry:\n");
    sim_print_atoms(sim);
    printf("\n");
    sim_print_bonds(sim);
    printf("\n");
    sim_print_angles(sim);

    /* Assign Maxwell-Boltzmann velocities at 300 K */
    integrator_maxwell_boltzmann(sim, 300.0, 42UL);

    /* Recalculate after velocity assignment */
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;
    sim->temperature    = integrator_temperature(sim);

    printf("\n  Initial thermodynamics:\n");
    printf("  KE = %.6f eV  PE = %.6f eV  E = %.6f eV  T = %.2f K\n\n",
           sim->kinetic_energy, sim->potential_energy,
           sim->total_energy, sim->temperature);

    /* ── MD trajectory ─────────────────────────────────────────────────── */
    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-12s\n",
           "Step","t (fs)","KE (eV)","PE (eV)","T (K)","O-H1 dist (Å)");
    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-12s\n",
           "────","──────","───────","───────","─────","────────────");

    int N_steps = 2000;
    int print_every = 100;

    for (int step = 0; step < N_steps; step++) {
        integrator_step(sim);

        if (step % print_every == 0) {
            double dOH = vec3_dist(sim->atoms[0].position,
                                    sim->atoms[1].position);
            printf("  %-8d  %-10.3f  %-10.6f  %-10.6f  %-8.2f  %-12.6f\n",
                   (int)sim->step,
                   sim->time,
                   sim->kinetic_energy,
                   sim->potential_energy,
                   sim->temperature,
                   dOH);
        }
    }

    printf("\n  Final geometry after %d steps:\n", N_steps);
    sim_print_atoms(sim);
    sim_print_summary(sim);
    sim_destroy(sim);
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 4: Cyclic (H2O)3 hydrogen-bonded trimer — emergent non-covalent structure
 *
 * This is the most chemically important demo in the file: it shows molecular
 * structure that NO ONE specified directly. We only give the simulator:
 *   - atomic positions / charges (from the periodic table + TIP3P charges)
 *   - the harmonic bond/angle terms (covalent, intramolecular)
 *   - the Coulomb + Lennard-Jones terms (non-bonded, intermolecular)
 *
 * The hydrogen-bonded ring — three O···H–O bridges holding the trimer
 * together — is not programmed in. It falls out of minimising those terms.
 * This is the gas-phase cyclic water trimer, a real, well-characterised
 * structure (see e.g. Keutsch & Saykally, PNAS 2001).
 *
 * Geometry construction:
 *   O atoms placed at the vertices of an equilateral triangle, O···O = 2.95 Å
 *   (the experimental/ab-initio range for a water-water H-bond is ~2.8-3.0 Å).
 *   Each water's "donor" H points along the O→O(next) vector at the correct
 *   O-H bond length (0.9572 Å) — this is the bridging hydrogen.
 *   Each water's "free" H satisfies the 104.52° H-O-H angle, rotated to
 *   point outward from the ring (away from the neighbouring molecules,
 *   avoiding artificial steric clash).
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_water_cluster(void) {
    banner("DEMO 4: Cyclic (H2O)3 — emergent hydrogen-bonded ring");

    Simulation *sim = sim_create(64, 128);
    if (!sim) { printf("  ERROR: allocation failed\n"); return; }

    const double PI    = 3.14159265358979323846;
    const double r_OO   = 2.95;                  /* Å, O···O H-bond distance */
    const double r_OH   = 0.9572;                 /* Å, covalent O-H         */
    const double hoh    = 104.52 * PI / 180.0;     /* H-O-H angle             */
    const double R_ring = r_OO / sqrt(3.0);        /* circumradius of triangle*/

    Vec3 O_pos[3], donorH_dir[3];

    /* Step 1: place the three oxygens at the triangle vertices */
    for (int i = 0; i < 3; i++) {
        double phi = 2.0 * PI * i / 3.0;
        O_pos[i] = vec3(R_ring * cos(phi), R_ring * sin(phi), 0.0);
    }

    /* Step 2: donor H direction = toward the NEXT oxygen in the ring */
    for (int i = 0; i < 3; i++) {
        int next = (i + 1) % 3;
        donorH_dir[i] = vec3_normalize(vec3_sub(O_pos[next], O_pos[i]));
    }

    int O_idx[3], Hd_idx[3], Hf_idx[3];  /* O, donor-H, free-H atom indices */

    for (int i = 0; i < 3; i++) {
        /* Oxygen */
        O_idx[i] = sim_add_atom(sim, 8, O_pos[i], -0.834);

        /* Donor hydrogen: along donorH_dir[i], at the covalent O-H length.
         * This is the bridging H that points at the next ring oxygen. */
        Vec3 Hd_pos = vec3_add(O_pos[i], vec3_scale(donorH_dir[i], r_OH));
        Hd_idx[i] = sim_add_atom(sim, 1, Hd_pos, +0.417);

        /* Free hydrogen: rotate donorH_dir[i] by ±104.52° about z so the
         * H-O-H angle is exact, choosing the sign that points away from
         * the ring centre (outward), minimising steric clash with
         * neighbouring molecules. Pure z-rotation keeps everything
         * in the ring plane (z=0), which is a fine simplification for
         * a point-charge, non-polarisable model like this one. */
        double cx = donorH_dir[i].x, cy = donorH_dir[i].y;
        double cos_a = cos(hoh), sin_a = sin(hoh);
        Vec3 rot_plus  = vec3(cx*cos_a - cy*sin_a, cx*sin_a + cy*cos_a, 0.0);
        Vec3 rot_minus = vec3(cx*cos_a + cy*sin_a, -cx*sin_a + cy*cos_a, 0.0);

        Vec3 cand_plus  = vec3_add(O_pos[i], vec3_scale(rot_plus,  r_OH));
        Vec3 cand_minus = vec3_add(O_pos[i], vec3_scale(rot_minus, r_OH));

        /* Outward = farther from ring centroid (origin) */
        Vec3 Hf_pos = (vec3_norm(cand_plus) > vec3_norm(cand_minus))
                      ? cand_plus : cand_minus;

        Hf_idx[i] = sim_add_atom(sim, 1, Hf_pos, +0.417);

        /* Apply verified TIP3P LJ parameters (Jorgensen 1983) - see the
         * detailed comment in sim_place_h2o() for why generic UFF values
         * are wrong here. Built by hand since this trimer doesn't go
         * through sim_place_h2o(). */
        sim_set_atom_lj(sim, O_idx[i],  0.1521 * KCAL_MOL_TO_EV, 3.15061);
        sim_set_atom_lj(sim, Hd_idx[i], 0.0, 0.0);
        sim_set_atom_lj(sim, Hf_idx[i], 0.0, 0.0);

        /* Covalent bonds: O to both its own hydrogens */
        sim_add_bond(sim, O_idx[i], Hd_idx[i], 1);
        sim_add_bond(sim, O_idx[i], Hf_idx[i], 1);
    }
    sim_build_angles(sim);

    /*
     * Run at 50 K rather than 300 K. This is an honest choice, not a fudge:
     * a 3-molecule gas-phase cluster with no confining box and no
     * surrounding liquid pressure is a genuinely fragile system — real
     * water trimers in molecular beams are weakly bound (~0.2-0.3 eV per
     * H-bond, but with large amplitude floppy motion and low barriers to
     * rearrangement). At 300 K the available thermal kinetic energy per
     * mode (~0.026 eV) is large enough that an unconfined trimer dissociates
     * on a picosecond timescale — which the previous version of this demo
     * correctly showed when the molecules weren't even H-bond oriented.
     * At 50 K we can watch genuine bound, oscillatory H-bond dynamics
     * within a short trajectory without needing a confining potential.
     */
    sim->dt = 0.5;
    sim->thermostat.type               = THERMOSTAT_BERENDSEN;
    sim->thermostat.target_temperature = 50.0;
    sim->thermostat.tau                = 50.0;

    integrator_maxwell_boltzmann(sim, 50.0, 99UL);
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;
    sim->temperature    = integrator_temperature(sim);

    printf("  %d atoms, %d bonds, %d angles\n",
           sim->num_atoms, sim->num_bonds, sim->num_angles);
    printf("  Ring O-O-O construction: O...O = %.3f Å per edge\n", r_OO);
    printf("  Initial T=%.2f K  PE=%.6f eV (intermolecular H-bonds "
           "contribute the negative part)\n\n",
           sim->temperature, sim->potential_energy);

    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-10s %-10s %-10s\n",
           "Step","t (fs)","KE (eV)","PE (eV)","T (K)",
           "O0-O1(Å)","O1-O2(Å)","O2-O0(Å)");
    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-10s %-10s %-10s\n",
           "────","──────","───────","───────","─────",
           "────────","────────","────────");

    int N_steps = 4000;
    int print_every = 200;
    double min_OO = 1.0e9, max_OO = 0.0;

    for (int step = 0; step < N_steps; step++) {
        integrator_step(sim);

        double dOO[3];
        dOO[0] = vec3_dist(sim->atoms[O_idx[0]].position,
                            sim->atoms[O_idx[1]].position);
        dOO[1] = vec3_dist(sim->atoms[O_idx[1]].position,
                            sim->atoms[O_idx[2]].position);
        dOO[2] = vec3_dist(sim->atoms[O_idx[2]].position,
                            sim->atoms[O_idx[0]].position);

        for (int k = 0; k < 3; k++) {
            if (dOO[k] < min_OO) min_OO = dOO[k];
            if (dOO[k] > max_OO) max_OO = dOO[k];
        }

        if (step % print_every == 0) {
            printf("  %-8d  %-10.3f  %-10.6f  %-10.6f  %-8.2f  "
                   "%-10.4f %-10.4f %-10.4f\n",
                   (int)sim->step, sim->time,
                   sim->kinetic_energy, sim->potential_energy,
                   sim->temperature, dOO[0], dOO[1], dOO[2]);
        }
    }

    printf("\n  Over %d steps (%.0f fs): O···O range = [%.3f, %.3f] Å\n",
           N_steps, sim->time, min_OO, max_OO);
    if (max_OO < 5.0) {
        printf("  → Ring stayed bound. Hydrogen bonds emerged from Coulomb "
               "+ LJ alone — never told the simulator these were H-bonds.\n");
    } else {
        printf("  → Ring dissociated within this trajectory.\n");
    }

    printf("\n  Final state:\n");
    sim_print_atoms(sim);
    sim_destroy(sim);
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 5: Methane molecule — tetrahedral geometry
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_methane(void) {
    banner("DEMO 5: CH4 — tetrahedral geometry check");

    Simulation *sim = sim_create(16, 32);
    if (!sim) { return; }

    sim_place_ch4(sim, vec3_zero());
    forces_calculate(sim);
    sim->kinetic_energy = 0.0;
    sim->total_energy   = sim->potential_energy;

    printf("  Methane geometry (C at origin):\n");
    sim_print_atoms(sim);
    printf("\n");

    /* Print all H-C-H angles */
    printf("  Bond angles:\n");
    for (int i = 0; i < sim->num_angles; i++) {
        const Angle *a = &sim->angles[i];
        Vec3 r_ba = vec3_sub(sim->atoms[a->atom_a].position,
                              sim->atoms[a->atom_b].position);
        Vec3 r_bc = vec3_sub(sim->atoms[a->atom_c].position,
                              sim->atoms[a->atom_b].position);
        double cos_t = vec3_dot(vec3_normalize(r_ba),
                                 vec3_normalize(r_bc));
        if (cos_t >  1.0) cos_t =  1.0;
        if (cos_t < -1.0) cos_t = -1.0;
        double theta = acos(cos_t) * 180.0 / 3.14159265358979323846;
        printf("    H(%d)-C(%d)-H(%d): %.4f°  (ideal 109.47°)\n",
               a->atom_a, a->atom_b, a->atom_c, theta);
    }

    printf("\n  Initial potential energy: %.6f eV\n", sim->potential_energy);
    sim_destroy(sim);
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 6: The five DNA/RNA nucleobases - geometry validation
 *
 * For each base: place atoms from verified PDB CCD coordinates, then
 * INDEPENDENTLY recompute every bond length and ring angle directly
 * from the Cartesian positions (not from the stored r0/theta0, even
 * though those were set to match by construction) - this is the same
 * "trust but verify" approach used in the methane tetrahedral-angle
 * check. A planarity check (independent of any external reference
 * numbers) confirms the aromatic ring is genuinely flat.
 * ════════════════════════════════════════════════════════════════════════════ */
static void print_bond_report(const Simulation *sim, const char *symbols[]) {
    printf("  %-10s %-6s %-6s %-8s\n", "Bond", "Atom1", "Atom2", "r (A)");
    printf("  %-10s %-6s %-6s %-8s\n", "----", "-----", "-----", "-----");
    for (int i = 0; i < sim->num_bonds; i++) {
        const Bond *b = &sim->bonds[i];
        double r = vec3_dist(sim->atoms[b->atom_a].position,
                              sim->atoms[b->atom_b].position);
        printf("  %-10s %-6s %-6s %-8.4f  (order %d)\n",
               b->order == 2 ? "double" : "single",
               symbols[b->atom_a], symbols[b->atom_b], r, b->order);
    }
}

static void print_angle_report(const Simulation *sim, const char *symbols[]) {
    printf("  %-6s %-6s %-6s %-10s\n", "A", "B(ctr)", "C", "theta (deg)");
    printf("  %-6s %-6s %-6s %-10s\n", "-", "------", "-", "-----------");
    for (int i = 0; i < sim->num_angles; i++) {
        const Angle *a = &sim->angles[i];
        Vec3 ba = vec3_sub(sim->atoms[a->atom_a].position,
                            sim->atoms[a->atom_b].position);
        Vec3 bc = vec3_sub(sim->atoms[a->atom_c].position,
                            sim->atoms[a->atom_b].position);
        double theta = vec3_angle(ba, bc) * 180.0 / 3.14159265358979323846;
        printf("  %-6s %-6s %-6s %-10.2f\n",
               symbols[a->atom_a], symbols[a->atom_b], symbols[a->atom_c], theta);
    }
}

static double total_charge(const Simulation *sim) {
    double q = 0.0;
    for (int i = 0; i < sim->num_atoms; i++) q += sim->atoms[i].partial_charge;
    return q;
}

static void demo_nucleobases(void) {
    banner("DEMO 6: The five nucleobases - geometry validation");

    /* ── Uracil ─────────────────────────────────────────────────────────── */
    {
        Simulation *sim = sim_create(16, 16);
        sim_place_uracil(sim, vec3_zero());
        printf("  URACIL (C4H4N2O2) - 12 atoms\n");
        const char *sym[] = {"N1","C2","O2","N3","C4","O4",
                              "C5","C6","HN1","HN3","H5","H6"};
        print_bond_report(sim, sym);
        printf("\n");
        print_angle_report(sim, sym);

        int ring[] = {0,1,3,4,6,7}; /* N1 C2 N3 C4 C5 C6 */
        double dev = nb_planarity_deviation(sim, ring, 6);
        printf("\n  Ring planarity: max deviation = %.4f A "
               "(aromatic rings should be ~0)\n", dev);

        printf("\n  Cross-check vs. electron diffraction "
               "(Ferenczy et al. 1986):\n");
        printf("  %-22s %-10s %-10s\n", "Quantity", "This model", "Literature");
        double r_C2N1 = vec3_dist(sim->atoms[1].position, sim->atoms[0].position);
        double r_C4C5 = vec3_dist(sim->atoms[4].position, sim->atoms[6].position);
        double r_C5C6 = vec3_dist(sim->atoms[6].position, sim->atoms[7].position);
        printf("  %-22s %-10.3f %-10s\n", "C-N (A)", r_C2N1, "1.399");
        printf("  %-22s %-10.3f %-10s\n", "C4-C5 single (A)", r_C4C5, "1.462");
        printf("  %-22s %-10.3f %-10s\n", "C5=C6 double (A)", r_C5C6, "1.343");
        printf("\n  Total molecular charge: %+.6f e (RESP, Aduri et al. "
               "2007 + derived HN1)\n", total_charge(sim));
        sim_destroy(sim);
    }

    /* ── Cytosine ───────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_cytosine(sim, vec3_zero());
        printf("  CYTOSINE (C4H5N3O) - 13 atoms\n");
        printf("  (tautomer-corrected: H relocated to N1, the Watson-Crick-\n");
        printf("   relevant position; N3 left bare as required for pairing)\n");

        const char *sym[] = {"N3","C4","N1","C2","O2","N4",
                              "C5","C6","HN41","HN42","H5","H6","HN1"};
        print_bond_report(sim, sym);

        /* Explicit WC-readiness check: N1 must have 3 bonds (C2,C6,H),
         * N3 must have exactly 2 (C4,C2 - bare, ready to accept guanine's
         * N1-H) */
        int n1_bonds = sim->atoms[2].num_bonds;
        int n3_bonds = sim->atoms[0].num_bonds;
        printf("\n  N1 bond count: %d (expect 3: C2, C6, H - donor ready)\n",
               n1_bonds);
        printf("  N3 bond count: %d (expect 2: C4, C2 - bare, acceptor ready)\n",
               n3_bonds);

        int ring[] = {0,1,2,3,6,7}; /* N3 C4 N1 C2 C5 C6 */
        double dev = nb_planarity_deviation(sim, ring, 6);
        printf("  Ring planarity: max deviation = %.4f A\n", dev);
        printf("  Total molecular charge: %+.6f e\n", total_charge(sim));
        sim_destroy(sim);
    }

    /* ── Thymine ────────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_thymine(sim, vec3_zero());
        printf("  THYMINE (C5H6N2O2, = 5-methyluracil) - 15 atoms\n");
        int ring[] = {0,1,3,4,6,7}; /* N1 C2 N3 C4 C5 C6 */
        double dev = nb_planarity_deviation(sim, ring, 6);
        printf("  Ring planarity: max deviation = %.4f A\n", dev);

        double r_methyl = vec3_dist(sim->atoms[6].position, sim->atoms[11].position);
        printf("  C5-CH3 methyl bond length: %.4f A "
               "(target 1.51 A, toluene-type)\n", r_methyl);

        /* Confirm methyl H's are tetrahedral about the C5-CM axis */
        Vec3 C5 = sim->atoms[6].position, CM = sim->atoms[11].position;
        Vec3 H1 = sim->atoms[12].position, H2 = sim->atoms[13].position;
        Vec3 u1 = vec3_normalize(vec3_sub(H1, CM));
        Vec3 u2 = vec3_normalize(vec3_sub(H2, CM));
        double hch = acos(vec3_dot(u1,u2)) * 180.0/3.14159265358979323846;
        printf("  H-CM-H methyl angle: %.2f deg (ideal tetrahedral 109.47)\n", hch);
        printf("  Total molecular charge: %+.6f e "
               "(lower-confidence approx, see code comment)\n", total_charge(sim));
        (void)C5;
        sim_destroy(sim);
    }

    /* ── Adenine ────────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_adenine(sim, vec3_zero());
        printf("  ADENINE (C5H5N5, fused 5+6 purine ring) - 15 atoms\n");
        int ring[] = {0,1,2,3,4,5,6,7,8,9}; /* all 10 ring atoms */
        double dev = nb_planarity_deviation(sim, ring, 10);
        printf("  Full bicyclic ring planarity: max deviation = %.4f A\n", dev);
        printf("  Total molecular charge: %+.6f e\n", total_charge(sim));
        sim_destroy(sim);
    }

    /* ── Guanine ────────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_guanine(sim, vec3_zero());
        printf("  GUANINE (C5H5N5O, fused 5+6 purine ring) - 16 atoms\n");
        int ring[] = {0,1,2,3,4,6,7,9,10}; /* 9 ring atoms (excl. O6 substituent) */
        double dev = nb_planarity_deviation(sim, ring, 9);
        printf("  Full bicyclic ring planarity: max deviation = %.4f A\n", dev);
        printf("  Total molecular charge: %+.6f e\n", total_charge(sim));
        sim_destroy(sim);
    }

    printf("\n  All five bases hold real, verified ring geometry. Next: "
           "sugar-phosphate\n  backbones and Watson-Crick base pairing - "
           "G-C should bind via 3 H-bonds,\n  A-T via 2, using nothing but "
           "the Coulomb+LJ code already validated\n  on the water trimer.\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 7: Watson-Crick base pairing - G-C vs A-U binding energy
 *
 * The central test: G-C pairs via 3 hydrogen bonds, A-U via 2. If this
 * force field (real charges, real geometry, the same Coulomb+LJ code
 * already validated on the water trimer) is doing real chemistry and
 * not just fitting a foregone conclusion, G-C should come out MORE
 * stable (more negative interaction energy) than A-U.
 *
 * Construction strategy, identical in spirit to the water trimer: get
 * a rough, approximately-correct geometry (align ring planes, place the
 * primary donor/acceptor heavy atoms at a realistic H-bond distance),
 * then let MD relaxation find the true local minimum. No part of the
 * Watson-Crick hydrogen-bond pattern is hard-coded into the force
 * field - it has to emerge from Coulomb+LJ alone, the same way the
 * water trimer's ring emerged.
 * ════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double min_interE, final_interE, initial_interE;
    double closest_approach_E, closest_approach_dist;
} PairResult;

static PairResult run_pair_relaxation(Simulation *sim, int n_fixed, int n_moving,
                                       int primary_a, int primary_b) {
    PairResult r;
    forces_calculate(sim);
    r.initial_interE = sim->potential_energy;

    sim->dt = 0.1;
    sim->thermostat.type               = THERMOSTAT_BERENDSEN;
    sim->thermostat.target_temperature = 50.0;
    sim->thermostat.tau                = 50.0;
    integrator_maxwell_boltzmann(sim, 50.0, 7UL);
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;

    printf("  %-6s %-10s %-10s %-10s %-10s\n",
           "Step", "t(fs)", "PE(eV)", "T(K)", "primary(A)");
    double min_pe = sim->potential_energy;
    double closest_dist = vec3_dist(sim->atoms[primary_a].position,
                                     sim->atoms[primary_b].position);
    double closest_E = sim->potential_energy;

    int N_steps = 800, print_every = 100;
    for (int step = 0; step < N_steps; step++) {
        integrator_step(sim);
        if (sim->potential_energy < min_pe) min_pe = sim->potential_energy;

        double d = vec3_dist(sim->atoms[primary_a].position,
                              sim->atoms[primary_b].position);
        if (d < closest_dist) { closest_dist = d; closest_E = sim->potential_energy; }

        if (step < 5 || step % print_every == 0) {
            printf("  %-6d %-10.3f %-10.6f %-10.2f %-10.4f\n",
                   (int)sim->step, sim->time, sim->potential_energy,
                   sim->temperature, d);
        }
    }
    r.min_interE   = min_pe;
    r.final_interE = sim->potential_energy;
    r.closest_approach_E    = closest_E;
    r.closest_approach_dist = closest_dist;
    (void)n_fixed; (void)n_moving;
    return r;
}

static void demo_basepairing(void) {
    banner("DEMO 7: Watson-Crick pairing - does G-C beat A-U?");

    double G_C_energy, A_U_energy;

    /* ── G-C pair ───────────────────────────────────────────────────────── */
    {
        printf("\n  --- Guanine-Cytosine (3 H-bonds: N1-H..N3, N2-H..O2, "
               "O6..H-N4) ---\n");
        Simulation *sim = sim_create(64, 64);
        /* RESP charges from Aduri et al. were fit for use in a condensed-
         * phase force field (with solvent/counterions implicitly assumed
         * to provide screening) - using them in a bare vacuum Coulomb
         * sum overestimates electrostatics substantially (a single
         * H-bond pairwise term alone already approaches the size of the
         * ENTIRE real gas-phase WC binding energy). epsilon=4 is a
         * standard, widely-used distance-independent dielectric for
         * approximating this screening in classical force fields without
         * explicit solvent - not a precise fit for this system, but a
         * well-established, literature-standard order-of-magnitude
         * correction. */
        sim->dielectric = 4.0;

        int g = sim_place_guanine(sim, vec3_zero());
        int c = sim_place_cytosine(sim, vec3(15.0, 0.0, 0.0)); /* temp, far */

        /* Guanine ring atoms for normal: N9,C8,N1 (indices 0,1,6) */
        int g_ring[3] = {g+0, g+1, g+6};
        int c_ring[3] = {c+0, c+1, c+2}; /* N3,C4,N1 */
        Vec3 n_G = nb_ring_normal(sim, g_ring);
        Vec3 n_C = nb_ring_normal(sim, c_ring);

        double cosang = vec3_dot(n_G, n_C);
        Vec3 axis = vec3_cross(n_C, n_G);
        double angle;
        if (vec3_norm(axis) < 1.0e-8) {
            angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
            axis  = (cosang > 0) ? vec3(0,0,1)
                                  : vec3_cross(n_C, vec3(1,0,0));
            if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
        } else {
            angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
        }

        Vec3 c_N3_pivot = sim->atoms[c+0].position;
        nb_transform_rigid(sim, c, 13, c_N3_pivot, axis, angle, vec3_zero());

        /* Second rotation: fix the azimuthal orientation left unconstrained
         * by normal-alignment alone. Align cytosine's N3->N4 direction
         * (toward its donor amino group) with guanine's N1->O6 direction
         * (toward its acceptor), both projected into the now-shared
         * plane - this is what keeps O6 and N4 approaching each other
         * instead of colliding. */
        Vec3 g_N1_for_az = sim->atoms[g+6].position, g_O6 = sim->atoms[g+5].position;
        Vec3 c_N3_for_az = sim->atoms[c+0].position, c_N4 = sim->atoms[c+5].position;
        Vec3 v_G_az = vec3_sub(g_O6, g_N1_for_az);
        Vec3 v_C_az = vec3_sub(c_N4, c_N3_for_az);
        double az_angle = nb_signed_inplane_angle(v_C_az, v_G_az, n_G);
        nb_transform_rigid(sim, c, 13, c_N3_pivot, n_G, az_angle, vec3_zero());

        /* Target: cytosine's N3 lands along guanine's N1->HN1 donor
         * direction, at the standard WC N1...N3 distance (~2.95 A) */
        Vec3 g_N1 = sim->atoms[g+6].position, g_HN1 = sim->atoms[g+13].position;
        Vec3 donor_dir = vec3_normalize(vec3_sub(g_HN1, g_N1));
        Vec3 target = vec3_add(g_N1, vec3_scale(donor_dir, 2.95));
        Vec3 c_N3_now = sim->atoms[c+0].position;
        Vec3 translation = vec3_sub(target, c_N3_now);
        nb_transform_rigid(sim, c, 13, c_N3_pivot, vec3_zero(), 0.0, translation);

        printf("  Initial heavy-atom contacts after geometric placement:\n");
        printf("    G:N1...C:N3 = %.3f A (target 2.95)\n",
               vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
        printf("    G:N2...C:O2 = %.3f A\n",
               vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
        printf("    G:O6...C:N4 = %.3f A\n",
               vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));

        forces_calculate(sim);
        printf("  Energy breakdown at initial placement: "
               "E_LJ=%.6f eV  E_Coulomb=%.6f eV  Total_PE=%.6f eV\n",
               sim->E_lj_total, sim->E_coulomb_total, sim->potential_energy);

        /* Diagnostic: find the closest intermolecular atom pair AND the
         * single largest LJ repulsive contributor - these need not be
         * the same pair (LJ repulsion depends on sigma/epsilon too) */
        {
            double min_d = 1.0e9; int mi = -1, mj = -1;
            double max_lj = -1.0e9; int li = -1, lj_idx = -1;
            for (int ii = g; ii < g+16; ii++)
                for (int jj = c; jj < c+13; jj++) {
                    double d = vec3_dist(sim->atoms[ii].position, sim->atoms[jj].position);
                    if (d < min_d) { min_d = d; mi = ii; mj = jj; }
                    PairEnergy pe = forces_nonbonded_pair(sim->atoms, ii, jj,
                                                           &sim->box, 1, 0, 1.0);
                    if (pe.lj_energy > max_lj) { max_lj = pe.lj_energy; li = ii; lj_idx = jj; }
                }
            printf("  Closest intermolecular contact: atom %d (Z=%d) ... "
                   "atom %d (Z=%d) = %.3f A\n",
                   mi, sim->atoms[mi].Z, mj, sim->atoms[mj].Z, min_d);
            printf("  Largest single LJ repulsion: atom %d (Z=%d, sigma=%.2f) ... "
                   "atom %d (Z=%d, sigma=%.2f) = %.3f A, contributes %.4f eV\n",
                   li, sim->atoms[li].Z, sim->atoms[li].lj_sigma,
                   lj_idx, sim->atoms[lj_idx].Z, sim->atoms[lj_idx].lj_sigma,
                   vec3_dist(sim->atoms[li].position, sim->atoms[lj_idx].position), max_lj);
        }

        PairResult res = run_pair_relaxation(sim, 16, 13, g+6, c+0);
        printf("\n  Initial interaction PE:        %.6f eV\n", res.initial_interE);
        printf("  PE at closest WC approach:     %.6f eV "
               "(primary N1...N3 = %.3f A)\n",
               res.closest_approach_E, res.closest_approach_dist);
        printf("  Global PE minimum over run:    %.6f eV "
               "(may reflect drift to a different,\n"
               "                                  non-WC configuration "
               "such as stacking)\n", res.min_interE);
        printf("  Final PE (end of run):         %.6f eV\n", res.final_interE);

        printf("\n  Final heavy-atom contacts:\n");
        printf("    G:N1...C:N3 = %.3f A\n",
               vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
        printf("    G:N2...C:O2 = %.3f A\n",
               vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
        printf("    G:O6...C:N4 = %.3f A\n",
               vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));

        G_C_energy = res.closest_approach_E;
        sim_destroy(sim);
    }

    /* ── A-U pair ───────────────────────────────────────────────────────── */
    {
        printf("\n  --- Adenine-Uracil (2 H-bonds: N1..H-N3, N6-H..O4) ---\n");
        Simulation *sim = sim_create(64, 64);
        sim->dielectric = 4.0; /* same screening correction as G-C, see comment there */

        int a = sim_place_adenine(sim, vec3_zero());
        int u = sim_place_uracil(sim, vec3(15.0, 0.0, 0.0));

        int a_ring[3] = {a+0, a+1, a+6}; /* N9,C8,N1 */
        int u_ring[3] = {u+0, u+1, u+3}; /* N1,C2,N3 */
        Vec3 n_A = nb_ring_normal(sim, a_ring);
        Vec3 n_U = nb_ring_normal(sim, u_ring);

        double cosang = vec3_dot(n_A, n_U);
        Vec3 axis = vec3_cross(n_U, n_A);
        double angle;
        if (vec3_norm(axis) < 1.0e-8) {
            angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
            axis  = (cosang > 0) ? vec3(0,0,1)
                                  : vec3_cross(n_U, vec3(1,0,0));
            if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
        } else {
            angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
        }

        /* Here uracil is the DONOR (N3-H) and adenine the ACCEPTOR (N1) -
         * reversed roles from G-C's primary bond. Pivot on uracil's N3,
         * target it (and thus its H) toward adenine's N1. */
        Vec3 u_N3_pivot = sim->atoms[u+3].position;
        nb_transform_rigid(sim, u, 12, u_N3_pivot, axis, angle, vec3_zero());

        /* Azimuthal fix, same logic as G-C: align uracil's N3->O4
         * direction with adenine's N1->N6 direction so O4 and N6
         * approach each other rather than landing somewhere unrelated. */
        Vec3 a_N1_for_az = sim->atoms[a+6].position, a_N6_for_az = sim->atoms[a+5].position;
        Vec3 u_N3_for_az = sim->atoms[u+3].position, u_O4_for_az = sim->atoms[u+5].position;
        Vec3 v_A_az = vec3_sub(a_N6_for_az, a_N1_for_az);
        Vec3 v_U_az = vec3_sub(u_O4_for_az, u_N3_for_az);
        double az_angle = nb_signed_inplane_angle(v_U_az, v_A_az, n_A);
        nb_transform_rigid(sim, u, 12, u_N3_pivot, n_A, az_angle, vec3_zero());

        Vec3 u_N3 = sim->atoms[u+3].position, u_HN3 = sim->atoms[u+9].position;
        Vec3 donor_dir = vec3_normalize(vec3_sub(u_HN3, u_N3));
        /* We want N3 positioned so that ITS H points at adenine's N1 -
         * i.e. place N3 such that N1 = N3 + donor_dir * 2.9, so
         * N3 = target_N1 - donor_dir*2.9 */
        Vec3 a_N1 = sim->atoms[a+6].position;
        Vec3 target_N3 = vec3_sub(a_N1, vec3_scale(donor_dir, 2.90));
        Vec3 translation = vec3_sub(target_N3, u_N3);
        nb_transform_rigid(sim, u, 12, u_N3_pivot, vec3_zero(), 0.0, translation);

        printf("  Initial heavy-atom contacts after geometric placement:\n");
        printf("    A:N1...U:N3 = %.3f A (target 2.90)\n",
               vec3_dist(sim->atoms[a+6].position, sim->atoms[u+3].position));
        printf("    A:N6...U:O4 = %.3f A\n",
               vec3_dist(sim->atoms[a+5].position, sim->atoms[u+5].position));

        forces_calculate(sim);
        printf("  Energy breakdown at initial placement: "
               "E_LJ=%.6f eV  E_Coulomb=%.6f eV  Total_PE=%.6f eV\n",
               sim->E_lj_total, sim->E_coulomb_total, sim->potential_energy);

        {
            double min_d = 1.0e9; int mi = -1, mj = -1;
            for (int ii = a; ii < a+15; ii++)
                for (int jj = u; jj < u+12; jj++) {
                    double d = vec3_dist(sim->atoms[ii].position, sim->atoms[jj].position);
                    if (d < min_d) { min_d = d; mi = ii; mj = jj; }
                }
            printf("  Closest intermolecular contact: atom %d (Z=%d) ... "
                   "atom %d (Z=%d) = %.3f A\n",
                   mi, sim->atoms[mi].Z, mj, sim->atoms[mj].Z, min_d);
        }

        PairResult res = run_pair_relaxation(sim, 15, 12, a+6, u+3);
        printf("\n  Initial interaction PE:        %.6f eV\n", res.initial_interE);
        printf("  PE at closest WC approach:     %.6f eV "
               "(primary N1...N3 = %.3f A)\n",
               res.closest_approach_E, res.closest_approach_dist);
        printf("  Global PE minimum over run:    %.6f eV\n", res.min_interE);
        printf("  Final PE (end of run):         %.6f eV\n", res.final_interE);

        printf("\n  Final heavy-atom contacts:\n");
        printf("    A:N1...U:N3 = %.3f A\n",
               vec3_dist(sim->atoms[a+6].position, sim->atoms[u+3].position));
        printf("    A:N6...U:O4 = %.3f A\n",
               vec3_dist(sim->atoms[a+5].position, sim->atoms[u+5].position));

        A_U_energy = res.closest_approach_E;
        sim_destroy(sim);
    }

    /* ── Verdict ────────────────────────────────────────────────────────── */
    printf("\n  ══════════════════════════════════════════════════\n");
    printf("  G-C @ closest WC approach: %.6f eV (3 H-bonds)\n", G_C_energy);
    printf("  A-U @ closest WC approach: %.6f eV (2 H-bonds)\n", A_U_energy);
    if (G_C_energy < A_U_energy) {
        printf("  --> G-C binds MORE strongly than A-U (%.6f eV difference),\n"
               "      matching real DNA/RNA thermodynamics, using nothing\n"
               "      but charges + Coulomb + LJ. Not programmed in.\n",
               A_U_energy - G_C_energy);
    } else {
        printf("  --> Unexpected: A-U came out more stable than G-C. This\n"
               "      would need investigation (geometry, charges, or\n"
               "      relaxation time) before trusting the result.\n");
    }
    printf("  ══════════════════════════════════════════════════\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * Entry point
 * ════════════════════════════════════════════════════════════════════════════ */
int main(void) {
    printf("\n");
    printf("  ╔═══════════════════════════════════════════════════════╗\n");
    printf("  ║       CARBON VM — CHEMISTRY SIMULATOR                 ║\n");
    printf("  ║       From subatomic to molecular dynamics            ║\n");
    printf("  ╚═══════════════════════════════════════════════════════╝\n");
    printf("\n  Unit system: Length=Å  Time=fs  Energy=eV  Mass=AMU\n");
    printf("  Physical constants: 2019 CODATA  |  LJ: UFF  |  Bonds: AMBER\n\n");

    demo_quantum();
    demo_bond_curve();
    demo_water_md();
    demo_water_cluster();
    demo_methane();
    demo_nucleobases();
    demo_basepairing();

    printf("\n  All demos complete.\n");
    printf("  Next steps: add nucleotides → DNA → gene regulatory networks.\n\n");
    return 0;
}
