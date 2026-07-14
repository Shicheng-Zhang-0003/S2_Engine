#include <stdio.h>
#include <string.h>
#include "../include/periodic_table.h"
#include "../include/constants.h"

/*
 * periodic_table.c
 *
 * Full element data for Z=1..18 (H→Ar), stubbed beyond.
 * Extend the table by filling in the remaining entries.
 *
 * Madelung orbital filling order (used by pt_electron_config):
 * 1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p
 *
 * LJ ε values converted from UFF kcal/mol → eV (÷ 23.0605):
 * LJ σ values in Å (UFF x1 column, the collision diameter).
 */

/* ── Madelung sequence: pairs of (n, l) in filling order ─────────────────── */
static const int MADELUNG_N[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
static const int MADELUNG_L[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
static const int MADELUNG_LEN = 19;

/* ── Orbital capacity: 2*(2l+1) ──────────────────────────────────────────── */
static int orbital_cap(int l) { return 2*(2*l+1); }

/* ══════════════════════════════════════════════════════════════════════════
 * Periodic table data
 * ══════════════════════════════════════════════════════════════════════════
 * Fields: Z, symbol, name, mass[AMU], electronegativity,
 *         atomic_radius[Å], covalent_radius[Å], vdw_radius[Å],
 *         ionization_energy[eV], electron_affinity[eV],
 *         valence,
 *         ground_config (populated at runtime by pt_electron_config),
 *         lj_epsilon[eV], lj_sigma[Å]
 *
 * Noble gas electron affinities set to 0 (they are negative).
 * Electronegativity = -1 for noble gases (undefined on Pauling scale).
 */
const Element PERIODIC_TABLE[MAX_ELEMENTS + 1] = {
    /* [0] unused sentinel */
    {0},

    /* Z=1  H  Hydrogen */
    { 1,"H","Hydrogen",       1.008,  2.20, 1.20, 0.31, 1.10, 13.598, 0.754,  1,
      {0}, 0.044*KCAL_MOL_TO_EV, 2.886 },

    /* Z=2  He Helium */
    { 2,"He","Helium",        4.003, -1.00, 1.40, 0.28, 1.40, 24.587, 0.000,  0,
      {0}, 0.056*KCAL_MOL_TO_EV, 2.362 },

    /* Z=3  Li Lithium */
    { 3,"Li","Lithium",       6.941,  0.98, 1.82, 1.28, 1.82,  5.392, 0.618,  1,
      {0}, 0.025*KCAL_MOL_TO_EV, 2.451 },

    /* Z=4  Be Beryllium */
    { 4,"Be","Beryllium",     9.012,  1.57, 1.53, 0.96, 1.53,  9.323, 0.000,  2,
      {0}, 0.085*KCAL_MOL_TO_EV, 2.745 },

    /* Z=5  B  Boron */
    { 5,"B","Boron",         10.811,  2.04, 1.92, 0.84, 1.92,  8.298, 0.280,  3,
      {0}, 0.180*KCAL_MOL_TO_EV, 4.083 },

    /* Z=6  C  Carbon */
    { 6,"C","Carbon",        12.011,  2.55, 1.70, 0.77, 1.70, 11.260, 1.262,  4,
      {0}, 0.105*KCAL_MOL_TO_EV, 3.851 },

    /* Z=7  N  Nitrogen */
    { 7,"N","Nitrogen",      14.007,  3.04, 1.55, 0.71, 1.55, 14.534, 0.000,  3,
      {0}, 0.069*KCAL_MOL_TO_EV, 3.660 },

    /* Z=8  O  Oxygen */
    { 8,"O","Oxygen",        15.999,  3.44, 1.52, 0.66, 1.52, 13.618, 1.461,  2,
      {0}, 0.060*KCAL_MOL_TO_EV, 3.500 },

    /* Z=9  F  Fluorine */
    { 9,"F","Fluorine",      18.998,  3.98, 1.47, 0.64, 1.47, 17.423, 3.401,  1,
      {0}, 0.050*KCAL_MOL_TO_EV, 3.364 },

    /* Z=10 Ne Neon */
    {10,"Ne","Neon",         20.180, -1.00, 1.54, 0.58, 1.54, 21.565, 0.000,  0,
      {0}, 0.042*KCAL_MOL_TO_EV, 3.243 },

    /* Z=11 Na Sodium */
    {11,"Na","Sodium",       22.990,  0.93, 2.27, 1.66, 2.27,  5.139, 0.548,  1,
      {0}, 0.030*KCAL_MOL_TO_EV, 2.983 },

    /* Z=12 Mg Magnesium */
    {12,"Mg","Magnesium",    24.305,  1.31, 1.73, 1.41, 1.73,  7.646, 0.000,  2,
      {0}, 0.111*KCAL_MOL_TO_EV, 3.021 },

    /* Z=13 Al Aluminium */
    {13,"Al","Aluminium",    26.982,  1.61, 1.84, 1.21, 1.84,  5.986, 0.441,  3,
      {0}, 0.505*KCAL_MOL_TO_EV, 4.499 },

    /* Z=14 Si Silicon */
    {14,"Si","Silicon",      28.086,  1.90, 2.10, 1.11, 2.10,  8.151, 1.385,  4,
      {0}, 0.402*KCAL_MOL_TO_EV, 4.295 },

    /* Z=15 P  Phosphorus */
    {15,"P","Phosphorus",    30.974,  2.19, 1.80, 1.07, 1.80, 10.487, 0.747,  3,
      {0}, 0.305*KCAL_MOL_TO_EV, 4.147 },

    /* Z=16 S  Sulphur */
    {16,"S","Sulphur",       32.065,  2.58, 1.80, 1.05, 1.80, 10.360, 2.077,  2,
      {0}, 0.274*KCAL_MOL_TO_EV, 4.035 },

    /* Z=17 Cl Chlorine */
    {17,"Cl","Chlorine",     35.453,  3.16, 1.75, 1.02, 1.75, 12.968, 3.613,  1,
      {0}, 0.227*KCAL_MOL_TO_EV, 3.947 },

    /* Z=18 Ar Argon */
    {18,"Ar","Argon",        39.948, -1.00, 1.88, 1.06, 1.88, 15.760, 0.000,  0,
      {0}, 0.185*KCAL_MOL_TO_EV, 3.868 },

    /* Z=19 K  Potassium — mass only stub */
    {19,"K","Potassium",     39.098,  0.82, 2.75, 2.03, 2.75,  4.341, 0.501,  1,
      {0}, 0.035*KCAL_MOL_TO_EV, 3.812 },

    /* Z=20 Ca Calcium */
    {20,"Ca","Calcium",      40.078,  1.00, 2.31, 1.76, 2.31,  6.113, 0.018,  2,
      {0}, 0.238*KCAL_MOL_TO_EV, 3.399 },
    /* Z=21..35 stubs with mass only */
    {21,"Sc","Scandium",     44.956,  1.36, 2.11, 1.70, 2.11,  6.561, 0.188,  3, {0}, 0.019*KCAL_MOL_TO_EV, 3.295},
    {22,"Ti","Titanium",     47.867,  1.54, 2.00, 1.60, 2.00,  6.828, 0.079,  4, {0}, 0.017*KCAL_MOL_TO_EV, 3.175},
    {23,"V","Vanadium",      50.942,  1.63, 1.92, 1.53, 1.92,  6.746, 0.525,  5, {0}, 0.016*KCAL_MOL_TO_EV, 3.144},
    {24,"Cr","Chromium",     51.996,  1.66, 1.85, 1.39, 1.85,  6.767, 0.666,  3, {0}, 0.015*KCAL_MOL_TO_EV, 3.023},
    {25,"Mn","Manganese",    54.938,  1.55, 1.79, 1.61, 1.79,  7.434, 0.000,  2, {0}, 0.013*KCAL_MOL_TO_EV, 2.961},
    {26,"Fe","Iron",         55.845,  1.83, 1.72, 1.52, 1.72,  7.902, 0.151,  3, {0}, 0.013*KCAL_MOL_TO_EV, 2.912},
    {27,"Co","Cobalt",       58.933,  1.88, 1.67, 1.50, 1.67,  7.881, 0.662,  3, {0}, 0.014*KCAL_MOL_TO_EV, 2.872},
    {28,"Ni","Nickel",       58.693,  1.91, 1.63, 1.24, 1.63,  7.640, 1.156,  2, {0}, 0.015*KCAL_MOL_TO_EV, 2.834},
    {29,"Cu","Copper",       63.546,  1.90, 1.40, 1.32, 1.40,  7.726, 1.235,  2, {0}, 0.005*KCAL_MOL_TO_EV, 3.495},
    {30,"Zn","Zinc",         65.38,   1.65, 1.39, 1.22, 1.39,  9.394, 0.000,  2, {0}, 0.124*KCAL_MOL_TO_EV, 2.763},
    {31,"Ga","Gallium",      69.723,  1.81, 1.87, 1.22, 1.87,  5.999, 0.300,  3, {0}, 0.415*KCAL_MOL_TO_EV, 4.383},
    {32,"Ge","Germanium",    72.630,  2.01, 2.11, 1.20, 2.11,  7.900, 1.233,  4, {0}, 0.379*KCAL_MOL_TO_EV, 4.280},
    {33,"As","Arsenic",      74.922,  2.18, 1.85, 1.19, 1.85,  9.815, 0.814,  3, {0}, 0.309*KCAL_MOL_TO_EV, 4.230},
    {34,"Se","Selenium",     78.971,  2.55, 1.90, 1.20, 1.90,  9.752, 2.021,  2, {0}, 0.291*KCAL_MOL_TO_EV, 4.205},
    {35,"Br","Bromine",      79.904,  2.96, 1.85, 1.20, 1.85, 11.814, 3.365,  1, {0}, 0.251*KCAL_MOL_TO_EV, 4.189},
    {36,"Kr","Krypton",      83.798, -1.00, 2.02, 1.16, 2.02, 14.000, 0.000,  0, {0}, 0.220*KCAL_MOL_TO_EV, 4.141},
};

/* ══════════════════════════════════════════════════════════════════════════
 * Implementation
 * ══════════════════════════════════════════════════════════════════════════ */

const Element *pt_element(int Z) {
    if (Z < 1 || Z > MAX_ELEMENTS) return NULL;
    return &PERIODIC_TABLE[Z];
}

const Element *pt_by_symbol(const char *symbol) {
    for (int Z = 1; Z <= 36; Z++) {
        if (strcmp(PERIODIC_TABLE[Z].symbol, symbol) == 0)
            return &PERIODIC_TABLE[Z];
    }
    return NULL;
}

/*
 * Fill an ElectronConfig for atomic number Z using the Madelung (diagonal)
 * rule. Fills orbitals in Madelung order until all Z electrons are placed.
 *
 * Note: Chromium (Z=24) and Copper (Z=29) and a few others are exceptions
 * to the strict Madelung rule; we handle those here explicitly.
 */
void pt_electron_config(int Z, ElectronConfig *cfg) {
    memset(cfg, 0, sizeof(*cfg));
    cfg->total_electrons = Z;

    /* Handle the well-known Madelung exceptions for 3d transition metals */
    /* We store any exception as a delta to apply after Madelung filling.  */
    int electrons_left = Z;

    for (int i = 0; i < MADELUNG_LEN && electrons_left > 0; i++) {
        int n = MADELUNG_N[i];
        int l = MADELUNG_L[i];
        int cap = orbital_cap(l);
        int fill = electrons_left < cap ? electrons_left : cap;
        cfg->config[n-1][l] += fill;
        electrons_left -= fill;
    }

    /* Madelung exceptions: swap 4s→3d for half-filled and full-d stability */
    if (Z == 24) { /* Cr: [Ar] 3d5 4s1 instead of 3d4 4s2 */
        cfg->config[3][0] -= 1;  /* remove one 4s */
        cfg->config[2][2] += 1;  /* add one 3d    */
    }
    if (Z == 29) { /* Cu: [Ar] 3d10 4s1 instead of 3d9 4s2 */
        cfg->config[3][0] -= 1;
        cfg->config[2][2] += 1;
    }

    /* Compute valence electrons (outermost shell sum) */
    int max_shell = 0;
    for (int s = 0; s < MAX_SHELLS; s++) {
        int shell_count = 0;
        for (int sub = 0; sub < 4; sub++) shell_count += cfg->config[s][sub];
        if (shell_count > 0) max_shell = s;
    }
    for (int sub = 0; sub < 4; sub++)
        cfg->valence_electrons += cfg->config[max_shell][sub];
}

/*
 * Write a human-readable electron configuration string into buf.
 * Example output: "1s2 2s2 2p2"
 */
void pt_config_string(const ElectronConfig *cfg, char *buf, int buflen) {
    static const char sub_labels[] = "spdf";
    int pos = 0;
    for (int n = 0; n < MAX_SHELLS && pos < buflen - 1; n++) {
        for (int l = 0; l < 4 && pos < buflen - 1; l++) {
            int count = cfg->config[n][l];
            if (count == 0) continue;
            int written = snprintf(buf + pos, buflen - pos, "%d%c%d ",
                                   n+1, sub_labels[l], count);
            if (written > 0) pos += written;
        }
    }
    if (pos > 0 && buf[pos-1] == ' ') buf[pos-1] = '\0';
    else buf[pos] = '\0';
}

void pt_print_element(int Z) {
    const Element *e = pt_element(Z);
    if (!e) { printf("Element Z=%d not found.\n", Z); return; }

    ElectronConfig cfg;
    pt_electron_config(Z, &cfg);
    char cfgstr[128];
    pt_config_string(&cfg, cfgstr, sizeof(cfgstr));

    printf("══════════════════════════════════════════\n");
    printf("  %s (%s)  Z=%d\n", e->name, e->symbol, e->Z);
    printf("══════════════════════════════════════════\n");
    printf("  Mass              : %.4f AMU\n",  e->mass);
    printf("  Electronegativity : %.2f (Pauling)\n", e->electronegativity);
    printf("  Atomic radius     : %.3f Å (vdW)\n", e->atomic_radius);
    printf("  Covalent radius   : %.3f Å\n",    e->covalent_radius);
    printf("  Ionisation energy : %.3f eV\n",   e->ionization_energy);
    printf("  Electron affinity : %.3f eV\n",   e->electron_affinity);
    printf("  Common valence    : %d\n",         e->valence);
    printf("  LJ ε              : %.5f eV\n",   e->lj_epsilon);
    printf("  LJ σ              : %.4f Å\n",    e->lj_sigma);
    printf("  Config            : %s\n",         cfgstr);
    printf("══════════════════════════════════════════\n");
}
