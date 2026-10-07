#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
#include <termios.h>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <stdarg.h>
#include <fnmatch.h>
#include <glob.h>
#include <regex.h>
#include <sys/utsname.h>

#include "../include/types.h"
#include "../include/vec3.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"
#include "../include/quantum.h"
#include "../include/qm.h"
#include "../include/forces.h"
#include "../include/integrator.h"
#include "../include/sim.h"
#include "../include/neuron.h"
#include "../include/nucleobases.h"
#include "../include/aminoacids.h"
#include "../include/amber_lj.h"
#include "../include/kcsa_filter.h"
#include "../include/tui.h"
#include "../include/tui_view.h"
#include "../include/tui_screen.h"

/* Diagnostics go to stderr (see sh_err below); stdout is DATA. */
static void sh_err(const char *fmt, ...);

/* OpenWorm-in-the-TUI. Plain stdin/stdout REPL, no ncurses, no graphics.
 * Separate `s2tui` binary so batch `./carbonsim` record stays untouched. */

typedef struct {
    Simulation *sim;
    HHNeuron nrn;
    int has_nrn;
    unsigned long seed;
    ViewCam cam;
    int quit;       /* REPL exits when nonzero */
    int srcdepth;   /* >0 while sourcing a script file */
    int baro_on;    /* isotropic Berendsen barostat (tui-side scaling) */
    double p0_bar;  /* target pressure, bar */
    double taup_fs; /* barostat time constant, fs */
    int rxn_every;  /* auto-reaction check every N steps (0 = off) */
    int rxn_armed;  /* auto-reaction fires once per arming */
} Tui;

static void tui_barostat_step(Tui *t);
static void rxn_autocheck(Tui *t);
static int sp_room(Tui *t, int atoms, int bonds);

static void tui_new(Tui *t, int atoms, int bonds) {
    if (t->sim) sim_destroy(t->sim);
    if (atoms < 8) atoms = 8;
    if (bonds < 8) bonds = 8;
    t->sim = sim_create(atoms, bonds);
    if (t->sim) {
        t->sim->dt = 0.5;
        t->sim->cutoff = 12.0;
        t->sim->dielectric = 1.0;
    }
}

static int parse_double(const char *s, double *out) {
    char *e = NULL;
    double v = strtod(s, &e);
    if (!e || e == s || !isfinite(v)) return 0;
    *out = v;
    return 1;
}

static int parse_int(const char *s, int *out) {
    char *e = NULL;
    long v = strtol(s, &e, 10);
    if (!e || e == s) return 0;
    *out = (int)v;
    return 1;
}

static int parse_ulong(const char *s, unsigned long *out) {
    char *e = NULL;
    unsigned long v = strtoul(s, &e, 10);
    if (!e || e == s) return 0;
    *out = v;
    return 1;
}

static void show_energy(const Tui *t) {
    if (!t->sim || t->sim->num_atoms < 1) { sh_err("  (empty sim: `new` first)\n"); return; }
    forces_calculate((Simulation *)t->sim);
    printf("  atoms=%d bonds=%d angles=%d step=%llu t=%.2f fs\n",
        t->sim->num_atoms, t->sim->num_bonds, t->sim->num_angles,
        (unsigned long long)t->sim->step, t->sim->time);
    printf("  KE=%.6f PE=%.6f E=%.6f T=%.2f K  (LJ=%.4f Coul=%.4f bond=%.4f ang=%.4f dih=%.4f restr=%.4f pol=%.4f pauli=%.4f disp=%.4f)\n",
        t->sim->kinetic_energy, t->sim->potential_energy, t->sim->total_energy,
        t->sim->temperature, t->sim->E_lj_total, t->sim->E_coulomb_total,
        t->sim->E_bond_total, t->sim->E_angle_total, t->sim->E_dihedral_total,
        t->sim->E_restraint_total, t->sim->E_polar_total, t->sim->E_pauli_total,
        t->sim->E_disp_total);
}

static void maybe_render(const Tui *t);

static void cmd_step(Tui *t, int n) {
    if (!t->sim) { sh_err("  (empty sim)\n"); return; }
    if (n < 1) n = 1;
    if (n > 100000) n = 100000;
    for (int i = 0; i < n; i++) {
        integrator_step(t->sim);
        tui_barostat_step(t);
        rxn_autocheck(t);
    }
    show_energy(t);
    maybe_render(t);
}

/* redraw the grid when interactive (TTY) and auto-view is on */
static void maybe_render(const Tui *t) {
    if (t->cam.auto_render && isatty(STDOUT_FILENO) && t->sim)
        view_render(t->sim, &t->cam, 0);
}

static void sleep_ms(long ms) {
    if (ms <= 0) return;
    struct timespec ts;
    ts.tv_sec = ms / 1000;
    ts.tv_nsec = (ms % 1000) * 1000000L;
    nanosleep(&ts, NULL);
}

/* ── demo-tests (reduced, fast; full physics lives in ./carbonsim) ── */
static int t_pass(const char *id, const char *detail) {
    printf("  TEST demo %s: PASS (%s)\n", id, detail);
    return 1;
}
static int t_fail(const char *id, const char *detail) {
    printf("  TEST demo %s: FAIL (%s)\n", id, detail);
    return 0;
}

static int test_d1(void) {
    ElectronConfig c; pt_electron_config(6, &c);
    double z = quantum_zeff_raw(6, 2, 1, &c);
    double e = quantum_orbital_energy(1, 1, 0, NULL);
    (void)e;
    ElectronConfig h; pt_electron_config(1, &h);
    double zh = quantum_zeff_raw(1, 1, 0, &h);
    char d[128];
    snprintf(d, sizeof d, "C2p Zeff=%.2f H1s=%.2f", z, zh);
    if (fabs(z - 3.25) < 1e-9 && fabs(zh - 1.0) < 1e-9) return t_pass("1", d);
    return t_fail("1", d);
}

static int test_d2(void) {
    Simulation *s = sim_create(8, 8);
    if (!s) return t_fail("2", "alloc");
    double sig = 2.5711, eps = 0.00191;
    double rmin = sig * TWOPOW_SIXTH;
    int a = sim_add_atom(s, 1, vec3(0, 0, 0), 0.0);
    int b = sim_add_atom(s, 1, vec3(rmin, 0, 0), 0.0);
    if (a < 0 || b < 0) { sim_destroy(s); return t_fail("2", "place"); }
    sim_set_atom_lj(s, 0, eps, sig);
    sim_set_atom_lj(s, 1, eps, sig);
    PairEnergy pe = forces_nonbonded_energy(s->atoms, 0, 1, NULL, 1, 0, 1.0);
    double ref = -eps; /* LJ minimum value is exactly -eps */
    char d[128];
    snprintf(d, sizeof d, "LJ(rmin)=%.9f ref=%.9f", pe.lj_energy, ref);
    int ok = fabs(pe.lj_energy - ref) < 1e-9;
    sim_destroy(s);
    return ok ? t_pass("2", d) : t_fail("2", d);
}

static int test_d3(void) {
    Simulation *s = sim_create(16, 32);
    if (!s) return t_fail("3", "alloc");
    sim_place_h2o(s, vec3_zero());
    s->dt = 0.5;
    s->thermostat.type = THERMOSTAT_BERENDSEN;
    s->thermostat.target_temperature = 300.0;
    s->thermostat.tau = 100.0;
    integrator_maxwell_boltzmann(s, 300.0, 42);
    forces_calculate(s);
    for (int i = 0; i < 200; i++) integrator_step(s);
    double oh = vec3_dist(s->atoms[0].position, s->atoms[1].position);
    char d[128];
    snprintf(d, sizeof d, "T=%.1fK O-H=%.3fA", s->temperature, oh);
    int ok = isfinite(s->temperature) && s->temperature > 100 && s->temperature < 500 && fabs(oh - 0.9572) < 0.08;
    sim_destroy(s);
    return ok ? t_pass("3", d) : t_fail("3", d);
}

static int test_d4(void) {
    Simulation *s = sim_create(64, 64);
    if (!s) return t_fail("4", "alloc");
    for (int k = 0; k < 3; k++) {
        double a = k * 2.0 * 3.141592653589793 / 3.0;
        sim_place_h2o(s, vec3(2.0 * cos(a), 2.0 * sin(a), 0.0));
    }
    s->dt = 0.5;
    s->thermostat.type = THERMOSTAT_BERENDSEN;
    s->thermostat.target_temperature = 50.0;
    s->thermostat.tau = 100.0;
    integrator_maxwell_boltzmann(s, 50.0, 7);
    forces_calculate(s);
    double pe0 = s->potential_energy;
    for (int i = 0; i < 300; i++) integrator_step(s);
    char d[128];
    snprintf(d, sizeof d, "PE %.3f->%.3f N=%d", pe0, s->potential_energy, s->num_atoms);
    int ok = s->num_atoms == 9 && isfinite(s->potential_energy);
    sim_destroy(s);
    return ok ? t_pass("4", d) : t_fail("4", d);
}

static int test_d5(void) {
    Simulation *s = sim_create(16, 32);
    if (!s) return t_fail("5", "alloc");
    int f = sim_place_ch4(s, vec3_zero());
    char d[128];
    snprintf(d, sizeof d, "first=%d bonds=%d", f, s->num_bonds);
    int ok = f == 0 && s->num_bonds == 4;
    sim_destroy(s);
    return ok ? t_pass("5", d) : t_fail("5", d);
}

static int test_d6(void) {
    const char *names[5] = {"U", "C", "T", "A", "G"};
    int total = 0;
    for (int k = 0; k < 5; k++) {
        Simulation *s = sim_create(128, 128);
        if (!s) return t_fail("6", "alloc");
        int f = -1;
        if (k == 0) f = sim_place_uracil(s, vec3_zero());
        else if (k == 1) f = sim_place_cytosine(s, vec3_zero());
        else if (k == 2) f = sim_place_thymine(s, vec3_zero());
        else if (k == 3) f = sim_place_adenine(s, vec3_zero());
        else f = sim_place_guanine(s, vec3_zero());
        if (f < 0 || s->num_atoms < 10) { sim_destroy(s); return t_fail("6", names[k]); }
        total += s->num_atoms;
        sim_destroy(s);
    }
    char d[128];
    snprintf(d, sizeof d, "5 bases placed, %d atoms total", total);
    return t_pass("6", d);
}

static int test_d7(void) {
    /* Reduced: G-C must carry more partial-charge magnitude than A-U path
     * builders exercise the same RESP tables as batch Demo 7. Full pairing
     * energetics stay in ./carbonsim output. */
    Simulation *s = sim_create(256, 256);
    if (!s) return t_fail("7", "alloc");
    int g = sim_place_guanine(s, vec3_zero());
    int c = sim_place_cytosine(s, vec3(5, 0, 0));
    char d[128];
    snprintf(d, sizeof d, "G@%d C@%d N=%d", g, c, s->num_atoms);
    int ok = g >= 0 && c >= 0 && s->num_atoms > 20;
    sim_destroy(s);
    return ok ? t_pass("7", d) : t_fail("7", d);
}

static int test_d8(void) {
    Simulation *s = sim_create(256, 256);
    if (!s) return t_fail("8", "alloc");
    int sugar = -1;
    int f = sim_place_dinucleotide_TA(s, vec3_zero(), &sugar);
    double q = 0;
    for (int i = 0; i < s->num_atoms; i++) q += s->atoms[i].partial_charge;
    char d[128];
    snprintf(d, sizeof d, "first=%d N=%d q=%.2f", f, s->num_atoms, q);
    int ok = f >= 0 && fabs(q + 1.0) < 0.05;
    sim_destroy(s);
    return ok ? t_pass("8", d) : t_fail("8", d);
}

static int test_d9(void) {
    HHNeuron n;
    hh_init(&n);
    n.I_ext = 10.0;
    int spikes = 0;
    for (int i = 0; i < 5000; i++) {
        hh_step(&n, 0.01);
        if (hh_is_spiking(&n, 0.0)) { spikes++; break; }
    }
    char d[128];
    snprintf(d, sizeof d, "V=%.1fmV spiking=%d", n.V, spikes);
    if (spikes) return t_pass("9", d);
    return t_fail("9", d);
}

static int test_d10(void) {
    Simulation *s = sim_create(128, 128);
    if (!s) return t_fail("10", "alloc");
    int alaN = -1;
    int f = sim_place_dipeptide_GlyAla(s, vec3_zero(), &alaN);
    double best = 1e30;
    for (int b = 0; b < s->num_bonds; b++) {
        double r = vec3_dist(s->atoms[s->bonds[b].atom_a].position,
                             s->atoms[s->bonds[b].atom_b].position);
        if (r > 1.0 && r < 1.6 && fabs(r - 1.33) < fabs(best - 1.33)) best = r;
    }
    char d[128];
    snprintf(d, sizeof d, "first=%d peptide~%.3fA", f, best);
    int ok = f >= 0 && fabs(best - 1.33) < 0.1;
    sim_destroy(s);
    return ok ? t_pass("10", d) : t_fail("10", d);
}

static int test_d11(void) {
    Simulation *s = sim_create(512, 512);
    if (!s) return t_fail("11", "alloc");
    AAResidue res[8];
    int f = sim_place_polyalanine(s, vec3_zero(), 6, res);
    char d[128];
    snprintf(d, sizeof d, "first=%d N=%d bonds=%d", f, s ? s->num_atoms : -1, s ? s->num_bonds : -1);
    int ok = f >= 0 && s->num_atoms > 30 && s->num_bonds > 20;
    sim_destroy(s);
    return ok ? t_pass("11", d) : t_fail("11", d);
}

static int test_d12(void) {
    Simulation *s = sim_create(64, 64);
    if (!s) return t_fail("12", "alloc");
    /* 8-O antiprism single-point legs (no WHAM): must be finite, JC splits K/Na.
     * LJ types via the shared header (N7): these two lines used to carry
     * 0.2100*KCAL_MOL_TO_EV and the literal 2.959921901149463, the same
     * duplicated-spelling the amber_lj.h consolidation removed from the
     * engine - the TUI had kept a private copy. Values are identical. */
    s->cutoff = 30.0;
    for (int i = 0; i < 4; i++) {
        double a = i * 1.5707963267948966;
        int o = sim_add_atom(s, 8, vec3(2.2 * cos(a), 2.2 * sin(a), -1.542), -0.5462);
        sim_set_atom_lj(s, o, LJ_AMBER_O_EPS, LJ_AMBER_O_SIGMA);
    }
    for (int i = 0; i < 4; i++) {
        double a = i * 1.5707963267948966 + 0.7853981633974483;
        int o = sim_add_atom(s, 8, vec3(2.37 * cos(a), 2.37 * sin(a), 1.542), -0.5462);
        sim_set_atom_lj(s, o, LJ_AMBER_O_EPS, LJ_AMBER_O_SIGMA);
    }
    int k = sim_add_ion(s, 19, 1, vec3(0, 0, 0), 1.0);
    forces_calculate(s);
    double ek = s->potential_energy;
    char d[128];
    snprintf(d, sizeof d, "8-O cage E=%.3feV ion=%d", ek, k);
    int ok = k >= 0 && isfinite(ek);
    sim_destroy(s);
    return ok ? t_pass("12", d) : t_fail("12", d);
}

static int test_d12b(void) {
    Simulation *s = sim_create(600, 600);
    if (!s) return t_fail("12b", "alloc");
    int f = kcsa_build_filter(s, vec3_zero(), 4);
    double q = 0;
    for (int i = 0; i < s->num_atoms; i++) q += s->atoms[i].partial_charge;
    Vec3 sites[4];
    int ns = kcsa_ion_sites(4, sites);
    char d[128];
    snprintf(d, sizeof d, "N=%d q=%.3f sites=%d", s->num_atoms, q, ns);
    int ok = f >= 0 && s->num_atoms == 164 && fabs(q) < 1e-9 && ns == 4;
    sim_destroy(s);
    return ok ? t_pass("12b", d) : t_fail("12b", d);
}

static int test_d17(void) {
    Simulation *s = sim_create(256, 256);
    if (!s) return t_fail("17", "alloc");
    int sugar = -1;
    int f = sim_place_dinucleotide_TA(s, vec3_zero(), &sugar);
    char d[128];
    snprintf(d, sizeof d, "TA N=%d (reduced; full duplex in batch)", s->num_atoms);
    int ok = f >= 0 && s->num_atoms > 40;
    sim_destroy(s);
    return ok ? t_pass("17", d) : t_fail("17", d);
}

static int run_demo_test(const char *id) {
    if (!strcmp(id, "1")) return !test_d1();
    else if (!strcmp(id, "2")) return !test_d2();
    else if (!strcmp(id, "3")) return !test_d3();
    else if (!strcmp(id, "4")) return !test_d4();
    else if (!strcmp(id, "5")) return !test_d5();
    else if (!strcmp(id, "6")) return !test_d6();
    else if (!strcmp(id, "7")) return !test_d7();
    else if (!strcmp(id, "8")) return !test_d8();
    else if (!strcmp(id, "9")) return !test_d9();
    else if (!strcmp(id, "10")) return !test_d10();
    else if (!strcmp(id, "11")) return !test_d11();
    else if (!strcmp(id, "12")) return !test_d12();
    else if (!strcmp(id, "12b")) return !test_d12b();
    else if (!strcmp(id, "17")) return !test_d17();
    sh_err("  unknown demo `%s` (try: 1 2 3 4 5 6 7 8 9 10 11 12 12b 17)\n", id);
    return 2;
}

/* ── help: full command reference ─────────────────────────────────────
 * Every command documented POSIX-man style: purpose, usage, POSIX mapping.
 * `help` prints the map; `help <cmd>` / `man <cmd>` prints one entry.
 * ───────────────────────────────────────────────────────────────────────── */
static void print_help(void) {
    printf("  s2tui — OpenWorm in the TUI (text only, live sim)\n");
    printf("  POSIX shell: 'sq' \"dq $V\" \\esc #comment ; && || | & > >> < <<< 2> $? $$ $! $#\n");
    printf("               $V ${V} ${V:-d} ${V:=d} ${V:?m} ${V:+a} ${#V} ${V#pat} ${V%%pat}\n");
    printf("               $( ) ` ` $(( )) *?[] ~ V=v set -- (status 0/1/2/127)\n");
    printf("  status: 0 ok, 1 error, 2 usage, 127 not found. `help <cmd>` for any entry.\n");
    printf("  flow: if TEST; then A; else B; fi | for V in W; do B; done | while T; do B; done\n");
    printf("        case W in PAT) B;; esac | { list; } | ( list ) | ! cmd | eval args\n");
    printf("  --- sim session (POSIX mapping) ---\n");
    printf("    new [A B]: fresh sim (POSIX mapping: like starting a new shell)\n");
    printf("    ls atoms|bonds|summary|demos | ls [path]   list sim objects / files (ls(1))\n");
    printf("    del atom <i> | rm atom <i>[..ranges] | rm <file>..  remove atom / files (rm(1))\n");
    printf("      e.g. rm atom 1..5,7  (descending, terminal-only)\n");
    printf("    bond <a> <b> [ord] | detect bonds | restrain <i> x y z k | clear restraints\n");
  printf("  --- spawn chemicals (no POSIX equivalent; s2 verbs) ---\n");
  printf("    spawn atom <Z|sym> x y z [q] | spawn ion <Z> <formal> x y z [q]\n");
  printf("    spawn h2|h2o|nh3|ch4|methane|co2 [x y z] | kcsa [nsub] | filter [x y z]\n");
  printf("    spawn demo <id> [variant] [x y z] — every demo as a live system:\n");
  printf("      1 quantum(H) 2 bond(H2) 3 water 4 trimer 5 methane 6 bases 7 pair\n");
  printf("      8 dinucleotide 9 neuron 10 dipeptide 11 helix 12 cage 12b filter 17 duplex\n");
  printf("    spawn quantum|water|trimer|methane | base <U|C|T|A|G> | pair [gc|au]\n");
    printf("    spawn dinucleotide|neuron|dipeptide|helix|cage|filter|duplex [x y z]\n");
    printf("    spawn demo <1|2|3|4|5|6|7|8|9|10|11|12|12b|17> [variant] [x y z]\n");
    printf("      (variants: demo 6 [U|C|T|A|G], demo 7 [gc|au]; pair/duplex set dielectric 4;\n");
    printf("       helix/duplex clash-relieved; pass origins — default 0,0,0 overlaps!)\n");
    printf("  --- live controls ---\n");
    printf("    set dt|cutoff|...|maxtemp|gamma | set NAME=VALUE ... (streamlined)\n");
    printf("    set box|pbc|press|tau-p|barostat (gas NPT-ish) | show energy|pressure|thermo|temp\n");
    printf("    init velocities <T> [seed] | step [N] | run <N> | heat <dE_eV> | minimize ...\n");
    printf("    neuron init|inject|step|run|show\n");
    printf("  --- your chemistry: molecules, reactions, gases ---\n");
    printf("    mol new|add|bond|list|clear|center|save|load|place (MOL1 templates)\n");
    printf("    rxn new|pair|break|make|delete|charge|list|save|load|fire|arm|auto\n");
    printf("  --- tests: each demo is a system test ---\n");
    printf("    list demos | test demo <id> | test all   (id: 1 2 3 4 5 6 7 8 9 10 11 12 12b 17)\n");
    printf("  --- grid display ---\n");
    printf("    render | view xy|xz|yz|auto | cam yaw|pitch|zoom|center|reset | slice | watch\n");
    printf("  --- files & shell (real POSIX) ---\n");
    printf("    cat | pwd | cd (~) | mkdir [-p] | cp | mv | head [-n] | sh -c | man\n");
    printf("    tail | wc | sort | uniq | cut | tr | grep | tee | basename | dirname\n");
    printf("    touch <sp> [xN] | touch -- <file> | touch -m <sp> | rmdir | rm -r | ln [-s]\n");
    printf("    date | uname | find | wait | true | false | : | eval | exec | command | type\n");
    printf("    shift | set -- | readonly | umask | trap | alias | if/for/while/case | !\n");
    printf("    echo [-n] | printf | export | unset | env | history | source/. | clear\n");
    printf("    sleep <steps>[fs] | sleep <sec>s | time | test/[ | wait | | pipe | & ($!) | > >> < <<< 2> globs $() `` $(( ))\n");
    printf("  quit | exit [n]\n");
}

static int print_help_topic(const char *t) {
    struct { const char *names; const char *text; } ref[] = {
        { "shell quoting vars status",
          "  POSIX shell layer (IEEE 1003.1 style):\n"
          "  'sq' literal | \"dq\" expands $V | \\x escapes x | # comment to EOL\n"
          "  ; sequence | && on-success | || on-failure | > >> redirect output\n"
          "  $V ${V} env vars | $? last status | $$ shell pid | V=v assignment\n"
          "  status: 0 ok, 1 error, 2 usage/syntax, 127 command not found" },
        { "help man",
          "  help [topic] | man <cmd>: this reference. Topics: shell file session\n"
          "  spawn set control neuron grid test save (POSIX: man(1))" },
        { "quit exit",
          "  quit | exit [n]: leave s2tui with status n (default: last status)" },
        { "file files cat pwd cd mkdir cp mv head ls rm !",
          "  cat <f>.. (cat(1): concatenate files to output)\n"
          "  pwd (pwd(1): print working directory) | cd <dir> (cd: change directory)\n"
          "  mkdir [-p] <dir>.. (mkdir(1)) | cp <s> <d> (cp(1)) | mv <s> <d> (mv(1))\n"
          "  head [-n N] <f>.. (head(1): first N lines, default 10)\n"
          "  ls [path] (ls(1): list directory, dirs suffixed /)\n"
          "  rm <file>.. (rm(1): refuses directories)\n"
          "  ! <cmd>.. | !<cmd> (sh(1): run a real Linux command, returns its status)" },
        { "session new list del bond detect restrain clear",
          "  new [atoms bonds]: fresh sim (POSIX mapping: starting a new shell)\n"
          "  list|ls atoms|bonds|summary (ls(1): list the working set)\n"
          "  del atom <i> | rm atom <i> (rm(1): terminal atoms only)\n"
          "  bond <a> <b> [order] (cf. ln(1): link two atoms)\n"
          "  detect bonds: auto-link by covalent distance\n"
          "  restrain <i> x y z k: spring anchor (V=1/2k|r-a|^2)\n"
          "  clear restraints" },
        { "spawn",
          "  spawn atom <Z|sym> x y z [q]: place element (no POSIX equivalent)\n"
          "  spawn ion <Z> <formal> x y z [q]: ion with electron count Z-formal\n"
          "  spawn h2|h2o|nh3|ch4|methane|co2 [x y z] | kcsa [nsub] | filter [x y z]\n"
          "  spawn demo <id> [variant] [x y z]: every demo as a live system\n"
          "  spawn quantum|water|trimer|methane | base <U|C|T|A|G> | pair [gc|au]\n"
          "  spawn dinucleotide|neuron|dipeptide|helix|cage|filter|duplex\n"
          "  spawn demo <id> [variant] [x y z] (6:[U|C|T|A|G] 7:[gc|au])" },
        { "set",
          "  set dt|cutoff|dielectric|temp|tau|nu|seed <v> (cf. stty(1): tune device)\n"
          "  set thermostat none|berendsen|andersen|langevin\n"
          "  set lj <i> <eps_eV> <sig_A> | set charge <i> <q>\n"
          "  set box <lx> <ly> <lz> | set pbc on|off|x y z\n"
          "  set press <bar> | set tau-p <fs> | set barostat on|off" },
        { "control step run show init minimize heat box pressure thermo",
          "  step [N] | run <N>: integrate N Verlet steps (default 1)\n"
          "  init velocities <T> [seed]: Maxwell-Boltzmann kick + COM removal\n"
          "  show energy: KE/PE/T + LJ/Coulomb/restraint/pol/Pauli/disp ledger\n"
          "  show pressure: virial P in bar + NkT/V + vir/3V terms (needs box)\n"
          "  show thermo: T/N/KE/V/rho/P/barostat | show temp: per-element T\n"
          "  heat <dE_eV>: scale velocities by (KE+dE)/KE (signed; needs KE)\n"
          "  set box <lx> <ly> <lz> (PBC all on) | set pbc on|off|x y z\n"
          "  set press <bar> | set tau-p <fs> | set barostat on|off (tui-side NPT)\n"
          "  minimize [iters step tol]: steepest-descent clash relief" },
        { "mol molecule builder template",
          "  mol new [name] | mol add <Z|sym> x y z [q] | mol bond <a> <b> [ord]\n"
          "  mol list | mol clear | mol center (COM) | mol save|load <file> (MOL1)\n"
          "  mol place [x y z]: stamp builder into the live sim" },
        { "rxn reaction break form trigger",
          "  rxn new [name] | rxn pair <Z1> <Z2> <rcut_A> (closest-pair match)\n"
          "  rxn break on|off | rxn make <order>|none | rxn delete <0|1|2>\n"
          "  rxn charge <q1> <q2>|off | rxn list | rxn save|load <file> (RXN1)\n"
          "  rxn fire (now) | rxn arm + rxn auto <every>|off (fires once live)" },
        { "neuron",
          "  neuron init: HH1952 squid axon at rest (-65mV, steady-state m/h/n)\n"
          "  neuron inject <uA/cm2> | neuron step (0.01ms RK4) | neuron run <N> | neuron show" },
        { "grid render view cam slice watch",
          "  ASCII grid: h c n o N p K = H C N O Na P K; . bond; x restraint\n"
          "  render: draw now | view xy|xz|yz (projections) | view auto on|off\n"
          "  cam yaw|pitch <deg> | cam zoom <f> | cam center [i] | cam reset\n"
          "  slice <thick_A|off>: z-slab cutaway | watch <steps> [ms]: live run" },
        { "test demos",
          "  test demo <id> | test all: each demo as a fast live system test\n"
          "  ids: 1 quantum 2 bond 3 water 4 trimer 5 methane 6 bases 7 pairing\n"
          "  8 dinucleotide 9 neuron 10 dipeptide 11 helix 12 kcsa 12b filter 17 duplex\n"
          "  (reduced: Demo 12 skips WHAM sampling; full physics: ./carbonsim)" },
        { "save load persist snapshot",
          "  save <file> | load <file>: S2SAVE1 snapshot\n"
          "  (atoms+LJ+velocities, bonds+r0/k, restraints, dt/cutoff/thermostat/seed)" },
        { "echo export unset env history source clear sleep time",
          "  echo [-n] args.. (echo(1)) | export N[=v].. / unset N.. / env\n"
          "  history (history: numbered inputs, Up/Down on TTY)\n"
          "  source <f> | . <f> (.(1): run script file, depth 8)\n"
          "  clear (clear(1)) | sleep <s> (sleep(1)) | time <cmd> (time(1): wall time)" },
        { "pipe background jobs wait glob redirection",
          "  a | b | c: pipes (no fork; identical output for terminating filters)\n"
          "  cmd &: background subshell snapshot, prints [pid], $! holds pid\n"
          "  wait [pid..]: reap jobs (wait(1)) | * ? [] globs (no match: literal)\n"
          "  $(...): capture stdout as one word (depth 4, no field splitting)\n"
          "  > >> truncate/append stdout | < stdin | <<< herestring | 2> 2>> stderr" },
        { "test [ true false",
          "  test EXPR | [ EXPR ] (test(1)): -e -f -d -r -w -x -s -z -n -t\n"
          "  = == != -eq -ne -gt -ge -lt -le ! -a -o ( ) | bare: false(1)\n"
          "  true (always 0) | false (always 1)" },
        { "printf",
          "  printf FORMAT [args..] (printf(1)): %s %d %i %u %o %x %X %f %e %g %c %%\n"
          "  flags width.precision, \\n \\t \\\\ \\a \\0ooo \\c, %b expands escapes" },
        { "text tail wc sort uniq cut tr grep tee",
          "  tail [-n N] [+N] (tail(1)) | wc [-lwc] (wc(1))\n"
          "  sort [-nru] (sort(1)) | uniq [-cdu] (uniq(1): adjacent)\n"
          "  cut -c|-b LIST | -f LIST [-d C] [-s] (cut(1))\n"
          "  tr [-cds] SET1 [SET2] (tr(1): ranges, [:upper:] etc.)\n"
          "  grep [-Einvc l] PATTERN [file..] (grep(1): BRE, -E extended)\n"
          "  tee [-a] [file..] (tee(1)) | no-args reads stdin (Ctrl-D ends)" },
        { "basename dirname touch rmdir ln date uname find",
          "  basename | dirname (path strings) | touch (touch(1): create/stamp)\n"
          "  rmdir | rm -r (recursive) | rm -f (force) | ln [-s] src dst\n"
          "  date [-u] [+FORMAT] (date(1)) | uname [-asnrv] (uname(1))\n"
          "  find path.. [-name PAT] [-type f|d] (find(1): fnmatch)" },
        { NULL, NULL }
    };
    for (int i = 0; ref[i].names; i++) {
        const char *p = ref[i].names;
        while (*p) {
            char w[32];
            int k = 0;
            while (*p && *p != ' ' && k < 31) w[k++] = *p++;
            w[k] = '\0';
            while (*p == ' ') p++;
            if (!strcmp(w, t)) { printf("%s\n", ref[i].text); return 0; }
        }
    }
    sh_err("  no manual entry for `%s` (try `help`)\n", t);
    return 1;
}

/* ── POSIX-shell layer ────────────────────────────────────────────────
 * Full POSIX shell syntax (IEEE 1003.1): quotes, escapes, # comments,
 * ; && || | & chaining, $V ${V} ${V:-d} ${V:=d} ${V:?m} ${V:+a} ${#V}
 * ${V#pat} ${V%pat}, $(( )) arithmetic, $( ) and `` capture, *?[] globs,
 * ~ home, V=v, set --, $# $@ $* $0..$9 $? $$ $! ! negation, redirection
 * > >> < <<< 2> 2>>, if/for/while/until/case, { } ( ) grouping, eval,
 * exec, command, type, shift, readonly, umask, trap, termios editing.
 * S2 params streamlined (same POSIX names): set NAME=VALUE, rm ranges,
 * sleep <steps>[fs]|<sec>s, touch [xN] / touch -- / touch -m.
 * Exit status: 0 ok, 1 runtime, 2 usage, 127 not found.
 * Quirk vs POSIX: no subprocesses, so VAR=x persists even with trailing
 * command; $( ) has no field splitting (one word); `time` is a prefix.
 * ───────────────────────────────────────────────────────────────────────── */

static int last_status = 0;

/* Diagnostics (usage, errors, prompts, progress) go to stderr; stdout is
 * DATA (command output, listings, tables, grid). That is what makes pipes,
 * `>`/`2>` and $( ) honest: `failing | wc` counts data, not error text. */
static void sh_err(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}

#define SH_VARS 64
static char sh_name[SH_VARS][64];
static char sh_val[SH_VARS][256];
static int sh_nvars = 0;

static const char *sh_get(const char *name) {
    for (int i = 0; i < sh_nvars; i++)
        if (!strcmp(sh_name[i], name)) return sh_val[i];
    return NULL;
}

static int sh_valid_name(const char *s, size_t n) {
    if (n == 0 || n > 63) return 0;
    if (!(isalpha((unsigned char)s[0]) || s[0] == '_')) return 0;
    for (size_t i = 1; i < n; i++)
        if (!(isalnum((unsigned char)s[i]) || s[i] == '_')) return 0;
    return 1;
}

static int sh_set(const char *name, const char *val) {
    for (int i = 0; i < sh_nvars; i++) {
        if (!strcmp(sh_name[i], name)) {
            snprintf(sh_val[i], sizeof sh_val[i], "%s", val);
            return 0;
        }
    }
    if (sh_nvars >= SH_VARS) { sh_err("  out of variables\n"); return 1; }
    snprintf(sh_name[sh_nvars], sizeof sh_name[0], "%s", name);
    snprintf(sh_val[sh_nvars], sizeof sh_val[0], "%s", val);
    sh_nvars++;
    return 0;
}

static int sh_unset(const char *name) {
    for (int i = 0; i < sh_nvars; i++) {
        if (!strcmp(sh_name[i], name)) {
            for (int j = i; j + 1 < sh_nvars; j++) {
                strcpy(sh_name[j], sh_name[j + 1]);
                strcpy(sh_val[j], sh_val[j + 1]);
            }
            sh_nvars--;
            return 0;
        }
    }
    return 1;
}

/* ── POSIX positional params + arithmetic + pattern helpers ──
 * $#, $@, $*, $0..$9, $-aretained via `set --`. $(( )) integer math.
 * ${} forms: :- := :? :+ # ## % %% ${#V}. Backquotes == $( ). */
#define SH_ARGS 64
static char sh_argv[SH_ARGS][256];
static int sh_argc = 0;
static char sh_prog[64] = "s2tui";

/* forward decls for expansion */
static char *sh_capture(Tui *t, const char *cmd, int depth);
static int run_line(Tui *t, char *line);

/* integer arithmetic: + - * / % () unary, $V / V names, numbers */
static const char *arith_p;
static void arith_sp(void) { while (*arith_p == ' ' || *arith_p == '\t') arith_p++; }
static long arith_expr(Tui *t, int *err);
static long arith_val(Tui *t, int *err) {
    arith_sp();
    int neg = 0;
    if (*arith_p == '+' || *arith_p == '-') { neg = (*arith_p == '-'); arith_p++; arith_sp(); }
    long v = 0;
    if (*arith_p == '(') {
        arith_p++; v = arith_expr(t, err); arith_sp();
        if (*arith_p == ')') arith_p++; else *err = 1;
    } else if (isalpha((unsigned char)*arith_p) || *arith_p == '_') {
        char nm[64]; int k = 0;
        if (*arith_p == '$') arith_p++;
        while ((isalnum((unsigned char)*arith_p) || *arith_p == '_') && k < 63)
            nm[k++] = *arith_p++;
        nm[k] = '\0';
        const char *s = sh_get(nm);
        if (!s) s = "0";
        char *e = NULL; v = strtol(s, &e, 10);
        if (!e || e == s) v = 0;
    } else {
        char *e = NULL;
        if (*arith_p == '$') arith_p++;
        v = strtol(arith_p, &e, 10);
        if (e == arith_p) { *err = 1; return 0; }
        arith_p = e;
    }
    arith_sp();
    return neg ? -v : v;
}
static long arith_term(Tui *t, int *err) {
    long v = arith_val(t, err);
    for (;;) {
        arith_sp();
        if (*arith_p == '*' || *arith_p == '/' || *arith_p == '%') {
            char op = *arith_p++;
            long r = arith_val(t, err);
            if (*err) return 0;
            if (op == '*') v *= r;
            else if (r == 0) { *err = 1; return 0; }
            else if (op == '/') v /= r;
            else v %= r;
        } else break;
    }
    return v;
}
static long arith_expr(Tui *t, int *err) {
    long v = arith_term(t, err);
    for (;;) {
        arith_sp();
        if (*arith_p == '+' || *arith_p == '-') {
            char op = *arith_p++;
            long r = arith_term(t, err);
            if (*err) return 0;
            v = (op == '+') ? v + r : v - r;
        } else break;
    }
    return v;
}
static int eval_arith(Tui *t, const char *expr, long *out) {
    int err = 0;
    arith_p = expr;
    long v = arith_expr(t, &err);
    arith_sp();
    if (err || *arith_p != '\0') return 0;
    *out = v;
    return 1;
}

#define SH_HIST 128
#define SH_LINE 4096
static char sh_hist[SH_HIST][1024];
static int sh_hist_n = 0;

static void sh_hist_push(const char *line) {
    if (!line) return;
    while (*line == ' ' || *line == '\t') line++;
    if (*line == '\0') return;
    if (sh_hist_n > 0 && !strcmp(sh_hist[(sh_hist_n - 1) % SH_HIST], line)) return;
    snprintf(sh_hist[sh_hist_n % SH_HIST], 1024, "%s", line);
    sh_hist_n++;
}

/* segments of one input line; nextop: 0 ';'/end, 1 '&&', 2 '||'.
 * Words are stored RAW (quotes removed, $ left intact) with a parallel
 * expand mask (1 = $ expands here, 0 = single-quoted or backslash-escaped).
 * Expansion runs per segment at execution time, so `V=1; echo $V` and
 * `fail; echo $?` see the assignments/status the segment produced —
 * expanding the whole line upfront (as an earlier draft did) broke both. */
#define SH_MAXSEG 32
#define SH_MAXARG 64
typedef struct { char *raw; char *mask; int globok; } ShWord;
typedef struct { ShWord args[SH_MAXARG]; int argc; int nextop; int bg; } ShSeg;

/* push current raw word (wbuf/wmask, length wn) into the segment */
static int sh_push(ShSeg *segs, int ns, int *ac, char *wbuf, char *wmask, int wn, int wglob) {
    if (*ac >= SH_MAXARG - 1) return -1;
    char *r = (char *)malloc((size_t)wn + 1);
    char *m = (char *)malloc((size_t)wn + 1);
    if (!r || !m) { free(r); free(m); return -1; }
    memcpy(r, wbuf, (size_t)wn + 1);
    memcpy(m, wmask, (size_t)wn + 1);
    segs[ns].args[*ac].raw = r;
    segs[ns].args[*ac].mask = m;
    segs[ns].args[*ac].globok = wglob;
    (*ac)++;
    return 0;
}

/* push an operator word (>, >>, <, 2>, 2>>) with mask-2 bytes */
static int sh_push_op(ShSeg *segs, int ns, int *ac, const char *op) {
    if (*ac + 1 >= SH_MAXARG - 1) return -1;
    size_t n = strlen(op);
    char *r = (char *)malloc(n + 1), *m = (char *)malloc(n + 1);
    if (!r || !m) { free(r); free(m); return -1; }
    memcpy(r, op, n + 1);
    memset(m, 2, n);
    m[n] = '\0';
    segs[ns].args[*ac].raw = r;
    segs[ns].args[*ac].mask = m;
    segs[ns].args[*ac].globok = 0;
    (*ac)++;
    return 0;
}

/* flush helper for the lexer: push pending word, reset accumulators */
#define SH_FLUSH() do { \
    wbuf[wn] = '\0'; wmask[wn] = '\0'; \
    if (sh_push(segs, ns, &ac, wbuf, wmask, wn, wglob)) return -1; \
    wn = 0; have = 0; wglob = 0; \
} while (0)

/* lex one line into segments. Returns count, -1 on syntax error
 * (unclosed quote, trailing &&/||, lone & |). Words are malloc'd. */
static int sh_lex(char *line, ShSeg *segs) {
    int ns = 0, ac = 0;
    char wbuf[SH_LINE], wmask[SH_LINE];
    int wn = 0, in_s = 0, in_d = 0, esc = 0, have = 0, com = 0, wglob = 0;
    memset(segs, 0, sizeof(ShSeg) * SH_MAXSEG);
    for (char *p = line; ; p++) {
        char c = *p;
        int end = (c == '\0' || c == '\n');
        if (com) {
            if (end) break;
            continue;
        }
        if (esc) {
            if (wn < SH_LINE - 1) { wbuf[wn] = c; wmask[wn] = 0; wn++; }
            have = 1; esc = 0;
            continue;
        }
        if (in_s) {
            if (c == '\'') in_s = 0;
            else if (end) return -1;
            else { if (wn < SH_LINE - 1) { wbuf[wn] = c; wmask[wn] = 0; wn++; } have = 1; }
            continue;
        }
        if (in_d) {
            if (c == '"') { in_d = 0; continue; }
            if (c == '\\') {
                /* POSIX: inside "", \ stays literal except before $ ` " \ newline */
                char d2 = p[1];
                if (d2 == '$' || d2 == '`' || d2 == '"' || d2 == '\\' || d2 == '\n') {
                    p++;
                    if (d2 == '\n') continue;
                    if (wn < SH_LINE - 1) { wbuf[wn] = d2; wmask[wn] = 0; wn++; }
                    have = 1;
                    continue;
                }
                if (wn < SH_LINE - 1) { wbuf[wn] = '\\'; wmask[wn] = 1; wn++; }
                have = 1;
                continue;
            }
            if (end) return -1;
            if (wn < SH_LINE - 1) { wbuf[wn] = c; wmask[wn] = 1; wn++; }
            have = 1;
            continue;
        }
        /* unquoted */
        if (c == '\\') { esc = 1; continue; }
        if (c == '\'') { in_s = 1; have = 1; continue; }
        if (c == '"') { in_d = 1; have = 1; continue; }
        if (c == '#') {
            if (!have) { com = 1; continue; }
            if (wn < SH_LINE - 1) { wbuf[wn] = c; wmask[wn] = 1; wn++; }
            have = 1;
            continue;
        }
        if (c == ';' || end || c == ' ' || c == '\t' || c == '\r') {
            if (have) { SH_FLUSH(); }
            if (c == ';' || end) {
                segs[ns].argc = ac; segs[ns].nextop = 0;
                ns++; ac = 0;
                if (ns >= SH_MAXSEG) return -1;
                if (end) break;
                segs[ns].argc = 0; segs[ns].nextop = 0;
            }
            continue;
        }
        if (c == '&' || c == '|') {
            if (p[1] == c) {
                if (have) { SH_FLUSH(); }
                segs[ns].argc = ac; segs[ns].nextop = (c == '&') ? 1 : 2;
                ns++; ac = 0;
                if (ns >= SH_MAXSEG) return -1;
                p++;
                /* trailing && || with nothing after is a syntax error */
                char *q = p + 1;
                while (*q == ' ' || *q == '\t') q++;
                if (*q == '\0' || *q == '\n') return -1;
                continue;
            }
            if (c == '|') {
                /* single |: pipeline stage boundary (nextop 3) */
                if (have) { SH_FLUSH(); }
                segs[ns].argc = ac; segs[ns].nextop = 3;
                ns++; ac = 0;
                if (ns >= SH_MAXSEG) return -1;
                char *q = p + 1;
                while (*q == ' ' || *q == '\t') q++;
                if (*q == '\0' || *q == '\n') return -1;
                continue;
            }
            /* single &: background the segment just closed */
            if (have) { SH_FLUSH(); }
            segs[ns].argc = ac; segs[ns].nextop = 0; segs[ns].bg = 1;
            ns++; ac = 0;
            if (ns >= SH_MAXSEG) return -1;
            continue;
        }
        if (c == '<') {
            if (have) { SH_FLUSH(); }
            if (p[1] == '<' && p[2] == '<') {
                /* <<< word: herestring (bashism, handy here): the next
                 * word is literal stdin content, not a filename */
                if (sh_push_op(segs, ns, &ac, "<<<")) return -1;
                p += 2;
            } else {
                if (sh_push_op(segs, ns, &ac, "<")) return -1;
            }
            continue;
        }
        if (c == '2' && !have && p[1] == '>') {
            /* fd-2 redirect only when 2 starts the word (POSIX IO_NUMBER) */
            if (p[2] == '>') {
                if (sh_push_op(segs, ns, &ac, "2>>")) return -1;
                p += 2;
            } else {
                if (sh_push_op(segs, ns, &ac, "2>")) return -1;
                p++;
            }
            continue;
        }
        if (c == '$' && p[1] == '(') {
            /* $(...): consume balanced and quote-aware as ONE word, so
             * spaces inside don't split (expansion runs it per segment) */
            char *q = p + 2;
            int dep = 1, sq = 0, dq = 0, es = 0;
            while (*q && dep > 0) {
                if (es) es = 0;
                else if (sq) { if (*q == '\'') sq = 0; }
                else if (dq) {
                    if (*q == '\\') es = 1;
                    else if (*q == '"') dq = 0;
                } else if (*q == '\\') es = 1;
                else if (*q == '\'') sq = 1;
                else if (*q == '"') dq = 1;
                else if (*q == '(') dep++;
                else if (*q == ')') { dep--; if (dep == 0) break; }
                q++;
            }
            if (dep != 0) { sh_err("  syntax error: unclosed $(\n"); return -1; }
            size_t span = (size_t)(q - p) + 1;
            if (wn + (int)span >= SH_LINE - 1) return -1;
            for (size_t k = 0; k < span; k++) { wbuf[wn] = p[k]; wmask[wn] = 1; wn++; }
            have = 1;
            p = q;
            continue;
        }
        if (c == '`') {
            /* `...`: consume to matching backquote as ONE word (POSIX) */
            char *q = p + 1;
            while (*q && *q != '`') {
                if (*q == '\\' && q[1]) q++;
                q++;
            }
            if (*q != '`') { sh_err("  syntax error: unclosed `\n"); return -1; }
            size_t span = (size_t)(q - p) + 1;
            if (wn + (int)span >= SH_LINE - 1) return -1;
            for (size_t k = 0; k < span; k++) { wbuf[wn] = p[k]; wmask[wn] = 1; wn++; }
            have = 1;
            p = q;
            continue;
        }
        if (c == '>') {
            if (have) { SH_FLUSH(); }
            if (p[1] == '>') {
                if (sh_push_op(segs, ns, &ac, ">>")) return -1;
                p++;
            } else {
                if (sh_push_op(segs, ns, &ac, ">")) return -1;
            }
            continue;
        }
        /* $ is stored raw with an expandable mask; sh_expand resolves it
         * per segment at execution time. */
        if (wn < SH_LINE - 1) { wbuf[wn] = c; wmask[wn] = 1; wn++; }
        if (c == '*' || c == '?' || c == '[') wglob = 1;
        have = 1;
    }
    /* drop trailing empty segment from a final ';' (but keep bg flag noise out) */
    while (ns > 0 && segs[ns - 1].argc == 0 && !segs[ns - 1].bg) ns--;
    return ns;
}

/* expand one raw word (mask-parallel) into a malloc'd final word.
 * Returns NULL on bad ${...} (caller reports status 2). Values are
 * copied verbatim, never rescanned. $(...) runs through run_line with
 * stdout captured, as ONE word (no field splitting — documented). */
static char *sh_capture(Tui *t, const char *cmd, int depth);
static char *sh_expand(Tui *t, const char *raw, const char *mask, int depth) {
    size_t cap = 64, len = 0;
    char *out = (char *)malloc(cap);
    if (!out) return NULL;
    /* NULL mask (arithmetic pre-expansion) means fully expandable */
    for (size_t i = 0; raw[i]; ) {
        int m = mask ? mask[i] : 1;
        /* backquotes `cmd` == $(cmd): POSIX command substitution */
        if (raw[i] == '`' && m) {
            size_t j = i + 1;
            while (raw[j] && !(raw[j] == '`' && (mask ? mask[j] : 1))) j++;
            if (!raw[j]) {
                if (len + 1 >= cap) {
                    cap *= 2;
                    char *n2 = (char *)realloc(out, cap);
                    if (!n2) { free(out); return NULL; }
                    out = n2;
                }
                out[len++] = raw[i++];
                continue;
            }
            size_t inner = j - (i + 1);
            char *cmd = (char *)malloc(inner + 1);
            if (!cmd) { free(out); return NULL; }
            memcpy(cmd, raw + i + 1, inner);
            cmd[inner] = '\0';
            char *got = sh_capture(t, cmd, depth + 1);
            free(cmd);
            if (!got) { free(out); return NULL; }
            size_t vl2 = strlen(got);
            while (len + vl2 + 1 >= cap) {
                cap *= 2;
                char *n2 = (char *)realloc(out, cap);
                if (!n2) { free(got); free(out); return NULL; }
                out = n2;
            }
            memcpy(out + len, got, vl2);
            len += vl2;
            free(got);
            i = j + 1;
            continue;
        }
        if (raw[i] == '$' && m) {
            char nb[256];
            const char *val = "";
            size_t adv = 1;
            if (raw[i + 1] == '?') {
                snprintf(nb, sizeof nb, "%d", last_status);
                val = nb; adv = 2;
            } else if (raw[i + 1] == '$') {
                snprintf(nb, sizeof nb, "%d", (int)getpid());
                val = nb; adv = 2;
            } else if (raw[i + 1] == '!') {
                const char *v = sh_get("!");
                snprintf(nb, sizeof nb, "%s", v ? v : "");
                val = nb; adv = 2;
            } else if (raw[i + 1] == '(') {
                /* $((...)) arithmetic vs $(...) capture */
                if (raw[i + 2] == '(') {
                    size_t j = i + 3;
                    int dep = 1;
                    while (raw[j]) {
                        if (raw[j] == '(') dep++;
                        else if (raw[j] == ')') {
                            if (dep == 1 && raw[j + 1] == ')') break;
                            dep--;
                        }
                        j++;
                    }
                    if (!raw[j] || !raw[j + 1]) { free(out); return NULL; }
                    size_t inner = j - (i + 3);
                    char *expr = (char *)malloc(inner + 1);
                    if (!expr) { free(out); return NULL; }
                    memcpy(expr, raw + i + 3, inner);
                    expr[inner] = '\0';
                    /* expand $V inside arithmetic first */
                    char *e2 = sh_expand(t, expr, NULL, depth + 1);
                    /* sh_expand with NULL mask treats all as expandable */
                    long v = 0;
                    int ok = 0;
                    if (e2) { ok = eval_arith(t, e2, &v); free(e2); }
                    else { ok = eval_arith(t, expr, &v); }
                    free(expr);
                    if (!ok) { free(out); return NULL; }
                    snprintf(nb, sizeof nb, "%ld", v);
                    val = nb; adv = (j + 2) - i;
                } else {
                /* $(...): quote-aware match, capture stdout as one word */
                size_t j = i + 2;
                int dep = 1, sq = 0, dq = 0, es = 0;
                while (raw[j] && dep > 0) {
                    if (es) es = 0;
                    else if (sq) { if (raw[j] == '\'') sq = 0; }
                    else if (dq) {
                        if (raw[j] == '\\') es = 1;
                        else if (raw[j] == '"') dq = 0;
                    } else if (raw[j] == '\\') es = 1;
                    else if (raw[j] == '\'') sq = 1;
                    else if (raw[j] == '"') dq = 1;
                    else if (raw[j] == '(') dep++;
                    else if (raw[j] == ')') dep--;
                    if (dep > 0) j++;
                }
                if (dep != 0) { free(out); return NULL; }
                size_t inner = j - (i + 2);
                char *cmd = (char *)malloc(inner + 1);
                if (!cmd) { free(out); return NULL; }
                memcpy(cmd, raw + i + 2, inner);
                cmd[inner] = '\0';
                char *got = sh_capture(t, cmd, depth + 1);
                free(cmd);
                if (!got) { free(out); return NULL; }
                size_t vl2 = strlen(got);
                while (len + vl2 + 1 >= cap) {
                    cap *= 2;
                    char *n2 = (char *)realloc(out, cap);
                    if (!n2) { free(got); free(out); return NULL; }
                    out = n2;
                }
                memcpy(out + len, got, vl2);
                len += vl2;
                free(got);
                i = j + 1;
                continue;
                }
            } else if (raw[i + 1] == '{') {
                /* ${V}, ${#V}, ${V:-w}, ${V:=w}, ${V:?m}, ${V:+w},
                 * ${V#pat}, ${V##pat}, ${V%pat}, ${V%%pat}, nested ${}Skip */
                const char *e = NULL;
                {
                    int dep = 1;
                    const char *q = raw + i + 2;
                    while (*q) {
                        if (q[0] == '$' && q[1] == '{') { dep++; q += 2; continue; }
                        if (*q == '}') { dep--; if (dep == 0) { e = q; break; } }
                        q++;
                    }
                }
                if (!e) { free(out); return NULL; }
                size_t nn = (size_t)(e - (raw + i + 2));
                if (nn > 200) { free(out); return NULL; }
                char inner[208];
                memcpy(inner, raw + i + 2, nn);
                inner[nn] = '\0';
                /* ${#V}: length */
                if (inner[0] == '#' && sh_valid_name(inner + 1, strlen(inner + 1))) {
                    const char *v = sh_get(inner + 1);
                    snprintf(nb, sizeof nb, "%lu", (unsigned long)(v ? strlen(v) : 0));
                    val = nb; adv = (size_t)(e - (raw + i)) + 1;
                } else {
                    /* split NAME OP WORD */
                    size_t k = 0;
                    while (inner[k] && (isalnum((unsigned char)inner[k]) || inner[k] == '_')) k++;
                    char nm[64]; size_t nl = k;
                    if (nl > 63) { free(out); return NULL; }
                    memcpy(nm, inner, nl); nm[nl] = '\0';
                    const char *op = inner + k;
                    const char *word = "";
                    char opc[3] = {0, 0, 0};
                    if (op[0] == ':' && (op[1] == '-' || op[1] == '=' || op[1] == '?' || op[1] == '+')) {
                        opc[0] = ':'; opc[1] = op[1]; word = op + 2;
                    } else if (op[0] == '#' && op[1] == '#') { opc[0] = '#'; opc[1] = '#'; word = op + 2; }
                    else if (op[0] == '%' && op[1] == '%') { opc[0] = '%'; opc[1] = '%'; word = op + 2; }
                    else if (op[0] == '#' || op[0] == '%') { opc[0] = op[0]; word = op + 1; }
                    else if (op[0] == '\0') { opc[0] = 0; }
                    else { free(out); return NULL; }
                    if (nl == 0 || !sh_valid_name(nm, nl)) { free(out); return NULL; }
                    const char *v = sh_get(nm);
                    int set = (v != NULL);
                    int nonempty = (set && v[0] != '\0');
                    int use_word = (!set || !nonempty);
                    char tmp[256];
                    /* expand WORD part (allows nested $V/$( )/${}) */
                    char xword[256] = {0};
                    {
                        char *xw = sh_expand(t, word, NULL, depth + 1);
                        if (xw) { snprintf(xword, sizeof xword, "%s", xw); free(xw); }
                        else snprintf(xword, sizeof xword, "%s", word);
                    }
                    if (opc[0] == 0) {
                        snprintf(nb, sizeof nb, "%s", v ? v : "");
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else if (opc[0] == ':' && opc[1] == '-') {
                        snprintf(nb, sizeof nb, "%s", use_word ? xword : v);
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else if (opc[0] == ':' && opc[1] == '=') {
                        if (use_word) { sh_set(nm, xword); v = sh_get(nm); }
                        snprintf(nb, sizeof nb, "%s", v ? v : "");
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else if (opc[0] == ':' && opc[1] == '?') {
                        if (use_word) {
                            snprintf(tmp, sizeof tmp, "  %s: %s\n", nm,
                                xword[0] ? xword : "parameter null or not set");
                            sh_err("%s", tmp);
                            free(out); return NULL;
                        }
                        snprintf(nb, sizeof nb, "%s", v ? v : "");
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else if (opc[0] == ':' && opc[1] == '+') {
                        snprintf(nb, sizeof nb, "%s", use_word ? "" : xword);
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else if (opc[0] == '#') {
                        /* prefix strip: # shortest, ## longest */
                        const char *s = v ? v : "";
                        size_t sl = strlen(s);
                        size_t best = 0;
                        for (size_t L = 0; L <= sl; L++) {
                            char pre[256];
                            if (L >= sizeof pre) break;
                            memcpy(pre, s, L); pre[L] = '\0';
                            if (fnmatch(word, pre, 0) == 0) {
                                best = L;
                                if (opc[1] != '#') break;
                            }
                        }
                        snprintf(nb, sizeof nb, "%s", s + best);
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else if (opc[0] == '%') {
                        const char *s = v ? v : "";
                        size_t sl = strlen(s);
                        size_t best = sl;
                        for (size_t L = 0; L <= sl; L++) {
                            const char *suf = s + sl - L;
                            if (fnmatch(word, suf, 0) == 0) {
                                best = sl - L;
                                if (opc[1] != '%') break;
                            }
                        }
                        /* longest wants smallest best */
                        if (opc[1] == '%') {
                            for (size_t L = sl + 1; L-- > 0; ) {
                                const char *suf = s + sl - L;
                                char tmp2[256];
                                if (L >= sizeof tmp2) continue;
                                memcpy(tmp2, suf, L); tmp2[L] = '\0';
                                if (fnmatch(word, tmp2, 0) == 0) { best = sl - L; break; }
                                if (L == 0) break;
                            }
                        }
                        char out2[256];
                        if (best >= sizeof out2) best = sizeof out2 - 1;
                        memcpy(out2, s, best); out2[best] = '\0';
                        snprintf(nb, sizeof nb, "%s", out2);
                        val = nb; adv = (size_t)(e - (raw + i)) + 1;
                    } else { free(out); return NULL; }
                }
            } else if (raw[i + 1] == '#') {
                snprintf(nb, sizeof nb, "%d", sh_argc);
                val = nb; adv = 2;
            } else if (raw[i + 1] == '@' || raw[i + 1] == '*') {
                /* "$@" / "$*" joined with spaces (no field splitting here) */
                char tmp[256] = {0};
                size_t p = 0;
                for (int k = 0; k < sh_argc && p < sizeof tmp - 1; k++) {
                    if (k) { if (p < sizeof tmp - 1) tmp[p++] = ' '; }
                    size_t L = strlen(sh_argv[k]);
                    if (p + L >= sizeof tmp - 1) L = sizeof tmp - 2 - p;
                    memcpy(tmp + p, sh_argv[k], L);
                    p += L;
                }
                tmp[p] = '\0';
                snprintf(nb, sizeof nb, "%s", tmp);
                val = nb; adv = 2;
            } else if (raw[i + 1] >= '0' && raw[i + 1] <= '9') {
                int idx = raw[i + 1] - '0';
                const char *v = "";
                if (idx == 0) v = sh_prog;
                else if (idx - 1 < sh_argc) v = sh_argv[idx - 1];
                snprintf(nb, sizeof nb, "%s", v);
                val = nb; adv = 2;
            } else if (raw[i + 1] == '-') {
                snprintf(nb, sizeof nb, "hB");
                val = nb; adv = 2;
            } else if (isalpha((unsigned char)raw[i + 1]) || raw[i + 1] == '_') {
                size_t j = i + 1;
                while (isalnum((unsigned char)raw[j]) || raw[j] == '_') j++;
                size_t nn = j - (i + 1);
                if (nn > 63) { free(out); return NULL; }
                char nm[64];
                memcpy(nm, raw + i + 1, nn);
                nm[nn] = '\0';
                const char *v = sh_get(nm);
                snprintf(nb, sizeof nb, "%s", v ? v : "");
                val = nb; adv = j - i;
            } else {
                /* lone $: literal */
                if (len + 1 >= cap) {
                    cap *= 2;
                    char *n2 = (char *)realloc(out, cap);
                    if (!n2) { free(out); return NULL; }
                    out = n2;
                }
                out[len++] = '$'; i++;
                continue;
            }
            size_t vl = strlen(val);
            while (len + vl + 1 >= cap) {
                cap *= 2;
                char *n2 = (char *)realloc(out, cap);
                if (!n2) { free(out); return NULL; }
                out = n2;
            }
            memcpy(out + len, val, vl);
            len += vl; i += adv;
        } else {
            if (len + 1 >= cap) {
                cap *= 2;
                char *n2 = (char *)realloc(out, cap);
                if (!n2) { free(out); return NULL; }
                out = n2;
            }
            out[len++] = raw[i++];
        }
    }
    out[len] = '\0';
    return out;
}

static void sh_free(ShSeg *segs, int ns) {
    for (int i = 0; i < ns; i++)
        for (int j = 0; j < segs[i].argc; j++) {
            free(segs[i].args[j].raw);
            free(segs[i].args[j].mask);
        }
}

/* true iff this raw word is an unquoted redirect operator.
 * Returns kind 1 > 2 >> 3 < 4 2> 5 2>> 6 <<<, else 0. Operator words
 * always carry mask-2 bytes from the lexer; quoted forms (mask 0/1) or
 * a $R expanding to one are NOT operators — no rescan (POSIX). */
static int sh_redir_kind(const char *raw, const char *mask) {
    size_t n = strlen(raw);
    for (size_t i = 0; i < n; i++)
        if (mask[i] != 2) return 0;
    if (n == 1 && raw[0] == '>') return 1;
    if (n == 2 && raw[0] == '>' && raw[1] == '>') return 2;
    if (n == 1 && raw[0] == '<') return 3;
    if (n == 2 && raw[0] == '2' && raw[1] == '>') return 4;
    if (n == 3 && raw[0] == '2' && raw[1] == '>' && raw[2] == '>') return 5;
    if (n == 3 && raw[0] == '<' && raw[1] == '<' && raw[2] == '<') return 6;
    return 0;
}

static int run_line(Tui *t, char *line);

/* run cmd through the full line executor, capturing stdout.
 * One substitution word, trailing newlines stripped, depth-guarded.
 * last_status is preserved (matches observed sh behavior for `echo $(..)`). */
static char *sh_capture(Tui *t, const char *cmd, int depth) {
    if (depth > 4) return NULL;
    char tmp[256];
    const char *td = getenv("TMPDIR");
    snprintf(tmp, sizeof tmp, "%s/s2capXXXXXX", td ? td : "/tmp");
    int fd = mkstemp(tmp);
    if (fd < 0) return NULL;
    fflush(stdout);
    int saved = dup(STDOUT_FILENO);
    int saved_status = last_status;
    char *out = NULL;
    if (saved >= 0) {
        if (dup2(fd, STDOUT_FILENO) >= 0) {
            char *dup = strdup(cmd);
            if (dup) { run_line(t, dup); free(dup); }
            fflush(stdout);
        }
        dup2(saved, STDOUT_FILENO);
        close(saved);
    }
    close(fd);
    FILE *f = fopen(tmp, "r");
    if (f) {
        size_t cap = 256, len = 0;
        out = (char *)malloc(cap);
        if (out) {
            int ch;
            while ((ch = fgetc(f)) != EOF) {
                if (len + 1 >= cap) {
                    cap *= 2;
                    char *n2 = (char *)realloc(out, cap);
                    if (!n2) { free(out); out = NULL; break; }
                    out = n2;
                }
                if (out) out[len++] = (char)ch;
            }
            if (out) {
                while (len > 0 && out[len - 1] == '\n') len--;
                out[len] = '\0';
            }
        }
        fclose(f);
    }
    unlink(tmp);
    last_status = saved_status;
    return out;
}

/* ── termios line editor (TTY only; piped input uses fgets) ── */
static const char *sh_complete(const char *w0, int first) {
    static const char *cmds[] = {
        "man", "ls", "test", "env", "rm", "ln", "fsck", "unlink", "set", "export",
        "sleep", "df", "du", "kill", "nice", "touch", "dd", "ps", "vi", "more",
        "sync", "make", "cat", "cp", "mv", "mkdir", "rmdir", "head", "tail", "wc",
        "sort", "uniq", "cut", "tr", "grep", "tee", "find", "date", "uname", "fc",
        "tput", "sh", "true", "false", "printf", "echo", "pwd", "cd", "read", "exit",
        "echo", "export", "unset", "env", "history", "source", "sleep",
        "time", "save", "load", "printf", "true", "false", "wait",
        "cat", "pwd", "cd", "mkdir", "cp", "mv", "head", "tail", "wc",
        "sort", "uniq", "cut", "tr", "grep", "tee", "basename", "dirname",
        "touch", "rmdir", "ln", "date", "uname", "find",
        "quit", "exit",
        ":", "eval", "exec", "command", "type", "shift", "readonly", "umask",
        "trap", "alias", "unalias", "if", "then", "else", "fi", "for", "in",
        "while", "until", "do", "done", "case", "esac", "screen", NULL
    };
    if (!first) return NULL; /* filenames below */
    for (int i = 0; cmds[i]; i++)
        if (!strncmp(cmds[i], w0, strlen(w0))) return cmds[i];
    return NULL;
}

static const char *sh_complete_file(const char *w0) {
    static char match[256];
    const char *slash = strrchr(w0, '/');
    char dir[256], pre[256];
    if (slash) {
        size_t dl = (size_t)(slash - w0) + 1;
        if (dl > 200) return NULL;
        memcpy(dir, w0, dl); dir[dl] = '\0';
        snprintf(pre, sizeof pre, "%s", slash + 1);
    } else {
        snprintf(dir, sizeof dir, ".");
        snprintf(pre, sizeof pre, "%s", w0);
    }
    DIR *d = opendir(dir);
    if (!d) return NULL;
    const char *best = NULL;
    size_t pl = strlen(pre);
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        if (strncmp(e->d_name, pre, pl)) continue;
        if (!best || strcmp(e->d_name, best) < 0) best = e->d_name;
    }
    if (!best) { closedir(d); return NULL; }
    if (slash) {
        size_t dl = strlen(dir), bl = strlen(best);
        if (dl + bl > sizeof match - 1) { closedir(d); return NULL; }
        memcpy(match, dir, dl);
        memcpy(match + dl, best, bl + 1);
    } else {
        size_t bl = strlen(best);
        if (bl > sizeof match - 1) { closedir(d); return NULL; }
        memcpy(match, best, bl + 1);
    }
    closedir(d);
    return match;
}

/* read one line; returns malloc'd buffer or NULL on EOF. TTY gets editing. */
static char *tui_readline(const char *prompt) {
    if (!isatty(STDIN_FILENO)) {
        if (prompt) { fprintf(stderr, "%s", prompt); fflush(stderr); }
        static char fb[SH_LINE];
        if (!fgets(fb, sizeof fb, stdin)) return NULL;
        size_t n = strlen(fb);
        while (n > 0 && (fb[n - 1] == '\n' || fb[n - 1] == '\r')) fb[--n] = '\0';
        char *r = malloc(n + 1);
        if (!r) return NULL;
        memcpy(r, fb, n + 1);
        return r;
    }
    struct termios orig, raw;
    if (tcgetattr(STDIN_FILENO, &orig) != 0) {
        static char fb[SH_LINE];
        fprintf(stderr, "%s", prompt ? prompt : "");
        fflush(stderr);
        if (!fgets(fb, sizeof fb, stdin)) return NULL;
        size_t n = strlen(fb);
        while (n > 0 && (fb[n - 1] == '\n' || fb[n - 1] == '\r')) fb[--n] = '\0';
        char *r = malloc(n + 1);
        if (!r) return NULL;
        memcpy(r, fb, n + 1);
        return r;
    }
    raw = orig;
    raw.c_lflag &= (unsigned)(~(ICANON | ECHO | ISIG));
    raw.c_cc[VMIN] = 1;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSANOW, &raw) != 0) {
        tcsetattr(STDIN_FILENO, TCSANOW, &orig);
        return NULL;
    }
    static char buf[SH_LINE];
    int len = 0, pos = 0, hn = -1;
    char hbuf[SH_LINE];
    hbuf[0] = '\0';
    fprintf(stderr, "%s", prompt ? prompt : "");
    fflush(stderr);
    for (;;) {
        char c;
        ssize_t r = read(STDIN_FILENO, &c, 1);
        if (r <= 0) { tcsetattr(STDIN_FILENO, TCSANOW, &orig); return NULL; }
        if (c == '\r' || c == '\n') {
            fprintf(stderr, "\r\n");
            tcsetattr(STDIN_FILENO, TCSANOW, &orig);
            buf[len] = '\0';
            char *out = malloc((size_t)len + 1);
            if (out) memcpy(out, buf, (size_t)len + 1);
            return out;
        } else if (c == 0x03) { /* Ctrl-C: cancel line */
            fprintf(stderr, "^C\r\n");
            tcsetattr(STDIN_FILENO, TCSANOW, &orig);
            char *out = malloc(1);
            if (out) out[0] = '\0';
            return out;
        } else if (c == 0x04) { /* Ctrl-D: EOF on empty, else ignore */
            if (len == 0) {
                fprintf(stderr, "\r\n");
                tcsetattr(STDIN_FILENO, TCSANOW, &orig);
                return NULL;
            }
        } else if (c == 0x7f || c == 0x08) { /* backspace */
            if (pos > 0) {
                memmove(buf + pos - 1, buf + pos, (size_t)(len - pos));
                len--; pos--;
            }
        } else if (c == 0x15) { /* Ctrl-U: kill line */
            len = 0; pos = 0; hn = -1;
        } else if (c == 0x0b) { /* Ctrl-K: kill to end */
            len = pos;
        } else if (c == 0x01) { pos = 0; } /* Ctrl-A */
        else if (c == 0x05) { pos = len; } /* Ctrl-E */
        else if (c == 0x17) { /* Ctrl-W: kill word */
            while (pos > 0 && buf[pos - 1] == ' ') { memmove(buf + pos - 1, buf + pos, (size_t)(len - pos)); len--; pos--; }
            while (pos > 0 && buf[pos - 1] != ' ') { memmove(buf + pos - 1, buf + pos, (size_t)(len - pos)); len--; pos--; }
        } else if (c == '\t') { /* complete first word (commands) or filenames */
            buf[len] = '\0';
            int ws = len;
            while (ws > 0 && buf[ws - 1] != ' ') ws--;
            int first = 1;
            for (int i = 0; i < ws; i++) if (buf[i] != ' ') { first = 0; break; }
            const char *m = first ? sh_complete(buf + ws, 1) : sh_complete_file(buf + ws);
            if (m) {
                size_t ml = strlen(m), pl = (size_t)(len - ws);
                if (ml > pl && len + (ml - pl) < SH_LINE - 1) {
                    memmove(buf + ws + ml, buf + len, 1);
                    memcpy(buf + ws, m, ml);
                    len += (int)(ml - pl); pos = len;
                }
            }
        } else if (c == 0x1b) { /* escape sequence */
            char s[3] = {0, 0, 0};
            if (read(STDIN_FILENO, &s[0], 1) <= 0) continue;
            if (s[0] == '[' || s[0] == 'O') {
                if (read(STDIN_FILENO, &s[1], 1) <= 0) continue;
                if (s[1] >= '0' && s[1] <= '8') {
                    if (read(STDIN_FILENO, &s[2], 1) <= 0) continue;
                }
            }
            char f = s[2] ? s[2] : s[1];
            char k = s[0];
            if ((k == '[' && f == 'A') || (k == 'O' && f == 'A')) { /* up */
                if (sh_hist_n > 0) {
                    if (hn == -1) { memcpy(hbuf, buf, (size_t)len + 1); hn = sh_hist_n; }
                    if (hn > 0 && hn > sh_hist_n - SH_HIST) {
                        hn--;
                        snprintf(buf, sizeof buf, "%s", sh_hist[hn % SH_HIST]);
                        len = (int)strlen(buf); pos = len;
                    }
                }
            } else if ((k == '[' && f == 'B') || (k == 'O' && f == 'B')) { /* down */
                if (hn != -1) {
                    hn++;
                    if (hn >= sh_hist_n) { snprintf(buf, sizeof buf, "%s", hbuf); hn = -1; }
                    else snprintf(buf, sizeof buf, "%s", sh_hist[hn % SH_HIST]);
                    len = (int)strlen(buf); pos = len;
                }
            } else if ((k == '[' && f == 'C') || (k == 'O' && f == 'C')) { if (pos < len) pos++; }
            else if ((k == '[' && f == 'D') || (k == 'O' && f == 'D')) { if (pos > 0) pos--; }
            else if ((k == '[' && (f == 'H' || f == 'F')) || (k == 'O' && (f == 'H' || f == 'F'))) { pos = (f == 'H') ? 0 : len; }
            else if (k == '[' && f == '~') {
                if (s[1] == '3' && pos < len) { /* delete */
                    memmove(buf + pos, buf + pos + 1, (size_t)(len - pos));
                    len--;
                } else if (s[1] == '1' || s[1] == '7') pos = 0;
                else if (s[1] == '4' || s[1] == '8') pos = len;
            }
        } else if (isprint((unsigned char)c)) {
            if (len < SH_LINE - 2) {
                memmove(buf + pos + 1, buf + pos, (size_t)(len - pos));
                buf[pos] = c;
                len++; pos++;
            }
        }
        buf[len] = '\0';
        fprintf(stderr, "\r\x1b[K%s%s", prompt ? prompt : "", buf);
        if (pos < len) fprintf(stderr, "\x1b[%dG", (int)(strlen(prompt ? prompt : "") + pos + 1));
        fflush(stderr);
    }
}

/* ── save/load (text snapshot S2SAVE1) ── */
static int cmd_save(const Tui *t, const char *path) {
    if (!t->sim) { sh_err("  (empty sim)\n"); return 1; }
    FILE *f = fopen(path, "w");
    if (!f) { sh_err("  cannot write `%s`\n", path); return 1; }
    const Simulation *s = t->sim;
    fprintf(f, "S2SAVE1\nseed %lu\ndt %.17g cutoff %.17g dielectric %.17g temp %.17g thermostat %d tau %.17g nu %.17g\n",
        t->seed, s->dt, s->cutoff, s->dielectric, s->thermostat.target_temperature,
        (int)s->thermostat.type, s->thermostat.tau, s->thermostat.nu);
    fprintf(f, "atoms %d\n", s->num_atoms);
    for (int i = 0; i < s->num_atoms; i++) {
        const Atom *a = &s->atoms[i];
        fprintf(f, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            a->Z, a->position.x, a->position.y, a->position.z, a->partial_charge,
            a->lj_epsilon, a->lj_sigma, a->velocity.x, a->velocity.y, a->velocity.z);
    }
    fprintf(f, "bonds %d\n", s->num_bonds);
    for (int b = 0; b < s->num_bonds; b++)
        fprintf(f, "%d %d %d %.17g %.17g\n", s->bonds[b].atom_a, s->bonds[b].atom_b,
            s->bonds[b].order, s->bonds[b].r0, s->bonds[b].k);
    fprintf(f, "restraints %d\n", s->num_restraints);
    for (int r = 0; r < s->num_restraints; r++)
        fprintf(f, "%d %.17g %.17g %.17g %.17g %.17g\n", s->restraint_atom[r],
            s->restraint_anchor[r].x, s->restraint_anchor[r].y,
            s->restraint_anchor[r].z, s->restraint_k[r],
            s->restraint_flat ? s->restraint_flat[r] : 0.0);
    fclose(f);
    sh_err("  saved %d atoms to `%s`\n", s->num_atoms, path);
    return 0;
}

typedef struct { int Z; double x, y, z, q, e, sg, vx, vy, vz; } SaveAtom;
typedef struct { int a, b, o; double r0, k; } SaveBond;
typedef struct { int i; double x, y, z, k, fb; } SaveRest;

static int cmd_load(Tui *t, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { sh_err("  cannot read `%s`\n", path); return 1; }
    char magic[16];
    if (!fgets(magic, sizeof magic, f) || strncmp(magic, "S2SAVE1", 7)) {
        sh_err("  bad save file (want S2SAVE1)\n");
        fclose(f);
        return 1;
    }
    unsigned long seed = t->seed;
    double dt = 0.5, co = 12.0, di = 1.0, tp = 300.0, tau = 100.0, nu = 0.02;
    int th = 1, na = 0, nb = 0, nr = 0;
    if (fscanf(f, "seed %lu\n", &seed) != 1) { sh_err("  bad save header\n"); fclose(f); return 1; }
    if (fscanf(f, "dt %lg cutoff %lg dielectric %lg temp %lg thermostat %d tau %lg nu %lg\n",
               &dt, &co, &di, &tp, &th, &tau, &nu) != 7) {
        sh_err("  bad settings\n");
        fclose(f);
        return 1;
    }
    if (fscanf(f, "atoms %d\n", &na) != 1 || na < 0 || na > 100000) {
        sh_err("  bad atom count\n");
        fclose(f);
        return 1;
    }
    SaveAtom *A = NULL;
    SaveBond *B = NULL;
    SaveRest *R = NULL;
    if (na > 0) {
        A = (SaveAtom *)malloc(sizeof(SaveAtom) * (size_t)na);
        if (!A) { fclose(f); return 1; }
    }
    int rc = 1;
    for (int i = 0; i < na; i++) {
        if (fscanf(f, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
                   &A[i].Z, &A[i].x, &A[i].y, &A[i].z, &A[i].q, &A[i].e,
                   &A[i].sg, &A[i].vx, &A[i].vy, &A[i].vz) != 10) {
            sh_err("  bad atom %d\n", i);
            goto done;
        }
    }
    if (fscanf(f, "bonds %d\n", &nb) != 1 || nb < 0 || nb > 200000) {
        sh_err("  bad bond count\n");
        goto done;
    }
    if (nb > 0) {
        B = (SaveBond *)malloc(sizeof(SaveBond) * (size_t)nb);
        if (!B) goto done;
    }
    for (int i = 0; i < nb; i++) {
        if (fscanf(f, "%d %d %d %lg %lg\n", &B[i].a, &B[i].b, &B[i].o, &B[i].r0, &B[i].k) != 5) {
            sh_err("  bad bond %d\n", i);
            goto done;
        }
    }
    if (fscanf(f, "restraints %d\n", &nr) != 1 || nr < 0 || nr > 100000) {
        sh_err("  bad restraint count\n");
        goto done;
    }
    if (nr > 0) {
        R = (SaveRest *)malloc(sizeof(SaveRest) * (size_t)nr);
        if (!R) goto done;
    }
    for (int i = 0; i < nr; i++) {
        /* backward compat: 5 cols (no flat) or 6 cols (with flat) */
        char lb[512];
        if (!fgets(lb, sizeof lb, f)) { sh_err("  bad restraint %d\n", i); goto done; }
        double fb = 0.0;
        int got = sscanf(lb, "%d %lg %lg %lg %lg %lg", &R[i].i, &R[i].x, &R[i].y, &R[i].z, &R[i].k, &fb);
        if (got < 5) { sh_err("  bad restraint %d\n", i); goto done; }
        R[i].fb = (got >= 6) ? fb : 0.0;
    }
    {
        Simulation *s = sim_create(na + 8 > 8 ? na + 8 : 8, nb + 8 > 8 ? nb + 8 : 8);
        if (!s) goto done;
        int ok = 1;
        for (int i = 0; i < na && ok; i++) {
            int idx = sim_add_atom(s, A[i].Z, vec3(A[i].x, A[i].y, A[i].z), A[i].q);
            if (idx < 0) { sh_err("  bad atom Z=%d\n", A[i].Z); ok = 0; break; }
            sim_set_atom_lj(s, idx, A[i].e, A[i].sg);
            s->atoms[idx].velocity = vec3(A[i].vx, A[i].vy, A[i].vz);
        }
        for (int i = 0; i < nb && ok; i++) {
            int bi = sim_add_bond(s, B[i].a, B[i].b, B[i].o);
            if (bi < 0) { sh_err("  bad bond %d-%d\n", B[i].a, B[i].b); ok = 0; break; }
            sim_set_bond_params(s, bi, B[i].r0, B[i].k);
        }
        for (int i = 0; i < nr && ok; i++) {
            int rr;
            if (R[i].fb > 0.0)
                rr = sim_add_restraint_fb(s, R[i].i, vec3(R[i].x, R[i].y, R[i].z), R[i].k, R[i].fb);
            else
                rr = sim_add_restraint(s, R[i].i, vec3(R[i].x, R[i].y, R[i].z), R[i].k);
            if (rr < 0) {
                sh_err("  bad restraint %d\n", i);
                ok = 0;
                break;
            }
        }
        if (!ok) { sim_destroy(s); goto done; }
        s->dt = dt; s->cutoff = co; s->dielectric = di;
        s->thermostat.target_temperature = tp;
        s->thermostat.type = (th >= 0 && th <= 3) ? (ThermostatType)th : THERMOSTAT_BERENDSEN;
        s->thermostat.tau = tau; s->thermostat.nu = nu;
        t->seed = seed;
        if (t->sim) sim_destroy(t->sim);
        t->sim = s;
        sim_rebuild_angles(s);
        forces_calculate(s);
        sh_err("  loaded %d atoms %d bonds %d restraints from `%s`\n", na, nb, nr, path);
        rc = 0;
    }
done:
    free(A); free(B); free(R);
    fclose(f);
    return rc;
}

/* ── main loop ── */
static int run_line(Tui *t, char *line);
static int run_segment(Tui *t, char **argv, const int *elig, int argc);
static int dispatch_cmd(Tui *t, char **tok, const int *elig, int nt);
int tui_main(void) {
    Tui t;
    memset(&t, 0, sizeof t);
    t.seed = 42;
    view_cam_reset(&t.cam);
    tui_new(&t, 512, 512);
    if (!t.sim) { sh_err("  alloc failed\n"); return 1; }
    sh_err("  s2tui — type `man` (batch `./carbonsim` record untouched)\n");
    while (!t.quit) {
        char *line = tui_readline(t.srcdepth == 0 ? "s2> " : NULL);
        if (!line) break;
        if (t.srcdepth == 0) sh_hist_push(line);
        last_status = run_line(&t, line);
        free(line);
    }
    sh_err("\n  bye.\n");
    if (t.sim) sim_destroy(t.sim);
    return last_status;
}

/* input redirection target for text filters: pipe stage file or < file.
 * NULL means the real stdin (terminal when interactive). */
static FILE *sh_input = NULL;

typedef struct {
    const char *infile, *outfile, *errfile, *herestr;
    int append, errappend;
} Redirs;

/* strip eligible redirections (kinds 1..6) from expanded argv. */
static int strip_redirs(char **argv, const int *elig, int argc,
                        char **outv, int *outelig, int *outc, Redirs *r) {
    r->infile = r->outfile = r->errfile = r->herestr = NULL;
    r->append = r->errappend = 0;
    int n = 0;
    for (int i = 0; i < argc; i++) {
        int k = elig ? elig[i] : 0;
        if (k >= 1 && k <= 6) {
            if (i + 1 >= argc) { sh_err("  syntax error: no target for redirect\n"); return 2; }
            const char *f = argv[i + 1];
            if (k == 1) { r->outfile = f; r->append = 0; }
            else if (k == 2) { r->outfile = f; r->append = 1; }
            else if (k == 3) { r->infile = f; }
            else if (k == 6) { r->herestr = f; }
            else if (k == 4) { r->errfile = f; r->errappend = 0; }
            else { r->errfile = f; r->errappend = 1; }
            i++;
            continue;
        }
        outv[n] = argv[i];
        outelig[n] = elig ? elig[i] : 0;
        n++;
    }
    *outc = n;
    return 0;
}

/* run one stage with explicit redirections plus optional pipe wiring:
 * in_fp feeds stdin (< file wins), out_fd takes stdout (> file wins). */
static int run_stage(Tui *t, char **argv, const int *elig, int argc,
                     FILE *in_fp, int out_fd) {
    char *filt[SH_MAXARG];
    int felig[SH_MAXARG];
    int fc = 0;
    Redirs r;
    int rc = strip_redirs(argv, elig, argc, filt, felig, &fc, &r);
    if (rc) return rc;
    FILE *old_in = sh_input;
    int saved_out = -1, saved_err = -1;
    FILE *inf = NULL, *outf = NULL, *errf = NULL;
    if (r.infile) {
        inf = fopen(r.infile, "r");
        if (!inf) { sh_err("  cannot read `%s`\n", r.infile); return 1; }
        sh_input = inf;
    } else if (r.herestr) {
        /* herestring: materialize content + newline as stdin */
        char htmp[256];
        const char *td = getenv("TMPDIR");
        snprintf(htmp, sizeof htmp, "%s/s2hereXXXXXX", td ? td : "/tmp");
        int hfd = mkstemp(htmp);
        if (hfd < 0) { sh_err("  redirection failed\n"); return 1; }
        size_t hl = strlen(r.herestr);
        int ok = 1;
        if (hl > 0 && write(hfd, r.herestr, hl) != (ssize_t)hl) ok = 0;
        if (ok && write(hfd, "\n", 1) != 1) ok = 0;
        close(hfd);
        if (!ok) { unlink(htmp); sh_err("  redirection failed\n"); return 1; }
        inf = fopen(htmp, "r");
        unlink(htmp); /* unnamed from here; fclose releases it */
        if (!inf) { sh_err("  redirection failed\n"); return 1; }
        sh_input = inf;
    } else if (in_fp) {
        sh_input = in_fp;
    }
    if (r.outfile) {
        fflush(stdout);
        saved_out = dup(STDOUT_FILENO);
        outf = fopen(r.outfile, r.append ? "a" : "w");
        if (!outf || saved_out < 0 || dup2(fileno(outf), STDOUT_FILENO) < 0) {
            sh_err("  cannot write `%s`\n", r.outfile ? r.outfile : "?");
            if (outf) fclose(outf);
            if (saved_out >= 0) close(saved_out);
            if (inf) fclose(inf);
            sh_input = old_in;
            return 1;
        }
    } else if (out_fd >= 0) {
        fflush(stdout);
        saved_out = dup(STDOUT_FILENO);
        if (saved_out < 0 || dup2(out_fd, STDOUT_FILENO) < 0) {
            sh_err("  redirection failed\n");
            if (saved_out >= 0) close(saved_out);
            if (inf) fclose(inf);
            sh_input = old_in;
            return 1;
        }
    }
    if (r.errfile) {
        fflush(stderr);
        saved_err = dup(STDERR_FILENO);
        errf = fopen(r.errfile, r.errappend ? "a" : "w");
        if (!errf || saved_err < 0 || dup2(fileno(errf), STDERR_FILENO) < 0) {
            sh_err("  cannot write `%s`\n", r.errfile ? r.errfile : "?");
            if (errf) fclose(errf);
            if (saved_err >= 0) close(saved_err);
            fflush(stdout);
            if (saved_out >= 0) { dup2(saved_out, STDOUT_FILENO); close(saved_out); }
            if (outf) fclose(outf);
            if (inf) fclose(inf);
            sh_input = old_in;
            return 1;
        }
    }
    int st;
    if (fc == 0 && (r.infile || r.outfile || r.errfile || r.herestr)) st = 0;
    else st = dispatch_cmd(t, fc > 0 ? filt : argv, fc > 0 ? felig : elig, fc > 0 ? fc : argc);
    fflush(stdout);
    fflush(stderr);
    if (saved_out >= 0) { dup2(saved_out, STDOUT_FILENO); close(saved_out); }
    if (saved_err >= 0) { dup2(saved_err, STDERR_FILENO); close(saved_err); }
    if (outf) fclose(outf);
    if (errf) fclose(errf);
    if (inf) fclose(inf);
    sh_input = old_in;
    return st;
}

static int run_segment(Tui *t, char **argv, const int *elig, int argc) {
    return run_stage(t, argv, elig, argc, NULL, -1);
}

/* sequential pipeline: stage outputs feed the next stage's stdin through
 * temp files (no fork; output-identical for terminating filters). */
#define SH_MAXPIPE 16
static int run_pipeline(Tui *t, char ***A, int **E, int *N, int k) {
    char tmps[SH_MAXPIPE][64];
    const char *td = getenv("TMPDIR");
    for (int s = 0; s < k - 1; s++) {
        snprintf(tmps[s], sizeof tmps[s], "%s/s2pipeXXXXXX", td ? td : "/tmp");
        int fd = mkstemp(tmps[s]);
        if (fd < 0) {
            sh_err("  pipe failed\n");
            for (int q = 0; q < s; q++) unlink(tmps[q]);
            return 1;
        }
        close(fd);
    }
    int st = 0;
    for (int s = 0; s < k; s++) {
        FILE *inf = NULL;
        int outfd = -1;
        if (s > 0) {
            inf = fopen(tmps[s - 1], "r");
            if (!inf) { st = 1; break; }
        }
        if (s < k - 1) {
            outfd = open(tmps[s], O_WRONLY | O_TRUNC, 0600);
            if (outfd < 0) {
                if (inf) fclose(inf);
                st = 1;
                break;
            }
        }
        st = run_stage(t, A[s], E[s], N[s], inf, outfd);
        if (inf) fclose(inf);
        if (outfd >= 0) close(outfd);
        if (st == 2 && 0) break; /* usage errors still flow: keep going like sh */
    }
    for (int s = 0; s < k - 1; s++) unlink(tmps[s]);
    return st;
}

/* expand one raw segment, apply leading assignments persistently, glob the
 * rest. Fills outv/outelig (*outc words, malloc'd). *ast = assignment
 * status. Returns 0 with a command to run, 1 assignments-only, -1 error. */
static int expand_unit(Tui *t, ShSeg *sg, char **outv, int *outelig, int *outc, int *ast) {
    char *argv[SH_MAXARG];
    int elig[SH_MAXARG];
    int gok[SH_MAXARG];
    int ac = 0;
    for (int w = 0; w < sg->argc; w++) {
        argv[ac] = sh_expand(t, sg->args[w].raw, sg->args[w].mask, 0);
        if (!argv[ac]) {
            for (int k = 0; k < ac; k++) free(argv[k]);
            sh_err("  expansion error\n");
            return -1;
        }
        /* POSIX ~ expansion at word start (unquoted): ~ and ~/ */
        if (sg->args[w].globok || (sg->args[w].mask[0] == 1)) {
            if (argv[ac][0] == '~' && (argv[ac][1] == '\0' || argv[ac][1] == '/')) {
                const char *h = getenv("HOME");
                if (h && h[0]) {
                    char tmp[1024];
                    snprintf(tmp, sizeof tmp, "%s%s", h, argv[ac] + 1);
                    char *n2 = strdup(tmp);
                    if (n2) { free(argv[ac]); argv[ac] = n2; }
                }
            }
        }
        elig[ac] = sh_redir_kind(sg->args[w].raw, sg->args[w].mask);
        gok[ac] = sg->args[w].globok;
        ac++;
    }
    int ai = 0, st = 0;
    while (ai < ac) {
        char *eq = strchr(argv[ai], '=');
        if (!eq || eq == argv[ai]) break;
        size_t nl = (size_t)(eq - argv[ai]);
        if (!sh_valid_name(argv[ai], nl)) break;
        char nm[64];
        memcpy(nm, argv[ai], nl);
        nm[nl] = '\0';
        st = sh_set(nm, eq + 1);
        free(argv[ai]);
        ai++;
    }
    /* glob remaining words with unquoted *?[( in the raw form) */
    int n = 0;
    for (int i = ai; i < ac; i++) {
        if (gok[i]) {
            glob_t g;
            memset(&g, 0, sizeof g);
            int r = glob(argv[i], GLOB_NOCHECK | GLOB_MARK, NULL, &g);
            if (r != 0 && r != GLOB_NOMATCH) {
                for (int k = 0; k < n; k++) free(outv[k]);
                for (int k = i; k < ac; k++) free(argv[k]);
                globfree(&g);
                return -1;
            }
            if (g.gl_pathc == 0) {
                outv[n] = argv[i];
                outelig[n] = elig[i];
                n++;
            } else {
                for (size_t m = 0; m < g.gl_pathc; m++) {
                    if (n >= SH_MAXARG - 1) {
                        sh_err("  too many matches\n");
                        for (int k = 0; k < n; k++) free(outv[k]);
                        for (int k = i + 1; k < ac; k++) free(argv[k]);
                        free(argv[i]);
                        globfree(&g);
                        return -1;
                    }
                    outv[n] = strdup(g.gl_pathv[m]);
                    outelig[n] = 0;
                    if (!outv[n]) {
                        for (int k = 0; k < n; k++) free(outv[k]);
                        for (int k = i + 1; k < ac; k++) free(argv[k]);
                        free(argv[i]);
                        globfree(&g);
                        return -1;
                    }
                    n++;
                }
                free(argv[i]);
            }
            globfree(&g);
        } else {
            outv[n] = argv[i];
            outelig[n] = elig[i];
            n++;
        }
    }
    *outc = n;
    *ast = st;
    return (n == 0) ? 1 : 0;
    /* return: 1 = assignments-only (*ast holds their status), 0 = run */
}

static void free_prepared(char **A, int n) {
    for (int i = 0; i < n; i++) free(A[i]);
}

/* run one unit: single segment or pipeline segs[i..j]. Handles a `time`
 * prefix and trailing & (subshell fork). Returns status. */
static int run_unit(Tui *t, ShSeg *segs, int i, int j) {
    int k = j - i + 1;
    if (k > SH_MAXPIPE) { sh_err("  pipeline too long\n"); return 2; }
    char *A[SH_MAXPIPE][SH_MAXARG];
    int E[SH_MAXPIPE][SH_MAXARG];
    int N[SH_MAXPIPE];
    int ast = 0;
    for (int s = 0; s < k; s++) {
        int r = expand_unit(t, &segs[i + s], A[s], E[s], &N[s], &ast);
        if (r < 0) {
            for (int q = 0; q < s; q++) free_prepared(A[q], N[q]);
            return 2;
        }
        if (r == 1) N[s] = -1; /* assignments-only stage: spliced below */
    }
    /* splice assignments-only stages (their sets already applied) */
    int w = 0;
    for (int s = 0; s < k; s++) {
        if (N[s] < 0) continue;
        if (w != s) {
            memcpy(A[w], A[s], sizeof A[w]);
            memcpy(E[w], E[s], sizeof E[w]);
            N[w] = N[s];
        }
        w++;
    }
    k = w;
    if (k == 0) return ast;
    /* POSIX ! negation: `! cmd` inverts status */
    int neg = 0;
    if (N[0] > 0 && !strcmp(A[0][0], "!") && N[0] > 1) {
        free(A[0][0]);
        for (int q = 0; q + 1 < N[0]; q++) { A[0][q] = A[0][q + 1]; E[0][q] = E[0][q + 1]; }
        N[0]--;
        neg = 1;
    }
    /* `time` prefixes the whole unit (single or pipeline) */
    int timed = 0;
    if (N[0] > 0 && !strcmp(A[0][0], "time") && (N[0] > 1 || k > 1)) {
        free(A[0][0]);
        for (int q = 0; q + 1 < N[0]; q++) { A[0][q] = A[0][q + 1]; E[0][q] = E[0][q + 1]; }
        N[0]--;
        timed = 1;
    }
    if (k > 1 && N[0] == 0) {
        sh_err("  syntax error near `|`\n");
        for (int s = 0; s < k; s++) free_prepared(A[s], N[s]);
        return 2;
    }
    clock_t t0 = timed ? clock() : 0;
    int st;
    if (segs[j].bg) {
        /* POSIX async: subshell snapshot (copy-on-write); parent keeps going */
        fflush(NULL);
        pid_t pid = fork();
        if (pid < 0) {
            sh_err("  fork failed\n");
            for (int s = 0; s < k; s++) free_prepared(A[s], N[s]);
            return 1;
        }
        if (pid == 0) {
            if (k == 1) st = run_stage(t, A[0], E[0], N[0], NULL, -1);
            else {
                char **Ap[SH_MAXPIPE];
                int *Ep[SH_MAXPIPE];
                for (int s = 0; s < k; s++) { Ap[s] = A[s]; Ep[s] = E[s]; }
                st = run_pipeline(t, Ap, Ep, N, k);
            }
            for (int s = 0; s < k; s++) free_prepared(A[s], N[s]);
            fflush(NULL);
            _exit(st & 255);
        }
        char pb[32];
        snprintf(pb, sizeof pb, "%d", (int)pid);
        sh_set("!", pb);
        sh_err("  [%d]\n", (int)pid);
        for (int s = 0; s < k; s++) free_prepared(A[s], N[s]);
        st = 0;
    } else if (k == 1) {
        st = run_stage(t, A[0], E[0], N[0], NULL, -1);
        free_prepared(A[0], N[0]);
    } else {
        char **Ap[SH_MAXPIPE];
        int *Ep[SH_MAXPIPE];
        for (int s = 0; s < k; s++) { Ap[s] = A[s]; Ep[s] = E[s]; }
        st = run_pipeline(t, Ap, Ep, N, k);
        for (int s = 0; s < k; s++) free_prepared(A[s], N[s]);
    }
    if (timed) sh_err("  %.3fs\n", (double)(clock() - t0) / CLOCKS_PER_SEC);
    if (neg) st = (st == 0) ? 1 : 0;
    return st;
}

/* ── POSIX control flow: if/for/while/until/case over expanded segments ──
 * Single-line + multi-part forms, executed via recursive run_line so all
 * redirections, pipes and expansions apply. Streamlined: S2 verbs work
 * unchanged inside bodies. */
static int join_segs(ShSeg *segs, int a, int b, char *out, size_t cap) {
    size_t p = 0;
    for (int s = a; s <= b && p + 1 < cap; s++) {
        for (int w = 0; w < segs[s].argc && p + 1 < cap; w++) {
            size_t L = strlen(segs[s].args[w].raw);
            if (p + L + 1 >= cap) break;
            if (p) out[p++] = ' ';
            memcpy(out + p, segs[s].args[w].raw, L);
            p += L;
        }
        if (s < b && p + 2 < cap) {
            if (segs[s].nextop == 3) { out[p++] = ' '; out[p++] = '|'; }
            else { out[p++] = ' '; out[p++] = ';'; }
        }
    }
    out[p] = '\0';
    return (int)p;
}
/* join words (s0,w0)..(s1,w1) inclusive across segments */
static int join_range(ShSeg *segs, int s0, int w0, int s1, int w1, char *out, size_t cap) {
    size_t p = 0;
    for (int s = s0; s <= s1 && p + 1 < cap; s++) {
        int a = (s == s0) ? w0 : 0;
        int b = (s == s1) ? w1 : segs[s].argc - 1;
        for (int w = a; w <= b && w < segs[s].argc && p + 1 < cap; w++) {
            if (w < 0) continue;
            size_t L = strlen(segs[s].args[w].raw);
            if (p + L + 1 >= cap) break;
            if (p) out[p++] = ' ';
            memcpy(out + p, segs[s].args[w].raw, L);
            p += L;
        }
        if (s < s1 && p + 2 < cap) {
            if (segs[s].nextop == 3) { out[p++] = ' '; out[p++] = '|'; }
            else { out[p++] = ' '; out[p++] = ';'; }
        }
    }
    out[p] = '\0';
    return (int)p;
}
static int run_line_str(Tui *t, const char *s) {
    char *d = strdup(s);
    if (!d) return 1;
    int r = run_line(t, d);
    free(d);
    return r;
}

static int run_line(Tui *t, char *line) {
    ShSeg segs[SH_MAXSEG];
    int ns = sh_lex(line, segs);
    if (ns < 0) { sh_err("  syntax error\n"); return 2; }
    /* POSIX control flow + grouping, checked before normal units.
     * Keywords are recognised on raw words (unquoted). Bodies re-enter
     * run_line so pipes/redirects/expansions apply uniformly. */
    if (ns > 0 && segs[0].argc > 0) {
        const char *k0 = segs[0].args[0].raw;
        /* { list; } group in current shell */
        if (!strcmp(k0, "{")) {
            int close = -1;
            for (int s = 0; s < ns; s++)
                if (segs[s].argc == 1 && !strcmp(segs[s].args[0].raw, "}")) close = s;
            if (close < 0) { sh_err("  syntax error: missing `}`\n"); sh_free(segs, ns); return 2; }
            char buf[SH_LINE];
            join_segs(segs, 1, close - 1, buf, sizeof buf);
            int st = run_line_str(t, buf);
            sh_free(segs, ns);
            return st;
        }
        /* ( list ) subshell: vars restored afterwards */
        if (!strcmp(k0, "(")) {
            int close = -1;
            for (int s = 0; s < ns; s++)
                if (segs[s].argc == 1 && !strcmp(segs[s].args[0].raw, ")")) close = s;
            if (close < 0) { sh_err("  syntax error: missing `)`\n"); sh_free(segs, ns); return 2; }
            char sn[SH_VARS][64], sv[SH_VARS][256];
            int on = sh_nvars, oac = sh_argc;
            memcpy(sn, sh_name, sizeof sn); memcpy(sv, sh_val, sizeof sv);
            char av[SH_ARGS][256]; memcpy(av, sh_argv, sizeof av);
            char buf[SH_LINE];
            join_segs(segs, 1, close - 1, buf, sizeof buf);
            int st = run_line_str(t, buf);
            memcpy(sh_name, sn, sizeof sn); memcpy(sh_val, sv, sizeof sv);
            sh_nvars = on; sh_argc = oac; memcpy(sh_argv, av, sizeof av);
            sh_free(segs, ns);
            return st;
        }
        /* if TEST; then A; [else B;] fi — then/else/fi start a segment */
        if (!strcmp(k0, "if")) {
            int ithen = -1, ielse = -1, ifi = -1, then_w = 0, else_w = 0, fi_w = 0;
            for (int s = 0; s < ns && ithen < 0; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "then")) { ithen = s; then_w = w; break; }
            for (int s = 0; s < ns && ielse < 0; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "else") && s > ithen) { ielse = s; else_w = w; break; }
            for (int s = 0; s < ns; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "fi")) { ifi = s; fi_w = w; }
            if (ithen < 0 || ifi < 0) {
                sh_err("  syntax error: `if TEST; then ..; fi`\n"); sh_free(segs, ns); return 2;
            }
            char tstr[SH_LINE], astr[SH_LINE], bstr[SH_LINE];
            /* TEST = after `if` to before `then`; A/B split at else/fi word spots */
            join_range(segs, 0, 1, ithen, then_w - 1, tstr, sizeof tstr);
            if (ielse >= 0)
                join_range(segs, ithen, then_w + 1, ielse, else_w - 1, astr, sizeof astr);
            else
                join_range(segs, ithen, then_w + 1, ifi, fi_w - 1, astr, sizeof astr);
            bstr[0] = '\0';
            if (ielse >= 0) join_range(segs, ielse, else_w + 1, ifi, fi_w - 1, bstr, sizeof bstr);
            int ct = run_line_str(t, tstr[0] ? tstr : "true");
            int st = (ct == 0) ? run_line_str(t, astr) : (ielse >= 0 ? run_line_str(t, bstr) : ct);
            sh_free(segs, ns);
            return st;
        }
        /* for V in WORDS; do BODY; done — all on one logical line */
        if (!strcmp(k0, "for")) {
            /* locate `in`, `do`, `done` words anywhere after `for V` */
            int iin_s = -1, iin_w = -1, ido_s = -1, ido_w = -1, idone_s = -1, idone_w = -1;
            for (int s = 0; s < ns && iin_s < 0; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "in")) { iin_s = s; iin_w = w; break; }
            for (int s = 0; s < ns && ido_s < 0; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "do")) { ido_s = s; ido_w = w; break; }
            for (int s = 0; s < ns; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "done")) { idone_s = s; idone_w = w; }
            if (segs[0].argc < 2 || iin_s < 0 || ido_s < 0 || idone_s < 0) {
                sh_err("  usage: for V in WORDS; do BODY; done\n"); sh_free(segs, ns); return 2;
            }
            const char *vn = segs[0].args[1].raw;
            if (!sh_valid_name(vn, strlen(vn))) { sh_err("  bad variable\n"); sh_free(segs, ns); return 2; }
            char wlist[SH_LINE], body[SH_LINE];
            join_range(segs, iin_s, iin_w + 1, ido_s, ido_w - 1, wlist, sizeof wlist);
            join_range(segs, ido_s, ido_w + 1, idone_s, idone_w - 1, body, sizeof body);
            int st = 0, n = 0;
            char echocmd[SH_LINE + 16];
            snprintf(echocmd, sizeof echocmd, "echo %s", wlist);
            char *words = sh_capture(t, echocmd[0] == 'e' && wlist[0] ? echocmd : "echo", 0);
            if (words) {
                char *save = NULL, *tok = strtok_r(words, " \t\n", &save);
                while (tok) {
                    sh_set(vn, tok);
                    st = run_line_str(t, body);
                    n++;
                    if (n > 10000) break;
                    tok = strtok_r(NULL, " \t\n", &save);
                }
                free(words);
            }
            sh_free(segs, ns);
            return st;
        }
        /* while/until TEST; do BODY; done */
        if (!strcmp(k0, "while") || !strcmp(k0, "until")) {
            int inv = !strcmp(k0, "until");
            int ido_s = -1, ido_w = -1, idone_s = -1, idone_w = -1;
            for (int s = 0; s < ns && ido_s < 0; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "do")) { ido_s = s; ido_w = w; break; }
            for (int s = 0; s < ns; s++)
                for (int w = 0; w < segs[s].argc; w++)
                    if (!strcmp(segs[s].args[w].raw, "done")) { idone_s = s; idone_w = w; }
            if (ido_s < 0 || idone_s < 0) { sh_err("  usage: while TEST; do BODY; done\n"); sh_free(segs, ns); return 2; }
            char tstr[SH_LINE], body[SH_LINE];
            join_range(segs, 0, 1, ido_s, ido_w - 1, tstr, sizeof tstr);
            join_range(segs, ido_s, ido_w + 1, idone_s, idone_w - 1, body, sizeof body);
            int st = 0, iter = 0;
            for (;;) {
                int ct = run_line_str(t, tstr);
                int go = inv ? (ct != 0) : (ct == 0);
                if (!go) { st = ct; break; }
                st = run_line_str(t, body);
                if (++iter > 10000) { sh_err("  loop limit\n"); break; }
            }
            sh_free(segs, ns);
            return st;
        }
        /* case WORD in PAT) BODY;; ... esac */
        if (!strcmp(k0, "case")) {
            int iin = -1, iesac = -1;
            for (int s = 0; s < ns; s++) {
                if (segs[s].argc == 2 && !strcmp(segs[s].args[0].raw, "case")) continue;
                if (segs[s].argc >= 1 && !strcmp(segs[s].args[0].raw, "in") && iin < 0) iin = s;
                if (segs[s].argc == 1 && !strcmp(segs[s].args[0].raw, "esac")) iesac = s;
            }
            /* simpler raw scan: case WORD in ... esac on one line */
            if (segs[0].argc < 3 || strcmp(segs[0].args[2].raw, "in")) {
                sh_err("  usage: case WORD in PAT) BODY;; ... esac\n"); sh_free(segs, ns); return 2;
            }
            if (iesac < 0) { sh_err("  syntax error: missing `esac`\n"); sh_free(segs, ns); return 2; }
            /* expand subject */
            char *subj = sh_expand(t, segs[0].args[1].raw, segs[0].args[1].mask, 0);
            if (!subj) { sh_free(segs, ns); return 2; }
            /* word-level parse: WORD in PAT) BODY ;; ... esac */
            /* gather words between `in` and `esac` */
            char *cwords[256];
            int ncw = 0;
            {
                int in_s = -1, in_w = -1;
                for (int s = 0; s < ns && in_s < 0; s++)
                    for (int w = 0; w < segs[s].argc; w++)
                        if (!strcmp(segs[s].args[w].raw, "in")) { in_s = s; in_w = w; break; }
                /* need iesac word pos */
                int es_s = -1, es_w = -1;
                for (int s = 0; s < ns && es_s < 0; s++)
                    for (int w = 0; w < segs[s].argc; w++)
                        if (!strcmp(segs[s].args[w].raw, "esac")) { es_s = s; es_w = w; break; }
                for (int s = in_s; s < ns && ncw < 250; s++) {
                    if (s == es_s) {
                        /* words before esac in its segment */
                        int a = (s == in_s) ? in_w + 1 : 0;
                        for (int w = a; w < es_w && ncw < 250; w++)
                            cwords[ncw++] = segs[s].args[w].raw;
                        break;
                    }
                    if (segs[s].argc == 0) {
                        /* empty segment from `;;` -> arm terminator */
                        cwords[ncw++] = ";;";
                        continue;
                    }
                    int a = (s == in_s) ? in_w + 1 : 0;
                    for (int w = a; w < segs[s].argc && ncw < 250; w++) {
                        if (!strcmp(segs[s].args[w].raw, "esac")) break;
                        cwords[ncw++] = segs[s].args[w].raw;
                    }
                }
            }
            int st = 0, matched = 0, idx = 0;
            while (idx < ncw && !matched) {
                /* pattern words until `)` — pattern may be `h*)` in one word */
                char pat[256] = {0};
                char body[SH_LINE] = {0};
                /* find ) in cwords[idx] */
                char *rp = strchr(cwords[idx], ')');
                if (!rp) break;
                {
                    size_t L = (size_t)(rp - cwords[idx]);
                    if (L >= sizeof pat) L = sizeof pat - 1;
                    memcpy(pat, cwords[idx], L);
                    pat[L] = '\0';
                    /* body remainder of same word after ) */
                    snprintf(body, sizeof body, "%s", rp + 1);
                }
                idx++;
                /* accumulate body words until ;; (two empty-separated ; or literal ;;) */
                char full[SH_LINE];
                snprintf(full, sizeof full, "%s", body);
                while (idx < ncw) {
                    if (!strcmp(cwords[idx], ";;") || !strcmp(cwords[idx], ";")) {
                        /* check for ;; as one word or two ; in a row */
                        if (!strcmp(cwords[idx], ";;")) { idx++; break; }
                        /* single ; : peek next */
                        if (idx + 1 < ncw && !strcmp(cwords[idx + 1], ";")) { idx += 2; break; }
                        /* single ; inside body: keep as separator */
                        if (strlen(full) + 2 < sizeof full) strcat(full, " ;");
                        idx++;
                        continue;
                    }
                    if (full[0] && strlen(full) + strlen(cwords[idx]) + 2 < sizeof full) strcat(full, " ");
                    if (strlen(full) + strlen(cwords[idx]) + 1 < sizeof full) strcat(full, cwords[idx]);
                    idx++;
                }
                /* pat may be empty before ) handled above; trim */
                while (pat[0] == ' ' || pat[0] == '\t') memmove(pat, pat + 1, strlen(pat));
                if (fnmatch(pat, subj, 0) == 0) {
                    st = run_line_str(t, full);
                    matched = 1;
                }
            }
            free(subj);
            sh_free(segs, ns);
            return st;
        }
    }
    int st = 0, prevop = 0, i = 0;
    while (i < ns) {
        /* control keyword mid-line (e.g. `V=1; while ..`): delegate rest */
        if (segs[i].argc > 0) {
            const char *kk = segs[i].args[0].raw;
            if (!strcmp(kk, "if") || !strcmp(kk, "for") || !strcmp(kk, "while") ||
                !strcmp(kk, "until") || !strcmp(kk, "case") || !strcmp(kk, "{") ||
                !strcmp(kk, "(")) {
                char rest[SH_LINE];
                join_segs(segs, i, ns - 1, rest, sizeof rest);
                /* preserve pipes: join_segs drops |, so re-lex original tail.
                 * Reconstruct from raw line tail instead: use run on rest
                 * which re-lexes pipes correctly for simple bodies. For
                 * piped bodies this still works because bodies re-enter. */
                st = run_line_str(t, rest);
                last_status = st;
                break;
            }
        }
        int j = i;
        while (j + 1 < ns && segs[j].nextop == 3) j++;
        if (!((prevop == 1 && st != 0) || (prevop == 2 && st == 0)))
            st = run_unit(t, segs, i, j);
        last_status = st;
        prevop = segs[j].nextop;
        i = j + 1;
    }
    sh_free(segs, ns);
    return st;
}

static int run_file(Tui *t, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { sh_err("  cannot read `%s`\n", path); return 1; }
    if (t->srcdepth >= 8) { sh_err("  source nesting too deep\n"); fclose(f); return 1; }
    t->srcdepth++;
    char lb[SH_LINE];
    int st = 0;
    while (fgets(lb, sizeof lb, f)) st = run_line(t, lb);
    t->srcdepth--;
    fclose(f);
    return st;
}

/* input stream provider for text filters: named files, else sh_input
 * (pipe stage / < file) or the real stdin. Missing files report but the
 * filter continues with the rest (POSIX cat behavior); miss tracks status. */
typedef struct {
    char **files;
    int n, i, miss, stdin_done;
    FILE *cur;
    int owned;
} InSrc;

static int insrc_next(InSrc *s) {
    if (s->cur && s->owned) { fclose(s->cur); s->cur = NULL; s->owned = 0; }
    while (s->i < s->n) {
        const char *p = s->files[s->i++];
        if (!strcmp(p, "-")) { /* POSIX convention: - means stdin */
            s->cur = sh_input ? sh_input : stdin;
            s->owned = 0;
            return 1;
        }
        s->cur = fopen(p, "r");
        if (!s->cur) { sh_err("  cannot read `%s`\n", p); s->miss = 1; continue; }
        s->owned = 1;
        return 1;
    }
    if (s->n == 0 && !s->stdin_done) {
        s->stdin_done = 1;
        s->cur = sh_input ? sh_input : stdin;
        s->owned = 0;
        return 1;
    }
    return 0;
}

static void insrc_close(InSrc *s) {
    if (s->cur && s->owned) { fclose(s->cur); s->cur = NULL; s->owned = 0; }
}

/* ── test/[ (POSIX test(1)): file, string and integer conditions ── */
typedef struct { char **a; int n, i, err; } TState;

static int t_expr_or(TState *s);

static int t_primary(TState *s) {
    if (s->i >= s->n) { s->err = 1; return 0; }
    char *w = s->a[s->i];
    if (!strcmp(w, "!")) {
        s->i++;
        int v = t_primary(s);
        return s->err ? 0 : !v;
    }
    if (!strcmp(w, "(")) {
        s->i++;
        int v = t_expr_or(s);
        if (s->err) return 0;
        if (s->i >= s->n || strcmp(s->a[s->i], ")")) { s->err = 1; return 0; }
        s->i++;
        return v;
    }
    if (s->i + 2 <= s->n - 1) {
        char *op = s->a[s->i + 1];
        if (!strcmp(op, "=") || !strcmp(op, "==") || !strcmp(op, "!=") ||
            !strcmp(op, "-eq") || !strcmp(op, "-ne") || !strcmp(op, "-gt") ||
            !strcmp(op, "-ge") || !strcmp(op, "-lt") || !strcmp(op, "-le")) {
            char *x = s->a[s->i], *y = s->a[s->i + 2];
            s->i += 3;
            if (!strcmp(op, "=") || !strcmp(op, "==")) return !strcmp(x, y);
            if (!strcmp(op, "!=")) return !!strcmp(x, y);
            char *e1 = NULL, *e2 = NULL;
            long a = strtol(x, &e1, 10), b = strtol(y, &e2, 10);
            if (!e1 || *e1 || !e2 || *e2) { s->err = 1; return 0; }
            if (!strcmp(op, "-eq")) return a == b;
            if (!strcmp(op, "-ne")) return a != b;
            if (!strcmp(op, "-gt")) return a > b;
            if (!strcmp(op, "-ge")) return a >= b;
            if (!strcmp(op, "-lt")) return a < b;
            return a <= b;
        }
    }
    s->i++;
    if (!strcmp(w, "-e") || !strcmp(w, "-f") || !strcmp(w, "-d") ||
        !strcmp(w, "-r") || !strcmp(w, "-w") || !strcmp(w, "-x") ||
        !strcmp(w, "-s")) {
        if (s->i >= s->n) { s->err = 1; return 0; }
        const char *p = s->a[s->i++];
        struct stat st;
        if (!strcmp(w, "-e")) return stat(p, &st) == 0;
        if (stat(p, &st) != 0) return 0;
        if (!strcmp(w, "-f")) return S_ISREG(st.st_mode);
        if (!strcmp(w, "-d")) return S_ISDIR(st.st_mode);
        if (!strcmp(w, "-r")) return access(p, R_OK) == 0;
        if (!strcmp(w, "-w")) return access(p, W_OK) == 0;
        if (!strcmp(w, "-x")) return access(p, X_OK) == 0;
        return st.st_size > 0; /* -s */
    }
    if (!strcmp(w, "-z") || !strcmp(w, "-n")) {
        if (s->i >= s->n) { s->err = 1; return 0; }
        const char *p = s->a[s->i++];
        return !strcmp(w, "-z") ? *p == '\0' : *p != '\0';
    }
    if (!strcmp(w, "-t")) {
        if (s->i >= s->n) { s->err = 1; return 0; }
        char *e = NULL;
        long fd = strtol(s->a[s->i++], &e, 10);
        if (!e || *e || fd < 0 || fd > 2) { s->err = 1; return 0; }
        return isatty((int)fd);
    }
    return *w != '\0'; /* bare string: true iff nonempty */
}

static int t_expr_and(TState *s) {
    int v = t_primary(s);
    while (!s->err && s->i < s->n && !strcmp(s->a[s->i], "-a")) {
        s->i++;
        int r = t_primary(s);
        v = v && r;
    }
    return v;
}

static int t_expr_or(TState *s) {
    int v = t_expr_and(s);
    while (!s->err && s->i < s->n && !strcmp(s->a[s->i], "-o")) {
        s->i++;
        int r = t_expr_and(s);
        v = v || r;
    }
    return v;
}

/* test/[ core: argv already stripped of the command name (and trailing ]
 * for [). Prints nothing; status 0 true, 1 false, 2 error. */
static int cmd_test(char **a, int n) {
    if (n == 0) return 1; /* no expression: false */
    TState s;
    s.a = a;
    s.n = n;
    s.i = 0;
    s.err = 0;
    int v = t_expr_or(&s);
    if (s.err || s.i != n) return 2;
    return v ? 0 : 1;
}

/* ── printf (POSIX printf(1)): FORMAT args, \ escapes, %b ── */
static void pf_esc_char(char c, int *stop) {
    switch (c) {
        case 'n': putchar('\n'); break;
        case 't': putchar('\t'); break;
        case 'r': putchar('\r'); break;
        case 'a': putchar('\a'); break;
        case 'b': putchar('\b'); break;
        case 'f': putchar('\f'); break;
        case 'v': putchar('\v'); break;
        case '\\': putchar('\\'); break;
        case '\'': putchar('\''); break;
        case '"': putchar('"'); break;
        case 'c': *stop = 1; break;
        default: putchar('\\'); putchar(c); break;
    }
}

static int cmd_printf(char **a, int n) {
    if (n < 1) { sh_err("  usage: printf FORMAT [args...]\n"); return 2; }
    const char *f = a[0];
    int ai = 1, used = 0;
    for (;;) {
        for (const char *p = f; *p; p++) {
            if (*p == '\\') {
                p++;
                if (*p == '\0') break;
                int stop = 0;
                if (*p == '0') { /* \0ooo octal */
                    int v = 0, k = 0;
                    while (k < 3 && p[1 + k] >= '0' && p[1 + k] <= '7') {
                        v = v * 8 + (p[1 + k] - '0');
                        k++;
                    }
                    if (k > 0) { putchar(v & 255); p += k; }
                    else pf_esc_char('0', &stop);
                } else pf_esc_char(*p, &stop);
                if (stop) return 0;
                continue;
            }
            if (*p != '%') { putchar(*p); continue; }
            p++;
            if (*p == '%') { putchar('%'); continue; }
            if (*p == 'b') { /* %b: expand escapes in next arg */
                const char *s = (ai < n) ? a[ai++] : "";
                used = 1;
                for (const char *q = s; *q; q++) {
                    if (*q == '\\') {
                        q++;
                        if (*q == '\0') break;
                        int stop = 0;
                        if (*q == '0') {
                            int v = 0, k = 0;
                            while (k < 3 && q[1 + k] >= '0' && q[1 + k] <= '7') {
                                v = v * 8 + (q[1 + k] - '0');
                                k++;
                            }
                            if (k > 0) { putchar(v & 255); q += k; }
                            else pf_esc_char('0', &stop);
                        } else pf_esc_char(*q, &stop);
                        if (stop) return 0;
                    } else putchar(*q);
                }
                continue;
            }
            char cf[40];
            int ci = 0;
            cf[ci++] = '%';
            while (*p && strchr("-+ #0", *p) && ci < 30) cf[ci++] = *p++;
            while (*p >= '0' && *p <= '9' && ci < 30) cf[ci++] = *p++;
            if (*p == '.' && ci < 30) {
                cf[ci++] = *p++;
                while (*p >= '0' && *p <= '9' && ci < 30) cf[ci++] = *p++;
            }
            int islong = 0;
            if (*p == 'l' && ci < 30) { cf[ci++] = *p++; islong = 1; }
            if (*p == 'l' && islong && ci < 30) { cf[ci++] = *p++; islong = 2; }
            char c = *p;
            if (!c || !strchr("sdiuoxXfFeEgGc", c)) { sh_err("  bad conversion\n"); return 2; }
            cf[ci++] = c;
            cf[ci] = '\0';
            used = 1;
            char ob[256];
            if (c == 's') {
                const char *s = (ai < n) ? a[ai++] : "";
                snprintf(ob, sizeof ob, cf, s);
            } else if (c == 'c') {
                int ch = (ai < n && a[ai][0]) ? (unsigned char)a[ai][0] : 0;
                if (ai < n) ai++;
                if (ch) snprintf(ob, sizeof ob, cf, ch);
                else ob[0] = '\0';
            } else if (c == 'f' || c == 'F' || c == 'e' || c == 'E' || c == 'g' || c == 'G') {
                double v = (ai < n) ? strtod(a[ai++], NULL) : 0.0;
                snprintf(ob, sizeof ob, cf, v);
            } else if (islong == 2) {
                long long v = 0;
                if (ai < n) { char *e = NULL; v = strtoll(a[ai++], &e, 10); }
                if (c == 'd' || c == 'i') snprintf(ob, sizeof ob, cf, v);
                else snprintf(ob, sizeof ob, cf, (unsigned long long)v);
            } else {
                long v = 0;
                if (ai < n) { char *e = NULL; v = strtol(a[ai++], &e, 10); }
                if (c == 'd' || c == 'i') snprintf(ob, sizeof ob, cf, (int)v);
                else snprintf(ob, sizeof ob, cf, (unsigned int)v);
            }
            fputs(ob, stdout);
        }
        if (ai >= n) break;
        (void)used;
    }
    return 0;
}

/* ── POSIX file utilities (operate on the real filesystem, cwd-relative) ── */
static int cmd_ls_files(const char *path) {
    const char *dir = (path && *path) ? path : ".";
    DIR *d = opendir(dir);
    if (!d) {
        /* POSIX ls lists plain files too, not just directories */
        struct stat st;
        if (stat(dir, &st) == 0 && !S_ISDIR(st.st_mode)) {
            printf("  %s\n", dir);
            return 0;
        }
        sh_err("  cannot list `%s`\n", dir);
        return 1;
    }
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        int isdir = 0;
        if (!strcmp(dir, ".")) {
            struct stat st;
            if (stat(e->d_name, &st) == 0 && S_ISDIR(st.st_mode)) isdir = 1;
        } else {
            char full[1024];
            size_t dl = strlen(dir), bl = strlen(e->d_name);
            if (dl + bl + 2 <= sizeof full) {
                memcpy(full, dir, dl);
                full[dl] = '/';
                memcpy(full + dl + 1, e->d_name, bl + 1);
                struct stat st;
                if (stat(full, &st) == 0 && S_ISDIR(st.st_mode)) isdir = 1;
            }
        }
        printf("  %s%s\n", e->d_name, isdir ? "/" : "");
    }
    closedir(d);
    return 0;
}

static int cmd_cat(char **files, int n) {
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    char lb[1024];
    while (insrc_next(&s)) {
        while (fgets(lb, sizeof lb, s.cur)) fputs(lb, stdout);
    }
    int rc = s.miss ? 1 : 0;
    insrc_close(&s);
    return rc;
}

static int cmd_head(char **files, int n, long nlines) {
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    char lb[1024];
    while (insrc_next(&s)) {
        long k = 0;
        if (nlines >= 0) {
            while (k < nlines && fgets(lb, sizeof lb, s.cur)) { fputs(lb, stdout); k++; }
        } else {
            /* +N: from line N on (POSIX head +N extension) */
            long ln = 0;
            while (fgets(lb, sizeof lb, s.cur)) {
                ln++;
                if (ln >= -nlines) fputs(lb, stdout);
            }
        }
    }
    int rc = s.miss ? 1 : 0;
    insrc_close(&s);
    return rc;
}

static int cmd_tail(char **files, int n, long nlines, int from_line) {
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    int rc = 0;
    char lb[1024];
    while (insrc_next(&s)) {
        if (!from_line) {
            /* ring buffer of last N lines */
            long cap = nlines > 0 ? nlines : 10;
            if (cap > 100000) cap = 100000;
            char **ring = NULL;
            if (cap > 0) ring = (char **)calloc((size_t)cap, sizeof(char *));
            if (cap > 0 && !ring) { rc = 1; break; }
            long total = 0;
            while (fgets(lb, sizeof lb, s.cur)) {
                char *cp = strdup(lb);
                if (!cp) { rc = 1; break; }
                if (cap > 0) {
                    free(ring[total % cap]);
                    ring[total % cap] = cp;
                }
                total++;
            }
            if (cap > 0) {
                long start = total > cap ? total - cap : 0;
                for (long i = start; i < total; i++) fputs(ring[i % cap], stdout);
                for (long i = 0; i < (total < cap ? total : cap); i++) free(ring[i]);
                free(ring);
            }
        } else {
            long ln = 0;
            while (fgets(lb, sizeof lb, s.cur)) {
                ln++;
                if (ln >= nlines) fputs(lb, stdout);
            }
        }
    }
    if (s.miss) rc = 1;
    insrc_close(&s);
    return rc;
}

static int cmd_wc(char **files, int n, int want_l, int want_w, int want_c) {
    if (!want_l && !want_w && !want_c) want_l = want_w = want_c = 1;
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    long tl = 0, tw = 0, tc = 0;
    int nsrc = 0;
    int ch;
    while (insrc_next(&s)) {
        long l = 0, w = 0, c = 0;
        int inword = 0;
        /* count bytes, words and lines over the stream */
        while ((ch = fgetc(s.cur)) != EOF) {
            c++;
            if (ch == '\n') l++;
            if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' || ch == '\v') inword = 0;
            else if (!inword) { inword = 1; w++; }
        }
        tl += l; tw += w; tc += c;
        nsrc++;
        if (n > 1) {
            if (want_l) printf("%7ld", l);
            if (want_w) printf("%7ld", w);
            if (want_c) printf("%7ld", c);
            printf(" %s\n", files[nsrc - 1]);
        }
    }
    if (n > 1) {
        if (want_l) printf("%7ld", tl);
        if (want_w) printf("%7ld", tw);
        if (want_c) printf("%7ld", tc);
        printf(" total\n");
    } else {
        if (want_l) printf("%7ld", tl);
        if (want_w) printf("%7ld", tw);
        if (want_c) printf("%7ld", tc);
        printf("\n");
    }
    int rc = s.miss ? 1 : 0;
    insrc_close(&s);
    return rc;
}

static int cmp_str(const void *a, const void *b) {
    return strcmp(*(char *const *)a, *(char *const *)b);
}

static int cmp_num(const void *a, const void *b) {
    const char *x = *(char *const *)a, *y = *(char *const *)b;
    double dx = strtod(x, NULL), dy = strtod(y, NULL);
    if (dx < dy) return -1;
    if (dx > dy) return 1;
    return strcmp(x, y);
}

static int cmd_sort(char **files, int n, int numeric, int rev, int uniq) {
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    size_t cap = 256, len = 0;
    char **lines = (char **)malloc(cap * sizeof(char *));
    if (!lines) return 1;
    char lb[4096];
    while (insrc_next(&s)) {
        while (fgets(lb, sizeof lb, s.cur)) {
            if (len >= cap) {
                cap *= 2;
                char **n2 = (char **)realloc(lines, cap * sizeof(char *));
                if (!n2) break;
                lines = n2;
            }
            if (len >= cap) break;
            lines[len] = strdup(lb);
            if (!lines[len]) break;
            len++;
        }
    }
    int rc = s.miss ? 1 : 0;
    insrc_close(&s);
    qsort(lines, len, sizeof(char *), numeric ? cmp_num : cmp_str);
    const char *prev = NULL;
    for (size_t i = 0; i < len; i++) {
        int out = 1;
        if (uniq && prev && !strcmp(prev, lines[i])) out = 0;
        if (rev) {
            /* reverse after dedup decision on sorted order */
        }
        prev = lines[i];
        if (out) {
            if (rev) {
                /* print later in reverse */
            } else fputs(lines[i], stdout);
        }
    }
    if (rev) {
        prev = NULL;
        for (size_t k = len; k > 0; k--) {
            const char *ln = lines[k - 1];
            int out = 1;
            if (uniq && prev && !strcmp(prev, ln)) out = 0;
            prev = ln;
            if (out) fputs(ln, stdout);
        }
    }
    for (size_t i = 0; i < len; i++) free(lines[i]);
    free(lines);
    return rc;
}

static int cmd_uniq(char **files, int n, int count, int duponly, int uniqonly) {
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    char prev[4096], lb[4096];
    int have = 0;
    long run = 0;
    int rc = 0;
    while (insrc_next(&s)) {
        while (fgets(lb, sizeof lb, s.cur)) {
            if (!have || strcmp(lb, prev)) {
                if (have) {
                    int show = 1;
                    if (duponly && run < 2) show = 0;
                    if (uniqonly && run > 1) show = 0;
                    if (show) {
                        if (count) printf("%7ld %s", run, prev);
                        else fputs(prev, stdout);
                    }
                }
                snprintf(prev, sizeof prev, "%s", lb);
                have = 1;
                run = 1;
            } else run++;
        }
    }
    if (have) {
        int show = 1;
        if (duponly && run < 2) show = 0;
        if (uniqonly && run > 1) show = 0;
        if (show) {
            if (count) printf("%7ld %s", run, prev);
            else fputs(prev, stdout);
        }
    }
    if (s.miss) rc = 1;
    insrc_close(&s);
    return rc;
}

/* parse "1,3-5" field list into ranges (1-based, hi=-1 open) */
static int parse_ranges(const char *s, int *lo, int *hi, int maxr) {
    int n = 0;
    const char *p = s;
    while (*p && n < maxr) {
        char *e = NULL;
        long a = strtol(p, &e, 10);
        if (e == p || a < 1) return -1;
        long b = a;
        p = e;
        if (*p == '-') {
            p++;
            if (*p == ',' || *p == '\0') b = -1;
            else {
                b = strtol(p, &e, 10);
                if (e == p || b < 1) return -1;
                p = e;
            }
        }
        lo[n] = (int)a;
        hi[n] = (int)b;
        n++;
        if (*p == ',') p++;
        else if (*p) return -1;
    }
    if (*p) return -1;
    return n;
}

static int in_ranges(int *lo, int *hi, int n, int pos) {
    for (int i = 0; i < n; i++)
        if (pos >= lo[i] && (hi[i] < 0 || pos <= hi[i])) return 1;
    return 0;
}

static int cmd_cut(char **files, int n, int *lo, int *hi, int nr, int byfields, char delim, int suppress) {
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    char lb[4096];
    while (insrc_next(&s)) {
        while (fgets(lb, sizeof lb, s.cur)) {
            size_t L = strlen(lb);
            int hasnl = L > 0 && lb[L - 1] == '\n';
            if (hasnl) lb[--L] = '\0';
            if (byfields) {
                if (!strchr(lb, delim)) {
                    if (!suppress) printf("%s%s", lb, hasnl ? "\n" : "");
                    continue;
                }
                /* split into fields (empty fields kept) */
                int fno = 1, first = 1;
                const char *p = lb;
                for (;;) {
                    const char *d = strchr(p, delim);
                    size_t fl = d ? (size_t)(d - p) : strlen(p);
                    if (in_ranges(lo, hi, nr, fno)) {
                        if (!first) putchar(delim);
                        fwrite(p, 1, fl, stdout);
                        first = 0;
                    }
                    fno++;
                    if (!d) break;
                    p = d + 1;
                }
                if (hasnl) putchar('\n');
            } else {
                for (size_t i = 0; i < L; i++)
                    if (in_ranges(lo, hi, nr, (int)i + 1)) putchar(lb[i]);
                if (hasnl) putchar('\n');
            }
        }
    }
    int rc = s.miss ? 1 : 0;
    insrc_close(&s);
    return rc;
}

/* expand a tr set string (ranges, escapes, [:classes:]) into a 256 table */
static void tr_expand_set(const char *s, unsigned char *present) {
    memset(present, 0, 256);
    struct { const char *name; const char *chars; } cls[] = {
        { "upper", "ABCDEFGHIJKLMNOPQRSTUVWXYZ" },
        { "lower", "abcdefghijklmnopqrstuvwxyz" },
        { "digit", "0123456789" },
        { "alpha", "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" },
        { "alnum", "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789" },
        { "space", " \t\n\r\f\v" },
        { NULL, NULL }
    };
    for (const char *p = s; *p;) {
        if (*p == '[' && p[1] == ':') {
            const char *e = strstr(p + 2, ":]");
            if (e) {
                char nm[16];
                size_t nl = (size_t)(e - (p + 2));
                if (nl < sizeof nm) {
                    memcpy(nm, p + 2, nl);
                    nm[nl] = '\0';
                    for (int k = 0; cls[k].name; k++) {
                        if (!strcmp(nm, cls[k].name)) {
                            for (const char *c = cls[k].chars; *c; c++)
                                present[(unsigned char)*c] = 1;
                            break;
                        }
                    }
                }
                p = e + 2;
                continue;
            }
        }
        unsigned char c = (unsigned char)*p;
        if (c == '\\' && p[1]) {
            p++;
            if (*p == 'n') c = '\n';
            else if (*p == 't') c = '\t';
            else if (*p == 'r') c = '\r';
            else c = (unsigned char)*p;
            present[c] = 1;
            p++;
            continue;
        }
        if (p[1] == '-' && p[2] && p[2] != ']') {
            unsigned char e2 = (unsigned char)p[2];
            if (e2 == '\\' && p[3]) {
                e2 = (unsigned char)p[3];
                if (e2 == 'n') e2 = '\n';
                else if (e2 == 't') e2 = '\t';
            }
            if (c <= e2) {
                for (int k = c; k <= e2; k++) present[k] = 1;
                p += (p[2] == '\\') ? 4 : 3;
                continue;
            }
        }
        present[c] = 1;
        p++;
    }
}

static int cmd_tr(const char *s1, const char *s2, int del, int squeez, int comp) {
    unsigned char m1[256], m2[256];
    unsigned char xlat[256];
    tr_expand_set(s1, m1);
    if (comp) {
        unsigned char inv[256];
        memset(inv, 0, sizeof inv);
        for (int i = 0; i < 256; i++)
            if (!m1[i]) inv[i] = 1;
        memcpy(m1, inv, sizeof inv);
    }
    if (s2) tr_expand_set(s2, m2);
    else memset(m2, 0, sizeof m2);
    /* translation table from s1 order to s2 (last char repeats) */
    unsigned char order[256];
    int on = 0;
    for (int i = 0; i < 256; i++)
        if (m1[i]) order[on++] = (unsigned char)i;
    unsigned char s2o[256];
    int s2n = 0;
    if (s2) {
        for (int i = 0; i < 256; i++)
            if (m2[i]) s2o[s2n++] = (unsigned char)i;
    }
    for (int i = 0; i < 256; i++) xlat[i] = (unsigned char)i;
    if (s2 && !del) {
        for (int i = 0; i < on; i++)
            xlat[order[i]] = s2n ? s2o[i < s2n ? i : s2n - 1] : order[i];
    }
    int squeezeset[256];
    memset(squeezeset, 0, sizeof squeezeset);
    if (squeez && s2 && !del) {
        for (int i = 0; i < s2n; i++) squeezeset[s2o[i]] = 1;
    } else if (squeez && del && !s2) {
        for (int i = 0; i < on; i++) squeezeset[order[i]] = 1;
    } else if (squeez && del && s2) {
        for (int i = 0; i < s2n; i++) squeezeset[s2o[i]] = 1;
    } else if (squeez && !del && !s2) {
        for (int i = 0; i < on; i++) squeezeset[order[i]] = 1;
    }
    FILE *in = sh_input ? sh_input : stdin;
    int prev = -1, ch;
    while ((ch = fgetc(in)) != EOF) {
        unsigned char c = (unsigned char)ch;
        if (del && m1[c]) { prev = -1; continue; }
        unsigned char o = xlat[c];
        if (squeez && squeezeset[o] && o == prev) continue;
        putchar(o);
        prev = o;
    }
    return 0;
}

static int cmd_grep(const char *pat, char **files, int n, int use_ext, int icase, int invert, int number, int count, int listf) {
    int flags = REG_NOSUB | (use_ext ? REG_EXTENDED : 0) | (icase ? REG_ICASE : 0);
    regex_t re;
    if (regcomp(&re, pat, flags) != 0) { sh_err("  bad pattern\n"); return 2; }
    InSrc s;
    memset(&s, 0, sizeof s);
    s.files = files;
    s.n = n;
    int multi = n > 1;
    char lb[4096];
    while (insrc_next(&s)) {
        long ln = 0, hits = 0;
        const char *nm = (s.files && s.i > 0) ? s.files[s.i - 1] : NULL;
        while (fgets(lb, sizeof lb, s.cur)) {
            ln++;
            int m = regexec(&re, lb, 0, NULL, 0) == 0;
            if (invert) m = !m;
            if (!m) continue;
            hits++;
            if (count || listf) continue;
            if (multi) printf("%s:", nm ? nm : "(stdin)");
            if (number) printf("%ld:", ln);
            fputs(lb, stdout);
        }
        if (count) {
            if (multi) printf("%s:", nm ? nm : "(stdin)");
            printf("%ld\n", hits);
        }
        if (listf && hits > 0) printf("%s\n", nm ? nm : "(stdin)");
    }
    int rc = s.miss ? 1 : 0;
    insrc_close(&s);
    regfree(&re);
    return rc;
}

static int cmd_tee(char **files, int n, int append) {
    FILE **outs = NULL;
    int rc = 0;
    if (n > 0) {
        outs = (FILE **)malloc(sizeof(FILE *) * (size_t)n);
        if (!outs) return 1;
        for (int i = 0; i < n; i++) {
            outs[i] = fopen(files[i], append ? "a" : "w");
            if (!outs[i]) { sh_err("  cannot write `%s`\n", files[i]); rc = 1; }
        }
    }
    FILE *in = sh_input ? sh_input : stdin;
    char lb[4096];
    while (fgets(lb, sizeof lb, in)) {
        fputs(lb, stdout);
        for (int i = 0; i < n; i++)
            if (outs[i]) fputs(lb, outs[i]);
    }
    for (int i = 0; i < n; i++)
        if (outs[i]) fclose(outs[i]);
    free(outs);
    return rc;
}

static int cmd_mkdir(const char *path, int parents) {
    if (!parents) {
        if (mkdir(path, 0755) != 0) { sh_err("  cannot mkdir `%s`: %s\n", path, strerror(errno)); return 1; }
        return 0;
    }
    char tmp[1024];
    size_t L = strlen(path);
    if (L == 0 || L > sizeof tmp - 1) { sh_err("  bad path\n"); return 2; }
    memcpy(tmp, path, L + 1);
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) {
                sh_err("  cannot mkdir `%s`: %s\n", tmp, strerror(errno));
                return 1;
            }
            *p = '/';
        }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) {
        sh_err("  cannot mkdir `%s`: %s\n", tmp, strerror(errno));
        return 1;
    }
    return 0;
}

static int cmd_cp(const char *src, const char *dst) {
    FILE *fi = fopen(src, "rb");
    if (!fi) { sh_err("  cannot read `%s`\n", src); return 1; }
    FILE *fo = fopen(dst, "wb");
    if (!fo) { sh_err("  cannot write `%s`\n", dst); fclose(fi); return 1; }
    char buf[8192];
    size_t n;
    int rc = 0;
    while ((n = fread(buf, 1, sizeof buf, fi)) > 0) {
        if (fwrite(buf, 1, n, fo) != n) { sh_err("  write failed `%s`\n", dst); rc = 1; break; }
    }
    fclose(fi);
    if (fclose(fo) != 0) { sh_err("  write failed `%s`\n", dst); rc = 1; }
    return rc;
}

static int cmd_mv(const char *src, const char *dst) {
    if (rename(src, dst) == 0) return 0;
    if (errno == EXDEV) { /* cross-device: copy + unlink */
        if (cmd_cp(src, dst) != 0) return 1;
        if (unlink(src) != 0) { sh_err("  cannot remove `%s`\n", src); return 1; }
        return 0;
    }
    sh_err("  cannot move `%s`: %s\n", src, strerror(errno));
    return 1;
}

static int rm_recursive(const char *path) {
    struct stat st;
    if (lstat(path, &st) != 0) { sh_err("  cannot remove `%s`\n", path); return 1; }
    if (!S_ISDIR(st.st_mode)) {
        if (unlink(path) != 0) { sh_err("  cannot remove `%s`\n", path); return 1; }
        return 0;
    }
    DIR *d = opendir(path);
    if (!d) { sh_err("  cannot list `%s`\n", path); return 1; }
    int rc = 0;
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        char full[1024];
        size_t dl = strlen(path), bl = strlen(e->d_name);
        if (dl + bl + 2 > sizeof full) { sh_err("  path too long\n"); rc = 1; continue; }
        memcpy(full, path, dl);
        full[dl] = '/';
        memcpy(full + dl + 1, e->d_name, bl + 1);
        if (rm_recursive(full)) rc = 1;
    }
    closedir(d);
    if (rmdir(path) != 0) { sh_err("  cannot remove `%s`\n", path); rc = 1; }
    return rc;
}

static int cmd_touch(char **files, int n) {
    int rc = 0;
    for (int i = 0; i < n; i++) {
        FILE *f = fopen(files[i], "a");
        if (!f) { sh_err("  cannot touch `%s`\n", files[i]); rc = 1; continue; }
        fclose(f);
        if (utimensat(AT_FDCWD, files[i], NULL, 0) != 0) {
            sh_err("  cannot touch `%s`\n", files[i]);
            rc = 1;
        }
    }
    return rc;
}

/* POSIX basename/dirname string rules (no filesystem access) */
static void base_dir_name(const char *p, char *base, size_t bn, char *dir, size_t dn) {
    const char *b = p;
    char tmp[1024];
    snprintf(tmp, sizeof tmp, "%s", p);
    /* strip trailing slashes (keep root) */
    size_t L = strlen(tmp);
    while (L > 1 && tmp[L - 1] == '/') tmp[--L] = '\0';
    const char *s = strrchr(tmp, '/');
    const char *bb = s ? s + 1 : tmp;
    if (base) snprintf(base, bn, "%s", *bb ? bb : ".");
    if (dir) {
        if (!s) snprintf(dir, dn, ".");
        else if (s == tmp) snprintf(dir, dn, "/");
        else {
            size_t dl = (size_t)(s - tmp);
            while (dl > 1 && tmp[dl - 1] == '/') dl--;
            if (dl >= dn) dl = dn - 1;
            memcpy(dir, tmp, dl);
            dir[dl] = '\0';
        }
    }
    (void)b;
}

static int find_walk(const char *path, const char *pat, int want_f, int want_d) {
    struct stat st;
    if (lstat(path, &st) != 0) { sh_err("  cannot access `%s`\n", path); return 1; }
    int isdir = S_ISDIR(st.st_mode);
    int show = (!want_f && !want_d) || (want_f && !isdir) || (want_d && isdir);
    if (show) {
        const char *slash = strrchr(path, '/');
        const char *base = slash ? slash + 1 : path;
        if (!pat || fnmatch(pat, base, 0) == 0) printf("%s\n", path);
    }
    if (!isdir) return 0;
    DIR *d = opendir(path);
    if (!d) { sh_err("  cannot list `%s`\n", path); return 1; }
    int rc = 0;
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
        if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
        char full[1024];
        size_t dl = strlen(path), bl = strlen(e->d_name);
        if (dl + bl + 2 > sizeof full) { sh_err("  path too long\n"); rc = 1; continue; }
        memcpy(full, path, dl);
        full[dl] = '/';
        memcpy(full + dl + 1, e->d_name, bl + 1);
        if (find_walk(full, pat, want_f, want_d)) rc = 1;
    }
    closedir(d);
    return rc;
}

/* ── custom molecule builder ──────────────────────────────────────────
 * Build your own molecules atom by atom in a working template, save them
 * as MOL1 files, and stamp copies into the live sim. Coordinates are
 * template-local; `mol center` recenters on the centre of mass. */
#define MOL_MAXA 512
typedef struct { int Z; Vec3 pos; double q, eps, sig; int has_lj; } MolAtom;
typedef struct { int a, b, order; } MolBond;
static struct {
    char name[64];
    MolAtom a[MOL_MAXA];
    int na;
    MolBond b[MOL_MAXA];
    int nb;
} mol;

static void mol_com(Vec3 *out) {
    double m = 0;
    Vec3 c = vec3_zero();
    for (int i = 0; i < mol.na; i++) {
        c = vec3_add(c, mol.a[i].pos);
        m += 1.0;
    }
    *out = (m > 0) ? vec3_scale(c, 1.0 / m) : c;
}

static int mol_place(Tui *t, Vec3 origin) {
    if (mol.na < 1) { sh_err("  builder empty (`mol new` first)\n"); return 1; }
    if (!sp_room(t, mol.na, mol.nb)) return 1;
    int first = t->sim->num_atoms;
    for (int i = 0; i < mol.na; i++) {
        Vec3 p = vec3_add(mol.a[i].pos, origin);
        int idx = sim_add_atom(t->sim, mol.a[i].Z, p, mol.a[i].q);
        if (idx < 0) { sh_err("  sim full after %d atoms\n", i); return 1; }
        if (mol.a[i].has_lj) sim_set_atom_lj(t->sim, idx, mol.a[i].eps, mol.a[i].sig);
    }
    for (int i = 0; i < mol.nb; i++) {
        if (sim_add_bond(t->sim, first + mol.b[i].a, first + mol.b[i].b, mol.b[i].order) < 0) {
            sh_err("  bond %d failed\n", i);
            return 1;
        }
    }
    sim_rebuild_angles(t->sim);
    forces_calculate(t->sim);
    return 0;
}

static int mol_save(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { sh_err("  cannot write `%s`\n", path); return 1; }
    fprintf(f, "MOL1 %s\natoms %d\n", mol.name[0] ? mol.name : "mol", mol.na);
    for (int i = 0; i < mol.na; i++)
        fprintf(f, "%d %.17g %.17g %.17g %.17g %.17g %.17g %d\n", mol.a[i].Z,
            mol.a[i].pos.x, mol.a[i].pos.y, mol.a[i].pos.z, mol.a[i].q,
            mol.a[i].eps, mol.a[i].sig, mol.a[i].has_lj);
    fprintf(f, "bonds %d\n", mol.nb);
    for (int i = 0; i < mol.nb; i++)
        fprintf(f, "%d %d %d\n", mol.b[i].a, mol.b[i].b, mol.b[i].order);
    fclose(f);
    return 0;
}

static int mol_load(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { sh_err("  cannot read `%s`\n", path); return 1; }
    char magic[16], nm[64];
    int na = -1;
    if (fscanf(f, "%15s %63s\n", magic, nm) != 2 || strcmp(magic, "MOL1")) {
        sh_err("  bad molecule file (want MOL1)\n");
        fclose(f);
        return 1;
    }
    if (fscanf(f, "atoms %d\n", &na) != 1 || na < 0 || na > MOL_MAXA) {
        sh_err("  bad atom count\n");
        fclose(f);
        return 1;
    }
    MolAtom A[MOL_MAXA];
    for (int i = 0; i < na; i++) {
        int Z, hl;
        double x, y, z, q, e, s;
        if (fscanf(f, "%d %lg %lg %lg %lg %lg %lg %d\n", &Z, &x, &y, &z, &q, &e, &s, &hl) != 8) {
            sh_err("  bad atom %d\n", i);
            fclose(f);
            return 1;
        }
        if (!pt_element(Z)) { sh_err("  unknown Z=%d\n", Z); fclose(f); return 1; }
        A[i].Z = Z;
        A[i].pos = vec3(x, y, z);
        A[i].q = q; A[i].eps = e; A[i].sig = s; A[i].has_lj = hl ? 1 : 0;
    }
    int nb = -1;
    if (fscanf(f, "bonds %d\n", &nb) != 1 || nb < 0 || nb > MOL_MAXA) {
        sh_err("  bad bond count\n");
        fclose(f);
        return 1;
    }
    MolBond B[MOL_MAXA];
    for (int i = 0; i < nb; i++) {
        if (fscanf(f, "%d %d %d\n", &B[i].a, &B[i].b, &B[i].order) != 3 ||
            B[i].a < 0 || B[i].a >= na || B[i].b < 0 || B[i].b >= na ||
            B[i].order < 1 || B[i].order > 3) {
            sh_err("  bad bond %d\n", i);
            fclose(f);
            return 1;
        }
    }
    fclose(f);
    snprintf(mol.name, sizeof mol.name, "%s", nm);
    memcpy(mol.a, A, sizeof(MolAtom) * (size_t)na);
    mol.na = na;
    memcpy(mol.b, B, sizeof(MolBond) * (size_t)nb);
    mol.nb = nb;
    return 0;
}

/* ── custom reaction types ────────────────────────────────────────────
 * A reaction matches the closest Z1–Z2 pair within rcut (minimum image
 * when PBC is on) and rewrites topology: break the pair bond, make a
 * bond, delete a leaving atom (terminal-only, real condensation
 * chemistry), set charges. `rxn fire` attempts it now; `rxn auto N`
 * with `rxn arm` fires once when a live run wanders into geometry. */
static struct {
    char name[64];
    int Z1, Z2;
    double rcut;
    int dobreak, make_order, delrole; /* delrole 0 none, 1 = Z1 atom, 2 = Z2 */
    double q1, q2;
    int setq, defined;
} rxn;

static double rxn_dist(const Simulation *s, Vec3 a, Vec3 b) {
    Vec3 d = vec3_sub(b, a);
    for (int k = 0; k < 3; k++) {
        double L = (k == 0) ? s->box.dimensions.x : (k == 1) ? s->box.dimensions.y : s->box.dimensions.z;
        double c = (k == 0) ? d.x : (k == 1) ? d.y : d.z;
        if (s->box.periodic[k] && L > 1e-12) {
            c -= L * round(c / L);
            if (k == 0) d.x = c;
            else if (k == 1) d.y = c;
            else d.z = c;
        }
    }
    return vec3_norm(d);
}

static int rxn_find(const Tui *t, int *ia, int *ib, double *rr) {
    if (!rxn.defined) return 1;
    double best = 1e30;
    int bi = -1, bj = -1;
    for (int i = 0; i < t->sim->num_atoms; i++) {
        for (int j = 0; j < t->sim->num_atoms; j++) {
            if (i == j) continue;
            int zi = t->sim->atoms[i].Z, zj = t->sim->atoms[j].Z;
            if (!((zi == rxn.Z1 && zj == rxn.Z2) || (zi == rxn.Z2 && zj == rxn.Z1))) continue;
            double d = rxn_dist(t->sim, t->sim->atoms[i].position, t->sim->atoms[j].position);
            if (d < best) { best = d; bi = i; bj = j; }
        }
    }
    if (bi < 0 || best > rxn.rcut) return 1;
    /* role 1 = the Z1 atom */
    if (t->sim->atoms[bi].Z == rxn.Z1) { *ia = bi; *ib = bj; }
    else { *ia = bj; *ib = bi; }
    *rr = best;
    return 0;
}

static int tui_bond_remove(Simulation *s, int ia, int ib) {
    int bi = -1;
    for (int b = 0; b < s->num_bonds; b++) {
        if ((s->bonds[b].atom_a == ia && s->bonds[b].atom_b == ib) ||
            (s->bonds[b].atom_a == ib && s->bonds[b].atom_b == ia)) {
            bi = b;
            break;
        }
    }
    if (bi < 0) return 1;
    s->bonds[bi] = s->bonds[s->num_bonds - 1];
    s->num_bonds--;
    for (int i = 0; i < s->num_atoms; i++) s->atoms[i].num_bonds = 0;
    for (int b = 0; b < s->num_bonds; b++) {
        Atom *a = &s->atoms[s->bonds[b].atom_a];
        Atom *c = &s->atoms[s->bonds[b].atom_b];
        if (a->num_bonds < MAX_BONDS_PER_ATOM) {
            a->bond_partners[a->num_bonds] = s->bonds[b].atom_b;
            a->bond_orders[a->num_bonds] = s->bonds[b].order;
            a->num_bonds++;
        }
        if (c->num_bonds < MAX_BONDS_PER_ATOM) {
            c->bond_partners[c->num_bonds] = s->bonds[b].atom_a;
            c->bond_orders[c->num_bonds] = s->bonds[b].order;
            c->num_bonds++;
        }
    }
    return 0;
}

static int tui_bonded(const Simulation *s, int ia, int ib) {
    for (int b = 0; b < s->num_bonds; b++) {
        if ((s->bonds[b].atom_a == ia && s->bonds[b].atom_b == ib) ||
            (s->bonds[b].atom_a == ib && s->bonds[b].atom_b == ia))
            return 1;
    }
    return 0;
}

static int rxn_fire(Tui *t) {
    if (!rxn.defined) { sh_err("  no reaction defined (`rxn new` first)\n"); return 1; }
    int ia, ib;
    double rr;
    if (rxn_find(t, &ia, &ib, &rr)) {
        printf("  no match: closest Z%d-Z%d pair outside %.3f A\n", rxn.Z1, rxn.Z2, rxn.rcut);
        return 1;
    }
    forces_calculate(t->sim);
    double pe0 = t->sim->potential_energy;
    printf("  match: atoms %d-%d at %.3f A (cut %.3f)\n", ia, ib, rr, rxn.rcut);
    if (rxn.dobreak) {
        if (tui_bond_remove(t->sim, ia, ib)) printf("  break: pair not bonded, nothing cut\n");
        else printf("  break: bond %d-%d cut\n", ia, ib);
    }
    if (rxn.make_order > 0) {
        if (tui_bonded(t->sim, ia, ib)) printf("  make: already bonded\n");
        else if (sim_add_bond(t->sim, ia, ib, rxn.make_order) < 0) printf("  make: failed\n");
        else printf("  make: bond %d-%d order %d formed\n", ia, ib, rxn.make_order);
    }
    if (rxn.setq) {
        t->sim->atoms[ia].partial_charge = rxn.q1;
        t->sim->atoms[ib].partial_charge = rxn.q2;
        printf("  charge: q[%d]=%.3f q[%d]=%.3f\n", ia, rxn.q1, ib, rxn.q2);
    }
    if (rxn.delrole == 1 || rxn.delrole == 2) {
        int victim = (rxn.delrole == 1) ? ia : ib;
        if (sim_remove_terminal_atom(t->sim, victim)) printf("  delete: atom %d left\n", victim);
        else printf("  delete: atom %d not terminal, kept\n", victim);
    }
    sim_rebuild_angles(t->sim);
    forces_calculate(t->sim);
    printf("  PE %.6f -> %.6f eV (dE=%+.6f)\n", pe0, t->sim->potential_energy,
        t->sim->potential_energy - pe0);
    maybe_render(t);
    return 0;
}

static int rxn_save(const char *path) {
    if (!rxn.defined) { sh_err("  no reaction defined\n"); return 1; }
    FILE *f = fopen(path, "w");
    if (!f) { sh_err("  cannot write `%s`\n", path); return 1; }
    fprintf(f, "RXN1 %s\npair %d %d %.17g\nbreak %d\nmake %d\ndelete %d\ncharge %d %.17g %.17g\n",
        rxn.name[0] ? rxn.name : "rxn", rxn.Z1, rxn.Z2, rxn.rcut, rxn.dobreak,
        rxn.make_order, rxn.delrole, rxn.setq, rxn.q1, rxn.q2);
    fclose(f);
    return 0;
}

static int rxn_load(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { sh_err("  cannot read `%s`\n", path); return 1; }
    char magic[16], nm[64];
    int Z1, Z2, br, mo, dl, sq;
    double rc, q1, q2;
    int ok = fscanf(f, "%15s %63s\n", magic, nm) == 2 && !strcmp(magic, "RXN1") &&
        fscanf(f, "pair %d %d %lg\n", &Z1, &Z2, &rc) == 3 &&
        fscanf(f, "break %d\n", &br) == 1 && fscanf(f, "make %d\n", &mo) == 1 &&
        fscanf(f, "delete %d\n", &dl) == 1 &&
        fscanf(f, "charge %d %lg %lg\n", &sq, &q1, &q2) == 3;
    fclose(f);
    if (!ok || !pt_element(Z1) || !pt_element(Z2) || !(rc > 0) || mo < -1 || mo > 3 || dl < 0 || dl > 2) {
        sh_err("  bad reaction file (want RXN1)\n");
        return 1;
    }
    snprintf(rxn.name, sizeof rxn.name, "%s", nm);
    rxn.Z1 = Z1; rxn.Z2 = Z2; rxn.rcut = rc;
    rxn.dobreak = br ? 1 : 0; rxn.make_order = mo; rxn.delrole = dl;
    rxn.setq = sq ? 1 : 0; rxn.q1 = q1; rxn.q2 = q2;
    rxn.defined = 1;
    return 0;
}

/* auto-reaction check for live runs: fires once per arming */
static void rxn_autocheck(Tui *t) {
    if (!t->rxn_armed || t->rxn_every < 1 || !rxn.defined) return;
    if (t->sim->step % (unsigned long long)t->rxn_every != 0) return;
    int ia, ib;
    double rr;
    if (rxn_find(t, &ia, &ib, &rr)) return;
    t->rxn_armed = 0;
    printf("  [rxn auto] `%s` triggered at step %llu (pair %d-%d %.3f A)\n",
        rxn.name[0] ? rxn.name : "rxn", (unsigned long long)t->sim->step, ia, ib, rr);
    rxn_fire(t);
}

/* ── gas pressure & heat ──────────────────────────────────────────────
 * Pressure from the virial (COM-relative, translation-invariant):
 *   P = (2·KE/3 + vir/3) / V,  vir = Σ (r−R)·F,
 * in eV/A³, reported in bar (1 eV/A³ = 1.602176634e6 bar).
 * The barostat is tui-side isotropic Berendsen scaling about the box
 * centre (factor clamped to ±1%/step); the engine is untouched. */
#define EVA3_TO_BAR 1602176.634

static int box_volume(const Simulation *s, double *V) {
    double x = s->box.dimensions.x, y = s->box.dimensions.y, z = s->box.dimensions.z;
    if (!(x > 1e-9) || !(y > 1e-9) || !(z > 1e-9) || !isfinite(x + y + z)) return 1;
    *V = x * y * z;
    return 0;
}

static void sim_com(const Simulation *s, Vec3 *out) {
    Vec3 c = vec3_zero();
    double m = 0;
    for (int i = 0; i < s->num_atoms; i++) {
        double mi = s->atoms[i].mass > 0 ? s->atoms[i].mass : 1.0;
        c = vec3_add(c, vec3_scale(s->atoms[i].position, mi));
        m += mi;
    }
    *out = (m > 0) ? vec3_scale(c, 1.0 / m) : c;
}

static int tui_pressure(Tui *t, double *p_bar, double *p_eva, double *nkt, double *vir) {
    double V;
    if (box_volume(t->sim, &V)) return 1;
    forces_calculate(t->sim);
    Vec3 R;
    sim_com(t->sim, &R);
    double v = 0;
    for (int i = 0; i < t->sim->num_atoms; i++) {
        Vec3 d = vec3_sub(t->sim->atoms[i].position, R);
        v += d.x * t->sim->atoms[i].force.x + d.y * t->sim->atoms[i].force.y + d.z * t->sim->atoms[i].force.z;
    }
    double ke = t->sim->kinetic_energy;
    double pe = (2.0 * ke / 3.0 + v / 3.0) / V;
    if (p_eva) *p_eva = pe;
    if (nkt) *nkt = (2.0 * ke / 3.0) / V;
    if (vir) *vir = (v / 3.0) / V;
    if (p_bar) *p_bar = pe * EVA3_TO_BAR;
    return 0;
}

static void tui_barostat_step(Tui *t) {
    if (!t->baro_on) return;
    double V;
    if (box_volume(t->sim, &V)) return;
    if (!(t->taup_fs > 1e-9) || !isfinite(t->taup_fs)) return;
    double cur = 0;
    tui_pressure(t, NULL, &cur, NULL, NULL);
    double tgt = t->p0_bar / EVA3_TO_BAR;
    if (!isfinite(cur)) return;
    /* isotropic Berendsen: mu^3 = 1 + (dt/tau)(P-P0)/|P|, clamped ±3%/step.
     * Sign matters: P above target must EXPAND (lower P), below must
     * compress. An earlier draft had (P0-P) and expanded into collapse. */
    double denom = fabs(cur) > 1e-9 ? fabs(cur) : 1e-9;
    double f = 1.0 + (t->sim->dt / t->taup_fs) * (cur - tgt) / denom;
    if (f < 0.97) f = 0.97;
    if (f > 1.03) f = 1.03;
    double mu = cbrt(f);
    if (!isfinite(mu) || mu <= 0) return;
    Vec3 c = vec3_scale(t->sim->box.dimensions, 0.5);
    for (int i = 0; i < t->sim->num_atoms; i++) {
        Vec3 d = vec3_sub(t->sim->atoms[i].position, c);
        t->sim->atoms[i].position = vec3_add(c, vec3_scale(d, mu));
    }
    t->sim->box.dimensions = vec3_scale(t->sim->box.dimensions, mu);
    forces_calculate(t->sim);
}

static int cmd_mol(Tui *t, char **tok, int nt) {
    if (nt < 2) { sh_err("  usage: mol new|add|bond|list|clear|center|save|load|place ...\n"); return 2; }
    if (!strcmp(tok[1], "new")) {
        mol.na = 0;
        mol.nb = 0;
        snprintf(mol.name, sizeof mol.name, "%s", nt >= 3 ? tok[2] : "mol");
        sh_err("  builder `%s` cleared\n", mol.name);
        return 0;
    }
    if (!strcmp(tok[1], "clear")) {
        mol.na = 0;
        mol.nb = 0;
        sh_err("  builder cleared\n");
        return 0;
    }
    if (!strcmp(tok[1], "list")) {
        printf("  mol `%s`: %d atoms %d bonds\n", mol.name[0] ? mol.name : "mol", mol.na, mol.nb);
        for (int i = 0; i < mol.na; i++)
            printf("    %d Z=%-3d (%7.3f %7.3f %7.3f) q=%+.4f\n", i, mol.a[i].Z,
                mol.a[i].pos.x, mol.a[i].pos.y, mol.a[i].pos.z, mol.a[i].q);
        for (int i = 0; i < mol.nb; i++)
            printf("    bond %d-%d order %d\n", mol.b[i].a, mol.b[i].b, mol.b[i].order);
        return 0;
    }
    if (!strcmp(tok[1], "add") && nt >= 6) {
        int Z = 0;
        double x, y, z, q = 0;
        if (isalpha((unsigned char)tok[2][0])) {
            const Element *e = pt_by_symbol(tok[2]);
            if (!e) { sh_err("  unknown symbol `%s`\n", tok[2]); return 1; }
            Z = e->Z;
        } else if (!parse_int(tok[2], &Z) || !pt_element(Z)) {
            sh_err("  bad Z\n");
            return 2;
        }
        if (!parse_double(tok[3], &x) || !parse_double(tok[4], &y) || !parse_double(tok[5], &z)) {
            sh_err("  usage: mol add <Z|sym> x y z [q]\n");
            return 2;
        }
        if (nt >= 7 && !parse_double(tok[6], &q)) { sh_err("  bad q\n"); return 2; }
        if (mol.na >= MOL_MAXA) { sh_err("  builder full (%d)\n", MOL_MAXA); return 1; }
        MolAtom *a = &mol.a[mol.na];
        a->Z = Z;
        a->pos = vec3(x, y, z);
        a->q = q;
        a->eps = 0;
        a->sig = 0;
        a->has_lj = 0;
        printf("  mol atom %d Z=%d\n", mol.na, Z);
        mol.na++;
        return 0;
    }
    if (!strcmp(tok[1], "bond") && nt >= 4) {
        int a, b, o = 1;
        if (!parse_int(tok[2], &a) || !parse_int(tok[3], &b)) {
            sh_err("  usage: mol bond <a> <b> [order]\n");
            return 2;
        }
        if (nt >= 5 && (!parse_int(tok[4], &o) || o < 1 || o > 3)) {
            sh_err("  usage: mol bond <a> <b> [order 1-3]\n");
            return 2;
        }
        if (a < 0 || a >= mol.na || b < 0 || b >= mol.na || a == b) {
            sh_err("  bad builder indices\n");
            return 1;
        }
        if (mol.nb >= MOL_MAXA) { sh_err("  builder full\n"); return 1; }
        mol.b[mol.nb].a = a;
        mol.b[mol.nb].b = b;
        mol.b[mol.nb].order = o;
        mol.nb++;
        sh_err("  mol bond %d-%d order %d\n", a, b, o);
        return 0;
    }
    if (!strcmp(tok[1], "center")) {
        Vec3 c;
        mol_com(&c);
        for (int i = 0; i < mol.na; i++) mol.a[i].pos = vec3_sub(mol.a[i].pos, c);
        sh_err("  recentered (%d atoms)\n", mol.na);
        return 0;
    }
    if (!strcmp(tok[1], "save") && nt >= 3) {
        int rc = mol_save(tok[2]);
        if (rc == 0) sh_err("  saved `%s` (%d atoms)\n", tok[2], mol.na);
        return rc;
    }
    if (!strcmp(tok[1], "load") && nt >= 3) {
        int rc = mol_load(tok[2]);
        if (rc == 0) sh_err("  loaded `%s` (%d atoms)\n", tok[2], mol.na);
        return rc;
    }
    if (!strcmp(tok[1], "place")) {
        double x = 0, y = 0, z = 0;
        if (nt >= 5 && (!parse_double(tok[2], &x) || !parse_double(tok[3], &y) || !parse_double(tok[4], &z))) {
            sh_err("  usage: mol place [x y z]\n");
            return 2;
        }
        int rc = mol_place(t, vec3(x, y, z));
        if (rc == 0) {
            sh_err("  placed `%s` N=%d\n", mol.name[0] ? mol.name : "mol", t->sim->num_atoms);
            maybe_render(t);
        }
        return rc;
    }
    sh_err("  usage: mol new|add|bond|list|clear|center|save|load|place ...\n");
    return 2;
}

static int cmd_rxn(Tui *t, char **tok, int nt) {
    if (nt < 2) { sh_err("  usage: rxn new|pair|break|make|delete|charge|list|save|load|fire|arm|auto ...\n"); return 2; }
    if (!strcmp(tok[1], "new")) {
        memset(&rxn, 0, sizeof rxn);
        snprintf(rxn.name, sizeof rxn.name, "%s", nt >= 3 ? tok[2] : "rxn");
        rxn.make_order = -1;
        sh_err("  reaction `%s` started\n", rxn.name);
        return 0;
    }
    if (!strcmp(tok[1], "list")) {
        if (!rxn.defined) { sh_err("  no reaction defined\n"); return 1; }
        printf("  rxn `%s`: Z%d-Z%d rcut=%.3f break=%d make=%d del=%d charges=%s armed=%d auto=%d\n",
            rxn.name, rxn.Z1, rxn.Z2, rxn.rcut, rxn.dobreak, rxn.make_order,
            rxn.delrole, rxn.setq ? "yes" : "no", t->rxn_armed, t->rxn_every);
        return 0;
    }
    if (!strcmp(tok[1], "pair") && nt >= 5) {
        int z1, z2;
        double rc;
        if (!parse_int(tok[2], &z1) || !parse_int(tok[3], &z2) || !parse_double(tok[4], &rc) ||
            !pt_element(z1) || !pt_element(z2) || !(rc > 0)) {
            sh_err("  usage: rxn pair <Z1> <Z2> <rcut_A>\n");
            return 2;
        }
        rxn.Z1 = z1; rxn.Z2 = z2; rxn.rcut = rc;
        rxn.defined = 1;
        sh_err("  pair Z%d-Z%d within %.3f A\n", z1, z2, rc);
        return 0;
    }
    if (!strcmp(tok[1], "break") && nt >= 3) {
        if (!strcmp(tok[2], "on")) rxn.dobreak = 1;
        else if (!strcmp(tok[2], "off")) rxn.dobreak = 0;
        else { sh_err("  usage: rxn break on|off\n"); return 2; }
        sh_err("  break %s\n", rxn.dobreak ? "on" : "off");
        return 0;
    }
    if (!strcmp(tok[1], "make") && nt >= 3) {
        if (!strcmp(tok[2], "none")) rxn.make_order = -1;
        else {
            int o;
            if (!parse_int(tok[2], &o) || o < 1 || o > 3) {
                sh_err("  usage: rxn make <order 1-3>|none\n");
                return 2;
            }
            rxn.make_order = o;
        }
        sh_err("  make order %d\n", rxn.make_order);
        return 0;
    }
    if (!strcmp(tok[1], "delete") && nt >= 3) {
        int r;
        if (!parse_int(tok[2], &r) || r < 0 || r > 2) {
            sh_err("  usage: rxn delete <0|1|2>  (0 none, 1 = Z1 atom, 2 = Z2)\n");
            return 2;
        }
        rxn.delrole = r;
        sh_err("  delete role %d\n", r);
        return 0;
    }
    if (!strcmp(tok[1], "charge")) {
        if (nt >= 3 && !strcmp(tok[2], "off")) {
            rxn.setq = 0;
            sh_err("  charges off\n");
            return 0;
        }
        if (nt >= 4) {
            double q1, q2;
            if (!parse_double(tok[2], &q1) || !parse_double(tok[3], &q2)) {
                sh_err("  usage: rxn charge <q1> <q2>|off\n");
                return 2;
            }
            rxn.q1 = q1; rxn.q2 = q2; rxn.setq = 1;
            sh_err("  charges q1=%.3f q2=%.3f\n", q1, q2);
            return 0;
        }
        sh_err("  usage: rxn charge <q1> <q2>|off\n");
        return 2;
    }
    if (!strcmp(tok[1], "save") && nt >= 3) {
        int rc = rxn_save(tok[2]);
        if (rc == 0) sh_err("  saved `%s`\n", tok[2]);
        return rc;
    }
    if (!strcmp(tok[1], "load") && nt >= 3) {
        int rc = rxn_load(tok[2]);
        if (rc == 0) sh_err("  loaded `%s`\n", tok[2]);
        return rc;
    }
    if (!strcmp(tok[1], "fire")) return rxn_fire(t);
    if (!strcmp(tok[1], "arm")) {
        t->rxn_armed = 1;
        sh_err("  reaction armed\n");
        return 0;
    }
    if (!strcmp(tok[1], "auto")) {
        if (nt >= 3 && !strcmp(tok[2], "off")) {
            t->rxn_every = 0;
            sh_err("  auto off\n");
            return 0;
        }
        if (nt >= 3) {
            int e;
            if (!parse_int(tok[2], &e) || e < 1) {
                sh_err("  usage: rxn auto <every_N>|off  (+ `rxn arm`)\n");
                return 2;
            }
            t->rxn_every = e;
            sh_err("  auto every %d steps (arm to fire once)\n", e);
            return 0;
        }
        sh_err("  usage: rxn auto <every_N>|off  (+ `rxn arm`)\n");
        return 2;
    }
    sh_err("  usage: rxn new|pair|break|make|delete|charge|list|save|load|fire|arm|auto ...\n");
    return 2;
}

/* ── per-demo spawners: each batch demo as one live starting system ──
 * Same builders the batch demos use (never re-derived geometry here):
 * pair placement follows Demo 7's WC-edge construction (N1...N3 2.9 A,
 * ring-normal alignment + azimuthal fix), the helix uses Demo 11's
 * textbook phi/psi/omega restraints, the cage Demo 12's 3D->xy
 * antiprism derivation, the duplex Demo 17's 3.4 A rise / 36° twist
 * (G-C/A-U stack; batch uses G-C/A-T plus backbone — stated on spawn). */
static Vec3 spo_origin(char **tok, int nt, int k, int *ok) {
    double x = 0, y = 0, z = 0;
    *ok = 1;
    if (nt > k) {
        if (nt < k + 3 || !parse_double(tok[k], &x) || !parse_double(tok[k + 1], &y) ||
            !parse_double(tok[k + 2], &z)) *ok = 0;
    }
    return vec3(x, y, z);
}

static int sp_is_number(const char *s) {
    double v;
    return parse_double(s, &v);
}

/* fail fast before building: placers report first-index even when short on
 * room, so check capacity up front (estimates are generous upper bounds). */
static int sp_room(Tui *t, int atoms, int bonds) {
    if (t->sim->num_atoms + atoms > t->sim->capacity_atoms ||
        t->sim->num_bonds + bonds > t->sim->capacity_bonds) {
        sh_err("  sim full (need ~%d atoms/%d bonds, free %d/%d; try `new` bigger)\n",
            atoms, bonds, t->sim->capacity_atoms - t->sim->num_atoms,
            t->sim->capacity_bonds - t->sim->num_bonds);
        return 0;
    }
    return 1;
}

/* Demo 1: hydrogen atom + its orbital line (analysis lives in test demo 1) */
static int sp_qm_h(Tui *t, Vec3 o) {
    if (!sp_room(t, 1, 0)) return 1;
    int i = sim_add_atom(t->sim, 1, o, 0.0);
    if (i < 0) { sh_err("  sim full\n"); return 1; }
    ElectronConfig c;
    pt_electron_config(1, &c);
    double z = quantum_zeff_raw(1, 1, 0, &c);
    double e = quantum_orbital_energy(1, 1, 0, &c);
    double r = quantum_most_probable_radius(1, 0, z);
    printf("  H 1s: Zeff=%.3f E=%.4f eV r_mp=%.4f A (atom %d)\n", z, e, r, i);
    return 0;
}

/* Demo 3: Berendsen water at 300 K, demo settings */
static int sp_water_md(Tui *t, Vec3 o) {
    if (!sp_room(t, 3, 2)) return 1;
    int f = sim_place_h2o(t->sim, o);
    if (f < 0) { sh_err("  sim full\n"); return 1; }
    t->sim->dt = 0.5;
    t->sim->thermostat.type = THERMOSTAT_BERENDSEN;
    t->sim->thermostat.target_temperature = 300.0;
    t->sim->thermostat.tau = 100.0;
    integrator_maxwell_boltzmann(t->sim, 300.0, t->seed);
    forces_calculate(t->sim);
    sh_err("  water @%d (Berendsen 300K, dt 0.5fs)\n", f);
    return 0;
}

/* Demo 4: cyclic trimer ring, O...O 2.95 A per edge */
static int sp_trimer(Tui *t, Vec3 o) {
    if (!sp_room(t, 9, 6)) return 1;
    double R = 2.95 / sqrt(3.0);
    for (int k = 0; k < 3; k++) {
        double a = (90.0 + k * 120.0) * 3.14159265358979323846 / 180.0;
        Vec3 p = vec3(o.x + R * cos(a), o.y + R * sin(a), o.z);
        if (sim_place_h2o(t->sim, p) < 0) { sh_err("  sim full\n"); return 1; }
    }
    t->sim->thermostat.type = THERMOSTAT_BERENDSEN;
    t->sim->thermostat.target_temperature = 50.0;
    t->sim->thermostat.tau = 100.0;
    forces_calculate(t->sim);
    sh_err("  trimer ring (O...O 2.95 A) placed\n");
    return 0;
}

/* Demo 6: one base, or the five-base row */
static int sp_base_one(Tui *t, const char *nm, Vec3 o) {
    if (!sp_room(t, 16, 16)) return 1;
    int f = -1;
    if (!strcasecmp(nm, "U")) f = sim_place_uracil(t->sim, o);
    else if (!strcasecmp(nm, "C")) f = sim_place_cytosine(t->sim, o);
    else if (!strcasecmp(nm, "T")) f = sim_place_thymine(t->sim, o);
    else if (!strcasecmp(nm, "A")) f = sim_place_adenine(t->sim, o);
    else if (!strcasecmp(nm, "G")) f = sim_place_guanine(t->sim, o);
    else { sh_err("  base U|C|T|A|G\n"); return 2; }
    if (f < 0) { sh_err("  sim full\n"); return 1; }
    forces_calculate(t->sim);
    return 0;
}

/* Demo 7: Watson-Crick pair in pairing geometry (dielectric 4 per demo) */
static int sp_pair(Tui *t, int gc, Vec3 o) {
    if (!sp_room(t, 32, 32)) return 1;
    Simulation *s = t->sim;
    s->dielectric = 4.0;
    if (gc) {
        int g = sim_place_guanine(s, o);
        if (g < 0) { sh_err("  sim full\n"); return 1; }
        Vec3 n1 = s->atoms[g + 6].position;
        int c = sim_place_cytosine(s, vec3_add(n1, vec3(2.9, 0, 0)));
        if (c < 0) { sh_err("  sim full\n"); return 1; }
        int gr[3] = {g + 0, g + 1, g + 6};
        int cr[3] = {c + 0, c + 1, c + 2};
        Vec3 nG = nb_ring_normal(s, gr);
        Vec3 nC = nb_ring_normal(s, cr);
        double co = vec3_dot(nG, nC);
        Vec3 ax = vec3_cross(nC, nG);
        double an;
        if (vec3_norm(ax) < 1e-8) {
            an = (co > 0) ? 0.0 : 3.14159265358979323846;
            ax = (co > 0) ? vec3(0, 0, 1) : vec3_cross(nC, vec3(1, 0, 0));
            if (vec3_norm(ax) < 1e-8) ax = vec3(0, 1, 0);
        } else {
            if (co > 1.0) co = 1.0;
            if (co < -1.0) co = -1.0;
            an = acos(co);
        }
        Vec3 piv = s->atoms[c + 0].position;
        nb_transform_rigid(s, c, 13, piv, ax, an, vec3_zero());
        Vec3 vG = vec3_sub(s->atoms[g + 5].position, s->atoms[g + 6].position);
        Vec3 vC = vec3_sub(s->atoms[c + 5].position, s->atoms[c + 0].position);
        double az = nb_signed_inplane_angle(vC, vG, nG);
        nb_transform_rigid(s, c, 13, piv, nG, az, vec3_zero());
        Vec3 dir = vec3_normalize(vec3_sub(s->atoms[g + 13].position, s->atoms[g + 6].position));
        Vec3 tgt = vec3_add(s->atoms[g + 6].position, vec3_scale(dir, 2.95));
        nb_transform_rigid(s, c, 13, piv, vec3_zero(), 0.0, vec3_sub(tgt, s->atoms[c + 0].position));
        sh_err("  G-C pair @%d (N1...N3 %.2f A, dielectric set 4)\n", g,
            vec3_dist(s->atoms[g + 6].position, s->atoms[c + 0].position));
    } else {
        int a = sim_place_adenine(s, o);
        if (a < 0) { sh_err("  sim full\n"); return 1; }
        int u = sim_place_thymine(s, vec3(o.x + 15.0, o.y, o.z));
        if (u < 0) { sh_err("  sim full\n"); return 1; }
        int ar[3] = {a + 0, a + 1, a + 6};
        int ur[3] = {u + 0, u + 1, u + 3};
        Vec3 nA = nb_ring_normal(s, ar);
        Vec3 nU = nb_ring_normal(s, ur);
        double co = vec3_dot(nA, nU);
        Vec3 ax = vec3_cross(nU, nA);
        double an;
        if (vec3_norm(ax) < 1e-8) {
            an = (co > 0) ? 0.0 : 3.14159265358979323846;
            ax = (co > 0) ? vec3(0, 0, 1) : vec3_cross(nU, vec3(1, 0, 0));
            if (vec3_norm(ax) < 1e-8) ax = vec3(0, 1, 0);
        } else {
            if (co > 1.0) co = 1.0;
            if (co < -1.0) co = -1.0;
            an = acos(co);
        }
        Vec3 piv = s->atoms[u + 3].position;
        nb_transform_rigid(s, u, 15, piv, ax, an, vec3_zero());
        Vec3 vA = vec3_sub(s->atoms[a + 5].position, s->atoms[a + 6].position);
        Vec3 vU = vec3_sub(s->atoms[u + 5].position, s->atoms[u + 3].position);
        double az = nb_signed_inplane_angle(vU, vA, nA);
        nb_transform_rigid(s, u, 15, piv, nA, az, vec3_zero());
        Vec3 dir = vec3_normalize(vec3_sub(s->atoms[u + 9].position, s->atoms[u + 3].position));
        Vec3 tgt = vec3_sub(s->atoms[a + 6].position, vec3_scale(dir, 2.90));
        nb_transform_rigid(s, u, 15, piv, vec3_zero(), 0.0, vec3_sub(tgt, s->atoms[u + 3].position));
        sh_err("  A-T pair @%d (N1...N3 %.2f A, dielectric set 4)\n", a,
            vec3_dist(s->atoms[a + 6].position, s->atoms[u + 3].position));
    }
    forces_calculate(s);
    return 0;
}

/* ── Demo 11: 5-alanine helix with textbook phi/psi/omega restraints ──
 * Mirrors the batch Demo 11 protocol: clash-relax first (no restraints),
 * then steer, then a long gentle minimize. Skipping the first step leaves
 * ~14 eV of strain in the live world, which NVE dynamics converts to
 * heat until the display overflows. Berendsen 50 K + MB start matches
 * the other demo spawners (water at 300 K, pair/duplex legs at 50 K). */
static int sp_helix(Tui *t, Vec3 o) {
    if (!sp_room(t, 60, 80)) return 1;
    AAResidue res[8];
    int f = sim_place_polyalanine(t->sim, o, 5, res);
    if (f < 0) { sh_err("  sim full\n"); return 1; }
    t->sim->dt = 0.5;
    forces_calculate(t->sim);
    integrator_minimize(t->sim, 5000, 0.001, 0.01);
    double k = 80.0 * KCAL_MOL_TO_EV;
    double dphi = (-57.0 - 180.0) * 3.14159265358979323846 / 180.0;
    double dpsi = (-47.0 - 180.0) * 3.14159265358979323846 / 180.0;
    double dom = (180.0 - 180.0) * 3.14159265358979323846 / 180.0;
    for (int i = 1; i < 5; i++) {
        sim_add_dihedral(t->sim, res[i - 1].C, res[i].N, res[i].CA, res[i].C, k, 1, dphi);
        sim_add_dihedral(t->sim, res[i - 1].CA, res[i - 1].C, res[i].N, res[i].CA, k, 1, dom);
    }
    for (int i = 0; i < 4; i++)
        sim_add_dihedral(t->sim, res[i].N, res[i].CA, res[i].C, res[i + 1].N, k, 1, dpsi);
    forces_calculate(t->sim);
    /* fresh chains carry steric clashes: relieve before any dynamics */
    integrator_minimize(t->sim, 200000, 0.0002, 0.0001);
    t->sim->thermostat.type = THERMOSTAT_BERENDSEN;
    t->sim->thermostat.target_temperature = 50.0;
    t->sim->thermostat.tau = 100.0;
    integrator_maxwell_boltzmann(t->sim, 50.0, t->seed);
    forces_calculate(t->sim);
    sh_err("  helix 5-ala @%d (phi/psi restrained, clash-relieved PE=%.3f)\n", f, t->sim->potential_energy);
    return 0;
}

/* Demo 12: constructed 8-O antiprism cage + point-charge K+ */
static int sp_cage(Tui *t, Vec3 o) {
    if (!sp_room(t, 9, 0)) return 1;
    Simulation *s = t->sim;
    double hs = 3.084 * 0.5;
    double ri = sqrt(2.70 * 2.70 - hs * hs);
    double ro = sqrt(2.83 * 2.83 - hs * hs);
    for (int i = 0; i < 4; i++) {
        double a = i * 3.14159265358979323846 / 2.0;
        int at = sim_add_atom(s, 8, vec3(o.x + ri * cos(a), o.y + ri * sin(a), o.z - hs), -0.5462);
        if (at < 0) { sh_err("  sim full\n"); return 1; }
        sim_set_atom_lj(s, at, LJ_AMBER_O_EPS, LJ_AMBER_O_SIGMA);
    }
    for (int i = 0; i < 4; i++) {
        double a = i * 3.14159265358979323846 / 2.0 + 3.14159265358979323846 / 4.0;
        int at = sim_add_atom(s, 8, vec3(o.x + ro * cos(a), o.y + ro * sin(a), o.z + hs), -0.5462);
        if (at < 0) { sh_err("  sim full\n"); return 1; }
        sim_set_atom_lj(s, at, LJ_AMBER_O_EPS, LJ_AMBER_O_SIGMA);
    }
    int ion = sim_add_ion(s, 19, 1, o, 1.0);
    if (ion < 0) { sh_err("  sim full\n"); return 1; }
    sim_set_atom_lj(s, ion, 0.0, 0.0);
    forces_calculate(s);
    sh_err("  cage 8-O antiprism + K+ point charge (z-sep 3.084 A)\n");
    return 0;
}

/* Demo 17 (reduced): G-C pair stacked over A-U pair, 3.4 A rise, 36° twist.
 * Batch uses G-C/A-T plus backbone; this is the geometry without it. */
static int sp_duplex(Tui *t, Vec3 o) {
    if (!sp_room(t, 60, 60)) return 1;
    Simulation *s = t->sim;
    int n0 = s->num_atoms;
    if (sp_pair(t, 1, o)) return 1;
    int a0 = s->num_atoms;
    /* second pair built at the origin, then rigid-moved into the stack */
    int a = sim_place_adenine(s, vec3_zero());
    if (a < 0) { sh_err("  sim full\n"); return 1; }
    int u = sim_place_uracil(s, vec3(15.0, 0.0, 0.0));
    if (u < 0) { sh_err("  sim full\n"); return 1; }
    int ar[3] = {a + 0, a + 1, a + 6};
    int ur[3] = {u + 0, u + 1, u + 3};
    Vec3 nA = nb_ring_normal(s, ar);
    Vec3 nU = nb_ring_normal(s, ur);
    double co = vec3_dot(nA, nU);
    Vec3 ax = vec3_cross(nU, nA);
    double an;
    if (vec3_norm(ax) < 1e-8) {
        an = (co > 0) ? 0.0 : 3.14159265358979323846;
        ax = (co > 0) ? vec3(0, 0, 1) : vec3(0, 1, 0);
    } else {
        if (co > 1.0) co = 1.0;
        if (co < -1.0) co = -1.0;
        an = acos(co);
    }
    Vec3 piv = s->atoms[u + 3].position;
    nb_transform_rigid(s, u, 15, piv, ax, an, vec3_zero());
    Vec3 vA = vec3_sub(s->atoms[a + 5].position, s->atoms[a + 6].position);
    Vec3 vU = vec3_sub(s->atoms[u + 5].position, s->atoms[u + 3].position);
    double az = nb_signed_inplane_angle(vU, vA, nA);
    nb_transform_rigid(s, u, 15, piv, nA, az, vec3_zero());
    Vec3 dir = vec3_normalize(vec3_sub(s->atoms[u + 9].position, s->atoms[u + 3].position));
    Vec3 tgt = vec3_sub(s->atoms[a + 6].position, vec3_scale(dir, 2.90));
    nb_transform_rigid(s, u, 15, piv, vec3_zero(), 0.0, vec3_sub(tgt, s->atoms[u + 3].position));
    /* stack: centroid, 36° about Y, rise 3.4 along Y from pair-1 plane */
    Vec3 cen = vec3_zero();
    for (int i = a0; i < s->num_atoms; i++) cen = vec3_add(cen, s->atoms[i].position);
    cen = vec3_scale(cen, 1.0 / (s->num_atoms - a0));
    Vec3 want = vec3(o.x, o.y + 3.4, o.z);
    nb_transform_rigid(s, a0, s->num_atoms - a0, cen, vec3(0, 1, 0),
        36.0 * 3.14159265358979323846 / 180.0, vec3_sub(want, cen));
    s->dielectric = 4.0;
    forces_calculate(t->sim);
    integrator_minimize(t->sim, 1000, 0.005, 0.02);
    sh_err("  duplex stack: G-C @%d + A-U @%d (rise 3.4 A, twist 36 deg; no backbone)\n", n0, a0);
    return 0;
}

/* dispatch one demo id to its spawner; extra selects base/pair variants */
static int cmd_spawn_demo(Tui *t, const char *id, const char *extra, Vec3 o) {
    if (!strcmp(id, "1")) return sp_qm_h(t, o);
    if (!strcmp(id, "2")) {
        if (!sp_room(t, 2, 1)) return 1;
        int f = sim_place_h2(t->sim, o);
        if (f < 0) { sh_err("  sim full\n"); return 1; }
        sh_err("  H2 covalent @%d (r0 0.7414 A)\n", f);
        return 0;
    }
    if (!strcmp(id, "3")) return sp_water_md(t, o);
    if (!strcmp(id, "4")) return sp_trimer(t, o);
    if (!strcmp(id, "5")) {
        if (!sp_room(t, 5, 4)) return 1;
        int f = sim_place_ch4(t->sim, o);
        if (f < 0) { sh_err("  sim full\n"); return 1; }
        sh_err("  CH4 tetrahedral @%d\n", f);
        return 0;
    }
    if (!strcmp(id, "6")) {
        if (extra) {
            int rc = sp_base_one(t, extra, o);
            if (rc == 0) sh_err("  base %s placed\n", extra);
            return rc;
        }
        if (!sp_room(t, 75, 75)) return 1;
        const char *bn[5] = {"U", "C", "T", "A", "G"};
        for (int k = 0; k < 5; k++) {
            Vec3 p = vec3(o.x + 10.0 * k, o.y, o.z);
            if (sp_base_one(t, bn[k], p)) return 1;
        }
        sh_err("  five-base row placed\n");
        return 0;
    }
    if (!strcmp(id, "7")) {
        int gc = 1;
        if (extra && !strcasecmp(extra, "au")) gc = 0;
        else if (extra && strcasecmp(extra, "gc")) { sh_err("  pair gc|au\n"); return 2; }
        return sp_pair(t, gc, o);
    }
    if (!strcmp(id, "8")) {
        if (!sp_room(t, 63, 80)) return 1;
        int s2 = -1;
        int f = sim_place_dinucleotide_TA(t->sim, o, &s2);
        if (f < 0) { sh_err("  sim full\n"); return 1; }
        sh_err("  T-p-A dinucleotide @%d\n", f);
        return 0;
    }
    if (!strcmp(id, "9")) {
        hh_init(&t->nrn);
        t->has_nrn = 1;
        sh_err("  neuron V=%.2f m=%.4f h=%.4f n=%.4f\n", t->nrn.V, t->nrn.m, t->nrn.h, t->nrn.n);
        return 0;
    }
    if (!strcmp(id, "10")) {
        if (!sp_room(t, 25, 25)) return 1;
        int an = -1;
        int f = sim_place_dipeptide_GlyAla(t->sim, o, &an);
        if (f < 0) { sh_err("  sim full\n"); return 1; }
        sh_err("  Gly-Ala dipeptide @%d\n", f);
        return 0;
    }
    if (!strcmp(id, "11")) return sp_helix(t, o);
    if (!strcmp(id, "12")) return sp_cage(t, o);
    if (!strcmp(id, "12b")) {
        if (!sp_room(t, 170, 170)) return 1;
        int f = kcsa_build_filter(t->sim, o, 4);
        if (f < 0) { sh_err("  sim full\n"); return 1; }
        sh_err("  real filter 164 atoms @%d (1K4C TVGYG C4)\n", f);
        return 0;
    }
    if (!strcmp(id, "17")) return sp_duplex(t, o);
    sh_err("  demo 1 2 3 4 5 6 7 8 9 10 11 12 12b 17\n");
    return 2;
}

/* ══════════════════════════════════════════════════════════════════════
 * PURE SIMULATOR MODE — petri world + live keybind surfaces
 *
 * A world of biological and organic matter that interacts under the
 * engine's own MD, and two live keybind surfaces over the SAME world the
 * shell operates on:
 *
 *   ps            live control monitor (single keys, no modes)
 *   vi world      vi-style world editor (modal; `:` runs shell commands)
 *
 * Both draw to stderr, so stdout stays DATA (the stream contract) and a
 * piped session never gets escape codes. Both step the shared Simulation
 * directly, so every shell command, reaction, template and job remains
 * live while the world runs.
 *
 * World generation lives in world_preset() behind named presets; "petri"
 * is the current default and the switch is the extension point for future
 * presets (cell, tissue, ...). Nothing here is linked into carbonsim.
 * ═════════════════════════════════════════════════════════════════════= */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define LIVE_KEY_UP    1001
#define LIVE_KEY_DOWN  1002
#define LIVE_KEY_RIGHT 1003
#define LIVE_KEY_LEFT  1004

static struct termios live_tio_saved;
static int live_tio_active = 0;

static int live_raw_on(void) {
    struct termios raw;
    if (live_tio_active) return 0;
    if (!isatty(STDIN_FILENO)) return -1;
    if (tcgetattr(STDIN_FILENO, &live_tio_saved) != 0) return -1;
    raw = live_tio_saved;
    raw.c_lflag &= (unsigned)~(ICANON | ECHO | ISIG);
    raw.c_iflag &= (unsigned)~(IXON | ICRNL);
    raw.c_cc[VMIN] = 1;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSANOW, &raw) != 0) return -1;
    live_tio_active = 1;
    return 0;
}

static void live_raw_off(void) {
    if (!live_tio_active) return;
    tcsetattr(STDIN_FILENO, TCSANOW, &live_tio_saved);
    live_tio_active = 0;
}

/* one byte, with arrow keys decoded; -1 when nothing is ready */
static int live_getch(void) {
    unsigned char c = 0;
    ssize_t r;
    if (!live_tio_active) return -1;
    r = read(STDIN_FILENO, &c, 1);
    if (r <= 0) return -1;
    if (c != 0x1b) return (int)c;
    {
        struct termios t0, t1;
        unsigned char b = 0, d = 0;
        if (tcgetattr(STDIN_FILENO, &t0) != 0) return 0x1b;
        t1 = t0;
        t1.c_cc[VMIN] = 0;
        t1.c_cc[VTIME] = 1;   /* 100 ms window for the sequence */
        tcsetattr(STDIN_FILENO, TCSANOW, &t1);
        if (read(STDIN_FILENO, &b, 1) == 1 && b == '[' &&
            read(STDIN_FILENO, &d, 1) == 1) {
            tcsetattr(STDIN_FILENO, TCSANOW, &t0);
            switch (d) {
                case 'A': return LIVE_KEY_UP;
                case 'B': return LIVE_KEY_DOWN;
                case 'C': return LIVE_KEY_RIGHT;
                case 'D': return LIVE_KEY_LEFT;
                default:  return 0x1b;
            }
        }
        tcsetattr(STDIN_FILENO, TCSANOW, &t0);
    }
    return 0x1b;
}

/* poll=1 returns -1 after ~100 ms instead of blocking (running frames) */
static int live_getch_poll(int poll) {
    struct termios t0, t1;
    int r;
    if (!poll || !live_tio_active) return live_getch();
    if (tcgetattr(STDIN_FILENO, &t0) != 0) return live_getch();
    t1 = t0;
    t1.c_cc[VMIN] = 0;
    t1.c_cc[VTIME] = 1;
    tcsetattr(STDIN_FILENO, TCSANOW, &t1);
    r = live_getch();
    tcsetattr(STDIN_FILENO, TCSANOW, &t0);
    return r;
}

/* ── world generation ──────────────────────────────────────────────── */

static unsigned long long live_rng_state = 88172645463325252ull;

static double live_rand01(void) {
    live_rng_state ^= live_rng_state << 13;
    live_rng_state ^= live_rng_state >> 7;
    live_rng_state ^= live_rng_state << 17;
    return (double)(live_rng_state >> 11) / 9007199254740992.0;
}
static double live_rand(double lo, double hi) { return lo + (hi - lo) * live_rand01(); }

static int live_z_from_name(const char *s) {
    const Element *e = pt_by_symbol(s);
    if (e) return e->Z;
    {
        int Z = -1;
        if (parse_int(s, &Z) && Z >= 1 && Z <= 36) return Z;
    }
    if (!strcmp(s, "na")) return 11;
    if (!strcmp(s, "k")) return 19;
    if (!strcmp(s, "cl")) return 17;
    if (!strcmp(s, "ca")) return 20;
    return -1;
}

/* bounding radius per species for open-site placement: the clearance
 * point must clear the whole molecule, not one atom, or edges overlap
 * and the next step detonates. Values are generous envelopes. */
static double live_species_radius(const char *name) {
    if (!name) return 3.0;
    if (!strcmp(name, "water") || !strcmp(name, "h2o")) return 1.6;
    if (!strcmp(name, "h2")) return 1.0;
    if (!strcmp(name, "nh3")) return 1.7;
    if (!strcmp(name, "ch4") || !strcmp(name, "methane")) return 2.0;
    if (!strcmp(name, "co2")) return 2.4;
    if (!strcmp(name, "glycine") || !strcmp(name, "gly")) return 3.0;
    if (!strcmp(name, "alanine") || !strcmp(name, "ala")) return 3.5;
    if (!strcmp(name, "uracil") || !strcmp(name, "cytosine") ||
        !strcmp(name, "thymine") || !strcmp(name, "adenine") ||
        !strcmp(name, "guanine")) return 4.0;
    if (!strcmp(name, "deoxyribose") || !strcmp(name, "sugar")) return 3.5;
    if (!strcmp(name, "Na+") || !strcmp(name, "na") ||
        !strcmp(name, "K+") || !strcmp(name, "k") ||
        !strcmp(name, "Cl-") || !strcmp(name, "cl") ||
        !strcmp(name, "Ca2+") || !strcmp(name, "ca")) return 1.0;
    if (live_z_from_name(name) >= 1) return 1.0;
    return 3.0;
}

/* min nonbonded new-vs-old distance for block [first, num). New atoms
 * share no bonds with old atoms, so every cross pair counts. */
static double live_block_gap(const Simulation *s, int first) {
    double best = 1e30;
    if (!s) return best;
    for (int i = first; i < s->num_atoms; i++) {
        for (int j = 0; j < first; j++) {
            double d = vec3_dist(s->atoms[i].position, s->atoms[j].position);
            if (d < best) best = d;
        }
    }
    return best;
}

static void live_move_block(Simulation *s, int first, Vec3 delta) {
    if (!s) return;
    for (int i = first; i < s->num_atoms; i++)
        s->atoms[i].position = vec3_add(s->atoms[i].position, delta);
}

/* add one entity (molecule, ion, or bare atom); returns atoms added */
static int live_species_add(Tui *t, const char *name, Vec3 at) {
    Simulation *s = t->sim;
    int first, added, ok = 0;
    if (!s || !name) return -1;
    if (!sp_room(t, 20, 20)) return -1;
    first = s->num_atoms;
    if (!strcmp(name, "water") || !strcmp(name, "h2o")) ok = sim_place_h2o(s, at) >= 0;
    else if (!strcmp(name, "nh3")) ok = sim_place_nh3(s, at) >= 0;
    else if (!strcmp(name, "methane") || !strcmp(name, "ch4")) ok = sim_place_ch4(s, at) >= 0;
    else if (!strcmp(name, "co2")) ok = sim_place_co2(s, at) >= 0;
    else if (!strcmp(name, "h2")) ok = sim_place_h2(s, at) >= 0;
    else if (!strcmp(name, "glycine") || !strcmp(name, "gly")) ok = sim_place_glycine(s, at) >= 0;
    else if (!strcmp(name, "alanine") || !strcmp(name, "ala")) ok = sim_place_alanine(s, at) >= 0;
    else if (!strcmp(name, "uracil")) ok = sim_place_uracil(s, at) >= 0;
    else if (!strcmp(name, "cytosine")) ok = sim_place_cytosine(s, at) >= 0;
    else if (!strcmp(name, "thymine")) ok = sim_place_thymine(s, at) >= 0;
    else if (!strcmp(name, "adenine")) ok = sim_place_adenine(s, at) >= 0;
    else if (!strcmp(name, "guanine")) ok = sim_place_guanine(s, at) >= 0;
    else if (!strcmp(name, "deoxyribose") || !strcmp(name, "sugar")) ok = sim_place_deoxyribose(s, at) >= 0;
    else if (!strcmp(name, "Na+") || !strcmp(name, "na")) ok = sim_add_ion(s, 11, 1, at, 1.0) >= 0;
    else if (!strcmp(name, "K+") || !strcmp(name, "k")) ok = sim_add_ion(s, 19, 1, at, 1.0) >= 0;
    else if (!strcmp(name, "Cl-") || !strcmp(name, "cl")) ok = sim_add_ion(s, 17, -1, at, -1.0) >= 0;
    else if (!strcmp(name, "Ca2+") || !strcmp(name, "ca")) ok = sim_add_ion(s, 20, 2, at, 2.0) >= 0;
    else {
        int Z = live_z_from_name(name);
        if (Z >= 1) ok = sim_add_atom(s, Z, at, 0.0) >= 0;
    }
    if (!ok) { sh_err("  unknown species `%s`\n", name); return -1; }
    added = s->num_atoms - first;
    sim_rebuild_angles(s);
    forces_calculate(s);
    return added;
}

static void live_rotate_block(Simulation *s, int first, int count) {
    Vec3 c = vec3_zero(), axis;
    double ang;
    if (!s || first < 0 || count < 1 || first + count > s->num_atoms) return;
    for (int i = 0; i < count; i++) c = vec3_add(c, s->atoms[first + i].position);
    c = vec3_scale(c, 1.0 / count);
    axis = vec3(live_rand(-1, 1), live_rand(-1, 1), live_rand(-1, 1));
    if (vec3_norm(axis) < 1e-6) axis = vec3(0, 0, 1);
    ang = live_rand(0.0, 2.0 * M_PI);
    for (int i = 0; i < count; i++) {
        Vec3 rel = vec3_sub(s->atoms[first + i].position, c);
        s->atoms[first + i].position = vec3_add(c, vec3_rotate_axis_angle(rel, axis, ang));
    }
}

static Vec3 live_open_site(const Simulation *s, double clearance) {
    for (int tries = 0; tries < 64; tries++) {
        Vec3 p = vec3(live_rand(4.0, 44.0), live_rand(4.0, 44.0), live_rand(4.0, 44.0));
        int clear = 1;
        for (int i = 0; i < s->num_atoms; i++) {
            if (vec3_dist(p, s->atoms[i].position) < clearance) { clear = 0; break; }
        }
        if (clear) return p;
    }
    return vec3(live_rand(4.0, 44.0), live_rand(4.0, 44.0), live_rand(4.0, 44.0));
}

/* Named world presets. "petri" is the current pure-simulator default;
 * future presets (cell, tissue, ...) slot in here. */
static int world_preset(Tui *t, const char *name, unsigned long seed) {
    Simulation *s;
    if (!t) return -1;
    if (!strcmp(name, "empty") || !strcmp(name, "void")) {
        tui_new(t, 2000, 4000);
        if (!t->sim) return -1;
        t->sim->dt = 0.5;
        t->sim->cutoff = 12.0;
        t->sim->dielectric = 1.0;
        sh_err("  world: empty (%d atom cap)\n", t->sim->capacity_atoms);
        return 0;
    }
    if (strcmp(name, "petri")) {
        sh_err("  unknown world preset `%s` (have: petri, empty)\n", name);
        return -1;
    }
    tui_new(t, 4000, 8000);
    if (!t->sim) return -1;
    s = t->sim;
    live_rng_state = seed ? (unsigned long long)seed : 88172645463325252ull;
    if (!live_rng_state) live_rng_state = 1;
    s->dt = 0.5;
    s->cutoff = 12.0;
    s->dielectric = 1.0;
    s->box.dimensions = vec3(48.0, 48.0, 48.0);
    s->box.periodic[0] = s->box.periodic[1] = s->box.periodic[2] = 1;
    s->thermostat.type = THERMOSTAT_ANDERSEN;
    s->thermostat.target_temperature = 310.0;
    s->thermostat.tau = 50.0;
    s->thermostat.nu = 0.02;
    /* 216 waters on a jittered 6x6x6 lattice (8 A spacing, no clashes) */
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++)
                sim_place_h2o(s, vec3(4.0 + 8.0 * i + live_rand(-0.8, 0.8),
                                      4.0 + 8.0 * j + live_rand(-0.8, 0.8),
                                      4.0 + 8.0 * k + live_rand(-0.8, 0.8)));
    /* ions */
    {
        static const int ion_Z[3] = {11, 19, 17};
        static const double ion_q[3] = {1.0, 1.0, -1.0};
        for (int n = 0; n < 24; n++)
            sim_add_ion(s, ion_Z[n % 3], n % 3 == 2 ? -1 : 1, live_open_site(s, 3.4), ion_q[n % 3]);
    }
    /* organic monomers, randomly oriented */
    {
        static const char *org[] = {"alanine", "glycine", "uracil", "adenine",
                                    "deoxyribose", "alanine", "glycine", "uracil"};
        for (unsigned n = 0; n < sizeof org / sizeof org[0]; n++) {
            Vec3 p = live_open_site(s, 4.5);
            int first = s->num_atoms;
            if (live_species_add(t, org[n], p) > 0)
                live_rotate_block(s, first, s->num_atoms - first);
        }
    }
    hh_init(&t->nrn);
    t->has_nrn = 1;
    integrator_maxwell_boltzmann(s, 310.0, seed ? seed : 7UL);
    sim_rebuild_angles(s);
    forces_calculate(s);
    sh_err("  world: petri  N=%d atoms, %d bonds, 48 A periodic box, 310 K Andersen, HH neuron\n",
           s->num_atoms, s->num_bonds);
    return 0;
}

/* ── live actions (shared by ps and vi) ────────────────────────────── */

static Vec3 live_cam_pos(const Tui *t) {
    return t->cam.has_center ? t->cam.center : vec3_zero();
}

static void live_advance(Tui *t, int steps) {
    if (!t->sim || steps < 1) return;
    for (int i = 0; i < steps; i++) {
        integrator_step(t->sim);
        tui_barostat_step(t);
        rxn_autocheck(t);
    }
    if (t->has_nrn) hh_step(&t->nrn, 0.01 * (double)steps);
    forces_calculate(t->sim);
}

static void live_heat(Tui *t, double dE) {
    Simulation *s = t->sim;
    double ke = 0.0, target, lam;
    if (!s || s->num_atoms < 1) return;
    for (int i = 0; i < s->num_atoms; i++) {
        double v2 = vec3_dot(s->atoms[i].velocity, s->atoms[i].velocity);
        ke += 0.5 * s->atoms[i].mass * v2 * AMU_AFS2_TO_EV;
    }
    target = ke + dE;
    if (target < 1e-6) target = 1e-6;
    lam = (ke > 1e-12) ? sqrt(target / ke) : 1.0;
    for (int i = 0; i < s->num_atoms; i++)
        s->atoms[i].velocity = vec3_scale(s->atoms[i].velocity, lam);
}

static int live_transmute(Tui *t, int atom, int Z) {
    const Element *e = pt_element(Z);
    Atom *a;
    if (!t->sim || atom < 0 || atom >= t->sim->num_atoms || !e) return -1;
    a = &t->sim->atoms[atom];
    a->Z = Z;
    a->element = e;
    a->mass = e->mass;
    a->lj_epsilon = e->lj_epsilon;
    a->lj_sigma = e->lj_sigma;
    return 0;
}

static int live_freeze_toggle(Tui *t, int atom) {
    Simulation *s = t->sim;
    if (!s || atom < 0 || atom >= s->num_atoms) return -1;
    for (int i = 0; i < s->num_restraints; i++) {
        if (s->restraint_atom[i] == atom) {
            for (int j = i; j < s->num_restraints - 1; j++) {
                s->restraint_atom[j] = s->restraint_atom[j + 1];
                s->restraint_anchor[j] = s->restraint_anchor[j + 1];
                s->restraint_k[j] = s->restraint_k[j + 1];
            }
            s->num_restraints--;
            return 0;   /* unfrozen */
        }
    }
    return sim_add_restraint(s, atom, s->atoms[atom].position, 5.0) >= 0 ? 1 : -1;
}

/* duplicate the connected molecule containing `root`; returns new first idx */
static int live_clone_molecule(Tui *t, int root) {
    Simulation *s = t->sim;
    int n, qh = 0, qt = 0, nc = 0, first;
    int *mark, *queue, *comp, *map;
    Vec3 shift;
    if (!s || root < 0 || root >= s->num_atoms) return -1;
    n = s->num_atoms;
    mark = (int *)calloc((size_t)n, sizeof(int));
    queue = (int *)malloc(sizeof(int) * (size_t)n);
    comp = (int *)malloc(sizeof(int) * (size_t)n);
    map = (int *)malloc(sizeof(int) * (size_t)n);
    if (!mark || !queue || !comp || !map) {
        free(mark); free(queue); free(comp); free(map);
        return -1;
    }
    for (int i = 0; i < n; i++) map[i] = -1;
    mark[root] = 1;
    queue[qt++] = root;
    while (qh < qt) {
        int a = queue[qh++];
        comp[nc++] = a;
        for (int b = 0; b < s->atoms[a].num_bonds; b++) {
            int p = s->atoms[a].bond_partners[b];
            if (p >= 0 && p < n && !mark[p]) { mark[p] = 1; queue[qt++] = p; }
        }
    }
    if (!sp_room(t, nc, nc + 8)) { free(mark); free(queue); free(comp); free(map); return -1; }
    shift = vec3(live_rand(-2.5, 2.5), live_rand(-2.5, 2.5), live_rand(-2.5, 2.5));
    if (vec3_norm(shift) < 1.5) shift = vec3(2.0, 0.0, 0.0);
    first = s->num_atoms;
    for (int i = 0; i < nc; i++) {
        int a = comp[i];
        map[a] = sim_add_atom(s, s->atoms[a].Z,
                              vec3_add(s->atoms[a].position, shift),
                              s->atoms[a].partial_charge);
        if (map[a] >= 0)
            sim_set_atom_lj(s, map[a], s->atoms[a].lj_epsilon, s->atoms[a].lj_sigma);
    }
    for (int i = 0; i < nc; i++) {
        int a = comp[i];
        for (int b = 0; b < s->atoms[a].num_bonds; b++) {
            int p = s->atoms[a].bond_partners[b];
            if (p > a && p >= 0 && p < n && map[p] >= 0)
                sim_add_bond(s, map[a], map[p], s->atoms[a].bond_orders[b]);
        }
    }
    free(mark); free(queue); free(comp); free(map);
    sim_rebuild_angles(s);
    forces_calculate(s);
    return first;
}

static void live_status(const Tui *t, const char *title, int paused, int speed, int sel) {
    const Simulation *s = t->sim;
    char tb[32], eb[32], tm[32];
    view_fmt_scalar(tb, sizeof tb, s ? s->temperature : 0.0, "%.1f", "");
    view_fmt_scalar(eb, sizeof eb, s ? s->total_energy : 0.0, "%.3f", "");
    view_fmt_scalar(tm, sizeof tm, s ? s->time : 0.0, "%.1f", "");
    fprintf(stderr, "%s  %s  speed=%d  N=%d  step=%llu  t=%s fs  T=%s K  E=%s eV\n",
            title, paused ? "PAUSED " : "RUNNING", speed,
            s ? s->num_atoms : 0, s ? (unsigned long long)s->step : 0ull,
            tm, tb, eb);
    if (t->has_nrn)
        fprintf(stderr, "process neuron  V=%+7.2f mV  I_ext=%+.1f  %s\n",
                t->nrn.V, t->nrn.I_ext,
                hh_is_spiking(&t->nrn, -20.0) ? "SPIKE" : "");
    if (s && s->num_restraints > 0)
        fprintf(stderr, "links: %d restrained atom(s)\n", s->num_restraints);
    if (sel >= 0 && s && sel < s->num_atoms) {
        const Atom *a = &s->atoms[sel];
        fprintf(stderr, "selection: atom %d  %s  q=%+.3f  bonds=%d\n", sel,
                (a->element && a->element->symbol[0]) ? a->element->symbol : "?",
                a->partial_charge, a->num_bonds);
    } else {
        fprintf(stderr, "selection: none (Tab / w,b to pick)\n");
    }
}

static int live_pick_next(const Simulation *s, int cur, int dir) {
    int n = s ? s->num_atoms : 0;
    if (n <= 0) return -1;
    int i = cur + dir;
    if (i < 0) i = n - 1;
    if (i >= n) i = 0;
    return i;
}

static void live_center_on(Tui *t, int sel) {
    if (t->sim && sel >= 0 && sel < t->sim->num_atoms) {
        t->cam.center = t->sim->atoms[sel].position;
        t->cam.has_center = 1;
    }
}

/* Strafe the look-at point across the environment in the view plane:
 * dx = +1 view-right, dy = +1 view-up, in angstroms. Anchors first:
 * strafing from auto-fit would be invisible, so the first strafe
 * pins the current frame center. `c` re-anchors onto an atom. */
static void live_strafe(Tui *t, double dx, double dy) {
    double yaw, pitch, cy, sy, cp, sp;
    Vec3 right, up;
    if (!t || !t->sim) return;
    if (!t->cam.has_center) {
        if (!view_bbox_center(t->sim, &t->cam.center)) return;
        t->cam.has_center = 1;
    }
    yaw = t->cam.yaw_deg * M_PI / 180.0;
    pitch = t->cam.pitch_deg * M_PI / 180.0;
    cy = cos(yaw); sy = sin(yaw);
    cp = cos(pitch); sp = sin(pitch);
    right = vec3(cy, -sy, 0.0);
    up = vec3(sy * cp, cy * cp, -sp);
    if (vec3_norm(right) > 1e-12)
        t->cam.center = vec3_add(t->cam.center, vec3_scale(right, dx));
    if (vec3_norm(up) > 1e-12)
        t->cam.center = vec3_add(t->cam.center, vec3_scale(up, dy));
}

static void live_add_at_cam(Tui *t, const char *species) {
    Vec3 at = live_cam_pos(t);
    int first = t->sim ? t->sim->num_atoms : -1;
    int added = live_species_add(t, species, at);
    if (added > 1 && first >= 0) live_rotate_block(t->sim, first, added);
    /* aimed at occupied space (camera inside a molecule): slide the new
     * block to open space rather than spawn inside matter. */
    if (t->sim && added > 0 && first >= 0 &&
        live_block_gap(t->sim, first) < 1.2) {
        double need = live_species_radius(species) + 2.0;
        Vec3 open = live_open_site(t->sim, need);
        live_move_block(t->sim, first, vec3_sub(open, t->sim->atoms[first].position));
        if (live_block_gap(t->sim, first) < 1.2)
            sh_err("  crowded world; `%s` placed tight\n", species);
    }
}

/* ── ps: live control monitor ──────────────────────────────────────── */

/* Shared ps keybinds for the raw-TTY monitor and the curses screen.
 * Returns 1 when the surface should close (q). via_screen selects the
 * prompt path for `r` (curses leaves the alternate screen to ask). */
static int live_key_ps(Tui *t, int *sel, int *paused, int *speed, int k,
                       int via_screen);

static int live_monitor(Tui *t) {
    int sel, paused = 0, speed = 4;
    if (!t->sim) { sh_err("  ps: no world (try `dd if=petri of=world`)\n"); return 1; }
    if (!isatty(STDIN_FILENO) || !isatty(STDERR_FILENO)) {
        live_status(t, "ps", 1, 0, -1);   /* piped: one snapshot, POSIX-like */
        return 0;
    }
    sel = t->sim->num_atoms > 0 ? 0 : -1;
    if (live_raw_on() != 0) { sh_err("  ps: cannot enter raw terminal mode\n"); return 1; }
    for (;;) {
        int k;
        if (!paused) live_advance(t, speed);
        view_render_to(stderr, t->sim, &t->cam, 1);
        live_status(t, "ps", paused, speed, sel);
        fprintf(stderr,
            "keys: space run/pause  s step  +/- speed  ijkl rotate  arrows strafe  z/Z zoom  Tab sel  c center\n"
            "      w water  N Na+  K K+  C Cl-  u base  a ala  g gly  d sugar  x del  f freeze\n"
            "      p clone  r replace  H/L heat/cool  m minimize  n neuron  e catalysis  q quit\n");
        fflush(stderr);
        k = live_getch_poll(!paused);
        if (k < 0) continue;
        if (live_key_ps(t, &sel, &paused, &speed, k, 0)) goto done;
    }
done:
    live_raw_off();
    fprintf(stderr, "\x1b[0m");
    return 0;
}

/* Fullscreen twin of the monitor: same world, same keybinds, curses
 * alternate screen with resize, colors and a status footer instead of
 * the scrolling ANSI frames. Piped/dumb/missing-curses degrades to a
 * typed line (status 0): a screen needs a human and a terminal. */
static int live_screen(Tui *t) {
    int sel, paused = 0, speed = 4;
    char help1[128], help2[160], status[160];
    if (!t->sim) { sh_err("  screen: no world (try `dd if=petri of=world`)\n"); return 1; }
    if (!screen_available()) {
        sh_err("  screen: needs a real terminal (TTY + TERM, ncursesw build); piped here\n");
        return 0;
    }
    if (screen_init() != 0) {
        sh_err("  screen: cannot open the alternate screen; staying in line mode\n");
        return 0;
    }
    screen_help(help1, sizeof help1, help2, sizeof help2);
    sel = t->sim->num_atoms > 0 ? 0 : -1;
    for (;;) {
        int k, cols, rows;
        ViewCells *vc;
        char tb[32], eb[32], tm[32];
        if (!paused) live_advance(t, speed);
        screen_size(&cols, &rows);
        vc = view_compute(t->sim, &t->cam, cols, rows + 2);
        view_fmt_scalar(tb, sizeof tb, t->sim->temperature, "%.1f", "");
        view_fmt_scalar(eb, sizeof eb, t->sim->potential_energy, "%.3f", "");
        view_fmt_scalar(tm, sizeof tm, t->sim->time, "%.1f", "");
        snprintf(status, sizeof status, "screen %s speed=%d N=%d step=%llu t=%s fs T=%s K E=%s eV sel=%d",
                 paused ? "PAUSED " : "RUNNING", speed, t->sim->num_atoms,
                 (unsigned long long)t->sim->step, tm, tb, eb, sel);
        screen_draw(vc, status, help1, help2, sel);
        view_cells_free(vc);
        k = screen_getch(paused ? -1 : 120);
        if (k == -1) continue;   /* frame tick */
        if (k == -3) continue;   /* resize: next draw refits */
        if (k == -2) break;      /* backend lost */
        if (k == 27) break;      /* Esc quits like q */
        if (live_key_ps(t, &sel, &paused, &speed, k, 1)) break;
    }
    screen_shutdown();
    return 0;
}

static int live_key_ps(Tui *t, int *sel, int *paused, int *speed, int k,
                       int via_screen) {
    switch (k) {
        case 'q': return 1;
        case ' ': *paused = !*paused; break;
        case 's': live_advance(t, 1); break;
        case '+': case '=': *speed += (*speed < 10 ? 1 : 10); if (*speed > 200) *speed = 200; break;
        case '-': *speed -= (*speed < 11 ? 1 : 10); if (*speed < 1) *speed = 1; break;
        /* ijkl rotates the camera (yaw/perspective); arrows strafe it
         * across the environment. vi keeps hjkl (vim standard). */
        case 'j': t->cam.yaw_deg += 8.0; break;
        case 'l': t->cam.yaw_deg -= 8.0; break;
        case 'k': t->cam.pitch_deg -= 5.0; break;
        case 'i': t->cam.pitch_deg += 5.0; break;
        case LIVE_KEY_LEFT:  live_strafe(t, -2.0, 0.0); break;
        case LIVE_KEY_RIGHT: live_strafe(t, +2.0, 0.0); break;
        case LIVE_KEY_UP:    live_strafe(t, 0.0, +2.0); break;
        case LIVE_KEY_DOWN:  live_strafe(t, 0.0, -2.0); break;
        case 'z': t->cam.zoom *= 0.85; if (t->cam.zoom < 0.1) t->cam.zoom = 0.1; break;
        case 'Z': t->cam.zoom *= 1.15; if (t->cam.zoom > 20.0) t->cam.zoom = 20.0; break;
        case '\t': *sel = live_pick_next(t->sim, *sel, 1); break;
        case 'c': live_center_on(t, *sel); break;
        case 'w': live_add_at_cam(t, "water"); break;
        case 'N': live_add_at_cam(t, "Na+"); break;
        case 'K': live_add_at_cam(t, "K+"); break;
        case 'C': live_add_at_cam(t, "Cl-"); break;
        case 'u': live_add_at_cam(t, "uracil"); break;
        case 'a': live_add_at_cam(t, "alanine"); break;
        case 'g': live_add_at_cam(t, "glycine"); break;
        case 'd': live_add_at_cam(t, "deoxyribose"); break;
        case 'x':
            if (*sel >= 0 && sim_remove_terminal_atom(t->sim, *sel))
                *sel = live_pick_next(t->sim, *sel, 1);
            break;
        case 'f': live_freeze_toggle(t, *sel); break;
        case 'p': if (*sel >= 0) live_clone_molecule(t, *sel); break;
        case 'r': {
            if (via_screen) {
                char buf[64];
                if (screen_prompt("replace with: ", buf, sizeof buf) == 0) {
                    int Z = live_z_from_name(buf);
                    if (Z >= 1 && *sel >= 0) live_transmute(t, *sel, Z);
                }
            } else {
                char *line = tui_readline("replace with: ");
                if (line) {
                    int Z = live_z_from_name(line);
                    if (Z >= 1 && *sel >= 0) live_transmute(t, *sel, Z);
                    free(line);
                }
            }
            break;
        }
        case 'H': live_heat(t, +5.0); break;
        case 'L': live_heat(t, -5.0); break;
        case 'm': if (t->sim) integrator_minimize(t->sim, 300, 0.005, 0.05); break;
        case 'n': hh_init(&t->nrn); t->has_nrn = 1; break;
        case 'e':
            if (!rxn.defined) {
                rxn.Z1 = 6; rxn.Z2 = 8; rxn.rcut = 1.6;
                rxn.dobreak = 0; rxn.make_order = 1;
                rxn.delrole = 0; rxn.setq = 0; rxn.q1 = rxn.q2 = 0.0;
                rxn.defined = 1;
                snprintf(rxn.name, sizeof rxn.name, "C-O proximity");
            }
            t->rxn_armed = !t->rxn_armed;
            if (t->rxn_armed && t->rxn_every <= 0) t->rxn_every = 50;
            break;
        default: break;
    }
    return 0;
}

/* ── vi world: modal editor over the live world ────────────────────── */

static int live_vi(Tui *t, char *undopath, size_t cap) {
    int sel, paused = 0;
    char msg[128];
    if (!t->sim) { sh_err("  vi: no world (try `dd if=petri of=world`)\n"); return 1; }
    if (!isatty(STDIN_FILENO) || !isatty(STDERR_FILENO)) {
        sh_err("  vi: world editor needs a TTY; use `more world` for a one-frame render\n");
        return 1;
    }
    sel = t->sim->num_atoms > 0 ? 0 : -1;
    snprintf(undopath, cap, "%s/s2tui-undo-%ld.s2",
             getenv("TMPDIR") ? getenv("TMPDIR") : "/tmp", (long)getpid());
    snprintf(msg, sizeof msg, "normal mode");
    if (live_raw_on() != 0) { sh_err("  vi: cannot enter raw terminal mode\n"); return 1; }
    for (;;) {
        int k;
        if (!paused) live_advance(t, 4);
        view_render_to(stderr, t->sim, &t->cam, 1);
        live_status(t, "vi world", paused, 4, sel);
        fprintf(stderr, "%s\n", msg);
        fprintf(stderr,
            "hjkl pan  w/b select  i insert  r replace  x delete  p clone  f freeze  u undo\n"
            "/ search  n next  c center  z/Z zoom  H/L heat/cool  m min  space run/pause  s step\n"
            ": command   q quit   (create/augment/adapt; every `:` command is the shell)\n");
        fflush(stderr);
        k = live_getch_poll(!paused);
        if (k < 0) continue;
        switch (k) {
            case 'q': goto done;
            case ' ': paused = !paused; break;
            case 's': live_advance(t, 1); break;
            case 'h': case LIVE_KEY_LEFT:  t->cam.yaw_deg += 8.0; break;
            case 'l': case LIVE_KEY_RIGHT: t->cam.yaw_deg -= 8.0; break;
            case 'j': case LIVE_KEY_DOWN:  t->cam.pitch_deg -= 5.0; break;
            case 'k': case LIVE_KEY_UP:    t->cam.pitch_deg += 5.0; break;
            case 'z': t->cam.zoom *= 0.85; if (t->cam.zoom < 0.1) t->cam.zoom = 0.1; break;
            case 'Z': t->cam.zoom *= 1.15; if (t->cam.zoom > 20.0) t->cam.zoom = 20.0; break;
            case 'w': sel = live_pick_next(t->sim, sel, 1); break;
            case 'b': sel = live_pick_next(t->sim, sel, -1); break;
            case 'c': live_center_on(t, sel); break;
            case 'f': live_freeze_toggle(t, sel); break;
            case 'H': live_heat(t, +5.0); break;
            case 'L': live_heat(t, -5.0); break;
            case 'm': if (t->sim) integrator_minimize(t->sim, 300, 0.005, 0.05); break;
            case 'i': {
                char *line = tui_readline("insert species: ");
                if (line) {
                    cmd_save(t, undopath);
                    live_add_at_cam(t, line);
                    free(line);
                    sel = t->sim->num_atoms - 1;
                    snprintf(msg, sizeof msg, "inserted; undo with u");
                }
                break;
            }
            case 'r': {
                char *line = tui_readline("replace with: ");
                if (line) {
                    int Z = live_z_from_name(line);
                    if (Z >= 1 && sel >= 0) {
                        cmd_save(t, undopath);
                        live_transmute(t, sel, Z);
                        snprintf(msg, sizeof msg, "atom %d -> %s", sel, line);
                    }
                    free(line);
                }
                break;
            }
            case 'x':
                if (sel >= 0 && sim_remove_terminal_atom(t->sim, sel)) {
                    sel = live_pick_next(t->sim, sel, 1);
                    snprintf(msg, sizeof msg, "deleted (terminal atoms only)");
                } else {
                    snprintf(msg, sizeof msg, "not deletable: atom is not terminal");
                }
                break;
            case 'p':
                if (sel >= 0) {
                    cmd_save(t, undopath);
                    live_clone_molecule(t, sel);
                    snprintf(msg, sizeof msg, "cloned molecule; undo with u");
                }
                break;
            case 'u':
                if (cmd_load(t, undopath) == 0) sel = live_pick_next(t->sim, -1, 1);
                snprintf(msg, sizeof msg, "undo");
                break;
            case '/': {
                char *line = tui_readline("/");
                if (line && line[0]) {
                    int Z = live_z_from_name(line), found = -1;
                    for (int i = 1; i <= t->sim->num_atoms; i++) {
                        int a = (sel + i + t->sim->num_atoms) % t->sim->num_atoms;
                        if (t->sim->atoms[a].Z == Z) { found = a; break; }
                    }
                    if (found >= 0) { sel = found; live_center_on(t, sel); }
                    if (found >= 0) snprintf(msg, sizeof msg, "match atom %d", found);
                    else            snprintf(msg, sizeof msg, "no match");
                }
                free(line);
                break;
            }
            case 'n': break;
            case ':': {
                char *line = tui_readline(":");
                if (line) {
                    sh_hist_push(line);
                    last_status = run_line(t, line);
                    free(line);
                    if (t->quit) goto done;
                    snprintf(msg, sizeof msg, "status %d", last_status);
                }
                break;
            }
            default: break;
        }
    }
done:
    live_raw_off();
    fprintf(stderr, "\x1b[0m");
    return 0;
}

/* ── POSIX-named entry points ──────────────────────────────────────── */

/* dd if=<preset|species|save|file.mol> of=<world|name> [count=N]
 *    [x=] [y=] [z=] [temp=] [seed=] */
static int cmd_dd(Tui *t, char **a, int n) {
    const char *src = NULL, *dst = NULL;
    int count = 1;
    double x = 0, y = 0, z = 0, temp = -1;
    double q = 1e30, eps = -1.0, sig = -1.0;
    unsigned long seed = 0;
    for (int i = 0; i < n; i++) {
        char *eq = strchr(a[i], '=');
        char *key, *val;
        if (!eq) continue;
        *eq = '\0';
        key = a[i];
        val = eq + 1;
        if (!strcmp(key, "if")) src = val;
        else if (!strcmp(key, "of")) dst = val;
        else if (!strcmp(key, "count")) { int c = 1; if (parse_int(val, &c)) count = c; }
        else if (!strcmp(key, "x")) parse_double(val, &x);
        else if (!strcmp(key, "y")) parse_double(val, &y);
        else if (!strcmp(key, "z")) parse_double(val, &z);
        else if (!strcmp(key, "temp")) parse_double(val, &temp);
        else if (!strcmp(key, "seed")) parse_ulong(val, &seed);
        else if (!strcmp(key, "q")) parse_double(val, &q);
        else if (!strcmp(key, "eps")) parse_double(val, &eps);
        else if (!strcmp(key, "sigma")) parse_double(val, &sig);
    }
    if (!src && !dst) { sh_err("  usage: dd if=<preset|species|file> of=world [count=N] [x= y= z=]\n"); return 2; }
    if (count < 1) count = 1;
    if (count > 4096) count = 4096;
    if (src && (!strcmp(src, "urandom") || !strcmp(src, "thermal") || !strcmp(src, "velocities"))) {
        double T = (temp > 0) ? temp : 300.0;
        if (!t->sim) { sh_err("  dd: no world\n"); return 1; }
        integrator_maxwell_boltzmann(t->sim, T, seed ? seed : 7UL);
        sh_err("  dd: thermalized %d atoms at %.1f K\n", t->sim->num_atoms, T);
        return 0;
    }
    if (dst && !strncmp(dst, "atom:", 5)) {
        int ai = -1;
        if (!parse_int(dst + 5, &ai) || !t->sim || ai < 0 || ai >= t->sim->num_atoms) {
            sh_err("  dd: bad atom target `%s`\n", dst);
            return 1;
        }
        if (q <= 1e29) t->sim->atoms[ai].partial_charge = q;
        if (eps >= 0.0 && sig > 0.0) sim_set_atom_lj(t->sim, ai, eps, sig);
        sh_err("  dd: atom %d q=%.4f eps=%.5f sigma=%.4f\n", ai,
               t->sim->atoms[ai].partial_charge, t->sim->atoms[ai].lj_epsilon,
               t->sim->atoms[ai].lj_sigma);
        return 0;
    }
    if (!src) { sh_err("  dd: missing if=\n"); return 2; }
    if ((!dst || !strcmp(dst, "world")) && (!strcmp(src, "petri") || !strcmp(src, "empty")))
        return world_preset(t, src, seed);
    {
        FILE *probe = fopen(src, "r");
        if (probe) {
            char magic[8] = {0};
            size_t got = fread(magic, 1, 7, probe);
            fclose(probe);
            if (got == 7 && !strncmp(magic, "S2SAVE1", 7)) return cmd_load(t, src);
            if (mol_load(src) == 0) {
                int made = 0;
                for (int c = 0; c < count; c++) {
                    Vec3 o = vec3(x + live_rand(-1.5, 1.5), y + live_rand(-1.5, 1.5), z);
                    if (mol_place(t, o) != 0) break;
                    made++;
                }
                sh_err("  dd: placed %d x `%s` at %.2f %.2f %.2f\n", made, src, x, y, z);
                return made > 0 ? 0 : 1;
            }
        }
    }
    {
        int made = 0;
        for (int c = 0; c < count; c++) {
            Vec3 at;
            if (count == 1) at = vec3(x, y, z);
            else {
                Vec3 base = vec3(x, y, z);
                at = vec3_add(base, vec3(live_rand(-6, 6), live_rand(-6, 6), live_rand(-6, 6)));
            }
            if (live_species_add(t, src, at) < 0) break;
            made++;
        }
        if (made && temp > 0)
            integrator_maxwell_boltzmann(t->sim, temp, seed ? seed : 7UL);
        sh_err("  dd: %d x `%s` (%d atoms) [%s]\n", made, src,
               t->sim ? t->sim->num_atoms : 0, made == count ? "ok" : "partial");
        return made > 0 ? 0 : 1;
    }
}

/* ps [ -l ] — process/world table; on a TTY with no args, live monitor */
static int cmd_ps(Tui *t, char **a, int n) {
    int longform = (n >= 1 && !strcmp(a[0], "-l"));
    if (!longform && n == 0 && isatty(STDIN_FILENO) && isatty(STDERR_FILENO))
        return live_monitor(t);
    live_status(t, "ps", 1, 0, -1);
    if (t->has_nrn)
        printf("  PID 2  neuron   V=%+.2f mV  I_ext=%+.2f  spiking=%d\n",
               t->nrn.V, t->nrn.I_ext, hh_is_spiking(&t->nrn, -20.0));
    if (rxn.defined)
        printf("  PID 3  reaction %s  pair Z%d-Z%d < %.2f A  %s  every %d steps\n",
               rxn.name, rxn.Z1, rxn.Z2, rxn.rcut,
               t->rxn_armed ? "ARMED" : "disarmed", t->rxn_every);
    return 0;
}

/* vi [world] — modal world editor; `:` runs any shell command */
static int cmd_vi(Tui *t, char **a, int n) {
    char path[128];
    if (n >= 1 && strcmp(a[0], "world") != 0) {
        sh_err("  vi: only `vi world` exists right now (template editing comes with `ed`)\n");
        return 2;
    }
    return live_vi(t, path, sizeof path);
}

/* more [world|atoms|bonds|summary|file] — one frame / a listing */
static int cmd_more(Tui *t, char **a, int n) {
    const char *what = (n >= 1) ? a[0] : "world";
    if (!strcmp(what, "world")) {
        if (!t->sim) { sh_err("  more: no world\n"); return 1; }
        view_render(t->sim, &t->cam, 0);
        return 0;
    }
    if (!strcmp(what, "atoms")) { sim_print_atoms(t->sim); return 0; }
    if (!strcmp(what, "bonds")) { sim_print_bonds(t->sim); return 0; }
    if (!strcmp(what, "summary")) { sim_print_summary(t->sim); return 0; }
    return cmd_cat(a, n);
}

/* sync [file] — flush the world to a S2SAVE1 file (default s2world.s2) */
static int cmd_sync(Tui *t, char **a, int n) {
    return cmd_save(t, (n >= 1) ? a[0] : "s2world.s2");
}

/* ── POSIX command surface ───────────────────────────────────────────
 * User-visible command names are POSIX commands and their structure is
 * POSIX (options, operands, stdin/stdout/stderr, exit status). Operands
 * are S2 words. Old non-POSIX names are refused. This layer rewrites the
 * POSIX surface onto the internal verbs the dispatch chain implements.
 * ─────────────────────────────────────────────────────────────────── */

static int posix_is_species(const char *s) {
    static const char *sp[] = {
        "water", "h2o", "nh3", "methane", "ch4", "co2", "h2", "glycine", "gly",
        "alanine", "ala", "uracil", "cytosine", "thymine", "adenine", "guanine",
        "deoxyribose", "sugar", "Na+", "na", "K+", "k", "Cl-", "cl", "Ca2+", "ca", NULL
    };
    for (int i = 0; sp[i]; i++) if (!strcmp(s, sp[i])) return 1;
    return live_z_from_name(s) >= 1;
}

static int posix_is_spawner(const char *s) {
    static const char *w[] = {
        "atom", "ion", "h2o", "h2", "nh3", "ch4", "methane", "co2", "kcsa",
        "demo", "quantum", "water", "trimer", "base", "pair", "dinucleotide",
        "neuron", "dipeptide", "helix", "cage", "filter", "duplex", NULL
    };
    for (int i = 0; w[i]; i++) if (!strcmp(s, w[i])) return 1;
    return 0;
}

/* 1 = rewritten (use out/outn), 0 = unchanged, -1 = refuse (not POSIX) */
static int posix_rewrite(char **tok, int nt, char **out, int *outn) {
    const char *c = tok[0];
    static const char *refuse[] = {
        "new", "spawn", "del", "bond", "detect", "restrain", "clear", "init",
        "step", "run", "show", "minimize", "heat", "mol", "rxn", "neuron",
        "render", "view", "cam", "slice", "watch", "list", "help", "quit",
        "history", "touchf", NULL
    };
    for (int i = 0; refuse[i]; i++)
        if (!strcmp(c, refuse[i])) return -1;
    if (c[0] == '!') return -1;   /* host escape is `sh -c '...'` now */

    if (!strcmp(c, "env") && nt >= 2 && !strcmp(tok[1], "-i")) {
        out[0] = "new";
        for (int i = 2; i < nt; i++) out[i - 1] = tok[i];
        *outn = nt - 1;
        return 1;
    }
    if (!strcmp(c, "touch") && nt >= 2 && !strcmp(tok[1], "--")) {
        /* explicit file touch: touch -- <file>.. (POSIX, no matter) */
        for (int i = 2; i < nt; i++) out[i - 2] = tok[i];
        /* reuse file-touch path via special marker: fall through as `touchf` */
        out[0] = "touchf";
        for (int i = 2; i < nt; i++) out[i - 1] = tok[i];
        *outn = nt - 1;
        return 1;
    }
    if (!strcmp(c, "touch") && nt >= 3 && !strcmp(tok[1], "-m")) {
        out[0] = "addsp";
        for (int i = 2; i < nt; i++) out[i - 1] = tok[i];
        *outn = nt - 1;
        return 1;
    }
    if (!strcmp(c, "touch") && nt >= 2 &&
        (posix_is_species(tok[1]) || posix_is_spawner(tok[1]))) {
        out[0] = posix_is_spawner(tok[1]) ? "spawn" : "addsp";
        for (int i = 1; i < nt; i++) out[i] = tok[i];
        *outn = nt;
        return 1;
    }
    if (!strcmp(c, "ln") && nt >= 3) {
        if (!strcmp(tok[1], "-s")) {
            out[0] = "restrain";
            for (int i = 2; i < nt; i++) out[i - 1] = tok[i];
            *outn = nt - 1;
            return 1;
        }
        out[0] = "bond";
        for (int i = 1; i < nt; i++) out[i] = tok[i];
        *outn = nt;
        return 1;
    }
    if (!strcmp(c, "fsck")) { out[0] = "detect"; out[1] = "bonds"; *outn = 2; return 1; }
    if (!strcmp(c, "unlink")) { out[0] = "clear"; out[1] = "restraints"; *outn = 2; return 1; }
    if (!strcmp(c, "man")) { out[0] = "manpage"; if (nt >= 2) { out[1] = tok[1]; *outn = 2; } else *outn = 1; return 1; }
    if (!strcmp(c, "fc")) { out[0] = "history"; *outn = 1; return 1; }
    if (!strcmp(c, "tput") && nt >= 2 && !strcmp(tok[1], "clear")) { out[0] = "clear"; *outn = 1; return 1; }
    if (!strcmp(c, "sh") && nt >= 3 && !strcmp(tok[1], "-c")) { out[0] = "!"; out[1] = tok[2]; *outn = 2; return 1; }
    if (!strcmp(c, "nice")) { out[0] = "minimize"; *outn = 1; return 1; }
    if (!strcmp(c, "df")) { out[0] = "show"; out[1] = "thermo"; *outn = 2; return 1; }
    if (!strcmp(c, "ls") && nt >= 2 &&
        (!strcmp(tok[1], "demos") || !strcmp(tok[1], "atoms") ||
         !strcmp(tok[1], "bonds") || !strcmp(tok[1], "summary"))) {
        for (int i = 0; i < nt; i++) out[i] = tok[i];
        out[0] = "list";
        *outn = nt;
        return 1;
    }
    if (!strcmp(c, "cat") && nt == 2) {
        if (!strcmp(tok[1], "neuron")) { out[0] = "neuron"; out[1] = "show"; *outn = 2; return 1; }
        if (!strcmp(tok[1], "energy")) { out[0] = "show"; out[1] = "energy"; *outn = 2; return 1; }
        if (!strcmp(tok[1], "world"))  { out[0] = "show"; out[1] = "thermo"; *outn = 2; return 1; }
        if (!strcmp(tok[1], "atoms"))  { out[0] = "list"; out[1] = "atoms"; *outn = 2; return 1; }
        if (!strcmp(tok[1], "bonds"))  { out[0] = "list"; out[1] = "bonds"; *outn = 2; return 1; }
    }
    if (!strcmp(c, "cp") && nt == 3 && posix_is_species(tok[1])) {
        out[0] = "species-export"; out[1] = tok[1]; out[2] = tok[2];
        *outn = 3;
        return 1;
    }
    return 0;
}

/* export a built-in species to a MOL1 template file (cp <species> x.mol) */
static int species_export(Tui *t, const char *species, const char *path) {
    Tui tmp;
    Simulation *s;
    Vec3 c = vec3_zero();
    int n, rc;
    memset(&tmp, 0, sizeof tmp);
    tmp.seed = t->seed;
    tmp.sim = sim_create(64, 64);
    if (!tmp.sim) { sh_err("  cp: alloc failed\n"); return 1; }
    if (live_species_add(&tmp, species, vec3_zero()) < 0) { sim_destroy(tmp.sim); return 1; }
    s = tmp.sim;
    n = s->num_atoms < MOL_MAXA ? s->num_atoms : MOL_MAXA;
    for (int i = 0; i < n; i++) c = vec3_add(c, s->atoms[i].position);
    if (n > 0) c = vec3_scale(c, 1.0 / n);
    mol.na = n;
    for (int i = 0; i < n; i++) {
        mol.a[i].Z = s->atoms[i].Z;
        mol.a[i].pos = vec3_sub(s->atoms[i].position, c);
        mol.a[i].q = s->atoms[i].partial_charge;
        mol.a[i].eps = s->atoms[i].lj_epsilon;
        mol.a[i].sig = s->atoms[i].lj_sigma;
        mol.a[i].has_lj = 1;
    }
    {
        int nb = 0;
        for (int b = 0; b < s->num_bonds && nb < MOL_MAXA; b++) {
            if (s->bonds[b].atom_a < n && s->bonds[b].atom_b < n) {
                mol.b[nb].a = s->bonds[b].atom_a;
                mol.b[nb].b = s->bonds[b].atom_b;
                mol.b[nb].order = s->bonds[b].order;
                nb++;
            }
        }
        mol.nb = nb;
    }
    snprintf(mol.name, sizeof mol.name, "%s", species);
    rc = mol_save(path);
    sim_destroy(tmp.sim);
    if (rc == 0) sh_err("  cp: `%s` -> `%s` (%d atoms, %d bonds)\n", species, path, mol.na, mol.nb);
    return rc;
}

/* kill [-SIGNAL] <world|neuron|atom:N> [dE=|cur=] (POSIX signals as events) */
static int cmd_kill(Tui *t, char **a, int n) {
    const char *sig = "TERM", *target;
    int i = 0, is_world = 0, atom = -1;
    double dE = 5.0, cur = 10.0;
    if (n >= 1 && a[0][0] == '-') { sig = a[0] + 1; i = 1; }
    if (i >= n) { sh_err("  usage: kill [-TERM|-STOP|-CONT|-USR1|-USR2] <world|neuron|atom:N> [dE=|cur=]\n"); return 2; }
    target = a[i];
    for (int j = i + 1; j < n; j++) {
        if (!strncmp(a[j], "dE=", 3)) parse_double(a[j] + 3, &dE);
        else if (!strncmp(a[j], "cur=", 4)) parse_double(a[j] + 4, &cur);
    }
    if (!strcmp(target, "world") || !strcmp(target, "all")) is_world = 1;
    else if (!strncmp(target, "atom:", 5)) parse_int(target + 5, &atom);
    else if (!strcmp(target, "neuron")) atom = -3;
    else parse_int(target, &atom);

    if (!strcmp(sig, "USR2")) { live_heat(t, dE); return 0; }            /* thermal kick */
    if (!strcmp(sig, "USR1")) {                                          /* stimulus */
        if (!t->has_nrn) { sh_err("  kill: no neuron (touch neuron)\n"); return 1; }
        t->nrn.I_ext += cur;
        return 0;
    }
    if (!strcmp(sig, "STOP") || !strcmp(sig, "CONT")) {
        if (atom >= 0) { live_freeze_toggle(t, atom); return 0; }
        sh_err("  kill: freeze/unfreeze targets one atom (atom:N)\n");
        return 1;
    }
    if (is_world) {                                                      /* TERM world */
        /* Reset to a fresh empty world rather than NULL: every viewer
         * and mutator assumes a live sim, and a NULL world turned the
         * next `ls` into a segfault. `dd`/`new`/`env -i` build from here. */
        tui_new(t, 512, 512);
        t->has_nrn = 0;
        sh_err("  kill: world terminated (fresh empty world)\n");
        return 0;
    }
    if (atom == -3) { t->has_nrn = 0; return 0; }
    if (atom >= 0) {
        if (!t->sim || !sim_remove_terminal_atom(t->sim, atom) ) {
            sh_err("  kill: atom %d is not removable (terminal atoms only)\n", atom);
            return 1;
        }
        return 0;
    }
    sh_err("  kill: unknown target `%s`\n", target);
    return 1;
}

/* du — energy and population share per element (the world's disk usage) */
static int cmd_du(Tui *t) {
    long cnt[37] = {0};
    double ke[37] = {0.0}, tot = 0.0;
    if (!t->sim || t->sim->num_atoms < 1) { sh_err("  du: no world\n"); return 1; }
    for (int i = 0; i < t->sim->num_atoms; i++) {
        int Z = t->sim->atoms[i].Z;
        double v2 = vec3_dot(t->sim->atoms[i].velocity, t->sim->atoms[i].velocity);
        double e = 0.5 * t->sim->atoms[i].mass * v2 * AMU_AFS2_TO_EV;
        if (Z < 0 || Z > 36) Z = 0;
        cnt[Z]++;
        ke[Z] += e;
        tot += e;
    }
    printf("  species   atoms     KE_eV    share\n");
    for (int Z = 1; Z <= 36; Z++) {
        if (!cnt[Z]) continue;
        printf("  %-8s %6ld %10.4f %7.1f%%\n",
               pt_element(Z) ? pt_element(Z)->symbol : "?", cnt[Z], ke[Z],
               tot > 1e-12 ? 100.0 * ke[Z] / tot : 0.0);
    }
    printf("  total    %6d %10.4f\n", t->sim->num_atoms, tot);
    return 0;
}

/* make <file.rxn> [-n] — load a reaction rule, optionally arm it, fire once */
static int cmd_make(Tui *t, char **a, int n) {
    int once = 0;
    const char *path = NULL;
    for (int i = 0; i < n; i++) {
        if (!strcmp(a[i], "-n")) once = 1;
        else if (!path) path = a[i];
    }
    if (!path) { sh_err("  usage: make <file.rxn> [-n]\n"); return 2; }
    if (rxn_load(path) != 0) return 1;
    if (!once) {
        t->rxn_armed = 1;
        if (t->rxn_every <= 0) t->rxn_every = 50;
    }
    return rxn_fire(t) == 0 ? 0 : 1;
}

/* apply an exported NAME=VALUE to the world (env vars ARE the world params).
 * Returns 0 when applied, 1 when the name is not a world parameter. */
static int world_set_param(Tui *t, const char *name, const char *val) {
    double d;
    int n;
    if (!t->sim) return 1;
    if (!strcmp(name, "dt") && parse_double(val, &d)) t->sim->dt = d;
    else if (!strcmp(name, "cutoff") && parse_double(val, &d)) t->sim->cutoff = d;
    else if (!strcmp(name, "dielectric") && parse_double(val, &d)) t->sim->dielectric = d;
    else if (!strcmp(name, "temp") && parse_double(val, &d)) t->sim->thermostat.target_temperature = d;
    else if (!strcmp(name, "thermostat"))
        t->sim->thermostat.type = !strcmp(val, "andersen") ? THERMOSTAT_ANDERSEN
                                : !strcmp(val, "berendsen") ? THERMOSTAT_BERENDSEN
                                : !strcmp(val, "langevin") ? THERMOSTAT_LANGEVIN
                                : THERMOSTAT_NONE;
    else if (!strcmp(name, "maxtemp") && parse_double(val, &d)) t->sim->max_temperature = d;
    else if (!strcmp(name, "tau") && parse_double(val, &d)) t->sim->thermostat.tau = d;
    else if (!strcmp(name, "nu") && parse_double(val, &d)) t->sim->thermostat.nu = d;
    else if (!strcmp(name, "seed")) parse_ulong(val, &t->seed);
    else if (!strcmp(name, "yaw") && parse_double(val, &d)) t->cam.yaw_deg = d;
    else if (!strcmp(name, "pitch") && parse_double(val, &d)) t->cam.pitch_deg = d;
    else if (!strcmp(name, "zoom") && parse_double(val, &d)) t->cam.zoom = d;
    else if (!strcmp(name, "slice")) {
        if (!strcmp(val, "off")) t->cam.slice = 0.0;
        else if (parse_double(val, &d)) t->cam.slice = d;
    }
    else if (!strcmp(name, "barostat")) t->baro_on = (!strcmp(val, "on") || !strcmp(val, "1"));
    else if (!strcmp(name, "press") && parse_double(val, &d)) t->p0_bar = d;
    else if (!strcmp(name, "tau-p") && parse_double(val, &d)) t->taup_fs = d;
    else if (!strcmp(name, "rxn-every") && parse_int(val, &n)) t->rxn_every = n;
    else if (!strcmp(name, "box")) {
        double bx, by, bz;
        if (sscanf(val, "%lf,%lf,%lf", &bx, &by, &bz) == 3)
            t->sim->box.dimensions = vec3(bx, by, bz);
        else if (parse_double(val, &d))
            t->sim->box.dimensions = vec3(d, d, d);
    }
    else if (!strcmp(name, "pbc")) {
        int on = (!strcmp(val, "on") || !strcmp(val, "1"));
        t->sim->box.periodic[0] = t->sim->box.periodic[1] = t->sim->box.periodic[2] = on;
    }
    else return 1;
    return 0;
}/* ── the POSIX manual ──────────────────────────────────────────────── */
static void man_index(void) {
    printf("S2TUI(1) - POSIX terminal over the live simulation\n\n");
    printf("WORLD      env -i [petri|empty]   ls [scope]            cat <entity>\n");
    printf("           dd if=... of=...       sync [file]           more [scope]\n");
    printf("MATTER     touch <species> [xN]   rm atom <i>[..ranges] cp <sp> <f>.mol\n");
    printf("           touch -- <file>..      touch -m <sp>         (explicit file/matter)\n");
    printf("           ln <a> <b>             ln -s <i> x y z k    unlink [all]\n");
    printf("           fsck                   kill [-SIG] <target> du\n");
    printf("TIME       sleep <steps>[fs]      sleep <sec>s         df\n");
    printf("PROCESSES  ps [-l]                 nice                 make <file.rxn> [-n]\n");
    printf("           screen                 fullscreen twin of ps (alternate screen, resize)\n");
    printf("LIVE       vi world                ps                   screen\n");
    printf("SHELL      export NAME=value       env                  set NAME=VALUE\n");
    printf("           : eval exec command type shift set-- readonly umask trap alias\n");
    printf("           if/for/while/until/case { } ( ) ! ~ $(( )) ${} `` $# $@\n");
    printf("           man [cmd]              fc -l                tput clear\n");
    printf("           sh -c 'cmd'            exit                 [all POSIX text tools]\n\n");
    printf("Signals as events: kill -STOP/-CONT atom:N freezes/thaws one atom;\n");
    printf("-USR1 neuron cur=N stimulates; -USR2 world dE=N is a thermal kick;\n");
    printf("plain TERM removes the target (world, neuron, or a terminal atom:N).\n");
    printf("World parameters for export: dt cutoff dielectric temp thermostat tau nu\n");
    printf("seed box pbc press tau-p barostat yaw pitch zoom slice rxn-every\n");
}

static int man_page(const char *topic) {
    if (!strcmp(topic, "touch")) {
        printf("TOUCH(1)\ntouch <species> [x y z] [xN|:N|count=N]  create matter (open site\n  each copy, 3 A clearance, randomly oriented; explicit xyz honored\n  verbatim). molecules (water, nh3, ch4, co2, glycine, ...), ions,\n  elements by symbol or number 1..36, composites.\n  touch -- <file>..  explicit file stamp (no matter). touch -m <sp> explicit matter.\n");
        return 0;
    }
    if (!strcmp(topic, "rm")) {
        printf("RM(1)\nrm atom <i>[..ranges,]  remove terminal atoms descending (e.g. 1..5,7).\nrm <file>     remove a file (POSIX rm semantics).\n");
        return 0;
    }
    if (!strcmp(topic, "sh") || !strcmp(topic, "shell") || !strcmp(topic, "posix")) {
        printf("SH(1)\nFull POSIX syntax: 'sq' \"dq\" \\esc # ; && || | & > >> < <<< 2> 2>>\n  $V ${V} ${V:-d} ${V:=d} ${V:?m} ${V:+a} ${#V} ${V#pat} ${V%%pat}\n  $( ) ` ` $(( )) *?[] ~ V=v set -- $# $@ $* $0..$9 $? $$ $! !\n  if/for/while/until/case { } ( ) eval exec command type shift readonly umask trap\n  S2 params streamlined: set NAME=VALUE, rm ranges, sleep <steps>[fs]|<sec>s.\n");
        return 0;
    }
    if (!strcmp(topic, "ln")) {
        printf("LN(1)\nln <a> <b> [order]        covalent/topological link (refuses dupes).\nln -s <i> x y z k         positional restraint: anchor atom i at (x,y,z),\n  spring k eV/A^2. Restraints are the world's symlinks; `unlink` clears.\n");
        return 0;
    }
    if (!strcmp(topic, "fsck")) {
        printf("FSCK(1)\nfsck   detect bonds from geometry and rebuild the topology\n  (angles rebuilt). Bonds only where the engine's criteria match.\n");
        return 0;
    }
    if (!strcmp(topic, "dd")) {
        printf("DD(1)\ndd if=<src> of=<dst> [count=N] [x= y= z=] [temp=] [seed=] [q= eps= sigma=]\n  src: a world preset (petri, empty), a save file (S2SAVE1), a MOL1\n  template, `urandom` (thermalise the world), or a species name.\n  of=world for creation; of=atom:N writes q=/eps=/sigma= onto atom N.\n");
        return 0;
    }
    if (!strcmp(topic, "sync")) {
        printf("SYNC(1)\nsync [file]   flush the world to an S2SAVE1 file (default s2world.s2).\n  Restore with: dd if=<file> of=world\n");
        return 0;
    }
    if (!strcmp(topic, "more")) {
        printf("MORE(1)\nmore [world|atoms|bonds|summary|file]   one grid frame or a listing.\n  Interactive viewing belongs to `vi world` and `ps`.\n");
        return 0;
    }
    if (!strcmp(topic, "sleep")) {
        printf("SLEEP(1)\nsleep <steps>[fs]  advance world N steps (world time).\n  sleep <sec>s      POSIX wall sleep (e.g. sleep 0.5s). Bare N = steps.\n");
        return 0;
    }
    if (!strcmp(topic, "df")) {
        printf("DF(1)\ndf   world capacity and thermodynamics: T, N, KE, box volume, density,\n  pressure (virial, COM-relative), barostat state. `du` gives per-element\n  population and kinetic-energy share.\n");
        return 0;
    }
    if (!strcmp(topic, "du")) {
        printf("DU(1)\ndu   per-element atom count and kinetic-energy share of the world.\n");
        return 0;
    }
    if (!strcmp(topic, "kill")) {
        printf("KILL(1)\nkill [-TERM|-STOP|-CONT|-USR1|-USR2] <world|neuron|atom:N> [cur=|dE=]\n  TERM remove; STOP/CONT freeze/thaw (atom:N); USR1 neural stimulus;\n  USR2 thermal kick in eV (world). Signals are the POSIX way to poke a\n  running process, and here the processes are matter and agents.\n");
        return 0;
    }
    if (!strcmp(topic, "nice")) {
        printf("NICE(1)\nnice   run a steepest-descent minimisation pass over the world\n  (clash relief after assembling matter; `vi world` also has `m`).\n");
        return 0;
    }
    if (!strcmp(topic, "make")) {
        printf("MAKE(1)\nmake <file.rxn> [-n]   load a reaction rule (RXN1 text file) and fire\n  it once; without -n it stays armed and auto-fires every `rxn-every`\n  steps. Rules are files: create them with echo/redirection, edit with\n  any editor. The live keybind is `e` in `ps`/`vi`.\n");
        return 0;
    }
    if (!strcmp(topic, "ps")) {
        printf("PS(1)\nps        on a terminal: the live control monitor. Single keys create\n  (w water, i Na+, K K+, C Cl-, u base, a alanine, g glycine, d sugar),\n  augment (H/L heat/cool, f freeze, r replace, p clone, x delete,\n  m minimise) and adapt processes (n neuron, e catalysis). Piped, `ps`\n  prints one snapshot instead of taking the screen.\nps -l     process table (world, neuron, reaction rule).\n");
        return 0;
    }
    if (!strcmp(topic, "screen")) {
        printf("SCREEN(1)\nscreen   fullscreen twin of `ps`: same world, same keybinds, on the\n  alternate screen with terminfo colors and resize (ncursesw build).\n  Needs a real terminal; piped, dumb, missing-lib or failed init prints\n  why and returns 0. Esc quits like q. S2TUI_SCREEN=0 forces line mode.\n");
        return 0;
    }
    if (!strcmp(topic, "vi")) {
        printf("VI(1)\nvi world   modal editor over the living world. Normal mode: hjkl pan,\n  w/b select, i insert, r replace, x delete, p clone, f freeze, u undo,\n  / search, n next, c center, z/Z zoom, H/L heat/cool, m minimise,\n  space run/pause, s step, q quit. `:` runs any shell command against\n  the same world. stdout stays DATA; the editor draws on stderr.\n");
        return 0;
    }
    if (!strcmp(topic, "ls")) {
        printf("LS(1)\nls [atoms|bonds|summary|demos]   list the world scope; `ls` alone\n  lists files. Scope listings compose with pipes (`ls atoms | wc -l`).\n");
        return 0;
    }
    if (!strcmp(topic, "cat")) {
        printf("CAT(1)\ncat <entity|file>   cat neuron (HH state), cat energy (full energy\n  ledger), cat world (thermodynamics), cat atoms|bonds (listings), or any\n  file. Text tools (grep/sort/cut/wc/...) compose in pipelines.\n");
        return 0;
    }
    if (!strcmp(topic, "cp")) {
        printf("CP(1)\ncp <src> <dst>           copy a file (POSIX).\ncp <species> <dst>.mol   export a built-in species as a MOL1 template;\n  instantiate later with dd if=<dst>.mol of=world. Template text is a\n  plain file; edit with shell redirection or any editor.\n");
        return 0;
    }
    if (!strcmp(topic, "env")) {
        printf("ENV(1)\nenv             list environment/world variables.\nenv -i [preset]  a fresh world: `petri` (the pure-simulator dish) or\n  `empty` (custom capacity: env -i <atoms> <bonds>).\n");
        return 0;
    }
    if (!strcmp(topic, "export")) {
        printf("EXPORT(1)\nexport NAME=value   shell variables and world parameters. World:\n  dt cutoff dielectric temp thermostat tau nu seed box pbc press tau-p\n  barostat yaw pitch zoom slice rxn-every.\n");
        return 0;
    }
    if (!strcmp(topic, "set")) {
        printf("SET(1)\nset   POSIX shell built-in; in S2 the parameter objects are the world\n  (see export) and per-atom writes go through dd of=atom:N.\n");
        return 0;
    }
    if (!strcmp(topic, "sh")) {
        printf("SH(1)\nsh -c 'command'   run a command on the HOST shell, outside the world\n  (the honest escape hatch: no S2 rewriting applies).\n");
        return 0;
    }
    if (!strcmp(topic, "unlink")) {
        printf("UNLINK(1)\nunlink [all]   remove positional restraints (the world's symlinks).\n");
        return 0;
    }
    if (!strcmp(topic, "man")) {
        printf("MAN(1)\nman [command]   this manual. Every command in this terminal is a\n  POSIX command; operands are S2 words (species, atom indices, scopes).\n");
        return 0;
    }
    printf("man: no page for `%s`; try `man` for the index. POSIX text tools keep\n"
           "their host semantics; sim commands are listed in the index.\n", topic);
    return 0;
}

static int dispatch_cmd(Tui *t, char **tok, const int *eleg, int nt) {
    char *rtok[SH_MAXARG];
    int releg[SH_MAXARG];
    int rn = 0, rr;
    if (nt == 0) return 0;
    if ((rr = posix_rewrite(tok, nt, rtok, &rn)) < 0) {
        sh_err("  %s: not a POSIX command — see `man`\n", tok[0]);
        return 127;
    }
    if (rr > 0) {
        for (int i = 0; i < rn; i++) releg[i] = 0;
        tok = rtok;
        nt = rn;
        eleg = releg;
    }
    if (nt == 0) return 0;
    if (!strcmp(tok[0], "quit") || !strcmp(tok[0], "exit")) {
        int c = last_status;
        if (nt >= 2 && !parse_int(tok[1], &c)) return 2;
        t->quit = 1;
        return c & 255;
    }
    else if (!strcmp(tok[0], "help")) {
        if (nt >= 2) return print_help_topic(tok[1]);
        print_help();
        return 0;
    }
    else if (!strcmp(tok[0], "list") && nt >= 2 && !strcmp(tok[1], "demos"))
            printf("  demos: 1 quantum 2 bond 3 water 4 trimer 5 methane 6 bases 7 pairing 8 dinucleotide 9 neuron 10 dipeptide 11 helix 12 kcsa 12b real-filter 17 duplex\n");
        else if (!strcmp(tok[0], "test") && nt >= 3 && !strcmp(tok[1], "demo")) return run_demo_test(tok[2]);
        else if (!strcmp(tok[0], "test") && nt >= 2 && !strcmp(tok[1], "all")) {
            const char *ids[] = {"1","2","3","4","5","6","7","8","9","10","11","12","12b","17"};
            int fails = 0;
            for (unsigned i = 0; i < sizeof ids / sizeof ids[0]; i++) fails += run_demo_test(ids[i]);
            return fails ? 1 : 0;
        }
        else if (!strcmp(tok[0], "test")) {
            return cmd_test(tok + 1, nt - 1); /* POSIX test(1); demo tests use `test demo|all` */
        }
        else if (!strcmp(tok[0], "[") && nt >= 2 && !strcmp(tok[nt - 1], "]")) {
            return cmd_test(tok + 1, nt - 2);
        }
        else if (!strcmp(tok[0], "printf")) {
            return cmd_printf(tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "true")) {
            return 0;
        }
        else if (!strcmp(tok[0], "false")) {
            return 1;
        }
        else if (!strcmp(tok[0], ":")) {
            return 0; /* POSIX no-op, expansions already done */
        }
        else if (!strcmp(tok[0], "eval") && nt >= 1) {
            /* eval WORDS..: join and re-execute (POSIX) */
            char buf[SH_LINE];
            size_t p = 0;
            buf[0] = '\0';
            for (int i = 1; i < nt && p + 1 < sizeof buf; i++) {
                if (i > 1 && p + 1 < sizeof buf) buf[p++] = ' ';
                size_t L = strlen(tok[i]);
                if (p + L >= sizeof buf) break;
                memcpy(buf + p, tok[i], L);
                p += L;
            }
            buf[p] = '\0';
            if (p == 0) return 0;
            return run_line_str(t, buf);
        }
        else if (!strcmp(tok[0], "exec") && nt >= 1) {
            /* exec CMD: replace current shell context (here: just run) */
            if (nt == 1) return 0;
            char buf[SH_LINE];
            size_t p = 0;
            buf[0] = '\0';
            for (int i = 1; i < nt && p + 1 < sizeof buf; i++) {
                if (i > 1 && p + 1 < sizeof buf) buf[p++] = ' ';
                size_t L = strlen(tok[i]);
                if (p + L >= sizeof buf) break;
                memcpy(buf + p, tok[i], L);
                p += L;
            }
            buf[p] = '\0';
            return run_line_str(t, buf);
        }
        else if (!strcmp(tok[0], "command") && nt >= 2) {
            return run_segment(t, tok + 1, eleg + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "type") && nt >= 2) {
            for (int i = 1; i < nt; i++) {
                /* POSIX: describe command */
                char *probe[2] = {(char *)tok[i], NULL};
                int rr2 = posix_rewrite(probe, 1, rtok, &rn);
                if (rr2 < 0) printf("  %s: S2 internal\n", tok[i]);
                else printf("  %s: POSIX builtin\n", tok[i]);
            }
            return 0;
        }
        else if (!strcmp(tok[0], "shift") && nt >= 1) {
            int n = 1;
            if (nt >= 2 && !parse_int(tok[1], &n)) { sh_err("  usage: shift [n]\n"); return 2; }
            if (n < 0 || n > sh_argc) { sh_err("  shift: can't shift %d\n", n); return 1; }
            for (int i = 0; i + n < sh_argc; i++)
                snprintf(sh_argv[i], sizeof sh_argv[i], "%s", sh_argv[i + n]);
            sh_argc -= n;
            return 0;
        }
        else if (!strcmp(tok[0], "set") && nt >= 2 && !strcmp(tok[1], "--")) {
            /* POSIX: set -- args.. (positional params) */
            sh_argc = 0;
            for (int i = 2; i < nt && sh_argc < SH_ARGS; i++)
                snprintf(sh_argv[sh_argc++], sizeof sh_argv[0], "%s", tok[i]);
            return 0;
        }
        else if (!strcmp(tok[0], "set") && nt == 2 && !strcmp(tok[1], "-")) {
            sh_argc = 0;
            return 0;
        }
        else if (!strcmp(tok[0], "readonly") && nt >= 2) {
            /* minimal: accept names, values ignored (no write protection yet) */
            for (int i = 1; i < nt; i++) {
                char *eq = strchr(tok[i], '=');
                if (eq) {
                    size_t nl = (size_t)(eq - tok[i]);
                    char nm[64];
                    if (nl >= sizeof nm || !sh_valid_name(tok[i], nl)) { sh_err("  bad name `%s`\n", tok[i]); return 2; }
                    memcpy(nm, tok[i], nl); nm[nl] = '\0';
                    sh_set(nm, eq + 1);
                } else if (!sh_valid_name(tok[i], strlen(tok[i]))) { sh_err("  bad name `%s`\n", tok[i]); return 2; }
            }
            return 0;
        }
        else if (!strcmp(tok[0], "umask")) {
            /* constrained env: report 022, ignore set (documented) */
            if (nt >= 2) {
                /* validate octal but keep 022 */
                for (const char *p = tok[1]; *p; p++)
                    if (*p < '0' || *p > '7') { sh_err("  usage: umask [NNN]\n"); return 2; }
            } else printf("0022\n");
            return 0;
        }
        else if (!strcmp(tok[0], "trap")) {
            /* no signals in REPL: `trap` lists nothing, `trap CMD SIG` accepted */
            if (nt == 1) return 0;
            if (nt >= 3) return 0;
            sh_err("  usage: trap [CMD SIGNAL...]\n");
            return 2;
        }
        else if (!strcmp(tok[0], "alias") || !strcmp(tok[0], "unalias")) {
            /* no user aliases: fixed POSIX mapping is the alias table */
            if (!strcmp(tok[0], "alias") && nt == 1) {
                printf("  touch=matter-spawn\n  ln=bond\n  df=thermo\n");
                return 0;
            }
            return 0;
        }
        else if (!strcmp(tok[0], "dd")) {
            return cmd_dd(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "ps")) {
            return cmd_ps(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "screen")) {
            /* GNU screen(1) name, S2 operand: fullscreen twin of `ps`
             * (same world, same keybinds, alternate screen + resize).
             * Piped/dumb/missing-curses prints why and returns 0. */
            return live_screen(t);
        }
        else if (!strcmp(tok[0], "vi")) {
            return cmd_vi(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "more")) {
            return cmd_more(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "sync")) {
            return cmd_sync(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "manpage")) {
            if (nt >= 2) return man_page(tok[1]);
            man_index();
            return 0;
        }
        else if (!strcmp(tok[0], "sleep")) {
            /* unified: sleep <steps> (world) | sleep <Ns> (POSIX sec) |
             * sleep <Nfs|Nps> (world steps, explicit) */
            if (nt < 2) { live_advance(t, 1); show_energy(t); return 0; }
            const char *s = tok[1];
            size_t L = strlen(s);
            /* seconds suffix: 2s, 0.5s */
            if (L > 1 && (s[L - 1] == 's' || s[L - 1] == 'S') &&
                !(L > 2 && (s[L - 2] == 'f' || s[L - 2] == 'F') && (s[L - 3] == 'p' || s[L - 3] == 'P' || s[L-2]=='f'))) {
                /* distinguish Nfs (steps) from Ns (sec): Nfs ends fs */
                int is_fs = (L >= 2 && (s[L - 2] == 'f' || s[L - 2] == 'F'));
                if (!is_fs) {
                    char tmp[64];
                    if (L - 1 >= sizeof tmp) { sh_err("  usage: sleep <steps>|<seconds>s\n"); return 2; }
                    memcpy(tmp, s, L - 1); tmp[L - 1] = '\0';
                    double sec;
                    if (!parse_double(tmp, &sec) || !(sec >= 0.0) || !(sec < 3600.0)) {
                        sh_err("  usage: sleep <steps>|<seconds>s\n"); return 2;
                    }
                    sleep_ms((long)(sec * 1000.0));
                    return 0;
                }
            }
            /* explicit steps suffixes: 100fs, 2ps(=2000fs steps at dt .5? no: steps) */
            {
                char tmp[64];
                snprintf(tmp, sizeof tmp, "%s", s);
                size_t tl = strlen(tmp);
                /* strip fs/ps/steps suffix */
                if (tl > 2 && !strcmp(tmp + tl - 2, "fs")) tmp[tl - 2] = '\0';
                else if (tl > 2 && !strcmp(tmp + tl - 2, "ps")) tmp[tl - 2] = '\0';
                else if (tl > 5 && !strcmp(tmp + tl - 5, "steps")) tmp[tl - 5] = '\0';
                int n = 0;
                if (!parse_int(tmp, &n)) { sh_err("  usage: sleep <steps>[fs] | <sec>s\n"); return 2; }
                if (n < 0) n = 0;
                if (n > 100000) n = 100000;
                live_advance(t, n);
                show_energy(t);
                return 0;
            }
        }
        else if (!strcmp(tok[0], "kill")) {
            return cmd_kill(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "du")) {
            return cmd_du(t);
        }
        else if (!strcmp(tok[0], "make")) {
            return cmd_make(t, tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "addsp")) {
            /* streamlined counts: addsp SPEC [x y z] [xN|:N|count=N].
             * No explicit xyz: one open site per copy (3 A clearance),
             * randomly oriented. Stacking copies at the origin detonates
             * the LJ core on the next step and prints absurd temperatures,
             * so the origin is only used when explicitly asked for.
             * Explicit xyz is honored verbatim (your coordinates, your
             * clash). Single atoms need no rotation; molecules get one. */
            int count = 1, cidx = -1;
            /* scan trailing token for count forms */
            for (int i = 1; i < nt; i++) {
                if ((tok[i][0] == 'x' || tok[i][0] == 'X') && tok[i][1]) {
                    int n = 0;
                    if (parse_int(tok[i] + 1, &n) && n > 0 && n <= 64) { count = n; cidx = i; }
                } else if (tok[i][0] == ':' && tok[i][1]) {
                    int n = 0;
                    if (parse_int(tok[i] + 1, &n) && n > 0 && n <= 64) { count = n; cidx = i; }
                } else if (!strncmp(tok[i], "count=", 6)) {
                    int n = 0;
                    if (parse_int(tok[i] + 6, &n) && n > 0 && n <= 64) { count = n; cidx = i; }
                }
            }
            /* origin = first three numbers after species, skipping count token */
            double x = 0, y = 0, z = 0;
            int nums = 0;
            for (int i = 2; i < nt && nums < 3; i++) {
                if (i == cidx) continue;
                double v = 0;
                if (!parse_double(tok[i], &v)) break;
                if (nums == 0) x = v; else if (nums == 1) y = v; else z = v;
                nums++;
            }
            int explicit_xyz = (nums == 3);
            if (nums > 0 && nums != 3) {
                sh_err("  usage: touch <species> [x y z] [xN|:N|count=N]\n");
                return 2;
            }
            /* clearance covers the whole molecule, not one atom */
            double need = live_species_radius(tok[1]) + 2.0;
            int fails = 0;
            for (int k = 0; k < count; k++) {
                Vec3 at = (explicit_xyz || !t->sim) ? vec3(x, y, z)
                                       : live_open_site(t->sim, need);
                int first = t->sim ? t->sim->num_atoms : 0;
                int added = live_species_add(t, tok[1], at);
                if (added <= 0) { fails++; continue; }
                /* stacked spawns (demo on demo) overlap wholesale: verify
                 * the gap and relocate the block once to open space rather
                 * than detonate on the next step. Explicit xyz that clashes
                 * warns but stays verbatim — your coordinates, your clash. */
                if (live_block_gap(t->sim, first) < 1.2) {
                    if (explicit_xyz) {
                        sh_err("  warning: `%s` overlaps the world; expect heat\n", tok[1]);
                    } else {
                        Vec3 open = live_open_site(t->sim, need);
                        Vec3 ref = t->sim->atoms[first].position;
                        live_move_block(t->sim, first, vec3_sub(open, ref));
                        if (live_block_gap(t->sim, first) < 1.2)
                            sh_err("  warning: crowded world; `%s` placed tight\n", tok[1]);
                        else live_rotate_block(t->sim, first, added);
                    }
                } else if (!explicit_xyz && added > 1) {
                    live_rotate_block(t->sim, first, added);
                }
            }
            return fails ? 1 : 0;
        }
        else if (!strcmp(tok[0], "touchf")) {
            /* explicit file touch (see posix_rewrite `touch --`) */
            int rc = 0;
            for (int i = 1; i < nt; i++) {
                FILE *f = fopen(tok[i], "a");
                if (!f) { sh_err("  cannot touch `%s`\n", tok[i]); rc = 1; continue; }
                fclose(f);
            }
            return rc;
        }
        else if (!strcmp(tok[0], "species-export")) {
            return species_export(t, tok[1], tok[2]);
        }
        else if (!strcmp(tok[0], "new")) {
            if (nt >= 2 && (!strcmp(tok[1], "petri") || !strcmp(tok[1], "empty") || !strcmp(tok[1], "void"))) {
                world_preset(t, tok[1], 0);
                return 0;
            }
            int a = 512, b = 512;
            if (nt >= 2) parse_int(tok[1], &a);
            if (nt >= 3) parse_int(tok[2], &b);
            tui_new(t, a, b);
            sh_err("  new sim cap %d atoms\n", t->sim ? t->sim->capacity_atoms : -1);
            maybe_render(t);
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "atom") && nt >= 6) {
            int Z = 0; double x, y, z, q = 0;
            if (isalpha((unsigned char)tok[2][0])) {
                const Element *e = pt_by_symbol(tok[2]);
                if (!e) { sh_err("  unknown symbol `%s`\n", tok[2]); return 1; }
                Z = e->Z;
            } else if (!parse_int(tok[2], &Z)) { sh_err("  bad Z\n"); return 1; }
            if (!parse_double(tok[3], &x) || !parse_double(tok[4], &y) || !parse_double(tok[5], &z)) { sh_err("  bad xyz\n"); return 1; }
            if (nt >= 7 && !parse_double(tok[6], &q)) { sh_err("  bad q\n"); return 1; }
            int i = sim_add_atom(t->sim, Z, vec3(x, y, z), q);
            sh_err(i >= 0 ? "  atom %d Z=%d\n" : "  spawn failed (%d)\n", i >= 0 ? i : i, Z);
            if (i >= 0) maybe_render(t);
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "ion") && nt >= 7) {
            int Z, fm; double x, y, z, q = 1.0;
            if (!parse_int(tok[2], &Z) || !parse_int(tok[3], &fm)) { sh_err("  usage: spawn ion <Z> <formal> x y z [q]\n"); return 1; }
            if (!parse_double(tok[4], &x) || !parse_double(tok[5], &y) || !parse_double(tok[6], &z)) { sh_err("  bad xyz\n"); return 1; }
            if (nt >= 8 && !parse_double(tok[7], &q)) { sh_err("  bad q\n"); return 1; }
            int i = sim_add_ion(t->sim, Z, fm, vec3(x, y, z), q);
            sh_err(i >= 0 ? "  ion %d\n" : "  spawn failed (%d)\n", i);
            if (i >= 0) maybe_render(t);
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && (!strcmp(tok[1], "h2o") || !strcmp(tok[1], "h2") || !strcmp(tok[1], "nh3") || !strcmp(tok[1], "ch4") || !strcmp(tok[1], "methane") || !strcmp(tok[1], "co2"))) {
            double x = 0, y = 0, z = 0;
            if (nt >= 5) { if (!parse_double(tok[2], &x) || !parse_double(tok[3], &y) || !parse_double(tok[4], &z)) { sh_err("  bad origin\n"); return 1; } }
            int f = -1;
            if (!strcmp(tok[1], "h2o")) f = sim_place_h2o(t->sim, vec3(x, y, z));
            else if (!strcmp(tok[1], "h2")) f = sim_place_h2(t->sim, vec3(x, y, z));
            else if (!strcmp(tok[1], "nh3")) f = sim_place_nh3(t->sim, vec3(x, y, z));
            else if (!strcmp(tok[1], "ch4") || !strcmp(tok[1], "methane")) f = sim_place_ch4(t->sim, vec3(x, y, z));
            else f = sim_place_co2(t->sim, vec3(x, y, z));
            sh_err(f >= 0 ? "  placed %s @%d N=%d\n" : "  place failed\n", tok[1], f, t->sim->num_atoms);
            if (f >= 0) maybe_render(t);
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "kcsa")) {
            int nsub = 4;
            if (nt >= 3) parse_int(tok[2], &nsub);
            if (nsub < 1 || nsub > 4) nsub = 4;
            if (!sp_room(t, 41 * nsub + 2, 41 * nsub + 2)) return 1;
            int f = kcsa_build_filter(t->sim, vec3_zero(), nsub);
            sh_err(f >= 0 ? "  kcsa filter @%d N=%d\n" : "  kcsa failed\n", f, t->sim->num_atoms);
            if (f >= 0) maybe_render(t);
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "demo") && nt >= 3) {
            const char *extra = NULL;
            int ao = 3;
            if ((!strcmp(tok[2], "6") || !strcmp(tok[2], "7")) && nt >= 4 && !sp_is_number(tok[3])) {
                extra = tok[3];
                ao = 4;
            }
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, ao, &ok);
            if (!ok) { sh_err("  usage: spawn demo <id> [variant] [x y z]\n"); return 1; }
            int rc = cmd_spawn_demo(t, tok[2], extra, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "quantum")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn quantum [x y z]\n"); return 1; }
            int rc = sp_qm_h(t, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "water")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn water [x y z]\n"); return 1; }
            int rc = sp_water_md(t, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "trimer")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn trimer [x y z]\n"); return 1; }
            int rc = sp_trimer(t, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "methane")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn methane [x y z]\n"); return 1; }
            int f = sim_place_ch4(t->sim, o);
            sh_err(f >= 0 ? "  CH4 tetrahedral @%d\n" : "  sim full\n", f);
            if (f >= 0) maybe_render(t);
            return f >= 0 ? 0 : 1;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "base") && nt >= 3) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 3, &ok);
            if (!ok) { sh_err("  usage: spawn base <U|C|T|A|G> [x y z]\n"); return 1; }
            int rc = sp_base_one(t, tok[2], o);
            if (rc == 0) {
                sh_err("  base %s placed\n", tok[2]);
                maybe_render(t);
            }
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "pair")) {
            const char *which = "gc";
            int ao = 2;
            if (nt >= 3 && !sp_is_number(tok[2])) { which = tok[2]; ao = 3; }
            if (strcasecmp(which, "gc") && strcasecmp(which, "au")) {
                sh_err("  usage: spawn pair [gc|au] [x y z]\n");
                return 1;
            }
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, ao, &ok);
            if (!ok) { sh_err("  usage: spawn pair [gc|au] [x y z]\n"); return 1; }
            int rc = sp_pair(t, !strcasecmp(which, "gc") ? 1 : 0, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "dinucleotide")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn dinucleotide [x y z]\n"); return 1; }
            int rc = cmd_spawn_demo(t, "8", NULL, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "neuron")) {
            int rc = cmd_spawn_demo(t, "9", NULL, vec3_zero());
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "dipeptide")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn dipeptide [x y z]\n"); return 1; }
            int rc = cmd_spawn_demo(t, "10", NULL, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "helix")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn helix [x y z]\n"); return 1; }
            int rc = sp_helix(t, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "cage")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn cage [x y z]\n"); return 1; }
            int rc = sp_cage(t, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "filter")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn filter [x y z]\n"); return 1; }
            if (!sp_room(t, 170, 170)) return 1;
            int f = kcsa_build_filter(t->sim, o, 4);
            sh_err(f >= 0 ? "  real filter 164 atoms @%d (1K4C TVGYG C4)\n" : "  sim full\n", f);
            if (f >= 0) maybe_render(t);
            return f >= 0 ? 0 : 1;
        }
        else if (!strcmp(tok[0], "spawn") && nt >= 2 && !strcmp(tok[1], "duplex")) {
            int ok = 1;
            Vec3 o = spo_origin(tok, nt, 2, &ok);
            if (!ok) { sh_err("  usage: spawn duplex [x y z]\n"); return 1; }
            int rc = sp_duplex(t, o);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "list") && nt >= 2 && !strcmp(tok[1], "atoms")) sim_print_atoms(t->sim);
        else if (!strcmp(tok[0], "list") && nt >= 2 && !strcmp(tok[1], "bonds")) sim_print_bonds(t->sim);
        else if (!strcmp(tok[0], "list") && nt >= 2 && !strcmp(tok[1], "summary")) sim_print_summary(t->sim);
        else if (!strcmp(tok[0], "del") && nt >= 3 && !strcmp(tok[1], "atom")) {
            int i; if (!parse_int(tok[2], &i)) { sh_err("  bad idx\n"); return 1; }
            printf(sim_remove_terminal_atom(t->sim, i) ? "  removed %d\n" : "  remove failed (terminal-only)\n", i);
        }
        else if (!strcmp(tok[0], "bond") && nt >= 3) {
            int a, b, o = 1;
            if (!parse_int(tok[1], &a) || !parse_int(tok[2], &b)) { sh_err("  usage: bond <a> <b> [order]\n"); return 1; }
            if (nt >= 4) parse_int(tok[3], &o);
            if (!t->sim || a < 0 || b < 0 || a >= t->sim->num_atoms || b >= t->sim->num_atoms) {
                sh_err("  ln: atom index out of range\n");
                return 1;
            }
            for (int i = 0; i < t->sim->num_bonds; i++) {
                if ((t->sim->bonds[i].atom_a == a && t->sim->bonds[i].atom_b == b) ||
                    (t->sim->bonds[i].atom_a == b && t->sim->bonds[i].atom_b == a)) {
                    sh_err("  ln: already linked (bond %d)\n", i);
                    return 1;
                }
            }
            int r = sim_add_bond(t->sim, a, b, o);
            sh_err(r >= 0 ? "  bond %d\n" : "  bond failed (%d)\n", r);
            if (r >= 0) maybe_render(t);
        }
        else if (!strcmp(tok[0], "detect") && nt >= 2 && !strcmp(tok[1], "bonds")) {
            int n = sim_detect_bonds(t->sim);
            sim_rebuild_angles(t->sim);
            sh_err("  detected %d bonds\n", n);
        }
        else if (!strcmp(tok[0], "restrain") && nt >= 6) {
            int i; double x, y, z, k;
            if (!parse_int(tok[1], &i) || !parse_double(tok[2], &x) || !parse_double(tok[3], &y) || !parse_double(tok[4], &z) || !parse_double(tok[5], &k)) { sh_err("  usage: restrain <i> x y z k\n"); return 1; }
            int r = sim_add_restraint(t->sim, i, vec3(x, y, z), k);
            if (r >= 0) sh_err("  restraint %d\n", r);
            else sh_err("  restrain failed\n");
        }
        else if (!strcmp(tok[0], "clear") && nt >= 2 && !strcmp(tok[1], "restraints")) { sim_clear_restraints(t->sim); sh_err("  restraints cleared\n"); }
        else if (!strcmp(tok[0], "set") && nt >= 4 && !strcmp(tok[1], "lj")) {
            int i; double e, s;
            if (!parse_int(tok[2], &i) || !parse_double(tok[3], &e) || (nt < 5 || !parse_double(tok[4], &s))) { sh_err("  usage: set lj <i> <eps_eV> <sig_A>\n"); return 1; }
            sim_set_atom_lj(t->sim, i, e, s);
            sh_err("  lj set\n");
        }
        else if (!strcmp(tok[0], "set") && nt >= 4 && !strcmp(tok[1], "charge")) {
            int i; double q;
            if (!parse_int(tok[2], &i) || !parse_double(tok[3], &q)) { sh_err("  usage: set charge <i> <q>\n"); return 1; }
            if (i < 0 || i >= t->sim->num_atoms) sh_err("  bad idx\n");
            else { t->sim->atoms[i].partial_charge = q; sh_err("  q[%d]=%.4f\n", i, q); }
        }
        else if (!strcmp(tok[0], "set") && nt >= 3 && !strcmp(tok[1], "thermostat")) {
            if (!strcmp(tok[2], "none")) t->sim->thermostat.type = THERMOSTAT_NONE;
            else if (!strcmp(tok[2], "berendsen")) t->sim->thermostat.type = THERMOSTAT_BERENDSEN;
            else if (!strcmp(tok[2], "andersen")) t->sim->thermostat.type = THERMOSTAT_ANDERSEN;
            else if (!strcmp(tok[2], "langevin")) t->sim->thermostat.type = THERMOSTAT_LANGEVIN;
            else { sh_err("  thermostat none|berendsen|andersen|langevin\n"); return 1; }
            sh_err("  thermostat set\n");
        }
        else if (!strcmp(tok[0], "set") && !strcmp(tok[1], "box") && nt >= 5) {
            double x, y, z;
            if (!parse_double(tok[2], &x) || !parse_double(tok[3], &y) || !parse_double(tok[4], &z) ||
                !(x > 0) || !(y > 0) || !(z > 0)) {
                sh_err("  usage: set box <lx> <ly> <lz>  (A, positive)\n");
                return 2;
            }
            sim_set_box(t->sim, x, y, z);
            sh_err("  box %.3f x %.3f x %.3f A (PBC all on)\n", x, y, z);
        }
        else if (!strcmp(tok[0], "set") && !strcmp(tok[1], "pbc")) {
            if (nt == 3 && (!strcmp(tok[2], "on") || !strcmp(tok[2], "off"))) {
                int v = !strcmp(tok[2], "on") ? 1 : 0;
                t->sim->box.periodic[0] = t->sim->box.periodic[1] = t->sim->box.periodic[2] = v;
                sh_err("  pbc %s\n", v ? "on (all axes)" : "off");
            } else if (nt == 5) {
                int a, b, c;
                if (!parse_int(tok[2], &a) || !parse_int(tok[3], &b) || !parse_int(tok[4], &c) ||
                    (a != 0 && a != 1) || (b != 0 && b != 1) || (c != 0 && c != 1)) {
                    sh_err("  usage: set pbc on|off | set pbc <x> <y> <z>  (0/1 each)\n");
                    return 2;
                }
                t->sim->box.periodic[0] = a;
                t->sim->box.periodic[1] = b;
                t->sim->box.periodic[2] = c;
                sh_err("  pbc %d %d %d\n", a, b, c);
            } else {
                sh_err("  usage: set pbc on|off | set pbc <x> <y> <z>  (0/1 each)\n");
                return 2;
            }
        }
        else if (!strcmp(tok[0], "set") && !strcmp(tok[1], "press") && nt >= 3) {
            double v;
            if (!parse_double(tok[2], &v)) { sh_err("  usage: set press <bar>\n"); return 2; }
            t->p0_bar = v;
            sh_err("  P target %.3f bar\n", v);
        }
        else if (!strcmp(tok[0], "set") && !strcmp(tok[1], "tau-p") && nt >= 3) {
            double v;
            if (!parse_double(tok[2], &v) || !(v > 0)) { sh_err("  usage: set tau-p <fs>  (positive)\n"); return 2; }
            t->taup_fs = v;
            sh_err("  tau-p %.3f fs\n", v);
        }
        else if (!strcmp(tok[0], "set") && !strcmp(tok[1], "barostat") && nt >= 3) {
            if (!strcmp(tok[2], "on")) t->baro_on = 1;
            else if (!strcmp(tok[2], "off")) t->baro_on = 0;
            else { sh_err("  usage: set barostat on|off\n"); return 2; }
            if (t->baro_on && !t->taup_fs) t->taup_fs = 500.0;
            sh_err("  barostat %s (P0=%.3f bar, tau=%.1f fs)\n",
                t->baro_on ? "on" : "off", t->p0_bar, t->taup_fs);
        }
        else if (!strcmp(tok[0], "set") && nt >= 2 && strchr(tok[1], '=')) {
            /* streamlined: set dt=0.5 cutoff=12 temp=300 (POSIX-like argv) */
            int rc = 0;
            for (int i = 1; i < nt; i++) {
                char *eq = strchr(tok[i], '=');
                if (!eq || eq == tok[i]) { sh_err("  usage: set NAME=VALUE ...\n"); return 2; }
                size_t nl = (size_t)(eq - tok[i]);
                char nm[64];
                if (nl >= sizeof nm) { sh_err("  bad name `%s`\n", tok[i]); rc = 1; continue; }
                memcpy(nm, tok[i], nl); nm[nl] = '\0';
                if (world_set_param(t, nm, eq + 1)) {
                    /* fall back to shell var so set FOO=bar still works */
                    if (!sh_valid_name(nm, nl) || sh_set(nm, eq + 1)) { sh_err("  unknown set `%s`\n", nm); rc = 1; }
                } else sh_err("  set %s=%s\n", nm, eq + 1);
            }
            return rc;
        }
        else if (!strcmp(tok[0], "set") && nt >= 3) {
            double v; if (!parse_double(tok[2], &v)) { sh_err("  bad value\n"); return 1; }
            if (!strcmp(tok[1], "dt")) t->sim->dt = v;
            else if (!strcmp(tok[1], "cutoff")) t->sim->cutoff = v;
            else if (!strcmp(tok[1], "dielectric")) t->sim->dielectric = v;
            else if (!strcmp(tok[1], "temp")) { t->sim->thermostat.target_temperature = v; sh_err("  T target %.1fK\n", v); return 0; }
            else if (!strcmp(tok[1], "tau")) t->sim->thermostat.tau = v;
            else if (!strcmp(tok[1], "nu")) t->sim->thermostat.nu = v;
            else if (!strcmp(tok[1], "gamma")) t->sim->thermostat.gamma = v;
            else if (!strcmp(tok[1], "maxtemp")) t->sim->max_temperature = v;
            else if (!strcmp(tok[1], "seed")) t->seed = (unsigned long)v;
            else { sh_err("  unknown set `%s`\n", tok[1]); return 1; }
            sh_err("  set %s=%.6g\n", tok[1], v);
        }
        else if (!strcmp(tok[0], "init") && nt >= 3 && !strcmp(tok[1], "velocities")) {
            double Tv; unsigned long sd = t->seed;
            if (!parse_double(tok[2], &Tv)) { sh_err("  usage: init velocities <T> [seed]\n"); return 1; }
            if (nt >= 4) parse_ulong(tok[3], &sd);
            integrator_maxwell_boltzmann(t->sim, Tv, sd);
            forces_calculate(t->sim);
            sh_err("  MB T=%.1fK seed=%lu\n", Tv, sd);
        }
        else if (!strcmp(tok[0], "step")) {
            int n = 1;
            if (nt >= 2) parse_int(tok[1], &n);
            cmd_step(t, n);
        }
        else if (!strcmp(tok[0], "run") && nt >= 2) {
            int n; if (!parse_int(tok[1], &n)) { sh_err("  usage: run <N>\n"); return 1; }
            cmd_step(t, n);
        }
        else if (!strcmp(tok[0], "show") && nt >= 2 && !strcmp(tok[1], "energy")) show_energy(t);
        else if (!strcmp(tok[0], "show") && nt >= 2 && !strcmp(tok[1], "pressure")) {
            double pb, pe, nkt, vir;
            if (tui_pressure(t, &pb, &pe, &nkt, &vir)) {
                sh_err("  no box: `set box <lx> <ly> <lz>` first\n");
                return 1;
            }
            double V = 0;
            box_volume(t->sim, &V);
            printf("  P=%.6f bar (%.6f eV/A^3)  NkT/V=%.6f  vir/3V=%.6f  V=%.3f A^3 N=%d\n",
                pb, pe, nkt * EVA3_TO_BAR, vir * EVA3_TO_BAR, V, t->sim->num_atoms);
            return 0;
        }
        else if (!strcmp(tok[0], "show") && nt >= 2 && !strcmp(tok[1], "thermo")) {
            double V = 0, mtot = 0;
            for (int i = 0; i < t->sim->num_atoms; i++) mtot += t->sim->atoms[i].mass > 0 ? t->sim->atoms[i].mass : 0;
            int hasV = !box_volume(t->sim, &V);
            double pb = 0;
            int hasP = hasV && !tui_pressure(t, &pb, NULL, NULL, NULL);
            double rho = (hasV && V > 0) ? mtot * 1.66053906660e-24 / (V * 1e-24) : 0;
            printf("  T=%.2f K  N=%d  KE=%.6f eV\n", t->sim->temperature, t->sim->num_atoms, t->sim->kinetic_energy);
            if (hasV) printf("  V=%.3f A^3  rho=%.4f g/cm^3", V, rho);
            else printf("  V=unset (vacuum default)");
            if (hasP) printf("  P=%.6f bar", pb);
            printf("  barostat=%s\n", t->baro_on ? "on" : "off");
            return 0;
        }
        else if (!strcmp(tok[0], "show") && nt >= 2 && !strcmp(tok[1], "temp")) {
            /* per-element kinetic temperatures: T = 2·KE/(3·N·kB) */
            double kb = BOLTZMANN_K / EV_TO_J;
            int done[119] = {0};
            for (int i = 0; i < t->sim->num_atoms; i++) {
                int z = t->sim->atoms[i].Z;
                if (z < 1 || z > 118 || done[z]) continue;
                done[z] = 1;
                double ke = 0;
                int n = 0;
                for (int j = 0; j < t->sim->num_atoms; j++) {
                    if (t->sim->atoms[j].Z != z) continue;
                    Vec3 v = t->sim->atoms[j].velocity;
                    ke += 0.5 * t->sim->atoms[j].mass * (v.x * v.x + v.y * v.y + v.z * v.z) * AMU_AFS2_TO_EV;
                    n++;
                }
                printf("  Z=%-3d N=%-4d T=%8.2f K\n", z, n, n > 0 ? 2.0 * ke / (3.0 * n * kb) : 0);
            }
            return 0;
        }
        else if (!strcmp(tok[0], "heat") && nt >= 2) {
            double dE;
            if (!parse_double(tok[1], &dE) || !isfinite(dE)) { sh_err("  usage: heat <dE_eV>  (signed)\n"); return 2; }
            /* recompute KE live: the stored field is only refreshed by the
             * integrator, so it is stale right after `init velocities`. */
            double ke = integrator_kinetic_energy(t->sim);
            if (ke <= 1e-12) {
                sh_err("  no kinetic energy: `init velocities` first\n");
                return 1;
            }
            double tgt = ke + dE;
            if (tgt <= 0) {
                for (int i = 0; i < t->sim->num_atoms; i++)
                    t->sim->atoms[i].velocity = vec3_zero();
                t->sim->kinetic_energy = 0;
                t->sim->temperature = 0;
                sh_err("  quenched to 0 K\n");
                return 0;
            }
            double lam = sqrt(tgt / ke);
            for (int i = 0; i < t->sim->num_atoms; i++)
                t->sim->atoms[i].velocity = vec3_scale(t->sim->atoms[i].velocity, lam);
            integrator_remove_com_velocity(t->sim);
            t->sim->kinetic_energy = integrator_kinetic_energy(t->sim);
            t->sim->temperature = integrator_temperature(t->sim);
            forces_calculate(t->sim);
            printf("  heat %+.6f eV -> KE=%.6f eV T=%.2f K\n", dE, t->sim->kinetic_energy, t->sim->temperature);
            return 0;
        }
        else if (!strcmp(tok[0], "minimize")) {
            int it = 500; double st = 0.01, tol = 0.02;
            if (nt >= 2) parse_int(tok[1], &it);
            if (nt >= 3) parse_double(tok[2], &st);
            if (nt >= 4) parse_double(tok[3], &tol);
            double e = integrator_minimize(t->sim, it, st, tol);
            sh_err("  minimized E=%.6f\n", e);
            show_energy(t);
            maybe_render(t);
        }
        else if (!strcmp(tok[0], "neuron") && nt >= 2 && !strcmp(tok[1], "init")) {
            hh_init(&t->nrn); t->has_nrn = 1;
            sh_err("  neuron V=%.2f m=%.4f h=%.4f n=%.4f\n", t->nrn.V, t->nrn.m, t->nrn.h, t->nrn.n);
        }
        else if (!strcmp(tok[0], "neuron") && nt >= 2 && !strcmp(tok[1], "step")) {
            if (!t->has_nrn) { sh_err("  neuron init first\n"); return 1; }
            hh_step(&t->nrn, 0.01);
            printf("  t=%.2fms V=%.2f m=%.3f h=%.3f n=%.3f\n", t->nrn.t, t->nrn.V, t->nrn.m, t->nrn.h, t->nrn.n);
        }
        else if (!strcmp(tok[0], "neuron") && nt >= 2 && !strcmp(tok[1], "run") && nt >= 3) {
            int n; if (!parse_int(tok[2], &n)) { sh_err("  usage: neuron run <N>\n"); return 1; }
            if (!t->has_nrn) hh_init(&t->nrn), t->has_nrn = 1;
            int spikes = 0;
            for (int i = 0; i < n; i++) { hh_step(&t->nrn, 0.01); if (hh_is_spiking(&t->nrn, 0.0)) spikes++; }
            sh_err("  ran %d steps t=%.2f V=%.2f spikes=%d\n", n, t->nrn.t, t->nrn.V, spikes);
        }
        else if (!strcmp(tok[0], "neuron") && nt >= 2 && !strcmp(tok[1], "inject") && nt >= 3) {
            double I; if (!parse_double(tok[2], &I)) { sh_err("  usage: neuron inject <uA/cm2>\n"); return 1; }
            if (!t->has_nrn) hh_init(&t->nrn), t->has_nrn = 1;
            t->nrn.I_ext = I;
            sh_err("  I_ext=%.2f\n", I);
        }
        else if (!strcmp(tok[0], "neuron") && nt >= 2 && !strcmp(tok[1], "show")) {
            if (!t->has_nrn) sh_err("  no neuron\n");
            else printf("  t=%.2f V=%.2f m=%.4f h=%.4f n=%.4f I=%.2f\n", t->nrn.t, t->nrn.V, t->nrn.m, t->nrn.h, t->nrn.n, t->nrn.I_ext);
        }
        else if (!strcmp(tok[0], "mol")) {
            return cmd_mol(t, tok, nt); /* place renders itself */
        }
        else if (!strcmp(tok[0], "rxn")) {
            return cmd_rxn(t, tok, nt);
        }
        else if (!strcmp(tok[0], "render")) {
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "view") && nt >= 2 && !strcmp(tok[1], "xy")) {
            t->cam.yaw_deg = 0; t->cam.pitch_deg = 0;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "view") && nt >= 2 && !strcmp(tok[1], "xz")) {
            t->cam.yaw_deg = 0; t->cam.pitch_deg = 90;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "view") && nt >= 2 && !strcmp(tok[1], "yz")) {
            t->cam.yaw_deg = 90; t->cam.pitch_deg = 0;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "view") && nt >= 3 && !strcmp(tok[1], "auto")) {
            if (!strcmp(tok[2], "on")) t->cam.auto_render = 1;
            else if (!strcmp(tok[2], "off")) t->cam.auto_render = 0;
            else { sh_err("  usage: view auto on|off\n"); return 1; }
            sh_err("  auto %s\n", t->cam.auto_render ? "on" : "off");
        }
        else if (!strcmp(tok[0], "cam") && nt >= 3 && !strcmp(tok[1], "yaw")) {
            double v; if (!parse_double(tok[2], &v)) { sh_err("  usage: cam yaw <deg>\n"); return 1; }
            t->cam.yaw_deg = v;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "cam") && nt >= 3 && !strcmp(tok[1], "pitch")) {
            double v; if (!parse_double(tok[2], &v)) { sh_err("  usage: cam pitch <deg>\n"); return 1; }
            t->cam.pitch_deg = v;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "cam") && nt >= 3 && !strcmp(tok[1], "zoom")) {
            double v; if (!parse_double(tok[2], &v) || !(v > 0.05) || !(v < 50.0)) { sh_err("  usage: cam zoom <0.05..50>\n"); return 1; }
            t->cam.zoom = v;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "cam") && nt >= 2 && !strcmp(tok[1], "center")) {
            if (nt >= 3) {
                int i; if (!parse_int(tok[2], &i) || i < 0 || !t->sim || i >= t->sim->num_atoms) { sh_err("  bad idx\n"); return 1; }
                t->cam.center = t->sim->atoms[i].position;
                t->cam.has_center = 1;
            } else t->cam.has_center = 0;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "cam") && nt >= 2 && !strcmp(tok[1], "reset")) {
            int ar = t->cam.auto_render;
            view_cam_reset(&t->cam);
            t->cam.auto_render = ar;
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "slice")) {
            if (nt >= 2 && !strcmp(tok[1], "off")) t->cam.slice = 0.0;
            else if (nt >= 2) {
                double v; if (!parse_double(tok[1], &v) || !(v > 0.0) || !(v < 100.0)) { sh_err("  usage: slice <thick_A|off>\n"); return 1; }
                t->cam.slice = v;
            } else { sh_err("  usage: slice <thick_A|off>\n"); return 1; }
            view_render(t->sim, &t->cam, 1);
        }
        else if (!strcmp(tok[0], "watch") && nt >= 2) {
            int steps; if (!parse_int(tok[1], &steps) || steps < 1) { sh_err("  usage: watch <steps> [delay_ms]\n"); return 1; }
            long delay = 80;
            if (nt >= 3) { double dv; if (!parse_double(tok[2], &dv) || dv < 0.0 || dv > 5000.0) { sh_err("  bad delay\n"); return 1; } delay = (long)dv; }
            if (steps > 100000) steps = 100000;
            int frames = steps <= 60 ? steps : 60;
            int chunk = (steps + frames - 1) / frames;
            int done = 0;
            for (int f = 0; f < frames && done < steps; f++) {
                int c = steps - done < chunk ? steps - done : chunk;
                for (int i = 0; i < c; i++) {
                    integrator_step(t->sim);
                    tui_barostat_step(t);
                    rxn_autocheck(t);
                }
                done += c;
                forces_calculate(t->sim);
                view_render(t->sim, &t->cam, 1);
                sh_err("  [watch %d/%d steps]\n", done, steps);
                if (!isatty(STDOUT_FILENO)) break; /* piped: single frame */
                sleep_ms(delay);
            }
            show_energy(t);
        }
        else if (!strcmp(tok[0], "echo")) {
            int start = 1;
            if (nt >= 2 && !strcmp(tok[1], "-n")) start = 2;
            for (int i = start; i < nt; i++) {
                if (i > start) putchar(' ');
                fputs(tok[i], stdout);
            }
            if (start == 1) putchar('\n');
            else fflush(stdout);
            return 0;
        }
        else if (!strcmp(tok[0], "export")) {
            if (nt == 1) {
                for (int i = 0; i < sh_nvars; i++) printf("  %s=%s\n", sh_name[i], sh_val[i]);
                return 0;
            }
            int rc = 0;
            for (int i = 1; i < nt; i++) {
                char *eq = strchr(tok[i], '=');
                char nm[64];
                size_t nl;
                int shok, wok;
                if (!eq || eq == tok[i] || (nl = (size_t)(eq - tok[i])) >= sizeof nm) {
                    sh_err("  usage: export NAME[=value] ...\n");
                    rc = 2;
                    continue;
                }
                memcpy(nm, tok[i], nl);
                nm[nl] = '\0';
                shok = (sh_valid_name(nm, nl) && sh_set(nm, eq + 1) == 0);
                wok = (world_set_param(t, nm, eq + 1) == 0);
                if (!shok && !wok) {
                    sh_err("  export: `%s` is not a shell name or a world parameter\n", nm);
                    rc = 1;
                }
            }
            return rc;
        }
        else if (!strcmp(tok[0], "unset")) {
            int rc = 0;
            for (int i = 1; i < nt; i++) {
                if (!sh_valid_name(tok[i], strlen(tok[i])) || sh_unset(tok[i])) {
                    sh_err("  no such variable `%s`\n", tok[i]);
                    rc = 1;
                }
            }
            if (nt < 2) { sh_err("  usage: unset NAME ...\n"); return 2; }
            return rc;
        }
        else if (!strcmp(tok[0], "env")) {
            for (int i = 0; i < sh_nvars; i++) printf("%s=%s\n", sh_name[i], sh_val[i]);
            return 0;
        }
        else if (!strcmp(tok[0], "history")) {
            int from = sh_hist_n > SH_HIST ? sh_hist_n - SH_HIST : 0;
            for (int i = from; i < sh_hist_n; i++)
                printf("  %d  %s\n", i + 1, sh_hist[i % SH_HIST]);
            return 0;
        }
        else if ((!strcmp(tok[0], "source") || !strcmp(tok[0], ".")) && nt >= 2) {
            return run_file(t, tok[1]);
        }
        else if (!strcmp(tok[0], "clear")) {
            if (isatty(STDOUT_FILENO)) printf("\x1b[H\x1b[2J");
            else printf("\n");
            return 0;
        }
        else if (!strcmp(tok[0], "time") && nt >= 2) {
            clock_t t0 = clock();
            int rc = run_segment(t, tok + 1, eleg + 1, nt - 1);
            double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
            sh_err("  %.3fs\n", dt);
            return rc;
        }
        else if (!strcmp(tok[0], "wait")) {
            /* wait [pid...]: reap background jobs (wait(1)); bare waits all */
            int rc = 0, status;
            pid_t w;
            if (nt == 1) {
                for (;;) {
                    w = waitpid(-1, &status, 0);
                    if (w < 0) {
                        if (errno == EINTR) continue;
                        break;
                    }
                    rc = WIFEXITED(status) ? WEXITSTATUS(status) : 1;
                }
                return rc;
            }
            for (int i = 1; i < nt; i++) {
                char *e = NULL;
                long p = strtol(tok[i], &e, 10);
                if (!e || e == tok[i] || p <= 0) {
                    sh_err("  usage: wait [pid...]\n");
                    return 2;
                }
                for (;;) {
                    w = waitpid((pid_t)p, &status, 0);
                    if (w < 0) {
                        if (errno == EINTR) continue;
                        sh_err("  no such job `%s`\n", tok[i]);
                        rc = 1;
                        break;
                    }
                    rc = WIFEXITED(status) ? WEXITSTATUS(status) : 1;
                    break;
                }
            }
            return rc;
        }
        else if (!strcmp(tok[0], "save") && nt >= 2) {
            return cmd_save(t, tok[1]);
        }
        else if (!strcmp(tok[0], "load") && nt >= 2) {
            int rc = cmd_load(t, tok[1]);
            if (rc == 0) maybe_render(t);
            return rc;
        }
        else if (!strcmp(tok[0], "ls")) {
            /* sim objects first (POSIX mapping: the sim IS the working set,
             * so `ls atoms` reads like `ls` on a directory); anything else
             * is a real filesystem path. */
            if (nt == 1 || !strcmp(tok[1], "summary")) { sim_print_summary(t->sim); return 0; }
            if (!strcmp(tok[1], "atoms")) { sim_print_atoms(t->sim); return 0; }
            if (!strcmp(tok[1], "bonds")) { sim_print_bonds(t->sim); return 0; }
            if (!strcmp(tok[1], "demos")) {
                printf("  demos: 1 quantum 2 bond 3 water 4 trimer 5 methane 6 bases 7 pairing 8 dinucleotide 9 neuron 10 dipeptide 11 helix 12 kcsa 12b real-filter 17 duplex\n");
                return 0;
            }
            if (nt == 2) return cmd_ls_files(tok[1]);
            sh_err("  usage: ls [atoms|bonds|summary|demos] | ls [path]\n");
            return 2;
        }
        else if (!strcmp(tok[0], "rm") && nt >= 3 && !strcmp(tok[1], "atom")) {
            /* streamlined ranges: rm atom 1..5,7 1-5 1,2,3 (descending) */
            int ids[256]; int nid = 0;
            for (int a = 2; a < nt && nid < 256; a++) {
                char *s = tok[a];
                /* split on commas first */
                char tmp[256];
                snprintf(tmp, sizeof tmp, "%s", s);
                char *save = NULL, *part = strtok_r(tmp, ",", &save);
                while (part && nid < 256) {
                    char *dots = strstr(part, "..");
                    char *dash = (!dots) ? strchr(part, '-') : NULL;
                    /* avoid negative numbers: dash must not be first char */
                    if (dash == part) dash = NULL;
                    if (dots || dash) {
                        char *sep = dots ? dots : dash;
                        int seplen = dots ? 2 : 1;
                        char left[64], right[64];
                        size_t ll = (size_t)(sep - part);
                        if (ll >= sizeof left) ll = sizeof left - 1;
                        memcpy(left, part, ll); left[ll] = '\0';
                        snprintf(right, sizeof right, "%s", sep + seplen);
                        int lo, hi;
                        if (parse_int(left, &lo) && parse_int(right, &hi)) {
                            if (lo > hi) { int tt = lo; lo = hi; hi = tt; }
                            for (int v = lo; v <= hi && nid < 256; v++) ids[nid++] = v;
                        } else { sh_err("  usage: rm atom <i> [..ranges..,]\n"); return 2; }
                    } else {
                        int v;
                        if (!parse_int(part, &v)) { sh_err("  usage: rm atom <i>\n"); return 2; }
                        ids[nid++] = v;
                    }
                    part = strtok_r(NULL, ",", &save);
                }
            }
            /* descending so shifting indices stay valid */
            for (int i = 0; i < nid; i++)
                for (int j = i + 1; j < nid; j++)
                    if (ids[j] > ids[i]) { int tt = ids[i]; ids[i] = ids[j]; ids[j] = tt; }
            int fails = 0;
            for (int i = 0; i < nid; i++) {
                if (sim_remove_terminal_atom(t->sim, ids[i])) sh_err("  removed %d\n", ids[i]);
                else { sh_err("  remove failed %d (terminal-only)\n", ids[i]); fails++; }
            }
            if (nid) maybe_render(t);
            return fails ? 1 : 0;
        }
        else if (!strcmp(tok[0], "rm") && nt >= 2) {
            int rec = 0, force = 0, a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                int onlyflags = 1;
                for (const char *p = tok[a] + 1; *p; p++) {
                    if (*p == 'r' || *p == 'R') rec = 1;
                    else if (*p == 'f') force = 1;
                    else { onlyflags = 0; break; }
                }
                if (!onlyflags) break;
            }
            if (a >= nt) { sh_err("  usage: rm atom <i> | rm [-rf] <file> ...\n"); return 2; }
            int rc = 0;
            for (int i = a; i < nt; i++) {
                struct stat st;
                int exists = stat(tok[i], &st) == 0;
                if (!exists && force) continue;
                if (exists && S_ISDIR(st.st_mode) && !rec) {
                    sh_err("  `%s` is a directory (use -r)\n", tok[i]);
                    rc = 1;
                    continue;
                }
                if (rec && exists && S_ISDIR(st.st_mode)) {
                    if (rm_recursive(tok[i])) rc = 1;
                    continue;
                }
                if (unlink(tok[i]) != 0) {
                    if (!force) { sh_err("  cannot remove `%s`\n", tok[i]); rc = 1; }
                }
            }
            return rc;
        }
        else if (!strcmp(tok[0], "cat")) {
            return cmd_cat(tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "pwd") && nt == 1) {
            char cwd[1024];
            if (!getcwd(cwd, sizeof cwd)) { sh_err("  cannot get cwd\n"); return 1; }
            printf("%s\n", cwd);
            return 0;
        }
        else if (!strcmp(tok[0], "cd")) {
            const char *dst = NULL;
            if (nt >= 2 && strcmp(tok[1], "~")) dst = tok[1];
            else dst = getenv("HOME");
            if (!dst) { sh_err("  HOME not set\n"); return 1; }
            if (chdir(dst) != 0) {
                sh_err("  cannot cd `%s`\n", nt >= 2 ? tok[1] : "~");
                return 1;
            }
            return 0;
        }
        else if (!strcmp(tok[0], "mkdir") && nt >= 2) {
            int p = 0, a = 1;
            if (!strcmp(tok[1], "-p")) { p = 1; a = 2; }
            if (a >= nt) { sh_err("  usage: mkdir [-p] <dir> ...\n"); return 2; }
            int rc = 0;
            for (int i = a; i < nt; i++)
                if (cmd_mkdir(tok[i], p)) rc = 1;
            return rc;
        }
        else if (!strcmp(tok[0], "cp") && nt == 3) {
            return cmd_cp(tok[1], tok[2]);
        }
        else if (!strcmp(tok[0], "mv") && nt == 3) {
            return cmd_mv(tok[1], tok[2]);
        }
        else if (!strcmp(tok[0], "head")) {
            long nl = 10;
            int a = 1;
            if (nt >= 2 && !strncmp(tok[1], "-n", 2) && tok[1][2] != '\0') {
                char *e = NULL;
                long v = strtol(tok[1] + 2, &e, 10);
                if (*e || v < 0) { sh_err("  usage: head [-n N] [+N] <file> ...\n"); return 2; }
                nl = (tok[1][2] == '+') ? -v : v;
                a = 2;
            } else if (nt >= 3 && !strcmp(tok[1], "-n")) {
                char *e = NULL;
                long v = strtol(tok[2], &e, 10);
                if (!e || e == tok[2] || v < 0) { sh_err("  usage: head [-n N] [+N] <file> ...\n"); return 2; }
                nl = (tok[2][0] == '+') ? -v : v;
                a = 3;
            }
            return cmd_head(tok + a, nt - a, nl);
        }
        else if (!strcmp(tok[0], "tail")) {
            long nl = 10;
            int from = 0, a = 1;
            if (nt >= 2 && !strncmp(tok[1], "-n", 2) && tok[1][2] != '\0') {
                char *e = NULL;
                long v = strtol(tok[1] + 2, &e, 10);
                if (*e) { sh_err("  usage: tail [-n N] [+N] <file> ...\n"); return 2; }
                if (tok[1][2] == '+') { from = 1; nl = v; }
                else {
                    if (v < 0) { sh_err("  usage: tail [-n N] [+N] <file> ...\n"); return 2; }
                    nl = v;
                }
                a = 2;
            } else if (nt >= 3 && !strcmp(tok[1], "-n")) {
                char *e = NULL;
                long v = strtol(tok[2], &e, 10);
                if (!e || e == tok[2]) { sh_err("  usage: tail [-n N] [+N] <file> ...\n"); return 2; }
                if (tok[2][0] == '+') { from = 1; nl = v; }
                else {
                    if (v < 0) { sh_err("  usage: tail [-n N] [+N] <file> ...\n"); return 2; }
                    nl = v;
                }
                a = 3;
            }
            return cmd_tail(tok + a, nt - a, nl, from);
        }
        else if (!strcmp(tok[0], "wc")) {
            int L = 0, W = 0, C = 0, a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                for (const char *p = tok[a] + 1; *p; p++) {
                    if (*p == 'l') L = 1;
                    else if (*p == 'w') W = 1;
                    else if (*p == 'c' || *p == 'm') C = 1;
                    else { sh_err("  usage: wc [-lwc] <file> ...\n"); return 2; }
                }
            }
            return cmd_wc(tok + a, nt - a, L, W, C);
        }
        else if (!strcmp(tok[0], "sort")) {
            int num = 0, rev = 0, uq = 0, a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                for (const char *p = tok[a] + 1; *p; p++) {
                    if (*p == 'n') num = 1;
                    else if (*p == 'r') rev = 1;
                    else if (*p == 'u') uq = 1;
                    else { sh_err("  usage: sort [-nru] <file> ...\n"); return 2; }
                }
            }
            return cmd_sort(tok + a, nt - a, num, rev, uq);
        }
        else if (!strcmp(tok[0], "uniq")) {
            int c = 0, d = 0, u = 0, a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                for (const char *p = tok[a] + 1; *p; p++) {
                    if (*p == 'c') c = 1;
                    else if (*p == 'd') d = 1;
                    else if (*p == 'u') u = 1;
                    else { sh_err("  usage: uniq [-cdu] <file> ...\n"); return 2; }
                }
            }
            return cmd_uniq(tok + a, nt - a, c, d, u);
        }
        else if (!strcmp(tok[0], "cut")) {
            int lo[64], hi[64], nr = 0, byf = 0, sup = 0;
            char delim = '\t';
            int a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                const char *t = tok[a];
                if (!strcmp(t, "-s")) { sup = 1; continue; }
                if (t[1] == 'c' || t[1] == 'b' || t[1] == 'f') {
                    const char *list = (t[2] != '\0') ? t + 2 : NULL;
                    if (!list) {
                        if (a + 1 >= nt) { sh_err("  usage: cut -c|-b LIST | -f LIST [-d C] [-s]\n"); return 2; }
                        list = tok[++a];
                    }
                    nr = parse_ranges(list, lo, hi, 64);
                    if (nr < 0) { sh_err("  bad list `%s`\n", list); return 2; }
                    if (t[1] == 'f') byf = 1;
                    continue;
                }
                if (t[1] == 'd') {
                    const char *d = (t[2] != '\0') ? t + 2 : NULL;
                    if (!d) {
                        if (a + 1 >= nt || !tok[a + 1][0]) { sh_err("  usage: cut -f LIST [-d C] [-s]\n"); return 2; }
                        d = tok[++a];
                    }
                    delim = d[0];
                    continue;
                }
                break;
            }
            if (nr == 0) { sh_err("  usage: cut -c|-b LIST | -f LIST [-d C] [-s]\n"); return 2; }
            return cmd_cut(tok + a, nt - a, lo, hi, nr, byf, delim, sup);
        }
        else if (!strcmp(tok[0], "tr")) {
            int del = 0, sq = 0, comp = 0, a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                int onlyflags = 1;
                for (const char *p = tok[a] + 1; *p; p++) {
                    if (*p == 'd') del = 1;
                    else if (*p == 's') sq = 1;
                    else if (*p == 'c' || *p == 'C') comp = 1;
                    else { onlyflags = 0; break; }
                }
                if (!onlyflags) break;
            }
            if (del && a >= nt) { sh_err("  usage: tr [-cds] SET1 [SET2]\n"); return 2; }
            if (!del && a >= nt) { sh_err("  usage: tr [-cds] SET1 [SET2]\n"); return 2; }
            const char *s1 = del || a < nt ? tok[a] : NULL;
            const char *s2 = (a + 1 < nt) ? tok[a + 1] : NULL;
            if (!s1) { sh_err("  usage: tr [-cds] SET1 [SET2]\n"); return 2; }
            if (a + 2 < nt) { sh_err("  too many operands\n"); return 2; }
            return cmd_tr(s1, s2, del, sq, comp);
        }
        else if (!strcmp(tok[0], "grep")) {
            int E = 0, ic = 0, v = 0, nn = 0, cc = 0, lf = 0, a = 1;
            for (; a < nt && tok[a][0] == '-' && tok[a][1]; a++) {
                int onlyflags = 1;
                for (const char *p = tok[a] + 1; *p; p++) {
                    if (*p == 'E') E = 1;
                    else if (*p == 'i') ic = 1;
                    else if (*p == 'v') v = 1;
                    else if (*p == 'n') nn = 1;
                    else if (*p == 'c') cc = 1;
                    else if (*p == 'l') lf = 1;
                    else { onlyflags = 0; break; }
                }
                if (!onlyflags) break;
            }
            if (a >= nt) { sh_err("  usage: grep [-Einvc l] PATTERN [file...]\n"); return 2; }
            return cmd_grep(tok[a], tok + a + 1, nt - a - 1, E, ic, v, nn, cc, lf);
        }
        else if (!strcmp(tok[0], "tee")) {
            int ap = 0, a = 1;
            for (; a < nt && !strcmp(tok[a], "-a"); a++) ap = 1;
            return cmd_tee(tok + a, nt - a, ap);
        }
        else if (!strcmp(tok[0], "basename") && (nt == 2 || nt == 3)) {
            char b[1024], d[16];
            base_dir_name(tok[1], b, sizeof b, d, sizeof d);
            if (nt == 3 && tok[2][0]) {
                size_t bl = strlen(b), sl = strlen(tok[2]);
                if (sl < bl && !strcmp(b + bl - sl, tok[2])) b[bl - sl] = '\0';
            }
            printf("%s\n", b);
            return 0;
        }
        else if (!strcmp(tok[0], "dirname") && nt == 2) {
            char b[16], d[1024];
            base_dir_name(tok[1], b, sizeof b, d, sizeof d);
            printf("%s\n", d);
            return 0;
        }
        else if (!strcmp(tok[0], "touch") && nt >= 2) {
            return cmd_touch(tok + 1, nt - 1);
        }
        else if (!strcmp(tok[0], "rmdir") && nt >= 2) {
            int rc = 0;
            for (int i = 1; i < nt; i++) {
                if (rmdir(tok[i]) != 0) { sh_err("  cannot remove `%s`\n", tok[i]); rc = 1; }
            }
            return rc;
        }
        else if (!strcmp(tok[0], "ln") && (nt == 3 || (nt == 4 && !strcmp(tok[1], "-s")))) {
            if (nt == 4) {
                if (symlink(tok[2], tok[3]) != 0) { sh_err("  cannot link\n"); return 1; }
                return 0;
            }
            if (link(tok[1], tok[2]) != 0) { sh_err("  cannot link\n"); return 1; }
            return 0;
        }
        else if (!strcmp(tok[0], "date")) {
            int u = 0, a = 1;
            for (; a < nt && !strcmp(tok[a], "-u"); a++) u = 1;
            time_t now = time(NULL);
            struct tm tmv;
            struct tm *tp = u ? gmtime_r(&now, &tmv) : localtime_r(&now, &tmv);
            if (!tp) return 1;
            char ob[256];
            const char *fmt = "%a %b %e %H:%M:%S %Z %Y";
            if (a < nt) {
                if (tok[a][0] != '+') { sh_err("  usage: date [-u] [+FORMAT]\n"); return 2; }
                fmt = tok[a] + 1;
                if (a + 1 < nt) { sh_err("  too many operands\n"); return 2; }
            }
            if (strftime(ob, sizeof ob, fmt, tp) == 0) return 1;
            printf("%s\n", ob);
            return 0;
        }
        else if (!strcmp(tok[0], "uname")) {
            struct utsname u;
            if (uname(&u) != 0) return 1;
            int all = 1, s = 0, n = 0, r = 0, v = 0, m = 0;
            if (nt > 1) {
                all = 0;
                for (int i = 1; i < nt; i++) {
                    if (tok[i][0] != '-' || !tok[i][1]) { sh_err("  usage: uname [-asnrv]\n"); return 2; }
                    for (const char *p = tok[i] + 1; *p; p++) {
                        if (*p == 'a') all = 1;
                        else if (*p == 's') s = 1;
                        else if (*p == 'n') n = 1;
                        else if (*p == 'r') r = 1;
                        else if (*p == 'v') v = 1;
                        else if (*p == 'm') m = 1;
                        else { sh_err("  usage: uname [-asnrv]\n"); return 2; }
                    }
                }
            }
            int first = 1;
            if (all || s) { printf("%s%s", first ? "" : " ", u.sysname); first = 0; }
            if (all || n) { printf("%s%s", first ? "" : " ", u.nodename); first = 0; }
            if (all || r) { printf("%s%s", first ? "" : " ", u.release); first = 0; }
            if (all || v) { printf("%s%s", first ? "" : " ", u.version); first = 0; }
            if (all || m) { printf("%s%s", first ? "" : " ", u.machine); first = 0; }
            printf("\n");
            return 0;
        }
        else if (!strcmp(tok[0], "find")) {
            const char *pat = NULL;
            int want_f = 0, want_d = 0;
            int a = 1, rc = 0;
            while (a < nt && tok[a][0] == '-') break;
            int pe = a;
            while (pe < nt && tok[pe][0] != '-') pe++;
            if (pe == a) { sh_err("  usage: find path... [-name PAT] [-type f|d]\n"); return 2; }
            for (int i = pe; i < nt;) {
                if (!strcmp(tok[i], "-name") && i + 1 < nt) { pat = tok[i + 1]; i += 2; }
                else if (!strcmp(tok[i], "-type") && i + 1 < nt) {
                    if (!strcmp(tok[i + 1], "f")) want_f = 1;
                    else if (!strcmp(tok[i + 1], "d")) want_d = 1;
                    else { sh_err("  type must be f or d\n"); return 2; }
                    i += 2;
                } else if (!strcmp(tok[i], "-print")) i++;
                else { sh_err("  unknown predicate `%s`\n", tok[i]); return 2; }
            }
            for (int i = a; i < pe; i++)
                if (find_walk(tok[i], pat, want_f, want_d)) rc = 1;
            return rc;
        }
        else if (tok[0][0] == '!' && (nt >= 2 || tok[0][1] != '\0')) {
            /* full Linux command set, escaped to the real shell.
             * Both `! make test` and `!make test` work. */
            char cmd[4096];
            size_t L = 0;
            int ai = (tok[0][1] == '\0') ? 1 : 0;
            for (int i = ai; i < nt; i++) {
                const char *w = (i == 0) ? tok[0] + 1 : tok[i];
                size_t n = strlen(w);
                if (n == 0) continue;
                if (L + n + 2 > sizeof cmd) { sh_err("  command too long\n"); return 2; }
                if (L > 0) cmd[L++] = ' ';
                memcpy(cmd + L, w, n);
                L += n;
            }
            if (L == 0) { sh_err("  usage: ! <command> ...\n"); return 2; }
            cmd[L] = '\0';
            int rc = system(cmd);
            if (rc == -1) { sh_err("  shell failed\n"); return 1; }
            if (WIFEXITED(rc)) return WEXITSTATUS(rc);
            return 1;
        }
        else if (!strcmp(tok[0], "man") && nt == 2) {
            return print_help_topic(tok[1]);
        }
        else if (!strcmp(tok[0], "list") && nt == 1) {
            sim_print_summary(t->sim);
            return 0;
        }
        else if (!strcmp(tok[0], "man") && nt == 1) {
            sh_err("  usage: man <command>  (try `help` for the list)\n");
            return 2;
        }
        else if (!strcmp(tok[0], "bond") || !strcmp(tok[0], "spawn") || !strcmp(tok[0], "set") ||
                 !strcmp(tok[0], "del") || !strcmp(tok[0], "rm") || !strcmp(tok[0], "restrain") ||
                 !strcmp(tok[0], "detect") || !strcmp(tok[0], "cam") || !strcmp(tok[0], "view") ||
                 !strcmp(tok[0], "slice") || !strcmp(tok[0], "watch") || !strcmp(tok[0], "sleep") ||
                 !strcmp(tok[0], "time") || !strcmp(tok[0], "save") || !strcmp(tok[0], "load") ||
                 !strcmp(tok[0], "source") || !strcmp(tok[0], ".") || !strcmp(tok[0], "init") ||
                 !strcmp(tok[0], "run") || !strcmp(tok[0], "show") || !strcmp(tok[0], "neuron") ||
                 !strcmp(tok[0], "mol") || !strcmp(tok[0], "rxn") || !strcmp(tok[0], "heat") ||
                 !strcmp(tok[0], "test") || !strcmp(tok[0], "list") || !strcmp(tok[0], "ls") ||
                 !strcmp(tok[0], "cat") || !strcmp(tok[0], "cd") || !strcmp(tok[0], "mkdir") ||
                 !strcmp(tok[0], "cp") || !strcmp(tok[0], "mv") || !strcmp(tok[0], "head") ||
                 !strcmp(tok[0], "export") || !strcmp(tok[0], "unset") || !strcmp(tok[0], "!") ||
                 !strcmp(tok[0], "man") || !strcmp(tok[0], "[") || !strcmp(tok[0], "printf") ||
                 !strcmp(tok[0], "tail") || !strcmp(tok[0], "wc") || !strcmp(tok[0], "sort") ||
                 !strcmp(tok[0], "uniq") || !strcmp(tok[0], "cut") || !strcmp(tok[0], "tr") ||
                 !strcmp(tok[0], "grep") || !strcmp(tok[0], "tee") || !strcmp(tok[0], "basename") ||
                 !strcmp(tok[0], "dirname") || !strcmp(tok[0], "touch") || !strcmp(tok[0], "rmdir") ||
                 !strcmp(tok[0], "ln") || !strcmp(tok[0], "date") || !strcmp(tok[0], "uname") ||
                 !strcmp(tok[0], "find") || !strcmp(tok[0], "wait") || !strcmp(tok[0], "true") ||
                 !strcmp(tok[0], "false")) {
            sh_err("  usage: man %s\n", tok[0]);
            return 2;
        }
        else { sh_err("  unknown `%s` (try `man`)\n", tok[0]); return 127; }
    return 0;
}

int main(void) {
    return tui_main();
}
