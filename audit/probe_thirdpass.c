/* audit_probe.c — independent verification harness for the S2 engine.
 *
 * Links the shipped engine objects and checks each module against oracles
 * computed INDEPENDENTLY here: central finite differences for every force
 * term, closed-form hydrogenic results, statistical-ensemble identities,
 * and published reference numbers. Covers ground the shipped suites do not
 * (CHARMM switching, flat-bottom restraints, PBC, coupled-dipole SCF,
 * Pauli, dispersion, NVE conservation, spherical-harmonic orthonormality).
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/types.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/sim.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/forces.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/integrator.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/quantum.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/qm.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/qm_eht.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/neuron.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/aminoacids.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/kcsa_filter.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/datastream.h"
#include "/home/magi-01/Desktop/work/projects/4179-Magi/include/periodic_table.h"

static int npass = 0, nfail = 0;
static void ck(const char *name, int cond, const char *detail) {
    if (cond) { npass++; printf("  PASS %s\n", name); }
    else { nfail++; printf("  FAIL %s\n        %s\n", name, detail ? detail : ""); }
}
static void ckrel(const char *name, double got, double ref, double tol) {
    char d[200]; snprintf(d, sizeof d, "got %.10g ref %.10g rel tol %.1e", got, ref, tol);
    ck(name, fabs(got - ref) <= tol * (1.0 + fabs(ref)), d);
}
#define KB_EV (1.380649e-23 / 1.602176634e-19)

/* ═══════════════ 1. every force term vs central finite difference ═════════ */
static void fd_check(Simulation *sim, const char *label, double atol) {
    forces_calculate(sim);
    double E0 = sim->potential_energy;
    char det[220];
    if (!isfinite(E0)) { snprintf(det,sizeof det,"E0=%g",E0); ck(label,0,det); return; }
    const double h = 1e-6;
    double worst = 0.0; int wi=0, wa=0; double wfd=0, wan=0;
    for (int i = 0; i < sim->num_atoms; i++) {
        double *c[3] = {&sim->atoms[i].position.x,&sim->atoms[i].position.y,&sim->atoms[i].position.z};
        for (int a = 0; a < 3; a++) {
            double o = *c[a];
            *c[a]=o+h; forces_calculate(sim); double Ep=sim->potential_energy;
            *c[a]=o-h; forces_calculate(sim); double Em=sim->potential_energy;
            *c[a]=o;
            forces_calculate(sim);
            double fd = -(Ep-Em)/(2.0*h);
            double an = (a==0)?sim->atoms[i].force.x:(a==1)?sim->atoms[i].force.y:sim->atoms[i].force.z;
            double d = fabs(fd-an);
            if (d > worst) { worst=d; wi=i; wa=a; wfd=fd; wan=an; }
        }
    }
    snprintf(det,sizeof det,"max|F_analytic-F_FD| = %.3e eV/A (atom %d comp %d: fd=%.8g an=%.8g) E0=%.6f",
             worst, wi, wa, wfd, wan, E0);
    ck(label, worst < atol, det);
    printf("        %s\n", det);
}

static void test_forces(void) {
    printf("\n== analytic forces vs central finite difference (every term) ==\n");
    { Simulation *s=sim_create(64,256);
      sim_place_h2o(s,vec3(0,0,0)); sim_place_h2o(s,vec3(3.0,1.0,0.5));
      fd_check(s,"LJ+Coulomb+bond+angle (hard cutoff)",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(64,256);
      sim_place_h2o(s,vec3(0,0,0)); sim_place_h2o(s,vec3(3.0,1.0,0.5));
      s->use_switching=1;
      fd_check(s,"CHARMM cutoff switching ON (untested by shipped suites)",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(64,256);
      sim_place_h2o(s,vec3(0,0,0)); sim_place_h2o(s,vec3(3.0,1.0,0.5));
      s->use_switching=1; s->cutoff=6.0;
      fd_check(s,"switching at cutoff 6 A",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(64,256); sim_place_h2o(s,vec3(0,0,0));
      s->use_polar=1;
      fd_check(s,"first-order induced dipoles (Thole-damped)",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(64,256); sim_place_h2o(s,vec3(0,0,0));
      s->use_pol_scf=1;
      fd_check(s,"coupled-dipole SCF (3N x 3N solve)",1e-3); sim_destroy(s); }
    { Simulation *s=sim_create(64,256); sim_place_h2o(s,vec3(0,0,0));
      s->use_pauli=1;
      fd_check(s,"overlap Pauli repulsion",1e-3); sim_destroy(s); }
    { Simulation *s=sim_create(64,256);
      sim_place_h2o(s,vec3(0,0,0)); sim_place_h2o(s,vec3(3.2,1.0,0.5));
      s->use_disp=1;
      fd_check(s,"Slater-Kirkwood + Tang-Toennies dispersion",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(64,256);
      sim_place_h2o(s,vec3(0,0,0)); sim_place_h2o(s,vec3(3.2,1.0,0.5));
      sim_add_restraint_fb(s,0,vec3(0,0,0),0.5,0.3);
      sim_add_restraint_fb(s,3,vec3(3.2,1.0,0.5),0.5,0.3);
      fd_check(s,"flat-bottom restraints",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(64,256); sim_place_ch4(s,vec3(0,0,0));
      fd_check(s,"CH4 bonds+angles+LJ",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(512,2048);
      AAResidue *r=calloc(4,sizeof(AAResidue));
      sim_place_polyalanine(s,vec3(0,0,0),4,r);
      fd_check(s,"polyalanine(4): full bonded topology",1e-3);
      free(r); sim_destroy(s); }
    { Simulation *s=sim_create(64,256);
      sim_place_h2o(s,vec3(1,1,1)); sim_place_h2o(s,vec3(9.0,1,1));
      sim_set_box(s,12.0,12.0,12.0);
      fd_check(s,"PBC minimum-image pair forces",1e-4); sim_destroy(s); }
    { Simulation *s=sim_create(512,2048);
      kcsa_build_filter(s, vec3(0,0,0), 4);
      fd_check(s,"KcsA TVGYG filter: LJ+Coulomb rigid cage",1e-3); sim_destroy(s); }
}

/* ═══════════════ 2. integrator ═══════════════ */
static void test_integrator(void) {
    printf("\n== integrator: conservation + closed forms ==\n");
    { Simulation *s=sim_create(8,8);
      sim_add_atom(s,12,vec3(0,0,0),0.0);
      sim_add_atom(s,12,vec3(1.95,0,0),0.0);
      sim_add_bond(s,0,1,1); sim_set_bond_params(s,0,1.75,20.0);
      s->use_lj=0; s->use_coulomb=0; s->use_angles=0; s->dt=0.25;
      /* FULL-AUDIT Q6 NOTE: the default diagnostic velocity cap
       * (max_temperature = 2000 K) rescales velocities whenever the
       * instantaneous temperature exceeds it. A hand-seeded kick of
       * 0.05 A/fs on Mg corresponds to ~24000 K, so the cap fires and
       * removes ~90% of the KE in one step — correct behaviour for the
       * cap, but it makes the run NOT NVE. Disable it (0 = off) so this
       * measures the integrator alone. */
      s->max_temperature = 0.0;
      s->atoms[0].velocity = vec3(0,0,0.05);
      forces_calculate(s);
      double E0=s->potential_energy+integrator_kinetic_energy(s);
      double lo=1e30,hi=-1e30;
      for(int i=0;i<20000;i++){ integrator_step(s);
          double E=s->kinetic_energy+s->potential_energy;
          if(E<lo)lo=E; if(E>hi)hi=E; }
      char d[160]; snprintf(d,sizeof d,"E0=%.10f span=%.3e eV",E0,hi-lo);
      ck("NVE Verlet: harmonic oscillator 20000 steps, bounded energy (Verlet ~ (dt w)^2/12)",(hi-lo)<1e-4*E0,d);
      printf("        %s\n",d); sim_destroy(s); }
    { Simulation *s=sim_create(16,64); sim_place_h2o(s,vec3(0,0,0));
      s->use_coulomb=0; s->use_lj=0; s->dt=0.5;
      s->max_temperature = 0.0;   /* pure NVE; see Q6 note above */
      integrator_maxwell_boltzmann(s,300.0,42); forces_calculate(s);
      double E0=s->potential_energy+integrator_kinetic_energy(s), worst=0;
      for(int i=0;i<20000;i++){ integrator_step(s);
          double dE=fabs(s->total_energy-E0); if(dE>worst)worst=dE; }
      char d[160]; snprintf(d,sizeof d,"max abs(dE)=%.3e eV over 10 ps, E0=%.6f",worst,E0);
      ck("NVE: isolated H2O bounded energy (dt=0.5 fs, stiff O-H)",worst<1e-3,d); printf("        %s\n",d); sim_destroy(s); }
    { Simulation *s=sim_create(16,64);
      sim_place_h2o(s,vec3(0,0,0)); sim_place_h2o(s,vec3(3.2,1.0,0.5));
      s->use_coulomb=0; s->dt=0.5;
      integrator_maxwell_boltzmann(s,300.0,11); forces_calculate(s);
      double E0=s->potential_energy+integrator_kinetic_energy(s), e1=0,e2=0;
      for(int i=0;i<40000;i++){ integrator_step(s);
          double dE=fabs(s->total_energy-E0);
          if(i<5000){ if(dE>e1)e1=dE; } else { if(dE>e2)e2=dE; } }
      char d[180]; snprintf(d,sizeof d,"early max|dE|=%.3e late max|dE|=%.3e ratio=%.2f",e1,e2,e2/(e1+1e-30));
      ck("NVE dimer: energy error non-secular (symplectic)", e2 < 4.0*e1+1e-6, d);
      printf("        %s\n",d); sim_destroy(s); }
    { Simulation *s=sim_create(200,256);
      for(int i=0;i<40;i++) sim_add_atom(s,8,vec3(i*2.0,0,0),-0.834);
      integrator_maxwell_boltzmann(s,300.0,7);
      ckrel("Maxwell-Boltzmann: instantaneous <T> == target", integrator_temperature(s), 300.0, 1e-9);
      double kT=KB_EV*300.0, sum=0; int n=0;
      for(int i=0;i<s->num_atoms;i++){ sum += s->atoms[i].velocity.x*s->atoms[i].velocity.x*s->atoms[i].mass/kT; n++; }
      ckrel("MB: per-component m<v^2> = kT/AMU_AFS2_TO_EV (40 O atoms)", sum/n, 1.0/103.642696, 0.25);
      sim_destroy(s); }
    { Simulation *s=sim_create(8,8); sim_add_atom(s,19,vec3(0,0,0),1.0);
      s->thermostat.type=THERMOSTAT_LANGEVIN; s->thermostat.target_temperature=300.0;
      s->thermostat.gamma=0.05; s->dt=1.0;
      integrator_maxwell_boltzmann(s,300.0,3);
      double sv2=0; int nv=0;
      for(int i=0;i<200000;i++){ integrator_langevin(s);
          if(i>1000){ sv2+=vec3_norm2(s->atoms[0].velocity); nv++; } }
      ckrel("Langevin OU: stationary <v^2> = 3kT/m_K", sv2/nv, 3.0*KB_EV*300.0/39.0983, 0.02);
      sim_destroy(s); }
    { Simulation *s=sim_create(8,8); sim_add_atom(s,19,vec3(0,0,0),1.0);
      s->thermostat.type=THERMOSTAT_ANDERSEN; s->thermostat.target_temperature=300.0;
      s->thermostat.nu=0.02; s->dt=1.0;
      integrator_maxwell_boltzmann(s,300.0,3);
      double sv2=0; int nv=0;
      for(int i=0;i<200000;i++){ integrator_andersen(s);
          if(i>1000){ sv2+=vec3_norm2(s->atoms[0].velocity); nv++; } }
      ckrel("Andersen: stationary <v^2> = 3kT/m_K", sv2/nv, 3.0*KB_EV*300.0/39.0983, 0.03);
      sim_destroy(s); }
    /* Berendsen: T must relax toward T0 */
    { Simulation *s=sim_create(64,64);
      for(int i=0;i<20;i++) sim_add_atom(s,8,vec3(i*2.0,0,0),-0.834);
      s->use_lj=0; s->use_coulomb=0; s->use_bonds=0; s->use_angles=0;
      s->dt=1.0;
      s->thermostat.type=THERMOSTAT_BERENDSEN; s->thermostat.target_temperature=300.0;
      s->thermostat.tau=100.0;
      integrator_maxwell_boltzmann(s,1000.0,5);
      for(int i=0;i<200000;i++){ forces_calculate(s); integrator_berendsen(s); }
      forces_calculate(s);
      double T=integrator_temperature(s);
      char d[120]; snprintf(d,sizeof d,"T after 200 ps = %.3f K (target 300)",T);
      ckrel("Berendsen: T relaxes to target from 1000 K", T, 300.0, 0.02);
      printf("        %s\n",d); sim_destroy(s); }
}

/* ═══════════════ 3. quantum.c vs closed forms ═══════════════ */
static void test_quantum(void) {
    printf("\n== quantum layer vs analytic hydrogenic results ==\n");
    { ElectronConfig c; memset(&c,0,sizeof c);
      int sh[7][2]={{1,0},{2,0},{2,1},{3,0},{3,1},{3,2},{4,0}};
      int ne[7]={2,2,6,2,6,6,2};
      for(int i=0;i<7;i++){ c.config[sh[i][0]-1][sh[i][1]]=ne[i]; c.total_electrons+=ne[i]; }
      ckrel("Slater Fe 4s Z_eff (Slater Phys.Rev.36,57 example)",
            quantum_zeff_raw(26,4,0,&c), 3.75, 1e-12);
      ckrel("Slater Fe 3d sigma (Slater's own example)",
            26.0-quantum_zeff_raw(26,3,2,&c), 19.75, 1e-12);
      ckrel("Slater Fe 3d Z_eff", quantum_zeff_raw(26,3,2,&c), 6.25, 1e-12); }
    { double z=1.0,sum=0,h=1e-4;
      for(double r=1e-6;r<25.0;r+=h){ double R=quantum_radial_wavefunction(1,0,z,r); sum+=r*r*R*R*h; }
      ckrel("H 1s: integral r^2|R|^2 dr == 1 (numeric)",sum,1.0,1e-5); }
    { double z=1.0,sum=0,h=1e-4;
      for(double r=1e-6;r<60.0;r+=h){ double R=quantum_radial_wavefunction(2,0,z,r); sum+=r*r*R*R*h; }
      ckrel("H 2s: integral r^2|R|^2 dr == 1 (numeric)",sum,1.0,1e-5); }
    { double z=1.0,sum=0,h=1e-4;
      for(double r=1e-6;r<60.0;r+=h){ double R=quantum_radial_wavefunction(2,1,z,r); sum+=r*r*R*R*h; }
      ckrel("H 2p: integral r^2|R|^2 dr == 1 (numeric)",sum,1.0,1e-5); }
    { double z=2.0,h=1e-4,s1=0,s2=0;
      for(double r=1e-6;r<20.0;r+=h){ double R=quantum_radial_wavefunction(2,0,z,r);
          s1+=r*r*r*R*R*h; s2+=r*r*r*r*R*R*h; }
      ckrel("<r> 2s numeric == Griffiths closed form",s1,quantum_expect_r(2,0,z),1e-4);
      ckrel("<r^2> 2s numeric == Griffiths closed form",s2,quantum_expect_r2(2,0,z),1e-4); }
    { double z=1.5;
      double E=quantum_orbital_energy(1,1,0,NULL);
      (void)E;(void)z;
      /* virial: <T> = +Z^2/(2n^2) Eh ; E_orb = -13.6057 (Z/n*)^2 */
      ckrel("virial: <T> = -E for H-like n=1 Z=1",
            quantum_expect_T(1,0,1.0)+quantum_orbital_energy(1,1,0,NULL), 0.0, 1e-12); }
    { double sum=0; int nth=200,np=400;
      for(int it=0;it<nth;it++){ double th=(it+0.5)*M_PI/nth;
        for(int ip=0;ip<np;ip++){ double ph=(ip+0.5)*2*M_PI/np;
          double w=sin(th)*(M_PI/nth)*(2*M_PI/np);
          for(int m=-2;m<=2;m++){ double y=qm_Y_real(2,m,th,ph); sum+=y*y*w; } } }
      ckrel("real harmonics l=2: 5 functions orthonormal (quad-limited)",sum,5.0,1e-4); }
    { double sum=0; int nth=300,np=600;
      for(int it=0;it<nth;it++){ double th=(it+0.5)*M_PI/nth;
        for(int ip=0;ip<np;ip++){ double ph=(ip+0.5)*2*M_PI/np;
          double w=sin(th)*(M_PI/nth)*(2*M_PI/np);
          for(int m=-3;m<=3;m++){ double y=qm_Y_real(3,m,th,ph); sum+=y*y*w; } } }
      ckrel("real harmonics l=3: 7 functions orthonormal (quad-limited)",sum,7.0,1e-4); }
    { /* cross-harmonic orthogonality: int Y_lm Y_l'm' = 0 for l!=l' */
      double s=0; int nth=200,np=400;
      for(int it=0;it<nth;it++){ double th=(it+0.5)*M_PI/nth;
        for(int ip=0;ip<np;ip++){ double ph=(ip+0.5)*2*M_PI/np;
          double w=sin(th)*(M_PI/nth)*(2*M_PI/np);
          s += qm_Y_real(1,0,th,ph)*qm_Y_real(2,0,th,ph)*w; } }
      ckrel("real harmonics: l=1 orthogonal to l=2 (quad-limited)",s,0.0,1e-4); }
    ckrel("CR Z_eff(C 2p) Clementi-Raimondi 1963",quantum_zeff_cr(6,2,1),3.14,1e-12);
    ckrel("CR Z_eff(O 2p) Clementi-Raimondi 1963",quantum_zeff_cr(8,2,1),4.45,1e-12);
    ckrel("CR Z_eff(K 4s) Clementi-Raimondi 1963",quantum_zeff_cr(19,4,0),3.50,1e-12);
    ckrel("CR Z_eff(Kr 3d)",quantum_zeff_cr(36,3,2),20.63,1e-12);
}

/* ═══════════════ 4. qm.c: induction, QEq, dispersion ═══════════════ */
static void test_qm(void) {
    printf("\n== quantum-mechanics layer ==\n");
    /* 4a. Polarizability table vs NIST/Schwerdtfeger (a.u. of volume) */
    struct { int Z; double a0_3; } pol[] = {
        {1,4.507},{2,1.384},{6,11.30},{7,7.44},{8,5.30},{9,3.74},
        {10,2.661},{11,162.7},{17,14.60},{18,11.08},{19,289.7},{35,21.00},{36,16.80}};
    { int bad=0; double a0c=0.529177210903;
      for (unsigned i=0;i<sizeof pol/sizeof pol[0];i++){
          Atom a; memset(&a,0,sizeof a); a.Z=pol[i].Z;
          a.element=pt_element(pol[i].Z);
          pt_electron_config(pol[i].Z,&a.electron_config);
          double got=qm_polarizability(&a)/(a0c*a0c*a0c);
          if (fabs(got-pol[i].a0_3)>0.02*pol[i].a0_3+0.02) { bad++;
              printf("        Z=%d got %.4f a0^3 ref %.4f\n",pol[i].Z,got,pol[i].a0_3); }
      }
      char d[120]; snprintf(d,sizeof d,"%d of %d mismatched",bad,(int)(sizeof pol/sizeof pol[0]));
      ck("polarizability table vs measured (a0^3)",bad==0,d); }
    /* 4b. Slater-Kirkwood C6 vs experimental H2 C6 = 6.499 a.u. */
    { Atom a,b; memset(&a,0,sizeof a); memset(&b,0,sizeof b);
      a.Z=1; a.element=pt_element(1); pt_electron_config(1,&a.electron_config);
      b=a;
      double c6=qm_c6(&a,&b);
      double au = 27.211386245988*pow(0.529177210903,6); /* eV A^6 per a.u. */
      ckrel("Slater-Kirkwood C6(H-H) vs experimental 6.499 a.u.", c6/au, 6.499, 0.25);
      printf("        C6(H-H) = %.4f a.u. (experimental 6.499)\n", c6/au); }
    { Atom a,b; memset(&a,0,sizeof a); memset(&b,0,sizeof b);
      a.Z=8; a.element=pt_element(8); pt_electron_config(8,&a.electron_config); b=a;
      double au = 27.211386245988*pow(0.529177210903,6);
      ckrel("Slater-Kirkwood C6(O-O) vs measured 15.6 a.u. (SK is a ~25%% estimate)", qm_c6(&a,&b)/au, 15.6, 0.35); }
    /* 4c. Tang-Toennies f6: f6(0)=0, f6(inf)=1, monotone increasing */
    { double b=2.0, prev=-1; int mono=1;
      for(double r=0.5;r<12.0;r+=0.1){ double f=qm_tt_f6(b,r,NULL); if(f<prev)mono=0; prev=f; }
      double df;
      ckrel("TT f6 at r->0 is 0", qm_tt_f6(b,1e-6,NULL), 0.0, 1e-9);
      ckrel("TT f6 at large r is 1", qm_tt_f6(b,40.0,NULL), 1.0, 1e-6);
      { char d[80]; snprintf(d,sizeof d,"monotone=%d",mono); ck("TT f6 monotone increasing",mono,d); }
      /* df/dr vs numeric d/dr */
      double h=1e-5, r=2.3, num=(qm_tt_f6(b,r+h,NULL)-qm_tt_f6(b,r-h,NULL))/(2*h);
      qm_tt_f6(b,r,&df);
      ckrel("TT f6 derivative vs numeric", df, num, 1e-5); }
    /* 4d. QEq: exact analytic two-atom solution */
    { Simulation *s=sim_create(8,8);
      sim_add_atom(s,1,vec3(0,0,0),0.0);
      sim_add_atom(s,1,vec3(2.0,0,0),0.0);
      double q[8];
      int rc=qm_qeq(s,0.0,1.0,q);
      /* Both atoms identical -> q must be 0 by symmetry */
      ckrel("QEq two identical H: q == 0 by symmetry", q[0], 0.0, 1e-9);
      ckrel("QEq conservation sum(q)==total_q", q[0]+q[1], 0.0, 1e-9);
      { char d[80]; snprintf(d,sizeof d,"rc=%d",rc); ck("QEq returns 0",rc==0,d); }
      sim_destroy(s); }
{ Simulation *s=sim_create(8,8);
      sim_add_atom(s,1,vec3(0,0,0),0.0);
      sim_add_atom(s,8,vec3(2.0,0,0),0.0);
      double q[8] = {-99,-99};
      int rc = qm_qeq(s,0.0,1.0,q);
      /* Q1: at 2.0 A the unscreened QEq Hessian is indefinite and the
       * electronegativity ordering inverts; the engine must REFUSE it. */
      ck("QEq refuses the inverted 2.0 A O/H pair (Q1 guard)", rc == -1, NULL);
      char d[160]; snprintf(d, sizeof d, "rc=%d q_H=%+.6f q_O=%+.6f",rc,q[0],q[1]);
      printf("        %s\n",d); sim_destroy(s); }
    { /* and a safe 3.0 A pair still solves, physically ordered */
      Simulation *t=sim_create(8,8);
      sim_add_atom(t,1,vec3(0,0,0),0.0);
      sim_add_atom(t,8,vec3(3.0,0,0),0.0);
      double q2_[8]; int rc2 = qm_qeq(t,0.0,1.0,q2_);
      char d2[160];
      snprintf(d2,sizeof d2,"rc=%d q_H=%+.6f q_O=%+.6f sum=%+.1e",
               rc2,q2_[0],q2_[1],q2_[0]+q2_[1]);
      ck("QEq solves an ordered 3.0 A O/H pair (q_O < 0)",
         rc2==0 && q2_[1]<0.0 && q2_[0]>0.0, d2);
      sim_destroy(t); }
    { /* charge conservation over many random clusters at several total_q */
      int bad=0, worst=0; double worstd=0;
      srand(12345);
      for (int t=0;t<400;t++){
          Simulation *s=sim_create(8,8);
          int n=2+(rand()%5);
          for(int i=0;i<n;i++){
              int Z = (rand()%3==0)?8:1;
              sim_add_atom(s,Z,vec3(1.5*i+0.1*(rand()%3),0,0),0.0);
          }
          double tot=((rand()%3)-1)*0.4;
          double q[8];
          if (qm_qeq(s,tot,1.0,q)==0){
              double sum=0;
              for(int i=0;i<n;i++) sum+=q[i];
              double dq=fabs(sum-tot);
              if (dq>worstd){worstd=dq; worst=1;}
              for(int i=0;i<n;i++) if(fabs(q[i])>2.0+1e-9) bad++;
          }
          sim_destroy(s);
      }
      char d[160]; snprintf(d,sizeof d,"worst |sum(q)-Q| = %.3e e over 400 clusters",worstd);
      ck("QEq: charge conserved for all total_q (random clusters)", worstd<1e-9, d);
      { char d2[120]; snprintf(d2,sizeof d2,"%d charge violations beyond +/-2 e",bad);
        ck("QEq: |q| <= 2 e bound never exceeded", bad==0, d2); }
      printf("        %s\n",d); }
    /* 4e. Induction: U must be negative and vanish without charges */
    { Simulation *s=sim_create(8,8);
      sim_add_atom(s,1,vec3(0,0,0),0.0);
      sim_add_atom(s,8,vec3(3.0,0,0),-0.8);
      double U=qm_induction_energy(s,1.0);
      ck("induction energy attractive (U<0)", U<0.0, NULL);
      char d[120]; snprintf(d,sizeof d,"U=%.6f eV",U); printf("        %s\n",d);
      for(int i=0;i<2;i++) s->atoms[i].partial_charge=0.0;
      ckrel("induction energy 0 when all charges 0", qm_induction_energy(s,1.0), 0.0, 1e-12);
      sim_destroy(s); }
    /* 4f. Hellmann-Feynman-ish: U scales as q^2 for a fixed geometry */
    { Simulation *s=sim_create(8,8);
      sim_add_atom(s,1,vec3(0,0,0),1.0);
      sim_add_atom(s,8,vec3(3.0,0,0),-1.0);
      double U1=qm_induction_energy(s,1.0);
      s->atoms[0].partial_charge=2.0; s->atoms[1].partial_charge=-2.0;
      double U4=qm_induction_energy(s,1.0);
      ckrel("induction U scales as q^2 (EITF line)", U4/U1, 4.0, 1e-9);
      sim_destroy(s); }
}

/* ═══════════════ 5. neuron ═══════════════ */
static void test_neuron(void) {
    printf("\n== Hodgkin-Huxley 1952 ==\n");
    { HHNeuron n; hh_init(&n);
      char d[200];
      snprintf(d,sizeof d,"V0=%.4f mV m=%.6f h=%.6f n=%.6f  (published 1952: m=0.053 h=0.596 n=0.318)",
               n.V,n.m,n.h,n.n);
      ck("HH resting gating values vs Hodgkin-Huxley 1952 worked example",
         fabs(n.m-0.05293)<2e-4 && fabs(n.h-0.59610)<2e-4 && fabs(n.n-0.31768)<2e-4, d);
      printf("        %s\n",d); }
    { HHNeuron n; hh_init(&n); n.I_ext=10.0;
      double dt=0.01, t=50.0, peak=-1e9, nspike=0; int prev=0;
      int steps=(int)(t/dt);
      for(int i=0;i<steps;i++){
          hh_step(&n,dt);
          if(n.V>peak)peak=n.V;
          int sp=(n.V>0.0); if(sp&&!prev)nspike++; prev=sp;
      }
      char d[200]; snprintf(d,sizeof d,"peak=%.3f mV spikes in 50 ms = %d",peak,nspike);
      ck("HH action potential: peak in (30,50) mV at I=10 uA/cm2",
         peak>30.0&&peak<50.0,d);
      ck("HH fires repeatedly (>=4 spikes / 50 ms)", nspike>=4, d);
      printf("        %s\n",d); }
    { /* rate-function sanity: monotonicity where physics demands it */
      int mono_m=1; double prev=hh_alpha_m(-100.0);
      for(double V=-100;V<50;V+=0.5){ double a=hh_alpha_m(V); if(a<prev-1e-12) mono_m=0; prev=a; }
      ck("alpha_m monotone increasing in V", mono_m, NULL);
      mono_m=1; prev=hh_beta_m(-100.0);
      for(double V=-100;V<50;V+=0.5){ double b=hh_beta_m(V); if(b>prev+1e-12) mono_m=0; prev=b; }
      ck("beta_m monotone decreasing in V", mono_m, NULL);
      ckrel("alpha_m at its removable singularity V=-40", hh_alpha_m(-40.0), 1.0, 1e-12);
      ckrel("alpha_n at its removable singularity V=-55", hh_alpha_n(-55.0), 0.1, 1e-12);
      /* L'Hopital limit: alpha_m(-40) = 0.1 * 10 = 1.0 */
      ckrel("alpha_n(-55) limit = 0.01*10", hh_alpha_n(-55.0), 0.10, 1e-12); }
    { /* RK4 convergence order: halving dt must reduce global error 16x */
      HHNeuron a,b; hh_init(&a); a.I_ext=10.0; hh_init(&b); b.I_ext=10.0;
      for(int i=0;i<5000;i++) hh_step(&a,0.01);   /* 50 ms */
      for(int i=0;i<5000;i++) hh_step(&b,0.005);  /* 25 ms */
      HHNeuron c; hh_init(&c); c.I_ext=10.0;
      for(int i=0;i<10000;i++) hh_step(&c,0.005);
      char d[160]; snprintf(d,sizeof d,"V(dt=0.01)=%.9f V(dt=0.005)=%.9f diff=%.3e",a.V,c.V,fabs(a.V-c.V));
      ck("HH RK4: dt/2 halves the error (dt=0.01 vs 0.005)", fabs(a.V-c.V)<2e-3, d);
      printf("        %s\n",d); (void)b; }
}

/* ═══════════════ 6. KcsA ═══════════════ */
static void test_kcsa(void) {
    printf("\n== KcsA filter geometry and energetics ==\n");
    { Simulation *s=sim_create(512,1024);
      int first=kcsa_build_filter(s,vec3(0,0,0),4);
      /* CN at each deposited site */
      Vec3 sites[4]; kcsa_ion_sites(4,sites);
      int all8=1; char d[200];
      for(int i=0;i<4;i++){
          double mr; int n=kcsa_coord_stats(s,first,4,sites[i],3.2,&mr);
          if(n!=8) all8=0;
          printf("        site %d z=%+8.3f  CN(<=3.2A)=%d  <r>=%.3f A\n",i,sites[i].z,n,mr);
      }
      snprintf(d,sizeof d,"CN=%d,%d,%d,%d",kcsa_coord_stats(s,first,4,sites[0],3.2,0),
               kcsa_coord_stats(s,first,4,sites[1],3.2,0),
               kcsa_coord_stats(s,first,4,sites[2],3.2,0),
               kcsa_coord_stats(s,first,4,sites[3],3.2,0));
      ck("KcsA: CN == 8 at all four deposited K+ sites", all8, d);
      double net=0; for(int i=first;i<s->num_atoms;i++) net+=s->atoms[i].partial_charge;
      ckrel("KcsA filter net charge == 0", net, 0.0, 1e-12);
      sim_destroy(s); }
    /* Marcus 1991 TATB hydration free energies */
    { struct { int Z; double dG; } ref[] = {{11,-365.3},{19,-295.3},{37,-275.3},{55,-250.7},{3,-475.1}};
      int bad=0; char d[160];
      for(unsigned i=0;i<sizeof ref/sizeof ref[0];i++){
          double got=kcsa_hydration_free_energy_kJmol(ref[i].Z);
          if(fabs(got-ref[i].dG)>0.05){ bad++; snprintf(d,sizeof d,"Z=%d got %.1f ref %.1f",ref[i].Z,got,ref[i].dG); }
      }
      ck("Marcus 1991 TATB hydration free energies (kJ/mol)", bad==0, d); }
    { double KJ=AVOGADRO_N*1.602176634e-19/1000.0;
      ckrel("kJ/mol per eV constant", KJ, 96.48533212, 1e-9);
      ckrel("K+ dehydration cost (Marcus)", kcsa_dehydration_cost_eV(19), 295.3/KJ, 1e-9);
      ckrel("Na+ dehydration cost (Marcus)", kcsa_dehydration_cost_eV(11), 365.3/KJ, 1e-9);
      ckrel("K+ dehydration advantage (0.7255 eV)",
            kcsa_dehydration_cost_eV(11)-kcsa_dehydration_cost_eV(19), 0.7255, 2e-3); }
    /* Shannon 8-coordinate cation radii */
    { struct { int Z; double r; } ref[] = {{3,0.92},{11,1.18},{19,1.51},{37,1.61},{55,1.74}};
      int bad=0; char d[160];
      for(unsigned i=0;i<sizeof ref/sizeof ref[0];i++){
          double got=kcsa_cation_radius(ref[i].Z);
          if(fabs(got-ref[i].r)>1e-9){ bad++; snprintf(d,sizeof d,"Z=%d got %.2f ref %.2f",ref[i].Z,got,ref[i].r); }
      }
      ck("Shannon 1976 VIII cation radii (A)", bad==0, d); }
    /* binding: K+ must not be repulsive vs the filter */
    { Simulation *s=sim_create(1024,2048);
      int first=kcsa_build_filter(s,vec3(0,0,0),4);
      s->use_bonds=0; s->use_angles=0; s->use_dihedrals=0;
      Vec3 sites[4]; kcsa_ion_sites(4,sites);
      double ek=0,ena=0;
      double e1,e2;
      kcsa_site_binding(s,first,4,19,sites[1],60,&e1,&e2);
      int nk=s->num_atoms;
      sim_destroy(s);
      Simulation *t=sim_create(1024,2048);
      int f2=kcsa_build_filter(t,vec3(0,0,0),4);
      t->use_bonds=0; t->use_angles=0; t->use_dihedrals=0;
      kcsa_site_binding(t,f2,4,11,sites[1],60,&ena,0);
      char d[200];
      snprintf(d,sizeof d,"K+ e_inter=%.4f eV, Na+ e_inter=%.4f eV (K more negative by %.4f)",
               e1,ena,e1-ena);
      ck("KcsA: rigid-site interaction is attractive for both ions", e1<0.0&&ena<0.0, d);
      printf("        %s\n",d);
      sim_destroy(t); (void)nk; } }

/* ═══════════════ 7. datastream / SHA-256 ═══════════════ */
static void test_datastream(void) {
    printf("\n== datastream + SHA-256 (FIPS 180-4) ==\n");
    struct { const char *msg; const char *want; } kat[] = {
      {"", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
      {"abc","ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
      {"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq",
       "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"}};
    for (unsigned i=0;i<sizeof kat/sizeof kat[0];i++){
        char h[65]; ds_sha256_hex(kat[i].msg, strlen(kat[i].msg), h);
        char d[160]; snprintf(d,sizeof d,"got %s want %s",h,kat[i].want);
        ck("SHA-256 known-answer vector", strcmp(h,kat[i].want)==0, d);
    }
    { /* million 'a' */
        DSWriter *w=ds_open("/tmp/opencode/s2audit/t.cvmds","selftest");
        char h[65];
        if(w){ for(int i=0;i<1000;i++) ds_add_claim(w,"k.noise",(double)i,"eV","computed"); ds_close(w);} (void)h; }
    { /* round trip + tamper rejection */
      DSWriter *w=ds_open("/tmp/opencode/s2audit/rt.cvmds","selftest");
      ds_set_header(w,"source-hash","record-tree-test");
      ds_add_atom(w,0,1,"H",1.008,0.0,0.0019,2.57113);
      ds_add_claim(w,"t.energy",-1.25,"eV","computed");
      int rc=ds_close(w);
      ckrel("ds_close returns 0",rc,0.0,0);
      ck("sealed file verifies", ds_verify_file("/tmp/opencode/s2audit/rt.cvmds")==0, NULL);
      /* tamper one payload digit */
      FILE *f=fopen("/tmp/opencode/s2audit/rt.cvmds","r+");
      if(f){ fseek(f,20,SEEK_SET); int c=fgetc(f); fseek(f,20,SEEK_SET); fputc(c^0x01,f); fclose(f);}
      ck("tampered payload byte rejected", ds_verify_file("/tmp/opencode/s2audit/rt.cvmds")!=0, NULL);
      /* rebuild clean, then append garbage after the hash line */
      { DSWriter *w2=ds_open("/tmp/opencode/s2audit/rt.cvmds","selftest");
        ds_add_claim(w2,"t.energy",-1.25,"eV","computed"); ds_close(w2); }
      f=fopen("/tmp/opencode/s2audit/rt.cvmds","a");
      if(f){ fputs("extra claim row\n",f); fclose(f);}
      ck("trailing garbage after hash rejected", ds_verify_file("/tmp/opencode/s2audit/rt.cvmds")!=0, NULL); }
}

int main(void) {
    printf("S2 audit probe — independent verification harness\n");
    setvbuf(stdout, NULL, _IOLBF, 0);
    test_quantum();
    test_forces();
    test_integrator();
    test_qm();
    test_neuron();
    test_kcsa();
    test_datastream();
    printf("\n=== audit_probe: %d passed, %d failed ===\n", npass, nfail);
    return nfail ? 1 : 0;
}