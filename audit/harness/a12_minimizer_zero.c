/* A12 — minimiser with force_tolerance == 0 and zero net force (finding M3). */
#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/integrator.h"
#include "../../include/forces.h"
static Simulation *converged(void){
  Simulation *s=sim_create(8,8);
  sim_add_atom(s,6,vec3(0,0,0),0);
  s->use_lj=0; s->use_coulomb=0;
  s->use_bonds=0; s->use_angles=0; s->use_dihedrals=0;
  return s;
}
int main(void){
  int bad=0;
  { Simulation *s=converged();
    double E=integrator_minimize(s,10,0.01,0.0);
    int ok=isfinite(E)&&isfinite(s->atoms[0].position.x);
    printf("%s  integrator_minimize   tol=0: E=%g pos.x=%g\n",ok?"PASS":"FAIL",E,s->atoms[0].position.x);
    if(!ok)bad++; sim_destroy(s); }
  { Simulation *s=converged();
    int frozen[1]={0};
    double E=integrator_minimize_frozen(s,frozen,10,0.01,0.0);
    int ok=isfinite(E)&&isfinite(s->atoms[0].position.x);
    printf("%s  minimize_frozen      tol=0: E=%g pos.x=%g\n",ok?"PASS":"FAIL",E,s->atoms[0].position.x);
    if(!ok)bad++; sim_destroy(s); }
  { Simulation *s=converged();
    double E=integrator_fire(s,10,0.5,0.0);
    int ok=isfinite(E)&&isfinite(s->atoms[0].position.x);
    printf("%s  integrator_fire      tol=0: E=%g pos.x=%g\n",ok?"PASS":"FAIL",E,s->atoms[0].position.x);
    if(!ok)bad++; sim_destroy(s); }
  return bad?1:0;
}
