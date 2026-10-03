#include <stdio.h>
#include <string.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
int main(void){
  double worstE=0, worstF=0;
  printf("  r(A)     engine E(eV)     closed form            rel diff\n");
  double sig=3.0, eps=0.01, q=0.8;
  for(double r=1.0; r<=12.0; r+=0.37){
    Atom at[2]; memset(at,0,sizeof at); Atom *a=&at[0], *b=&at[1];
    a->element=8; b->element=19;
    a->position=vec3(0,0,0); b->position=vec3(r,0,0);
    a->lj_sigma=sig; a->lj_epsilon=eps; a->partial_charge=q;
    b->lj_sigma=sig; b->lj_epsilon=eps; b->partial_charge=-q;
    PairEnergy pe=forces_nonbonded_energy(at,0,1,NULL,1,1,1.0);
    double sr=sig/r;
    double ref=4*eps*(pow(sr,12)-pow(sr,6)) - COULOMB_MD*q*q/r;
    double rel=fabs((pe.lj_energy+pe.coulomb_energy)-ref)/fabs(ref);
    if(rel>worstE) worstE=rel;
    printf("  %5.2f  %14.8f  %14.8f   %.2e\n",r,pe.lj_energy+pe.coulomb_energy,ref,rel);
  }
  printf("\n  worst relative deviation of the pair energy from the closed form: %.3e\n",worstE);
  /* now the FORCE, against the analytic derivative */
  printf("\n  force check against -dV/dr with the engine's documented convention\n");
  printf("  (F_a = (dV/dr) * rhat_ab, i.e. the potential gradient, not its negation)\n");
  worstF=0;
  for(double r=1.2; r<=11.0; r+=0.53){
    Atom at[2]; memset(at,0,sizeof at); Atom *a=&at[0], *b=&at[1];
    a->element=8; b->element=19;
    a->position=vec3(0,0,0); b->position=vec3(r,0,0);
    a->lj_sigma=sig; a->lj_epsilon=eps; a->partial_charge=q;
    b->lj_sigma=sig; b->lj_epsilon=eps; b->partial_charge=-q;
    forces_nonbonded_pair(at,0,1,NULL,1,1,1.0);
    double sr=sig/r;
    double dVdr = 4*eps*(-12*pow(sr,12)/r + 6*pow(sr,6)/r) + COULOMB_MD*q*q/(r*r);
    /* d_ij = a - b = -r xhat, so F_a = (dV/dr) * (d_ij/r) -> x comp = -dVdr */
    double rel=fabs(a->force.x + dVdr)/fabs(dVdr);
    if(rel>worstF) worstF=rel;
    if(r<4.0) printf("  r=%5.2f  engine Fx=%12.6f   -dV/dr=%12.6f   rel %.2e\n",r,a->force.x,-dVdr,rel);
  }
  printf("\n  worst relative deviation of the force from -dV/dr: %.3e\n",worstF);
  return (worstE>1e-12||worstF>1e-12);
}
