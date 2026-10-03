#include <stdio.h>
#include <string.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
static double V(Atom*at){ return forces_nonbonded_energy(at,0,1,NULL,1,1,1.0).lj_energy
                              + forces_nonbonded_energy(at,0,1,NULL,1,1,1.0).coulomb_energy; }
int main(void){
  double sig=3.0, eps=0.01, q=0.8, h=1e-6;
  double worst=0;
  printf("  Convention-free: engine force vs -d/dx of the ENGINE'S OWN energy\n");
  printf("     r(A)     engine Fx        central-difference -dV/dx     rel\n");
  for(double r=1.2; r<=11.0; r+=0.97){
    double F[3];
    for(int k=0;k<3;k++){
      Vec3 ax = (k==0)?vec3(h,0,0):(k==1)?vec3(0,h,0):vec3(0,0,h);
      Vec3 bx = (k==0)?vec3(-h,0,0):(k==1)?vec3(0,-h,0):vec3(0,0,-h);
      Atom p[2],m[2];
      for(int i=0;i<2;i++){ memset(&p[i],0,sizeof(Atom)); memset(&m[i],0,sizeof(Atom)); }
      p[0].element=8; p[1].element=19;
      p[0].position=vec3_add(vec3(0,0,0),ax); p[1].position=vec3(r,0,0);
      m[0].element=8; m[1].element=19;
      m[0].position=vec3_sub(vec3(0,0,0),ax); m[1].position=vec3(r,0,0);
      for(int i=0;i<2;i++){ p[i].lj_sigma=sig; p[i].lj_epsilon=eps; p[i].partial_charge=(i? -q:q);
                             m[i].lj_sigma=sig; m[i].lj_epsilon=eps; m[i].partial_charge=(i? -q:q); }
      F[k] = -(V(p)-V(m))/(2*h);
    }
    Atom at[2]; memset(at,0,sizeof at);
    at[0].element=8; at[1].element=19;
    at[0].position=vec3(0,0,0); at[1].position=vec3(r,0,0);
    at[0].lj_sigma=at[1].lj_sigma=sig;
    at[0].lj_epsilon=at[1].lj_epsilon=eps;
    at[0].partial_charge=q; at[1].partial_charge=-q;
    forces_nonbonded_pair(at,0,1,NULL,1,1,1.0);
    double err=0,ref=0;
    double fc[3]={at[0].force.x,at[0].force.y,at[0].force.z};
    for(int k=0;k<3;k++){ err+=pow(fc[k]-F[k],2); ref+=F[k]*F[k]; }
    double rel = sqrt(err)/(sqrt(ref)>1e-30?sqrt(ref):1.0);
    if(rel>worst) worst=rel;
    printf("  %5.2f  %13.6f  %19.6f   %.2e\n",r,at[0].force.x,F[0],rel);
  }
  printf("\n  worst relative deviation force vs -dV/dx: %.3e\n",worst);
  return worst>1e-6;
}
