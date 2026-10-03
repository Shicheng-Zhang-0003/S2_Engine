#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/kcsa_filter.h"
int main(void){
  for (int k=0;k<2;k++){
    double rion = k? 1.02 : 1.38, Z = k? 11:19;
    Simulation *s=sim_create(8,8);
    int ion = sim_add_ion(s,(int)Z,1,vec3(0,0,0),1.0);
    kcsa_set_ion_radius(s,ion,(int)Z);
    double ri = s->atoms[ion].lj_sigma, ei = s->atoms[ion].lj_epsilon;
    /* one O at the contact distance */
    sim_add_atom(s,8,vec3(2*rion+1.40,0,0),0);
    double ro = s->atoms[ion+1].lj_sigma, eo = s->atoms[ion+1].lj_epsilon;
    double sig = 0.5*(ri+ro), eps = sqrt(ei*eo);
    /* locate the ion-O minimum for an LJ pair (o=0), then the contact */
    double rm = pow(2.0*pow(sig,12.0)/pow(sig,6.0),1.0/6.0);
    printf("  %s  ion_sigma=%.4f A  O_sigma=%.4f A  LB_sigma=%.4f A  eps=%.4f\n",
           k?"Na+":"K+ ", ri, ro, sig, eps);
    printf("      LJ minimum r_m = 2^(1/6) sigma = %.4f A   (tabulated contact %.2f A)\n", 2.0*pow(sig,2.0)/pow(sig,2.0)*pow(2.0,1.0/6.0)*sig, 2*rion+1.40);
    printf("      check: 2^(1/6)*%.4f = %.4f\n", sig, pow(2.0,1.0/6.0)*sig);
    sim_destroy(s);
  }
  return 0;
}
