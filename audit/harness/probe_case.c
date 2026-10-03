#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/qm.h"
static void t(const char*tag,Simulation*s,int n,int ion){
  double q[128];
  int rc = qm_qeq(s, ion>=0?1.0:0.0, 1.0, q);
  double sum=0,mx=0;
  for(int i=0;i<n;i++){ if(i==ion)continue; sum+=q[i]; if(fabs(q[i])>mx)mx=fabs(q[i]); }
  printf("  %-34s rc=%d shell_sum=%+.4e max|q|=%.4f\n",tag,rc,sum,mx);
  sim_destroy(s);
}
int main(void){
  { Simulation *s=sim_create(16,32);
    for(int i=0;i<6;i++) sim_add_atom(s,8,vec3(0.95*i,0,0),0);
    t("6 unbonded O, 0.95 A",s,6,-1); }
  { Simulation *s=sim_create(16,32);
    for(int i=0;i<6;i++) sim_add_atom(s,8,vec3(0.80*i,0,0),0);
    t("6 unbonded O, 0.80 A",s,6,-1); }
  { Simulation *s=sim_create(16,32);
    for(int i=0;i<8;i++) sim_add_atom(s,1,vec3(0.55*i,0.1,0),0);
    t("8 unbonded H, 0.55 A",s,8,-1); }
  { Simulation *s=sim_create(16,32);
    for(int i=0;i<6;i++) sim_add_atom(s,8,vec3(0.95*i,0.1*i,0.05*i),0);
    int ion=sim_add_ion(s,19,1,vec3(0.4,0.4,0.4),1.0);
    t("SCF: 6 O + K+ (tight)",s,7,ion); }
  return 0;
}
