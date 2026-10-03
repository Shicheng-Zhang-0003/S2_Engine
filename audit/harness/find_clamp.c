#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/qm.h"
/* Find a DETERMINISTIC case where the unbounded QEq solution exceeds 2 e. */
int main(void){
  int bestN=0; double bestGap=0; static double bestP[32][3]; static int bestZ[32];
  for(int N=2;N<=8;N++){
    for(double sp=0.55; sp<2.0; sp+=0.05){
      Simulation *s=sim_create(16,32);
      /* alternating electronegativity, tight linear chain, unbonded */
      for(int i=0;i<N;i++) sim_add_atom(s, (i%2)?6:8, vec3(sp*i,0.0,0.0), 0);
      double q[128];
      if(qm_qeq(s,0.0,1.0,q)==0){
        double mx=0; for(int i=0;i<N;i++) if(fabs(q[i])>mx) mx=fabs(q[i]);
        if(mx>bestGap){ bestGap=mx; bestN=N;
          for(int i=0;i<N;i++){bestP[i][0]=sp*i;bestP[i][1]=0;bestP[i][2]=0;bestZ[i]=(i%2)?6:8;} }
      }
      sim_destroy(s);
    }
  }
  printf("  worst |q| found: %.4f at N=%d\n",bestGap,bestN);
  if(bestN>0){
    Simulation *s=sim_create(16,32);
    for(int i=0;i<bestN;i++) sim_add_atom(s,bestZ[i],vec3(bestP[i][0],0,0),0);
    double q[128];
    int rc=qm_qeq(s,0.0,1.0,q);
    double sum=0,mx=0;
    for(int i=0;i<bestN;i++){sum+=q[i]; if(fabs(q[i])>mx)mx=fabs(q[i]);}
    printf("  case: N=%d alternating C/O chain at %.2f A spacing, unbonded\n",bestN,bestP[1][0]);
    printf("  rc=%d  sum(q)=%+.4e  max|q|=%.4f  %s\n",rc,sum,mx,
           (fabs(sum)<1e-12&&mx<=2.0+1e-12)?"OK":"VIOLATES a stated invariant");
    for(int i=0;i<bestN;i++) printf("    q[%d] Z=%d = %+.6f\n",i,bestZ[i],q[i]);
    sim_destroy(s);
  }
  return 0;
}
