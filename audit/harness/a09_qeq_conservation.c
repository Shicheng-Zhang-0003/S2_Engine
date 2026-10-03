/* A09 — QEq clamp vs charge conservation (audit finding C2). */
#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/qm.h"
int main(void){
  int breaks=0,tested=0; double worst=0, worstmaxq=0;
  unsigned seed=987654321u;
  for(int t=0;t<4000;t++){
    int N=2+((seed>>16)%6); seed=seed*1103515245u+12345u;
    Simulation *s=sim_create(16,32);
    for(int i=0;i<N;i++){
      seed=seed*1103515245u+12345u; int pk=(int)((seed>>8)%4);
      int Z = pk==0?6:pk==1?8:pk==2?7:1;
      seed=seed*1103515245u+12345u; double x=0.9+3.0*((double)((seed>>8)%1000)/1000.0);
      seed=seed*1103515245u+12345u; double y=3.0*((double)((seed>>8)%1000)/1000.0);
      seed=seed*1103515245u+12345u; double z=3.0*((double)((seed>>8)%1000)/1000.0);
      sim_add_atom(s,Z,vec3(x,y,z),0);
    }
    double q[128];
    if(qm_qeq(s,0.0,1.0,q)!=0){ sim_destroy(s); continue; }
    tested++;
    double sum=0,maxq=0;
    for(int i=0;i<N;i++){ sum+=q[i]; if(fabs(q[i])>maxq)maxq=fabs(q[i]); }
    double e=fabs(sum); if(e>worst) worst=e;
    if(maxq>worstmaxq) worstmaxq=maxq;
    if(e>1e-9){ breaks++; if(breaks<=4) printf("  N=%d sum=%+.3e max|q|=%.4f\n",N,sum,maxq); }
    sim_destroy(s);
  }
  printf("\n  %d/%d random clusters violate sum(q)=0 at 1e-9\n",breaks,tested);
  printf("  worst |sum(q)|          = %.3e\n",worst);
  printf("  worst max|q| over all   = %.6f  (bound is 2.0)\n",worstmaxq);
  int ok = (breaks==0) && (worstmaxq<=2.0+1e-12);
  printf("%s  charge conservation holds AND the 2 e bound is respected\n", ok?"PASS":"FAIL");
  return ok?0:1;
}
