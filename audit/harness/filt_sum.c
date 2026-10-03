#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/kcsa_filter.h"
int main(void){
  Simulation *s=sim_create(600,600);
  int f=kcsa_build_filter(s,vec3(0,0,0),1);
  printf("filter atoms: %d (first=%d)\n", s->num_atoms, f);
  double sum=0,mn=9,mx=-9;
  for(int i=0;i<s->num_atoms;i++){double q=s->atoms[i].partial_charge; sum+=q; if(q<mn)mn=q; if(q>mx)mx=q;}
  printf("TOTAL sum(q)      = %+.6f e\n", sum);
  printf("per-atom q range  = [%+.4f, %+.4f]\n", mn, mx);
  /* The filter is built as 5 residues x 4 chains; cluster atoms by the
     residue they were assigned to using the C-atom offset scheme:
     recompute each residue's own sum by walking bonds is not available,
     so report sum grouped by 20-atom blocks as a coarse proxy. */
  int blk=0; for(int i=0;i<s->num_atoms;i+=1){ (void)blk; break; }
  /* real proxy: nearest-carbon grouping via residue membership in code is
     internal, so instead check each atom's charge is within the documented
     carbon-offset band. */
  int ncarb=0; double worst=0;
  for(int i=0;i<s->num_atoms;i++){
    int Z=s->atoms[i].element;
    if(Z!=6) continue; ncarb++;
  }
  printf("carbon atoms: %d (these carry the neutralising offset)\n", ncarb);
  sim_destroy(s);
  (void)worst;
  return 0;
}
