#include <stdio.h>
#include "../../include/sim.h"
#include "../../include/aminoacids.h"
static void rep(const char*tag, Simulation*s, int first){
  double sum=0; int n=0;
  for(int i=first;i<s->num_atoms;i++){ sum+=s->atoms[i].partial_charge; n++; }
  printf("  %-34s %2d atoms  sum(q) = %+.6f e\n",tag,n,sum);
}
int main(void){
  { Simulation *s=sim_create(64,64); sim_place_glycine(s,vec3(0,0,0));       rep("glycine (free)",s,0); sim_destroy(s); }
  { Simulation *s=sim_create(64,64); sim_place_alanine(s,vec3(0,0,0));      rep("alanine (free)",s,0); sim_destroy(s); }
  { Simulation *s=sim_create(128,128); int a; int n=sim_place_dipeptide_GlyAla(s,vec3(0,0,0),&a); rep("dipeptide Gly-Ala",s,0); sim_destroy(s); }
  { Simulation *s=sim_create(512,512); sim_place_polyalanine(s,vec3(0,0,0),4,0); rep("polyalanine x4",s,0); sim_destroy(s); }
  return 0;
}
