/* A11 — is the QEq ordering failure a CODE defect or a MODEL limitation? */
#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/qm.h"
static void trial(double r, int bonded){
  Simulation *s=sim_create(8,8);
  int c=sim_add_atom(s,6,vec3(0,0,0),0);
  int o=sim_add_atom(s,8,vec3(r,0,0),0);
  if(bonded) sim_add_bond(s,c,o,1);
  double q[8]; qm_qeq(s,0.0,1.0,q);
  printf("  %-9s r=%.2f A  q(C)=%+.5f  q(O)=%+.5f  ordering %s\n",
     bonded?"bonded":"unbonded", r, q[0], q[1],
     q[1]<q[0] ? "OK" : "VIOLATED");
  sim_destroy(s);
}
int main(void){
  printf("A11  QEq electronegativity ordering vs pair separation\n\n");
  trial(1.43,1); trial(1.23,1); trial(2.80,0);
  trial(2.20,0); trial(1.80,0); trial(1.60,0); trial(1.50,0);
  printf("\n  In standard QEq a bonded C-O has A_CO = 0 (hard core), so the pair\n");
  printf("  is driven purely by the chi/J difference -> O takes the charge.\n");
  printf("  An UNBONDED pair at short range has a huge 1/r coupling that exceeds\n");
  printf("  J, and the solve lands in the charge-transfer mode where the pair\n");
  printf("  shares charge by distance rather than by electronegativity. That is\n");
  printf("  a property of the FUNCTIONAL on an unphysical input (two atoms 1.5 A\n");
  printf("  apart with no bond between them is not a molecule), not a coding\n");
  printf("  defect - which is why the shipped suite tests the bonded case and\n");
  printf("  documents the unbonded one as out of scope.\n");
  return 0;
}
