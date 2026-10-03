/* A12b — characterise the zero-force minimiser path precisely. */
#include <stdio.h>
#include <math.h>
#include "../../include/sim.h"
#include "../../include/integrator.h"
#include "../../include/forces.h"
int main(void){
  Simulation *s=sim_create(8,8);
  sim_add_atom(s,6,vec3(0,0,0),0);
  s->use_lj=0;s->use_coulomb=0;s->use_bonds=0;s->use_angles=0;s->use_dihedrals=0;
  forces_calculate(s);
  printf("  single atom, all terms off: |F| = %.17g  (exactly zero: %s)\n",
     vec3_norm(s->atoms[0].force), vec3_norm(s->atoms[0].force)==0.0?"yes":"no");
  double E=integrator_minimize(s,10,0.01,0.0);
  printf("  after minimize: E=%g pos.x=%s\n", E, isfinite(s->atoms[0].position.x)?"finite":"NaN");
  printf("\n  The NaN rollback guard (isnan(E_new)) is what rescues this, not\n");
  printf("  the tolerance test: `max_force < force_tolerance` is 0 < 0 = FALSE,\n");
  printf("  so the loop body RUNS, scale = 0.01/0 = inf, and every displacement\n");
  printf("  becomes 0*inf = NaN. The next forces_calculate returns NaN, the\n");
  printf("  isnan branch rolls the positions back and returns early.\n");
  printf("  Observable result is correct but the function exits after ONE wasted\n");
  printf("  iteration instead of recognizing convergence.\n");
  return 0;
}
