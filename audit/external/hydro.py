import math,subprocess,os
# The engine tree is two levels up: audit/external -> audit -> v9R4.
import pathlib
V=str(pathlib.Path(__file__).resolve().parents[2])
open('/tmp/opencode/h.c','w').write(r'''
#include <stdio.h>
#include "quantum.h"
int main(void){
  int ns[]={1,2,3,4};
  for(int k=0;k<4;k++){int n=ns[k];
    for(int l=0;l<n;l++){
      double Z=3.7;
      printf("%d %d %.15e %.15e %.15e %.15e\n",n,l,
        quantum_expect_r(n,l,Z),quantum_expect_r2(n,l,Z),
        quantum_expect_invr(n,l,Z),quantum_expect_T(n,l,Z));
    }}
  return 0;}
''')
subprocess.run(['gcc','-O2','-std=c11','-I'+V+'/include','-o','/tmp/opencode/h','/tmp/opencode/h.c',
                V+'/build/quantum.o',V+'/build/periodic_table.o','-lm'],check=True)
out=subprocess.run(['/tmp/opencode/h'],capture_output=True,text=True).stdout
a0=0.529177210903
print("  n  l   <r>            <r^2>           <1/r>           <T>(eV)         vs Griffiths closed form")
worst=0
for line in out.strip().split('\n'):
    n,l,sr,sr2,sinvr,sT=[int(x) if i<2 else float(x) for i,x in enumerate(line.split())]
    Z=3.7
    er  = a0/(2*Z)*(3*n*n - l*(l+1))
    er2 = a0*a0*n*n/(2*Z*Z)*(5*n*n + 1 - 3*l*(l+1))
    einvr = Z/(a0*n*n)
    eT = (Z*Z/(2*n*n))*27.211386245988
    rel=[abs(sr/er-1),abs(sr2/er2-1),abs(sinvr/einvr-1),abs(sT/eT-1)]
    worst=max(worst,max(rel))
    flag="" if max(rel)<1e-13 else "  <-- MISMATCH"
    print("  %d  %d  %.10e  %.10e  %.10e  %.10e  max rel %.2e%s"%(n,l,sr,sr2,sinvr,sT,max(rel),flag))
print("\n  worst relative deviation from the closed forms: %.3e"%worst)
# Full-audit O11: failable (was print-only, exit 0 always).
import sys
sys.exit(0 if worst < 1e-12 else 1)
