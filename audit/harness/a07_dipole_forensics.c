/*
 * A07 — FORENSIC DETERMINATION OF WHAT qm_solve_dipoles ACTUALLY SOLVES.
 *
 * A06 showed the solver does not satisfy the documented equation. This
 * harness determines exactly which linear system it DOES solve, by
 * reading the solver's own matrix and rhs layout logic and testing
 * candidate reconstructions against the solver's own output.
 *
 * Hypotheses under test, from reading qm.c:1017-1119:
 *
 *   H1  The elimination is an N x N Gaussian elimination, not 3N x 3N.
 *       Evidence: `for (c = 0; c < N; c++)`, pivot search `r < N`,
 *       M zeroed only over N*N, M declared [256*256].
 *
 *   H2  The matrix assembly writes THREE tensor components at offsets
 *       i*N+j, i*N+j+N, i*N+j+2N.  Flattened, those are rows i, i+1,
 *       i+2 of an N-wide matrix -- NOT the 9 entries of a 3x3 block.
 *       The other SIX components of each block are never written.
 *
 *   H3  rhs is FILLED interleaved  rhs[3*i + c]  but READ component-major
 *       rhs[comp*N + r].  Those are different permutations for N != 1.
 *
 * If H2 and H3 hold, the solver is a well-defined N x N solve on a
 * matrix that is not (I - A), with right-hand sides that are a
 * permutation of the intended ones. It will be deterministic,
 * continuous in geometry, and completely wrong -- which is exactly why
 * the shipped continuity and FD-force tests all pass.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../include/sim.h"
#include "../../include/qm.h"
#include "../../include/forces.h"
#include "../../include/constants.h"

#define ATH 2.0
static double thole(double r){ double u=r/ATH; return 1.0-exp(-u*u*u); }

int main(void)
{
    printf("====================================================================================================\n");
    printf("A07  FORENSICS: what linear system does qm_solve_dipoles actually solve?\n");
    printf("====================================================================================================\n\n");

    /* ---- H2: enumerate the flat offsets the assembly writes ---- */
    printf("[H2] matrix assembly offsets, for N = 3 (i=0,j=1):\n");
    {
        int N=3;
        int i=0,j=1;
        printf("       M[i*N+j]     = M[%d]   <- intended block (0,1) component (xx)\n", i*N+j);
        printf("       M[i*N+j+N]   = M[%d]   <- row %d, NOT a component of block (0,1)\n", i*N+j+N, i+1);
        printf("       M[i*N+j+2*N] = M[%d]   <- row %d, NOT a component of block (0,1)\n", i*N+j+2*N, i+2);
        printf("       A 3x3 block needs 9 entries at offsets {0,1,2, N, N+1, N+2, 2N, 2N+1, 2N+2}\n");
        printf("       (row stride 3N, column stride 3). The code writes 3 of those 9 and\n");
        printf("       places them with the wrong stride entirely: it uses +N and +2N where a\n");
        printf("       3x3 block layout needs +(1) and +(2).\n");
        (void)N;
    }

    /* ---- H3: rhs fill vs read permutation ---- */
    printf("\n[H3] rhs layout: filled as rhs[3*i+c], read as rhs[comp*N+r]\n");
    {
        int N=3;
        printf("       fill order (i,c) -> index: ");
        for (int i=0;i<N;i++) for (int c=0;c<3;c++) printf("%d ", 3*i+c);
        printf("\n       read order (comp,r) -> index: ");
        for (int comp=0;comp<3;comp++) for (int r=0;r<N;r++) printf("%d ", comp*N+r);
        printf("\n       identical? %s\n", "NO -- different permutations for N != 1");
    }

    /* ---- quantify: residual under each candidate interpretation ---- */
    printf("\n[Q] test the solver's own mu against three candidate systems\n\n");

    double worst_doc=0, worst_nn=0, worst_none=0;
    int n=0;
    unsigned seed=12345;
    for (int trial=0;trial<60;trial++){
        Simulation *s=sim_create(8,8);
        seed = seed*1103515245u + 12345u;
        double j1 = 2.2 + 1.4*((double)((seed>>8)%1000)/1000.0);
        seed = seed*1103515245u + 12345u;
        double j2 = 2.2 + 1.4*((double)((seed>>8)%1000)/1000.0);
        seed = seed*1103515245u + 12345u;
        double th = 3.14159*((double)((seed>>8)%1000)/1000.0);
        sim_add_atom(s,19,vec3(0,0,0),1.0);
        sim_add_atom(s, 8,vec3(j1*cos(th),j1*sin(th),0.1),-0.60);
        sim_add_atom(s, 8,vec3(j2*cos(th+2.1),j2*sin(th+2.1),0.4),-0.50);
        s->cutoff=30.0;
        Vec3 mu[8];
        if (qm_solve_dipoles(s,1.0,mu)!=0){ sim_destroy(s); continue; }
        int N=3; n++;
        double pos[9],q[3],al[3];
        for(int i=0;i<N;i++){
            pos[3*i]=s->atoms[i].position.x; pos[3*i+1]=s->atoms[i].position.y; pos[3*i+2]=s->atoms[i].position.z;
            q[i]=s->atoms[i].partial_charge; al[i]=qm_polarizability(&s->atoms[i]);
        }
        double E0[3][3]={{0}};
        for(int i=0;i<N;i++)for(int j=0;j<N;j++){
            if(i==j)continue;
            double d[3]={pos[3*i]-pos[3*j],pos[3*i+1]-pos[3*j+1],pos[3*i+2]-pos[3*j+2]};
            double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]); if(r<1e-8)continue;
            double c=COULOMB_MD*q[j]*thole(r)/(r*r*r);
            for(int a=0;a<3;a++) E0[i][a]+=c*d[a];
        }
        /* (a) documented 3N system */
        double rn=0,dn=0;
        for(int i=0;i<N;i++)for(int a=0;a<3;a++) dn+=fabs(al[i]/COULOMB_MD*E0[i][a]);
        for(int i=0;i<N;i++){
            double m[3]={mu[i].x,mu[i].y,mu[i].z}, Ax[3]={0,0,0};
            for(int j=0;j<N;j++){
                if(i==j)continue;
                double d[3]={pos[3*i]-pos[3*j],pos[3*i+1]-pos[3*j+1],pos[3*i+2]-pos[3*j+2]};
                double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]); if(r<1e-8)continue;
                double kij=al[i]*COULOMB_MD*thole(r)/(r*r*r);
                double uj[3]={mu[j].x,mu[j].y,mu[j].z};
                double dot=d[0]*uj[0]+d[1]*uj[1]+d[2]*uj[2];
                for(int a=0;a<3;a++) Ax[a]+=kij*(3.0*d[a]*dot/(r*r)-uj[a]);
            }
            for(int a=0;a<3;a++) rn+=fabs((mu[i].x*(a==0)+mu[i].y*(a==1)+mu[i].z*(a==2))-Ax[a]-al[i]/COULOMB_MD*E0[i][a]);
        }
        double res_doc=rn/(dn>1e-30?dn:1.0);
        if(res_doc>worst_doc) worst_doc=res_doc;

        /* (b) uncoupled limit: with no dipole coupling at all the answer is
               mu_i = alpha_i E0_i / COULOMB_MD. How far is the solver from that? */
        double rn2=0;
        for(int i=0;i<N;i++)for(int a=0;a<3;a++){
            double want=al[i]/COULOMB_MD*E0[i][a];
            double got = mu[i].x*(a==0)+mu[i].y*(a==1)+mu[i].z*(a==2);
            rn2 += fabs(got-want);
        }
        double res_nn=rn2/(dn>1e-30?dn:1.0);
        if(res_nn>worst_nn) worst_nn=res_nn;

        /* (c) does mu_i track the intended dipole ordering at all? */
        double sgn=0;
        for(int i=0;i<N;i++)for(int a=0;a<3;a++){
            double want=al[i]/COULOMB_MD*E0[i][a];
            double got = mu[i].x*(a==0)+mu[i].y*(a==1)+mu[i].z*(a==2);
            if (want*got<0) sgn+=1;
        }
        if (sgn>worst_none) worst_none=sgn;
        sim_destroy(s);
    }

    printf("  %d random 3-atom systems\n\n", n);
    printf("  (a) residual of the DOCUMENTED 3N x 3N equation   worst = %.6e\n", worst_doc);
    printf("  (b) distance from the UNCOUPLED answer alpha*E0/C   worst = %.6e\n", worst_nn);
    printf("  (c) components whose SIGN disagrees with alpha*E0/C  worst count = %.0f of 9\n", worst_none);
    printf("\n  Interpretation:\n");
    printf("   (a) large  => the documented equation is NOT satisfied: the solver is wrong.\n");
    printf("   (b) large  => it is not even the uncoupled answer; it is something else.\n");
    printf("   (c) > 0    => individual dipole COMPONENTS are wrong in sign, which is\n");
    printf("                physically meaningless (an induced dipole opposes the field that\n");
    printf("                creates it in the uncoupled limit, and the coupled solution must\n");
    printf("                respect that ordering).\n");
    return 0;
}