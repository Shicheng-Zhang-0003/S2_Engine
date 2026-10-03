/*
 * A06 — THE COUPLED-DIPOLE SOLVER, called directly.
 *
 * qm.c:1002 documents the system it claims to solve as
 *
 *     (I - A) mu = alpha (*) E0,   A_ij = alpha_i k f(r) [3 d d - I] / r^3 / eps
 *
 * a 3N x 3N linear system. This harness:
 *   1. builds a small Simulation,
 *   2. calls qm_solve_dipoles(),
 *   3. computes the residual of the DOCUMENTED equation,
 *   4. compares against a correctly-assembled reference solve.
 *
 * If the residual is large the solver is not solving the equation it
 * documents, whatever it is solving.
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

/* Correct 3N x 3N assembly and solve, numpy-free. */
static int ref_solve(int N, const double *pos, const double *q,
                     const double *al, double *mu_out)
{
    static double E0[64][3];
    for (int i=0;i<N;i++){ E0[i][0]=E0[i][1]=E0[i][2]=0.0; }
    for (int i=0;i<N;i++) for (int j=0;j<N;j++){
        if (i==j) continue;
        double dx=pos[3*i]-pos[3*3*j+0];
        double dy=pos[3*i+1]-pos[3*j+1];
        double dz=pos[3*i+2]-pos[3*j+2];
        (void)dx;(void)dy;(void)dz;
    }
    for (int i=0;i<N;i++) for (int j=0;j<N;j++){
        if (i==j) continue;
        double d[3]={pos[3*i]-pos[3*j], pos[3*i+1]-pos[3*j+1], pos[3*i+2]-pos[3*j+2]};
        double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if (r<1e-8) continue;
        double f=thole(r), c=COULOMB_MD*q[j]*f/(r*r*r);
        for (int a=0;a<3;a++) E0[i][a]+=c*d[a];
    }
    int M3=3*N;
    static double M[64*64], rhs[64], sol[64];
    for (int i=0;i<M3*M3;i++) M[i]=0.0;
    for (int i=0;i<M3;i++) M[i*M3+i]=1.0;
    for (int i=0;i<N;i++) for (int j=0;j<N;j++){
        if (i==j) continue;
        double d[3]={pos[3*i]-pos[3*j], pos[3*i+1]-pos[3*j+1], pos[3*i+2]-pos[3*j+2]};
        double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if (r<1e-8) continue;
        double f=thole(r);
        if (f==0.0) continue;
        double r3=r*r*r;
        double u[3]={d[0]/r,d[1]/r,d[2]/r};
        double kij=al[i]*COULOMB_MD*f/r3;
        for (int a=0;a<3;a++) for (int b=0;b<3;b++)
            M[(3*i+a)*M3+(3*j+b)] -= kij*(3.0*u[a]*u[b] - (a==b?1.0:0.0));
    }
    for (int i=0;i<N;i++) for (int a=0;a<3;a++)
        rhs[3*i+a] = al[i]/COULOMB_MD*E0[i][a];
    for (int c=0;c<M3;c++){
        int piv=c; double best=fabs(M[c*M3+c]);
        for (int r=c+1;r<M3;r++) if (fabs(M[r*M3+c])>best){best=fabs(M[r*M3+c]);piv=r;}
        if (!(best>1e-12)) return -1;
        if (piv!=c){ for(int k=0;k<M3;k++){double t=M[c*M3+k];M[c*M3+k]=M[piv*M3+k];M[piv*M3+k]=t;}
                     double t=rhs[c];rhs[c]=rhs[piv];rhs[piv]=t; }
        double dg=M[c*M3+c];
        for (int r=c+1;r<M3;r++){
            double fc=M[r*M3+c]/dg; if (fc==0.0) continue;
            for(int k=0;k<M3;k++) M[r*M3+k]-=fc*M[c*M3+k];
            rhs[r]-=fc*rhs[c];
        }
    }
    for (int r=M3-1;r>=0;r--){
        double acc=rhs[r];
        for(int k=r+1;k<M3;k++) acc-=M[r*M3+k]*sol[k];
        sol[r]=acc/M[r*M3+r];
    }
    for (int i=0;i<N;i++) for(int a=0;a<3;a++) mu_out[3*i+a]=sol[3*i+a];
    return 0;
}

int main(void)
{
    printf("====================================================================================================\n");
    printf("A06  qm_solve_dipoles — residual of the DOCUMENTED equation  (I - A) mu = alpha (*) E0\n");
    printf("====================================================================================================\n\n");

    struct { const char *name; int Z; double q; } spec[] = {
        {"K+ (19)",  19,  1.0},
        {"O  (8)",    8, -0.60},
        {"H  (1)",    1,  0.20},
    };
    double worst_res=0.0, worst_rel=0.0;
    int ncase=0, nsing=0;

    for (int which=0; which<3; which++)
    for (int trial=0; trial<40; trial++)
    {
        Simulation *s = sim_create(8, 8);
        if (!s) continue;
        /* 3 atoms, deliberately NOT bonded so qm_pair_excluded is false */
        int i0 = sim_add_atom(s, 19, vec3(0,0,0), 1.0);
        int i1 = sim_add_atom(s,  8, vec3(2.6,0.2,0.1), -0.60);
        int i2 = sim_add_atom(s,  1, vec3(1.2,2.4,-0.3), 0.20);
        s->cutoff = 30.0;

        Vec3 mu[8];
        int rc = qm_solve_dipoles(s, 1.0, mu);
        ncase++;
        if (rc != 0) { nsing++; sim_destroy(s); continue; }

        int N = 3;
        double pos[9], q[3], al[3];
        for (int i=0;i<N;i++){
            pos[3*i+0]=s->atoms[i].position.x;
            pos[3*i+1]=s->atoms[i].position.y;
            pos[3*i+2]=s->atoms[i].position.z;
            q[i]=s->atoms[i].partial_charge;
            al[i]=qm_polarizability(&s->atoms[i]);
        }
        /* residual of the documented equation, built independently */
        double E0[3][3];
        for (int i=0;i<N;i++) E0[i][0]=E0[i][1]=E0[i][2]=0.0;
        for (int i=0;i<N;i++) for (int j=0;j<N;j++){
            if (i==j) continue;
            double d[3]={pos[3*i]-pos[3*j],pos[3*i+1]-pos[3*j+1],pos[3*i+2]-pos[3*j+2]};
            double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
            if (r<1e-8) continue;
            double c=COULOMB_MD*q[j]*thole(r)/(r*r*r);
            for (int a=0;a<3;a++) E0[i][a]+=c*d[a];
        }
        double rn=0.0, dn=0.0;
        for (int i=0;i<N;i++) for (int a=0;a<3;a++) dn += fabs(al[i]/COULOMB_MD*E0[i][a]);
        for (int i=0;i<N;i++){
            double mui[3]={mu[i].x,mu[i].y,mu[i].z};
            double Ax[3]={0,0,0};
            for (int j=0;j<N;j++){
                if (i==j) continue;
                double d[3]={pos[3*i]-pos[3*j],pos[3*i+1]-pos[3*j+1],pos[3*i+2]-pos[3*j+2]};
                double r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
                if (r<1e-8) continue;
                double f=thole(r); if (f==0.0) continue;
                double kij=al[i]*COULOMB_MD*f/(r*r*r);
                double uj[3]={mu[j].x,mu[j].y,mu[j].z};
                /* dot does NOT depend on the component index b; looping b
                 * around the accumulation applied the same term three
                 * times. The harness was wrong, not the engine. */
                double dot = d[0]*uj[0]+d[1]*uj[1]+d[2]*uj[2];
                for (int a=0;a<3;a++)
                    Ax[a] += kij*(3.0*d[a]*dot/(r*r) - uj[a]);
            }
            for (int a=0;a<3;a++){
                double want = al[i]/COULOMB_MD*E0[i][a];
                double got  = mui[a] - Ax[a];
                rn += fabs(got-want);
            }
        }
        double res = rn/(dn>1e-30?dn:1.0);
        if (res>worst_res) worst_res=res;
        double muref[9];
        if (ref_solve(N,pos,q,al,muref)==0){
            double num=0, den=0;
            for (int i=0;i<3*N;i++){
                double a1 = (&mu[0].x)[i], a2 = muref[i];
                num += (a1-a2)*(a1-a2); den += a2*a2;
            }
            double rel = sqrt(num/(den>1e-30?den:1.0));
            if (rel>worst_rel) worst_rel=rel;
        }
        sim_destroy(s);
    }

    printf("  %d systems (%d reported singular)\n", ncase, nsing);
    printf("\n  WORST residual ||(I-A)mu - alpha(*)E0|| / ||alpha(*)E0||   = %.6e\n", worst_res);
    printf("  WORST relative difference from an INDEPENDENT correct solve    = %.6e\n", worst_rel);
    printf("\n  If the solver solves the documented 3N x 3N system, the first number must be\n");
    printf("  at solver precision (~1e-13) and the second ~1e-13 too.\n");

    if (worst_res < 1e-10 && worst_rel < 1e-10) {
        printf("\nPASS  qm_solve_dipoles solves the documented equation.\n");
    } else {
        printf("\nFAIL  qm_solve_dipoles does NOT solve the documented equation.\n");
    }
    return (worst_res < 1e-10 && worst_rel < 1e-10) ? 0 : 1;
}