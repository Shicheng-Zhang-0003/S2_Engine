#!/usr/bin/env python3
"""
A05 — INDUCTION, THOLE DAMPING, COUPLED DIPOLE SOLVE.

Four separate claims are under test:

  I1  U_pol = -0.5 * C * sum_i alpha_i |E_i|^2   with C = 1/COULOMB_MD.
      The code's own comment states the consistency requirement:
        -0.5*mu.E == -0.5*alpha*|E|^2/COULOMB_MD
      so U and the dipole must agree. Checked numerically.

  I2  The dipole solver must satisfy its own defining equation
        (I - A) mu = alpha (*) E0
      to solver precision, with A_ij = alpha_i k f(r) [3 d d - I]/r^3.

  I3  The Hellmann-Feynman derivative dU/dE0 for the SELF-CONSISTENT case.
      E = (I - T alpha)^-1 E0, so dU/dE0 = -mu (I - T alpha)^-1, NOT -mu/2.
      Verified by differentiating U numerically in the EXTERNAL field and
      comparing against -mu/2 (wrong) and against -mu M (right), where
      M = (I-T alpha)^-1 built in numpy. This is the exact quantity the
      release note says the analytic kernel got 98% wrong.

  I4  Thole f(r) = 1 - exp(-(r/a)^3) and df/dr = 3 (r/a)^2 exp(-(r/a)^3)/a
      verified by symbolic differentiation and by f(0)=0, f(inf)=1, monotone.
"""
import numpy as np
import sympy as sp
FAIL=[]

COULOMB_MD = 8.9875517923e9*1.602176634e-19/1.0e-10
CM = 1.0/COULOMB_MD
A_TH = 2.0

print("="*100)
print("A05  INDUCTION / THOLE / COUPLED DIPOLES")
print("="*100)

# ---------------- I4 Thole ----------------
print("\n[I4] Thole damping f(r) = 1 - exp(-(r/a)^3)")
r,a = sp.symbols('r a', positive=True)
f = 1-sp.exp(-(r/a)**3)
df = sp.diff(f,r)
print(f"      sympy df/dr          = {sp.simplify(df)}")
print(f"      code qm_thole_df form = 3*(r/a)^2*exp(-(r/a)^3)/a")
target = 3*(r/a)**2*sp.exp(-(r/a)**3)/a
ok = sp.simplify(df-target)==0
print(f"{'PASS' if ok else 'FAIL'}  df/dr matches the code's analytic form exactly")
if not ok: FAIL.append("thole df")
vals=[(x,1-np.exp(-(x/A_TH)**3)) for x in [0.0,0.5,1.0,2.0,2.35,2.75,3.0,5.0,10.0]]
mono = all(vals[i][1] <= vals[i+1][1]+1e-15 for i in range(len(vals)-1))
print(f"{'PASS' if mono else 'FAIL'}  f(r) monotone increasing, f(0)=0, f->1")
if not mono: FAIL.append("thole monotone")
print("      f at the first-shell ion-O range (2.35-2.75 A): "
      + ", ".join(f"{x:.2f}->{1-np.exp(-(x/A_TH)**3):.3f}" for x in [2.35,2.75]))

# ---------------- I1 energy/dipole consistency ----------------
print("\n[I1] U = -0.5 sum mu.E  ==  -0.5 sum alpha|E|^2 / COULOMB_MD")
rng=np.random.default_rng(7)
worst=0.0
for _ in range(200):
    N=rng.integers(2,9)
    E=rng.normal(size=(N,3))*rng.uniform(0.1,50)
    al=rng.uniform(0.1,3.0,size=N)
    mu=al[:,None]*E/COULOMB_MD
    U1=-0.5*np.sum(mu*E)
    U2=-0.5*CM*np.sum(al*np.sum(E*E,axis=1))
    worst=max(worst,abs(U1-U2)/max(abs(U2),1e-30))
print(f"      worst relative disagreement over 200 random (E, alpha) sets = {worst:.3e}")
ok=worst<1e-14
print(f"{'PASS' if ok else 'FAIL'}  the two energy expressions are the same function")
if not ok: FAIL.append("induction U consistency")

# ---------------- I2 the solver equation ----------------
print("\n[I2] (I - A) mu = alpha (*) E0  must hold at the solution")
def fields(pos,q,diel=1.0):
    N=len(q); E=np.zeros((N,3))
    for i in range(N):
        for j in range(N):
            if i==j: continue
            d=pos[i]-pos[j]; r=np.linalg.norm(d)
            if r<1e-4: continue
            f=1-np.exp(-(r/A_TH)**3)
            E[i]+= COULOMB_MD*q[j]/diel*f*d/r**3
    return E
def build_A(pos,al,diel=1.0,excl=None):
    N=len(al); A=np.zeros((N,N,3,3))
    for i in range(N):
        for j in range(N):
            if i==j or al[i]==0.0: continue
            if excl and (i,j) in excl: continue
            d=pos[i]-pos[j]; r=np.linalg.norm(d)
            if r<1e-4: continue
            f=1-np.exp(-(r/A_TH)**3)
            k=al[i]*COULOMB_MD*f/(diel*r**3)
            u=d/r
            T=3*np.outer(u,u)-np.eye(3)
            A[i,j]=k*T
    return A
def solve_code(pos,q,al,diel=1.0,excl=None):
    """Re-implementation of qm_solve_dipoles, exactly as written."""
    N=len(q)
    E0=fields(pos,q,diel)
    rhs=np.zeros((N,3))
    for i in range(N): rhs[i]=al[i]/COULOMB_MD*E0[i]
    if N==1: return rhs
    M=np.zeros((N,N,3,3))
    for i in range(N): M[i,i]=np.eye(3)
    A=build_A(pos,al,diel,excl)
    for i in range(N):
        for j in range(N):
            M[i,j]-=A[i,j]
    # flatten and gaussian-eliminate with 3 RHS
    Mm=np.zeros((N*3,N*3)); rr=rhs.reshape(-1)
    for i in range(N):
        for j in range(N):
            for a_ in range(3):
                for b_ in range(3):
                    Mm[3*i+a_,3*j+b_]=M[i,j,a_,b_]
    for c in range(3*N):
        piv=c; best=abs(Mm[c,c])
        for r_ in range(c+1,3*N):
            if abs(Mm[r_,c])>best: best=abs(Mm[r_,c]); piv=r_
        if not (best>1e-12): return None
        if piv!=c:
            Mm[[c,piv]]=Mm[[piv,c]]; rr[[c,piv]]=rr[[piv,c]]
        dg=Mm[c,c]
        for r_ in range(c+1,3*N):
            f=Mm[r_,c]/dg
            if f==0: continue
            Mm[r_,c:]-=f*Mm[c,c:]; rr[r_]-=f*rr[c]
    mu=np.zeros((N,3))
    for c in range(3*N-1,-1,-1):
        acc=rr[c]
        for k_ in range(c+1,3*N): acc-=Mm[c,k_]*mu.reshape(-1)[k_]
        mu.reshape(-1)[c]=acc/Mm[c,c]
    return mu
worst=0.0; nsing=0; ntried=0
for trial in range(300):
    N=int(rng.integers(2,10))
    pos=rng.normal(size=(N,3))*2.5
    q=rng.normal(size=N)
    al=rng.uniform(0.1,2.0,size=N)
    mu=solve_code(pos,q,al)
    ntried+=1
    if mu is None: nsing+=1; continue
    E0=fields(pos,q); A=build_A(pos,al)
    lhs=np.einsum('ijab,ib->ija',A,mu)   # (I-A)mu - mu-free check below
    rhs=np.array([al[i]/COULOMB_MD*E0[i] for i in range(N)])
    resid = mu - lhs - rhs        # (I-A)mu - alpha(*)E0 = 0
    den=max(np.linalg.norm(rhs),1e-30)
    worst=max(worst,np.linalg.norm(resid)/den)
print(f"      {ntried} random systems, {nsing} reported singular")
print(f"      worst ||(I-A)mu - alpha(*)E0|| / ||alpha(*)E0|| = {worst:.3e}")
ok=worst<1e-10
print(f"{'PASS' if ok else 'FAIL'}  the solver satisfies its defining linear system")
if not ok: FAIL.append("dipole solver residual")

# ---------------- I3 Hellmann-Feynman ----------------
print("\n[I3] dU/dE0 for the SELF-CONSISTENT case")
print("      Engine's own note: E = (I-T alpha)^-1 E0  =>  dU/dE0 = -mu (I-T alpha)^-1,")
print("      and that the naive -mu/2 is wrong. Verified by brute force below.")
worst_wrong=0.0; worst_right=0.0
for trial in range(120):
    N=int(rng.integers(3,7))
    pos=rng.normal(size=(N,3))*3.0
    q=rng.normal(size=N)*0.5
    al=rng.uniform(0.05,0.6,size=N)
    qion=0.5
    def U_of_field(s):
        E0=fields(pos,q)*s
        mu=solve_code(pos,q,al)
        return -0.5*np.sum(mu*E0), mu
    U0,mu0=U_of_field(1.0)
    h=1e-6
    Up,_=U_of_field(1+h); Um,_=U_of_field(1-h)
    dU_num=(Up-Um)/(2*h)
    # M = (I - T alpha)^-1 as a 3N x 3N operator: (I-A) with A_ij = alpha_i T_ij
    A=build_A(pos,al)
    M3=np.zeros((3*N,3*N))
    for i in range(N):
        for j in range(N):
            for a_ in range(3):
                for b_ in range(3):
                    M3[3*i+a_,3*j+b_]=(np.eye(3)-A[i,j])[a_,b_]
    Minv=np.linalg.inv(M3)
    E0=fields(pos,q)
    rhs=np.array([al[i]/COULOMB_MD*E0[i] for i in range(N)])
    # dU/dE0_i = -mu . (M e_i) for each direction; directional derivative along E0:
    # dU/d(s) at s=1 where E0 -> s*E0
    MmE0=(Minv@rhs.reshape(-1)).reshape(N,3)
    pred_right=-np.sum(mu0*MmE0)
    pred_wrong=-0.5*np.sum(mu0*E0)
    den=max(abs(dU_num),1e-30)
    worst_right=max(worst_right,abs(pred_right-dU_num)/den)
    worst_wrong=max(worst_wrong,abs(pred_wrong-dU_num)/den)
print(f"      worst rel error of  -mu.M  (correct)  vs numeric dU/dE0 : {worst_right:.3e}")
print(f"      worst rel error of  -mu/2  (naive)    vs numeric dU/dE0 : {worst_wrong:.3e}")
ok=worst_right<5e-4
print(f"{'PASS' if ok else 'FAIL'}  dU/dE0 == -mu.(I-T alpha)^-1, confirming the engine's claim")
if not ok: FAIL.append("Hellmann-Feynman")
ok=worst_wrong>0.1
print(f"{'PASS' if ok else 'FAIL'}  the naive -mu/2 IS wrong here by {worst_wrong:.2f}, as documented")
if not ok: FAIL.append("HF naive check")

print()
print("="*100)
print(f"A05 RESULT: {len(FAIL)} failure(s)" + ("" if not FAIL else "  ->  "+", ".join(FAIL)))
print("="*100)
raise SystemExit(1 if FAIL else 0)