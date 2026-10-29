# audit-checks.py -- independent arithmetic checks for the Conway 99-graph papers.
# These are AUDIT checks by an independent review, NOT the authors' repository and NOT a
# verification of the papers' proofs or hypotheses. They check: the spectrum multiplicities,
# the congruence -> candidate-set arithmetic, the triangle count, and the quotient identity
# (M^2)_ij = (k-mu)delta_ij + (lambda-mu)M_ij + mu|O_j| on Paley(9)=srg(9,4,1,2) for all
# automorphisms. Requires: python3, numpy.  Run: python3 audit-checks.py
import itertools, numpy as np
print("=== 1. spectrum multiplicities (99,14,1,2) ===")
# f3+f-4=98 ; 14+3f3-4f-4=0
for f3 in range(99):
    f4=98-f3
    if 14+3*f3-4*f4==0: print("  f3,f-4 =",f3,f4)
print("=== 2. congruences -> candidate sets ===")
T=[x for x in range(0,34) if x%21==6]; print("  T<=33, T=6 mod21:",T)
t=[x for x in range(0,21) if x%21==2]; print("  t<=20, t=2 mod21:",t)
e=[x for x in range(0,50) if x%7==0 and x%2==1]; print("  e<=49, 7|e, odd:",e)
print("=== 3. triangle count ===")
print("  edges=99*14//2 =",99*14//2," triangles=",99*14//2//3)
print("=== 4. quotient identity on Paley(9) for all automorphisms ===")
def add(x,y): return ((x[0]+y[0])%3,(x[1]+y[1])%3)
def mul(x,y): return ((x[0]*y[0]-x[1]*y[1])%3,(x[0]*y[1]+x[1]*y[0])%3)
def neg(x): return ((-x[0])%3,(-x[1])%3)
def sub(x,y): return add(x,neg(y))
def pw(x,k):
    r=(1,0)
    for _ in range(k): r=mul(r,x)
    return r
els=[(a,b) for a in range(3) for b in range(3)]; idx={e:i for i,e in enumerate(els)}
sq=sorted({pw(x,2) for x in els if x!=(0,0)})
A=np.zeros((9,9),dtype=int)
for x in els:
    for y in els:
        if x!=y and sub(x,y) in sq: A[idx[x],idx[y]]=1
cands=[]
for a in sq:
    for b in els:
        for p in (1,3):
            perm=[idx[add(mul(a,pw(x,p)),b)] for x in els]
            if len(set(perm))==9 and all(A[perm[i],perm[j]]==A[i,j] for i in range(9) for j in range(9)):
                cands.append(tuple(perm))
cands=sorted(set(cands))
def orbs(p):
    seen=set(); o=[]
    for v in range(9):
        if v in seen: continue
        c=[]; u=v
        while u not in c: c.append(u); seen.add(u); u=p[u]
        o.append(sorted(c))
    return o
allok=True
for p in cands:
    O=orbs(p); OI={}
    for i,o in enumerate(O):
        for v in o: OI[v]=i
    sz=[len(o) for o in O]; M=np.zeros((len(O),len(O)),dtype=int)
    for i,o in enumerate(O):
        for y in range(9):
            if A[o[0],y]: M[i,OI[y]]+=1
    for i in range(len(O)):
        for j in range(len(O)):
            if (M@M)[i,j]!=(2 if i==j else 0)-M[i,j]+2*sz[j]: allok=False
print("  automorphisms:",len(cands)," quotient identity holds for all:",allok)
