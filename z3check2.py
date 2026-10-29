from z3 import *
import numpy as np
exec(open('paley.py').read().split('# find order-3')[0])

def order(p):
    q=list(p); o=1
    while q!=list(range(9)):
        q=[p[q[i]] for i in range(9)]; o+=1
    return o
def orbits(p):
    seen=set(); orbs=[]
    for v in range(9):
        if v in seen: continue
        o=[]; u=v
        while u not in o:
            o.append(u); seen.add(u); u=p[u]
        orbs.append(sorted(o))
    return orbs
def true_M(A,p):
    orbs=orbits(p); OI={}
    for i,o in enumerate(orbs):
        for v in o: OI[v]=i
    M=np.zeros((len(orbs),len(orbs)),dtype=int)
    for i,o in enumerate(orbs):
        x=o[0]
        for y in range(9):
            if A[x,y]: M[i,OI[y]]+=1
    return orbs,M

def feasible(p):
    orbs,Mt=true_M(A,p); k=len(orbs); sz=[len(o) for o in orbs]
    s=Solver(); m=[[Int(f"m{i}_{j}") for j in range(k)] for i in range(k)]
    for i in range(k):
        for j in range(k):
            s.add(m[i][j]>=0, m[i][j]<=sz[j])   # #neighbors of a vertex of orbit i inside orbit j <= |j|
            if sz[i]==sz[j] and i!=j: s.add(m[i][j]==m[j][i])
    for i in range(k): s.add(Sum([m[i][j] for j in range(k)])==4)
    for i in range(k):
        for j in range(k):
            s.add(Sum([m[i][x]*m[x][j] for x in range(k)])==(2 if i==j else 0)-m[i][j]+2*sz[j])
    return sz, (s.check()==sat), Mt

for p in cands:
    o=order(p)
    if o==2:
        orbs,Mt=true_M(A,p); sz=[len(o) for o in orbs]
        good=all(sum(Mt[i])==4 for i in range(len(sz))) and all((Mt@Mt)[i,j]==(2 if i==j else 0)-Mt[i,j]+2*sz[j] for i in range(len(sz)) for j in range(len(sz)))
        print(f"order=2 sizes={sz} quotient_feasible={feas} trueM_satisfies={good}")
        break
