from z3 import *
import numpy as np
exec(open('paley.py').read().split('# find order-3')[0])  # A, cands, els, idx

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

# GENERIC quotient feasibility check for a given orbit-size profile of srg(9,4,1,2): A^2 = 2I - A + 2J
def feasible_paley(p):
    orbs,Mtrue=true_M(A,p)
    k=len(orbs); sizes=[len(o) for o in orbs]
    s=Solver()
    m=[[None]*k for _ in range(k)]
    for i in range(k):
        for j in range(i,k):
            lo,hi=(0,0) if sizes[i]==1 and sizes[j]==1 else (0,min(sizes[i],sizes[j]) if sizes[i]!=sizes[j] else 0)
            v=Int(f"m{i}_{j}")
            if i==j: s.add(v>=0, v<= (0 if sizes[i]==1 else 2))
            else: s.add(v>=0, v<= min(sizes[i],sizes[j]))
            m[i][j]=v
            if i!=j: m[j][i]=v
    for i in range(k):
        s.add(Sum([m[i][j] for j in range(k)])==4)   # row sum = k = 4
    for i in range(k):
        for j in range(k):
            s.add(Sum([m[i][x]*m[x][j] for x in range(k)])==(2 if i==j else 0)-m[i][j]+2*sizes[j])
    r=s.check()
    ok = (r==sat)
    # does a true model exist among solutions?
    return sizes, ok, Mtrue

for p in cands[:200]:
    orbs,Mt=true_M(A,p)
    sz,feas,_=feasible_paley(p)
    o=order(p)
    if o in (1,2,3) :
        # verify true M satisfies constraints
        good=True
        for i in range(len(sz)):
            if sum(Mt[i])!=4: good=False
        for i in range(len(sz)):
            for j in range(len(sz)):
                if (Mt@Mt)[i,j]!=(2 if i==j else 0)-Mt[i,j]+2*sz[j]: good=False
        print(f"order={o} sizes={sz} quotient_feasible={feas} trueM_satisfies={good}")
        break
