import numpy as np
exec(open('paley.py').read().split('# find order-3')[0])  # reuse graph + cands

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
def buildM(A,p):
    orbs=orbits(p); OI={}
    for i,o in enumerate(orbs):
        for v in o: OI[v]=i
    M=np.zeros((len(orbs),len(orbs)),dtype=int)
    for i,o in enumerate(orbs):
        x=o[0]
        for y in range(9):
            if A[x,y]: M[i,OI[y]]+=1
    return orbs,M

print("Testing CORRECTED identity (M^2)_PQ = 2*delta_PQ - M_PQ + 2*|Q| on Paley(9):")
seen=set()
for p in cands:
    key=tuple(sorted(map(tuple,orbits(p))))
    orbs,M=buildM(A,p)
    sizes=[len(o) for o in orbs]
    ok=all((M@M)[P,Q]==((2 if P==Q else 0)-M[P,Q]+2*sizes[Q]) for P in range(len(orbs)) for Q in range(len(orbs)))
    sig=(order(p),tuple(sizes))
    if sig in seen: continue
    seen.add(sig)
    print(f"  order={order(p)} sizes={sizes} identity_holds={ok}")
