import itertools, numpy as np

# F_9 = F_3[i], i^2 = -1
def add(x,y): return ((x[0]+y[0])%3,(x[1]+y[1])%3)
def mul(x,y): return ((x[0]*y[0]-x[1]*y[1])%3,(x[0]*y[1]+x[1]*y[0])%3)
def neg(x): return ((-x[0])%3,(-x[1])%3)
def sub(x,y): return add(x,neg(y))
els=[(a,b) for a in range(3) for b in range(3)]
idx={e:i for i,e in enumerate(els)}
def pw(x,k):
    r=(1,0)
    for _ in range(k): r=mul(r,x)
    return r
sq=sorted({pw(x,2) for x in els if x!=(0,0)})
n=9
A=np.zeros((n,n),dtype=int)
for x in els:
    for y in els:
        if x!=y and sub(x,y) in sq:
            A[idx[x],idx[y]]=1
print("Paley(9) adjacency:")
print(A)
print("srg(9,4,1,2)?", np.array_equal(A@A, 2*np.eye(9,dtype=int)-A+2*np.ones((9,9),dtype=int)))

# candidate automorphisms x -> a*x + b  and  x -> a*x^3 + b
cands=[]
for a in sq:
    for b in els:
        for power in (1,3):
            perm=[idx[add(mul(a,pw(x,power)),b)] for x in els]
            if len(set(perm))==9 and all(A[perm[i],perm[j]]==A[i,j] for i in range(9) for j in range(9)):
                cands.append(perm)
cands=[list(p) for p in {tuple(c) for c in cands}]
print("num automorphisms found:", len(cands))

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
    orbs=orbits(p)
    OI={}
    for i,o in enumerate(orbs):
        for v in o: OI[v]=i
    m=len(orbs)
    M=np.zeros((m,m),dtype=int)
    for i,o in enumerate(orbs):
        x=o[0]
        for y in range(9):
            if A[x,y]: M[i,OI[y]]+=1
    return orbs,M

# find order-3 automorphism with 3 fixed vertices
print("\n--- testing framework identity (M^2)_PQ = 2*delta - M_PQ + 2|P| on Paley(9) ---")
for p in cands:
    orbs,M=buildM(A,p)
    sizes=[len(o) for o in orbs]
    ok=True
    for P in range(len(orbs)):
        for Q in range(len(orbs)):
            rhs = (2 if P==Q else 0) - M[P,Q] + 2*sizes[P]
            if (M@M)[P,Q]!=rhs: ok=False
    o=order(p)
    nfix=sum(1 for oo in orbs if len(oo)==1)
    print(f"order={o} fixedpts={nfix} orbit_sizes={sizes} identity_holds={ok}")
