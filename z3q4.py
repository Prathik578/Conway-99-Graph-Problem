from z3 import *
import sys
TARGET=int(sys.argv[1])
# 32 orbits: 0-3 adj z0, 4-7 adj z1, 8-11 adj z2 (type-1); 12-31 free
NT=12; NR=20; N=32
grp=[0]*4+[1]*4+[2]*4+[-1]*20
isT=[i<12 for i in range(N)]
s=Solver()
m=[[None]*N for _ in range(N)]
for i in range(N):
    for j in range(i,N):
        v=Int(f"m{i}_{j}"); s.add(v>=0,v<=2); m[i][j]=v; m[j][i]=v
for P in range(N):
    if isT[P]: s.add(m[P][P]==0)
    else: s.add(Or(m[P][P]==0,m[P][P]==2))
def sameFixed(P,Q): return 1 if (isT[P] and isT[Q] and grp[P]==grp[Q] and P!=Q) else 0
# row sums over size-3 orbits
for P in range(N):
    s.add(Sum([m[P][Q] for Q in range(N)])==(13 if isT[P] else 14))
# quotient identity folded
for P in range(N):
    for Q in range(N):
        s.add(Sum([m[P][k]*m[k][Q] for k in range(N)])==(12 if P==Q else 0)-m[P][Q]+6-3*sameFixed(P,Q))
# hints
for P in range(12,32): s.add(Sum([m[P][R] for R in range(12)])==6)
for R in range(12): s.add(Sum([m[P][R] for P in range(12,32)])==10)
tri=[If(m[P][P]==2,1,0) for P in range(12,32)]
s.add(Sum(tri)==TARGET)
s.set("timeout",260000)
r=s.check(); print("t=",TARGET,"RESULT:",r)
if r==sat:
    M=s.model(); val=lambda i,j: M.evaluate(m[i][j]).as_long()
    S=sum(val(i,j) for i in range(12,32) for j in range(12,32) if i<j)
    S2=sum(val(i,j)**2 for i in range(12,32) for j in range(12,32) if i<j)
    N2f=sum(1 for P in range(12,32) for R in range(12) if val(P,R)==2)
    print("S=",S,"S2=",S2,"N2f=",N2f,"| 80-t=",80-TARGET,"108-3t=",108-3*TARGET)
