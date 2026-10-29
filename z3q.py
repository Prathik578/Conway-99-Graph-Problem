from z3 import *
import itertools, sys

# nodes: 0,1,2 = fixed vertices; 3..14 = 12 type-1 orbits; 15..34 = 20 free orbits
NF=3; NT=12; NR=20; N=NF+NT+NR
size=[1 if i<NF else 3 for i in range(N)]
def typ(i):
    if i<3: return 'F'
    if i<15: return 'T'
    return 'R'

s=Solver()
m=[[None]*N for _ in range(N)]
# symmetric integer entries for size-3 pairs
for i in range(N):
    for j in range(i,N):
        if i<3 and j<3:
            if i==j: v=Int(f"m_{i}_{j}"); s.add(v==0)
            else: v=Int(f"m_{i}_{j}"); s.add(v==1)
        elif i<3 and j>=3:
            v=Int(f"m_{i}_{j}"); s.add(v>=0, v<=3)
        else:
            v=Int(f"m_{i}_{j}"); s.add(v>=0, v<=2)
        m[i][j]=v
        if i!=j: m[j][i]=v
# cross fixed<->size3: m[z][P]=3*y, m[P][z]=y
for z in range(3):
    for P in range(3,N):
        y=Int(f"y_{z}_{P}"); s.add(y>=0,y<=1)
        s.add(m[z][P]==3*y)
        s.add(m[P][z]==y)
# each size-3 orbit adjacent to exactly one fixed vertex iff type-1
for P in range(3,N):
    tot=Sum([m[P][z] for z in range(3)])
    if typ(P)=='T': s.add(tot==1)
    else: s.add(tot==0)
for z in range(3):
    s.add(Sum([m[P][z] for P in range(3,15)])==4)   # 4 type-1 orbits per fixed vertex
# row sums = 14
for i in range(N):
    s.add(Sum([m[i][j] for j in range(N)])==14)
# quotient identity (M^2)_ij = 12 delta_ij - M_ij + 2|j|
for i in range(N):
    for j in range(N):
        s.add(Sum([m[i][k]*m[k][j] for k in range(N)]) == (12 if i==j else 0) - m[i][j] + 2*size[j])

s.set("timeout", 220000)
r=s.check()
print("result:", r)
if r==sat:
    M=s.model()
    def val(i,j): return M.evaluate(m[i][j]).as_long()
    # read off sums over free orbits
    S=sum(val(i,j) for i in range(15,35) for j in range(15,35) if i<j)
    S2=sum(val(i,j)**2 for i in range(15,35) for j in range(15,35) if i<j)
    aP=[val(i,i) for i in range(15,35)]
    t=aP.count(2)
    print("t(meta)=",t,"S=",S,"S2=",S2, "80-t=",80-t,"108-3t=",108-3*t)
