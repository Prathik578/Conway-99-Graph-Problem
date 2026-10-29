from z3 import *
NF=3; NT=12; NR=20; N=NF+NT+NR   # 3 fixed, 12 type-1 (size3), 20 free (size3)
size=[1 if i<NF else 3 for i in range(N)]
def typ(i): return 'F' if i<3 else ('T' if i<15 else 'R')

s=Solver()
m={}
# size-3 x size-3 entries (symmetric), fixed x fixed
for i in range(N):
    for j in range(N):
        if i<NF and j<NF:
            m[(i,j)]=0 if i==j else 1
        elif i>=NF and j>=NF:
            if (j,i) in m: m[(i,j)]=m[(j,i)]
            else:
                v=Int(f"m_{i}_{j}"); s.add(v>=0,v<=2); m[(i,j)]=v
# fixed<->size3: m[z][P]=3*y, m[P][z]=y
for z in range(NF):
    for P in range(NF,N):
        y=Int(f"y_{z}_{P}"); s.add(y>=0,y<=1)
        m[(z,P)]=3*y; m[(P,z)]=y
# structure: each size3 orbit adj to exactly one fixed vertex iff type-1
for P in range(NF,N):
    tot=Sum([m[(P,z)] for z in range(NF)])
    s.add(tot==(1 if typ(P)=='T' else 0))
for z in range(NF):
    s.add(Sum([m[(P,z)] for P in range(NF,NF+NT)])==4)
# row sums = 14
for i in range(N):
    s.add(Sum([m[(i,j)] for j in range(N)])==14)
# quotient identity
for i in range(N):
    for j in range(N):
        s.add(Sum([m[(i,k)]*m[(k,j)] for k in range(N)])==(12 if i==j else 0)-m[(i,j)]+2*size[j])
s.set("timeout",200000)
r=s.check(); print("RESULT:",r)
if r==sat:
    M=s.model(); val=lambda i,j: M.evaluate(m[(i,j)]).as_long()
    S=sum(val(i,j) for i in range(15,35) for j in range(15,35) if i<j)
    S2=sum(val(i,j)**2 for i in range(15,35) for j in range(15,35) if i<j)
    t=sum(1 for i in range(15,35) if val(i,i)==2)
    print("t=",t,"S=",S,"S2=",S2,"| 80-t=",80-t,"108-3t=",108-3*t)
