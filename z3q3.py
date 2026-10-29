from z3 import *
import sys
TARGET=int(sys.argv[1])
NF=3; NT=12; NR=20; N=NF+NT+NR
size=[1 if i<NF else 3 for i in range(N)]
def typ(i): return 'F' if i<3 else ('T' if i<15 else 'R')
s=Solver()
m={}
for i in range(N):
    for j in range(N):
        if i<NF and j<NF: m[(i,j)]=0 if i==j else 1
        elif i>=NF and j>=NF:
            if (j,i) in m: m[(i,j)]=m[(j,i)]
            else:
                v=Int(f"m_{i}_{j}"); s.add(v>=0,v<=2); m[(i,j)]=v
for z in range(NF):
    for P in range(NF,N):
        y=Int(f"y_{z}_{P}"); s.add(y>=0,y<=1); m[(z,P)]=3*y; m[(P,z)]=y
for P in range(NF,N):
    s.add(Sum([m[(P,z)] for z in range(NF)])==(1 if typ(P)=='T' else 0))
for z in range(NF):
    s.add(Sum([m[(P,z)] for P in range(NF,NF+NT)])==4)
for i in range(N): s.add(Sum([m[(i,j)] for j in range(N)])==14)
# derived hints
for P in range(NF,NF+NT): s.add(m[(P,P)]==0)          # no type-1 triangle
for P in range(NF+NT,N): s.add(Or(m[(P,P)]==0,m[(P,P)]==2))  # free: independent or triangle
for P in range(NF+NT,N): s.add(Sum([m[(P,R)] for R in range(NF,NF+NT)])==6)
for R in range(NF,NF+NT): s.add(Sum([m[(P,R)] for P in range(NF+NT,N)])==10)
for P in range(NF+NT,N):
    for z in range(NF): s.add(Sum([m[(P,R)] for R in range(NF,NF+NT) if True])==6)  # (redundant)
# t count
tri=[If(m[(P,P)]==2,1,0) for P in range(NF+NT,N)]
s.add(Sum(tri)==TARGET)
# quotient identity
for i in range(N):
    for j in range(N):
        s.add(Sum([m[(i,k)]*m[(k,j)] for k in range(N)])==(12 if i==j else 0)-m[(i,j)]+2*size[j])
s.set("timeout",260000)
r=s.check(); print("t=",TARGET,"RESULT:",r)
if r==sat:
    M=s.model(); val=lambda i,j: M.evaluate(m[(i,j)]).as_long()
    S=sum(val(i,j) for i in range(15,35) for j in range(15,35) if i<j)
    S2=sum(val(i,j)**2 for i in range(15,35) for j in range(15,35) if i<j)
    print("S=",S,"S2=",S2,"| 80-t=",80-TARGET,"108-3t=",108-3*TARGET)
