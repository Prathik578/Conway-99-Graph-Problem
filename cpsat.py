from ortools.sat.python import cp_model
import sys
TARGET=int(sys.argv[1]); TL=int(sys.argv[2]) if len(sys.argv)>2 else 200
NT=12; NR=20; N=32
grp=[0]*4+[1]*4+[2]*4+[-1]*20
isT=[i<NT for i in range(N)]
M=cp_model.CpModel()
m={}
for P in range(N):
    for Q in range(P,N):
        v=M.NewIntVar(0,2,f"m{P}_{Q}"); m[(P,Q)]=v; m[(Q,P)]=v
for P in range(N):
    if isT[P]: M.Add(m[(P,P)]==0)
    else: M.AddAllowedAssignments([m[(P,P)]],[(0,),(2,)])
# row sums
for P in range(N):
    M.Add(sum(m[(P,Q)] for Q in range(N))==(13 if isT[P] else 14))
# folded quotient identity
for P in range(N):
    for Q in range(N):
        sf = 1 if (isT[P] and isT[Q] and grp[P]==grp[Q]) else 0
        terms=[]
        for k in range(N):
            z=M.NewIntVar(0,4,f"z{P}_{k}_{Q}"); M.AddMultiplicationEquality(z,[m[(P,k)],m[(k,Q)]]); terms.append(z)
        M.Add(sum(terms)==(12 if P==Q else 0)-m[(P,Q)]+6-3*sf)
# hints
for P in range(NT,N): M.Add(sum(m[(P,R)] for R in range(NT))==6)
for R in range(NT): M.Add(sum(m[(P,R)] for P in range(NT,N))==10)
# t count
# t fixed by symmetry breaking below
# symmetry breaking: fix free orbits 12..(12+TARGET-1) to be triangles (WLOG relabel)
for P in range(NT,NT+TARGET): M.Add(m[(P,P)]==2)
for P in range(NT+TARGET,N): M.Add(m[(P,P)]==0)
solver=cp_model.CpSolver(); solver.parameters.max_time_in_seconds=TL
solver.parameters.num_search_workers=8
st=solver.Solve(M)
print("t=",TARGET,"status:",solver.StatusName(st))
if st in (cp_model.OPTIMAL,cp_model.FEASIBLE):
    val=lambda P,Q: solver.Value(m[(P,Q)])
    S=sum(val(P,Q) for P in range(NT,N) for Q in range(NT,N) if P<Q)
    S2=sum(val(P,Q)**2 for P in range(NT,N) for Q in range(NT,N) if P<Q)
    N2f=sum(1 for P in range(NT,N) for R in range(NT) if val(P,R)==2)
    print("S=",S,"S2=",S2,"N2f=",N2f,"| 80-t=",80-TARGET,"108-3t=",108-3*TARGET)
