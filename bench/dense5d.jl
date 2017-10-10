using AlgebraicSolvers

include("solver_dense.jl")

X = @Ring x1 x2 x3 x4 x5

init(X)
    
d = [2,3] 
t = [solver_dense(k,X) for k in d]

save("dense5d.res", length(X),d,t)





