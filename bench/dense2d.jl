using AlgebraicSolvers

include("solver_dense.jl")

X = @Ring x1 x2

init(X)
    
d = [2,3,4,5] 
t = [solver_dense(k,X) for k in d]

save("dense2d.res", length(X),d,t)
