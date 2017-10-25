using AlgebraicSolvers

include("solver_dense.jl")

X = @Ring x1 x2 x3 x4

init(X)
    
d = [2,4,6,8,10,12,14,16] 
t = [solver_dense(k,X) for k in d]

save("dense4d.res", length(X),d,t)
