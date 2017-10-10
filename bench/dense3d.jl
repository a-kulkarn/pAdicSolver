using AlgebraicSolvers
include("solver_dense.jl")

X = @Ring x1 x2 x3

init(X)

dg = [2,3,5,7,9,11,13] 
t  = [solver_dense(d,X) for d in dg]

save("dense3d.res", length(X), dg, t)




