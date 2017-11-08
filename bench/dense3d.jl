using AlgebraicSolvers
include("solver_dense.jl")

X = @Ring x1 x2 x3

init(X)

dg = [2,4,6,8,10,12,14,16,18,20] 
t  = [solver_dense(d,X) for d in dg]

save("dense3d.res", length(X), dg, t)




