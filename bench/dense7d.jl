using AlgebraicSolvers

include("solver_dense.jl")

X = @Ring x1 x2 x3 x4 x5 x6 x7

solver_dense(2,X)
println()
#solver_dense(3,X)
#println()
#solver_dense(4,X)
#println()



