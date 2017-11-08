using AlgebraicSolvers

include("solver_dense.jl")

X = @Ring x1 x2

init(X)
    
d = [2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100] 
t = [solver_dense(k,X) for k in d]

save("dense2d.res", length(X),d,t)
