#using AlgebraicSolvers
#using LinearAlgebra

X = AS.@Ring x0 x1 x2
n = length(X)-1

d  = 5
M  = AS.monomials(X,d)
s  = length(M)

P = (2*rand(n,s)-fill(1.0,n,s))*M

Xi = AS.solve_macaulay(P,X)
Er = AS.rel_error(P,Xi)
println("-- Rel. error: ", norm(Er,Inf))
