
X = AS.@Ring x1 x2 x3
n = length(X)-1

d = 1
M = AS.monomials(X,d)
s = length(M)

#P = [x1^2+1.0, x2^2-2.0]
P = (2*rand(n,s)-fill(1.0,n,s))*M

Xi = AS.solve_macaulay(P,X)
println("-- sol ", Xi)
Er = AS.rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

