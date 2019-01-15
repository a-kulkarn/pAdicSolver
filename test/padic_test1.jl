
X = @Ring x1 x2
n = length(X)

d = 2
M = monomials(X,0:d)
s = length(M)

#P = [x1^2+1.0, x2^2-2.0]
P = (2*rand(n,s)-fill(1.0,n,s))*M

Xi = solve_macaulay(P,X)
println("-- sol ", Xi)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

