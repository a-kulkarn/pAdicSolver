
X = @Ring x1 x2 x3
n = length(X)

d = 2
M = monomials(X,0:d)
s = length(M)

P = (2*rand(n,s)-fill(1.0,n,s))*M
Xi = solve_macaulay(P,X)

Er = rel_error(P,Xi)
println("-- Rel error: ", norm(Er,Inf))
