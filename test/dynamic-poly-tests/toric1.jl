using AlgebraicSolvers

X = @Ring x1 x2
n = length(X)

A1 = [one(X),x1,x2,x1*x2,x1^2,x1^2*x2]
A2 = [one(X),x1,x1^2]

p1= rand(length(A1))'*A1
p2= rand(length(A2))'*A2
P= [p1,p2]

Xi = solve_toric(P,X)

println("-- sol ", Xi)
Er = rel_error(P,Xi,X)
println("-- Rel error: ", norm(Er,Inf));

