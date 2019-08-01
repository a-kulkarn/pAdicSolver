
using Dory
using pAdicSolver

Qp = PadicField(491,10)

R, (x1,x2,x3) = PolynomialRing(Qp, ["x1", "x2","x3"])
n = length(gens(R))-1

d = 1
M = monomials_of_degree(gens(R), d)
s = length(M)

P = hcat( [ [ 2*rand_padic_int(Qp)- 1 for i in 1:n] for j in 1:s]... )*M

sol = solve_macaulay(P)

println("\n-- sol ")
println(sol,"\n")

Er = rel_error(P,sol)
println("-- Rel error: ")
display(Er)
println()
