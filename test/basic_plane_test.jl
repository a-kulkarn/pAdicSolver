
using Test
using Dory
using pAdicSolver

# For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
Qp = PadicField(491,3)


R, (x1,x2) = PolynomialRing(Qp, ["x1", "x2"])
n = length(gens(R))

d = 10
M = monomials_of_degree(gens(R),0:d)
s = length(M)

failed_test_count = 0

global P = hcat( [ [ 2*rand_padic_int(Qp)- 1 for i in 1:n] for j in 1:s]... )*M

sol = solve_macaulay(P)

println("\n-- sol ")
println(sol,"\n")

Er = rel_error(P,sol)
println("-- Rel error: ")
display(Er)
println()

@test iszero(Er)
