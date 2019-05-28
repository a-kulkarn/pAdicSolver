


using Hecke
using Dory
using AlgebraicSolvers


# For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
#Qp = PadicField(89,30)
Qp = PadicField(41,30)

R, (p12, p13, p14, p23, p24) = PolynomialRing(Qp, ["p12", "p13", "p14", "p23", "p24"])
n = length(gens(R))

include("_cubic_surface_eqns.jl")

sol = solve_macaulay(P, gens(R), rho=5)
#matlist, F, B, N, Nr = solve_macaulay(P, gens(R), rho=4, test_mode=true);

#A = matrix(matlist[1]);
#E = eigspaces(A);

# sol = AS.solve_macaulay(P,X)

# println("\n-- sol ")
# println(sol,"\n")

Er = rel_error(P,sol)
println("-- Rel error: ")
display(Er)
println()


