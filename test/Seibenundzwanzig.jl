

using Dory
using pAdicSolver


# For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
Qp = PadicField(89,30)
#Qp = PadicField(41,30)

R, (p12, p13, p14, p23, p24) = PolynomialRing(Qp, ["p12", "p13", "p14", "p23", "p24"])
n = length(gens(R))

include(".cubic_surface_eqns.jl")

sol = solve_macaulay(P, rho=5, groebner=false, eigenvector_method="tropical")
#matlist, F, B, N, Nr = solve_macaulay(P, gens(R), rho=4, test_mode=true);

#A = matrix(matlist[1]);
#E = eigspaces(A);

# println("\n-- sol ")
# println(sol,"\n")

