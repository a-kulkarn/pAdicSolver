
using Dory
using pAdicSolver
using Hecke
using Test

Qp = PadicField(89,30)

R, (p12, p13, p14, p23, p24) = Hecke.PolynomialRing(Qp, ["p12", "p13", "p14", "p23", "p24"])
n = length(gens(R))

include(".cubic_surface_eqns.jl")
include(".cubic_surface_eqns_tropical_sol.jl")

#sol = solve_macaulay(P, rho=5, eigenvector_method="tropical")
sol = pAdicSolver.solve_macaulay_II(P, rho=5, eigenvector_method="tropical", verbose=false)

@test true_sol_27_lines == sol
