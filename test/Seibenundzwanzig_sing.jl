

using Dory
using pAdicSolver
using Singular
using Test

# We need more precision when using the Groebner basis methods.
Qp = PadicField(89,70)

R, (p12, p13, p14, p23, p24) = Singular.PolynomialRing(Singular.QQ,
                                                       ["p12", "p13", "p14", "p23", "p24"])
n = length(gens(R))

include(".cubic_surface_eqns.jl")
include(".cubic_surface_eqns_tropical_sol.jl")

I = Singular.Ideal(R, P)
sol = padic_solutions(I, Qp, eigenvector_method="tropical")

@test true_sol_27_lines == sol
