
using Dory
using pAdicSolver

Qp = PadicField(89,30)

R, (p12, p13, p14, p23, p24) = PolynomialRing(Qp, ["p12", "p13", "p14", "p23", "p24"])
n = length(gens(R))

include(".cubic_surface_eqns.jl")

rho=4;
X=gens(R);

println()
println("-- Degrees ", map(p->total_degree(p),P))

ish = false
println("-- Homogeneity ", ish)

t0 = time()
R, L = macaulay_mat(P, X, rho, ish)

println( valuation.(singular_values(matrix(R))) )

# R = matrix(R)

# F = padic_qr(transpose(R))

println("nullspace: ")
@time nu,N = nullspace(R)

size(N);

