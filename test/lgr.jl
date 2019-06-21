
using Dory
using pAdicSolver

Qp = PadicField(89,6)

R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40) =
    PolynomialRing(Qp,
                   ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12",
                    "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23",
                    "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", "x32", "x33", "x34",
                    "x35", "x36", "x37", "x38", "x39", "x40"])

# System of polynomial equations for the Lagrangian grassmanian in variables (x1, ..., x40)
include(".lgr_eqns.jl")

# We also need the random lines as well.
function random_linear()
    M = matrix( hcat([[ R(p)^rand(0:15) for i in 1:12] for j in 1:2]...))
    return (var_vec*M)[1,1]
end


#sol = solve_macaulay(P,X)

# println("\n-- sol ")
# println(sol,"\n")

try
    Er = rel_error(P,sol)
    println("-- Rel error: ")
    display(Er)
    println()
catch e
end
