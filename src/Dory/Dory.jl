module Dory

export broadcast, iterate, collect, matrix, rectangular_solve, my_nullspace, eigen, eigvecs, eigspaces, charpoly, MyEigen

export /, valuation, abs, modp, test_rings, rand, rand_padic_int, random_test_matrix, padic_qr, inverse_iteration, iseigenvector, singular_values

include("dory_matrix.jl")
include("padic_util.jl")

println("\nYou've found DORY! Welcome!")
println()

end
