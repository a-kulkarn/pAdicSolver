module HeckeExt

export broadcast, iterate, collect, matrix, rectangular_solve, my_nullspace, eigen, eigvecs, eigspaces, charpoly, MyEigen

export /, valuation, abs, modp, test_rings, rand, rand_padic_int, random_test_matrix, padic_qr, inverse_iteration, iseigenvector, singular_values

include("matrix_util_ext.jl")
include("padic_util.jl")

println("HeckeExt loaded!")

end
