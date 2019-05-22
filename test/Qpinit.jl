



Qp = PadicField(7,5)
A = random_test_matrix(Qp,5)


# TODO: improve the matrix constructor so that it is actually intuitive.

# entarr =( vcat( [4 + 5*7^1  2 + 3*7^1  1 + 3*7^1  2 + 2*7^1 ],
#        [3 + 1*7^1  5  1 + 4*7^1  6 + 4*7^1 ],
#        [2*7^1  4*7^1  2  3*7^1 ],
#        [4 + 4*7^1  3 + 3*7^1  4 + 3*7^1  4 + 2*7^1 ]))

# entarr = [Qp(1)*a for a in entarr]

# A = matrix(Qp, entarr)


@time X = Dory.block_schur_form(A)
#H = Dory.hessenberg(A)

vals, spaces = Dory.power_iteration_decomposition(A, modp.(A))


display(vals[1].entries)

display(spaces[1].entries)
