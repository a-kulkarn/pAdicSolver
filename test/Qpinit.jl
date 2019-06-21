
using Dory
Qp = PadicField(7,30)

# TODO: improve the matrix constructor so that it is actually intuitive.

# entarr =( vcat( [4 + 5*7^1  2 + 3*7^1  1 + 3*7^1  2 + 2*7^1 ],
#        [3 + 1*7^1  5  1 + 4*7^1  6 + 4*7^1 ],
#        [2*7^1  4*7^1  2  3*7^1 ],
#        [4 + 4*7^1  3 + 3*7^1  4 + 3*7^1  4 + 2*7^1 ]))

# entarr = [Qp(1)*a for a in entarr]

# A = matrix(Qp, entarr)


# for i=1:10

     A = random_test_matrix(Qp,50)
#     @time X,V = Dory.block_schur_form(A)
#     @assert iszero(X*V - V*A)

     display(factor(charpoly(modp.(A))))
#     #display(valuation.(X))
#     println()

# end 

    
# @assert iszero(inv(V)*X*V - A)

vals, spaces = Dory.power_iteration_decomposition(A, modp.(A))

# E = Dory._eigenspaces_by_power_iteration(A)

# if length(vals) > 0
#     display(vals[1].entries)
#     display(spaces[1].entries)
# end
    
# display(E)

nothing
