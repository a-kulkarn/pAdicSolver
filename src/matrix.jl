################################################################################################
#
#  Monomial divisibility
#
################################################################################################

function is_divisible_by_x0(m)
    return total_degree(gcd(m, gens(parent(m))[1])) > 0
end

function monomials_divisible_by_x0(L,ish)
    if ish
        return filter(kvp->is_divisible_by_x0(first(kvp)), L) 
    else
        d  = maximum([total_degree(m) for m in keys(L)])
        return filter(kvp->total_degree(first(kvp))<d,L)
    end
end

function permute_and_divide_by_x0(L0,F,ish)
    B = Dict()
    m = size(F.Q,1) # The dimension of the quotient algebra.

    # Extract the column to monomial correspondence.    
    key_array = first.(sort(collect(L0), by=x->x[2]))

    if ish
        for i in 1:m
            m = copy( key_array[F.p[i]]  )
            push!(B, Dory.divexact(m, gens(parent(m))[1])=>i )
            # should test if the diag. coeff. is not small
        end
    else
        for i in 1:m
            push!(B,  key_array[F.p[i]]=>i )
            # should test if the diag. coeff. is not small
        end
    end
 
    return B
end
    


################################################################################################
#
#  Macaulay matrix and multiplication matrices.
#
################################################################################################


# INPUTS:
#  B -- basis for QR factorization
#  X -- variables for the polynomial ring
#  K -- the "R" from the QR factorization. More precisely, the transpose of R is a submatrix of K.
#  L -- Monomials of polynomial system
#  ish- the "is_homogeneous" boolian.
# OUTPUT:
# a list of matrices whose eigenvalues are the solution coordinates.
# The matrices represent multiplication-by-xi maps  ** in Q-coordinates **
function mult_matrices(B, X, K, L, ish = false)

    # For an affine system, '1' is needed as a monomial as well.
    if !ish Y = vcat([parent(X[1])(1)], X) else Y = X end

    # TODO: Best to turn this into a for loop.
    rng = base_ring(K)

    Z = fill(zero(rng), length(B), size(K,2))
    monomial_mats = fill(Z, length(Y))

    for j = 1:length(Y)
        v = Y[j]
        M = fill(zero(rng), length(B), size(K,2))

        for (m,i) in collect(B)
            k = L[m*v]
            M[i, :] = K.entries[k, :] # Might be nice to make this compatible with matrix type.
        end
        monomial_mats[j] = M
    end
    
    return monomial_mats
end


## Present issues.
# Performance is not good, and garbage collector runs way too much.
# bizzare reversal should be made more robust.
# homogeneity is not handled.
#
@doc Markdown.doc"""
    macaulay_mat(P::Array{Hecke.Generic.MPoly{T},1},
                      X::Array{Hecke.Generic.MPoly{T},1}, rho, ish) where T <: Hecke.RingElem

Constructs the sparse macaulay matrix defined by the polynomials `P` and degree bound `rho`. The argument `X`
is a list of variables of the ambient polynomial ring used to construct the multiplier monomials. 
"""
function macaulay_mat(P::Vector{T},
                      X::Vector{T}, rho, ish) where T <: Hecke.Generic.MPolyElem

    degrees = unique!(map(p->total_degree(p),P))
    monomial_set    = Set{T}()
    mult_monomials  = Array{Array{T}}(undef, maximum(degrees))
    
    for d in degrees
        if ish
            mult_monomials[d] = monomials_of_degree(X, rho-d)
        else
            mult_monomials[d] = monomials_of_degree(X, 0:rho-d)
        end
    end
    for p in P
        for m in mult_monomials[total_degree(p)]
            push!(monomial_set, monomials(m*p)...)
        end
    end

    # The method "isless" is defined in AbstractAlgebra. By default Julia will use this to sort.
    monomial_set = collect(monomial_set)
    sort!(monomial_set, rev=true)
    monomial_dict = Dict(monomial_set[i]=>i for i=1:length(monomial_set))

    # Create sparse rows for each m*p, with m a mulitplier monomial and p a polynomial.
    R = base_ring(parent(P[1]))    
    macaulay_matrix = sparse_matrix(R)
    for p in P
        for m in mult_monomials[total_degree(p)]

            srow = sparse_row(R, [monomial_dict[mon] for mon in monomials(m*p)],
                              collect(coefficients(p)))
            push!(macaulay_matrix, srow)
        end
    end

    return macaulay_matrix, monomial_dict
end
