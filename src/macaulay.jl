export macaulay_mat, solve_macaulay, coefficient_matrix

using LinearAlgebra

function is_not_homogeneous(p)
    L = [total_degree(t) for t in Hecke.terms(p)]
    return maximum(L) != minimum(L)
end

#
# Try to only generate the list of monomials once. 
#

# Creates the macaulay matrix of the polynomial system P.
#
# Should return
#    --- macaulay matrix
#    --- monomial "array"  (could be stored as a dictionary to enable indexing D[m]).

## Let's make a different implementation later


export macaulay_mat

## Present issues.
# Performance is not good, and garbage collector runs way too much.
# bizzare reversal should be made more robust.
# homogeneity is not handled.
#
function macaulay_mat(P::Array{Hecke.Generic.MPoly{T},1},
                      X::Array{Hecke.Generic.MPoly{T},1}, rho, ish) where T <: Hecke.RingElem

    degrees = unique!(map(p->total_degree(p),P))
    monomial_set    = Set{Hecke.Generic.MPoly{T}}()
    mult_monomials  = Array{Array{Hecke.Generic.MPoly{T}}}(undef, maximum(degrees))
    
    for d in degrees
        if ish
            mult_monomials[d] = monomials_of_degree(X, rho-d)
        else
            mult_monomials[d] = monomials_of_degree(X, 0:rho-d)
        end
    end
    @time for p in P
        for m in mult_monomials[total_degree(p)]
            push!(monomial_set, monomials(m*p)...)
        end
    end

    # The method "isless" is defined in AbstractAlgebra. By default Julia will use this to sort.
    monomial_set = collect(monomial_set)
    sort!(monomial_set, rev=true)
    monomial_dict = Dict(monomial_set[i]=>i for i=1:length(monomial_set))
    
    # Compute Macaulay stuff later
    # , store the result as a sparse_row.
    R = base_ring(parent(P[1]))
    
    macaulay_matrix = sparse_matrix(R)
    @time for p in P
        for m in mult_monomials[total_degree(p)]

            srow = sparse_row( R, [monomial_dict[mon] for mon in monomials(m*p)],
                                collect(coeffs(p)) )
            push!(macaulay_matrix, srow)
        end
    end

    return macaulay_matrix, monomial_dict
end

# Takes a list of polynomials and a basis of monomials and
# returns a matrix of coefficients corresponding to the
# monomial basis.
#
# L -- list of monomials
#import MultivariatePolynomials.coefficients
function coefficient_matrix(P::Vector, L)
    return Array(transpose(hcat([coeff(p, L) for p in P]...)))
end


## ***************************************************************************************    
# Main solver function

    # calls to:
    # macaulay_mat
    # mult_matrix
    # eigdiag
## ***************************************************************************************

"""
    solve_macaulay(P, X;
                        rho =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                        eigenvector_method="power",
                        test_mode=false )

Solve a 0-dimensional system of polynomial equations. (Presently, only over Qp)

#-------------------

INPUTS:
P   -- polynomial system, a Vector of AbstractAlgebra polynomials.
X   -- variables in the polynomial system
rho -- monomial degree of the system. Default is the macaulay degree.
eigenvector_method -- Strategy to solve for eigenvectors. Default is power iteration.

"""
function solve_macaulay(P, X;
                        rho =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                        eigenvector_method="power",
                        test_mode=false )

    # This solve function could be made to work with polynomials with FlintRR coefficients
    # as well, though this requires managing the type dispatch a bit and remodeling the
    # old DynamicPolynomials based subfunctions.
    
    println()
    println("-- Degrees ", map(p->total_degree(p),P))
    
    ish = !any(is_not_homogeneous, P)
    println("-- Homogeneity ", ish)

    t0 = time()
    R, L = macaulay_mat(P, X, rho, ish)
    # R = matrix(R)
    
    # Not efficient. This is intermediate code to check tests.
    L0 = monomials_divisible_by_x0(L, ish)
    
    println("-- Macaulay matrix ", size(R,1),"x",size(R,2),  "   ",time()-t0, "(s)"); t0 = time()
    @time N = nullspace(R)[2]

    println("-- -- rank of Macaulay matrix ", size(R,2) - size(N,2))
    println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    # The idea of the QR step is two-fold:
    # 1: Choose a well-conditioned *monomial* basis for the algebra from a given spanning 
    #    set (here, IdL0).
    #    This is accomplished by pivoting. The columns corresponding to F.p[1:size(N,2)] form
    #    a well-conditioned submatrix.
    #
    # 2: Present the algebra in Q-coordinates, which has many zeroes. Note that the choice of
    #    coordinates is not important in the final step, when the eigenvalues are calulated.
    #
    F, Nr = iwasawa_step(N, L0)
    B = permute_and_divide_by_x0(L0, F, ish)

    println("-- Qr basis ",  length(B), "   ",time()-t0, "(s)"); t0 = time()

    
    M = mult_matrices(B, X, Nr, L, ish)
    println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()

    if test_mode
        println("TESTING MODE: Computation incomplete. Returning partial result.")
        return M, F, B, N, Nr, R, IdL0, Idx
    end

    Xi = normalized_simultaneous_eigenvalues(M,ish, eigenvector_method)
    println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()

    # In the affine system, the distinguished monomial (i.e, "1" for that case) does not correspond
    # to a coordinate.
    if ish return Xi else return  Xi[:,2:size(Xi,2)] end
end


"""
    iwasawa_step

    Return the QR-factorization object (For PN₀P' = QR, return <inv(P)Q, R, P'>)
    together with  Nr = N*inv(inv(P)Q)^T.
"""

function iwasawa_step(N :: Array{T,2} where T <: Number, L0)
    F = qr( Array(transpose(N[IdL0,:])) , Val(true))
    return F, N*F.Q
end

function iwasawa_step(N :: Array{padic,2} , IdL0)
    
    F = padic_qr( transpose(matrix(N[IdL0,:])) , col_pivot=Val(true))
    Qinv = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv= Dory.inverse_permutation(F.p)

    X = Qinv[Fpinv,:].entries    
    Farr = QRPadicArrayPivoted( (F.Q.entries)[Fpinv,:], F.R.entries, F.q)
    
    return Farr, N*X
end


function iwasawa_step(N :: Hecke.Generic.MatSpaceElem{padic} , L0)

    F = padic_qr( transpose(N[ sort(collect(values(L0))),:]) , col_pivot=Val(true))
    Qinv = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv= Dory.inverse_permutation(F.p)

    X = Qinv[Fpinv,:]
    Farr = QRPadicArrayPivoted( (F.Q.entries)[Fpinv,:], F.R.entries, F.q)
    
    return Farr, N*X
end
