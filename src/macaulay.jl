export macaulay_mat, qr_basis, solve_macaulay, coefficient_matrix

using LinearAlgebra
using DynamicPolynomials

function is_not_homogeneous(p)
    L = [total_degree(t) for t in Hecke.terms(p)]
    return maximum(L) != minimum(L)
end


# Creates the macaulay matrix of the polynomial system P.
function macaulay_mat(P, L::AbstractVector, X, ish = false )
    d = maximum([total_degree(m) for m in L])
    if ish
        Q = [monomials_of_degree(X,d-total_degree(P[i])) for i in 1:length(P)]
    else
        Q = [monomials_of_degree(X,0:d-total_degree(P[i])) for i in 1:length(P)]
    end

    ### this looks like it can be optimized a bit.
    M = []
    for i in 1:length(P)
        for m in Q[i]
            push!(M,P[i]*m)
        end
    end
    ###    
    return coefficient_matrix(M, L)
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



# Main solver function, for us anyways.
# calls to:
# macaulay_mat
# qr_basis
# mult_matrix
# eigdiag
#-------------------
# INPUTS:
# P   -- polynomial system
# X   -- variables in the polynomial system
# rho -- monomial degree of the system. Default is the macaulay degree.
# eigenvector_method -- Strategy to solve for eigenvectors. Default is power iteration.
#
function solve_macaulay(P, X;
                        rho =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                        eigenvector_method="power",
                        test_mode=false )
    
    println()
    println("-- Degrees ", map(p->total_degree(p),P))
    
    ish = !any(is_not_homogeneous, P)
    println("-- Homogeneity ", ish)
    if ish
        L = [m for m in monomials_of_degree(X, rho)]
    else
        L = [m for m in monomials_of_degree(X, 0:rho)]
    end
    # We also specifically designate the "monomial" of x0 in the computations.
    # in the affine case, the monomial x0 is just "1", in which case we mean take the
    # monomials whose degrees are not maximal.
    #
    # The KEY property of this monomial basis L0 is that for any element b, xi*b remains inside the
    # larger monomial basis L. That is,
    #                                         X ⋅ L0 ⊂ L
    Idx = idx(L)
    L0 = monomials_divisible_by_x0(L, ish)
    IdL0 = [get(Idx, m,0) for m in L0]
    
    # START MAIN SOLVER
    t0 = time()
    println("-- Monomials ", length(L), " degree ", rho,"   ",time()-t0, "(s)"); t0 = time()

    R = macaulay_mat(P, L, X, ish)
    println("-- Macaulay matrix ", size(R,1),"x",size(R,2),  "   ",time()-t0, "(s)"); t0 = time()
    N = nullspace(R)

    println("-- -- rank of Macaulay matrix ", size(R,2) - size(N,2))
    println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    # The idea of the QR step is two-fold:
    # 1: Choose a well-conditioned *monomial* basis for the algebra from a given spanning set (here, IdL0).
    #    This is accomplished by pivoting. The columns corresponding to F.p[1:size(N,2)] form a well-conditioned
    #    submatrix.
    #
    # 2: Present the algebra in Q-coordinates, which has many zeroes. Note that the choice of coordinates
    #    is not important in the final step, when the eigenvalues are calulated.
    #
    F, Nr = iwasawa_step(N, IdL0)
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

function iwasawa_step(N :: Array{T,2} where T <: Number, IdL0)
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
