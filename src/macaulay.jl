export macaulay_mat, qr_basis, solve_macaulay, coefficient_matrix

using LinearAlgebra
using DynamicPolynomials

function is_not_homogeneous(p)
    L = [degree(t) for t in p]
    maximum(L) != minimum(L)
end


# Creates the macaulay matrix of the polynomial system P.
function macaulay_mat(P, L::AbstractVector, X, ish = false )
    d = maximum([degree(m) for m in L])
    if ish
        Q = [monomials(X,d-degree(P[i])) for i in 1:length(P)]
    else
        Q = [monomials(X,0:d-degree(P[i])) for i in 1:length(P)]
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

# AVI:
# Takes a list of polynomials and a basis of monomials and
# returns a matrix of coefficients corresponding to the
# monomial basis.
#
# AVI: Multivariate polynomials has a "coefficients" function.
# Should simplify this function a bit.
#
# AVI: The t.\alpha is the coefficient of the term.
#
# L -- list of monomials
import MultivariatePolynomials.coefficients
function coefficient_matrix(P::Vector, L)
    return Array(transpose(hcat([MultivariatePolynomials.coefficients(p, L) for p in P]...)))
end



# AVI: Main solver function, for us anyways.
# calls to:
# macaulay_mat
# qr_basis
# mult_matrix
# eigdiag
#------------
# call to the qr intrinsic uses float (not p-adic) arithmetic.
#-------------------
# INPUTS:
# P   -- polynomial system
# X   -- variables in the polynomial system
# rho -- monomial degree of the system. Default is the macaulay degree.
#
function solve_macaulay(P, X, rho =  sum(degree(P[i])-1 for i in 1:length(P)) + 1, test_mode=false )
    println()
    println("-- Degrees ", map(p->degree(p),P))
    ish = !any(is_not_homogeneous, P)
    println("-- Homogeneity ", ish)
    if ish
        L = [m for m in monomials(X, rho)]
    else
        L = [m for m in monomials(X, 0:rho)]
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
    println("-- -- rank of Macaulay matrix ", rank(matrix(R)))

    N = nullspace(R)
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
        return M, F, B, N, Nr, R
    end

    Xi = normalized_simultaneous_eigenvalues(M,ish)
    println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()

    # In the affine system, the distinguished monomial (i.e, "1" for that case) does not correspond
    # to a coordinate.
    if ish return Xi else return  Xi[:,2:size(Xi,2)] end
end


# Dispatcher for doing the right iwasawa algorithm
function iwasawa_step(N :: Array{T,2} where T <: Number, IdL0)
    F = qr( Array(transpose(N[IdL0,:])) , Val(true))
    return F, N*F.Q
end


# Dispatcher for doing the right iwasawa algorithm
function iwasawa_step(N :: Array{padic,2} , IdL0)
    F = qr( Array(transpose(N[IdL0,:])) , Val(true))
    return F, N*inv(matrix(parent(N[1,1]), Array(transpose(F.Q)))).entries
end
