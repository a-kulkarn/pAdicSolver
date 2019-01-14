export macaulay_mat, qr_basis, solve_macaulay

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

    ### this looks like it can be inlined
    M = []
    for i in 1:length(P)
        for m in Q[i]
            push!(M,P[i]*m)
        end
    end
    ###
    
    matrix(M,idx(L))
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
function solve_macaulay(P, X, rho =  sum(degree(P[i])-1 for i in 1:length(P)) + 1 )
    println()
    println("-- Degrees ", map(p->degree(p),P))
    ish = !any(is_not_homogeneous, P)
    println("-- Homogeneity ", ish)
    if ish
        L = [m for m in monomials(X, rho)]
    else
        L = [m for m in monomials(X, 0:rho)]
    end
    t0 = time()
    println("-- Monomials ", length(L), " degree ", rho,"   ",time()-t0, "(s)"); t0 = time()

    R = macaulay_mat(P, L, X, ish)
    println("-- Macaulay matrix ", size(R,1),"x",size(R,2),  "   ",time()-t0, "(s)"); t0 = time()

    N = nullspace(R)
    println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    B, Nr = qr_basis(N, L, ish)
    println("-- Qr basis ",  length(B), "   ",time()-t0, "(s)"); t0 = time()

    M = mult_matrix(B, X, Nr, L, ish)
    println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()

    Xi = eigdiag(M)
    println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()

    # This should not be at top-level. Fix needed.
    if (!ish)
        for i in 1:size(Xi,1) Xi[i,:]/=Xi[i,1] end
        Xi = Xi[:,2:size(Xi,2)]
    else
        for i in 1:size(Xi,1) Xi[i,:]/=norm(Xi[i,:]) end
    end
    Xi
end
