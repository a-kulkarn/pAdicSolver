

using Hecke
# Needs the new matrix utilities as well.

##############################################################################################
#                                                                                            #
#                          Basic extension to padic numbers                                  #
#                                                                                            #
##############################################################################################

import Base: /, +, abs 
import Hecke.valuation

function +(x::padic) return x end
function /(x::padic,y::padic) return x//y end


# Potential fix for the bug with the valuation function
# Note: there may be issues if Hecke depends on the
# valuation being integer type
function valuation(x::padic)
    if iszero(x)
        return Inf
    end
    return x.v
end

function abs(x::padic)
    p = parent(x).p
    return Float64(p)^(-valuation(x))
end

function modp(x::padic)
    Fp = ResidueRing(FlintZZ,parent(x).p)
    return Fp(lift(x))
end

## Test utilities
function test_rings()
    return Qp = FlintPadicField(7,20), ResidueRing(FlintZZ,7)
end

import Base.rand
function rand(Qp::FlintPadicField)
    p = Qp.p
    N = Qp.prec_max
    return Qp(rand(1:p^N))
end

function random_test_matrix(Qp)
    return matrix(Qp, rand.(fill(1:7^20),4,4))
end


##############################################################################################
#                                                                                            #
#                               p-adic linear algebra                                        #
#                                                                                            #
##############################################################################################

# Compute a factorization of a padic matrix A = QR, where R is
# upper triangular (Borel) and Q lies in the maximally compact
# subgroup of SL(Qp) = SL(Zp).
#
# It turns out that the algorithm to do this is just the LU factorization
# with pivoting.

struct QRPadicPivoted
    Q::Hecke.Generic.MatElem{padic}
    R::Hecke.Generic.MatElem{padic}
    p::Array{Int64,1}
end

function padic_qr(A::Hecke.Generic.MatElem{padic})
    n = size(A,1)
    L= identity_matrix(A.base_ring,n)
    U= deepcopy(A)
    P= Array(1:n)
    #identity_matrix(A.base_ring,n)
    
    for k=1:size(A,2)
        norm_list = abs.((U.entries)[k:n,k])
        maxn, m = findmax( norm_list );
        if iszero(maxn) continue end
            
        m=m+k-1;
        if m!=k
            # Note: we actually want inv(P), so we should compute differently
            # interchange rows m and k in U
            temp=U[k,:];
            U[k,:]=U[m,:];
            U[m,:]=temp;
            
            # interchange rows m and k in P
            temp=P[k];
            P[k]=P[m];
            P[m]=temp;

            #swap columns corresponding to the row operations already done
            if k >= 2
                Lent = L.entries
                temp=Lent[k,1:k-1];
                Lent[k,1:k-1]=Lent[m,1:k-1];
                Lent[m,1:k-1]=temp;
            end
        end
        for j=k+1:n
            L[j,k]=U[j,k]*inv(U[k,k]);
            U[j,:]=U[j,:]-L[j,k]*U[k,:];
        end
    end
    return QRPadicPivoted(L,U,P)
end

# Needs to be more robust. Also applied to the situation A is square but not of rank 1.
#
# a slightly generalized version of solve
# WARNING: does not check if the top block is non-singular
function rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})

    m = rows(A)
    n = cols(A)
    if rows(b_input) != m
        error("`A` and `b` must have the same number of rows.")
    end
    b = deepcopy(b_input)

    if m < n
        error("Underdetermined systems not yet supported.")
    end

    F = padic_qr(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            #scale = A[i, j]
            #b.zip_row_op(i, j, lambda x, y: x - y * scale)
            b[i,:] = b[i,:] - b[j,:]* F.Q[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:cols(b)
                if !iszero(b[i, j])
                    error("The system is inconsistent.")
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # backward substitution
    #for i in range(n - 1, -1, -1):  # original python

    for i in n:-1:1
        for j in (i+1):n
            #scale = A[i, j]
            #b.zip_row_op(i, j, lambda x, y: x - y * scale)
            b[i,:] = b[i,:] - b[j,:]*F.R[i,j]
        end
        #scale = A[i, i]
        #b.row_op(i, lambda x, _: x / scale)

        if !iszero(b[i,:]) && iszero(F.R[i,i])
            error("The system is inconsistent.")
        elseif !iszero(F.R[i,i])
            b[i,:] *= inv(F.R[i,i])
        end
    end

    return b
end


# function old_version_of_rectangular_solve() 
#     if rows(M) < cols(M)
#         error("Not implemented when rows(M) < cols(M)")
#     end
#     # Extract top nxn block
#     A = M[1:cols(M),:]

    
#     x = solve(A,b[1:cols(M),:])
    
#     if iszero(M*x - b)
#         return x
#     else
#         error("Linear system does not have a solution")
#     end
# end


# Solve for an eigenvector using inverse iteration.
# Note that the algorithm will not converge to a particular vector in general, but the norm of
#
# A*w - λ*w converges to zero. Here, λ is the unique eigenvalue closest to `shift`, (if it is unique).
#
# I should really optimize and stabilize this later
function inverse_iteration!(A,shift,v)
    In = identity_matrix(A.base_ring, size(A,1))
    B = A - shift*In

    if iszero(det(B))
        println("Value `shift` is exact eigenvalue. Returning only first basis vector")
        return nullspace(B)[:,1]
    end
    
    pow = inv(B)
    for i=1:100
        v = A.base_ring.p*pow*v
    end
    return v
end

function inverse_iteration(A, shift, v)
    w = deepcopy(v)
    w = inverse_iteration!(A,shift,w)
    return w
end

"""
eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)

Compute the eigenvectors of a padic matrix iteratively.

NOTE: I should write a corresponding eigen function.
"""

import LinearAlgebra.eigvecs
function eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)

    Qp = A.base_ring
    # First, make everything in A a p-adic integer
    println("Skipping an important step because lazy...")

    # Solve the problem modulo p
    Amp = broadcast(modp, A)
    E = eigen(Amp)
    eig_pairs = [ (E.values[i], E.vectors[:,i]) for i in 1:size(E.values)[1]]
    
    println("Assuming that the roots of the characteristic polynomial modulo p are all distinct")
   
    return  hcat([ inverse_iteration(A, Qp(lift(e)), matrix(Qp,lift(v))) for (e,v) in eig_pairs]...)
end


# function for testing
"""
Computes if v is an eigenvector of A. If so, returns the eigenvalue as well.

TODO: This function needs some work
"""
function iseigenvector(A,v)
    i=1
    while i<=size(v,2)
        if !iszero(v[i,1])
            break
        end
        i+=1
    end
    if i>size(v,2)
        return false
    end
    e = (A*v)[i,1]/v[i,1]
    return iszero(A*v - (A*v)[i,1]/v[i,1]*v),e
end



#PA - eP = LU
# PAinv(P) - eI = LUinv(P)
# inv(L)PAinv(P)L = Uinv(P)L + eI
#inv(L)(PA - eP)inv(P)L = Uinv(P)L
#inv(L)PAinv(P)L  = Uinv(P)L + eI

function one_iteration(A,Q,shift)
    eI = Qp(shift)*(identity_matrix(Qp,size(A,1)))
    L,U,P = padic_qr(A-eI)
    return U*inv(P)*L + eI, Q*inv(P)*L
end


println("padic_util 'Package' loaded!")
