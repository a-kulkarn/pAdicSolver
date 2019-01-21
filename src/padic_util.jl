

using Hecke
# Needs the new matrix utilities as well.

##############################################################################################
#                                                                                            #
#                          Basic extension to padic numbers                                  #
#                                                                                            #
##############################################################################################

import Base: +, abs 
import Hecke.valuation

function +(x::padic) return x end
function /(x::padic,y::padic) return x//y end


# Potential fix for the bug with the valuation function
# Note: there may be issues if Hecke depends on the
# valuation being non-infinite.
#
function valuation(x::padic)
    if iszero(x)
        return Inf
    end
    return Int64(x.v)
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
    return matrix(Qp, rand.(fill(1:7^Qp.prec_max),4,4))
end

##############################################################################################
#                                                                                            #
#                          Polynomials over p-adic fields                                    #
#                                                                                            #
##############################################################################################

# Getting coefficients of a flint polynomial is not intuitive.
# Flint system crashes if coefficients are not integers.
# Flint system crashes if leading coefficients are divisible by p.

# Lift termwise to a polynomial over the Flintegers.
import Hecke.lift
function lift(f :: Hecke.Generic.Poly{padic})
    R,_ = PolynomialRing(FlintZZ)
    return R([lift(c) for c in f.coeffs])
end

# This function is...kind of a hack.
# It is also very buggy since FLINT can only handle a specific case
# (integer polynomial, non-vanishing leading coefficient mod p)
function factor(f :: Hecke.Generic.Poly{padic})
    QpX = f.parent
    Qp = QpX.base_ring
    N = Qp.prec_max
    
    f_int = lift(f)
    H = factor_mod_pk_init(fint,Qp.p)
    D = factor_mod_pk(H,N)

    return Dict( QpX(lift(k))=>D[k] for k in keys(D))   
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
    q::Array{Int64,1}
end

function padic_qr(A::Hecke.Generic.MatElem{padic}; col_pivot=Val{false})
    n = size(A,1)
    L= identity_matrix(A.base_ring,n)
    U= deepcopy(A)
    P= Array(1:n)
    Pcol=Array(1:n)
    #identity_matrix(A.base_ring,n)
    
    for k=1:min(size(A,1),size(A,2))

        if col_pivot==Val{true}

            norm_list = abs.((U.entries)[k:n,k:size(A,2)])
            maxn, m = findmax( norm_list );
            
            m=m[2]+k-1;
            if m!=k
                # interchange columns m and k in U
                temp=U[:,k];
                U[:,k]=U[:,m];
                U[:,m]=temp;
                
                # interchange rows m and k in P
                temp=Pcol[k];
                Pcol[k]=Pcol[m];
                Pcol[m]=temp;
            end
        end
        
        norm_list = abs.((U.entries)[k:n,k])
        maxn, m = findmax( norm_list );
        if iszero(maxn) continue end
            
        m=m+k-1;
        if m!=k
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
            # A totally unneccesary prec. loss can be avoided here.
            # Conducting a division gives an exact result in one matrix, but reduces precision in
            # the other. Instead, we attempt to avoid the operation unless it is well-conditioned.
            #L[j,k]=U[j,k]*inv(U[k,k]);
            #U[j,:]=U[j,:]-L[j,k]*U[k,:];

            L[j,k]= _precision_stable_division(U[j,k], U[k,k])
            U[j,:]=U[j,:]-L[j,k]*U[k,:];
        end
    end
    return QRPadicPivoted(L,U,P,Pcol)
end

# Assumes that |a| ≤ |b| ≠ 0. Computes a padic integer x such that |a - xb| ≤ p^N, where N is the ring precision.
# TODO: investigate the precision.
#
function _precision_stable_division(a::padic, b::padic)
    Qp = parent(b)
    #if iszero(b) error("DivideError: integer division error") end
    if iszero(a) return zero(Qp) end
    
    x = Qp(a.u) * inv(Qp(b.u))
    x.v = a.v - b.v
    # x.N = something...
    return x
end

# stable version of nullspace for padic matrices.
function rank(A::Hecke.MatElem{padic})
    n = rows(A)
    m = cols(A)
    F = padic_qr(A)

    rank=0
    for i=1:min(n,m)
        if !iszero(F.R[i,i])
            rank += 1
        end
    end
    return rank
end

# stable version of nullspace for padic matrices.
# TODO: pivoting strategy in padic_qr does not provide the correct guarantees for this
# algorithm. I should implement a full-pivoted version.
import Hecke.nullspace
function nullspace(A::Hecke.MatElem{padic})

    m = rows(A)
    n = cols(A)
    F = padic_qr(transpose(A), col_pivot=Val{true})

    col_list = Array{Int64,1}()
    for i=1:min(n,m)
        if iszero(F.R[i,:])
            push!(col_list, i)
        end
    end

    Pinv = fill(0,length(F.p))
    for i=1:length(F.p)
        Pinv[F.p[i]] = i
    end
    
    Q = F.Q
    inv_unit_lower_triangular!(Q)
    Qinvt = transpose(Q)[Pinv,:]
    
    return length(col_list) + max(0,n-m), hcat(Qinvt[:, col_list], Qinvt[:,(m+1):n])
end

function inv_unit_lower_triangular!(L)
    m = size(L,1)
    n = size(L,2)
    if n != m
        error("Square matrix required for inverse")
    end

    for k = 1:n
        for i = k+1:n
            L[i, k] = -(L[i:i, k:i-1] * L[k:i-1, k:k])[1,1]
        end
    end
    return
end

function inv_unit_lower_triangular(L)
    L2 = deepcopy(L)
    inv_unit_lower_triangular!(L2)
    return L2
end

# Needs to be more robust. Also applied to the situation A is square but not of rank 1.
#
# a slightly generalized version of solve
# If A,b have different precisions, some strange things happen.
# TODO: honestly, just call this solve.
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
                    println(b)
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
            println(b)
            println()
            println(F.R)
            error("The system is inconsistent.")
        elseif !iszero(F.R[i,i])
            b[i,:] *= inv(F.R[i,i])
        end
    end

    return b
end


# Solve for an eigenvector using inverse iteration.
# Note that the algorithm will not converge to a particular vector in general, but the norm of
#
# A*w - λ*w converges to zero. Here, λ is the unique eigenvalue closest to `shift`, (if it is unique).
#
# TODO: I should really optimize and stabilize this later
function inverse_iteration!(A,shift,v)
    In = identity_matrix(A.base_ring, size(A,1))
    B = A - shift*In
    
    if rank(B) < cols(B)
        println("Value `shift` is exact eigenvalue.")
        return nullspace(B)[2]
    end

    function normalize(v)
        maxn, m = findmax( abs.(v.entries) )
        if iszero(maxn)
            return v
        end
        return v / v[m]
    end
    
    pow = rectangular_solve(B,identity_matrix(B.base_ring,size(B,1)))

    println(pow)
    println("---")
    println()
    
    for i=1:10
        v = normalize(pow*v)
        println(v)
        println()
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

function eigspaces(A::Hecke.Generic.Mat{T} where T <: padic)
    
    if size(A)[1] != size(A)[2]
        error("Input matrix must be square.")
    end
    
    Qp = A.base_ring
    
    # First, make everything in A a p-adic integer
    vals_of_A = valuation.( A.entries )
    min_val = minimum(vals_of_A)

    if min_val==Inf
        # In this case, A is the zero matrix.
        return identity_matrix(Qp, size(A)[1])
    end

    scale_factor = Qp(Qp.p)^min(0,Int64(-min_val))
    Aint = scale_factor * A
    
    # Solve the problem modulo p
    Amp = broadcast(modp, Aint)
    E = eigspaces(Amp)

    values_lift = fill(zero(Qp), length(E.values))
    spaces_lift = fill(zero(parent(A)), length(E.values))
    for i in 1:length(E.values)

        w = inverse_iteration(A, Qp(lift( E.values[i])), matrix(Qp, lift(E.spaces[i])))
        try
            boo, nu = iseigenvector(A,w[:,1])
            values_lift[i] = nu
            spaces_lift[i] = w

        catch
            error("Failure of convergence in inverse iteration. Likely a stability issue.")
        end
    end
    
    return EigenSpaceDec(values_lift, spaces_lift)
end



import LinearAlgebra.eigvecs
function eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)
    E = eigspaces(A)
    return hcat(E.spaces...)
end


# function eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)

#     if size(A)[1] != size(A)[2]
#         error("Input matrix must be square.")
#     end
    
#     Qp = A.base_ring
    
#     # First, make everything in A a p-adic integer
#     vals_of_A = valuation.( A.entries )
#     min_val = minimum(vals_of_A)

#     if min_val==Inf
#         # In this case, A is the zero matrix.
#         return identity_matrix(Qp, size(A)[1])
#     end

#     scale_factor = Qp(Qp.p)^Int64(min(0,-min_val))
#     Aint = scale_factor * A
    
#     # Solve the problem modulo p
#     Amp = broadcast(modp, Aint)
#     E = eigen(Amp)
#     eig_pairs = [ (E.values[i], E.vectors[:,i]) for i in 1:size(E.values)[1]]
    
#     println("Assuming that the roots of the characteristic polynomial modulo p are all distinct")

#     return  hcat([ inverse_iteration(A, Qp(lift(e)), matrix(Qp,lift(v))) for (e,v) in eig_pairs]...)
# end


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

    if iszero(A*v - (A*v)[i,1]/v[i,1]*v)
        return true,e
    else
        println(A*v - (A*v)[i,1]/v[i,1]*v)
        return false
    end
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
