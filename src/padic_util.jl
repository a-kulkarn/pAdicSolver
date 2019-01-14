###################################
## Some new p-adic linear algebra

# Compute a factorization of a padic matrix A = QR, where R is
# upper triangular (Borel) and Q lies in the maximally compact
# subgroup of SL(Qp) = SL(Zp).
#
# It turns out that the algorithm to do this is just the LU factorization
# with pivoting.

# Extend Base to function padically
using Hecke
import Base./
import Base.abs
import Hecke.valuation

function /(x::padic,y::padic)
    return x//y
end

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


# Hopefully, LAPACK selects pivots to put the largest
# diagonal or subdiagonal element on the diagonal. Experimentally
# this appears to be the case.

# UPDATE: I don't think LAPACK does this...
#function padic_qr(A::Array{padic,2})
#    return lu(A,Val(true))
#end

function padic_qr(A)
    n = size(A,1);
    L= identity_matrix(A.base_ring,n);
    U=deepcopy(A);
    P=identity_matrix(A.base_ring,n);
    
    for k=1:n
        norm_list = abs.((U.entries)[k:n,k])
        _, m = findmax( norm_list );
        m=m+k-1;
        if m!=k
            # Note: we actually want inv(P), so we should compute differently
            # interchange rows m and k in U
            temp=U[k,:];
            U[k,:]=U[m,:];
            U[m,:]=temp;
            # interchange rows m and k in P
            temp=P[k,:];
            P[k,:]=P[m,:];
            P[m,:]=temp;

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
    return L,U,P
end

function modp(x::padic)
    Fp = ResidueRing(FlintZZ,parent(x).p)
    return Fp(lift(x))
end

## Test utilities
function test_rings()
    return Qp = FlintPadicField(7,20), ResidueRing(FlintZZ,7)
end

function random_test_matrix(Qp)
    return matrix(Qp, rand.(fill(1:7^20),4,4))
end
####


# Really, I should just define an iterator structure on A
function modp(A)
    Fp = ResidueRing(FlintZZ, A.base_ring.p)
    return matrix(Fp, modp.(A.entries))
end


## Make things a little more consistent with the other Julia types
function my_nullspace(A::nmod_mat)
    nu,N = nullspace(A)
    return nu, nu==0 ? matrix(A.base_ring, fill(0,size(A,2),0)) : N[:,1:nu]
end


# Extend Linear algebra
import LinearAlgebra: eigen, Eigen, eigvals, eigvecs

struct MyEigen{T}
    values::Array{T,1}
    vectors::Union{Hecke.Generic.Mat{T}, nmod_mat}
end

# Returns an eigen factorization structure like the default LinearAlgebra.eigen function.
#
# Fails when A mod p has no eigenvalues, because list has eltype Any
function eigen(A::nmod_mat)
    R,_ = PolynomialRing(A.base_ring)
    g = charpoly(R, A)
    rts = roots(g)
    if isempty(rts) then
        error("Not implemented if no roots of char. poly. over the finite field")
    end
    
    Imat = identity_matrix(A.base_ring, size(A,1))

    eig_list = [ (r,my_nullspace(A-r*Imat)...) for r in rts]

    eig_vals = vcat([fill( ep[1] , ep[2]) for ep in eig_list]...)
    eig_vecs = hcat([ ep[3] for ep in eig_list]...)    
    return MyEigen(eig_vals, eig_vecs)
end

# See usual eigvecs
function eigvecs(A::nmod_mat)
    return eigen(A).vectors    
end

# See usual eigvals
function eigvals(A::nmod_mat)
    return eigen(A).values
end

import Hecke.Generic.charpoly
# Sometimes I don't want to construct the ambient ring
function charpoly(A::nmod_mat)
    R,_ = PolynomialRing(A.base_ring)
    return charpoly(R, A)
end




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

#import StatsBase.countmap
function eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)

    Qp = A.base_ring
    # First, make everything in A a p-adic integer
    println("Skipping an important step because lazy...")

    # Solve the problem modulo p
    Amp = modp(A)
    E = eigen(Amp)
    eig_pairs = [ (E.values[i], E.vectors[:,i]) for i in 1:size(E.values)[1]]
    
    println("Assuming that the roots of the characteristic polynomial modulo p are all distinct")
   
    return  hcat([ inverse_iteration(A, Qp(lift(e)), matrix(Qp,lift(v))) for (e,v) in eig_pairs]...)
end


# function for testing
function is_eigenvector(A,v)
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
