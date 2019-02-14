

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

# typesafe version
function float64_valuation(x::padic)
    if iszero(x)
        return Inf
    end
    return Float64(x.v)
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
    A = matrix(Qp, fill(zero(Qp),4,4))
    for i=1:4
        for j=1:4
            A[i,j] = rand(Qp)
        end
    end
    return A
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

    # Set constants
    n = size(A,1)
    m = size(A,2)
    basezero = zero(A.base_ring)
    
    L= identity_matrix(A.base_ring,n)
    Lent = L.entries  
    Umat= deepcopy(A)
    U = Umat.entries
    
    P= Array(1:n)
    Pcol=Array(1:size(A,2))
    #identity_matrix(A.base_ring,n)

    # We cache the maximum value of the matrix at each step, so we save an iteration pass
    # through the matrix.
    norm_list = abs.(U[1:n,1:m])
    max_val, max_val_index = findmax( norm_list );

    
    # Allocate specific working memory for multiplications.
    container_for_swap = padic(U[1,1].N)
    container_for_product = padic(U[1,1].N)

    # Allocate a 2-element array to hold the index of the maximum valuation.
    max_val_index_mut = [x for x in max_val_index.I]

    # Allocate a function to zero consecutive entries of a column
    function zero_subdiagonal_of_column!(U,k::Int64)
        for j = k+1:n
            zero!(U[j,k])
        end
    end

    # The index of the diagonal point is (k,k)
    function swap_prefix_of_row!(Lent,k::Int64,i::Int64)
        for r=1:k-1
            container_for_swap = Lent[k,r]
            Lent[k,r] = Lent[i,r] 
            Lent[i,r] = container_for_swap
        end
        return
    end
    
    for k=1:(min(n,m)::Int64)

        if col_pivot==Val{true}
            #norm_list = abs.(U[k:n,k:size(A,2)])
            #maxn, m = findmax( norm_list );
            #println("max found")
            
            col_index=max_val_index_mut[2]+k-1;
            if col_index!=k
                # interchange columns m and k in U
                temp=U[:,k];
                U[:,k]=U[:,col_index];
                U[:,col_index]=temp;
                
                # interchange entries m and k in Pcol
                temp=Pcol[k];
                Pcol[k]=Pcol[col_index];
                Pcol[col_index]=temp;
            end
        end
        
        norm_list = abs.(U[k:n,k])
        maxn, row_pivot_index = findmax( norm_list );
        if iszero(maxn) continue end

        row_pivot_index=row_pivot_index+k-1;
        if row_pivot_index!=k
            # interchange rows `row_pivot_index` and `k` in U
            #temp=U[k,:];
            #U[k,:]=U[row_pivot_index,:];
            #U[row_pivot_index,:]=temp;
            for r=1:m
                U[k,r], U[row_pivot_index,r] = U[row_pivot_index,r], U[k,r]
            end               
            
            # interchange entries `row_pivot_index` and k in P
            P[k],P[row_pivot_index] = P[row_pivot_index],P[k]
            #temp=P[k];
            #P[k]=P[row_pivot_index];
            #P[row_pivot_index]=temp;

            # swap columns corresponding to the row operations already done.
            # NOTE: there is nothing to do if k<2.
            #for r=1:k-1
            #    container_for_swap = Lent[k,r]
            #    Lent[k,r] = Lent[row_pivot_index,r] 
            #    Lent[row_pivot_index,r] = Lent[k,r]
            #end
            swap_prefix_of_row!(Lent, k, row_pivot_index)
        end
        
        max_val = Inf

        # Note to self: with Julia, the optimal thing to do is split up the row operations and write a for loop.
        # The entries left of the k-th column are zero, so skip these.
        ## Cache the values of L[j,k] first.
        #     
        if iszero(U[k,k]) continue end # We have already pivoted so that abs(U[k,k]) is maximal.
        for j=k+1:n
            _unsafe_precision_stable_division!(L[j,k],U[j,k], U[k,k])
        end

        #j=k+1
        #while j<= size(A,1)
        #    zero!(U[j,k])
        #    j += 1
        #end
        zero_subdiagonal_of_column!(U,k)
        
        for r=k+1:m
            for j=k+1:n
                # Compute U[j,r] = U[j,r] - L[j,k]*U[k,r]                
                Hecke.mul!(container_for_product, L[j,k], U[k,r])
                _unsafe_minus!(U[j,r], container_for_product)
                
                # Update the biggest valuaton element
                if float64_valuation(U[j,r]) < max_val
                    max_val = float64_valuation(U[j,r])
                    max_val_index_mut[1] = j
                    max_val_index_mut[2] = r
                end
            end
        end
        #println("row operations performed: ",k)
    end
    return QRPadicPivoted(L,Umat,P,Pcol)
end


# Performs subtraction in-place, x-> x-y 
function _unsafe_minus!(x::padic, y::padic)
    x.N = min(x.N, y.N)
    ccall((:padic_sub, :libflint), Nothing,
          (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
          x, x, y, parent(x))
    return
end

# Performs multiplication and stores the result in a preexisting container
@inline function _unsafe_mult!(container::padic, x::padic, y::padic)
   container.N = min(x.N + y.v, y.N + x.v)
   ccall((:padic_mul, :libflint), Nothing,
         (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
               container, x, y, parent(x))
   return
end

# Assumes that |a| ≤ |b| ≠ 0. Computes a padic integer x such that |a - xb| ≤ p^N, where N is the ring precision.
# This prevents the division of small numbers by powers of p.
## Somehow, this is slower than the other function...
#=
function _unsafe_precision_stable_division!(container::padic, a::padic, b::padic)

    if iszero(a) return a end

    Hecke.mul!(container, a, inv(b))
    container.v = container.v - container.v
    container.N = container.N #This is wrong! fix after seminar.
    
    return
end
=#

# Try again, hope for more speed!
function _unsafe_precision_stable_division!(container::padic, a::padic, b::padic)

    if iszero(a) return a end
    # Because the division is guarenteed to be stable, manually set the precsion.
    container.N = min(a.N, b.N)
    ctx = container.parent

    ccall((:padic_div, :libflint), Cint,
          (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
          container, a, b, ctx)    
    return
end


# function divexact(a::padic, b::padic)
#    iszero(b) && throw(DivideError())
#    check_parent(a, b)
#    ctx = parent(a)
#    z = padic(min(a.N - b.v, b.N - 2*b.v + a.v))
#    z.parent = ctx
#    ccall((:padic_div, :libflint), Cint,
#          (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
#                z, a, b, ctx)
#    return z
# end


# TODO: investigate the precision.
#
# function _precision_stable_division(a::padic, b::padic)
#     Qp = parent(b)
#     #if iszero(b) error("DivideError: integer division error") end
#     if iszero(a) return zero(Qp) end
    
#     x = Qp(a.u) * inv(Qp(b.u))
#     x.v = a.v - b.v
#     # x.N = something...
#     return x
# end

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
const TESTFLAG=false
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

    if TESTFLAG
        println("---pow---")
        println(pow)
        println("---")
        println()
    end
    
    for i=1:10
        v = normalize(pow*v)
        if TESTFLAG
            println(v)
            println()
        end
    end

    if TESTFLAG
        println("---end inv iteration---")
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

    scale_factor = Qp(Qp.p)^max(0,Int64(-min_val))
    Aint = scale_factor * A
    
    # Solve the problem modulo p
    Amp = modp.(Aint)
    E = eigspaces(Amp)

    values_lift = fill(zero(Qp), length(E.values))
    spaces_lift = fill(zero(parent(A)), length(E.values))
    for i in 1:length(E.values)

        w = inverse_iteration(A, Qp(lift( E.values[i])), matrix(Qp, lift(E.spaces[i])))
        
        boo, nu = iseigenvector(A,w[:,1])
        if !boo || typeof(nu) == String
            println("-------error data-------")
            println(nu)            
            error("Failure of convergence in inverse iteration. Likely a stability issue.")
        end

        values_lift[i] = nu
        spaces_lift[i] = w

    end
    
    return EigenSpaceDec(Qp, values_lift, spaces_lift)
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
Computes if v is an eigenvector of A. If so, returns the eigenvalue as well. If not, return the error.

TODO: This function needs some work. Also the wacky return structure should be changed.
"""
function iseigenvector(A,v)
    i=1
    while i<=size(v,1)
        if !iszero(v[i,1])
            break
        end
        i+=1
    end
    if i>size(v,1)
        return false, "zero"
    end
    e = (A*v)[i,1]/v[i,1]

    if iszero(A*v - (A*v)[i,1]/v[i,1]*v)
        return true,e
    else
        return false, A*v - (A*v)[i,1]/v[i,1]*v
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
