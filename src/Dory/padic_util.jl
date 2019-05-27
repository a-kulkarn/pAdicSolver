

# simple function to invert a permutation array.
function inverse_permutation(A::Array{Int64,1})
    Pinv = fill(0,length(A))
    for i=1:length(A)
        Pinv[A[i]] = i
    end
    return Pinv
end

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
    return Qp(rand(1:BigInt(p)^N))*Qp(p)^(rand(-N:N))
end

function rand_padic_int(Qp::FlintPadicField)
    p = Qp.p
    N = Qp.prec_max
    return Qp(rand(1:BigInt(p)^N))
end


function random_test_matrix(Qp,n=4)
    A = matrix(Qp, fill(zero(Qp),n,n))
    for i=1:n
        for j=1:n
            A[i,j] = rand_padic_int(Qp)
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
    H = factor_mod_pk_init(f_int,Qp.p)
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
    n = size(A,1)::Int64
    m = size(A,2)::Int64
    basezero = zero(A.base_ring)
    
    L= identity_matrix(A.base_ring,n)
    Lent = L.entries::Array{padic,2}  
    Umat= deepcopy(A)
    U = Umat.entries
    
    P= Array(1:n)
    Pcol=Array(1:m)

    # We cache the maximum value of the matrix at each step, so we save an iteration pass
    # through the matrix.
    val_list = float64_valuation.(U)
    min_val, min_val_index = findmin( val_list );
    
    # Allocate specific working memory for multiplications.
    container_for_swap    = padic(U[1,1].N)
    container_for_product = padic(U[1,1].N)
    container_for_div     = padic(U[1,1].N)
    
    # Allocate a 2-element array to hold the index of the maximum valuation.
    min_val_index_mut = [x for x in min_val_index.I]

    # Allocate a function to zero consecutive entries of a column
    function zero_subdiagonal_of_column!(U,k::Int64)
        for j = k+1:n
            zero!(U[j,k])
        end
    end
    
    
    for k=1:(min(n,m)::Int64)

        if col_pivot==Val{true}
            col_index=min_val_index_mut[2]
            if col_index!=k
                # interchange columns m and k in U
                for r=1:n
                    U[r,k], U[r,col_index] = U[r,col_index], U[r,k]
                end
                
                # interchange entries m and k in Pcol
                temp=Pcol[k];
                Pcol[k]=Pcol[col_index];
                Pcol[col_index]=temp;
            end
        end
        
        val_list = float64_valuation.(U[k:n,k])
        minn, row_pivot_index = findmin( val_list );
        if minn==Inf continue end

        row_pivot_index=row_pivot_index+k-1;
        if row_pivot_index!=k

            # interchange rows `row_pivot_index` and `k` in U
            for r=1:m
                U[k,r], U[row_pivot_index,r] = U[row_pivot_index,r], U[k,r]
            end               
            
            # interchange entries `row_pivot_index` and k in P
            P[k],P[row_pivot_index] = P[row_pivot_index],P[k]

            # swap columns corresponding to the row operations already done.
            swap_prefix_of_row!(Lent, k, row_pivot_index)
        end

        # Reset min_valuation for selecting new pivot.
        min_val = Inf

        # Note to self: with Julia, the optimal thing to do is split up the row operations
        # and write a for loop.
        # The entries left of the k-th column are zero, so skip these.
        # Cache the values of L[j,k] first.
        #
        if iszero(U[k,k]) continue end 

        # The use of the inversion command preserves relative precision. By row-pivoting,
        # the extra powers of p cancel to give the correct leading term.
        # the "lost" digits of precision for L[j,k] can simply be set to 0.
        container_for_inv = inv(U[k,k]) 
        
        for j=k+1:n
            Hecke.mul!(L[j,k],U[j,k], container_for_inv)
            L[j,k].N = parent(L[j,k]).prec_max            # L[j,k] is really an integer.
        end

        zero_subdiagonal_of_column!(U,k)
        
        for r=k+1:m
            for j=k+1:n
                # Compute U[j,r] = U[j,r] - L[j,k]*U[k,r]                
                Hecke.mul!(container_for_product, L[j,k], U[k,r])
                _unsafe_minus!(U[j,r], container_for_product)
                
                # Update the smallest valuation element
                if float64_valuation(U[j,r]) < min_val
                    min_val = float64_valuation(U[j,r])
                    min_val_index_mut[1] = j
                    min_val_index_mut[2] = r
                end
            end
        end
    end

    # Perform a last column pivot in case the corner entry is 0.
    #
    # ...
    
    @assert iszero(A[P,Pcol] - L*Umat)
    
    return QRPadicPivoted(L,Umat,P,Pcol)
end

# The index of the diagonal point is (k,k)
function swap_prefix_of_row!(Lent, k::Int64, i::Int64)
    for r=1:(k-1)
        container_for_swap = Lent[k,r]
        Lent[k,r] = Lent[i,r] 
        Lent[i,r] = container_for_swap
    end
    return
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
# @inline function _unsafe_mult!(container::padic, x::padic, y::padic)
#    container.N = min(x.N + y.v, y.N + x.v)
#    ccall((:padic_mul, :libflint), Nothing,
#          (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
#                container, x, y, parent(x))
#    return
# end

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
# function _unsafe_precision_stable_division!(container::padic, a::padic, b::padic)

#     if iszero(a) return a end
#     # Because the division is guarenteed to be stable, manually set the precsion.
#     container.N = min(a.N, b.N)
#     ctx = container.parent

#     ccall((:padic_div, :libflint), Cint,
#           (Ref{padic}, Ref{padic}, Ref{padic}, Ref{FlintPadicField}),
#           container, a, b, ctx)    
#     return
# end


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

# IMPORTANT!
# We deviate slightly from LinearAlgebra's SVD structure by putting a diagonal matrix for S.
struct SVDPadic
    U::Hecke.Generic.MatElem{padic}
    S::Hecke.Generic.MatElem{padic}
    Vt::Hecke.Generic.MatElem{padic}
    p::Array{Int64,1}
    q::Array{Int64,1}
end

# A padic analogue for svd
import LinearAlgebra.svd
function svd(A::Hecke.Generic.MatElem{padic})

    F = padic_qr(A, col_pivot=Val{true})
    G = padic_qr(transpose(F.R))

    @assert G.p == [i for i=1:length(G.p)]

    U = deepcopy(F.Q)
    S = transpose(G.R)
    Vt= transpose(G.Q)
    
    @assert iszero( A[F.p,F.q] - U*S*Vt)

    return SVDPadic(U,S,Vt,F.p,F.q)
end

# stable version of nullspace for padic matrices.
function rank(A::Hecke.MatElem{padic})
    n = nrows(A)
    m = ncols(A)
    F = padic_qr(A)

    rank=0
    for i=1:min(n,m)
        if !iszero(F.R[i,i])
            rank += 1
        end
    end
    return rank
end

# Returns the p-adic singular values of a matrix
function singular_values(A::Hecke.MatElem{padic})
    F = padic_qr(A,col_pivot=Val{true})
    return [ F.R[i,i] for i=1:minimum(size(A)) ]
end

# stable version of nullspace for padic matrices.
import Hecke.nullspace
function nullspace(A::Hecke.MatElem{padic})

    m = nrows(A)
    n = ncols(A)
    F = padic_qr(transpose(A), col_pivot=Val{true})

    col_list = Array{Int64,1}()
    for i=1:min(n,m)
        if iszero(F.R[i,:])
            push!(col_list, i)
        end
    end

    Pinv = inverse_permutation(F.p)   
    
    Q = F.Q
    inv_unit_lower_triangular!(Q)
    Qinvt = transpose(Q)[Pinv,:]
    
    return length(col_list) + max(0,n-m), hcat(Qinvt[:, col_list], Qinvt[:,(m+1):n])
end

function inv_unit_lower_triangular!(L::Hecke.Generic.MatElem{T} where T)

    m = size(L,1)::Int64
    n = size(L,2)::Int64    
    #if !issquare(L)
    #    error("Square matrix required for inverse")
    #end
    Qp = parent(L[1,1])
    container_for_mul = Qp()
    container_for_result = Qp()
    
    for k = 1:n
        for i = k+1:n
            container_for_result=zero(Qp)
            for r=k:i-1
                Hecke.mul!(container_for_mul, L[i,r], L[r,k])
                addeq!(container_for_result,  container_for_mul)
            end
            L[i,k] = -container_for_result
        end
    end

    return
end

function inv_unit_lower_triangular(L)
    L2 = deepcopy(L)
    inv_unit_lower_triangular!(L2)
    return L2
end

# A slightly generalized version of solve
# If A,b have different precisions, some strange things happen.
# TODO: honestly, just call this solve.
#
# Parameter `stable` determines whether qr or svd method is used. Default is for qr (speed).
#
function rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic}; stable::Bool=false)
    if !stable
        return _lu_rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})
    else
        return _svd_rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})
    end
end

# Specialization to lu-solve
function _lu_rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})

    m = nrows(A)
    n = ncols(A)
    if nrows(b_input) != m
        error("`A` and `b` must have the same number of rows.")
    end
    b = deepcopy(b_input)

    if m < n
        error("System is underdetermined. Use `underdetermined_solve` instead.")
    end

    F = padic_qr(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            b[i,:] = b[i,:] - b[j,:]* F.Q[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:ncols(b)
                if !iszero(b[i, j])
                    println()
                    println("--- Error data: ---")
                    println("bad entry at ", i," ",j)
                    println("entries: ", b[i,j])
                    println()
                    error("Line 461: The system is inconsistent.")
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # backward substitution
    for i in n:-1:1
        for j in (i+1):n
            b[i,:] = b[i,:] - b[j,:]*F.R[i,j]
        end
        #scale = A[i, i]
        #b.row_op(i, lambda x, _: x / scale)

        if !iszero(b[i,:]) && iszero(F.R[i,i])
            println()
            println("--- Error data: ---")
            println("bad entry at row ", i)
            error("Line 480: The system is inconsistent.")
        elseif !iszero(F.R[i,i])
            b[i,:] *= inv(F.R[i,i])
        end
    end

    return b
end

# Specialization to svd-solve
function _svd_rectangular_solve(A::Hecke.MatElem{padic}, b_input::Hecke.MatElem{padic})

    m = nrows(A)
    n = ncols(A)
    if nrows(b_input) != m
        error("`A` and `b` must have the same number of rows.")
    end
    b = deepcopy(b_input)

    if m < n
        error("System is underdetermined. Use `underdetermined_solve` instead.")
    end

    F = svd(A)
    b = b[F.p,:]

    # forward substitution, all diag entries are scaled to 1
    for i in 1:m
        for j in 1:(i-1)
            b[i,:] = b[i,:] - b[j,:]* F.U[i,j]
        end
    end

    # consistency check for overdetermined systems
    if m > n
        for i in (n+1):m
            for j in 1:ncols(b)
                if !iszero(b[i, j])
                    println()
                    println("--- Error data: ---")
                    println("bad entry at ", i," ",j)
                    println("entries: ", b[i,j])
                    # println()
                    # println(b)
                    # println()
                    # println(A)
                    error("Line 533: The system is inconsistent.")
                end
            end
        end
    end
    b = b[1:n, :]   # truncate zero rows if consistent

    # Scaling step
    for i in 1:n
        if !iszero(b[i,:]) && iszero(F.S[i,i])
            println()
            println("--- Error data: ---")
            error("The system is inconsistent: singular value: ", i," is zero, while `b[i,:]` is nonzero.")
        elseif !iszero(F.S[i,i])
            b[i,:] *= inv(F.S[i,i])
        end
    end

    # backward substitution
    for i in n:-1:1
        for j in (i+1):n
            b[i,:] = b[i,:] - b[j,:]*F.Vt[i,j]
        end
    end

    return b[F.q,:]
end
    
function underdetermined_solve()
    error("Not implemented.")
    return
end


#************************************************************
#
# Eigenvector iteration methods.
#
#************************************************************


#************************************************
#  Basic inverse iteration
#************************************************

# Solve for an eigenvector using inverse iteration.
# Note that the algorithm will not converge to a particular vector in general, but the norm of
#
# A*w - λ*w converges to zero. Here, λ is the unique eigenvalue closest to `shift`, (if it is unique).
#
# TODO: Separate invariant subspaces at high valuation.
const TESTFLAG=false
function inverse_iteration!(A,shift,V)

    # Note: If A is not known to precision at least one, really bad things happen.
    Qp = A.base_ring
    In = identity_matrix(A.base_ring, size(A,1))
    B  = A - shift*In
    
    if rank(B) < ncols(B)
        println("Value `shift` is exact eigenvalue. `shift` = ", shift)
        return [nullspace(B)[2]], [shift]
    end

    function normalize(V)
        maxn, m = findmax( abs.(V.entries) )
        if iszero(maxn)
            return V
        end
        return V / V[m]
    end
    
    @time pow = rectangular_solve(B,identity_matrix(B.base_ring,size(B,1)),stable=true)

    if TESTFLAG
        println("---pow---")
        println(pow)
        println("---")
        println()
    end
    
    @time for i=1:(A.base_ring.prec_max)
        V = normalize(pow*V)
        if TESTFLAG
            println(V)
            println()
        end
    end
    
    if TESTFLAG
        println("---end inv iteration---")
        println()
    end

    # Test for convergence and calculate eigenvalues,
    # if the algorithm hasn't converged, check for whether the remaining subspace can be
    # split by further iteration.
    X = try
        rectangular_solve(V, A*V, stable=true)        
    catch e
        error("Error in inverse iteration. Likely a stability issue.")
    end

    nu= trace(X)//size(X,2)
    Y = X - nu*identity_matrix(Qp,size(X,2))
    
    if iszero(Y)
        # In this case, the eigenvectors are at their maximum refinement.
        return [V],[nu]
    end

    # Since only eigenvectors (mod p) are given as the initial data, the operator Y *must* be
    # zero mod p. We scale out the denominator to try again.
    
    vals_of_Y = valuation.( Y )
    min_val = minimum(vals_of_Y)

    if min_val <=0
        error("Failure of convergence in inverse iteration.")
    end

    println("Second level iteration.")
    scale_factor = Qp(Qp.p)^Int64(-min_val)
    inv_scale_factor = Qp(Qp.p)^Int64(min_val)
    Ynew = scale_factor * Y    
    E = eigspaces(Ynew)

    return Array{typeof(V),1}([V*Esp for Esp in E.spaces]),
    Array{padic,1}([inv_scale_factor*nu for nu in E.values])
end

# Given an approximate subspace to an invariant subspace, return the
# invariant subspace and the eigenvalue.
#
# Raises an error if there is a failure of convergence.
#
function inverse_iteration(A, shift, v)
    w = deepcopy(v)
    wlist,nulist = inverse_iteration!(A,shift,w)    
    return wlist,nulist
end

function inverse_iteration_decomposition(A, Amp)

    Qp = A.base_ring
    E = eigspaces(Amp)

    values_lift = fill(zero(Qp), 0)
    spaces_lift = fill(zero(parent(A)), 0)

    for i in 1:length(E.values)

        # Approximate input data
        appx_eval = Qp( lift(E.values[i]) )
        appx_espace =  matrix(Qp, lift(E.spaces[i]) )

        # Apply inverse iteration step.
        wlist,nulist = inverse_iteration(A, appx_eval, appx_espace)

        # Append refined data to the main list.
        values_lift = vcat(values_lift, nulist)
        spaces_lift = vcat(spaces_lift,  wlist)
    end

    return values_lift, spaces_lift
end

#************************************************
#  Power iteration with recovery
#************************************************

# function power_iteration_decomposition(A, chi)
#     if is_diagonal(A)
#         return 
#     end
# end


###############################################################################
#
#   Hessenberg form
#
###############################################################################

function hessenberg!(A::Hecke.Generic.Mat{T} where T <: padic)
    !issquare(A) && error("Dimensions don't match in hessenberg")
    R = base_ring(A)
    n = nrows(A)
    u = R()
    t = R()
    
    for m = 1:n - 2

        val_list = float64_valuation.(A.entries[ m+1 : n , m ])
        minn, row_pivot_index = findmin( val_list );
        if minn==Inf continue end

        i = row_pivot_index + m;
        
        # Perform a row/column swap to move the pivot to the subdiagonal
        if i > m+1
            for j = m:n
                A[i, j], A[m+1, j] = A[m+1, j], A[i, j]
            end
            for j = 1:n
                A[j, i], A[j, m+1] = A[j, m+1], A[j, i]
            end
        end

        # cache the inverted pivot.
        h = -inv(A[m+1, m])

        # Perform the elimination.
        for i = m + 2:n
            if iszero(A[i, m]) continue end
            
            u = Hecke.mul!(u, A[i, m], h)

            # Row operatons
            for j = m+1:n
                t = Hecke.mul!(t, u, A[m+1, j])
                A[i, j] = addeq!(A[i, j], t)
            end
            u = -u

            # Column eliminations
            for j = 1:n
                t = Hecke.mul!(t, u, A[j, i])
                A[j, m+1] = addeq!(A[j, m+1], t)
            end
            A[i, m] = R()            
        end        
    end
end

"""
    hessenberg(A::Generic.MatrixElem{T}) where {T <: padic}
> Returns the Hessenberg form of M, i.e. an upper Hessenberg matrix
> which is similar to M. The upper Hessenberg form has nonzero entries
> above and on the diagonal and in the diagonal line immediately below the
> diagonal.
> A padically stable form of the algorithm is used, where pivots are 
> selected carefully.
"""
function hessenberg(A::Hecke.Generic.Mat{T} where T <: padic)
   !issquare(A) && error("Dimensions don't match in hessenberg")
   M = deepcopy(A)
   hessenberg!(M)
   return M
end


#************************************************
#  QR-iteration 
#************************************************

"""
   blockschurform

    Computes the block schur form of a padic matrix A, where the
    blocks correspond to the different eigenvalues of A modulo p.
"""
function block_schur_form(A::Hecke.Generic.Mat{T} where T <: padic)

    Qp = A.base_ring
    N = Qp.prec_max
    
    # Extract data from the reduction modulo p
    Aint  = _normalize_matrix(A)
    Amp   = modp.(Aint)
    chiAp = charpoly(Amp)

    B = hessenberg(A)
    id= identity_matrix(Qp, size(B)[1])
    
    for (rt,m) in roots_with_multiplicities(chiAp)
        
        lambdaI = lift(rt)*id

        # Regarding convergence. It seems like it needs a little extra time to
        # sort the terms via permutation.
        for i in 1:N*m
            F = padic_qr(B - lambdaI)

            # Note about Julia's syntax. A[:,F.p] = A*inv(P), for a permutation P.
            B = F.R[:,F.p] * F.Q + lambdaI
        end        
    end
    
    return B
end

function roots_with_multiplicities(f)
    F = Hecke.factor(f)
    return [(-g(0), m) for (g,m) in F if Hecke.degree(g) == 1]
end


function qr_iteration_decomposition(A,Amp)
    
end



function _normalize_matrix(A)

    Qp = A.base_ring
    vals_of_A = valuation.( A.entries )
    min_val = minimum(vals_of_A)

    scale_factor = Qp(Qp.p)^max(0,Int64(-min_val))
    return scale_factor * A
end

"""
eigvecs(A::Hecke.Generic.Mat{T} where T <: padic)

Compute the eigenvectors of a padic matrix iteratively.

The `method` parameter selects the method to be used to compute the eigenvectors.
The intended options are:

-- "inverse"
-- "classical"
-- "power"
-- "qr"

The default is "inverse", since at the moment this is the one that is implemented.

"""

function eigspaces(A::Hecke.Generic.Mat{T} where T <: padic; method="inverse")

    ## Input sanitization    
    if size(A)[1] != size(A)[2]
        error("Input matrix must be square.")
    end

    ## Set constants
    Qp = A.base_ring

    ##### Begin main computations #####
    
    if iszero(A)        
        return EigenSpaceDec(Qp, [zero(Qp)] , [identity_matrix(Qp, size(A)[1])] )
    end

    if method == "classical"
        error("Not Implemented")
    end
        
    # Extract data from the reduction modulo p
    Aint  = _normalize_matrix(A)
    Amp   = modp.(Aint)
    chiAp = charpoly(Amp)
    factors_chiAp = Hecke.factor(chiAp)
    
    if isirreducible(chiAp)
        empty_array = Array{padic,1}()
        return EigenSpaceDec(Qp, empty_array , [matrix(Qp, size(A)[1], 0, empty_array)] )
    end

    
    # FAILSAFE DURING DEVELOPMENT...
    # Fail automatically if there are large invariant subspaces mod p
    if any( e >= 2 for (f,e) in factors_chiAp if degree(f)==1 )
        error("Not implemented when roots are not squarefree") 
    end

    
    ## Call decomposition method
    if method == "inverse"
        values_lift, spaces_lift = inverse_iteration_decomposition(A, Amp)
    else
        error("Not Implemented")
    end

    ## Post-processing
    
    return EigenSpaceDec(Qp, values_lift, spaces_lift)
end



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

