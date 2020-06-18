
#import LinearAlgebra: nullspace
#using SparseArrays

#=
function kernel(A::Matrix)
    U,S,V = svd(A)
    r=1;
    while r<=min(size(A,1),size(A,2)) && abs(S[r])> 1.e-4
        r+=1
    end
    r-=1
    #P = F[:p]
    #Pi = fill(0, length(P))
    #for i in 1:length(P) Pi[P[i]]= i end
    V[:,r+1:end]
end

function nullspace(A::AbstractSparseMatrix)
    # F = qrfact(A')
    # R = F[:R]
    # r=1;
    # while r<=min(size(A,1),size(A,2)) && abs(R[r,r])> 1.e-4
    #     r+=1
    # end
    # r-=1
    # P = F[:prow]
    # Pi = fill(0, length(P))
    # for i in 1:length(P) Pi[P[i]]= i end
    # F[:Q][Pi,r+1:end]

    F = lufact(A')
    U = F[:U]
    r = 1
    while r<=min(size(A,1),size(A,2)) && abs(U[r,r])> 1.e-4
         r+=1
    end
    r-= 1
    L = F[:L]'

    L0= L[1:r,1:r]
    K = L[1:r,r+1:end]
    N = cat(1, - L0\K, eye(size(A,2)-r))

    P = fill(0, size(A,2))
    for i in 1:size(A,2) P[F[:p][i]]= i end
    return (F[:Rs] .* N[P,:])

end

function matrix(P::Vector, M::MonomialIdx)
    A = fill(zero(coeftype(P[1])), length(P), length(M))
    j = 1
    for p in P
        for t in p
            i = get(M, t.x, 0)
            if i != 0
                A[j,i] = t.α
            end
        end
        j+=1
    end
    A
end

function smatrix(P::Vector, M::MonomialIdx)
    I = Int64[]
    J = Int64[]
    T = coeftype(P[1])
    V = T[]
    for (p,j) in zip(P,1:length(P))
        for t in p
            i = get(M, t.x, 0)
            if i != 0
                push!(I,i); push!(J,j), push!(V,t.α)
            end
        end
    end
    sparse(J,I,V)
end


function mult_basis(N, L::Vector{T}, X) where T
    Idx = idx(L)
    L0 = T[]
    for m in L
        I = map(t->get(Idx,t,0), [v*m for v in X])
        if minimum(I)!=0
            push!(L0,m)
        end
    end
    N0 = fill(zero(N[1,1]), size(N,2),length(L0))
    for i in 1:length(L0)
        for j in 1:size(N,2)
            N0[j,i]= N[get(Idx,L0[i],0),j]
        end
    end
    N0
    F = qrfact(N0, Val{true})
    B = []
    for i in 1:size(N,2)
        push!(B,L0[F[:p][i]])
    end
    B
end

=#


##################################################################
## AVI: this looks like the part we care about.
## macaulay_solve runs, though it looks like excluding the other
## things will kill the toric functionality.


# function matrix(P::Vector, M::MonomialIdx)
#     A = fill(zero(coeftype(P[1])), length(P), length(M))
#     j = 1
#     for p in P
#         for t in p
#             i = get(M, t.x, 0)
#             if i != 0
#                 A[j,i] = t.α
#             end
#         end
#         j+=1
#     end
#     A
# end

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
    
"""
INPUTS:
 B -- basis for QR factorization
 X -- variables for the polynomial ring
 K -- the "R" from the QR factorization. More precisely, the transpose of R is a submatrix of K.
 L -- Monomials of polynomial system
 ish- the "is_homogeneous" boolian.
OUTPUT:
a list of matrices whose eigenvalues are the solution coordinates.
The matrices represent multiplication-by-xi maps  ** in Q-coordinates **
"""
function mult_matrices(B, X, K, L, ish = false)

    # For an affine system, '1' is needed as a monomial as well.
    if !ish Y = vcat([parent(X[1])(1)],X) else Y = X end
    
    function construct_monomial_mult_matrix(v)
        M = fill( eltype(K)(0), length(B), size(K,2) )
        
        for (m,i) in collect(B)
            k = L[m*v]
            M[i,:] = K.entries[k,:] # Might be nice to make this compatible with matrix type.
        end
        
        return M
    end
    
    return [construct_monomial_mult_matrix(v) for v in Y]
end


# Function to compute the eigenvalues of a list of (commuting) matrices, normalized
# by the eigenvalues of the first matrix.
#
# INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
# Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M,
#          normalized by eigenvalues of the first matrix.
#          NOTE: The first matrix always represents multiplication-by-"1" or multiplication-by-x0.
#          for the reason of precision, we choose different normalizations in the affine or
#          projective cases.
function normalized_simultaneous_eigenvalues(M :: Array{Array{T,2},1} where T <: Number, ish::Bool)
    M0 = sum(A*rand() for A in M)
    #t0=time()
    I0 = inv(M0)
    #println("... inv   ", time()-t0, "(s)"); t0=time()
    Mg = I0*M[1]

    E  = eigvecs(Mg)
    #println("... eig   ", time()-t0, "(s)"); t0=time()
    Z  = E\I0

    #t0 = time()
    #F = schurfact(Mg)
    #println("... schur ", time()-t0, "(s)"); t0=time()
    # E = F[:vectors]
    # Z = E'

    X = fill(Complex{Float64}(0.0),size(M0,1),length(M))
    for j in 1:length(M)
        Yj = Z*M[j]*E
        # D = Y\Yj
        for i in 1:size(M0,1)
            X[i,j]= Yj[i,i] #(Y[:,i]\Yj[:,i])[1] #D[i,i]
        end
    end

    function normalize_solution!(Xi, ish)
        Sol = Xi
        if (!ish)
            for i in 1:size(Sol,1) Sol[i,:]/=Sol[i,1] end
        else
            for i in 1:size(Sol,1) Sol[i,:]/=norm(Sol[i,:]) end
        end
        return Sol
    end

    return normalize_solution!(X, ish)
end



