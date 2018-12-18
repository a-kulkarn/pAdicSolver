export matrix, smatrix, mult_basis, mult_matrix, eigdiag, kernel, rel_error

import LinearAlgebra: nullspace
using SparseArrays

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

# AVI:
# INPUTS:
#  B -- basis for QR factorization
#  X -- variables for the polynomial ring
#  K -- ????
#  KM-- ????
#  ish- the "is_homogeneous" boolian.
# OUTPUT:
# ????
function mult_matrix(B, X, K, KM, ish = false)
    R = []
    Idx = idx(B)
    if !ish
        M = fill(0.0, length(B), size(K,2) )
        for (m,i) in Idx.terms
            k = get(KM, m, 0)
            if k != 0
                for j in 1:size(K,2)
                    M[i,j] = K[k,j]
                end
            end
        end
        push!(R,M)
    end
    for v in X
        M = fill(0.0, length(B), size(K,2) )
        for (m,i) in Idx.terms
            k = get(KM, m*v, 0)
            if k != 0
                for j in 1:size(K,2)
                    M[i,j] = K[k,j]
                end
            end
        end
        push!(R,M)
    end
    R
end

## 
function eigdiag(M)
    M0 = sum(M[i]*rand() for i in 1:length(M))
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
    X
end

function (p::Polynomial{B,T})(x::Vector) where {B,T}
   r = zero(x[1]);
   for m in p
      t=m.α
      for i in 1:length(m.x.z)
      	 t*=x[i]^m.x.z[i]
      end
      r+=t
   end
   r
end

##################################################################

#=
"""
Vector of relative errors of P at the points X
"""
function rel_error(p, Xi::Matrix, X = variables(p))
    r = fill(0.0, size(Xi,1), length(p))
    n = size(Xi,2)
    for i in 1: size(Xi,1)
        for j in 1:length(p)
            V = Xi[i,:]
            r[i,j]= norm(subs(p[j],X=>V))
            s = 0.0
            for t in p[j]
                s+= norm(subs(t, X => V))
            end
            r[i,j]/=s
        end
    end
    r
end
=#
