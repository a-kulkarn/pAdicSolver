export macaulay_mat, qr_basis, solve_macaulay

function is_not_homogeneous(p)
    L = [degree(t) for t in p]
    maximum(L) != minimum(L)
end
    
function macaulay_mat(P, L::AbstractVector, X, ish = false )
    d = maximum([degree(m) for m in L])
    if ish
        Q = [monomials(X,d-degree(P[i])) for i in 1:length(P)]
    else
        Q = [monomials(X,0:d-degree(P[i])) for i in 1:length(P)]
    end
    M = []
    for i in 1:length(P)
        for m in Q[i]
            push!(M,P[i]*m)
        end
    end
    matrix(M,idx(L))
end

function qr_basis(N, L, ish = false)
    Idx= idx(L)
    if ish
        L0 = filter(m->(m.z[1]>0), L)
    else
        d  = maximum([degree(m) for m in L])
        L0 = filter(m->degree(m)<d,L)
    end
    N0 = fill(zero(N[1,1]), size(N,2),length(L0))
    for i in 1:length(L0)
        for j in 1:size(N,2)
            N0[j,i]= N[get(Idx,L0[i],0),j]
        end
    end
    N0
    F= qrfact(N0, Val(true))
    B = []
    for i in 1:size(N,2)
        push!(B,L0[F[:p][i]])
    end
    if ish
        X = variables(L[1])
        for i in 1:length(B)
            B[i] = B[i]*X[1]^(-1)
        end
    end
    B
end
    
solve_macaulay = function(P, X, rho =  sum(degree(P[i])-1 for i in 1:length(P)) + 1 )
    println()
    println("-- Degrees ", map(p->degree(p),P))
    ish = !any(is_not_homogeneous, P)
    println("-- Homogeneity ",ish)
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
    B = qr_basis(N, L, ish)
    println("-- Qr basis ",  length(B), "   ",time()-t0, "(s)"); t0 = time()
    M = mult_matrix(B, X, N, idx(L), ish)
    println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()
    Xi = eigdiag(M)
    println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()
    #println("-- Rel. error ", eigdiag_error(M,Xi,Y,Z))
    if (!ish)
        Xi = diagm([1/Xi[i,1] for i in 1:length(B)])*Xi
        Xi = Xi[:,2:size(Xi,2)]
    else
        Xi = diagm([1/norm(Xi[i,:]) for i in 1:length(B)])*Xi
    end
    Xi
end
