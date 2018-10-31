export support, toric_mat, solve_toric

function support(p::Polynomial)
    X = variables(p)
    r = zero(p)
    for t in p
        r+= Monomial(X,t.x.z)
    end
    r
end

function toric_mat(P, A)
    M = []
    for i in 1:length(P)
        m = one(P[i])
        for j in 1:length(A)
            if j!= i
                m*=A[j]
            end
        end
        for t in m
            push!(M,P[i]*t.x)
        end
    end
    m*= A[length(A)]
    L = [t.x for t in m]
    R = matrix(M,idx(L))
    R, L
end


function solve_toric(P, X)
    t0 = time()
    A = [support(p) for p in P]
    R, L = toric_mat(P, A)
    println("-- Toric matrix ", size(R,1),"x",size(R,2),  "   ",time()-t0, "(s)"); t0 = time()
    N = nullspace(R)
    println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    B = mult_basis(N, L, X)
    println("-- Basis ", B, "  ", time()-t0, "(s)"); t0 = time()

    M = mult_matrix(B, X, N, idx(L), false)
    println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()

    Xi = eigdiag(M)
    println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()

    for i in 1:size(Xi,1) Xi[i,:]/=Xi[i,1] end
    Xi = Xi[:,2:size(Xi,2)]
    Xi
end
