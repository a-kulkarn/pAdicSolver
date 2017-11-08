export der, alpha_beta, jacobian, newton_iter, newton_improve, rel_error

function der{C,T}(p::Polynomial{C,T}, v::PolyVar{C})
    r = zero(p)
    X = variables(p)
    i = findfirst(X, v)
    if i>0
        for t in p
            if t.x.z[i] >0
                e = copy(t.x.z);
                d = e[i]
                e[i] = d-1;
                r+= Monomial(X,e)*t.α*d
            end
        end
    end
    r
end

function jacobian(P,X)
    [der(P[i],X[j]) for i in 1:length(P), j in 1:length(X)]
end

function newton_iter(F, J, Xi)
    J0 = fill(zero(Xi[1]), size(J,1), size(J,2))
    for i in 1:size(J,1)
        for j in 1:size(J,2)
            if J[i,j] != zero(J[i,j])
                J0[i,j] = J[i,j](Xi)
            end
        end
    end
    F0 = fill(zero(Xi[1]), length(F))
    for i in 1:length(F)
        if F[i] != zero(F[i])
                F0[i] = F[i](Xi)
        end
    end
    err = norm(F0)
    Xi -= J0\F0
    Xi, err
end

    
function newton_improve(Xi::Matrix, P, X,  eps::Float64=1.e-12, Nit::Int64 = 20)
    J = jacobian(P,X)
    for j in 1: size(Xi,1)
        i = 1
        err = 1.0
        V = Xi[j,:]
        while err>eps && i< Nit
            V, err = newton_iter(P, J, V)
            i+=1
        end
        if i==Nit && err>eps
            println("err: ", err, "   ", V,  " nwt: no convergence ",)
        end
        Xi[j,:] = V
    end
    Xi
end


"""
alpha, beta quantities for Newton convergence to an approximate zero.

 - If alpha < 0.125, then the approximate zero is withing 2*beta from Xi and Newton methods converges to it from Xi quadratically.
 - If alpha < 0.02, then Newton method converges from all points in the ball of center Xi and radius 2*beta.

"""
function alpha_beta(P::Vector, Xi::Vector)
    X  = variables(P[1])
    for i in 2:length(P) X = union(X,variables(P[i])) end
    F0 = fill(zero(Xi[1]),length(P))
    for (p,i) in zip(P,1:length(P))
        f=subs(p,X=>Xi)
        if f != zero(f)
            F0[i]=(f[1]).α
        end
    end
    J0 = fill(zero(Xi[1]), length(P), length(X))
    for i in 1:length(P)
        for j in 1:length(X)
            dP = subs(der(P[i],X[j]),X=>Xi)
            if dP != zero(dP)
                J0[i,j] = (dP[1]).α
            end
        end
    end
    beta = norm(J0\F0)

    s  =  sqrt(1 + norm(Xi)^2)
    mu =  sqrt(sum(norm(p,degree(p))^2 for p in P))
    mu *= norm(J0\diagm([sqrt(degree(p))*s^(degree(p)-1) for p in P]))

    d  = maximum([degree(p) for p in P])    
    gamma = 0.5*d*sqrt(d)*mu

    beta*gamma, 2*beta
end
