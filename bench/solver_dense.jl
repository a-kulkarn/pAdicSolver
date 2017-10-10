function solver_dense(d,X)

    n = length(X)

    M = monomials(X,0:d)
    s = length(M)

    P = (2*rand(n,s)-fill(1.0,n,s))*M
    t0 =time()
    Xi = solve_macaulay(P,X)
    t1= time()-t0
    println("-- Number of solutions: ",size(Xi,1),"    ",t1,"(s)" )
    Er = rel_error(P,Xi)
    println("-- Rel error: ", norm(Er,Inf))
    println()
    t1
end

function init(X)
    solver_dense(2,[X[1],X[2]])
end

function save(file::String, n, dg, t)
    io=open(file,"w")
    println(io,"n=", length(X))
    println(io,"d=", dg)
    println(io,"t=", t)
    close(io)
    println("n=", length(X))
    println("d=", dg)
    println("t=", t)
end
