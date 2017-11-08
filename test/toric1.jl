using AlgebraicSolvers

X = @Ring x1 x2
n = length(X)

A1 = [one(X),x1,x2,x1*x2,x1^2,x1^2*x2]
A2 = [one(X),x1,x1^2]

p1= rand(length(A1))'*A1
p2= rand(length(A2))'*A2
P= [p1,p2]

Xi = solve_toric(P,X)

#Err = norm_eval(P,Xi,X)
#println("-- Eval err: ", norm(Err,Inf), "   ",time()-t0, "(s)"); t0 = time()

Xi
# R, L = toric_mat(P, map(p->support(p),P))

# N = nullspace(R)

# B = [one(X), x1]

# M = matmult(B,X,N,idx(L))
# Xi, Y, Z = eigdiag(M)

