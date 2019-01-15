
using Hecke
include("../src/HeckeExt.jl")
include("../src/AlgebraicSolvers.jl")
using Main.HeckeExt
using Main.AlgebraicSolvers

Qp = PadicField(7,6)

function Base.zero(X::Type{padic})
    return Qp(0)
end


X = @Ring x1 x2
n = length(X)

d = 2
M = monomials(X,0:d)
s = length(M)

#P = [ Qp(1)*x1^2 + Qp(1), Qp(1)*x2^2- Qp(2)*2]
P = hcat( [ [ 2*rand(Qp)- 1 for i in 1:n] for j in 1:s]... )*M

Xi = solve_macaulay(P,X)
#println("-- sol ", Xi)
#Er = rel_error(P,Xi,X)
#println("-- Rel error: ", norm(Er,Inf));

