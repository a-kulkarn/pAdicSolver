
using Hecke
include("../src/HeckeExt.jl")
include("../src/AlgebraicSolvers.jl")

# module aliasing
HKE = Main.HeckeExt
AS = Main.AlgebraicSolvers

Qp = PadicField(7,6)
function Base.zero(X::Type{padic})
    return Qp(0)
end


X = AS.@Ring x1 x2
n = length(X)

d = 2
M = AS.monomials(X,0:d)
s = length(M)

#P = [ Qp(1)*x1^2 + Qp(1), Qp(1)*x2^2- Qp(2)*2]
P = hcat( [ [ 2*rand(1:10)- 1 for i in 1:n] for j in 1:s]... )*M

matlist, F, B, N, Nr = AS.solve_macaulay(P,X)
#println("-- sol ", Xi)
#Er = rel_error(P,Xi,X)
#println("-- Rel error: ", norm(Er,Inf));

