
using Hecke
include("../src/HeckeExt.jl")
include("../src/AlgebraicSolvers.jl")

# module aliasing
HKE = Main.HeckeExt
AS = Main.AlgebraicSolvers


function raw_sol_test(P,sol)
    return [[p(X=>sol.entries[i,2:size(sol,2)]) for i in 1:size(sol,1)] for p in P]
end

function rel_error(P,sol)
    return [[p(X=>sol.entries[i,:]) for i in 1:size(sol,1)] for p in P]
end

# For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
Qp = PadicField(89,6)
function Base.zero(X::Type{padic})
    return Qp(0)
end


X = AS.@Ring x1 x2
n = length(X)

d = 2
M = AS.monomials(X,0:d)
s = length(M)

#P = [ Qp(1)*x1^2 + Qp(1), Qp(1)*x2^2- Qp(2)*2]
P = hcat( [ [ 2*rand(Qp)- 1 for i in 1:n] for j in 1:s]... )*M

#matlist, F, B, N, Nr = AS.solve_macaulay(P,X);
sol = AS.solve_macaulay(P,X)

println("\n-- sol ")
println(sol,"\n")

Er = rel_error(P,sol)
println("-- Rel error: ")
println(Er);

