
include("HeckeExt.jl")

using Hecke
using Main.HeckeExt
using Profile

Qp = PadicField(7,4)
B = matrix(Qp, hcat([[rand(Qp) for i in 1:400] for j in 1:400]...))

N = eigspaces(B);
Profile.clear_malloc_data()

B[1,1] = 0
N = eigspaces(B);

#@time N = nullspace(B);



println("Finished.")

