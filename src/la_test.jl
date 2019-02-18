
include("HeckeExt.jl")

using Hecke
using Main.HeckeExt
using Profile

Qp = PadicField(7,4)
B = matrix(Qp, hcat([[rand(Qp) for i in 1:800] for j in 1:200]...))


Profile.clear_malloc_data()

B[1,1] = 0
@time qr = padic_qr(B);

#@time N = nullspace(B);



println("Finished.")

