
include("HeckeExt.jl")

using Hecke
using Main.HeckeExt

Qp = PadicField(7,6)
B = matrix(Qp, hcat([[rand(Qp) for i in 1:1775] for j in 1:8800]...))

@time qr = padic_qr(B);

println("Finished.")

