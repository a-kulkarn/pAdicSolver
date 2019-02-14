
include("HeckeExt.jl")

using Hecke
using Main.HeckeExt
using Profile

Qp = PadicField(7,6)
B = matrix(Qp, hcat([[rand(Qp) for i in 1:700] for j in 1:400]...))

@time qr = padic_qr(B);


B[1,1] = zero(Qp);

Profile.clear_malloc_data()
@time qr = padic_qr(B)

println(iszero(B[qr.p,:] - qr.Q*qr.R))

println("Finished.")

