

using Dory, Random

# Enforce that tests are identical.
Random.seed!(123)

Qp = PadicField(31,10)

R, (x1,x2) = PolynomialRing(Qp, ["x1", "x2"])
n = length(gens(R))
num_samp = 10

failed_test_count = 0
for d in [2,4,6,10, 15]

    tot_time  = 0
    tot_error = 0
    M = monomials_of_degree(gens(R),0:d)
    s = length(M)

    for i=1:num_samp
        try
            global P = hcat([[rand_padic_int(Qp) for i in 1:n] for j in 1:s]...)*M

            t0 = time()
            sol = solve_macaulay(P)
            t1 = time()
            tot_time += t1-t0
            
            println("\n-- sol ")
            println(sol,"\n")

            Er = rel_error(P,sol)
            println("-- Absolute forward error: ")
            display(Er)
            println()
        catch
            global failed_test_count += 1
        end
    end

    avg_time = tot_time/num_samp

    @info "" avg_time
end

println("----------------------")
println("Tests failed due to stability issues (probably): ", failed_test_count)
