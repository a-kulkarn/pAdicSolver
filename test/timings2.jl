

using Dory, Random

# Enforce that tests are identical.
Random.seed!(123)

Qp = PadicField(31,10)



failed_test_count = 0
for n in [2,4,5]

    R, vars = PolynomialRing(Qp, n+1)
    d = 2
    M = monomials_of_degree(gens(R), d)
    s = length(M)

    tot_time  = 0
    tot_error = 0

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
