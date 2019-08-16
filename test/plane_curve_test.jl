
using Test
using Dory
using pAdicSolver

# For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
Qp = PadicField(491,3)


R, (x1,x2) = PolynomialRing(Qp, ["x1", "x2"])
n = length(gens(R))

d = 10
M = monomials_of_degree(gens(R),0:d)
s = length(M)

failed_test_count = 0
for i=1:10
    try
        global P = hcat( [ [ 2*rand_padic_int(Qp)- 1 for i in 1:n] for j in 1:s]... )*M

        #matlist, F, B, N, Nr, RR, IdL0, IdL = solve_macaulay(P,gens(R), test_mode=true);
        sol = solve_macaulay(P)

        println("\n-- sol ")
        println(sol,"\n")

        Er = rel_error(P,sol)
        println("-- Rel error: ")
        display(Er)
        println()
        
    catch
        global failed_test_count += 1
    end
end

println("----------------------")
println("Tests failed due to stability issues (probably): ", failed_test_count)

@test failed_test_count==0
