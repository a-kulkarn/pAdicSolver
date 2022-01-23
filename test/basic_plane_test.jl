
using Dory

@testset "Basic plane curve tests" begin

    # For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
    Qp = PadicField(491,3)


    R, (x1,x2) = PolynomialRing(Qp, ["x1", "x2"])
    n = length(gens(R))

    d = 10
    M = monomials_of_degree(gens(R),0:d)
    s = length(M)

    failed_test_count = 0

    P = hcat([[2*rand_padic_int(Qp) - 1 for i in 1:n] for j in 1:s]...)*M

    sol = solve_macaulay(P)

    println("\n-- sol ")
    println(sol,"\n")

    Er = rel_error(P,sol)
    println("-- Rel error: ")
    display(Er)
    println()
    @test iszero(Er)
end


#=
Some primes:
2      3      5      7     11     13     17     19     23     29 
31     37     41     43     47     53     59     61     67     71 
73     79     83     89     97    101    103    107    109    113 
127    131    137    139    149    151    157    163    167    173 
179    181    191    193    197    199    211    223    227    229 
233    239    241    251    257    263    269    271    277    281 
283    293    307    311    313    317    331    337    347    349 
353    359    367    373    379    383    389    397    401    409 
419    421    431    433    439    443    449    457    461    463 

=#
