
@testset "Interface test" begin
#if true
    Qp = PadicField(7, 20)
    L1, _ = LaurentSeriesField(GF(5), 15, "t")
    L2, _ = LaurentSeriesField(QQ, 5, "t")

    test_series_name = Dict(Qp=>"pAdic", L1=>"Laurent (Fp)", L2=>"Laurent (QQ)")

    for K in [Qp, L1, L2]
        if true
        #@testset "Interface: $(test_series_name[K])" begin
            R, (x1, x2) = PolynomialRing(K, 2)
            R3, (y1, y2, y3) = PolynomialRing(K, 3)
            t = uniformizer(K)

            P = [x1 - 1, x2^2 - 1^2]
            Q = [y1 - 1 * y3, y2^2 - y3^2]

            eins = one(ZZ)
            true_sol = matrix(K, [eins -eins; eins eins])
            true_sol_hom = matrix(K, [eins -eins eins; eins eins eins])


            # Test the main function with no keywords supplied.
            sol = solve_affine_system(P)
            @test solcmp(true_sol, sol)

            psol = solve_projective_system(Q)
            @test solcmp(true_sol_hom, psol)

            # Test groebner functions with no keywords.
            gsol = solve_affine_groebner_system(P)
            @test solcmp(true_sol, gsol)

            gpsol = solve_projective_groebner_system(Q)
            @test solcmp(true_sol_hom, gpsol)

            #######################
            # Test with keywords

            # Prepare the options dictionaries.
            method_opts = [:truncated_normal_form, :tnf, :macaulay]
            method_opts = [:given_groebner, :given_GB, :givenGB]

            evmethod_opts = [:power, :inverse, :schur]
            ordering_opts = [:lex, :deglex, :degrevlex]

            for meth in method_opts for evmeth in evmethod_opts for ord in ordering_opts

                test_series = "method = $meth, eigenvector_method = $evmeth, ordering = $ord"
                @testset "Options: $test_series" begin
                    sol = solve_affine_system(P, method = meth, eigenvector_method = evmeth,
                                              ordering = ord)
                    @test solcmp(true_sol, sol)
                    
                    psol = solve_projective_system(Q, method = meth,
                                                   eigenvector_method = evmeth,
                                                   ordering = ord)

                    #@info " " psol
                    @test solcmp(true_sol_hom, psol)
                end
            end end end


            # tropsol = solve_affine_system(P, eigenvector_method=:tropical)
            # gsol = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)

            # @test solcmp(true_sol, sol, gsol)

        end
    end
end

@testset "Basic plane curve tests" begin
    # For now, we need a fairly large prime. p=7 goes wrong fairly quickly.
    Qp = PadicField(491,3)

    R, (x1, x2) = PolynomialRing(Qp, ["x1", "x2"])
    d = 10
    P = random_square_system(R, d)
    sol = solve_affine_system(P)
    @test iszero(forward_error(P, sol))
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
