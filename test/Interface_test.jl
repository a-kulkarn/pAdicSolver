
@testset "Interface test" begin
#if true
    Qp = PadicField(7, 20)
    L1, _ = LaurentSeriesField(GF(5), 15, "t")
    L2, _ = LaurentSeriesField(QQ, 5, "t")

    test_series_name = Dict(Qp=>"pAdic", L1=>"Laurent (Fp)", L2=>"Laurent (QQ)")

    for K in [Qp, L1, L2]
        @testset "Interface: $(test_series_name[K])" begin
        #if true
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

            #######################
            # Test tropical

            # tropsol = solve_affine_system(P, eigenvector_method=:tropical)
            # gsol = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)

            # @test solcmp(true_sol, sol, gsol)

        end
    end
end

