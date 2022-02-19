

@testset "Projective solving tests" begin

    Qp = PadicField(491,10)
    R, (x1,x2,x3) = PolynomialRing(Qp, 3)

    for d in [1, 2, 5]
        P = random_square_system(R, d, homogeneous=true)
        sol = solve_projective_system(P)
        @test iszero(forward_error(P, sol))
    end

    # A large system of degree 2.
    R, vars = PolynomialRing(Qp, 6)
    d = 2
    
    P = random_square_system(R, d, homogeneous=true)
    sol = solve_projective_system(P)
    @test iszero(forward_error(P, sol))
end

