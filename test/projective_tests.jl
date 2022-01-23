
using Dory

@testset "Projective solving tests" begin

    Qp = PadicField(491,10)

    R, (x1,x2,x3) = PolynomialRing(Qp, ["x1", "x2","x3"])
    n = length(gens(R))-1


    for d in [1, 2, 5]
        M = monomials_of_degree(gens(R), d)
        s = length(M)

        P = hcat([[2*rand_padic_int(Qp) - 1 for i in 1:n] for j in 1:s]...)*M
        sol = solve_macaulay(P)
        Er = rel_error(P,sol)

        @test iszero(Er)
    end
    
    R, (x0,x1,x2,x3,x4,x5) = PolynomialRing(Qp, ["x0","x1", "x2","x3","x4","x5"])
    n = length(gens(R))-1

    d = 2
    M = monomials_of_degree(gens(R), d)
    s = length(M)

    P = hcat([[2*rand_padic_int(Qp) - 1 for i in 1:n] for j in 1:s]...)*M
    sol = solve_macaulay(P)
    Er = rel_error(P,sol)

    @test iszero(Er)
end

