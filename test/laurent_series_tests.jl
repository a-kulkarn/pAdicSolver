

using Hecke, pAdicSolver
TropicalInterval = pAdicSolver.TropicalInterval

@testset "Basic systems over Laurent series fields" begin

    K, t = LaurentSeriesRing(GF(5), 15, "t")
    R, (x1, x2) = PolynomialRing(K, 2)

    P = [x1 - t, x2^2 - 1^2]

    eins = one(ZZ)
    true_sol = matrix(K, [t -eins; t eins])

    sol = solve_affine_system(P)
    tropsol = solve_affine_system(P, eigenvector_method=:tropical)
    gsol = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)

    @test solcmp(true_sol, sol, gsol)

    # Base ring of QQ
    K, t = LaurentSeriesRing(QQ, 15, "t")
    R, (x1, x2) = PolynomialRing(K, 2)

    P = [x1 - t, x2^2 - 1^2]
    true_sol = matrix(K, [t -eins; t eins])


    sol  = solve_affine_system(P)
    solg = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)

    @test solcmp(true_sol, sol, solg)

    tropsol = solve_affine_system(P, eigenvector_method=:tropical)

    trop_eins = TropicalInterval(fmpq(1), fmpq(1))
    trop_null = TropicalInterval(fmpq(0), fmpq(0))
    @test tropsol == [trop_eins trop_null; trop_eins trop_null]

end
