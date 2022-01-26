
using Oscar
import pAdicSolver

#@testset "Groebner option tests" begin

function kbase_from_GB(I)
    Oscar.singular_assure(I)
    sI = I.gens.S
    sI.isGB = true

    return sing_kbase = Oscar.Singular.kbase(sI)
end


if true

    R, Rvars = PolynomialRing(QQ, ["x1", "x2", "x3"])
    (x1, x2, x3) = Rvars

    I = ideal(R, [x1^2 - x2^2 + x3^2, x1^2 + x1 * x2 + 3 * x2^2 + 5 * x3^2, x1^3 + x2^3 + x3^3 - 1])
    Oscar.singular_assure(I)
    
    gb = groebner_basis(I)
    gblex = groebner_basis(I, ordering=:lex)
    gbdegrevlex = groebner_basis(I, ordering=:degrevlex)

    # We test this because kbase more or less ignores the ordering of the AbstractAlgebra parent.
    @test gb == gbdegrevlex
    @test gb != gblex
    
    I1 = ideal(R, gb)
    I2 = ideal(R, gblex)

    # Check to ensure that kbase is IGNORING the ordering of the parent, and only works
    # with respect to a :degrevlex Groebner basis
    @test length(gens(kbase_from_GB(I1))) == 12
    @test length(gens(kbase_from_GB(I2))) == 1

    # Raise error if the input parent does not have :degrevlex ordering.
    
    # basis of quotient R/id

    # B = Singular.kbase(gb)
    nothing
end


if true

    sol = pAdicSolver.solve_affine_groebner_system(gblex, leading_monomial.(gblex))

end 
