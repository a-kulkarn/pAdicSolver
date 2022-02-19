###############################################################################################
#
#    Groebner Hook tests
#
###############################################################################################

# NOTE: I don't know how to test this functionality within a "pkg> test pAdicSolver" session.
#       One can run this functionality manually using `include("groebner_hook.jl")`.

# Check if OSCAR is installed, and test the groebner hook.
import Pkg
using Test, Hecke, pAdicSolver
OSCAR_UUID = Base.UUID("f1435218-dba5-11e9-1e4d-f1a5fab5fc13")

if OSCAR_UUID in keys(Pkg.dependencies())
    
    # Add the hook
    import pAdicSolver._solver_groebner_hook
    import Oscar: ideal, groebner_basis

    function _solver_groebner_hook(P; ordering=ordering(parent(P[1])))
        return groebner_basis(ideal(P), ordering=ordering)
    end

    @testset "Groebner hook" begin

        K = QQ
        R, (x1, x2) = PolynomialRing(K, 2)
        R3, (y1, y2, y3) = PolynomialRing(K, 3)

        P = [x1 - 1, x2^2 - 1^2]
        Q = [y1 - 1 * y3, y2^2 - y3^2]

        eins = one(ZZ)
        true_sol = matrix(K, [eins -eins; eins eins])
        true_sol_hom = matrix(K, [eins -eins eins; eins eins eins])

        evmethod_opts = [:power, :inverse, :schur]
        ordering_opts = [:lex, :revlex, :grevlex, :deglex]

        for evmeth in evmethod_opts for ord in ordering_opts

            test_series = "eigenvector_method = $evmeth, ordering = $ord"
            @testset "Options: $test_series" begin
                sol = solve_system(P, method=:groebner, eigenvector_method=evmeth)
                @test solcmp(true_sol, sol)
            end
        end end
    end
end

nothing
