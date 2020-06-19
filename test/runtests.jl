
using Test
using Dory, pAdicSolver, Hecke, Singular

function row_set(M)
    return Set(M[i,:] for i=1:size(M,1))
end


@testset "Basic solving (Singular)" begin

    Qp = PadicField(433,3)

    R, (x1,x2) = Singular.PolynomialRing(Singular.QQ, ["x1", "x2"])
    n = length(gens(R))

    # Create the simple polynomials.
    P = [x2 - (x1^2 - 1), x2]

    sol_mat = matrix(Qp, [[-1,1], [0,0]])
    low_prec_sol_mat = (x->setprecision(x,2)).(sol_mat)
    
    true_sol_set = row_set(sol_mat)
    sol_set_low_prec = row_set(low_prec_sol_mat)

    #Set(matrix(Qp,1,2,[setprecision(ai,2) for ai in a]) for a in true_sol_set)
    
    # Launch solver
    I = Singular.Ideal(R, P)
    
    @test true_sol_set == row_set(padic_solutions(I, Qp, eigenvector_method="qr"))
    @test true_sol_set == row_set(padic_solutions(I, Qp, eigenvector_method="power"))
    # @test sol_set_low_prec == row_set(padic_solutions(I, Qp, eigenvector_method="inverse"))
    # @test true_sol_set == row_set(padic_solutions(I, Qp, eigenvector_method="classical"))

end

@testset "Basic solving (TNF)" begin

    Qp = PadicField(433,3)

    R, (x1,x2) = Hecke.PolynomialRing(Qp, ["x1", "x2"])
    n = length(gens(R))

    # Create the simple polynomials.
    P = [x2 - (x1^2 - 1), x2-(x1+1)]

    sol = padic_solutions(P, Qp)
    @test rel_error(P,sol) == Hecke.zero_matrix(Qp,2,2).entries
end

@testset "Verbose printing" begin
    Qp = PadicField(433,3)

    R, (x1,x2) = Hecke.PolynomialRing(Qp, ["x1", "x2"])
    n = length(gens(R))

    # Create the simple polynomials.
    P = [x2 - (x1^2 - 1), x2-(x1+1)]

    sol = padic_solutions(P, Qp)
    @test rel_error(P,sol) == Hecke.zero_matrix(Qp,2,2).entries
end

@testset "TNF solving with a condensed Macaulay matrix (missing monomials)" begin

    Qp = PadicField(433,3)

    R, (x1,x2) = Hecke.PolynomialRing(Qp, ["x1", "x2"])

    # Create the simple polynomials.
    P = [x2 - (x1^2 - 1), x2]

    sol_mat = matrix(Qp, [[-1,1], [0,0]])
    low_prec_sol_mat = (x->setprecision(x,2)).(sol_mat)
    
    true_sol_set = row_set(sol_mat)
    sol_set_low_prec = row_set(low_prec_sol_mat)


    @info "" padic_solutions(P, Qp, eigenvector_method="qr")
        
    @test true_sol_set == row_set(padic_solutions(P, Qp, eigenvector_method="qr"))
    @test true_sol_set == row_set(padic_solutions(P, Qp, eigenvector_method="power"))
    @test sol_set_low_prec == row_set(padic_solutions(P, Qp, eigenvector_method="inverse"))
    # @test true_sol_set == row_set(padic_solutions(I, Qp, eigenvector_method="classical"))
end


@testset "Tropical solving (Singular)" begin
    include("Seibenundzwanzig_sing.jl")
end

@testset "Tropical solving (TNF)" begin
    include("Seibenundzwanzig.jl")
end

nothing
