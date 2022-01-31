

@testset "Tropical solver tests" begin

    include(".cubic_surface_eqns.jl")

    sol = solve_macaulay(P, rho=5, eigenvector_method = :tropical)
    
    ## We should get this output running the tropical solver.
    tropical_test_data = [
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        5//4  3//4  0//1   5//4  1//2;
        2//1  1//1  0//1   2//1  1//1;        
        2//1  1//1  0//1   2//1  1//1;
        6//1  3//1  0//1   6//1  3//1]

    @test min.(sol) == tropical_test_data

end 
