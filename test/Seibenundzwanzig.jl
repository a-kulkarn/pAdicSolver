

@testset "Tropical solver tests" begin

    include(".cubic_surface_eqns.jl")

    sol = solve_macaulay(P, rho=5, eigenvector_method = :tropical)
    #matlist, F, B, N, Nr = solve_macaulay(P, gens(R), rho=4, test_mode=true);

    #A = matrix(matlist[1]);
    #E = eigspaces(A);

    # println("\n-- sol ")
    # println(sol,"\n")


    #######
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

    @test sol == tropical_test_data

end 
