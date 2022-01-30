

@testset "Tropical solver tests (part 2)" begin


    p = 1000003
    N = 15
    K = Hecke.PadicField(p, N)
    R, (x1, x2, x3, x4, x5) = Hecke.PolynomialRing(K, 5)
    
    GB = [2*x3 - 39000117*x4 + 5*x5,
          1000003*x2 - 1,
          x1 - 1,
          21996131716197184*x4*x5 + 35100208458309461*x4 - 3158009474*x5^2 - 5232015125*x5 - 985,
          4182039634666882990820525906472*x4^2 + 1597064285781592634326475*x4 - 92060053845332998*x5^2 - 241723414460067851*x5 - 56121168503,
          (16324044970945243992612460993471986029743497621*x4 + 1635912215634320703145614369145738988672*x5^3 + 3145181581802732307984319762588946014782*x5^2
          - 799304976455501742146405393837628704021*x5 + 4579907243568677265846766872644839)]

    sol1 = pAdicSolver.solve_affine_system(GB, eigenvector_method=:tropical)
    sol2 = pAdicSolver.solve_affine_groebner_system(GB, eigenvector_method=:tropical, ordering=:degrevlex)
    
    # The result we are supposed to get
    true_solution = let
        mat = [0 -1 -1 -2 -1;
               0 -1 -1 -2 -1;
               0 -1 -1 -1 -1;
               0 -1 0 0 0]
        
        fmpq.(mat)
    end

    # M = pAdicSolver._multiplication_matrices(Val(:given_GB), curry.(GB), false, ordering=:degrevlex)
    # sol = pAdicSolver.simultaneous_eigenvalues([matrix(A) for A in M], method=:tropical)

    #@info " " min.(sol1) min.(sol2)

    @test min.(sol1) == min.(sol2) == true_solution


    Kx,(x1,x2,x3) = PolynomialRing(QQ,3)
    p = 32003
    p_adic_precision=29
    G =  [x2 + p * x3, x1 - p * x3 + 1, 1024128004*x3^2 - 1024224011*x3]

    Qp = PadicField(p,p_adic_precision)
    Gp = [change_base_ring(Qp,f) for f in G]
    sol = pAdicSolver.solve_affine_groebner_system(Gp, eigenvector_method=:tropical, ordering=:degrevlex)

    true_solution = let
        mat = [0  1  0;
               0 29 29]
        
        fmpq.(mat)
    end

    true_error = let
        O = fmpq(0)
        [O 1 O;
         O Inf Inf]
    end
    
    @test min.(sol) == true_solution
    @test max.(sol) == true_error
end
