

using Hecke, pAdicSolver
TropicalInterval = pAdicSolver.TropicalInterval

toList(A) = [A[i, :] for i=1:size(A,1)]

    
function solcmp(args...)
    length(args) < 2 && throw(ArgumentError("Must compare at least two solution matrices."))

    # We can't use dictionaries due to isequal(x) being sensitive to precision (which
    # we do not intend to check for here). Instead, we expect that the length of these solution
    # arrays is reasonable, so we use a quadratic-time solution.
    
    X = toList(args[1])
    for arg in args
        Y = toList(arg)
        length(Y) != length(X) && return false

        # Cross out entries from Y one at a time until we are out.        
        for x in X
            i = findfirst(y->y==x, Y) # Very important to not use isequal.
            i == nothing && return false
            deleteat!(Y, i)
        end

        # If some elements remain, the sets are not equal.
        isempty(Y) || return false
    end
    return true

    #all(y in X for y in Y) || return false
    #all(x in X for x in X) || return false
    #return true

    # NOTE: This didn't work. Methinks there is an issue with a hash function and precision.
    # toSet(A) = Set(A[i, :] for i=1:size(A,1))
    # X = toSet(args[1])
    # for Y in args
    #     !issetequal(X, toSet(Y)) && return false
    # end
    # return true
end

K, t = LaurentSeriesRing(GF(5), 15, "t")
R, (x1, x2) = PolynomialRing(K, 2)

P = [x1 - t, x2^2 - 1^2]

eins = one(ZZ)
true_sol = matrix(K, [t -eins; t eins])


sol = solve_affine_system(P)
tropsol = solve_affine_system(P, eigenvector_method=:tropical)
sol = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)

@info " " sol true_sol
!isempty(sol) && @info " " collect(toSet(sol))[1] == collect(toSet(true_sol))[1]
!isempty(sol) && @info " " collect(toSet(sol))[2] == collect(toSet(true_sol))[2]
!isempty(sol) && @info " " collect(toSet(sol))[2] in collect(toSet(true_sol))

@assert solcmp(true_sol, sol)

# Base ring of QQ
K, t = LaurentSeriesRing(QQ, 15, "t")
R, (x1, x2) = PolynomialRing(K, 2)

P = [x1 - t, x2^2 - 1^2]
true_sol = matrix(K, [t -eins; t eins])


sol  = solve_affine_system(P)
solg = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)

@info " " sol true_sol solg

@assert solcmp(true_sol, true_sol, solg)

tropsol = solve_affine_system(P, eigenvector_method=:tropical)

trop_eins = TropicalInterval(fmpq(1), fmpq(1))
trop_null = TropicalInterval(fmpq(0), fmpq(0))
@assert tropsol == [trop_eins trop_null; trop_eins trop_null]

