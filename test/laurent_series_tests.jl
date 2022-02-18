

using Hecke, pAdicSolver

K, t = LaurentSeriesRing(GF(5), 15, "t")
R, (x1, x2) = PolynomialRing(K, 2)

P = [x1 - t, x2^2 - 1^2]

sol = solve_affine_system(P)
tropsol = solve_affine_system(P, eigenvector_method=:tropical)
sol = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)


K, t = LaurentSeriesRing(QQ, 15, "t")
R, (x1, x2) = PolynomialRing(K, 2)

P = [x1 - t, x2^2 - 1^2]
sol = solve_affine_system(P)
tropsol = solve_affine_system(P, eigenvector_method=:tropical)
sol = solve_affine_groebner_system(P, ordering=:lex, eigenvector_method=:schur)
