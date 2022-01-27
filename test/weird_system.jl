
using Hecke


R, (x1, x2, x3, x4, x5) = PolynomialRing(QQ, 5)
#S, (x1, x2, x3, x4, x5) = PolynomialRing(QQ, 5, ordering=:degrevlex)
#R, (x4, x5) = PolynomialRing(QQ, 5)

prime = 1000003
#prime = 233
prec = 15

GB = [2*x3 - 39000117*x4 + 5*x5,
      1000003*x2 - 1,
      x1 - 1,
      21996131716197184*x4*x5 + 35100208458309461*x4 - 3158009474*x5^2 - 5232015125*x5 - 985,
      4182039634666882990820525906472*x4^2 + 1597064285781592634326475*x4 - 92060053845332998*x5^2 - 241723414460067851*x5 - 56121168503,
      16324044970945243992612460993471986029743497621*x4 + 1635912215634320703145614369145738988672*x5^3 + 3145181581802732307984319762588946014782*x5^2 - 799304976455501742146405393837628704021*x5 + 4579907243568677265846766872644839]


lp = [x3, x2,  x1,  x4 * x5,  x4^2, x5^3]
#lp = [x4 * x5,  x4^2, x5^3]

K = PadicField(prime, prec)


#steve = f->pAdicSolver.change_parent(S, f)
curry = f->change_base_ring(K,f)

nothing

#pAdicSolver.solve_affine_groebner_system(curry.(GB), curry.(lp))
#pAdicSolver.solve_affine_system(curry.(GB), eigenvector_method="tropical")


# The result we are supposed to get
true_solution = convert(Matrix{Rational{BigInt}}, hcat([ 0, -1, 0, 0, 0 ], [0, -1, -1, -1, -1 ], [ 0, -1, -1, -2, -1 ], [0, -1, -1, -2, -1 ]))
