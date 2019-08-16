
using Test
using Dory
using pAdicSolver

p = 101
Qp = PadicField(p, 70)

# vars_str = [ "x1",  "x2",  "x3",  "x4",  "x5",  "x6",  "x7",  "x8",  "x9",  "x10", "x11", "x12"]
# R, (x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8,  x9,  x10,  x11,  x12) = PolynomialRing(Qp, vars_str)


# # 3x4 matrix of variables.
# var_vec = matrix([x1  x2  x3  x4  x5  x6  x7  x8  x9  x10  x11  x12])
# L = [x1  x2  x3  x4;
#      x5  x6  x7  x8;
#      x9  x10  x11  x12]

function wrap(RX)
    vars = RX[2]
    for x in vars
        global :($v) = x
    end
    return RX[1]
end

R, (x1,x2,x3,x5,x6) = PolynomialRing(Qp, ["x1","x2","x3","x5","x6"])
var_vec = matrix(R, [1 x1 x2 x3 x5 x6])

L = [x1 x2 x3 2;
      1 x5 x6 1;
      1  1  1 1]


A = matrix( R,  L)

function random_linear()
    M = matrix( hcat([[ R(p)^rand(0:15) for i in 1:size(var_vec,2)] for j in 1:2]...))
    return (var_vec*M)[1,1]
end

function affine_random_linear()
    M = matrix( R, [ [rand_padic_int(Qp) for i in 1:size(var_vec,2)] for j in 1:2] )
    return (var_vec*M)[1,1]
end

h1 = affine_random_linear()
h2 = affine_random_linear()

eqns = vcat(minors(A, 3), [affine_random_linear() for i=1:1]);

# Presently, the inefficiency of the algorithm makes this inadvisable to run.
# The positive dimensional components are definitely causing an issue here.
# This is actually quite a bad example from the perspective of the solver.
sol = solve_macaulay(eqns, rho=4, eigenvector_method="tropical")


# Prevent printing.
nothing
