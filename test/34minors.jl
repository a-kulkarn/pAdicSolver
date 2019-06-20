

p = 13
Qp = PadicField(p, 50)

vars_str = [ "x1",  "x2",  "x3",  "x4",  "x5",  "x6",  "x7",  "x8",  "x9",  "x10", "x11", "x12"]
R, (x1,  x2,  x3,  x4,  x5,  x6,  x7,  x8,  x9,  x10,  x11,  x12) = PolynomialRing(Qp, vars_str)


# 3x4 matrix of variables.
var_vec = matrix([x1  x2  x3  x4  x5  x6  x7  x8  x9  x10  x11  x12])
L = [x1  x2  x3  x4;  x5  x6  x7  x8;  x9  x10  x11  x12]

A = matrix( R,  L)

function random_linear()
    M = matrix( hcat([[ R(p)^rand(0:15) for i in 1:12] for j in 1:2]...))
    return (var_vec*M)[1,1]
end

h1 = random_linear()
h2 = random_linear()

eqns = vcat(minors(A, 3), [h1, h2]);

# Prevent printing everything.
nothing


