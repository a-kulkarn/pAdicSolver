using Singular

function random_linear_equations(R::Singular.PolyRing)
    var_vec   = Vector(gens(R))
    rand_vec  = Vector([rand(-1000:1000) for i in 1:length(gens(R))])

    return (transpose(var_vec)*rand_vec)
end

R, Rvars =
    Singular.PolynomialRing(Singular.QQ,
    #PolynomialRing(Singular.N_ZpField(32003),
                   ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12",
                    "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20", "x21", "x22", "x23",
                    "x24", "x25", "x26", "x27", "x28", "x29", "x30", "x31", "x32", "x33", "x34",
                    "x35", "x36", "x37", "x38", "x39", "x40"], ordering=:degrevlex)

(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40) = Rvars

include(".lgr_eqns.jl")


sid   = Singular.Ideal(R, [eval(Meta.parse("$p")) for p in ps])
sidl  = Singular.Ideal(R, [random_linear_equations(R) for i = 1:11])
sidl2 = Singular.Ideal(R, [ Rvars[j] - R(rand(-1000:1000)) for j=vcat(26:28, 30, 34:40) ])

@time sid = sid + sidl2
@time G = Singular.std(sid);

# basis of quotient R/id
B = Singular.kbase(G)


# Convert to the right type.
using Hecke
using pAdicSolver

#P = map(f->rauls_change_base_ring(f,FlintQQ, Hecke.PolynomialRing(FlintQQ, 40, ordering=:degrevlex)[1]), gens(G))

#P = gens(G)

#sol = solve_macaulay(P, groebner=true, eigenvector_method="tropical")

sol = padic_solutions(G, PadicField(23,100))

nothing
