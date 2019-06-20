
using AbstractAlgebra

export @Ring
macro Ring(args...)
    X = DynamicPolynomials.PolyVar{true}[DynamicPolynomials.PolyVar{true}(string(arg)) for arg in args]
    V = [buildpolyvar(DynamicPolynomials.PolyVar{true}, args[i], X[i]) for i in 1:length(X)]
    push!(V, :(TMP = $X) )
    reduce((x,y) -> :($x; $y), V; init = :() )
end

function buildpolyvar(::Type{PV}, arg, var) where PV
    :($(esc(arg)) = $var)
end


DynPoly = DynamicPolynomials
GeneralDynPoly = DynPoly.Polynomial{B, S} where {B, S}

## Need new versions of the AA_* constructors for terms/monomials/etc...


export DPmonomials
DPmonomials = DynamicPolynomials.monomials

export AAPolynomialRing, AAPolynomial

function extract_variables(P)
    mons = p.m
    return unique(vcat( [m.vars for m in mons]...))    
end

"""
    Given a coefficient ring, which is a subtype of NCRing, return the AbstractAlgebra
"""
function AAPolynomialRing(coeff_ring::T where T <: AbstractAlgebra.NCRing,
                          X::Array{DynPoly.PolyVar{B},1} where B )
    
    return PolynomialRing(coeff_ring, [string(x) for x in X])
end

function AAPolynomialRing(P::T where T <:GeneralDynPoly )
    if DynPoly.degree(P)==0
        error("Cannot construct AbstractAlgebra parent from degree 0 polynomial")
    elseif !(P.a[1] <: AbstractAlgebra.NCRing || typeof(P.a[1]) != Int64)
        error("Coefficient type has no canonical parent")
    end

    if typeof(P.a[1]) == Int64
        R = FlintIntegerRing()
    else
        R = parent(a[1])
    end
    
    vars = extract_variables(P)

    return PolynomialRing(R, [string(v) for v in vars])
end

"""

    Converts a DynamicPolynomial to an AbstractAlgebra polynomial with parent `RX`. If no parent is
    provided, once is guessed using  "AAPolynomialRing".

    Default assignment of Dynamic Polynomial variables to the parent is greedy. However, one
    can provide a dictionary to make explicit the desired assignment.

"""

function AAPolynomial(P::T where T <:GeneralDynPoly )
    AAPolynomial(P, AAPolynomialRing(P))
end


# TODO: type assert the var_assignment input at some point.
function AAPolynomial(P::T where T <:GeneralDynPoly , RX; var_assignment=nothing)

    if var_assignment==nothing

        dvars = extract_variables(P)

        if size(dvars,1) > size(gens(RX),1)
            error("Input polynomial has more variables than adoptive parent ring")
        end

        varDic = Dict( zip(dvars,  gens(RX)[1:size(dvars)]) )
    else
        varDic = var_assignment
    end
    
    coeffs = [RX(c) for c in P.a]

    function coerce_monomial(m)
        e = m.z
        v = [ varDic[x] for x in m.vars]
        r = size(v,1)

        if r > 0
            return prod( v[i]^e[i] for i=1:r)
        else
            return RX(1)
        end
    end

    
    mons = [coerce_monomial(m) for m in P.x]
    r = size(mons,1)
    
    if iszero(P)
        return zero(RX)
    else
        return sum( mons[i]*coeffs[i] for i=1:r ) 
    end

end



    # if vars==nothing && parent==nothing
    #     if iszero(P) || !( typeof(P.a[1]) <: Hecke.NCRingElem )
    #         error("Cannot convert constant polynomial without an explicit parent for the coefficients.")
    #     end
    # end

    # if parent==nothing
    #     RX, RXvars = AAPolynomialRing(parent(P.a[1]), vars)
    # else
    #     vars = extract_variables(P):: T where T <: Array{DynPoly.PolyVar{true},1}
    #     RX, RXvars = AAPolynomialRing(parent(P.a[1]), vars)
    # end
