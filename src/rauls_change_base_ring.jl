
##
#  Code introduced by Raul in the new version of AbstractAlgebra.
##


export rauls_change_base_ring


function (R::Singular.Rationals)(a::Hecke.fmpq)
    num = convert(BigInt, numerator(a))
    den = convert(BigInt, denominator(a))
    return R( num // den)
end

function (R::Hecke.FlintRationalField)(a::Singular.n_Q)
    n = convert(BigInt, numerator(a))
    d = convert(BigInt, denominator(a))
    return R(n//d)
end

function (R::Hecke.FlintRationalField)(a::Singular.n_Z)
    n = convert(BigInt, a)
    return R(n)
end


@doc Markdown.doc"""
    change_base_ring(p::AbstractAlgebra.MPolyElem{T}, g, R::MPolyRing)
       where T <: RingElement
> Return the polynomial in R obtained by applying g to the coefficients of p.
"""
function rauls_change_base_ring(p::Hecke.MPolyElem{T}, g, R::Hecke.MPolyRing) where {T <: RingElement}

   z = g(zero(base_ring(p.parent)))
   base_ring(R) != parent(z) && error("Base rings do not match.")

   cvzip = zip(coeffs(p), exponent_vectors(p))
   M = MPolyBuildCtx(R)
   for (c, v) in cvzip
      push_term!(M, g(c), v)
   end

   return finish(M)
end

