import Base:
    zero, one, show, print, length, getindex, setindex!, copy, promote_rule, convert, eltype, get
    *, /, //, -, +, ==, ^, divrem, conj, rem, real, imag, diff

export idx, rev_idx, length

mutable struct MonomialIdx
    terms::Dict{DynamicPolynomials.Monomial{true},Int64}
end

MonomialIdx() =  MonomialIdx(Dict{DynamicPolynomials.Monomial{true},Int64}())

get(p::MonomialIdx, m::DynamicPolynomials.Monomial{true}, df:: Int64) = get(p.terms, m, df)

function Base.setindex!(p::MonomialIdx, v::Int64, m::DynamicPolynomials.Monomial{true})
    (p.terms)[m] = v
end

length(l::MonomialIdx) = length(l.terms)

function idx(M)
    I = MonomialIdx()
    i = 1;
    for m in M
        I[m] = i
        i+=1
    end
    return I
end

function rev_idx(M)
    I = MonomialIdx()
    i = length(M);
    for m in M
        I[m] = i
        i-=1
    end
    return I
end
