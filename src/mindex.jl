import Base:
    zero, one, show, showcompact, print, length, endof, getindex, setindex!, copy, promote_rule, convert, start, next, done, eltype, get
    *, /, //, -, +, ==, ^, divrem, conj, rem, real, imag, diff, norm

export idx, rev_idx, length

mutable struct MonomialIdx
    terms::Dict{DynamicPolynomials.Monomial{true},Int64}
end

MonomialIdx() =  MonomialIdx(Dict{DynamicPolynomials.Monomial{true},Int64}())

get(p::MonomialIdx, m::DynamicPolynomials.Monomial{true}, df:: Int64) = get(p.terms, m, df)

function setindex!(p::MonomialIdx, v, m::Monomial{true})
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
