import Base: zero, one, show, print, length, getindex, setindex!, copy, promote_rule, convert, eltype, get
import Base: *, /, //, -, +, ==, ^, divrem, conj, rem, real, imag, diff

mutable struct MonomialIdx
    terms::Dict{Dory.MPolyElem{T} where T, Int64} 
end

MonomialIdx() =  MonomialIdx(Dict{Dory.MPolyElem{T} where T,Int64}())

get(p::MonomialIdx, m::Dory.MPolyElem{T} where T, df:: Int64) = get(p.terms, m, df)

function Base.setindex!(p::MonomialIdx, v::Int64, m)
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
