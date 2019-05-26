
# NOTE: Probably should be moved to "Dory" at some point.

# Return the monomials in the variables specified by X of degrees
# given by itr. itr can either be a list, iterator, or a single number.
function monomials_of_degree(X,itr)
    return vcat([monomials_of_degree(X,d) for d in itr]...)
end

function monomials_of_degree(X,d::Int64)

    if isempty(X)
        error("Variables required in input")
    end

    if length(X) == 1
        return [X[1]^d]
    end

    if d<0
        return []
    end

    if d==0
        return [parent(X[1])(1)]
    end

    all_monomials = []

    n = length(X)
    Y = X[1:n-1]
    for j=0:d
        push!(all_monomials, [X[n]^j*m for m in monomials_of_degree(Y,d-j)])
    end
    return vcat(all_monomials...)
end

# TODO: Make a projective version of this
function dense_coefficients(f)
    R = parent(f)
    mons = monomials_of_degree(gens(R), 0:total_degree(f))
    return [coeff(f, m) for m in mons]
end

function coeff(f::AbstractAlgebra.MPolyElem{T}, L) where T <: RingElement
    return [AbstractAlgebra.coeff(f,m) for m in L]
end
