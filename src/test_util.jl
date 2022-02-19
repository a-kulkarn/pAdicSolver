
export abs_error, rel_error, forward_error, solcmp, random_square_system

@doc Markdown.doc"""
    solcmp(args...)

Compare the solution coordinates as sets with multiplicity. Coordinates are compared with
`==` rather than `isequal`.
"""
function solcmp(args...)
    length(args) < 2 && throw(ArgumentError("Must compare at least two solution matrices."))

    # We can't use dictionaries due to isequal(x) being sensitive to precision (which
    # we do not intend to check for here). Instead, we expect that the length of these 
    # solution arrays is reasonable, so we use a quadratic-time solution.
    
    X = [args[1][i, :] for i=1:size(args[1], 1)]
    for arg in args
        Y = [arg[i, :] for i=1:size(arg, 1)]
        length(Y) != length(X) && return false

        # Cross out entries from Y one at a time until we are out.        
        for x in X
            i = findfirst(y->y==x, Y) # Very important to not use isequal.
            i == nothing && return false
            deleteat!(Y, i)
        end

        # If some elements remain, the sets are not equal.
        isempty(Y) || return false
    end
    return true
end


################################################################################################
#
#  pNumerical error
#
################################################################################################

"""
Gives the list of evaluations of the polynomials in `P` at the points `sol`.
"""
function abs_error(P, sol)
    return [evaluate(p, sol.entries[i,:]) for i in 1:size(sol,1),  p in P]
end

# TODO: Supported for legacy reasons, but technically this is wrong.
rel_error = abs_error

forward_error = abs_error

################################################################################################
#
#  Random generation
#
################################################################################################

@doc Markdown.doc"""
    random_square_system(R, d; homogeneous = false)

Given a polynomial ring `R` and a positive integer `d`, return a random polynomial system
with `ngens(R)` equations of degree `d`. If `homogeneous = true`, return a random homogeneous
system with `ngens(R)-1` equations.
"""
function random_square_system(R::Hecke.Generic.MPolyRing, d::Int64; homogeneous = false)

    K = base_ring(R)

    if homogeneous
        n = length(gens(R))-1
        degs = d
    else
        n = length(gens(R))
        degs = 0:d
    end
    
    M = monomials_of_degree(gens(R), degs)
    
    return [randint(K) for i in 1:n, j in 1:length(M)] * M
end
