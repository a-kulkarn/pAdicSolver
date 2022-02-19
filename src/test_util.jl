
export rel_error, solcmp

"""
Gives the list of evaluations of the polynomials in `P` at the points `sol`.
"""
function rel_error(P,sol)
    return [evaluate(p, sol.entries[i,:]) for i in 1:size(sol,1),  p in P]
end


@doc Markdown.doc"""
    solcmp(args...)

Compare the solution coordinates as sets with multiplicity.
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

    #all(y in X for y in Y) || return false
    #all(x in X for x in X) || return false
    #return true

    # NOTE: This didn't work. Methinks there is an issue with a hash function and precision.
    # toSet(A) = Set(A[i, :] for i=1:size(A,1))
    # X = toSet(args[1])
    # for Y in args
    #     !issetequal(X, toSet(Y)) && return false
    # end
    # return true
end
