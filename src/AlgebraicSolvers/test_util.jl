
export rel_error

function rel_error(P,sol)
    return [evaluate(p, sol.entries[i,:]) for i in 1:size(sol,1),  p in P]
end
