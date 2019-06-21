
export rel_error

# Mourrain's relative error function
#
# function rel_error(P, Xi::Matrix, X = variables(P))
#     r = fill(0.0, size(Xi,1), length(P))
#     n = size(Xi,2)
#     for i in 1: size(Xi,1)
#         for j in 1:length(P)
#             V = Xi[i,:]
#             r[i,j]= norm(subs(P[j],X=>V))
#             s = 0.0
#             for t in P[j]
#                 s+= norm(subs(t, X => V))
#             end
#             r[i,j]/=s
#         end
#     end
#     return r
# end


"""
Gives the list of evaluations of the polynomials in `P` at the points `sol`.
"""
function rel_error(P,sol)
    return [evaluate(p, sol.entries[i,:]) for i in 1:size(sol,1),  p in P]
end
