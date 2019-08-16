
export rel_error, some_primes

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


@doc Markdown.doc"""
    rel_error(P,sol) --> ::Array{T,1} where T

Gives the list of evaluations of the polynomials in `P` at the points `sol`.
"""
function rel_error(P,sol)
    return [evaluate(p, sol.entries[i,:]) for i in 1:size(sol,1),  p in P]
end


@doc Markdown.doc"""
    some_primes() --> Array{Int64,2}

Return a 2-dimensional array of the primes up to 464.
"""
function some_primes()
    Array{Int64,2}(
   [  2      3      5      7     11     13     17     19     23     29; 
     31     37     41     43     47     53     59     61     67     71; 
     73     79     83     89     97    101    103    107    109    113; 
    127    131    137    139    149    151    157    163    167    173; 
    179    181    191    193    197    199    211    223    227    229; 
    233    239    241    251    257    263    269    271    277    281; 
    283    293    307    311    313    317    331    337    347    349; 
    353    359    367    373    379    383    389    397    401    409; 
    419    421    431    433    439    443    449    457    461    463 ] 
    )
end
