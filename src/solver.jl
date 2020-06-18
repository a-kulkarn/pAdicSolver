export padic_solutions

import Hecke.ishomogeneous
function ishomogeneous(I::Singular.sideal)
    G = Singular.gens(I)
    return all(ishomogeneous, G)
end

function ishomogeneous(f::Singular.spoly)
    d = total_degree(first(monomials(f)))
    return all(m->total_degree(m)==d, monomials(f))
end


function truncated_normal_form(I::Singular.sideal{<:Singular.spoly{<:T}} where T)
    return Singular.std(I)
end

"""
Returns the affine solutions to the polynomial equations defined by I.
The system must be 0 dimensional.
"""
function padic_solutions(I::Singular.sideal{<:Singular.spoly{<:T}} where T, R;
                         eigenvector_method="power")

    # Presently, the user is responsible for providing the padic ring R.
    #
    # Question: Is it possible to decide the right precision for the user at this stage?

    base_ring = Singular.base_ring(I)
    X  = gens(base_ring)
    GB = truncated_normal_form(I)
    ish = ishomogeneous(GB)

    ## Construct the multiplication matrices directly.
    B = Singular.gens(Singular.kbase(GB))

    xi_operators = []
    for v in vcat([base_ring(1)], X)
        push!(xi_operators, [Singular.reduce(v*b, GB) for b in B] )
    end

    
    M = [matrix(R, [[coeff(g, b) for b in B] for g in op]) for op in xi_operators]    
    M = [m.entries for m in M]

    # @info "" M[1] M[2] M[3]
    
    t0 = time()
    Xi = normalized_simultaneous_eigenvalues(M, ish, eigenvector_method)
    
    # In the affine system, the distinguished monomial (i.e, "1" for that case) does 
    # not correspond to a coordinate.
    if ish return Xi else return  Xi[:,2:size(Xi,2)] end
end

function padic_solutions(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1}, R;
                         eigenvector_method="power")
    
    @assert base_ring(parent(P[1])) == R
    return solve_macaulay_II(P, eigenvector_method=eigenvector_method)
end
