export padic_solutions

function truncated_normal_form(I::Singular.sideal{<:Singular.spoly{<:T}} where T)
    return Singular.std(I)
end

function padic_solutions(I::Singular.sideal{<:Singular.spoly{<:T}} where T, R;
                         eigenvector_method="power")

    @info "Computation started..."

    base_ring = Singular.base_ring(I)
    X  = gens(base_ring)
    GB = truncated_normal_form(I)
    ish = false

    ## Construct the multiplication matrices directly.
    B = Singular.gens(Singular.kbase(GB))

    xi_operators = []
    for v in vcat([base_ring(1)], X)
        push!(xi_operators, [Singular.reduce(v*b, GB) for b in B] )
    end

    
    M = [matrix(R, [[coeff(g, b) for b in B] for g in op]) for op in xi_operators]
    
    M = [m.entries for m in M]

    @info "" M[1] M[2] M[3]
    
    ## The prime will be specified by the user..
    ## The precsion should also be specified, or the user should request
    ## some feature to be invoked.

    # Question: How to decide the right precision for the user at this stage???        

    t0 = time()
    Xi = normalized_simultaneous_eigenvalues(M, ish, eigenvector_method)
    println("-- Eigen diag",  "   ",time()-t0, "(s)"); t0 = time()

    # In the affine system, the distinguished monomial (i.e, "1" for that case) does 
    # not correspond to a coordinate.
    if ish return Xi else return  Xi[:,2:size(Xi,2)] end
end

function padic_solutions(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1})
    return true
end
