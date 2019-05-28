

## Dispatcher for nullspace
function nullspace(A::Array{padic,2})
    M = matrix(A)
    return Hecke.nullspace(M)[2].entries
end

struct QRPadicArrayPivoted
    Q::Array{padic,2}
    R::Array{padic,2}
    p::Array{Int64,1}
end


# modify the norm function to make sense
import LinearAlgebra.norm
norm(A :: Array{padic,1}) = maximum( abs.(A))


# AVI:
# Function to compute the eigenvalues of a list of (commuting) matrices, in the
# specific case that the matrices are mult-by-coordinate-variable operators on R/I
#
# INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
# Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
function normalized_simultaneous_eigenvalues(
    inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)

    M = [ matrix(A) for A in inputM]
    Qp = base_ring(M[1])
    Mrand = sum(A*rand_padic_int(Qp) for A in M) # non-unit random causes problems

    println("Valuations of singular values: ")
    println(valuation.(singular_values(M[1])))

    # We will assume that M[1] is well-conditioned. for now.
    
    # The rectangular solve step is enough to kill off any helpful data mod p.
    @time I0 = rectangular_solve(M[1],identity_matrix(Mrand.base_ring,size(M[1],1)))
    @time Mg = I0*Mrand

    # eigen vectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
    eigen_subspaces  = eigspaces(Mg, method=method)
    
    #println("eigvalues: ", invariant_subspaces.values)
    #println()
    #println("eigspaces: ", length(invariant_subspaces.spaces))


    X = matrix( Qp, fill(zero(Qp), length(eigen_subspaces.spaces) ,length(M)))
    for j in 1:length(M)
        for i in 1:length(eigen_subspaces.spaces)
            
            V = eigen_subspaces.spaces[i]            
            Y = rectangular_solve(V,I0*M[j]*V)           
            X[i,j] = trace(Y)/Qp(size(Y,2))
        end
    end

    function normalize_solution!(Xi, ish)
        Sol = Xi
        
        if (!ish)
            i=1            
            while i <= size(Sol,1)
                scale_factor = Sol[i,1]
                if iszero(scale_factor)
                    println()
                    println("!-- Ignoring solution at infinity --!")
                    println()
                    
                    Sol = vcat(Sol[1:(i-1), :], Sol[i+1:size(Sol,1), :])
                else
                    Sol[i,:] *= inv(scale_factor)
                    i+=1
                end
            end
        #else
            # do nothing otherwise, for now.
        end
        return Sol
    end
    
    return normalize_solution!(X, ish)
end
