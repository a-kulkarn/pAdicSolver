
using Hecke

# Somehow enforce this dependency.
include("padic_util.jl")

####################

# dispatch to the p-adic qr utilities
function qr(A :: Array{padic,2}, pivot=Val(true))
    return padic_qr(matrix(A))
end


# AVI:
# Function to compute the eigenvalues of a list of (commuting) matrices, in the
# specific case that the matrices are mult-by-coordinate-variable operators on R/I
#
# INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
# Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
function normalized_simultaneous_eigenvalues(
    inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool)

    M = [ matrix(A) for A in inputM]
    Qp = base_ring(M[1])
    M0 = sum(A*rand(Qp) for A in M) # non-unit random causes problems

    I0 = inv(M0)
    Mg = I0*M[1]

    # eigen vectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
    E  = eigvecs(Mg) 

    X = matrix( Qp, fill(zero(Qp), size(E,2),length(M)))
    for j in 1:length(M)
        Yj = rectangular_solve(E, I0*M[j]*E)

        for i in 1:size(E,2)
            X[i,j]= Yj[i,i]
        end
    end
    return normalize_solution!(X)
end
