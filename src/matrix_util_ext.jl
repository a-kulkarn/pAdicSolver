
## Make some generic extensions to the matrix utilities in Nemo/Hecke

using Hecke, LinearAlgebra
import Hecke: Generic.Mat, Generic.MatElem, nmod_mat, Generic.charpoly

# I don't actually want LinearAlgebra as a dependency, but I do want to mimick the syntax
# to provide familiarity for the user.
#import LinearAlgebra: eigen, Eigen, eigvals, eigvecs

##############################################################################################
#                                                                                            #
#                             Basic interface                                                #
#                                                                                            #
##############################################################################################

# broadcast over the matrix and return a new matrix.
# note that the function must return an element x such that parent(x) is a ring
#
# NOTE: Sadly, this does not hook into Julia's dot-broadcast syntax. In theory it is possible
# to broadcast and return the same data type, but I don't know how to do that.
# 
function Base.broadcast(f, A::Hecke.Generic.Mat{T} where T)
    return matrix( parent(f(A[1,1])), f.(A.entries))
end

# In theory this works, but the matrix never ends up with the right shape.
function Base.iterate(A::Hecke.Generic.Mat{T}, state=1) where T
    return iterate(A.entries, state)
end

function Base.collect(A::Hecke.Generic.Mat{T}, state=1) where T
    return A.entries
end

function Hecke.matrix(A::Array{T,2} where T <: Hecke.Generic.Mat{S} where S)
    @assert reduce(==, [parent(x) for x in A]) 
    return matrix(parent(A[1,1], A))
end
                  
##############################################################################################
#                                                                                            #
#                          Cosmetic override to nullspace                                    #
#                                                                                            #
##############################################################################################


## Make things a little more consistent with the other Julia types
function my_nullspace(A::nmod_mat)
    nu,N = nullspace(A)
    return nu, nu==0 ? matrix(A.base_ring, fill(0,size(A,2),0)) : N[:,1:nu]
end


##############################################################################################
#                                                                                            #
#                          Generic Eigenvalue/Eigenvector functions                          #
#                                                                                            #
##############################################################################################


struct MyEigen{T}
    values::Array{T,1}
    vectors::Hecke.Generic.MatElem{T}
end

# Returns an eigen factorization structure like the default LinearAlgebra.eigen function.
#
# Fails when A mod p has no eigenvalues, because list has eltype Any
"""
eigen(A::nmod_mat)

Computes the Eigenvalue decomposition of A. Requires factorization of polynomials implemented
over the base ring. 
"""
function eigen(A::Hecke.Generic.MatElem{T}) where T
    R,_ = PolynomialRing(A.base_ring)
    g = charpoly(R, A)
    rts = roots(g)
    if isempty(rts)
        error("Not implemented if no roots of char. poly. over the finite field")
    end
    
    Imat = identity_matrix(A.base_ring, size(A,1))

    eig_list = [ (r,my_nullspace(A-r*Imat)...) for r in rts]

    eig_vals = vcat([fill( ep[1] , ep[2]) for ep in eig_list]...)
    eig_vecs = hcat([ ep[3] for ep in eig_list]...)    
    return MyEigen(eig_vals, eig_vecs)
end

# See usual eigvecs
import LinearAlgebra.eigvecs
function eigvecs(A::Hecke.Generic.MatElem{T}) where T
    return eigen(A).vectors    
end

# See usual eigvals
function eigvals(A::Hecke.Generic.MatElem{T}) where T
    return eigen(A).values
end

"""
charpoly(A::nmod_mat)

Returns the characteristic polynomial of A.
"""
function charpoly(A::Hecke.Generic.MatElem{T}) where T
    R,_ = PolynomialRing(A.base_ring)
    return charpoly(R, A)
end


# Needs to be more robust. Also applied to the situation A is square but not of rank 1.

# a slightly generalized version of solve
# WARNING: does not check if the top block is non-singular
function rectangular_solve(M::Hecke.MatElem{T}, b::Hecke.MatElem{T}) where T

    if rows(M) < cols(M)
        error("Not implemented when rows(M) < cols(M)")
    end
    # Extract top nxn block
    A = M[1:cols(M),:]

    
    x = solve(A,b[1:cols(M),:])
    
    if iszero(M*x - b)
        return x
    else
        error("Linear system does not have a solution")
    end
end


