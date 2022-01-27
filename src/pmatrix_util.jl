
const DEFAULT_VALUATION_OF_ZERO = BigInt(2^63-1)

######################################################################################################
#
#  Legacy code
#
######################################################################################################

"""
    nullspace(A::Array{padic,2})

Compute the nullspace of an array of padic numbers, viewed as a matrix.
"""
function nullspace(A::Array{padic,2})
    M = matrix(A)
    return Hecke.nullspace(M)[2].entries
end


struct QRPadicArrayPivoted
    Q::Array{padic,2}
    R::Array{padic,2}
    p::Array{Int64,1}
end

######################################################################################################
#
#  Tropical shifting choices in QR iteration
#
######################################################################################################

# Function to feed into Dory.block_schur_form for optimal performance
function tropical_shift(B)
    return zero(base_ring(B))
end


######################################################################################################
#
#  Normalized Simultaneous Eigenvalues.
#
######################################################################################################


@doc Markdown.doc"""
    normalized_simultaneous_eigenvalues(inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)
            --> eigenvalue_matrix :: Hecke.Generic.MatElem{T}

(Internal function) Compute the eigenvalues of a list of commuting matrices, given as an array. In the
specific case that the matrices are mult-by-coordinate-variable operators on R/I

INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
"""
function normalized_simultaneous_eigenvalues(inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)

    if method == "schur" || method == "qr" || method == "tropical"
        return nse_schur(inputM, ish, method)
    end
    
    M = [matrix(A) for A in inputM]
    Qp = base_ring(M[1])
    Mrand = sum(A*rand_padic_int(Qp) for A in M) # non-unit random causes problems

    @vprint :padic_solver 2 "Valuations of singular values:"
    @vprint :padic_solver 2 valuation.(singular_values(M[1]))

    # We will assume that M[1] is well-conditioned. for now.
    
    # The rectangular solve step is enough to kill off any helpful data mod p.
    @vtime :padic_solver 2 I0 = rectangular_solve(M[1], identity_matrix(Mrand.base_ring,size(M[1],1)))
    @vtime :padic_solver 2 Mg = I0 * Mrand

    # eigen vectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
    eigen_subspaces  = eigspaces(Mg, method=method)
    
    #println("eigvalues: ", invariant_subspaces.values)
    #println()
    #println("eigspaces: ", length(invariant_subspaces.spaces))


    X = matrix(Qp, fill(zero(Qp), length(eigen_subspaces.spaces), length(M)))
    for j in 1:length(M)
        for i in 1:length(eigen_subspaces.spaces)
            
            V = eigen_subspaces.spaces[i]            
            Y = rectangular_solve(V,I0*M[j]*V)           
            X[i,j] = trace(Y)/Qp(size(Y,2))
        end
    end

    function normalize_solution!(Xi, ish)
        Sol = Xi
        
        if !ish
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

@doc Markdown.doc"""
    nse_schur(inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)
            --> eigenvalue_matrix :: Hecke.Generic.MatElem{T}

(Internal function) Compute the eigenvalues of a list of commuting matrices, given as an array. In the
specific case that the matrices are mult-by-coordinate-variable operators on R/I. The method attempts to
use the Schur form decomposition.

INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
"""

function nse_schur(inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)
    
    M = [matrix(A) for A in inputM]
    Qp = base_ring(M[1])
    Mrand = sum(A * rand_padic_int(Qp) for A in M) # non-unit random causes problems

    @vprint :padic_solver 2 "Valuations of singular values: "
    @vprint :padic_solver 2 valuation.(singular_values(M[1]))

    # We will assume that M[1] is well-conditioned. for now.
    
    # The rectangular solve step is enough to kill off any helpful data mod p.
    @vtime :padic_solver 2 I0 = inv(M[1])

    
    #println("eigvalues: ", invariant_subspaces.values)
    #println()
    #println("eigspaces: ", length(invariant_subspaces.spaces))

    if method != "tropical"

        Mg = I0 * (M[2] + M[3] + M[5] + M[6])
        
        # eigenvectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
        X, V = Dory.block_schur_form(Mg)
        
        sol_array = Array{Array{padic,1},1}()
        for j in 1:length(M)

            # Put the other matrix into schur form
            Y = inv(V)*(I0*M[j])*V
            
            for i=1:size(X,2)
                if (i == 1 && iszero(X[i, i+1])) ||
                    (i==size(X, 2) && iszero(X[i-1, i])) ||
                    (iszero(X[i, i-1]) && iszero(X[i+1, i]))
                    
                    push!(sol_array[j], Y[i,i])
                end
            end
        end
        return normalize_solution(Xi, ish)
    end

    if method == "tropical"

        Mg = I0 * M[2]
        
        
        # eigen vectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
        X, V = Dory.block_schur_form(Mg, shift=tropical_shift)
        
        sol_array = Array{Array{Number,1},1}()
        #display(valuation.(X))
        
        for j in 1:length(M)

            # Put the other matrix into schur form
            ycoords = Array{Number,1}()
            Y = V * (I0 * M[j]) * inv(V)
            block_start_index = 1

            #display(valuation.(Y))
            
            for i=1:size(X,2)
                if (i == size(X,2) || iszero(X[i+1, i]))
                    
                    block_inds = block_start_index:i
                    block = Y[block_inds, block_inds]

                    sing_vals = singular_values(block)

                    @info valuation.(sing_vals)
                    
                    # In any particular block, the valuations of the eigenvalues are equal.
                    # We need to check if the block has a kernel, as a special default value needs to be assigned.
                    if zero(Qp) in sing_vals
                        val_of_eigenvalues = DEFAULT_VALUATION_OF_ZERO
                    else
                        sing_val_sizes = [BigInt(valuation(x)) for x in sing_vals]
                        val_of_eigenvalues = sum(sing_val_sizes) // length(sing_vals)
                    end
                    
                    push!(ycoords, [val_of_eigenvalues for r in block_inds]...)
                    block_start_index = i+1
                end                
            end
            push!(sol_array, ycoords)
        end

        if ish
            error("Tropical mode not implemented for projective systems.")
        end
        #return hcat( sol_array[2:length(sol_array)]...)
        return hcat(sol_array...)
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


function tropical_nse(inputM)
    
    M = [matrix(A) for A in inputM]
    Qp = base_ring(M[1])

    # The rectangular solve step is enough to kill off any helpful data mod p.
    I0 = inv(M[1])
    
    ########################################
    # Simultaneous diagonalization.

    X, V = Dory.block_schur_form(I0 * M[2], shift=tropical_shift)

    # Setup containers
    B = X
    C = deepcopy(X)
    zero_mat = zero_matrix(Qp, size(B,1), size(B, 2))
    
    for A in M[3:length(M)]
        B = mul!(B, I0, A)

        # Compute B = V * B * inv(B)
        B = let
            C = mul!(C, V, B)
            mul!(B, C, inv(V))
        end
        
        X, Vi = Dory.block_schur_form(B, shift=tropical_shift)

        # Compute V = Vi * V
        B = mul!(B, Vi, V)
        V = add!(V, B, zero_mat)
    end

    ## Cache the inverse
    IV = inv(V)
    
    ########################################
    # Extract eigenvalues after simultaneous diagonalization.
    sol_array = Array{Array{Number,1},1}()
    for j in 1:length(M)

        # Put the other matrix into schur form
        ycoords = Array{Number,1}()
        Y = V * (I0 * M[j]) * IV
        block_start_index = 1

        for i=1:size(X,2)
            if (i == size(X,2) || iszero(X[i+1, i]))
                
                block_inds = block_start_index:i
                block = Y[block_inds, block_inds]

                sing_vals = singular_values(block)

                #@info valuation.(sing_vals)
                
                # In any particular block, the valuations of the eigenvalues are equal.
                # We need to check if the block has a kernel, as a special default value needs to be assigned.
                if zero(Qp) in sing_vals
                    val_of_eigenvalues = DEFAULT_VALUATION_OF_ZERO
                else
                    sing_val_sizes = [BigInt(valuation(x)) for x in sing_vals]
                    val_of_eigenvalues = sum(sing_val_sizes) // length(sing_vals)
                end
                
                push!(ycoords, [val_of_eigenvalues for r in block_inds]...)
                block_start_index = i+1
            end                
        end
        push!(sol_array, ycoords)
    end

    return hcat(sol_array...)
end

