
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
#  Normalize solutions(?)
#
######################################################################################################

# Why did Bernard need this step? It seems totally pointless...
# function normalize_solution!(Xi, ish)
#     Sol = Xi
    
#     if !ish
#         i=1            
#         while i <= size(Sol,1)
#             scale_factor = Sol[i,1]
#             if iszero(scale_factor)
#                 println()
#                 println("!-- Ignoring solution at infinity --!")
#                 println()
                
#                 Sol = vcat(Sol[1:(i-1), :], Sol[i+1:size(Sol,1), :])
#             else
#                 Sol[i,:] *= inv(scale_factor)
#                 i+=1
#             end
#         end
#         #else
#         # do nothing otherwise, for now.
#     end
#     return Sol
# end


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
    simultaneous_eigenvalues(inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)
            --> eigenvalue_matrix :: Hecke.Generic.MatElem{T}

(Internal function) Compute the eigenvalues of a list of commuting matrices, given as an array. In the
specific case that the matrices are mult-by-coordinate-variable operators on R/I

INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
"""
function simultaneous_eigenvalues(M::Vector; method=:schur)

    if method in [:schur, :qr]
        return simultaneous_eigenvalues_schur(M)
        
    elseif method in [:tropical]
        return simultaneous_eigenvalues_tropical(M)
        
    elseif method in [:power]
        return simultaneous_eigenvalues_power(M)
    else
        @error "method = $method is not recognized by simultaneous_eigenvalues."
    end
end

function simultaneous_eigenvalues_power(M::Vector)
    
    Qp = base_ring(M[1])
    Mg = sum(A*rand_padic_int(Qp) for A in M) # non-unit random causes problems

    # TODO: Check well-conditionedness of Mg
    
    @vprint :padic_solver 2 "Valuations of singular values:"
    @vprint :padic_solver 2 valuation.(singular_values(M[1]))

    # eigen vectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
    eigen_subspaces  = eigspaces(Mg, method="power")
    
    #println("eigvalues: ", invariant_subspaces.values)
    #println()
    #println("eigspaces: ", length(invariant_subspaces.spaces))


    X = matrix(Qp, fill(zero(Qp), length(eigen_subspaces.spaces), length(M)))
    for j in 1:length(M)
        for i in 1:length(eigen_subspaces.spaces)
            
            V = eigen_subspaces.spaces[i]            
            Y = rectangular_solve(V, M[j] * V)           
            X[i,j] = trace(Y)/Qp(size(Y,2))
        end
    end
    
    return X
end

@doc Markdown.doc"""
    simultaneous_eigenvalues_schur(inputM :: Array{Array{T,2},1} where T <: FieldElem)

(Internal function) Compute the eigenvalues of a list of commuting matrices, given as an array. In the
specific case that the matrices are mult-by-coordinate-variable operators on R/I. The method attempts to
use the Schur form decomposition.

INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
"""

function simultaneous_eigenvalues_schur(M::Vector)
    
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
    return hcat(sol_array...)
end


# Function to extract information from a local field element.
elt_info(x) = (iszero(x), valuation(x), precision(x))


# TODO: Move to Dory?
function simultaneous_eigenvalues_tropical(M::Vector)

    # Basic setup
    Qp = base_ring(M[1])
    n = size(M[1], 1)
    base_zero = zero(Qp)

    # Setup containers
    # We also consider `V` to be one of the allocated containers.
    V = Dory.identity_matrix(Qp, n)
    B = deepcopy(V)
    C = deepcopy(V)
    zero_mat = zero_matrix(Qp, n, n)
    boolean_mats = [iszero.(A) for A in M]
    
    # Enforce good memory structure
    @assert !Dory._has_any_shared_refs(V)

    # Setup structure for recording blocks.
    block_ranges = [1:n]
    
    for i = 2:length(M)
        
        #@info "ROUND $i"

        # We want to perturb the original matrix a little to ensure that the blocks
        # detected in previous iterations stay intact.
        B = add!(B, M[i], zero_mat)
        
        # Set entries in the columns block_start:block_end under the diagonal block to zero.
        for block_ran in block_ranges
            for i = last(block_ran)+1:n
                for j = block_ran
                    B[i,j] = zero!(B[i,j])
                end
            end
        end
        
        #@info "Input of Block Schur Form:" elt_info.(B)
        X, Vi = Dory.block_schur_form(B, shift=tropical_shift, iter_bound=(m,N)->100)
        #@info "Output of Block Schur Form:" elt_info.(X)
        
        #@info "Output of Block Schur Form:" elt_info.(X)
        #@info "Transform:" elt_info.(V)
        #@info "Inverse precision" precision.(inv(V))

        # Update all the matrices
        for j=1:length(M)
            # Compute B = V * M[i] * inv(V)
            M[j] = let
                C = Hecke.mul!(C, Vi, M[j])
                Hecke.mul!(M[j], C, inv(Vi))
            end
        end

        # Update the predicted blocks up to pnumerical tolerance.

        # Broadcast isapprox to create boolean matrices.
        for j=1:length(M)
            A = M[j]
            norm_val = Dory.norm_valuation(A)
            min_prec = minimum(precision.(A))

            # NOTE: We need zero entries to show up as zero bits for the block form check.
            tol_check = x->!Dory.isapprox_zero(x, valuation_atol = min(norm_val, 0) + min_prec)
            boolean_mats[j] = tol_check.(A)

            if false && j == i
                @info " " norm_val min_prec
                @info " " elt_info.(A)
                @info boolean_mats[j]
            end

        end

        # Refine the list of ranges based on common blocks post-change of coordinates.
        all_block_ranges = [Dory.diagonal_block_ranges(A) for A in boolean_mats]
        block_ranges = let
            temp = Vector{UnitRange{Int}}()
            for rangs in all_block_ranges
                for r in rangs
                    if !any(issubset(r, S) for S in temp)
                        filter!(S->!issubset(S, r), temp)
                        push!(temp, r)
                    end
                end
            end
            temp
        end

        # Update V = Vi * V
        C = Hecke.mul!(C, Vi, V)
        V, C = C, V

        #for j = 1:length(M)
            #i == j && @info  "Round $(i),  Matrix $(j):" elt_info.(M[j])
        #end         
    end

    
    ########################################
    # Extract eigenvalues after simultaneous diagonalization.

    # Initialize the solution array.
    # sol_array = Hecke.MatrixSpace(Hecke.QQ, n, length(M))()

    sol_array = Matrix{TropicalInterval}(undef, n, length(M))
    
    for j = 1:length(M)
        for ran in block_ranges
            block = M[j][ran, ran]            
            sing_vals = singular_values(block)

            # NOTE: It might be faster to look at traces of powers of the matrix and use Newton's identities.
            
            #@info valuation.(sing_vals)
            
            # In any particular block, the valuations of the eigenvalues are equal.
            # We need to check if the block has a kernel, as a special default value needs to be assigned.


            if zero(Qp) in sing_vals
                val_of_eigenvalues = fmpq(precision(Qp))
                sol_interval = TropicalInterval(val_of_eigenvalues, Inf)
            else
                sing_val_sizes = (fmpq(valuation(x)) for x in sing_vals)
                val_of_eigenvalues = sum(sing_val_sizes) // length(sing_vals)
                sol_interval = TropicalInterval(val_of_eigenvalues, val_of_eigenvalues)
            end

            # Populate the solution array
            for i in ran
                sol_array[i,j] = sol_interval
            end
        end
    end
    
    return sol_array

    ####################
    # Junk (Hopefully)
    
    # sol_array = Array{Array{Number,1},1}()
    # for j in 1:length(M)

    #     # Put the other matrix into schur form
    #     ycoords = Array{Number,1}()
    #     #Y = V * M[j] * IV
    #     Y = M[j]
    #     block_start_index = 1

    #     #@info  "Matrix $(j):" elt_info.(Y)
        
    #     for i=1:size(X,2)
    #         if (i == size(X,2) || iszero(X[i+1, i]))
                
    #             block_inds = block_start_index:i
    #             block = Y[block_inds, block_inds]

    #             sing_vals = singular_values(block)

    #             #@info valuation.(sing_vals)
                
    #             # In any particular block, the valuations of the eigenvalues are equal.
    #             # We need to check if the block has a kernel, as a special default value needs to be assigned.
    #             if zero(Qp) in sing_vals
    #                 val_of_eigenvalues = DEFAULT_VALUATION_OF_ZERO
    #             else
    #                 sing_val_sizes = [BigInt(valuation(x)) for x in sing_vals]
    #                 val_of_eigenvalues = sum(sing_val_sizes) // length(sing_vals)
    #             end
                
    #             push!(ycoords, [val_of_eigenvalues for r in block_inds]...)
    #             block_start_index = i+1
    #         end                
    #     end
    #     push!(sol_array, ycoords)
    # end

    # return hcat(sol_array...)
end

