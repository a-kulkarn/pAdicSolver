
const DEFAULT_VALUATION_OF_ZERO = BigInt(2^63-1)

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
#  Simultaneous Eigenvalues.
#
######################################################################################################


@doc Markdown.doc"""
    simultaneous_eigenvalues(inputM :: Array{Array{T,2},1} where T <: FieldElem, ish::Bool, method)
            -> eigenvalue_matrix :: Hecke.Generic.MatElem{T}

(Internal function) Compute the eigenvalues of a list of commuting matrices, given as an array. In the
specific case that the matrices are mult-by-coordinate-variable operators on R/I

INPUTS: M -- list of commuting matrices corresponding to mult-by-xi operators
Outputs: A matrix whose j-th column are the eigenvalues of the j-th matrix in M
"""
function simultaneous_eigenvalues(M::Vector; method=:schur)

    if method isa String
        method = Symbol(method)
    end
    
    if method in [:schur, :qr]
        method = :schur
    end
    
    return simultaneous_eigenvalues(Val(method), M)
end

function simultaneous_eigenvalues(::Val{:power}, M::Vector)
    
    Qp = base_ring(M[1])
    Mg = sum(A*randunit(Qp) for A in M) # non-unit random causes problems

    # TODO: Check well-conditionedness of eigenvalue problem for Mg
    
    #@vprint :padic_solver 2 "Valuations of singular values:"
    #@vprint :padic_solver 2 valuation.(singular_values(M[1]))

    # eigen vectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
    eigen_subspaces = eigspaces(Mg, method = :power)
    
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


function simultaneous_eigenvalues(::Val{:schur}, M::Vector)
    
    Qp = base_ring(M[1])
    Mrand = sum(A * randunit(Qp) for A in M) # non-unit random causes problems

    if length(M) >= 6
        Mg = M[2] + M[3] + M[5] + M[6]
    else
        Mg = Mrand
    end
    
    # eigenvectors of inv(M0)*M[1], which are common eigenvectors of inv(M0)*M[i]
    X, V = Dory.block_schur_form(Mg)
    simple_blocks = filter(x->length(x)==1, Dory.diagonal_block_ranges(X))

    # Allocate and populate
    sol_array = matrix(Qp, fill(zero(Qp), length(simple_blocks), length(M)))
    
    for j = 1:length(M)
        Y = V * M[j] * inv(V)
        for ran in simple_blocks
            for i in ran
                sol_array[i, j] = Y[i, i]
            end
        end
    end
    return sol_array
end


function simultaneous_eigenvalues(::Val{:tropical}, M::Vector)

    # Basic setup
    Qp = base_ring(M[1])
    n = size(M[1], 1)
    base_zero = zero(Qp)

    # Setup containers
    # We also consider `V` to be one of the allocated containers.
    V = Dory.unaliased_identity_matrix(Qp, n)
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
            tol_check = x->!Dory.isapprox_zero(x, valuation_atol = min(norm_val, -5) + min_prec)
            boolean_mats[j] = tol_check.(A)

            if false && j == i
                @info " " norm_val min_prec
                @info " " elt_info.(A)
                @info " " boolean_mats[j]
            end

        end

        # Refine the list of ranges based on common blocks after change of coordinates.
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
end

