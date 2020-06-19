export macaulay_mat, solve_macaulay

using LinearAlgebra

function is_not_homogeneous(p)
    L = [total_degree(t) for t in Hecke.terms(p)]
    return maximum(L) != minimum(L)
end

export macaulay_mat

## Present issues.
# Performance is not good, and garbage collector runs way too much.
# bizzare reversal should be made more robust.
# homogeneity is not handled.
#
@doc Markdown.doc"""
    macaulay_mat(P::Array{Hecke.Generic.MPoly{T},1},
                      X::Array{Hecke.Generic.MPoly{T},1}, rho, ish) where T <: Hecke.RingElem

Constructs the sparse macaulay matrix defined by the polynomials `P` and degree bound `rho`. The argument `X`
is a list of variables of the ambient polynomial ring used to construct the multiplier monomials. 
"""
function macaulay_mat(P::Array{Hecke.Generic.MPoly{T},1},
                      X::Array{Hecke.Generic.MPoly{T},1}, rho, ish) where T <: Hecke.RingElem

    degrees = unique!(map(p->total_degree(p),P))
    monomial_set    = Set{Hecke.Generic.MPoly{T}}()
    mult_monomials  = Array{Array{Hecke.Generic.MPoly{T}}}(undef, maximum(degrees))
    
    for d in degrees
        if ish
            mult_monomials[d] = monomials_of_degree(X, rho-d)
        else
            mult_monomials[d] = monomials_of_degree(X, 0:rho-d)
        end
    end
    @time for p in P
        for m in mult_monomials[total_degree(p)]
            push!(monomial_set, monomials(m*p)...)
        end
    end

    # The method "isless" is defined in AbstractAlgebra. By default Julia will use this to sort.
    monomial_set = collect(monomial_set)
    sort!(monomial_set, rev=true)
    monomial_dict = Dict(monomial_set[i]=>i for i=1:length(monomial_set))

    # Create sparse rows for each m*p, with m a mulitplier monomial and p a polynomial.
    R = base_ring(parent(P[1]))    
    macaulay_matrix = sparse_matrix(R)
    @time for p in P
        for m in mult_monomials[total_degree(p)]

            srow = sparse_row(R, [monomial_dict[mon] for mon in monomials(m*p)],
                                collect(coeffs(p)))
            push!(macaulay_matrix, srow)
        end
    end

    return macaulay_matrix, monomial_dict
end

## ***************************************************************************************    
# Main solver function

    # calls to:
    # macaulay_mat
    # iwasawa_step
    # mult_matrix
    # eigdiag
## ***************************************************************************************

@doc Markdown.doc"""
    solve_macaulay(P :: Vector{Hecke.Generic.MPolyElem{T}} where T <: Hecke.RingElem;
                   rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                   eigenvector_method  :: String ="power",
                   test_mode  :: Bool =false )

Solve a 0-dimensional system of polynomial equations. (Presently, only over Qp.) More precisely,
compute the values of x in Qp^n such that 

    all([ iszero(p(x)) for p in P ]) == true

The options specify strategy parameters.

#-------------------

INPUTS:
- P        -- polynomial system, a Vector of AbstractAlgebra polynomials.
- rho      -- monomial degree of the system. Default is the macaulay degree.
- groebner -- Boolean indicating if the polynomial ring is already a groebner basis 
              w.r.t the ambient ring ordering.
- eigenvector_method -- Strategy to solve for eigenvectors. Default is power iteration.

"""
# Should have a verbose parameter.

function solve_macaulay(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1};
                        rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                        eigenvector_method :: String = "power",
                        test_mode :: Bool = false )

    # This solve function could be made to work with polynomials with FlintRR coefficients
    # as well, though this requires managing the type dispatch a bit and remodeling the
    # old DynamicPolynomials based subfunctions.
    
    the_ring = parent(P[1])
    X = gens(the_ring)
    
    println()
    println("\n-- Degrees ", map(p->total_degree(p),P))
    
    ish = !any(is_not_homogeneous, P)
    println("-- Homogeneity ", ish)

    t0 = time()

    R, L = macaulay_mat(P, X, rho, ish)           
    
    L0 = monomials_divisible_by_x0(L, ish)
    
    println("-- Macaulay matrix ", size(R,1),"x",size(R,2),  "   ",
            time()-t0, "(s)"); t0 = time()
    
    @time N = nullspace(R)[2]
    
    println("-- -- rank of Macaulay matrix ", size(R,2) - size(N,2))
    println("-- Null space ",size(N,1),"x",size(N,2), "   ",time()-t0, "(s)"); t0 = time()

    # The idea of the QR step is two-fold:
    # 1: Choose a well-conditioned *monomial* basis for the algebra from a given spanning 
    #    set (here, IdL0).
    #    This is accomplished by pivoting. The columns corresponding to F.p[1:size(N,2)] form
    #    a well-conditioned submatrix.
    #
    # 2: Present the algebra in Q-coordinates, which has many zeroes. Note that the choice of
    #    coordinates is not important in the final step, when the eigenvalues are calulated.
    #
    F, Nr = iwasawa_step(N, L0)
    B = permute_and_divide_by_x0(L0, F, ish)

    @info "" B
    @info "" L

    println("-- Qr basis ",  length(B), "   ",time()-t0, "(s)"); t0 = time()

    
    M = mult_matrices(B, X, Nr, L, ish)
    println("-- Mult matrices ",time()-t0, "(s)"); t0 = time()


    if test_mode
        println("TESTING MODE: Computation incomplete. Returning partial result.")
        return M, F, B, N, Nr, R, IdL0, Idx
    end

    Xi = normalized_simultaneous_eigenvalues(M, ish, eigenvector_method)
    println("-- Eigen diag", "   ", time()-t0, "(s)"); t0 = time()

    # In the affine system, the distinguished monomial (i.e, "1" for that case) does 
    # not correspond to a coordinate.
    if ish return Xi else return  Xi[:,2:size(Xi,2)] end
end

##############
# Improvement to the solve_macaulay function.


function solve_macaulay_II(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1};
                           rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                           eigenvector_method :: String = "power",
                           verbose::Bool = false,
                           test_mode :: Bool = false )

    # Start the clock.
    t0 = time()

    # This solve function could be made to work with polynomials with FlintRR coefficients
    # as well, though this requires managing the type dispatch a bit and remodeling the
    # old DynamicPolynomials based subfunctions.
    
    the_ring = parent(P[1])
    X = gens(the_ring)
    ish = !any(is_not_homogeneous, P)


    # Compute the truncated normal form
    #
    # Return the map N, together with
    # the sets of monomials (V,B), with N|B the identity.

    N, L, L0 = truncated_normal_form_map(P, rho=rho, verbose=verbose)

    # TODO: Need to figure out what to do with missing monomials.
    
    Nr, B = truncated_normal_form_section(N, L, L0, ish)
   
    verbose && println("-- Qr basis ",  length(B), "   ", time()-t0, "(s)");
    t0 = time()

    
    M = mult_matrices(B, X, Nr, L, ish)
    verbose && println("-- Mult matrices ",time()-t0, "(s)");
    t0 = time()

    if test_mode
        println("TESTING MODE: Computation incomplete. Returning partial result.")
        return M, F, B, N, Nr, R, IdL0, Idx
    end

    # Eigenvector step.
    Xi = normalized_simultaneous_eigenvalues(M, ish, eigenvector_method)

    verbose && println("-- Eigen diag", "   ", time()-t0, "(s)");
    t0 = time()

    # In the affine system, the distinguished monomial (i.e, "1" for that case) does 
    # not correspond to a coordinate.
    if ish return Xi else return  Xi[:,2:size(Xi,2)] end

end
    
##############

function truncated_normal_form_map(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1};
                                   rho::Integer= sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                                   verbose=false)

    the_ring = parent(P[1])
    K = base_ring(the_ring)
    
    X = gens(the_ring)
    ish = !any(is_not_homogeneous, P)

    t0 = time()
    R, L = macaulay_mat(P, X, rho, ish)
    
    verbose && println("-- Macaulay matrix $(size(R,1)) x $(size(R,2))   $(time()-t0) (s)");
    t0 = time()
    
    if verbose
        @time N = nullspace(R)[2]
    else
        N = nullspace(R)[2]
    end

    # Detect if there are missing monomials, and update N accordingly.
    # (Missing monomials automatically lie in the kernel).
    missing_monomials = let
        if ish
            np1 = length(X)
        else
            np1 = length(X)+1
        end
        !(length(L) == binomial(rho + np1 - 1, rho))
    end

    
    if missing_monomials
        V = ish ? monomials_of_degree(X, rho) : monomials_of_degree(X, 0:rho)

        missing_count = 0
        missing_mons = Array{typeof(X[1]), 1}()
        for m in V
            if !(m in keys(L))
                push!(missing_mons, m)
                missing_count += 1
            end
        end

        # Update indices in L
        for m in keys(L)
            L[m] += missing_count
        end

        # Add the new labels
        for i=1:length(missing_mons)
            L[missing_mons[i]] = i
        end
        
        # TODO: When AbstractAlgebra finally implements diagonal joining, use that.
        #
        # Update matrices via a diagonal join.
        newN = zero_matrix(K, missing_count + size(N,1), missing_count + size(N,2))
        for i=1:missing_count
            newN[i,i] = K(1)
        end
        for i=1:size(N,1)
            for j=1:size(N,2)
                newN[i+missing_count, j+missing_count] = N[i,j]
            end
        end        
        N = newN
    end

    # The monomials for the set such that [xi]B \in V.
    L0 = monomials_divisible_by_x0(L, ish) 
    
    verbose && println("-- -- rank of Macaulay matrix ", size(R,2) - size(N,2))
    verbose && println("-- Null space $(size(N,1)) x $(size(N,2))   $(time()-t0) (s)");
    t0 = time()

    return N, L, L0
end

function truncated_normal_form_section(N, L, L0, ish)
    # The idea of the QR step is two-fold:
    # 1: Choose a well-conditioned *monomial* basis for the algebra from a given spanning 
    #    set (here, IdL0).
    #    This is accomplished by pivoting. The columns corresponding to F.p[1:size(N,2)] form
    #    a well-conditioned submatrix.
    #
    # 2: Present the algebra in Q-coordinates, which has many zeroes. Note that the choice of
    #    coordinates is not important in the final step, when the eigenvalues are calulated.
    #

    # TODO: The result of this output should be mathematically meaningful.
    Nr, B = iwasawa_step(N, L0, ish)
    
    # Turn B into an algebra basis.
    if ish
        B = Dict(Dory.divexact(m, gens(parent(m))[1])=>i for (m,i) in B)
    end
 
    return Nr, B # Not actually a section, but whatever.
end

##############

@doc Markdown.doc"""
    iwasawa_step(N :: Hecke.Generic.MatSpaceElem{padic} , L0)

    Return the QR-factorization object (For PN_0P' = QR, return <inv(P)Q, R, P'>)
    together with  Nr = N*inv(inv(P)Q)^T.
"""
function iwasawa_step(N :: Hecke.Generic.MatSpaceElem{padic}, L0, ish)

    """
    New Goal for function.

    Given a matrix N and a subset of rows specified by L0, return

    N' -- A change of coordinates of N in the codomain.
    B' -- A subset of L0 specifying a stable square submatrix.
    """

    sorted_column_labels = sort(collect(values(L0)))
    
    F = padic_qr(transpose(N[sorted_column_labels,:]), col_pivot=Val(true))
    Qinv  = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv = Dory.inverse_permutation(F.p)

    # X = transpose(Qinv[Fpinv,:])
    X = transpose(Qinv)[Fpinv,:]
    
    #Farr = QRPadicArrayPivoted((F.Q.entries)[Fpinv,:], F.R.entries, F.q)

    
    # Next, extract the algebra basis.
    B = Dict()
    m = size(F.Q,1) # The dimension of the quotient algebra.

    # Extract the column to monomial correspondence.    
    key_array = first.(sort(collect(L0), by=x->x[2]))

    for i in 1:m
        push!(B,  key_array[F.q[i]]=>i)
        # should test if the diag. coeff. is not small
    end

    Nr = N*X

    #test_rows = sorted_column_labels[F.q[1:m]]
    #@info "" test_rows
    #@info "" valuation.(Array(Nr[test_rows, :]))
    
    return N*X, B
end


## Depreciate these...
function iwasawa_step(N :: Array{T,2} where T <: Number, L0)
    F = qr(Array(transpose(N[IdL0,:])) , Val(true))
    return F, N*F.Q
end

function iwasawa_step(N :: Array{padic,2} , IdL0)
    
    F = padic_qr( transpose(matrix(N[IdL0,:])) , col_pivot=Val(true))
    Qinv = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv= Dory.inverse_permutation(F.p)

    X = Qinv[Fpinv,:].entries    
    Farr = QRPadicArrayPivoted( (F.Q.entries)[Fpinv,:], F.R.entries, F.q)
    
    return Farr, N*X
end

function (R::FlintPadicField)(a::Singular.n_Q)
    return R(FlintQQ(a))
end

# function kbase_gens_from_GB(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1})

#     the_ring = parent(P[1])
#     @assert base_ring(the_ring) == FlintQQ
    
#     sing_R,sing_vars = Singular.PolynomialRing(Singular.QQ,
#                                                ["x$i" for i=1:nvars(the_ring)],
#                                                ordering=ordering(the_ring))
    
#     singular_B = kbase_gens_from_GB(map(f-> rauls_change_base_ring(f, Singular.QQ, sing_R), P))

#     return map(f-> rauls_change_base_ring(f, Hecke.FlintQQ, the_ring), singular_B)
# end

# function kbase_gens_from_GB(P::Array{<:Singular.spoly{<:Hecke.FieldElem},1})
#     I = Singular.Ideal(parent(P[1]), P)
#     I.isGB = true
#     return gens(Singular.kbase(I))
# end


# import Base.rem
# function rem(f::Hecke.Generic.MPolyElem{<:Hecke.FieldElem},
#              P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1})
#     return divrem(f,P)[2]
# end

# function rem(f::Singular.spoly{<:T where T},
#              P::Array{<:Singular.spoly{<:T where T},1})

#     I = Singular.Ideal(parent(P[1]), P)
#     I.isGB = true
#     return Singular.reduce(f, I)
# end

