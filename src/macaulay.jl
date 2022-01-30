export macaulay_mat, solve_macaulay
using LinearAlgebra

function is_not_homogeneous(p)
    L = [total_degree(t) for t in Hecke.terms(p)]
    return maximum(L) != minimum(L)
end

export macaulay_mat


######################################################################################################
#
# Creation of Linear algebra structures.
#
######################################################################################################

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
    for p in P
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
    for p in P
        for m in mult_monomials[total_degree(p)]

            srow = sparse_row(R, [monomial_dict[mon] for mon in monomials(m*p)],
                              collect(coefficients(p)))
            push!(macaulay_matrix, srow)
        end
    end

    return macaulay_matrix, monomial_dict
end


######################################################################################################
#
#  Solver (Interface)
#
######################################################################################################

@doc Markdown.doc"""
    solve_affine_groebner_system(P :: Vector{Hecke.Generic.MPolyElem{T}} where T <: Hecke.RingElem;
                   rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                   eigenvector_method  :: String ="power",
                   test_mode  :: Bool =false )

Given a system of polynomials `P` that is a Groebner basis, with leading monomials `LP`, return
the solutions to the polynomial system `P`.
"""

function solve_affine_system(P; kwds...)
    _solve_system_method_dispatch(P, false; kwds...)
end

function solve_projective_system(P; kwds...)
    @assert !any(is_not_homogeneous, P)
    _solve_system_method_dispatch(P, true; kwds...)
end


function _solve_system_method_dispatch(P, is_homogeneous; method = :truncated_normal_form, kwds...)

    if method in [:truncated_normal_form, :tnf, :macaulay]
        return _solver_engine(P, is_homogeneous; method = :tnf, kwds...)

    elseif method in [:given_groebner, :given_GB]
        return _solver_engine(P, is_homogeneous; method = :given_GB, kwds...)
        
    elseif method == :groebner
        # 1. Compute the groebner basis over the exact field.
        I = ideal(parent(P[1]), P)
        #gb = Oscar.groebner_basis(I)

        # Fetch the leading monomials somehow.
        @error "method == :groebner is not implemented."
        LP = "TODO"
        
        # 2. Call the Solver Engine with the groebner basis and leading monomials.
        return _solver_engine(P, is_homogeneous, LP; kwds...)
        
    elseif method isa String
        @error "Method must be of type Symbol."
    else
        @error "Specified method is not implemented"
    end
end

@doc Markdown.doc"""
    solve_affine_groebner_system(P :: Vector{Hecke.Generic.MPolyElem{T}} where T <: Hecke.RingElem;
                   rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                   eigenvector_method  :: String ="power",
                   test_mode  :: Bool =false )

Given a system of polynomials `P` that is a Groebner basis, with leading monomials `LP`, return
the solutions to the polynomial system `P`.
"""
function solve_affine_groebner_system(P; kwds...)
    return _solver_engine(P, false; method = :given_GB, kwds...)
end

function solve_projective_groebner_system(P; method = :given_GB, kwds...)
    return _solver_engine(P, true; kwds...)
end


@doc Markdown.doc"""
    solve_system(P :: Vector{Hecke.Generic.MPolyElem{T}} where T <: Hecke.RingElem;
                   rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                   eigenvector_method  :: String ="power",
                   test_mode  :: Bool =false )

Alias for `solve_affine_system`.
"""
solve_system = solve_affine_system


######################################################################################################
#
#  Legacy code
#
######################################################################################################

@doc Markdown.doc"""
    solve_macaulay(P :: Vector{Hecke.Generic.MPolyElem{T}} where T <: Hecke.RingElem;
                   rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                   eigenvector_method  :: String ="power",
                   test_mode  :: Bool =false )

Solve a 0-dimensional system of polynomial equations. (Presently, only over Qp.) More precisely,
compute the values of x in Qp^n such that 

    all([iszero(p(x)) for p in P]) == true

The options specify strategy parameters.

#-------------------

INPUTS:
- P        -- polynomial system, a Vector of AbstractAlgebra polynomials.
- rho      -- monomial degree of the system. Default is the macaulay degree.
- groebner -- Boolean indicating if the polynomial ring is already a groebner basis 
              w.r.t the ambient ring ordering.
- eigenvector_method -- Strategy to solve for eigenvectors. Default is power iteration.

"""
function solve_macaulay(P; kwds...)
    ish = !any(is_not_homogeneous, P)
    @vprint :padic_solver 2 "-- Homogeneity $(ish)"
    if !ish return solve_affine_system(P; kwds...) else return solve_projective_system(P; kwds...) end
end


######################################################################################################
#
#    Main solver functionality
#
#    calls to:
#        # macaulay_mat
#        # iwasawa_step
#        # mult_matrix
#        # eigdiag
#
######################################################################################################


function _solver_engine(P, is_homogeneous; method = :tnf, eigenvector_method = :power, kwds...)
    
    # NOTE: The first "multiplication matrix" either corrsponds to the operator [x0]B, or to [1]B,
    #       where B is some change of basis matrix.
    M = _multiplication_matrices(Val(method), P, is_homogeneous; kwds...)

    # Normalize the multiplication matrices.
    # We assume that the first multiplication matrix is well-conditioned.
    I0 = inv(matrix(M[1]))

    # Use the first operator to cancel out the weird change of basis from the truncated
    # normal form approach.
    M = [I0 * matrix(A) for A in M]
    
    # Simultaneous diagonalization.
    @vtime :padic_solver Xi = simultaneous_eigenvalues(M, method = eigenvector_method)

    # In the affine system, the distinguished_homogeneizing monomial (i.e, "1" for that case) does 
    # not correspond to a coordinate.
    if is_homogeneous return Xi else return Xi[:, 2:size(Xi,2)] end    
end

# function _solver_engine(P, is_homogeneous; eigenvector_method = "power", kwds...)
    
#     # NOTE: The first "multiplication matrix" either corrsponds to the operator [x0]B, or to [1]B,
#     #       where B is some change of basis matrix.
#     #
#     # Dispatch on the method argument.
#     M = _multiplication_matrices(Val(:given_GB), P, is_homogeneous, leading_mons_of_P; kwds...)

#     # Apply the Eigenvector method.
#     @vtime :padic_solver Xi = normalized_simultaneous_eigenvalues(M, is_homogeneous, eigenvector_method)
    
#     # In the affine system, the distinguished_homogeneizing monomial (i.e, "1" for that case) does 
#     # not correspond to a coordinate.
    
#     if is_homogeneous return Xi else return Xi[:, 2:size(Xi,2)] end
    
# end


@doc Markdown.doc"""
    _multiplication_matrices(method::Val{:tnf}, P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem}, 1}, is_homogeneous;
                            rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                            test_mode :: Bool = false )

Determines the multiplication-by-coordinate functions on the polynomial ring modulo the ideal generated by `P`.
The options specify strategy parameters.

#-------------------

INPUTS:
- P        -- polynomial system, a Vector of AbstractAlgebra polynomials.
- rho      -- monomial degree of the system. Default is the macaulay degree.
"""
function _multiplication_matrices(method::Val{:tnf}, P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem}, 1}, is_homogeneous::Bool;
                        rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
                        test_mode :: Bool = false)

    # This solve function could be made to work with polynomials with FlintRR coefficients
    # as well, though this requires managing the type dispatch a bit and remodeling the
    # old DynamicPolynomials based subfunctions.
    
    the_ring = parent(P[1])
    X = gens(the_ring)
    
    @vprint :padic_solver 2 "\n-- Degrees $(map(p->total_degree(p),P)))\n"
    
    t0 = time()
    R, L = macaulay_mat(P, X, rho, is_homogeneous)
    L0 = monomials_divisible_by_x0(L, is_homogeneous)

    msg = "-- Macaulay matrix $size(R,1) x $size(R,2) $(time()-t0) (s)\n"
    @vprint :padic_solver msg
    t0 = time()
    
    @vtime :padic_solver N = nullspace(R)[2]
    
    @vprint :padic_solver "-- -- rank of Macaulay matrix $(size(R,2) - size(N,2))\n"
    msg = string("-- Null space ", size(N,1), " x ", size(N,2), "   ", time()-t0, " (s)\n"); t0 = time()
    @vprint :padic_solver msg

    # The idea of the QR step is two-fold:
    # 1: Choose a well-conditioned *monomial* basis for the algebra from a given spanning 
    #    set (here, IdL0).
    #    This is accomplis_homogeneoused by pivoting. The columns corresponding to F.p[1:size(N,2)] form
    #    a well-conditioned submatrix.
    #
    # 2: Present the algebra in Q-coordinates, which has many zeroes. Note that the choice of
    #    coordinates is not important in the final step, when the eigenvalues are calulated.
    #
    F, Nr = iwasawa_step(N, L0)
    B = permute_and_divide_by_x0(L0, F, is_homogeneous)

    @vprint :padic_solver "-- Qr basis $length(B)   $(time()-t0) (s)\n"
    t0 = time()

    @vtime :padic_solver M = mult_matrices(B, X, Nr, L, is_homogeneous)
    t0 = time()

    if test_mode
        println("TESTING MODE: Computation incomplete. Returning partial result.")
        return M, F, B, N, Nr, R, IdL0, Idx
    end

    return M
end

function _multiplication_matrices(method::Val{:given_GB}, P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem}, 1}, is_homogeneous::Bool;
                                  ordering=ordering(P))

    the_ring = parent(P[1])
    K = base_ring(the_ring)
    X = gens(the_ring)
    insert!(X, 1, the_ring(1))
    
    # Ensure that the specified ordering is used for reduction/kbase.
    R, vars = PolynomialRing(K, length(gens(the_ring)), ordering=ordering)
    PR = [change_parent(R, g) for g in P]

    # Determine the basis for the quotient algebra.
    BR = kbase(PR, is_homogeneous)

    # Construct the multiplication operators as polynomial elements.
    XR = [change_parent(R, v) for v in X]
    xi_operators = [[rem(v*b, PR) for b in BR] for v in XR]

    if !(K isa Hecke.NonArchLocalField)
        # TODO: Figure out what the user should do in this case...
        K = PadicField(29,30)
    end

    #@info xi_operators
    
    M = [matrix(K, [coeff(g, b) for b in BR, g in op]) for op in xi_operators]
    M = [m.entries for m in M]

    return M 
end

# TODO: This function should replace the previous version if it tests well.
# function _multiplication_matrices_II(method::Val{:tnf}, P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem}, 1}, is_homogeneous;
#                         rho :: Integer =  sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
#                         test_mode :: Bool = false )
#     # Start the clock.
#     t0 = time()

#     # This solve function could be made to work with polynomials with FlintRR coefficients
#     # as well, though this requires managing the type dispatch a bit and remodeling the
#     # old DynamicPolynomials based subfunctions.
    
#     the_ring = parent(P[1])
#     X = gens(the_ring)

#     # Compute the truncated normal form
#     #
#     # Return the map N, together with
#     # the sets of monomials (V,B), with N|B the identity.

#     N, L, L0 = truncated_normal_form_map(P, is_homogeneous, rho=rho, verbose=verbose)

#     @vprint :padic_solver "-- -- rank of Macaulay matrix $(size(R,2) - size(N,2))"
#     @vprint :padic_solver "-- Null space $size(N,1) x $size(N,2)  $(time()-t0) (s)"
    
#     t0 = time()
#     Nr, B = truncated_normal_form_section(N, L, L0, is_homogeneous)
   
#     @vprint :padic_solver "-- Qr basis $length(B)   $(time()-t0) (s)"

#     @vtime :padic_solver M = mult_matrices(B, X, Nr, L, is_homogeneous)
#     return M
# end

######################################################################################################
#
# Truncated Normal Form logic
#
######################################################################################################

# function truncated_normal_form_map(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1}, is_homogeneous;
#                                    rho::Integer= sum(total_degree(P[i])-1 for i in 1:length(P)) + 1,
#                                    verbose=false)

#     the_ring = parent(P[1])
#     K = base_ring(the_ring)
    
#     X = gens(the_ring)
#     ish = !any(is_not_homogeneous, P)

#     t0 = time()
#     R, L = macaulay_mat(P, X, rho, ish)

#     @vprint :padic_solver "-- Macaulay matrix $size(R,1) x $size(R,2) $(time()-t0) (s)"
#     t0 = time()
    
#     @vtime :padic_solver N = nullspace(R)[2]

#     # Detect if there are missing monomials, and update N accordingly.
#     # (Missing monomials automatically lie in the kernel).
#     missing_monomials = let
#         if ish
#             np1 = length(X)
#         else
#             np1 = length(X)+1
#         end
#         !(length(L) == binomial(rho + np1 - 1, rho))
#     end

    
#     if missing_monomials
#         V = ish ? monomials_of_degree(X, rho) : monomials_of_degree(X, 0:rho)

#         missing_count = 0
#         missing_mons = Array{typeof(X[1]), 1}()
#         for m in V
#             if !(m in keys(L))
#                 push!(missing_mons, m)
#                 missing_count += 1
#             end
#         end

#         # Update indices in L
#         for m in keys(L)
#             L[m] += missing_count
#         end

#         # Add the new labels
#         for i=1:length(missing_mons)
#             L[missing_mons[i]] = i
#         end
        
#         # TODO: When AbstractAlgebra finally implements diagonal joining, use that.
#         #
#         # Update matrices via a diagonal join.
#         newN = zero_matrix(K, missing_count + size(N,1), missing_count + size(N,2))
#         for i=1:missing_count
#             newN[i,i] = K(1)
#         end
#         for i=1:size(N,1)
#             for j=1:size(N,2)
#                 newN[i+missing_count, j+missing_count] = N[i,j]
#             end
#         end        
#         N = newN
#     end

#     # The monomials for the set such that [xi]B \in V.
#     L0 = monomials_divisible_by_x0(L, ish) 

#     @vprint :padic_solver "-- -- rank of Macaulay matrix $(size(R,2) - size(N,2))"
#     @vprint :padic_solver "-- Null space $size(N,1) x $size(N,2)  $(time()-t0) (s)"

#     return N, L, L0
# end

# function truncated_normal_form_section(N, L, L0, ish)
#     # The idea of the QR step is two-fold:
#     # 1: Choose a well-conditioned *monomial* basis for the algebra from a given spanning 
#     #    set (here, IdL0).
#     #    This is accomplished by pivoting. The columns corresponding to F.p[1:size(N,2)] form
#     #    a well-conditioned submatrix.
#     #
#     # 2: Present the algebra in Q-coordinates, which has many zeroes. Note that the choice of
#     #    coordinates is not important in the final step, when the eigenvalues are calulated.
#     #

#     # TODO: The result of this output should be mathematically meaningful.
#     Nr, B = iwasawa_step(N, L0, ish)
    
#     # Turn B into an algebra basis.
#     if ish
#         B = Dict(Dory.divexact(m, gens(parent(m))[1])=>i for (m,i) in B)
#     end
 
#     return Nr, B # Not actually a section, but whatever.
# end


######################################################################################################
#
#  Singular Dependency:
#  Functions for picking monomials given a Groebner basis.
#
######################################################################################################

# function (R::FlintPadicField)(a::Singular.n_Q)
#     return R(FlintQQ(a))
# end

# function kbase_gens_from_GB(P::Array{<:Hecke.Generic.MPolyElem{<:Hecke.FieldElem},1})

#     the_ring = parent(P[1])
#     @assert base_ring(the_ring) == FlintQQ
    
#     sing_R, sing_vars = Singular.PolynomialRing(Singular.QQ,
#                                                 ["x$i" for i=1:nvars(the_ring)],
#                                                 ordering=ordering(the_ring))
    
#     singular_B = kbase_gens_from_GB(map(f-> change_base_ring(f, Singular.QQ, sing_R), P))

#     return map(f-> change_base_ring(f, Hecke.FlintQQ, the_ring), singular_B)
# end

# function kbase_gens_from_GB(P::Array{<:Singular.spoly{<:Hecke.FieldElem},1})
#     I = Singular.Ideal(parent(P[1]), P)
#     I.isGB = true
#     return gens(Singular.kbase(I))
# end

# function rem(f::Singular.spoly{<:Hecke.FieldElem},
#              P::Array{<:Singular.spoly{<:Hecke.FieldElem},1})

#     I = Singular.Ideal(parent(P[1]), P)
#     I.isGB = true
#     return Singular.reduce(f, I)
# end

# TODO: Move to Dory
import Base.rem
function rem(f::Hecke.Generic.MPolyElem{T}, P::Array{<:Hecke.Generic.MPolyElem{T},1}) where T
    return divrem(f,P)[2]
end

function rem(f::Hecke.Generic.PolyElem{T}, P::Array{<:Hecke.Generic.PolyElem{T},1}) where T
    return divrem(f,P)[2]
end

function rem(f::Hecke.Generic.MPolyElem{T}, P::Hecke.Generic.MPolyElem{T}) where T
    return divrem(f,P)[2]
end

function rem(f::Hecke.Generic.PolyElem{T}, P::Hecke.Generic.PolyElem{T}) where T
    return divrem(f,P)[2]
end


function change_parent(P, f)
    K  = base_ring(P)
    fP = MPolyBuildCtx(P)
    for (coeff, exp_vec) = zip(coefficients(f), exponent_vectors(f))
        push_term!(fP, K(coeff), exp_vec)
    end
    return finish(fP)
end 


@doc Markdown.doc"""
    kbase(P, ish)

Given a Groebner basis `P` for an ideal, and the corresponding leading monomials `LP`, return a basis
for the quotient algebra `k[x1, ..., xn]/P` as a k-vector space.

It is assumed that `P` is a Groebner basis with respect to the ordering of the parent.
"""
function kbase(P, ish)

    # TODO: This code is not optimal since it generates more monomials than is necessary.
    #       It is much more efficient to enumerate by degree and build up the basis via
    #       multiplying by monomials.

    LP = leading_monomial.(P)

    D = maximum(total_degree(f) for f in P)
    degree_range = ish ? D : 0:D
    mons = monomials_of_degree(parent(P[1]), degree_range)
    
    divisible_by_LP_elt = f->any(iszero(rem(f, lp)) for lp in LP)
    return filter(!divisible_by_LP_elt, mons)
end


######################################################################################################
# 
#  Oscar Dependency
#
######################################################################################################

# @doc Markdown.doc"""
#     Given a list of polynomials `P`, which is a Groebner basis for some ideal,
#     as well as the list of leading monomials LP, construct the generators for the quotient
#     ring `K[x]/P` as a K-vector space 
# """
# function kbase_from_GB(P, LP)

#     isempty(P) && @error "Input list of polynomials must not be empty."
    
#     R = parent(P[1])

#     for m in monomials_of_degree
        
#     end
        
#     @error "Not implemented."
# end
    
# function kbase_from_GB(I)
#     Oscar.singular_assure(I)
#     sI = I.gens.S
#     sI.isGB = true

#     return sing_kbase = Oscar.Singular.kbase(sI)
# end


######################################################################################################
# 
#  Computing a numerically stable basis to compute the multiplication-by-xi operators.
#
######################################################################################################

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
        
@doc Markdown.doc"""
    iwasawa_step(N :: Array{T,2} where T <: Number, L0)
    iwasawa_step(N :: Array{padic,2} , IdL0)
    iwasawa_step(N :: Hecke.Generic.MatSpaceElem{padic} , L0)

    Return the QR-factorization object (For PN₀P' = QR, return <inv(P)Q, R, P'>)
    together with  Nr = N*inv(inv(P)Q)^T.
"""

function iwasawa_step(N :: Array{T,2} where T <: Number, L0)
    F = qr(Array(transpose(N[IdL0,:])) , Val(true))
    return F, N*F.Q
end

function iwasawa_step(N :: Array{padic,2} , IdL0)
    
    F = padic_qr(transpose(matrix(N[IdL0,:])) , col_pivot=Val(true))
    Qinv = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv= Dory.inverse_permutation(F.p)

    X = Qinv[Fpinv,:].entries    
    Farr = QRPadicArrayPivoted((F.Q.entries)[Fpinv,:], F.R.entries, F.q)
    
    return Farr, N*X
end

function iwasawa_step(N :: Hecke.Generic.MatSpaceElem{padic}, L0)

    F = padic_qr(transpose(N[ sort(collect(values(L0))),:]), col_pivot=Val(true))
    Qinv = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv= Dory.inverse_permutation(F.p)

    X = Qinv[Fpinv,:]
    Farr = QRPadicArrayPivoted((F.Q.entries)[Fpinv,:], F.R.entries, F.q)
    
    return Farr, N*X
end
