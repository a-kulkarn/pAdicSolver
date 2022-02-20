export macaulay_mat, solve_macaulay, solve_affine_system, solve_projective_system
export solve_affine_groebner_system, solve_projective_groebner_system, solve_system


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
    @assert all(ishomogeneous, P)
    _solve_system_method_dispatch(P, true; kwds...)
end

@doc Markdown.doc"""
    _padic_solver_groebner_hook(P; ordering=ordering(parent(P[1])))

To access the `method = :groebner` functionality in the solver, one must import this function
and add a method which outputs a Groebner basis with respect to the specified ordering. This
method must accept the ordering keyword and get the ordering of the parent polynomial ring from
the input `P`.

EXAMPLE:
```
using pAdicSolver
using Oscar       # Has groebner basis functionality

import pAdicSolver._solver_groebner_hook as F

function F(P; ordering = ordering(parent(P[1])))
    return groebner_basis(ideal(P), ordering=ordering)
end

```
Then the following will work:
```
#Blah
```
"""
function _solver_groebner_hook(P; ordering=ordering(parent(P[1])))
    error("Groebner hook not added. See documentation for `_solver_groebner_hook`.")
end

function _solve_system_method_dispatch(P, is_homogeneous; method = :truncated_normal_form, kwds...)

    if method in [:truncated_normal_form, :tnf, :macaulay]
        return _solver_engine(P, is_homogeneous; method = :tnf, kwds...)

    elseif method in [:given_groebner, :given_GB, :givenGB]
        return _solver_engine(P, is_homogeneous; method = :given_GB, kwds...)
        
    elseif method == :groebner
        use_order = :ordering in keys(kwds) ? kwds[:ordering] : ordering(parent(P[1]))
        
        # Use the user supplied groebner basis function, or error if they have not given one.
        gb = _solver_groebner_hook(P, ordering = use_order)

        # 2. Call the Solver Engine with the groebner basis.
        return _solver_engine(gb, is_homogeneous; method=:given_GB, ordering=use_order, kwds...)
        
    elseif method isa String
        ArgumentError("Keyword 'method' must be of type Symbol.")
    else
        ArgumentError("method = $method  is not recognized by the solver.")
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

function solve_projective_groebner_system(P; kwds...)
    return _solver_engine(P, true; method = :given_GB, kwds...)
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
    ish = all(ishomogeneous, P)
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
    Mentries = _multiplication_matrices(Val(method), P, is_homogeneous; kwds...)
    M = [matrix(A) for A in Mentries]

    # Normalize the multiplication matrices.
    I0 = try
        # The first multiplication matrix is usually well-conditioned.
        inv(M[1])
    catch e
        !isa(e, Dory.InconsistentSystemError) && throw(e)
        rand_combo = sum(Dory.randint(parent(M[1][1,1])) * M[j] for j=1:length(M))

        # Try to invert one more time with a random linear combination, and if it fails,
        # there is a lack of precision to solve the polynomial system.
        try
            inv(rand_combo)
        catch e
            !isa(e, Dory.InconsistentSystemError) && throw(e)
            throw(Dory.InsufficientPrecisionError())
        end
    end

    # Use the first operator to cancel out the weird change of basis from the truncated
    # normal form approach.
    M = [I0 * A for A in M]
    
    # Simultaneous diagonalization.
    @vtime :padic_solver Xi = simultaneous_eigenvalues(M, method = eigenvector_method)

    # In the affine system, the distinguished_homogeneizing monomial (i.e, "1" for that case) does 
    # not correspond to a coordinate.
    if is_homogeneous return Xi else return Xi[:, 2:size(Xi,2)] end    
end



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


    missing_monomials = let
        if is_homogeneous
            np1 = length(X)
        else
            np1 = length(X)+1
        end
        !(length(L) == binomial(rho + np1 - 1, rho))
    end

    # TODO: In theory this causes a bug.
    @assert !missing_monomials
    
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
                                  ordering=ordering(parent(P[1])))

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

    if !(K isa Dory.DiscreteValuedField)
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
#     ish = !any(!ishomogeneous, P)

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
        
#         # TODO: Use AbstractAlgebra.diagonal_matrix to diagonal join.
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

######################################################################################################
#
#  Polynomial methods
#
######################################################################################################

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
#  Computing a numerically stable basis to compute the multiplication-by-xi operators.
#
######################################################################################################

struct QRPadicArrayPivoted{T}
    Q::Array{T,2}
    R::Array{T,2}
    p::Array{Int64,1}
end


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
    Fpinv = invperm(F.p)

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
    iwasawa_step(N :: Hecke.Generic.MatSpaceElem{<:Dory.DiscreteValuedFieldElem}, L0)

Return the QR-factorization object (For PN₀P' = QR, return <inv(P)Q, R, P'>)
together with  Nr = N*inv(inv(P)Q)^T.
"""
function iwasawa_step(N :: Hecke.Generic.MatSpaceElem{<:Dory.DiscreteValuedFieldElem}, L0)

    sorted_monomial_rows = sort(collect(values(L0)))
    Nkbase = N[sorted_monomial_rows, :]

    F = padic_qr(transpose(Nkbase), col_pivot=Val(true))
    Qinv = Dory.inv_unit_lower_triangular(F.Q)
    Fpinv= invperm(F.p)

    X = Qinv[Fpinv,:]
    Farr = QRPadicArrayPivoted((F.Q.entries)[Fpinv,:], F.R.entries, F.q)

    return Farr, N*X
end
