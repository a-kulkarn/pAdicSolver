
=== pAdic solver ===

Package for solving 0-dimensional padic polynomial systems. Examples for usage are contained in the "/test/" directory.

-- Dense polynomial system solver for p-adic polynomial systems.
-- Input polynomials are AbstractAlgebra.jl polynomials instead of "DynamicPolynomials". A conversion interface is included, but not particularly well-maintained.


REQUIRED PACKAGES:
-- Dory.jl and dependencies
-- MultivariatePolynomials  (not strictly needed, but used for legacy tests)
-- DynamicPolynomials       (not strictly needed, but used for legacy tests)

Dory.jl can be found at: https://github.com/a-kulkarn/Dory.git


This package is a fork of Bernard Mourrain's solver.
https://gitlab.inria.fr/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git


=== Installation  ===

Installation is with the standard Julia package manager.

    add https://github.com/a-kulkarn/Dory.git
    add https://github.com/a-kulkarn/pAdicSolver.git

-- At the moment, the package does not work with Julia's package manager.

=== Documentation ===

In progress. Presently, there are internal comments explaining functions.
