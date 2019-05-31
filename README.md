This repository contains two related packages. Examples for usage are contained in the "/test/" directory.


=== Dory.jl ===

Package to extend the functionality of Nemo/Hecke. Notable additions include

basic utilities
-- Allows Julia broadcasting for Nemo matrices.

padic linear algebra:
-- padic qr-factorization.
-- padic singular value decomposition.
-- padically stable solving of linear systems.
-- padically stable hessenburg form.
-- eigenvector solver (power and inverse iterations). [Only implemented for eigenvectors defined over Qp]
-- block schur form.

=== Algebraic solvers ===

This package is a fork of Bernard Mourrain's solver.
https://gitlab.inria.fr/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git

-- Dense polynomial system solver for p-adic polynomial systems.
-- Input changed to accept AbstractAlgebra.jl polynomials instead of "DynamicPolynomials". A conversion interface is included, but not particularly well-maintained.


=== Installation/Documentation ===

-- At the moment, the package does not work with Julia's package manager.
-- Plenty of interal comments explain what is going on, though documentation for the code is missing.