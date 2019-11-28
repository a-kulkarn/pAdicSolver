
# pAdic solver

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://a-kulkarn.github.io/pAdicSolvers/latest)

Package for solving 0-dimensional padic polynomial systems. Examples for usage are contained in the "/test/" directory.

- Dense polynomial system solver for p-adic polynomial systems.
- Input polynomials are AbstractAlgebra.jl polynomials instead of "DynamicPolynomials". A conversion interface is included, but not particularly well-maintained.


## REQUIRED PACKAGES:
- Dory.jl and dependencies
- MultivariatePolynomials  (not strictly needed, but used for legacy tests)
- DynamicPolynomials       (not strictly needed, but used for legacy tests)

Dory.jl can be found at: https://github.com/a-kulkarn/Dory.git


OPTIONAL PACKAGES:
-- ConvertDynamicPolynomials

ConvertDynamicPolynomials can be found at: https://github.com/a-kulkarn/ConvertDynamicPolynomials.git

This package is a fork of Bernard Mourrain's solver.
https://gitlab.inria.fr/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git

## Installation

Installation is with the standard Julia package manager.

    add https://github.com/a-kulkarn/Dory.git
    add https://github.com/a-kulkarn/pAdicSolver.git


# Documentation

In progress. Presently, there are internal comments explaining functions.
