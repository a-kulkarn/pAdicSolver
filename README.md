The package `AlgebraicSolvers.jl` provides some tools for solving polynomial equations

To install the package within julia:

```julia
using Pkg
Pkg.clone("https://gitlab.inria.fr/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git")
Pkg.build("AlgebraicSolvers")
```


To use it within julia:

```julia
using AlgebraicSolvers

X = @Ring x1 x2 x3
n = length(X)

d = 3
M = monomials(X,0:d)
s = length(M)

P = randn(n,s)*M

Xi = solve_macaulay(P,X)

```

