The package `AlgebraicSolvers.jl` provides tools for solving polynomial equations

To install the package within julia:

```julia
Pkg.clone("https://gitlab.inria.fr/AlgebraicGeometricModeling/AlgebraicSolvers.jl.git")
Pkg.build("AlgebraicSolvers")
```


To use it within julia:

```julia
using AlgebraicSolvers

X = @ring x1 x2 x3
n = length(X)

d = 3
M = monomials(X,0:d)
s = length(M)

P = (2*rand(n,s)-fill(1.0,n,s))*M

Xi = solve_macaulay(P,X)

```

