
module pAdicSolver

  using MultivariatePolynomials
  using DynamicPolynomials
  using Hecke
  using Dory

  #degree = DynamicPolynomials.maxdegree

  export coeftype
  coeftype(::Type{Polynomial{C, T}}) where {C, T} = T
  coeftype(p::Polynomial{C, T}) where {C, T} = T

  Base.one(X::Vector{PolyVar{true}}) = monomials(X,0)[1]

  include("convert_dynamic_polynomials.jl")
  include("mindex.jl")
  include("matrix.jl")
  include("pmatrix_util.jl")
  include("macaulay.jl")
  include("test_util.jl")
  #include("newton.jl")
  #include("toric.jl")



end
