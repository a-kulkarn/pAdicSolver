
module pAdicSolver

  using Hecke
  using Dory

  # Add a nullspace method for pAdic arrays.
  import Dory: nullspace

  include("mindex.jl")
  include("matrix.jl")
  include("pmatrix_util.jl")
  include("macaulay.jl")
  include("test_util.jl")
  #include("newton.jl")
  #include("toric.jl")

end
