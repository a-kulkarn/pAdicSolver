
module pAdicSolver

  using Markdown
  using Hecke
  using Dory

  # Add a nullspace method for pAdic arrays.
  import Dory: nullspace

  # TODO: check if Singular is installed.
  if true
      import Singular
      include("rauls_change_base_ring.jl") # TODO: Temporary
  end

  include("mindex.jl")
  include("matrix.jl")
  include("pmatrix_util.jl")
  include("macaulay.jl")
  include("test_util.jl")
  #include("newton.jl")
  #include("toric.jl")

  # Requires Singular.jl
  include("solver.jl")

end
