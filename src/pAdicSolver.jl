
module pAdicSolver

using Markdown
using Hecke
using Dory

# Add a nullspace method for pAdic arrays.
import Dory: nullspace

# We might choose to add some functionality depending on Singular later.
# import Singular

export TropicalInterval
export macaulay_mat, solve_macaulay, solve_affine_system, solve_projective_system
export solve_affine_groebner_system, solve_projective_groebner_system, solve_system
export abs_error, rel_error, forward_error, solcmp, random_square_system


#include("mindex.jl")
include("matrix.jl")
include("TropicalInterval.jl")
include("simultaneous_eigen.jl")
include("macaulay.jl")
include("test_util.jl")
#include("newton.jl")
#include("toric.jl")

###############################################################
#
# Init.
#
###############################################################

function __init__()
    # Attach the necessary verbose scopes to Hecke's global dicitionary.
    Hecke.add_verbose_scope(:padic_solver)
    return
end

end
