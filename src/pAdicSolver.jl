
module pAdicSolver

using Markdown
using Hecke
using Dory

# Add a nullspace method for pAdic arrays.
import Dory: nullspace

# TODO: check if Singular is installed.
if true
    import Singular
end

include("mindex.jl")
include("matrix.jl")
include("TropicalInterval.jl")
include("pmatrix_util.jl")
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
