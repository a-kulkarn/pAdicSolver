
using Test, Hecke, Dory, pAdicSolver

@testset verbose = true "pAdicSolver" begin
    include("basic_plane_test.jl")
    include("Seibenundzwanzig.jl")
    include("random_system_test.jl")
    include("weird_system.jl")
    include("insufficient_precision_test.jl")
    include("laurent_series_test.jl")
end
