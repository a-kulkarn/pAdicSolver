
using Test, Hecke, Dory, pAdicSolver

@testset verbose = true "pAdicSolver" begin
    include("basic_plane_test.jl")
    include("Seibenundzwanzig.jl")
    include("projective_tests.jl")
    include("weird_system.jl")
    include("low_precision_handling.jl")
    include("laurent_series_tests.jl")
end
