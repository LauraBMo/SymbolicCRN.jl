push!(LOAD_PATH, "$(homedir())/.julia/dev")

using Test
using SafeTestsets# , TestSetExtensions

@testset verbose = true "SymbolicCRN.jl" begin
    @safetestset "Cones.jl" begin include("Cones.jl") end;

    @safetestset "ConvexParameters.jl" begin include("ConvexParameters.jl") end;

    @safetestset "FindRoots.jl" begin end;

    @safetestset "HomotopyContinuation.jl" begin end;

    @safetestset "NemoUtils.jl" begin include("NemoUtils.jl") end;

    @safetestset "PosNegVertices.jl" begin end;

    @safetestset "Stability.jl" begin end;

    @safetestset "Utils.jl" begin end;

end
