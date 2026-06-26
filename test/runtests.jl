using Test

@testset "All tests" begin
    include("types.jl")
    include("synthesis.jl")
    include("properties.jl")
    include("averaging.jl")
    include("effective_medium.jl")
    include("conversion.jl")
    include("velocities.jl")
    include("io.jl")
end
