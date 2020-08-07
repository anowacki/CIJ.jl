using Test

@testset "All tests" begin
    include("types.jl")
    include("synthesis.jl")
    include("properties.jl")
    include("conversion.jl")
end
