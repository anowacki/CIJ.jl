using CIJ, Test

@testset "EC" begin
    @testset "Constructors" begin
        let
            # Check on size
            @test_throws ArgumentError EC([1 2; 3 4])
            @test_throws ArgumentError EC(rand(5,6))
            # No constructor
            @test_throws ArgumentError EC((1.0, 2.0))
            # Default precision
            @test EC(zeros(6,6)) isa EC{Float64}
            @test EC(zeros(Float32, 6, 6)) isa EC{Float32}
            @test EC(fill(1, 6, 6)) isa EC{float(Int)}
            @test EC(Array{Any,2}(rand(6,6))) isa EC{CIJ.DEFAULT_FLOAT}
            # Symmetrisation using upper half
            x = [10i + j for i in 1:6, j in 1:6]
            c = EC(x)
            @test c == [i<j ? 10i+j : 10j+i for i in 1:6, j in 1:6]
            # Warning for asymmetric matrices
            @test_logs (:warn, "input matrix not symmetrical: taking upper half"
                ) EC((x = rand(6,6); x[1,2] = -x[1,2]; x), warn=true)
            # Tuple and array construction
            c′ = EC(((10i + j for j in 1:6 for i in 1:6)...,))
            @test c == c′
            # Construction with another EC uses same underlying MMatrix
            c2 = EC(c)
            @test c == c2
            c[1,1] = rand()
            @test c2[1,1] == c[1,1]
            # Symmetrisation on setindex!
            c[1,2] = rand()
            @test c[2,1] == c[1,2]
            @test zero(EC) == zeros(CIJ.DEFAULT_FLOAT, 6, 6)
            @test zero(EC{CIJ.DEFAULT_FLOAT}) == zero(EC)
            @test zero(c) == zero(EC)
        end
    end

    @testset "Properties" begin
        let c = EC(rand(6,6))
            @test size(c) == (6, 6)
            @test length(c) == 36
        end
    end
end