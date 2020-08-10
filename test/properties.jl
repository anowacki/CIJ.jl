# Test calculations of properties of tensors

using Test, CIJ
using LinearAlgebra: Symmetric

@testset "Properties" begin
    @testset "Is iso" begin
        let ciso = CIJ.iso(vp=8000, vs=4500), caniso = CIJ.thom(8000, 4500, 0.1, 0.1, 0.1)
            @test is_iso(ciso)
            ciso[1,1] = nextfloat(ciso[1,1])
            @test !isapprox(ciso[1,1], ciso[2,2], atol=eps(eltype(ciso)))
            @test !is_iso(ciso)
            @test is_iso(ciso, atol=√eps(eltype(ciso)))
            @test !is_iso(caniso)
        end
    end

    @testset "Is stable" begin
        let c = CIJ.iso(vp=8000, vs=4500), c′ = CIJ.iso(vp=8000, vs=0)
            @test is_stable(c)
            @test is_stable(c.data)
            @test is_stable(Array(c))
            @test is_stable(Array(c.data))
            @test !is_stable(c′)
            @test !is_stable(c′.data)
            @test !is_stable(Array(c′))
            @test !is_stable(Array(c′.data))

            c[1,1] = -c[1,1]
            @test !is_stable(c)
            c[1,1] = c[2,2]
            c[1,3] = c[1,2] = c[2,3] = c[1,1]
            @test !is_stable(c)
        end
        @testset "Catch wrong exception" begin
            let c = CIJ.ol()[1]
                @test_logs (:warn, "matrix not symmetrical: using upper half only"
                    ) is_stable(rand(6,6))
                @test_throws ArgumentError is_stable(rand(5,5))
            end
        end
    end

    @testset "Is 6×6" begin
        @test CIJ.is_6x6(CIJ.ol()[1])
        @test CIJ.is_6x6(rand(6, 6))
        @test !CIJ.is_6x6(rand(5, 5))
    end

    @testset "Is symm" begin
        @test CIJ.is_symm(CIJ.ol()[1])
        @test CIJ.is_symm(Symmetric(rand(6, 6)))
        @test !CIJ.is_symm((x = Array(CIJ.ol()[1]); x[1,2] = -x[1,2]; x))
    end
end
