# Test conversion between representations

using Test, CIJ
import StaticArrays

@testset "Conversion" begin
    @testset "Stiffness-Compliance" begin
        let c = CIJ.ol()[1], atol=1e-12, c′ = copy(c), ca = Array(c), ca′ = copy(ca)
            # ECs
            s = C2S(c)
            @test s[1,1] ≈  1.16918e-8 atol=atol
            @test s[2,2] ≈  2.03838e-8 atol=atol
            @test s[3,3] ≈  1.70592e-8 atol=atol
            @test s[1,2] ≈ -3.04176e-9 atol=atol
            @test s[1,3] ≈ -2.58468e-9 atol=atol
            @test s[2,3] ≈ -5.77166e-9 atol=atol
            for i in 4:6
                @test s[i,i] ≈ 1/c[i,i]
            end
            @test all(isapprox.([vec(s[1:3,4:6]); s[4,5:6]; s[5,6]], 0.0, atol=atol))
            @test C2S(S2C(c)) ≈ c
            @test C2S!(c′) ≈ s
            @test typeof(c) == typeof(s)
            # Arrays
            ca = Array(c)
            sa = C2S(ca)
            @test typeof(ca) == typeof(sa)
            @test S2C(C2S(ca)) ≈ ca
            @test sa ≈ s
            @test C2S!(ca′) == sa
            @test S2C!(ca′) ≈ c
        end
    end
end

@testset "Cij/Cijkl" begin
    let vp = rand(), vs = 1.7*vp, ϵ = γ = δ = rand(), (α, β, γ) = 360 .* rand(3),
            c = CIJ.rot3(CIJ.thom(vp, vs, ϵ, γ, δ), α, β, γ), ca = Array(c)
        c4 = cijkl(c)
        @test c4 isa StaticArrays.MArray
        @test size(c4) == (3, 3, 3, 3)
        @test c4[1,1,1,1] == c[1,1]
        @test c4[2,2,2,2] == c[2,2]
        @test c4[3,3,3,3] == c[3,3]
        @test c4[2,3,2,3] == c[4,4]
        @test c4[1,3,1,3] == c[5,5]
        @test c4[1,2,1,2] == c[6,6]
        @test c4[1,1,2,2] == c[1,2]
        @test c4[1,1,3,3] == c[1,3]
        @test c4[2,2,3,3] == c[2,3]
        @test c4[1,1,1,3] == c[1,5]
        @test cij(cijkl(c)) ≈ c
        c4a = cijkl(ca)
        @test c4a isa StaticArrays.MArray
        @test cij(cijkl(ca)) ≈ ca
    end
end
