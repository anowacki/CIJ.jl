# Test creation of elastic tensors from parameters

using CIJ, Test

@testset "Synthesis" begin
    @testset "Iso" begin
        let vp=8000, vs=4500, λ=vp^2-2vs^2, μ=vs^2, K=vp^2-4/3*vs^2, G=μ
            @test_throws ArgumentError iso(vp=64e6, lam=23.5e6)
            # Vp, Vs
            c = iso(vp=vp, vs=vs)
            @test c isa EC{CIJ.DEFAULT_FLOAT}
            @test c[1,1] == c[2,2] == c[3,3] == 64_000_000
            @test c[4,4] == c[5,5] == c[6,6] == 20_250_000
            @test c[1,2] == c[1,3] == c[2,3] == 23_500_000
            # λ, μ
            @test iso(lam=λ, mu=μ) == c
            # K, G
            @test iso(K=K, G=G) == c
        end
    end

    @testset "thom_st" begin
        @testset "Invalid input" begin
            @test_throws ArgumentError thom_st(0, 1000, 0, 0, 0)
            @test_throws ArgumentError thom_st(4000, -1000, 0, 0, 0)
            @test_throws ArgumentError thom_st(4000, 2000, -1, 0, 0)
        end

        @testset "Known values" begin
            c = thom_st(4000, 2000, 0.01, -0.01, 0.03)
            @test c[1,1] == c[2,2] == 16_320_000
            @test c[3,3] == 16_000_000
            @test c[4,4] == c[5,5] == 4_000_000
            @test c[6,6] == 3_920_000
            @test c[1,2] == 8_480_000
            @test c[1,3] ≈ 8_393_546.708 atol=0.001
            @test c[2,3] == c[1,3]
            @test c[1,4] == c[1,5] == c[1,6] == c[2,4] == c[2,5] == c[2,6] ==
                c[3,4] == c[3,5] == c[3,6] == c[4,5] == c[4,6] == c[5,6] == 0
        end
    end
end
