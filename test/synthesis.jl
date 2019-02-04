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
end