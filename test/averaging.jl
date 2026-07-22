using CIJ
using Test

@testset "Averaging" begin
    @testset "VRH" begin
        C1, ρ1 = iso(vp=2000, vs=1000), 1800
        C2, ρ2 = iso(vp=1500, vs=100), 1200

        @testset "Same tensors" begin
            @test VRH(0.1, C1, ρ1, 0.2, C1, ρ1)[1] ≈ C1
            @test VRH(0.1, C1, ρ1, 0.2, C1, ρ1)[2] ≈ ρ1
        end

        @testset "Different tensors" begin
            C,  ρ = VRH(0.5, C1, ρ1, 1.5, C2, ρ2)

            @test ρ ≈ (1800*0.5 + 1200*1.5)/2
            @test C ≈ EC([
                 2.51792e6  2.24713e6  2.24713e6  0.0        0.0        0.0
                 2.24713e6  2.51792e6  2.24713e6  0.0        0.0        0.0
                 2.24713e6  2.24713e6  2.51792e6  0.0        0.0        0.0
                 0.0        0.0        0.0        1.35395e5  0.0        0.0
                 0.0        0.0        0.0        0.0        1.35395e5  0.0
                 0.0        0.0        0.0        0.0        0.0        1.35395e5
            ]) atol=10
            @test is_iso(C; atol=√eps(eltype(C)))
        end
    end
end