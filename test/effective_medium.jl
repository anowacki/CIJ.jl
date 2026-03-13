using CIJ
using Test

@testset "Effective media" begin
    @testset "pitl" begin
        C, ρ = CIJ.pitl(0.1, 2000, 1000, 1500, 0.01, 1500, 100, 1000)
        @test ρ ≈ (0.1*1500 + 0.01*1000)/0.11
        @test all(
            .≈(
                C,
                [ 3.87762e6  2.00137e6  1.95105e6      0.0      0.0       0.0
                  2.00137e6  3.87762e6  1.95105e6      0.0      0.0       0.0
                  1.95105e6  1.95105e6  3.58224e6      0.0      0.0       0.0
                  0.0        0.0        0.0        70898.4      0.0       0.0
                  0.0        0.0        0.0            0.0  70898.4       0.0
                  0.0        0.0        0.0            0.0      0.0  938125.0];
                  atol=100
            )
        )
    end

    @testset "tandon_and_weng" begin
        @testset "Water-filled cracks" begin
            C, ρ = tandon_and_weng(3240, 1800, 3800, 0.01, 0.01, 1500, 0, 1000)
            @test ρ ≈ 0.99*3800 + 0.01*1000 atol=0.001
            @test all(
                .≈(
                    C, 
                    [ 9.18418e6  3.57056e6  3.57056e6  0.0        0.0        0.0
                      3.57056e6  1.03229e7  3.86088e6  0.0        0.0        0.0
                      3.57056e6  3.86088e6  1.03229e7  0.0        0.0        0.0
                      0.0        0.0        0.0        3.23099e6  0.0        0.0
                      0.0        0.0        0.0        0.0        2.10493e6  0.0
                      0.0        0.0        0.0        0.0        0.0        2.10493e6];
                    atol=100
                )
            )
        end
    end

    @testset "thomsen_cracks" begin
        @testset "Invalid arguments" begin
            @test_throws ArgumentError thomsen_cracks(4000, 2000, -0.1, 0.01)
            @test_throws ArgumentError thomsen_cracks(4000, 2000, 1.1, 0.01)
        end

        @test all(
            .≈(
                thomsen_cracks(4850.0, 2770.0, 0.2, 0.06),
                [
                 2.50279e7  7.59037e6  7.20205e6  0.0       0.0       0.0
                 7.59037e6  2.50279e7  7.20205e6  0.0       0.0       0.0
                 7.20205e6  7.20205e6  2.35225e7  0.0       0.0       0.0
                 0.0        0.0        0.0        7.6729e6  0.0       0.0
                 0.0        0.0        0.0        0.0       7.6729e6  0.0
                 0.0        0.0        0.0        0.0       0.0       8.71878e6
                ],
                atol=100
            )
        )
    end
end
