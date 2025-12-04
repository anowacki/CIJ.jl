using CIJ
using Test

function components_equal(a, b; atol)
    all(≈(a[i,j], b[i,j]; atol) for i in 1:6 for j in 1:6)
end

@testset "Decomposition" begin
    @testset "B&S olivine" begin
        C = EC(
            [192  66  60   0   0   0
              66 160  56   0   0   0
              60  56 272   0   0   0
               0   0   0  60   0   0
               0   0   0   0  62   0
               0   0   0   0   0  49]
        )
        d = decompose(C)

        @testset "Proportions" begin
            @test d.p_iso ≈ 0.793 atol=0.001
            @test d.p_hex ≈ 0.152 atol=0.001
            @test d.p_tet + d.p_orth ≈ 0.055 atol=0.001
            @test d.p_mono ≈ 0 atol=0.001
            @test d.p_tri ≈ 0 atol=0.001
        end

        @testset "ECs" begin
            @testset "Isotropic" begin                
                @test d.Ciso[1,1] ≈ 194.7 atol=0.1
                @test d.Ciso[4,4] ≈ 63.7 atol=0.1
                @test d.Ciso[1,2] ≈ 67.3 atol=0.1
                @test d.Ciso[1,1] == d.Ciso[2,2] == d.Ciso[3,3]
                @test d.Ciso[4,4] == d.Ciso[5,5] == d.Ciso[6,6]
                @test d.Ciso[1,2] == d.Ciso[1,3] == d.Ciso[2,3]
                @test all(inds -> iszero(d.Ciso[inds...]),
                    [(i,j) for i in 1:6 for j in 1:6 if i !=j && i > 3 && j > 3]
                )
            end

            @testset "Hexagonal" begin                
                @test d.Chex[1,1] ≈ -21.7 atol=0.1
                @test d.Chex[2,2] ≈ -21.7 atol=0.1
                @test d.Chex[3,3] ≈ 77.3 atol=0.1
                @test d.Chex[4,4] ≈ -2.7 atol=0.1
                @test d.Chex[5,5] ≈ -2.7 atol=0.1
                @test d.Chex[6,6] ≈ -11.7 atol=0.1
                @test d.Chex[1,2] ≈ 1.7 atol=0.1
                @test d.Chex[1,3] ≈ -9.3 atol=0.1
                @test d.Chex[2,3] ≈ -9.3 atol=0.1
                @test all(inds -> iszero(d.Chex[inds...]),
                    [(i,j) for i in 1:6 for j in 1:6 if i !=j && i > 3 && j > 3]
                )
            end

            @testset "Tetragonal" begin
                @test d.Ctet[1,1] ≈ 3 atol=0.1
                @test d.Ctet[2,2] ≈ 3 atol=0.1
                @test d.Ctet[1,2] ≈ -3 atol=0.1
                @test d.Ctet[6,6] ≈ -3 atol=0.1
                @test all(inds -> iszero(d.Ctet[inds...]),
                    [
                        (i,j) for i in 1:6 for j in 1:6
                        if (i, j) ∉ [(1,1), (2,2), (1,2), (2,1), (6,6)]
                    ]
                )
            end

            @testset "Orthorhombic" begin
                @test d.Corth[1,1] ≈ 16 atol=0.1
                @test d.Corth[2,2] ≈ -16 atol=0.1
                @test d.Corth[3,3] ≈ 0 atol=0.1
                @test d.Corth[1,2] ≈ 0 atol=0.1
                @test d.Corth[1,3] ≈ 2 atol=0.1
                @test d.Corth[2,3] ≈ -2 atol=0.1
                @test d.Corth[4,4] ≈ -1 atol=0.1
                @test d.Corth[5,5] ≈ 1 atol=0.1
                @test d.Corth[6,6] ≈ 0 atol=0.1
                @test all(inds -> iszero(d.Corth[inds...]),
                    [(i,j) for i in 1:6 for j in 1:6 if i !=j && i > 3 && j > 3]
                )
            end

            @testset "Others" begin
                @test all(iszero, d.Cmono)
                @test all(iszero, d.Ctri)
            end
        end
    end

    @testset "B&S Enstatite" begin
        C = EC([
            225  54  72   0   0   0
             51 214  53   0   0   0
             72  53 178   0   0   0
              0   0   0  78   0   0
              0   0   0   0  82   0
              0   0   0   0   0  76
        ])
        d = decompose(C)

        @testset "Proportions" begin
            @test d.p_iso ≈ 0.908 atol=0.001
            @test d.p_hex ≈ 0.043 atol=0.001
            @test d.p_tet + d.p_orth ≈ 0.049 atol=0.001
        end

        @testset "Isotropic" begin
            @test components_equal(
                d.Ciso,
                [210.2  57.4  57.4     0     0     0
                  57.4 210.2  57.4     0     0     0
                  57.4  57.4 210.2     0     0     0
                     0     0     0  76.4     0     0
                     0     0     0     0  76.4     0
                     0     0     0     0     0  76.4]; atol=0.1
            )
        end

        @testset "Hexagonal" begin
            @test components_equal(
                d.Chex,
                [5.9   0   5.1   0   0   0
                   0 5.9   5.1   0   0   0
                 5.1 5.1 -32.2   0   0   0
                   0   0     0 3.6   0   0
                   0   0     0   0 3.6   0
                   0   0     0   0   0 3.0]; atol=0.1
            )
        end

        @testset "Tetragonal" begin
            @test components_equal(
                d.Ctet,
                [3.4 -3.4  0  0  0    0
                -3.4  3.4  0  0  0    0
                  0     0  0  0  0    0
                  0     0  0  0  0    0
                  0     0  0  0  0    0
                  0     0  0  0  0 -3.4]; atol=0.1
            )
        end

        @testset "Orthorhombic" begin
            @test components_equal(
                d.Corth,
                [5.5    0  9.5   0  0  0
                   0 -5.5 -9.5   0  0  0
                 9.5 -9.5    0   0  0  0
                   0    0    0  -2  0  0
                   0    0    0   0  2  0
                   0    0    0   0  0  0]; atol=0.1
            )
        end
    end
end
