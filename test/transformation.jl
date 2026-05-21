using CIJ
using Test

@testset "Transformations" begin
    @testset "rot3" begin
        c = CIJ.rot3(CIJ.ol()[1], 180rand(), 180rand(), 180rand())

        @testset "Identity" begin
            @test c ≈ rot3(c, 0, 0, 0)
            @test c ≈ rot3(c, 360, 0, 0)
            @test c ≈ rot3(c, 0, 360, 0)
            @test c ≈ rot3(c, 0, 0, 360)
        end

        rot_angle = 30

        @testset "x1" begin
            c′ = rot3(c, rot_angle, 0, 0)
            @test phase_vels(c, 0, 0).pol ≈ phase_vels(c′, 0, 0).pol + rot_angle atol=0.1
        end

        @testset "x2" begin
            c′ = rot3(c, 0, rot_angle, 0)
            @test phase_vels(c, -90, 0).pol ≈ phase_vels(c′, -90, 0).pol + rot_angle atol=0.1
        end

        @testset "x3" begin
            c′ = rot3(c, 0, 0, rot_angle)
            @test phase_vels(c, 0, 90).pol ≈ phase_vels(c′, -rot_angle, 90).pol atol=0.1
        end

        @testset "Allocations" begin
            @test Base.@allocations(rot3(c, rand(), rand(), rand())) == 0
        end
    end
end