using CIJ
using Test

datapath = joinpath(@__DIR__, "data")

@testset "IO" begin
    @testset "read" begin
        file = joinpath(datapath, "test.ecs")

        @testset "file" begin
            C, rho = CIJ.read(file)
            @test C == [11 12 13 14 15 16
                        12 22 23 24 25 26
                        13 23 33 34 35 36
                        14 24 34 44 45 46
                        15 25 35 45 55 56
                        16 26 36 46 56 66]
            @test rho == 1000
        end

        @testset "IOBuffer" begin
            @test CIJ.read(file) == CIJ.read(IOBuffer(read(file)))
        end

        @testset "Errors" begin
            @test_throws ErrorException CIJ.read("file_which_doesn't_exist.ecs")
            @test_throws ErrorException CIJ.read("nonexistent_directory/file.ecs")
            @test_throws ErrorException CIJ.read(IOBuffer("afd as ffs saf a"))
            @test_throws ErrorException CIJ.read(IOBuffer("1 7 1000"))
        end
    end

    @testset "write" begin
        mktempdir() do dir
            file = joinpath(dir, "test.ecs")

            @testset "Type $T" for T in (EC, Matrix)
                C, rho = CIJ.ol()
                C = T(C)
                @testset "file" begin
                    CIJ.write(file, C, rho, "X")
                    lines = readlines(file)
                    @test length(lines) == 23
                    @test lines[end] == "# X"
                    C′, rho′ = CIJ.read(file)
                    @test C ≈ C′/rho′
                    @test rho == rho′
                end

                @testset "IOBuffer" begin
                    io = IOBuffer()
                    CIJ.write(io, C, rho, "X")
                    seekstart(io)
                    lines = readlines(io)
                    @test length(lines) == 23
                    @test lines[end] == "# X"
                    seekstart(io)
                    C′, rho′ = CIJ.read(io)
                    @test C ≈ C′/rho′
                    @test rho == rho′
                end
            end
        end

        @testset "Default comment" begin
            io = IOBuffer()
            CIJ.write(io, CIJ.ol()...)
            lines = split(chomp(String(take!(io))), '\n')
            @test length(lines) == 23
            @test occursin(r"# Saved on .* using CIJ.jl on Julia .*",
                lines[end])
        end

        @testset "No comment" begin
            io = IOBuffer()
            CIJ.write(io, CIJ.ol()..., "")
            lines = split(chomp(String(take!(io))), '\n')
            @test length(lines) == 22
        end

        @testset "Unstable ECs" begin
            C = rand(5,5)
            rho = 1
            @test_throws ArgumentError CIJ.write(IOBuffer(), C, rho)
            C = rand(6,6)
            @test_logs (:warn,
                "elastic constants are not in the right form.  May be asymmetric or unstable."
                ) CIJ.write(IOBuffer(), C, rho)
        end
    end
end
