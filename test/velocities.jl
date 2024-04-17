using CIJ, Test

function are_approx(a::Tuple, b::Tuple; atol=1e-6)
    length(a) == length(b) || return false
    all(≈(aa, bb; atol=atol) for (aa, bb) in zip(a, b))
end

@testset "Velocities" begin
    @testset "incaz2cart" begin
        @test CIJ.incaz2cart(0, 0) ≈ [1, 0, 0]
        @test CIJ.incaz2cart(0, 90) ≈ [0, -1, 0]
        @test CIJ.incaz2cart(0, -90) ≈ [0, 1, 0]
        @test CIJ.incaz2cart(0, 45) ≈ [√2/2, -√2/2, 0]
        @test CIJ.incaz2cart(90, 0) ≈ [0, 0, 1]
        @test CIJ.incaz2cart(-45, 0) ≈ [√2/2, 0, -√2/2]
    end

    @testset "cart2incaz" begin
        @test are_approx(CIJ.cart2incaz(1, 0, 0), (0, 0))
        @test are_approx(CIJ.cart2incaz(0, 10, 0), (0, -90))
        @test are_approx(CIJ.cart2incaz(0, 0, -100), (-90, 0))
        @test are_approx(CIJ.cart2incaz(√2/2, √2/2, 1), (45, -45))
    end

    @testset "incaz2up" begin
        @test CIJ.incaz2up(0, 0) ≈ [0, 0, 1]
        @test CIJ.incaz2up(0, -90) ≈ [0, 0, 1]
        @test CIJ.incaz2up(45, 90) ≈ [0, √2/2, √2/2]
        @test CIJ.incaz2up(-45, -90) ≈ [0, √2/2, √2/2]
    end

    @testset "phase_vels" begin
        @testset "iso" begin
            let vp = 8000, vs = 4400, c = CIJ.iso(vp=8000, vs=4400),
                    azis = rand(0:360, 4), incs = rand(-90:90, 4)
                for (azi, inc) in zip(azis, incs)
                    v = CIJ.phase_vels(c, azi, inc)
                    @test v.vp ≈ vp
                    @test v.vs1 ≈ vs
                    @test v.vs2 ≈ vs
                    @test v.avs ≈ 0 atol=√eps()
                    @test v.xp ≈ CIJ.incaz2cart(inc, azi)
                end
            end
        end

        @testset "VTI" begin
            let vpv = 8800, vph = 8000, vsv = 4000, vsh = 4400,
                    η = 1.01,
                    c = CIJ.vti(vpv=vpv, vph=vph, vsh=vsh, vsv=vsv, eta=η)
                @testset "Horizontal x1" begin
                    v = CIJ.phase_vels(c, 0, 0)
                    @test v.vp ≈ vph
                    @test v.vs1 ≈ vsh
                    @test v.vs2 ≈ vsv
                    @test v.pol ≈ 90 || v.pol ≈ -90
                    @test v.avs ≈ 200*(v.vs1 - v.vs2)/(v.vs1 + v.vs2)
                    @test v.xp ≈ [1, 0, 0]
                    @test v.xs1 ≈ [0, 1, 0] || v.xs1 ≈ [0, -1, 0]
                    @test v.xs2 ≈ [0, 0, 1] || v.xs2 ≈ [0, 0, -1]
                end
                @testset "Horizontal x2" begin
                    v = CIJ.phase_vels(c, -90, 0)
                    @test v.vp ≈ vph
                    @test v.vs1 ≈ vsh
                    @test v.vs2 ≈ vsv
                    @test v.pol ≈ 90 || v.pol ≈ -90
                    @test v.avs ≈ 200*(v.vs1 - v.vs2)/(v.vs1 + v.vs2)
                    @test v.xp ≈ [0, 1, 0]
                    @test v.xs1 ≈ [1, 0, 0] || v.xs1 ≈ [-1, 0, 0]
                    @test v.xs2 ≈ [0, 0, 1] || v.xs2 ≈ [0, 0, -1]
                end
                @testset "Horizontal 45" begin
                    v = CIJ.phase_vels(c, -45, 0)
                    @test v.vp ≈ vph
                    @test v.vs1 ≈ vsh
                    @test v.vs2 ≈ vsv
                    @test v.pol ≈ 90 || v.pol ≈ -90
                    @test v.avs ≈ 200*(v.vs1 - v.vs2)/(v.vs1 + v.vs2)
                    @test v.xp ≈ [√2/2, √2/2, 0]
                    @test v.xs1 ≈ [-√2/2, √2/2, 0] || v.xs1 ≈ [√2/2, -√2/2, 0]
                    @test v.xs2 ≈ [0, 0, 1] || v.xs2 ≈ [0, 0, -1]
                end
                @testset "Vertical" begin
                    v = CIJ.phase_vels(c, 0, 90)
                    @test v.vp ≈ vpv
                    @test v.vs1 ≈ vsv
                    @test v.vs2 ≈ vsv
                    @test v.avs ≈ 0 atol=√eps()
                    @test v.xp ≈ [0, 0, 1]
                end
            end
        end

        @testset "Olivine" begin
            let (c, ρ) = CIJ.ol()
                # [100], corresponds to c1111 for P, or C11, and c1212 and c1313 for S,
                # or C66 and C55
                @testset "[100]" begin
                    v = CIJ.phase_vels(c, 0, 0)
                    @test v.vp ≈ √c[1,1]
                    @test v.vs1 ≈ √c[6,6]
                    @test v.vs2 ≈ √c[5,5]
                    @test v.pol ≈ 90 || v.pol ≈ -90
                    @test v.avs ≈ 200*(v.vs1 - v.vs2)/(v.vs1 + v.vs2)
                    @test v.xp ≈ [1, 0, 0]
                    @test v.xs1 ≈ [0, 1, 0] || v.xs1 ≈ [0, -1, 0]
                    @test v.xs2 ≈ [0, 0, 1] || v.xs2 ≈ [0, 0, -1]
                end

                # c2222 for P, c1212 and c2323 for S, or C66 and C44
                @testset "[010]" begin
                    v = CIJ.phase_vels(c, -90, 0)
                    @test v.vp ≈ √c[2,2]
                    @test v.vs1 ≈ √c[6,6]
                    @test v.vs2 ≈ √c[4,4]
                    @test v.pol ≈ 90 || v.pol ≈ -90
                    @test v.avs ≈ 200*(v.vs1 - v.vs2)/(v.vs1 + v.vs2)
                    @test v.xp ≈ [0, 1, 0]
                    @test v.xs1 ≈ [1, 0, 0] || v.xs1 ≈ [-1, 0, 0]
                    @test v.xs2 ≈ [0, 0, 1] || v.xs2 ≈ [0, 0, -1]
                end

                # Compared to CIJ_phasevels from
                # https://github.com/anowacki/elasticity
                @testset "(30, 40)" begin
                    v = CIJ.phase_vels(c, 30, 40)
                    @test v.vp ≈ 8557.0 atol=0.1
                    @test v.vs1 ≈ 5367.3 atol=0.1
                    @test v.vs2 ≈ 4625.1 atol=0.1
                    @test v.pol ≈ -32.7 atol=0.1
                    @test v.avs ≈ 14.8554 atol=0.0001
                    @test v.xp ≈ [0.768, -0.308, 0.560] atol=0.01
                    @test v.xs1 ≈ [0.634, 0.250, -0.731] atol=0.01
                    @test v.xs2 ≈ [0.085, 0.917, 0.388] atol=0.01
                end
            end
        end
    end
end