# Test calculations of properties of tensors

using Test, CIJ

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
            @test !is_stable(c′)
            c[1,1] = -c[1,1]
            @test !is_stable(c)
            c[1,1] = c[2,2]
            c[1,3] = c[1,2] = c[2,3] = c[1,1]
            @test !is_stable(c)
        end
    end
end