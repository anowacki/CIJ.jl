# Convert between representations of tensors

"""
    cijkl(C) -> c

Convert the 6x6 Voigt matrix `C` into the 3x3x3x3 tensor `c`
"""
function cijkl(C)
    a = VOIGT_CONTRACTION_MATRIX
    if ! is_6x6(C)
        error("CIJ.cijkl: Input must be a 6x6 array")
    end
    c = MArray{Tuple{3,3,3,3},eltype(C),4,81}(undef)
    @inbounds for i = 1:3, j = 1:3, k = 1:3, l = 1:3
        c[i,j,k,l] = C[a[i,j],a[k,l]]
    end
    return c
end

"""
    cij(c) -> C

Convert the 3x3x3x3 elasticity tensor `c` into the 6x6 Voigt matrix `C`
"""
function cij(c)
    a = VOIGT_CONTRACTION_MATRIX
    if size(c) != (3,3,3,3)
        error("CIJ.cij: Input must be a 3x3x3x3 tensor")
    end
    C = zero(EC{float(eltype(c))})
    @inbounds for i = 1:3, j = i:3, k = 1:3, l = k:3
        C[a[i,j],a[k,l]] = c[i,j,k,l]
    end
    return C
end

"""
    C2S(C) -> S

Return the inverse of the stiffness matrix `C`, the compliance matrix `S`.
"""
C2S(C) = LinearAlgebra.inv(C)
C2S(C::EC{T}) where T = EC{T}(C2S(C.data))

"""
    S2C(S) -> C

Return the inverse of the compliance matrix `S`, the stiffness matrix `C`.
"""
S2C(S) = C2S(S)

"""
    C2S!(C) -> S

Invert the stiffness matrix `C` in place, giving the compliance matrix `S`.
"""
C2S!(C) = C .= LinearAlgebra.inv!(LinearAlgebra.lu(C))
C2S!(C::EC) = C .= C2S(C.data)

"""
    S2C!(S) -> C

Invert the compliance matrix `S` in place, giving the siffness matrix `C`.
"""
S2C!(S) = C2S!(S)
