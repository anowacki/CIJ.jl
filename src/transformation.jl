# Transformation of tensors

"""
    rot3(C, a, b, x) -> Cr

Return a rotated tensor `Cr`, the result of rotating in turn the tensor `C`:

1. About the x1 axis by `a` degrees clockwise when looking down the x1 axis
2.   "    "  x2   "   " `b`    "        "       "     "      "   "  x2   "
3.   "    "  x3   "   " `b`    "        "       "     "      "   "  x3   "

Rotations are performed in that order,
"""
rot3(C, a, b, c) = transform(C, RotXYZ(deg2rad(a), deg2rad(b), deg2rad(c)))

"""
    transform(C, M) -> C′

Apply a transformation matrix `M` to a 6x6 matrix `C`, returning the transformed matrix `C′`.
"""
function transform(C, M)
    K = _transformation_matrix(M)
    return K * (C * transpose(K))
end

"""
    _transformation_matrix(M) -> K

Return the matrix `K` which represents the matrix with which to transform a Voigt
6x6 matrix by the 3x3 transformation matrix `M`.
"""
@inline function _transformation_matrix(M)
    K = MMatrix{6, 6, Float64}(undef)
    @inbounds for i = 1:3
        i1 = (i + 1)%3
        i1 = i1 == 0 ? 3 : i1
        i2 = (i + 2)%3
        i2 = i2 == 0 ? 3 : i2
        for j = 1:3
            j1 = (j + 1)%3
            j1 = j1 == 0 ? 3 : j1
            j2 = (j + 2)%3
            j2 = j2 == 0 ? 3 : j2
            K[i,j] = M[i,j]^2
            K[i,j+3] = 2*M[i,j1]*M[i,j2]
            K[i+3,j] = M[i1,j]*M[i2,j]
            K[i+3,j+3] = M[i1,j1]*M[i2,j2] + M[i1,j2]*M[i2,j1]
        end
    end
    K
end
