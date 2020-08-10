# Enquire about properties of tensors

"""
    is_iso(C) -> ::Bool

Return `true` if `C` is isotropic.

# Examples
```
julia> is_iso(CIJ.iso(vp=8000, vs=4400))
true

julia> is_iso(CIJ.thom(8000, 4400, 0.01, 0.01, 0.01))
false
```
"""
function is_iso(C; atol=eps(eltype(C)))
    all(x->isapprox(C[1,1], x, atol=atol), (C[2,2], C[3,3])) &&
    all(x->isapprox(C[4,4], x, atol=atol), (C[5,5], C[6,6])) &&
    all(x->isapprox(C[1,1]-2C[4,4], x, atol=atol), (C[1,2], C[1,3], C[2,3])) &&
    all(x->isapprox(0, x, atol=atol), C[i,j] for i in 1:3 for j in 4:6) &&
    all(x->isapprox(0, x, atol=atol), (C[4,5], C[4,6], C[5,6]))
end

"""
    is_stable(C) -> ::Bool

Return `false` if the input 6x6 matrix `C` is not positive definite (symmetric)
and hence not dynamically stable, and `true` otherwise.

# Examples
```
julia> c = CIJ.ol()[1]
6×6 EC{Float64}:
 9.55291e7  2.02981e7  2.13413e7  0.0       0.0        0.0
 2.02981e7  5.85693e7  2.28912e7  0.0       0.0        0.0
 2.13413e7  2.28912e7  6.95976e7  0.0       0.0        0.0
 0.0        0.0        0.0        1.9076e7  0.0        0.0
 0.0        0.0        0.0        0.0       2.29508e7  0.0
 0.0        0.0        0.0        0.0       0.0        2.34575e7

julia> is_stable(c)
true

julia> is_stable(EC(rand(6, 6)))
false
```
"""
is_stable(C) = _isposdef(C)

# It's faster to convert to an Array and then use `isposdef`,
# than it is to catch a PosDefException.
_isposdef(c::StaticArray{Tuple{6,6}}) = LinearAlgebra.isposdef(Array(c))
_isposdef(c::EC) = _isposdef(c.data)

function _isposdef(a::AbstractMatrix)
    is_6x6(a) || throw(ArgumentError("matrix must be 6×6"))
    is_symm(a) || @warn("matrix not symmetrical: using upper half only")
    LinearAlgebra.isposdef(LinearAlgebra.Hermitian(a))
end

"""
    is_symm(C) -> ::Bool

Return `true` if `C` is symmetrical.
"""
function is_symm(C)
    for i = 1:6, j = i+1:6
        if !isapprox(C[i,j], C[j,i]) return false end
    end
    return true
end
is_symm(::EC) = true

"""
    is_6x6(C) -> ::Bool

Return `true` if C is a 6x6 matrix.
"""
is_6x6(C) = size(C) == (6,6)
is_6x6(::EC) = true
