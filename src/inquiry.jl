# Enquire about properties of tensors

"""
    is_iso(C) -> ::Bool

Return `true` if `C` is isotropic.
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
and hence not dynamically stable.
"""
is_stable(C) = _is_posdef(C)
is_stable(C::EC) = _is_posdef(C.data)

function _is_posdef(a::AbstractMatrix)
    is_6x6(a) || throw(ArgumentError("matrix must be 6×6"))
    is_symm(a) || @warn("matrix not symmetrical: using upper half only")
    LinearAlgebra.isposdef(LinearAlgebra.Hermitian(a))
end
# Required since there is no working isposdef for StaticArrays
_is_posdef(a::StaticMatrix) = try
        LinearAlgebra.cholesky(a)
        true
    catch err
        if err isa LinearAlgebra.PosDefException
            false
        else
            rethrow(err)
        end
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
