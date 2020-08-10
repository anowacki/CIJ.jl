# Type definitions and interface definitions

"""
    EC{T} <: AbstractArray{T,2}

Wrapper type for elastic constants held as a 6×6 Voigt matrix.  The underlying storage
is a `StaticArrays.MArray`, meaning certain operations using this structure will be
quicker than using normal `Array`s.  For instance, elastic constant transformations
such as rotations, or conversion between stiffness and compliance, are each about
twice as fast.

    EC(data::AbstractArray{T}; warn=false) where T -> ec

Create a set of elastic constants `ec` from the array `data` which must have length 36.
Symmetry is enforced when converting, taking the upper half of the matrix.  No warning
is issued for non-symmetric input, unless `warn` is `true`.
"""
struct EC{T} <: AbstractArray{T,2}
    data::MArray{Tuple{6,6},T,2,36}

    # Arbitrary input requires symmetrisation and conversion to MMatrix if not already
    function EC{T}(data; warn=false) where T
        if data isa AbstractArray
            size(data) == (6,6) || throw(ArgumentError("ECs must have dimensions 6×6"))
        elseif data isa NTuple
            length(data) == 36 || throw(ArgumentError("ECs must have length 36"))
        end
        ec = new{T}(MMatrix{6,6,T}(data))
        _make_symmetric!(ec; warn=warn)
    end
    # Construction with another EC doesn't require symmetrisation
    EC{T}(ec::EC; warn=false) where T = new{T}(ec.data)
end

"Enforce symmetry in a 6×6 array.  No checks on shape performed."
function _make_symmetric!(x; warn=false)
    if warn
        atol = √eps(eltype(x))
        @inbounds for i in 1:6, j in i+1:6
            if !isapprox(x[i,j], x[j,i], atol=atol)
                @warn("input matrix not symmetrical: taking upper half")
                break
            end
        end
    end
    # Copy upper to lower half
    @inbounds for j in 1:6, i in j+1:6
        x[i,j] = x[j,i]
    end
    x
end

EC(x::AbstractArray{T}; kwargs...) where {T<:Number} = EC{float(T)}(x; kwargs...)
EC(x::NTuple{36,T}; kwargs...) where T = EC{float(T)}(x; kwargs...)
EC(x; kwargs...) = EC{DEFAULT_FLOAT}(x; kwargs...)

Base.size(x::EC) = (6,6)
Base.length(x::EC) = 36
Base.getindex(x::EC, i...) = x.data[i...]
Base.setindex!(x::EC, val, i, j) = x.data[i,j] = x.data[j,i] = val
Base.IndexStyle(::Type{EC}) = IndexStyle(MMatrix)

Base.zero(::Type{EC{T}}) where T = EC{T}(zero(MMatrix{6,6,T}))
Base.zero(::Type{EC}) = zero(EC{DEFAULT_FLOAT})
Base.zero(::EC{T}) where T = zero(EC{T})
