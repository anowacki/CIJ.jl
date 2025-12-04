# Type definitions and interface definitions

"""
    EC{T} <: AbstractArray{T,2}

Wrapper type for elastic constants held as a 6×6 Voigt matrix.  The underlying storage
is a `StaticArrays.MArray`, meaning certain operations using this structure will be
quicker than using normal `Array`s.  For instance, elastic constant transformations
such as rotations, or conversion between stiffness and compliance, are each about
twice as fast.

!!! note
    `EC`s are based on `StaticArrays.MArray`s, and are therefore mutable.
    However, if mutating an `EC` in-place, avoid broadcasting

---

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
    # Uninitialised
    EC{T}() where T = new{T}(MArray{Tuple{6,6},T,2,36}(undef))
    EC() = EC{DEFAULT_FLOAT}()
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
        x.data[i,j] = x[j,i]
    end
    x
end

EC(x::AbstractArray{T}; kwargs...) where {T<:Number} = EC{float(T)}(x; kwargs...)
EC(x::NTuple{36,T}; kwargs...) where T = EC{float(T)}(x; kwargs...)
EC(x; kwargs...) = EC{DEFAULT_FLOAT}(x; kwargs...)

Base.size(x::EC) = (6,6)
Base.length(x::EC) = 36
Base.iterate(x::EC) = Base.iterate(x.data)
Base.iterate(x::EC, i) = Base.iterate(x.data, i)
Base.getindex(x::EC, i...) = x.data[i...]
function Base.setindex!(x::EC, val, i, j)
    if i != j
        # @info "x[$j,$i] = $val (symmetry)"
        x.data[j,i] = val
    end
    # @info "x[$i,$j] = $val"
    x.data[i,j] = val
end
Base.LinearIndices(x::EC) = LinearIndices(x.data)
Base.IndexStyle(::Type{EC}) = IndexStyle(MMatrix)
Base.BroadcastStyle(::Type{<:EC}) = Base.BroadcastStyle(MMatrix)
# Base.BroadcastStyle(::Type{T}) where {T<:EC} = Base.BroadcastArrayStyle{T}()
# Base.similar(bc::Base.Broadcasted{})
Base.copyto!(x::EC, bc::Base.Broadcast.Broadcasted{Nothing}) = copyto!(x.data, bc)
Base.broadcasted(f::F, x::EC, arg2, args...) where F = Base.broadcasted(f, x.data, arg2, args...)
Base.broadcasted(::typeof(Base.identity), x::EC, arg) = Base.broadcasted(Base.identity, x.data, arg)
# Base.materialize!(x::EC, bc::Base.Broadcast.Broadcasted) = Base.materialize!(x.data, bc)

Base.zero(::Type{EC{T}}) where T = EC{T}(zero(MMatrix{6,6,T}))
Base.zero(::Type{EC}) = zero(EC{DEFAULT_FLOAT})
Base.zero(::EC{T}) where T = zero(EC{T})

Base.similar(::EC, ::Type{T}) where T = EC{T}()
Base.similar(::EC{T}) where T = EC{T}()
Base.similar(::Type{EC{T}}) where T = EC{T}()
Base.similar(::Type{<:EC}, ::Type{T}) where T = EC{T}()

Base.copy(x::EC) = EC(copy(x.data))

Base.isapprox(a::EC, b::EC) = Base.isapprox(a.data, b.data)
Base.isapprox(a::EC, b::AbstractArray) = Base.isapprox(a.data, b)
Base.isapprox(a::AbstractArray, b::EC) = Base.isapprox(a, b.data)
Base.isapprox(a, b::EC) = Base.isapprox(a, b.data)

# Conversion to StaticArrays types
Base.convert(::Type{SA}, x::EC) where {SA<:StaticArray} = SA(x.data)
