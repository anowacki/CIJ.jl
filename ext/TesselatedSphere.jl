"""
# TesselatedSphere

Create tesselations of the unit sphere based on repeated subdivision of
a regular icosahedron.

- To create a tesselation, use [`Tesselation`](@ref).
- To find out how many points or sides a tesselation has, use [`npts_in_level`](@ref)
  or [`faces_in_level`](@ref).
- To compute the inverse, use [`level_from_npts`](@ref) and [`level_from_nfaces`](@ref).
- [`triangles`](@ref) returns the coordinates of the faces, useful for plotting.
"""
module TesselatedSphere

using LinearAlgebra: normalize
using StaticArrays: StaticArrays, FieldVector, Size, SVector, @SVector

export
    Tesselation,
    faces_in_level,
    npts_in_level,
    level_from_npts,
    level_from_nfaces,
    triangles

"Total amount of memory available on the system in bytes"
const TOTAL_MEMORY = Int(Sys.total_memory())


#=
    Types
=#

"""
    Point{T}

A Cartesian point in 3D space.  For fields that relate to the Earth:

- `x` is defined to be the vector pointing from the centre of the Earth out
  through the Greenwich meridian and the equator,
- `y` passes through the equator and 90° east, and
- `z` passes through the north pole.

`Point`s making up a `Tesselation` should be normalised to the unit sphere.
"""
struct Point{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

# Ensures that the result of adding a different StaticArray to a Point is still a Point
StaticArrays.similar_type(::Type{<:Point}, ::Type{T}, s::Size{(3,)}) where T = Point{T}

"""
    Triangle

The set of point indices describing the points which make up a triangle.
The three indices are `a`, `b` and `c`, and refer to the position in the
`points` field of a [`Tesselation`](@ref).
"""
struct Triangle <: FieldVector{3, Int}
    a::Int
    b::Int
    c::Int
end

"""
    Tesselation{T}

The set of points and triangles which tesselate the unit sphere.  The type parameter
`T` defines the element type used for the coordinates.

# Fields
- `level`: Subdivision level, starting at 0
- `npoints`: Number of points
- `ntriangles`: Number of faces
- `points`: `Vector` of `StaticVector`s; each one contains the coordinates of the points
- `barycentres`: `Vector{StaticVector}` containing the coodinates of the face barycentres
- `triangles`: `Vector{StaticVector}` giving the three indices of `points` for each face

Note that `points` use the private type [`TesselatedFields.Point`](@ref) which is a
`StaticArrays.FieldVector{3}`.  You can access the coordinates of a `Point` as if
it were a normal 3-vector like `x = tess.points[1][1]`, or using its fields
`x`, `y` and `z` like `x = tess.points[1].x`.  Triangles behave similarly,
but with fields `a`, `b` and `c`.

You can also index `points` with a triangle, giving three points like so:

```
julia> t = Tesselation(0);

julia> tri = t.points[t.triangles[1]] # The set of points forming the first triangle
3-element StaticArrays.SVector{3, TERRA.TesselatedFields.Point{Float64}} with indices SOneTo(3):
 [0.7236067977499789, -0.5257311121191336, 0.4472135954999579]
 [0.7236067977499789, 0.5257311121191336, 0.4472135954999579]
 [0.0, 0.0, 1.0]

julia> tri[1].x # The first x coordinate of the first triangle's first vertex
0.7236067977499789
```
"""
struct Tesselation{T}
    level::Int
    npoints::Int
    ntriangles::Int
    points::Vector{Point{T}} # Coordinates of points
    barycentres::Vector{Point{T}} # Coordinates of barycentres
    triangles::Vector{Triangle} # Indices of points forming triangles around barycentres
end

#=
    Constructors
=#

"""
    Tesselation([T=Float64,] level; polar=true, check_memory=true) -> t

Create a `Tesselation{T}` `t` subdivided from the icosahedron `level` times.
If `polar` is true, then there are points lying on the poles.  By default,
this constructor tests whether the tesselation is likely to be too big to fit in the
current machine's memory and throws and `ArgumentError` if not.  This can be turned
off by passing `check_memory=false`.

The element type of the tesselation coordinates is given by `T` and defaults to
`Float64`.
"""
function Tesselation(T, level::Integer; polar=true, check_memory=true)
    level >= 0 || throw(ArgumentError("level cannot be negative"))
    ntriangles = faces_in_level(level)
    npts = npts_in_level(level)

    # Don't try and compute a tesselation which won't fit in memory unless we skip this check
    if check_memory
        # If T == Float64, the calculation of how much memory is needed itself overflows
        # when level >= 27, so we cannot possibly compute such refined tesselations
        level > 25 && throw(ArgumentError("levels > 25 are too large to be computed"))
        # We can probably calculate the number of bytes required without overflowing
        if _tesselation_memory(T, level) > TOTAL_MEMORY
            throw(ArgumentError(
                "memory required to hold tesselation exceeds machine's total memory"))
        end
    end

    # The base level-0 case
    tesselation = _icosahedron(T, polar)

    # Subdivide to the correct level
    for _ in 1:level
        tesselation = _subdivide(tesselation)
    end

    # Barycentres were not calculated in _icosahedron or _subdivide
    @assert isempty(tesselation.barycentres)
    append!(tesselation.barycentres,
        compute_barycentres(tesselation.points, tesselation.triangles))

    return tesselation
end

Tesselation(level::Integer; kwargs...) = Tesselation(Float64, level; kwargs...)


#=
    Helper functions
=#

"""
    _subdivide(tesselation) -> tesselation′

Subdivide a tesselation so that the number of triangles increases by a factor
of four, and the point spacing decreases by about 2.  This increases `tesselation.level`
by 1.

# Subdivision scheme
In the code, each triangle at `tesselation.level` has corners with indices `a`, `b`
and `c`.  New points with indices `d`, `e` and `f` are inserted and new triangles
formed.
```
      + a
     / `
    /   `
 d +-----+ f
  / `   / `
 /   ` /   `
+-----+-----+
b     e     c
```

This follows the explanation at http://richardssoftware.net/Home/Post/60
"""
function _subdivide(tess::Tesselation{T}) where T
    level = tess.level + 1
    npoints = npts_in_level(level)
    ntriangles = faces_in_level(level)
    points = Vector{Point{T}}(undef, npoints)
    points[1:tess.npoints] .= tess.points
    triangles = Vector{Triangle}(undef, ntriangles)
    newpoints = Dict{Tuple{Int,Int}, Int}()

    it_new = 1
    ip_new = tess.npoints + 1

    for triangle in tess.triangles
        a, b, c = triangle
        d, ip_new = _add_new_points!(newpoints, points, a, b, ip_new)
        e, ip_new = _add_new_points!(newpoints, points, b, c, ip_new)
        f, ip_new = _add_new_points!(newpoints, points, c, a, ip_new)

        triangles[it_new] = Triangle(a, d, f)
        triangles[it_new + 1] = Triangle(b, e, d)
        triangles[it_new + 2] = Triangle(c, f, e)
        triangles[it_new + 3] = Triangle(d, e, f)
        it_new += 4
    end

    # Don't compute barycentres here; wait until we have completed all levels.
    Tesselation{T}(level, npoints, ntriangles, points, Point{T}[], triangles)
end

"""
    _add_new_points!(newpoints, points, a, b, ip_new) -> index, ip_new′

Given a pair of point indices `a` and `b`, determine if this pair have been previously
used in this subdivision process to create a new halfway point between them, and
return the `index` of the existing or new halfway point.
Otherwise, calculate this new point and add it to the list of `points`.  `newpoints`
is a dictionary where the key `(min(a,b), max(a,b))` gives the index of the new
halfway point.  `ip_new` is the index in `points` at which any new points should
be added, and if necessary an incremented `ip_new′` is returned.
"""
function _add_new_points!(newpoints, points, a, b, ip_new)
    key = a < b ? (a, b) : (b, a)
    index_try = get(newpoints, key, 0)
    if index_try !== 0
        return index_try, ip_new
    else
        index = ip_new
        points[index] = normalize(sum(@SVector[points[a], points[b]])/2)
        newpoints[key] = index
        return index, ip_new + 1
    end
end


"""
    _tesselation_memory(T, level) -> size_in_bytes

Compute the total memory usage in bytes required for a `Tesselation{T}`.
Note that this will be a very slight underestimate, as the calculation ignores
the extra memory used by Julia arrays beyond the raw data.
"""
function _tesselation_memory(T, level)
    3*(2*sizeof(T)*npts_in_level(level) + sizeof(Int)*faces_in_level(level))
end

"""
    _icosahedron(T, polar) -> tesselation

Create a `Tesselation{T}` at level 0, where the points are simply the vertices
of an icosahedron.  If `polar` is `true`, create an icosahedron with points at
the poles; otherwise there is no point at the pole.

The `polar == true` case is the arrangement more common in global climate simulations
(e.g., Sadourny, Arakawa and Mintz, Monthly Weather Review, 1968).
The `polar == false` case may be better when looking at data at the poles
(Teanby, Computers & Geophysics, 2006).

!!! note
    Barycentres are not computed using this internal function.
"""
function _icosahedron(T, polar)
    level = 0
    np = 12
    nt = 20

    points, triangles = if polar
        _polar_icosahedron(T)
    else
        _nonpolar_icosahedron(T)
    end

    tess = Tesselation(level, np, nt, points, Point{T}[], triangles)
end

"Compute the barycentre of a set of points on the units sphere."
compute_barycentre(points) = normalize(sum(points))

"""
    compute_barycentres(points, triangles) -> barycentres

Compute the barycentres of `triangles`, whose indices give the vertex locations
in `points`.  The result is a vector of barycentres corresponding to each triangle.
"""
function compute_barycentres(points, triangles)
    [compute_barycentre(@SVector[points[t[i]] for i in 1:3]) for t in triangles]
end

"""
    _polar_icosahedron(T) -> points, triangles

Return the set of `points` and `triangles` defining an icosahedron with points placed
at the poles, with point type `Point{T}`.

Coordinates are taken from http://mathworld.wolfram.com/Icosahedron.html
"""
function _polar_icosahedron(T)
    # Two points at poles, then around small circles at z = ±a = 1/√5
    a = oneunit(T)/√5
    b = 2a
    points = Point{T}.(
        [
            (           0,             0,   1),
            (b*cospi(1//5), b*sinpi(1//5),  a),
            (b*cospi(3//5), b*sinpi(3//5),  a),
            (b*cospi(5//5), b*sinpi(5//5),  a),
            (b*cospi(7//5), b*sinpi(7//5),  a),
            (b*cospi(9//5), b*sinpi(9//5),  a),
            (            b,             0, -a),
            (b*cospi(2//5), b*sinpi(2//5), -a),
            (b*cospi(4//5), b*sinpi(4//5), -a),
            (b*cospi(6//5), b*sinpi(6//5), -a),
            (b*cospi(8//5), b*sinpi(8//5), -a),
            (            0,             0, -1)
        ])
    triangles = Triangle.(
        [
            ( 6,  2,  1),
            ( 2,  3,  1),
            ( 3,  4,  1),
            ( 4,  5,  1),
            ( 5,  6,  1),
            (11,  7,  6),
            ( 7,  2,  6),
            ( 7,  8,  2),
            ( 8,  3,  2),
            ( 8,  9,  3),
            ( 9,  4,  3),
            ( 9, 10,  4),
            (10,  5,  4),
            (10, 11,  5),
            (11,  6,  5),
            (12,  7, 11),
            (12,  8,  7),
            (12,  9,  8),
            (12, 10,  9),
            (12, 11, 10)
        ])
    points, triangles
end

"""
    _nonpolar_icosahedron(T) -> points, triangles

Return the set of `points` and `triangles` defining an icosahedron without points placed
at the poles, with point type `Point{T}`.
"""
function _nonpolar_icosahedron(T)
    # Golden ratio
    ϕ = Base.MathConstants.φ
    points = normalize.(Point{T}.([
            ( 0,  ϕ,  1),
            ( 0, -ϕ,  1),
            ( 0,  ϕ, -1),
            ( 0, -ϕ, -1),
            ( 1,  0,  ϕ),
            (-1,  0,  ϕ),
            ( 1,  0, -ϕ),
            (-1,  0, -ϕ),
            ( ϕ,  1,  0),
            (-ϕ,  1,  0),
            ( ϕ, -1,  0),
            (-ϕ, -1,  0)
        ]))
    triangles = Triangle.([
            ( 2,  4, 11),
            ( 5,  2, 11),
            ( 9,  5, 11),
            ( 7,  9, 11),
            (11,  7,  4),
            ( 4,  2, 12),
            ( 6, 12,  2),
            ( 2,  5,  6),
            ( 1,  6,  5),
            ( 5,  9,  1),
            ( 3,  1,  9),
            ( 9,  7,  3),
            ( 8,  3,  7),
            ( 7,  4,  8),
            (12,  8,  4),
            (12,  6, 10),
            ( 6,  1, 10),
            ( 1,  3, 10),
            ( 3,  8, 10),
            (10, 12,  8)
        ])
    points, triangles
end

"""
    triangles(t::Tesselation) -> x, y, z

Return `x`, `y` and `z`, where each is a `Vector{SVector{4}}` giving the location
of the three vertices of the triangles in `t`, plus the first again to close the
loop.
"""
function triangles(t::Tesselation{TC}) where {TC}
    tx = Vector{SVector{4,TC}}(undef, t.ntriangles)
    ty = similar(tx)
    tz = similar(tx)
    for (i, tri) in pairs(t.triangles)
        points = t.points[tri]
        # Add the first point back to the beginning
        points = points[@SVector[1,2,3,1]]
        tx[i] = first.(points)
        ty[i] = getindex.(points, 2)
        tz[i] = last.(points)
    end
    tx, ty, tz
end

"""
    faces_in_level(level) -> nfaces

Return the number of faces (triangles) `nfaces` in a tesselation of `level`.
"""
faces_in_level(level::Integer) = 20*2^(2*level)

"""
    npts_in_level(level) -> npts::Int

Return the number of points at a given tesselation level `level`.
"""
npts_in_level(level::Integer) = Int(10*2^(2*level) + 2)

"""
    level_from_nfaces(nfaces) -> level

Return the `level` of a tesselation which contains `nfaces` faces (triangles).
The routine only searches up to level 10.  An error is thrown for an `nfaces`
which is not a valid number.
"""
function level_from_nfaces(nfaces::Integer)
    level::Int = 0
    while level <= 10
        faces_in_level(level) == nfaces && return level
        level += 1
    end
    error("$nfaces does not match any tesselation level in range 0 < level <= 10")
end

"""
    level_from_npts(npts) -> level::Int

Return the level of a tesselation from the number of points `npts`.

If the number of points is not valid for any tesselation level up to 10,
throw an error, which may be caught.
"""
function level_from_npts(npts::Integer)
    level::Int = 0
    while level <= 10
        npts_in_level(level) == npts && return level
        level += 1
    end
    error("$npts does not match any tesselation in range 0 < level < 10")
end

end # module
