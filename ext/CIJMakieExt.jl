module CIJMakieExt

import CIJ
import GeometryBasics
import Makie

using StaticArrays: SVector

include("TesselatedSphere.jl")
import .TesselatedSphere as TS

"""
    hemisphere_axis(subfig; kwargs...) -> ::Makie.PolarAxis

Return a `Makie.PolarAxis` which represents an upper hemisphere into
which a call to `CIJ.plot_hemisphere!` can be made.

`kwargs` are keyword arguments passed to the `Makie.PolarAxis`
constructor.
"""
function CIJ.hemisphere_axis(subfig; kwargs...)
    ax = Makie.PolarAxis(subfig;
        direction=-1,
        rgridvisible=false,
        rlimits=(0, 1),
        rticklabelsvisible=false,
        thetagridvisible=false,
        theta_0=-π/2,
        thetaticks=([0, 3π/2], ["x₁", "x₂"]),
        kwargs...
    )
end

"""
    CIJ.plot_hemisphere(C, properties=(:vp, :avs); kwargs...) -> ::Makie.Figure

Return a plot of phase velocities of a Voigt elasticity matrix `C`
on an upper hemisphere.  One hemisphere is shown per entry in `properties`.

`properties` is a set (e.g., vector or tuple) or `Symbol`s, of which each
can be one of the following:

- `:vp`: P-wave velocity
- `:vs1`: Fast shear-wave velocity
- `:vs2`: Slow shear-wave velocity
- `:avs`: Shear-wave anisotropy (200⨯(Vₛ₁ – Vₛ₂)/(Vₛ₁ + Vₛ₂) %)

Velocities depend on the units of `C`; if `C` is normalised by density and
hence in m²/s², then velocities are in m/s.  The default axis titles
assume this is the case and if velocities are in, say, km/s, then
you can pass a different argument to the `units` keyword argument.

# Keyword arguments
- `fast_dirs = (false, true)`: If this is a single `Bool`, then plot fast
  shear-wave orientations on all axes.  If it is a set (e.g., vector or
  tuple) of `Bool`s, then the `i`th `property` has fast S-wave orientations
  shown if `fast_dirs[i]` is `true`.
- `height = 300`: Default hemisphere height in pixels.  Set to `nothing`
  to allow hemispheres to grow when resizing the figure.
- `projection = :infinity`: Projection to use when plotting.  Can be one of:
  - `:infinity`: Azimuthal projection at infinity; i.e., as if viewing a sphere
    from infinite distance.
  - `:lambert`: Lambert equal-area azimuthal projection.
- `resize = true`: Resize the figure before returning to fit the plot.
- `spacing = 2.5`: Spacing in ° for the grid sampling the hemisphere.
- `ticks = (color=:black, markersize=20, strokecolor=:black, strokewidth=1)`:
  `NamedTuple` or `Dict` containing keyword arguments passed to `Makie.scatter`
  when plotting the fast orientations (if `fast_dirs` is `true`).
- `units = "m/s"`: Assumed units of velocities when plotted.
"""
function CIJ.plot_hemisphere(
    C,
    properties=(:vp, :avs);
    fast_dirs=(properties .== :avs),
    height=300,
    projection=:infinity,
    resize::Bool=true,
    spacing=2.5,
    ticks=(),
    units="m/s",
)
    isempty(properties) && throw(ArgumentError("`properties` cannot be empty"))

    all(in((:vp, :vs1, :vs2, :avs)), properties) ||
        throw(ArgumentError(
            "`properties` must be an iterable containing only: `:vp`, `:vs1`, `:vs2` or `:avs`"
        ))

    # Plot coordinates and phase velocities
    (; θs, rs, values) = _phase_vels_and_hemisphere_coords(C, spacing; projection)

    fig = Makie.Figure()

    for (iprop, (property, fast)) in enumerate(zip(properties, Iterators.cycle(fast_dirs)))
        label = _label_from_property(property, units)

        icol = 2iprop - 1

        ax = CIJ.hemisphere_axis(fig[1,icol]; height=height, width=height, title=label)

        plot_vals = getproperty.(values, property)

        pl = Makie.contourf!(ax, θs, rs, plot_vals; colormap=Makie.Reverse(:turbo))

        # Update title with P-wave anisotropy for P-wave plots
        if property == :vp
            vp_min, vp_max = extrema(plot_vals)
            avp = 200*(vp_max - vp_min)/(vp_max + vp_min)
            ax.title = label * " (AVₚ = $(round(avp; sigdigits=3)) %)"
        end

        Makie.tightlimits!(ax)

        cb = Makie.Colorbar(fig[1,icol+1], pl;
            height=(isnothing(height) ? nothing : 0.7*height),
            tellheight=true,
        )
        
        Makie.colgap!(fig.layout, icol, 2)

        # Plot fast orientations
        if fast
            θs_fast, rs_fast, pols = _pols_and_hemisphere_coords(C, 2; projection)
            Makie.scatter!(ax, θs_fast, rs_fast;
                rotation=(-deg2rad.(pols) .- θs_fast),
                color=:black,
                marker=:vline,
                markersize=20,
                strokecolor=:white,
                strokewidth=1,
                ticks...
            )
        end
    end

    resize && Makie.resize_to_layout!(fig)

    fig
end

"""
    CIJ.plot_hemisphere(C, property::Symbol; kwargs...) -> ::Makie.Figure

Create a plot of a single `property`.
"""
CIJ.plot_hemisphere(C, properties::Symbol; kwargs...) =
    CIJ.plot_hemisphere(C, (properties,); kwargs...)


"""
    CIJ.plot_hemisphere!(ax::Makie.PolarAxis, C, property::Symbol=:vp; fast_dirs=(property == :avs), projection=:infinity, spacing=2.5, levels=10, ticks) -> ::Makie.Plot

Plot a single upper hemisphere plot of phase velocities into an existing
`Makie.PolarAxis` `ax` for a Voigt elastic constants matrix `C`.

See [`plot_hemisphere`](@ref CIJ.plot_hemisphere) for more details,
including keyword arguments.  Note that only those listed above can be used
for this function.

The following additional keyword arguments are only applicable to this function:
- `levels = 10`: `Int` containing the number of contours to plot, or a
  `AbstractVector{<:Real}` giving the values of the contours to plot.
"""
function CIJ.plot_hemisphere!(
    ax::Makie.PolarAxis,
    C,
    property::Symbol=:vp;
    fast_dirs=(property == :avs),
    levels=10,
    projection=:infinity,
    spacing=2.5,
    ticks=(),
)
    (; θs, rs, values) = _phase_vels_and_hemisphere_coords(C, spacing; projection)
    plot_vals = getproperty.(values, property)

    pl = Makie.contourf!(ax, θs, rs, plot_vals; colormap=Makie.Reverse(:turbo), levels)

    # Plot fast orientations
    if fast_dirs
        θs_fast, rs_fast, pols = _pols_and_hemisphere_coords(C, 2; projection)
        Makie.scatter!(ax, θs_fast, rs_fast;
            rotation=(-deg2rad.(pols) .- θs_fast),
            color=:black,
            marker=:vline,
            markersize=20,
            strokecolor=:white,
            strokewidth=1,
            ticks...
        )
    end

    pl
end

"""
    _phase_vels_and_hemisphere_coords(C, spacing; projection=:lambert) -> (; θs, rs, values)

Calculate the output of `CIJ.phase_vels(C, azi, inc)` on an azimuth-inclination
grid spaced `spacing`° apart.  These are returned in a named tuple as
`:values`.  Also return the plotting coordinates `θs` and `rs`.
`projection` is passed to [`_hemisphere_coords`](@ref).
"""
function _phase_vels_and_hemisphere_coords(C, spacing; projection=:lambert)
    # Sampling points
    azis = 0:spacing:360
    incs = 0:spacing:90

    # Matrix of phase velocities
    values = CIJ.phase_vels.(Ref(C), azis, incs')

    θs, rs = _hemisphere_coords(azis, incs; projection)

    (; θs, rs, values)
end

"""
    _pols_and_hemisphere_coords(C, level; projection=:lambert) -> θs, rs, pols

Calculate the fast shear-wave polariations `pols` at a set of evenly distributed
points about the upper hemisphere, defined by a tesselation of the sphere
at `level`.  Returns these as well as the azimuths `θs` and radii `rs` at
which to plot these points into a `Makie.PolarAxis`.  `projection` is
passed on to [`_hemisphere_coords`](@ref).
"""
function _pols_and_hemisphere_coords(C, level; projection=:lambert)
    t = TS.Tesselation(level)
    incs, azis = _cart2incaz(filter!(p -> p.z >= 0, t.points))
    pols = getproperty.(CIJ.phase_vels.(Ref(C), azis, incs), :pol)

    θs, rs = _hemisphere_coords(azis, incs; projection)

    θs, rs, pols
end

"""
    _hemisphere_coords(azis, incs; projection=:lambert) -> θs, rs

Calculate the polar plotting coordinates (azimuths `θs` and radii `rs`)
for the azimuths `azis` and inclinations `incs` (both degrees) in
the CIJ convention.  (See [`CIJ.phase_vels`](@ref).)

`projection` defines the convention used.
"""
function _hemisphere_coords(azis, incs; projection=:lambert)
    if projection == :lambert
        θs = deg2rad.(azis)
        # Maximum radius is √2, so we need to divide by √2 to get rs in the
        # range [0, 1].  The formula has a 2 in front, so the prefactor becomes
        # 2/√2 = √2
        rs = √2 .* sind.((90 .- incs)./2)
    elseif projection == :infinity
        θs = deg2rad.(azis)
        rs = sind.(90 .- incs)
    else
        throw(ArgumentError("unsupported projection '$projection'"))
    end

    θs, rs
end    

"""
    CIJ.plot_sphere(C::CIJ.ECs, property=:vp; kwargs...) -> figure, axis, plot

Create a plot of the variation of `property` over all directions for the
elastic constants in `C` and show this as a 3D sphere.  Return a
`Makie.FigureAxisPlot` object which can be iterated to give the whole
`figure` object, plus the sphere `axis` and `plot` objects which can then
be manipulated as usual Makie objects.

`property` may be one of:
- `:vp` (default): P-wave phase velocity
- `:vs1`: Fast S-wave phase velocity
- `:vs2`: Slow S-wave phase velocity
- `:avs`: Shear wave anisotropy (200⨯(Vₛ₁ – Vₛ₂)/(Vₛ₁ + Vₛ₂) %)

Velocities depend on the units of `C`; if `C` is normalised by density and
hence in m²/s², then velocities are in m/s.  The default colorbar label
assumes this is the case and if velocities are in, say, km/s, then
you can pass a different argument to the `units` keyword argument.

# Keyword arguments
- `ax_kwargs = ()`: `Dict`, named tuple or set of pairs giving extra keyword
  arguments to control the axis into which the sphere is plotted.  Passed
  to `Makie.Axis3`.  (Cannot be passed to `plot_sphere!`.)
- `fig_kwargs = ()`: `Dict`, named tuple or set of pairs giving extra keyword
  arguments to control the figure into which the sphere is plotted.  Passed
  to `Makie.Figure`.  (Cannot be passed to `plot_sphere!`.)
- `colorbar = true`: Whether (default) or not to show a colour scale on the
  right of the plot.  (Cannot be passed to `plot_sphere!`.)
- `directions = ()`: Tuple of two things: (1) vector of directions
  as tuples (azi°, inc°) and (2) text to plot at each of these directions.
- `fast_dirs`: Whether or not to plot ticks showing the orientation of the
  fast shear wave across the sphere.  Defaults to `true` unless plotting
  Vₚ.
- `p_dirs = false`: Whether or not to plot the P-wave particle motion orientation.
- `p_normals = false`: If `p_dirs` is `true`, then if `p_normals` is also `true`,
  plot the direction of the wave vector.  This can be useful to show where
  the P-wave particle motion differs significantly from the wave vector.
- `slow_dirs = false`: Whether or not to plot ticks showing the orientation of
  the slow shear wave.
- `level = 4`: Level of sphere refinement.  The larger the number, the finer
  the mesh of points sampling the sphere and the smoother the surface appears.
  Numbers above 6 or 7 may become slower and are unnecessarily fine.
- `units = m/s`: Units of velocity for colour scale label.  (Cannot be passed
  to `plot_sphere!`.)
"""
function CIJ.plot_sphere(
    C::CIJ.EC,
    property=:vp;
    fig_kwargs=(),
    ax_kwargs=(),
    colorbar=true,
    units="m/s",
    kwargs...
)
    # Throws if `property` is not something we can plot
    label = _label_from_property(property, units)

    fig = Makie.Figure(; fig_kwargs...)
    ax = Makie.Axis3(
        fig[1,1];
        aspect=:data,
        viewmode=:fit,
        limits=1.5.*(-1, 1, -1, 1, -1, 1),
        elevation=π/6,
        azimuth=π/4,
        ax_kwargs...
    )
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)

    pl = CIJ.plot_sphere!(ax, C, property; kwargs...)

    if colorbar
        Makie.Colorbar(fig[1,2], pl; label, height=Makie.Relative(0.6), width=10)
    end

    Makie.resize_to_layout!(fig)

    Makie.FigureAxisPlot(fig, ax, pl)
end

"""
    CIJ.plot_sphere!(ax::Makie.AbstractAxis, C::CIJ.ECs, property=:vp; kwargs...)

Like [`plot_sphere`](@ref CIJ.plot_sphere!), but insert the plot into an
existing `Makie.AbstractAxis` (like `Makie.Axis3`).  See
[`plot_sphere`](@ref CIJ.plot_sphere) for details of arguments and keyword
arguments.
"""
function CIJ.plot_sphere!(
    ax::Makie.AbstractAxis,
    C::CIJ.EC,
    property=:vp;
    fast_dirs::Bool=(property in (:vs1, :vs2, :avs)),
    slow_dirs::Bool=false,
    p_dirs::Bool=false,
    p_normals::Bool=false,
    level=4,
    directions=([], [])
)
    # Create even sampling of the sphere
    t = TS.Tesselation(level)

    # Convert to something Makie can plot
    points, mesh = _mesh_points_from_tesselation(t)

    # Calculate phase velocities at points
    incs, azis = _cart2incaz(points)
    output = CIJ.phase_vels.(Ref(C), azis, incs)

    # Get property of interest
    vals, color_limits = if property == :avs
        avs = getproperty.(output, :avs)
        clims = (0, maximum(avs))
        avs, clims
    elseif property in (:vp, :vs1, :vs2)
        vals_ = getproperty.(output, property)
        clims = extrema(vals_)
        vals_, clims
    else
        throw(ArgumentError("property must be one of :vp, :vs1, :vs2, :avs"))
    end

    # Plot surface
    pl = Makie.mesh!(ax, mesh; color=vals, colormap=Makie.Reverse(:turbo), colorrange=color_limits)

    # Fast orientations
    if fast_dirs || slow_dirs || p_dirs
        t_coarse = TS.Tesselation(2)
        p_coarse, m_coarse = _mesh_points_from_tesselation(t_coarse)
        incs_coarse, azis_coarse = _cart2incaz(p_coarse)
        output_coarse = CIJ.phase_vels.(Ref(C), azis_coarse, incs_coarse)

        Makie.scatter!(ax, 1.01 .* p_coarse, color=:black, markersize=0.05, markerspace=:data)

        if fast_dirs
            # Direction vector of S1
            xs1s = getproperty.(output_coarse, :xs1)
            fast_lines = _vector_lines(xs1s, p_coarse)
            Makie.lines!(ax, fast_lines, color=:black)
        end

        if slow_dirs
            xs2s = getproperty.(output_coarse, :xs2)
            slow_lines = _vector_lines(xs2s, p_coarse)
            Makie.lines!(ax, slow_lines, color=:gray)
        end

        if p_dirs
            # Direction vectors of P
            xps = getproperty.(output_coarse, :xp)
            p_lines = _vector_lines(xps, p_coarse)
            Makie.lines!(ax, p_lines, color=:red)

            # Ray directions (which may not align with `xps` for anisotropy)
            if p_normals
                normals = _vector_lines(p_coarse, p_coarse)
                Makie.lines!(ax, normals, color=:blue, linestyle=:dash)
            end
        end
    end

    # Primary directions
    Makie.arrows3d!(
        ax,
        fill(Makie.Point(0, 0, 0), 3),
        [Makie.Vec3(1.1.*v...) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))];
        taillength=0,
        shaftradius=0,
        tiplength=0.1,
        tipradius=0.05,
        color=:white
    )
    arrow_text_coords = ([1.2, 0, 0], [0, 1.2, 0], [0.1, 0.1, 1.2])
    Makie.text!(ax, arrow_text_coords...; text=["█▌", "█▌", "█▌"], align=(:center, :bottom), color=:white)
    Makie.text!(ax, arrow_text_coords...; text=["x₁", "x₂", "x₃"], align=(:center, :bottom))

    # Annotations at certain directions
    if !isempty(first(directions))
        ndirections = length(directions[1])
        if length(directions[2]) != ndirections
            throw(ArgumentError(
                "directions must contain a tuple of two vectors of the same length"
            ))
        end

        direction_vectors = map(directions[1]) do (azi, inc)
            Makie.Vec3((1.1 .* CIJ.incaz2cart(inc, azi))...)
        end
        text_positions = map(directions[1]) do (azi, inc)
            Makie.Point((1.2 .* CIJ.incaz2cart(inc, azi))...)
        end

        Makie.arrows3d!(
            ax,
            fill(Makie.Point(0, 0, 0), ndirections),
            direction_vectors;
            taillength=0,
            shaftradius=0,
            tiplength=0.1,
            tipradius=0.05,
            color=:white
        )
        Makie.text!(ax, text_positions; text=fill("█", length(direction_vectors)), color=:white)
        Makie.text!(ax, text_positions; text=directions[2])
    end

    return pl
end

"""
    _cart2incaz(rs::AbstractVector) -> incs::Vector, azis::Vector

Local version of [`CIJ.cart2incaz`](@ref) which takes vectors of vectors
and returns vectors of `azis` and `incs`.
"""
function _cart2incaz(rs::AbstractVector)
    inc_azs = CIJ.cart2incaz.(first.(rs), getindex.(rs, 2), last.(rs))
    first.(inc_azs), last.(inc_azs)
end

"""
    _label_from_property(property::Symbol, units) -> label::String

Return an appropriate label for the `property`, throwing an error
if `property` is not something we can plot.
"""
function _label_from_property(property::Symbol, units)
    label = if property == :vp
        "Vₚ ($units)"
    elseif property == :vs1
        "Vₛ₁ ($units)"
    elseif property == :vs2
        "Vₛ₂ ($units)"
    elseif property == :avs
        "AVₛ (%)"
    else
        throw(ArgumentError("property must be one of :vp, :vs1, :vs2, :avs"))
    end
end

"""
    _mesh_points_from_tesselation(t) -> points::Vector{<:GeometryBasics.Point3}, mesh::GeometryBasics.Mesh

Convert a `TesselatedSphere.Tesselation` `t` into a set of `points`
and a connected `mesh` for plotting.
"""
function _mesh_points_from_tesselation(t)
    points = reinterpret(GeometryBasics.Point3{eltype(eltype(t.points))}, t.points)
    triangles = reinterpret(
        GeometryBasics.TriangleFace{eltype(eltype(t.triangles))},
        t.triangles
    )
    mesh = GeometryBasics.Mesh(points, triangles)
    points, mesh
end

"""
    _vector_lines(vectors, points; point_scale=1.01, vector_scale=0.1) -> ::Vector{<:GeometryBasics.Point3}

Create a single vector of `GeometryBasics.Point3`s which consists of pairs
of points separated by a point containing `NaN`s.  The pairs of points
have each item in `points` ± each item of `vectors`.  Thus each line
segment is centred around `points` but has two ends either side, and its
orientation is determined by `vectors`.
"""
function _vector_lines(vectors, points; point_scale=1.01, vector_scale=0.1)
    # Work around https://github.com/MakieOrg/Makie.jl/issues/5298
    # and avoid creating a `Vector{Makie.Point3}`, since `Point3` is a
    # `UnionAll` and not a concrete type, and at present Makie can't
    # convert that automatically for plotting.
    V = eltype(first(vectors))
    P = eltype(first(vectors))
    T = promote_type(V, P, typeof(point_scale), typeof(vector_scale))
    # TODO: Replace with a version using `GeometryBasics.LineString`s?
    Iterators.flatten(
        (
            point_scale*p + vector_scale*x,
            point_scale*p - vector_scale*x,
            # Avoid creating a UnionAll
            GeometryBasics.Point3{T}(NaN32, NaN32, NaN32)
        )
        for (x, p) in zip(vectors, points)
    ) |> collect
end

end # module
