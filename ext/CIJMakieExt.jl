module CIJMakieExt

import CIJ
import GeometryBasics
import Makie

using StaticArrays: SVector

include("TesselatedSphere.jl")
import .TesselatedSphere as TS

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
- `:avs`: Shear wave anisotropy (200*(Vₛ₁ - Vₛ₂)/(Vₛ₁ + Vₛ₂) %)

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
- `directions = ()`: Tuple of two things: (1) vector of directions and (2)
  text to plot at each of these directions.
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

    fig = Makie.Figure(; fig_kwargs...)
    ax = Makie.Axis3(
        fig[1,1];
        aspect=:data,
        viewmode=:fit,
        limits=0.8.*(-1, 1, -1, 1, -1, 1),
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
)
    # Create even sampling of the sphere
    t = TS.Tesselation(level)

    # Convert to something Makie can plot
    points, mesh = _mesh_points_from_tesselation(t)

    # Calculate phase velocities at points
    incs, azis = _cart2incaz(points)
    output = CIJ.phase_vels.(Ref(C), azis, incs)

    # Get property of interest
    vals, color_limits = if property == :vp
        getproperty.(output, :vp), extrema(x -> x.vp, output)
    elseif property == :avs
        getproperty.(output, :avs), (0, maximum(x -> x.avs, output))
    elseif property == :vs1
        getproperty.(output, :vs1), extrema(x -> x.vs1, output)
    elseif property == :vs2
        getproperty.(output, :vs2), extrema(x -> x.vs2, output)
    else
        throw(ArgumentError("property must be one of :vp, :vs1, :vs2, :avs"))
    end

    # Plot surface
    pl = Makie.mesh!(ax, mesh, color=vals, colormap=Makie.Reverse(:turbo), colorrange=color_limits)

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
    Makie.arrows!(
        fill(Makie.Point(0, 0, 0), 3),
        [Makie.Vec3(1.1.*v...) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))];
        arrowsize=Makie.Vec3(0.15, 0.15, 0.2),
        color=:white
    )
    Makie.text!([1.2, 0, 0], [0, 1.2, 0], [0.1, 0.1, 1.2], text=["x₁", "x₂", "x₃"])

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
    # TODO: Replace with a version using `GeometryBasics.LineString`s?
    Iterators.flatten(
        (
            point_scale*p + vector_scale*x,
            point_scale*p - vector_scale*x,
            GeometryBasics.Point3(NaN32, NaN32, NaN32)
        )
        for (x, p) in zip(vectors, points)
    ) |> collect
end

end # module
