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
- `:avs`: Shear wave anisotropy

# Keyword arguments
- `ax_kwargs = ()`: `Dict`, named tuple or set of pairs giving extra keyword
  arguments to control the axis into which the sphere is plotted.  Passed
  to `Makie.Axis3`.
- `fig_kwargs = ()`: `Dict`, named tuple or set of pairs giving extra keyword
  arguments to control the figure into which the sphere is plotted.  Passed
  to `Makie.Figure`.
- `colorbar = true`: Whether (default) or not to show a colour scale on the
  right of the plot.
- `directions = ()`: Tuple of two things: (1) vector of directions and (2)
  text to plot at each of these directions.  Allows annotated
- `level = 3`: Level of sphere refinement.  The larger the number, the finer
  the mesh of points sampling the sphere and the smoother the surface appears.
  Numbers above 6 or 7 may become slower and are unnecessarily fine.
"""
function CIJ.plot_sphere(
    C::CIJ.EC,
    property=:vp;
    fig_kwargs=(),
    ax_kwargs=(),
    colorbar=true,
    kwargs...
)
    label = if property == :vp
        "Vₚ (m/s)"
    elseif property == :vs1
        "Vₛ₁ (m/s)"
    elseif property == :vs2
        "Vₛ₂ (m/s)"
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
        limits=1.0.*(-1, 1, -1, 1, -1, 1),
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
    CIJ.plot_sphere!(ax::Makie.AbstractAxis, C::CIJ.ECs, property=:vp; level=3)

Like [`plot_sphere`](@ref CIJ.plot_sphere!), but insert the plot into an
existing `Makie.AbstractAxis` (like `Makie.Axis3`).  See
[`plot_sphere`](@ref CIJ.plot_sphere) for details of arguments and keyword
arguments.
"""
function CIJ.plot_sphere!(
    ax::Makie.AbstractAxis,
    C::CIJ.EC,
    property=:vp;
    level=3,
)
    # Create even sampling of the sphere
    t = TS.Tesselation(level)

    # Convert to something Makie can plot
    points = reinterpret(GeometryBasics.Point3{eltype(eltype(t.points))}, t.points)
    mesh = GeometryBasics.Mesh(points, [GeometryBasics.TriangleFace(SVector(tt)...) for tt in t.triangles])

    # Calculate phase velocities at points
    azi = atand.(first.(points), -getindex.(points, 2))
    inc = asind.(last.(points))
    output = CIJ.phase_vels.(Ref(C), azi, inc)

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

    # Primary directions
    Makie.arrows!(
        fill(Makie.Point(0, 0, 0), 3),
        [Makie.Vec3(1.1.*v...) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))];
        arrowsize=Makie.Vec3(0.15, 0.15, 0.2),
        color=:white
    )
    Makie.text!([1.2, 0, 0], [0, 1.2, 0], [0.1, 0.1, 1.2], text=["x₁", "x₂", "x₃"])

    # Other directions
    # for (dir, label) in arrows
    #     Makie.arrows!()
    # end

    pl
end

function _cart2aziinc(x, y, z)
    azi = atand.(x, -y)
    inc = asind(z)
    azi, inc
end
_cart2aziinc(r::AbstractVector) = _cart2aziinc(r[1], r[2], r[3])


end # module
