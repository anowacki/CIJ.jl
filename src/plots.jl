"""
    plot_sphere!

Plot a spherical view of the phase velocities of an
elastic tensor into an existing plot axis.

To use this function, you must first load a backend for the
[Makie.jl](https://docs.makie.org/stable/) plotting package by doing e.g.
`using GLMakie`, `using CairoMakie`, `using WGLMakie`, etc.  Once you
have done this, the docstring will update.
"""
function plot_sphere! end

"""
    plot_sphere

Plot a spherical view of the phase velocities of an elastic tensor.

To use this function, you must first load a backend for the
[Makie.jl](https://docs.makie.org/stable/) plotting package by doing e.g.
`using GLMakie`, `using CairoMakie`, `using WGLMakie`, etc.  Once you
have done this, the docstring will update.
"""
function plot_sphere end

"""
    plot_hemisphere!

Plot an upper hemispherical view of the phase velocities of an
elastic tensor into an existing plot axis.

To use this function, you must first load a backend for the
[Makie.jl](https://docs.makie.org/stable/) plotting package by doing e.g.
`using GLMakie`, `using CairoMakie`, `using WGLMakie`, etc.  Once you
have done this, the docstring will update.
"""
function plot_hemisphere! end

"""
    plot_hemisphere

Plot an upper hemispherical view of the phase velocities of an elastic tensor.

To use this function, you must first load a backend for the
[Makie.jl](https://docs.makie.org/stable/) plotting package by doing e.g.
`using GLMakie`, `using CairoMakie`, `using WGLMakie`, etc.  Once you
have done this, the docstring will update.
"""
function plot_hemisphere end

"""
    hemisphere_axis

Return a `Makie.PolarAxis` into which upper hemisphere data can be plotted
by passing azimuth and inclination.

To use this function, you must first load a backend for the
[Makie.jl](https://docs.makie.org/stable/) plotting package by doing e.g.
`using GLMakie`, `using CairoMakie`, `using WGLMakie`, etc.  Once you
have done this, the docstring will update.
"""
function hemisphere_axis end
