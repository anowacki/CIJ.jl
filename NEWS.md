# CIJ.jl v0.2.2 release notes

## New features and non-breaking changes
- `CIJ.plot_hemisphere` plots a 2D view of an upper hemisphere showing
  phase velocities and fast orientations.
- `CIJ.plot_hemisphere!` does the same, but plots into an existing
  `Makie.PolarAxis`.


# CIJ.jl v0.2.1 release notes

## New features and non-breaking changes
- `CIJ.phase_vels` now also returns particle motions for P, S1 and
  S2.
- `CIJ.phase_vels` no longer allocates at all and is faster.
- If you load a backend for the [Makie.jl](https://docs.makie.org/stable/)
  plotting package by doing e.g. `using GLMakie`, `using CairoMakie`,
  `using WGLMakie`, etc., then you can plot phase velocities and
  anisotropy using the new `CIJ.plot_sphere` function.
  `CIJ.plot_sphere!` is also available to plot into an existing
  `Makie.Axis3`.


# CIJ.jl v0.2.0 release notes

## Breaking changes

- `CIJ.write` no longer extends `Base.write` and so must be called
  with the module name.
- `CIJ.write` also now places the file name or `IO` object first,
  for consistency with `Base.write`.

## New features and non-breaking changes
- Julia v1.6 is now the minimum required version.
- `CIJ.read` and `CIJ.write` can now read/write from/to an `IO`
  object, like an `IOBuffer`.

## Notable bug fixes
- Several incorrect calls to `warn` (now `@warn`) were updated.
