# CIJ.jl v0.2.0 release notes

## Breaking changes

- `CIJ.write` no longer extends `Base.write` and so must be called
  with the module name.
- `CIJ.write` also now places the file name or `IO` object first,
  for consistency with `Base.write`.

## New features
- `CIJ.read` and `CIJ.write` can now read/write from/to an `IO`
  object, like an `IOBuffer`.

## Notable bug fixes
- Several incorrect calls to `warn` (now `@warn`) were updated.
