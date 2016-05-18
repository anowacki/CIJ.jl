# CIJ.jl

## What is CIJ.jl?
A [Julia](http://julialang.org) package for dealing with linear elastic
constants, with particular applicability to geophysics problems.


## How to install
Although not registered as an official package, CIJ.jl can be added to your
Julia install like so:

```julia
Pkg.clone("https://github.com/anowacki/CIJ.jl")
```

You then need only do

```julia
import CIJ
```

and if that works, you're ready to go.


## How to use
CIJ.jl, to make your life as easy as possible, does not define a type for
elastic constants, but relies on them being arrays.  The basic kind is the
Voigt matrix, which is a 6 &times; 6 array.

Throughout CIJ.jl, we make no assumptions about the *units* of a matrix `C`,
but when it does matter (for instance, calculating phase velocities with
`phase_vels()`), then it is assumed you are dealing with *density-normalised*
constants (i.e., the units are m<sup>2</sup>&nbsp;s<sup>-2</sup>), sometimes called
*A<sub>ij</sub>* instead.  (This is the same as Pa/(kg&nbsp;m<sup>-3</sup>)),

Let's try out a simple example of what we can do:

```julia
C = zeroes(6, 6)
CIJ.is_stable(C) # -> false
```

Here, the `is_stable()` function tells one whether a set of elastic constants
`C` is dynamically possible.  It turns out, a material where all constants are
zero is not.  Olivine, however, should be, so

```julia
C, rho = CIJ.ol() # Return density-normalised constants, and density, for olivine
is_stable(C) # -> true
```

is not a surprise.

### Calculating phase velocities
One of the common uses for the package is to compute the *phase velocities* in
a given direction through some elastic constants.  This is done using the
`phase_vels()` function like so:

```julia
C, rho = CIJ.ol()
az, inc = 20, 45 # Directions
vp, vs1, vs2, pol, avs = CIJ.phase_vels(C, az, inc)
```

`vp`, `vs1` and `vs2` and so on, are velocities in m&nbsp;s<sup>-1</sup>,
`pol` is the orientation of the fast shear wave in degrees, and `avs` is the
percentage shear wave anisotropy along this direction.

### Converting to Voigt notation
If you have an 81-component tensor `c`, how do you get the Voigt matrix `C`?

```julia
C = CIJ.cij(c)
```

And the other way?

```julia
c = CIJ.cijkl(C)
```

## Getting help
Functions are documented, so at the REPL type `?` to get a `help?>` prompt,
and type the name of the function:

```julia
help?> CIJ.phase_vels

  phase_vels(C, az, inc) -> vp, vs1, vs2, pol, avs
  
  Calculate the phase velocities for the 6x6 elasticity matrix C
  along the direction (az, inc), in degrees, and return P-wave
  velocity vp, the fast and slow shear wave velocities, vs1 and
  vs2, the polarisation of the fast shear wave pol, and the shear
  wave velocity anisotropy, avs. Velocities are in m/s if the
  tensor C is in m^2/s^2 (i.e., is a density-normalised tensor,
  sometimes called A).
  
  az is the azimuth in degrees measured from the x1 towards to
  -x2 axis.
  
  inc is the inclination in degrees from the x1-x2 plane towards
  the x3 axis.
```

## Dependencies
None, apart from Julia itself.  I do not test CIJ.jl on versions older than
`v0.3`, however.


## Other software

* If you use MATLAB, then you should use [MSAT](https://github.com/andreww/MSAT/).
* If you use Fortran, then you should investigate the module
  `anisotropy_ajn` which is in the
  [seismo-fortran](https://github.com/anowacki/seismo-fortran) repo.


## Why the name?
Linear elastic constants are a fourth-rank tensor, relating the stress
_&sigma;_ in a material to the strain _&epsilon;_ with the relationship

*&sigma;<sub>ij</sub> = c<sub>ijkl</sub> &epsilon;<sub>kl</sub>*,

where _i_, _j_, _k_ and _l_ are indices taking values 1 to 3 and
representing the three cartesian directions in space.  This leads to there
being 3 &times; 3&times; 3 &times; 3 = 81 numbers in the tensor *c*, which
is somewhat unwieldy.

However, certain symmetries mean one can reduce this to 21, and represent
the 4-tensor with a symmetric 2-tensor or matrix instead; the so-called 'Voigt
notation'.  Typically, the lowercase 4-tensor *c<sub>ijkl</sub>* becomes the
uppercase matrix *C<sub>ij</sub>*, and thus the package is born.


## Acknowledgments
Credit goes to [James Wookey](http://www1.gly.bris.ac.uk/~wookey) and
[Mike Kendall](http://www1.gly.bris.ac.uk/~jmk/) for the original set of Fortran
routines on which the code is based, and which lives on in the
[seismo-fortran](https://github.com/anowacki/seismo-fortran) repo.

