# Calculate velocities for tensors

"""
    phase_vels(C::EC, az, inc) -> (; vp, vs1, vs2, pol, avs, xp, xs1, xs2)
    phase_vels(C::AbstractMatrix, az, inc) -> (; vp, vs1, vs2, pol, avs, xp, xs1, xs2)

Calculate the phase velocities for the 6⨯6 elasticity matrix `C` along the direction
(`az`, `inc`), in degrees, and return the following in a named tuple:

- `:vp`: P-wave velocity `vp`
- `:vs1`: fast shear-wave velocity
- `:vs2`: slow shear-wave velocity
- `:pol`: the polarisation of the fast shear wave
- `:avs`: the shear wave velocity anisotropy (200(Vₛ₁ - Vₛ₂)/(Vₛ₁ + Vₛ₂) %)
- `:xp`: The particle motion vector of the P-wave.  This is guaranteed to
  be the direction of the two possible which points close to the wave vector.
- `:xs1`: The particle motion vector of the fast S-wave
- `:xs2`: The particle motion vector of the slow S-wave

Velocities are in m/s if the tensor `C` is in m^2/s^2 (i.e., is a
density-normalised tensor, sometimes called ``A``).  Vectors are given
as unit vectors with components along the x1, x2 and x3 directions.

`az` is the azimuth in degrees measured from the x1 towards the -x2 axis.

`inc` is the inclination in degrees from the x1-x2 plane towards the x3 axis.

`pol` is measured when looking towards the (0,0,0) point along the ray (in
the negative propagation direction).  `pol` increases clockwise away from the vector
normal to the propagation direction which points towards to the x3-axis.

# Example
```
julia> C, ρ = CIJ.ol();

julia> CIJ.phase_vels(C, 20, 45)
(vp = 8590.639816324323, vs1 = 5422.96811629596, vs2 = 4602.7034882853395, pol = -20.68250375350948, avs = 16.36328538101716, xp = [0.7600919761839179, -0.19404857608999868, 0.6201655729386067], xs1 = [0.6469844914122812, 0.13700188972411798, -0.7500943607867029], xs2 = [-0.06059088720936796, -0.9713782128138683, -0.22968045641220122])
```
"""
function phase_vels(C, az, inc)
    x = incaz2cart(inc, az)

    # Create the Christoffel matrix
    T = make_T(C, x)

    # Find eigenvectors and values, which correspond to velocities^2 (val[i])
    # and polarisation vectors (vec[:,i]).
    # They seem to be sorted into descending order of eigenvalue, but we double
    # check anyway.
    vals, vecs = LinearAlgebra.eigen(T)
    ip = argmax(vals)
    is2 = argmin(vals)
    is1 = 6 - ip - is2
    # Velocities
    vp  = sqrt(vals[ip])
    vs1 = sqrt(vals[is1])
    vs2 = sqrt(vals[is2])
    # Polarisations
    xp  = vecs[:,ip]
    xs1 = vecs[:,is1]
    xs2 = vecs[:,is2]

    # Flip direction of xp if it is pointing in about the opposite direction
    # to the wavevector
    if dot(x, xp) < 0
        xp = -xp
    end

    # Calculate S1 polarisation and amount of shear wave anisotropy
    pol = get_pol(inc, az, x, xs1)
    avs = 200*(vs1 - vs2)/(vs1 + vs2)

    return (vp=vp, vs1=vs1, vs2=vs2, pol=pol, avs=avs, xp=xp, xs1=xs1, xs2=xs2)
end

function group_vels(C, az, inc)
    # Return the group velocities in the az, inc direction for a 6x6
    # Voigt elasticity matrix C, in an elastic medium.
    # INPUT:
    #    C(6,6)  : Elasticity matrix
    #    az      : Azimuth (degrees) away from x1 towards -x2 axis
    #    inc     : Inclination (degrees) away from x1-x2 plane towards x3

    x = incaz2cart(inc, az)
    # Create the Christoffel matrix
    T = make_T(C, x)
    # Find eigenvectors, giving phase velocities
    vals, vecs = LinearAlgebra.eigen(T)
    ip = indmax(vals)
    is2 = indmin(vals)
    is1 = 6 - ip - is2
    vp = sqrt(vals[ip])
    vs1 = sqrt(vals[is1])
    vs2 = sqrt(vals[is2])
    xp = vec(vecs[:,ip])
    xs1 = vec(vecs[:,is1])
    xs2 = vec(vecs[:,is2])
    # Calculate group velocities, which are phase velocities projected onto
    # group velocity directions
    pp = xp./vp
    ps1 = xs1./vs1
    ps2 = xs2./vs2
    p = [pp'; ps1'; ps2']
    vg = zero(MVector{3,float(eltype(C))})
    ijkl = VOIGT_CONTRACTION_MATRIX
    for i=1:3, j=1:3, k=1:3, l=1:3
        m = ijkl[i,j]
        n = ijkl[k,l]
        vg[i] = vg[i] + C[m,n]*p[i,k]*x[j]*x[l]
    end
    return vg[1], vg[2], vg[3]
end

"Projection of fast shear wave polarisation onto wavefront plane"
function get_pol(inc, az, x, xs1)    
    xs1p = cross(x, cross(x, xs1))
    xs1p /= norm(xs1p)
    # Local up vector
    u = incaz2up(inc, az)
    # Angle between the local up vector and the projected fast orientation
    v = cross(xs1p, u)
    pol = rad2deg(atan(norm(v), dot(xs1p, u)))
    # If v is codirectional with x, then the angle is correct; otherwise, we're
    # the wrong way round
    if dot(x, v) < 0
        pol = -pol
    end
    # Put in range -90° to 90°
    return mod(pol + 90.0, 180.0) - 90.0
end

"""
    incaz2up(inc, az) -> ::SVector{3}

For a vector pointing along a given azimuth `az` and inclination `inc`,
return the local 'up' direction, which is normal to the original direction and
within the plane defined by it and the z-axis.
"""
function incaz2up(inc, az)
    a = deg2rad(az) + pi
    i = deg2rad(inc)
    sini = sin(i)
    return @SVector [cos(a)*sini, -sin(a)*sini, cos(i)]
end

"""
    make_T(C, x) -> T

Create the Christoffel matrix, `T`, given the Voigt matrix `C` and unit vector `x`.
"""
function make_T(C, x)
    T = zero(MMatrix{3,3,float(eltype(C))})
    ijkl = VOIGT_CONTRACTION_MATRIX
    @inbounds for i = 1:3, j = 1:3, k = 1:3, l = 1:3
        m = ijkl[i,j]
        n = ijkl[k,l]
        T[i,k] = T[i,k] + C[m,n]*x[j]*x[l]
    end
    return LinearAlgebra.Hermitian(SMatrix(T)) # Real-symmetric
end

"""
    incaz2cart(inc, az) -> v

Return a 3-vector `v` which is the cartesian direction corresponding to
the azimuth `az` and inclination `inc` (degrees).
"""
function incaz2cart(inc, az)
    i = deg2rad(inc)
    a = deg2rad(az)
    cosi = cos(i)
    return @SVector [cos(a)*cosi, -sin(a)*cosi, sin(i)]
end

"""
    cart2incaz(x, y, z) -> inc, az

Return the inclination `inc` and azimuth `az` (both degrees) from
the cartesian direction defined by the vector `[x, y, z]`.
The vector need not be a unit vector as it is normalised inside
this function.
"""
function cart2incaz(x, y, z)
    x̂, ŷ, ẑ = (x, y, z)./hypot(x, y, z)
    inc = asind(ẑ)
    az = atand(-ŷ, x̂)
    return inc, az
end
