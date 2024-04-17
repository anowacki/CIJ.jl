# Calculate velocities for tensors

"""
    phase_vels(C, az, inc) -> vp, vs1, vs2, pol, avs

Calculate the phase velocities for the 6x6 elasticity matrix `C` along the direction
(`az`, `inc`), in degrees, and return P-wave velocity `vp`, the fast and slow shear
wave velocities, `vs1` and `vs2`, the polarisation of the fast shear wave `pol`,
and the shear wave velocity anisotropy, `avs`.  Velocities are in m/s if the tensor
`C` is in m^2/s^2 (i.e., is a density-normalised tensor, sometimes called A).

`az` is the azimuth in degrees measured from the x1 towards the -x2 axis.

`inc` is the inclination in degrees from the x1-x2 plane towards the x3 axis.

`pol` is measured when looking towards the (0,0,0) point along the ray (in
the negative propagation direction).  `pol` increases clockwise away from the vector
normal to the propagation direction which points towards to the x3-axis.
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
    vp  = sqrt(vals[ip])
    vs1 = sqrt(vals[is1])
    vs2 = sqrt(vals[is2])
    # xp  = vec(vecs[:,ip])
    xs1 = vecs[:,is1]
    # xs2 = vec(vecs[:,is2])
    # Calculate S1 polarisation and amount of shear wave anisotropy
    pol = get_pol(inc, az, x, xs1)
    avs = 200*(vs1 - vs2)/(vs1 + vs2)
    return (vp=vp, vs1=vs1, vs2=vs2, pol=pol, avs=avs)
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
