"""
`module CIJ`

CIJ contains routines for dealing with elastic constants.
"""
module CIJ

import Base.write

export
	Au,
	C2S,
	HillG,
	HillK,
	Panning_VTI,
	ReussG,
	ReussK,
	VRH,
	VoigtG,
	VoigtK,
	cij,
	cijkl,
	global_VTI,
	is_stable,
	iso,
	phase_vels,
	pitl,
	rot3,
    symm,
    symm!,
    tandon_and_weng,
	thom,
	thom_st

"""
`iso(;vp=nothing, vs=nothing, lam=nothing, mu=nothing, K=nothing, G=nothing) -> C`

Return an isotropic Voigt stiffness matrix `C` from a pair of the following:

`vp`, `vs` : Isotropic velocities in m/s

`lam`, `mu`: Lamé parameters divided by density in m^2/s^2

`K`, `G`   : Bulk and shear moduli dvided by density in m^2/s^2
"""
function iso(;vp=nothing, vs=nothing, lam=nothing, mu=nothing, K=nothing, G=nothing)
	C = zeros(6,6)
	if vp != nothing && vs != nothing
		C[1,1] = vp^2
		C[4,4] = vs^2
		C[1,2] = C[1,1] - 2C[4,4]
	elseif lam != nothing && mu != nothing
		C[1,1] = lam + 2mu
		C[4,4] = mu
		C[1,2] = lam
	elseif K != nothing && G != nothing
		C[1,1] = K + 4G/3.
		C[4,4] = G
		C[1,2] = C[1,1] - 2C[4,4]
	else
		error("iso: Must request a pair of vp,vs; lam,mu; K,G")
	end
	C[2,2] = C[3,3] = C[1,1]
	C[5,5] = C[6,6] = C[4,4]
	C[1,3] = C[2,3] = C[3,1] = C[3,2] = C[2,1] = C[1,2]
	return C
end

"""
`cijkl(C) -> c`

Convert the 6x6 Voigt matrix `C` into the 3x3x3x3 tensor `c`
"""
function cijkl(C)
	const a = [1 6 5
	           6 2 4
			   5 4 3]
	if ! is_6x6(C)
		error("CIJ.cijkl: Input must be a 6x6 array")
	end
	c = Array{eltype(C)}(3, 3, 3, 3)
	for i = 1:3, j = 1:3, k = 1:3, l = 1:3
		c[i,j,k,l] = C[a[i,j],a[k,l]]
	end
	return c
end

"""
`cij(c) -> C`

Convert the 3x3x3x3 elasticity tensor `c` into the 6x6 Voigt matrix `C`
"""
function cij(c)
	const a = [1 6 5
	           6 2 4
			   5 4 3]
	if size(c) != (3,3,3,3)
		error("CIJ.cij: Input must be a 3x3x3x3 tensor")
	end
	C = zeros(6,6)
	for i = 1:3, j = 1:3, k = 1:3, l = 1:3
		C[a[i,j],a[k,l]] = c[i,j,k,l]
	end
	return C
end

"""
`thom(vp, vs, eps, gam, del) -> C`

Return the 6x6 Voigt matrix defined `C` by the weak anisotropy parameters of
Thomsen (1986) Weak elastic anisotropy.  Geophysics, 51, 10, 1954-1966.

Output is density-normalised tensor.
"""
function thom(vp, vs, eps, gam, del)
	if vp <= 0.
		error("CIJ.thom: vp must be greater than 0")
	elseif vs < 0.
		error("CIJ.thom: vs must be greater than or equal to 0")
	end

	C = zeros(6,6)
	C[3,3] = vp^2
	C[1,1] = C[2,2] = C[3,3]*(2.*eps + 1.)
	C[4,4] = C[5,5] = vs^2
	C[6,6] = C[4,4]*(2.*gam + 1.)

	b = 2.*C[4,4]
	term = C[3,3] - C[4,4]
	c = C[4,4]^2 - (2.*del*C[3,3]*term + term^2)
	d = b^2 - 4.*c
	if (d < 0.)
		error("CIJ.thom: S-velocity too high or delta too negative")
	end
	C[1,3] = C[2,3] = C[3,1] = C[3,2] = -b/2. + sqrt(d)/2.
	C[1,2] = C[2,1] = C[1,1] - 2.*C[6,6]
	return C
end

"""
`thom_st(vp, vs, eps, gam, delst) -> C`

Return the 6x6 Voigt matrix `C` defined by the general anisotropy parameters of
Thomsen (1986).

Output has same units as input.
"""
function thom_st(vp, vs, eps, gam, delst)
	if vp <= 0.
		error("CIJ.thom: vp must be greater than 0")
	elseif vs < 0.
		error("CIJ.thom: vs must be greater than or equal to 0")
	end

	C = zeros(6,6)
	C[3,3] = vp^2
	C[1,1] = C[2,2] = C[3,3]*(2.*eps + 1.)
	C[4,4] = C[5,5] = vs^2
	C[6,6] = C[4,4]*(2.*gam + 1.)
	a = 2.
	b = 4.*C[4,4]
	c = C[4,4]^2 - 2.delst*C[3,3]^2 - (C[3,3] - C[4,4])*(C[1,1] + C[3,3] - 2.*C[4,4])
	if b^2 - 4.*a*c < 0.
		error("CIJ.thom_st: S velocity too high or delta too negative")
	end
	C[1,3] = C[3,1] = C[2,3] = C[3,2] = (-b + sqrt(b^2 - 4.*a*c))/(2.*a)
	C[1,2] = C[2,1] = C[1,1] - 2.*C[6,6]
	return C
end

"""
`global_VTI(vp, vs, xi, phi, eta) -> C`

Return the 6x6 Voigt matrix `C` defined by the radial anisotropy parameters
as used typically in global seismology.
"""
function global_VTI(vp, vs, xi, phi, eta)
	if vp <= 0.
		error("CIJ.global_VTI: vp must be greater than 0")
	elseif vs < 0.
		error("CIJ.global_VTI: vs must be greater than or equal to 0")
	elseif xi <= 0. || phi <= 0. || eta <= 0.
		error("CIJ.global_VTI: xi, phi and eta must be greater then 0")
	end

	L = 15.*((3.*phi + 8. + 4.*eta)*vs^2 - (phi + 1. - 2.*eta)*vp^2) /
		((6. + 4.*eta + 5.*xi)*(3.*phi + 8. +4.*eta) - 8.*(phi + 1. -2.*eta)*(1. - eta))
	A = (15.*vp^2 - 8.*(1. - eta)*L) / (3.*phi + 8. + 4.*eta)
	F = eta*(A - 2.*L)
	C = phi*A
	N = xi*L
	C12 = A - 2.*N
	return [ A  C12 F  0. 0. 0.
	        C12  A  F  0. 0. 0.
			 F   F  C  0. 0. 0.
			 0.  0. 0. L  0. 0.
			 0.  0. 0. 0. L  0.
			 0.  0. 0. 0. 0. N]
end

"""
`global_VTI(vp, vs, xi, phi, eta) -> C`

Return the 6x6 Voigt matrix `C` defined by the radial anisotropy parameters
defined by Panning and Romanowicz as an approximation to general VTI.
"""
function Panning_VTI(vp, vs, xi, phi)
	L = 3.*vs^2/(2. + xi)
	N = xi*L
	A = 5.*vp^2/(4. + phi)
	C = phi*A
	F = A - 2.*L
	C12 = A - 2.*N
	return [ A  C12 F  0. 0. 0.
	        C12  A  F  0. 0. 0.
			 F   F  C  0. 0. 0.
			 0.  0. 0. L  0. 0.
			 0.  0. 0. 0. L  0.
			 0.  0. 0. 0. 0. N]
end

"""
`pitl(d1, vp1, vs1, rho1, d2, vp2, vs2, rho2) -> C, rho`

Return the effective 6x6 Voigt matrix `C` and density `rho` for a medium defined
by two periodic layers where each layer `i` is defined by:

`di`:   The proportion of the total medium

`vpi`:  P-wave velocity (m/s)

`vsi`:  S-wave velocity (m/s)

`rhoi`: Density (kg/m^3)

`C` is density-normalised.
"""
function pitl(d1, vp1, vs1, rho1, d2, vp2, vs2, rho2)
	C = zeros(6,6)
	# Lamé parameters from velocities
	m1 = rho1*vs1^2
	m2 = rho2*vs2^2
	l1 = rho1*vp1^2 - 2.*m1
	l2 = rho2*vp2^2 - 2.*m2
	# Time-saving terms
    l1p2m1 = l1 + 2.*m1
    l2p2m2 = l2 + 2.*m2
    # D term, p. 785
    D = (d1 + d2)*(d1*l2p2m2 + d2*l1p2m1)
    # Eq. (7)
    C[1,1] = ((d1+d2)^2*l1p2m1*l2p2m2 + 4.*d1*d2*(m1 - m2)*((l1 + m1) - (l2 + m2)))/D
    C[2,2] = C[1,1]
    C[1,2] = C[2,1] = ((d1 + d2)^2*l1*l2 + 2.*(l1*d1 + l2*d2)*(m2*d1 + m1*d2))/D
    C[1,3] = ((d1 + d2)*(l1*d1*l2p2m2 + l2*d2*l1p2m1))/D
    C[3,1] = C[2,3] = C[1,3]
    C[3,2] = C[2,3]
    C[3,3] = ((d1 + d2)^2*l1p2m1*l2p2m2)/D
    C[4,4] = C[5,5] = (d1 + d2)*m1*m2/(d1*m2 + d2*m1)
    C[6,6] = (m1*d1 + m2*d2)/(d1 + d2)
    # Normalise back to average density
	rho = (d1*rho1 + d2*rho2)/(d1 + d2)
    C /= rho
	return (C, rho)
end

"""
`tandon_and_weng(vp, vs, rho, del, c, vpi, vsi, rhoi) -> C, rho`

Return the effective elastic constants `C` and density `rho` of a medium with matrix
velocities `vp` and `vs` (m/s) and density `rho` (kg/m^3), and inclusion
velocities `vpi` and `vsi`, density `rhoi`.

`del` is the aspect ration of spheroidal inclusions: <1 oblate, >1 prolate

`c` is the volume fraction of inclusions (0 <= c <= 1).
"""
function tandon_and_weng(vp, vs, rho, del, c, vpi, vsi, rhoi)
    #This implementation is based on the Fortran code by Mike Kendall
    del > 0 || error("CIJ.tandon_and_weng: `del` must be > 0.")
    del == 1 && error("CIJ.tandon_and_weng: Theory not valid for `del == 1`")
    0 <= c <= 1 || error("CIJ.tandon_and_weng: `c` must be in range 0 - 1")
    vp == vpi && vs == vsi && rho == rhoi && return CIJ.iso(vp=vp, vs=vs), rho

    C = zeros(6,6)
    #  weighted average density
    rho_out = (1.0 - c)*rho + c*rhoi

    amu  = vs^2*rho
    amui = vsi^2*rhoi
    alam = vp^2*rho - 2.0*amu
    alami = vpi^2*rhoi - 2.0*amui
    bmi = alami + amui*2.0/3.0
    bmps = alam + amu
    #  Young's modulus for matrix
    E0 = amu*(3.0*alam + 2.0*amu)/(alam + amu)
    #  Poisson's ratio of the matrix.
    anu = alam/(2.0*(alam + amu))

    #  Some time saving terms
    t1 = del^2 - 1.0
    t2 = 1.0 - anu
    t3 = 1.0 - 2.0*anu
    t4 = 3.0 * del^2
    t5 = 1.0 - del^2

    # D1, D2 and D3 from Tandon and Weng (1984) (just before equation (18)).
    D1 = 1.0 + 2.0*(amui - amu)/(alami - alam)
    D2 = (alam + 2.0*amu)/(alami - alam)
    D3 = alam/(alami - alam)

    # g and g' terms (appendix of Tandon and Weng 1984). g is for prolate spheroidal
    # inclusions (del>1), whilst g' is for disc-like (oblate) inclusions (del<1).
    if (del >= 1) then
        acshdel = log(del + sqrt(t1))
        g = (del*sqrt(t1) - acshdel)*del/sqrt(t1^3)
    else
        # g' below
        g = (acos(del) - del*sqrt(t5))*del/sqrt(t5^3) ;
    end

    # Eshelby's Sijkl tensor (appendix of Tandon and Weng 1984).
    s11 = (t3 + (t4-1.0)/t1 - (t3 + t4/t1)*g)/(2.0*t2)
    s22 = (t4/(t1*2.0) + (t3 - 9.0/(4.0*t1))*g)/(4.0*t2)
    s33 = s22
    s23 = (del^2/(2.0*t1) - (t3 + 3.0/(4.0*t1))*g)/(4.0*t2)
    s32 = s23
    s21 = (-2.0*del*del/t1 + (t4/t1 - t3)*g)/(4.0*t2)
    s31 = s21
    s12 = (-1.0*(t3 + 1.0/t1) + (t3 + 3.0/(2.0*t1))*g)/(2.0*t2)
    s13 = s12
    s44 = (del*del/(2.0*t1) + (t3 - 3.0/(4.0*t1))*g)/(4.0*t2)
    s66 = (t3 - (t1+2.0)/t1 - (t3 - 3.0*(t1+2.0)/t1)*g/2.0)/(4.0*t2)
    s55 = s66

    # Tandon and Weng's B terms (after equation 17).
    B1 = c*D1 + D2 + (1.0-c)*(D1*s11 + 2.0*s21)
    B2 = c + D3 + (1.0-c)*(D1*s12 + s22 + s23)
    B3 = c + D3 + (1.0-c)*(s11 + (1.0+D1)*s21)
    B4 = c*D1 + D2 + (1.0-c)*(s12 + D1*s22 + s23)
    B5 = c + D3 + (1.0-c)*(s12 + s22 + D1*s23)

    # Tandon and Weng's A terms (after equation 20).
    A1 = D1*(B4 + B5) - 2.0*B2
    A2 = (1.0 + D1)*B2 - (B4 + B5)
    A3 = B1 - D1*B3
    A4 = (1.0 + D1)*B1 - 2.0*B3
    A5 = (1.0 - D1)/(B4 - B5)
    A = 2.0*B2*B3 - B1*(B4+B5)

    # Tandon and Weng (1984) equations (25) (28) (31) (32)
    E11 = E0/(1.0+c*(A1+2.0*anu*A2)/A)
    E22 = E0/(1.0+c*(-2.0*anu*A3 + (1.0-anu)*A4 + (1.0+anu)*A5*A)/(2.0*A))
    amu12 = amu*(1.0 + c/(amu/(amui-amu) + 2.0*(1.0-c)*s66))
    amu23 = amu*(1.0 + c/(amu/(amui-amu) + 2.0*(1.0-c)*s44))

    # Sayers equation (36)
    anu31 = anu - c*(anu*(A1+2.0*anu*A2)+(A3-anu*A4)) / (A + c*(A1+2.0*anu*A2))

    # T&W equation (36)
    #     aK12 term; bmps=plane strain bulk modulus
    anum = (1.0+anu)*(1.0-2.0*anu)
    denom = 1.0 - anu*(1.0+2.0*anu31) + c*(2.0*(anu31-anu)*A3 + (1.0-anu*(1.0+2.0*anu31))*A4)/A
    aK23 = bmps*anum/denom
    anu12tst = E11/E22 - (1.0/amu23 + 1.0/aK23)*E11/4.0

    # Cij - Sayers' (1992) equations (24)-(29).
    # Conversion
    C[2,2] = amu23 + aK23
    C[3,3] = C[2,2]
    C[1,1] = E11 + 4.0*anu12tst*aK23
    C[2,3] = -amu23 + aK23
    C[1,2] = 2.0*anu31*aK23
    C[1,3] = C[1,2]
    C[5,5] = amu12
    C[6,6] = C[5,5]
    C[4,4] = (C[2,2] - C[2,3])/2.0

    # Fill out matrix by symmetry
    CIJ.symm!(C)

    # apply density normalisation
    C/rho_out, rho_out
end

"""
`phase_vels(C, az, inc) -> vp, vs1, vs2, pol, avs`

Calculate the phase velocities for the 6x6 elasticity matrix `C` along the direction
(`az`, `inc`), in degrees, and return P-wave velocity `vp`, the fast and slow shear
wave velocities, `vs1` and `vs2`, the polarisation of the fast shear wave `pol`,
and the shear wave velocity anisotropy, `avs`.  Velocities are in m/s if the tensor
`C` is in m^2/s^2 (i.e., is a density-normalised tensor, sometimes called A).

`az` is the azimuth in degrees measured from the x1 towards to -x2 axis.

`inc` is the inclination in degrees from the x1-x2 plane towards the x3 axis.
"""
function phase_vels(C, az, inc)
	# Return the phase velocities in the az, inc direction for a 6x6
	# Voigt elasticity matrix C.
	# INPUT:
	#	C(6,6)  : Elasticity matrix
	#	az      : Azimiuth (degrees) away from x1 towards -x2 axis
	#	inc     : Inclination (degrees) away from x1-x2 plane towards x3 axis

	x = incaz2cart(inc, az)
	# Create the Christoffel matrix
	T = make_T(C, x)
	# Find eigenvectors and values, which correspond to velocities^2 (val[i])
	# and polarisation vectors (vec[:,i]).
	# They seem to be sorted into descending order of eigenvalue, but we double
	# check anyway.
	eval, evec = eig(T)
	ip = indmax(eval)
	is2 = indmin(eval)
	is1 = 6 - ip - is2
	vp  = sqrt(eval[ip])
	vs1 = sqrt(eval[is1])
	vs2 = sqrt(eval[is2])
	# xp  = vec(evec[:,ip])
	xs1 = vec(evec[:,is1])
	# xs2 = vec(evec[:,is2])
	# Calculate S1 polarisation and amount of shear wave anisotropy
	pol = get_pol(inc, az, vec(x), vec(xs1))
	avs = 200.*(vs1 - vs2)/(vs1 + vs2)
	return (vp, vs1, vs2, pol, avs)
end

function group_vels(C, az, inc)
	# Return the group velocities in the az, inc direction for a 6x6
	# Voigt elasticity matrix C, in an elastic medium.
	# INPUT:
	#	C(6,6)  : Elasticity matrix
	#	az      : Azimuth (degrees) away from x1 towards -x2 axis
	#	inc     : Inclination (degrees) away from x1-x2 plane towards x3

	x = incaz2cart(inc, az)
	# Create the Christoffel matrix
	T = make_T(C, x)
	# Find eigenvectors, giving phase velocities
	eval, evec = eig(T)
	ip = indmax(eval)
	is2 = indmin(eval)
	is1 = 6 - ip - is2
	vp = sqrt(eval[ip])
	vs1 = sqrt(eval[is1])
	vs2 = sqrt(eval[is2])
	xp = vec(evec[:,ip])
	xs1 = vec(evec[:,is1])
	xs2 = vec(evec[:,is2])
	# Calculate group velocities, which are phase velocities projected onto
	# group velocity directions
	pp = xp/vp
	ps1 = xs1/vs1
	ps2 = xs2/vs2
	p = [pp'; ps1'; ps2']
	vg = zeros(3)
	const ijkl = [1 6 5; 6 2 4; 5 4 3]
	for i=1:3, j=1:3, k=1:3, l=1:3
		m = ijkl[i,j]
		n = ijkl[k,l]
		vg[i] = vg[i] + C[m,n]*p[i,k]*x[j]*x[l]
	end
	return vg[1], vg[2], vg[3]
end

function get_pol(inc, az, x, xs1)
	# Projection of fast shear wave polarisation onto wavefront plane
	xs1p = cross(x, cross(x, xs1))
	xs1p /= norm(xs1p)
	# Local up vector
	u = incaz2up(inc, az)
	# Angle between the local up vector and the projected fast orientation
	v = cross(xs1p, u)
	pol = rad2deg(atan2(norm(v), dot(xs1p, u)))
	# If v is codirectional with x, then the angle is correct; otherwise, we're
	# the wrong way round
	if dot(x, v) < 0.
		pol = -pol
	end
	# Put in range -90° to 90°
	return modulo(pol + 90., 180.) - 90.
end

function incaz2up(inc, az)
	a = deg2rad(az) + pi
	i = deg2rad(inc)
	sini = sin(i)
	return [cos(a)*sini, -sin(a)*sini, cos(i)]
end

function make_T(C, x)
	# Create the Christoffel matrix, T
	# INPUT:
	#	C(6,6) : Voigt elasticity matrix
	#	x(3)   : Direction of interest
	T = zeros(3,3)
	const ijkl = [1 6 5; 6 2 4; 5 4 3]
	for i = 1:3, j = 1:3, k = 1:3, l = 1:3
		m = ijkl[i,j]
		n = ijkl[k,l]
		T[i,k] = T[i,k] + C[m,n]*x[j]*x[l]
	end
	return T
end

function incaz2cart(inc, az)
	# Return a 3-vector which is the cartesian direction corresponding to
	# the direction az, inc (degrees)
	i = deg2rad(inc)
	a = deg2rad(az)
	cosi = cos(i)
	return [cos(a)*cosi, -sin(a)*cosi, sin(i)]
end

function VRH(VF1, C1, rh1, VF2, C2, rh2)
	# Return the VRH elasticity and density
	voigt = (C1.*VF1 + C2.*VF2)/(VF1 + VF2)
	reuss = (C2S(C1).*VF1 + C2S(C2).*VF2)/(VF1 + VF2)
	return ((voigt + S2C(reuss))/2., (rh1*VF1 + rh2*VF2)/(VF1 + VF2))
end

@doc doc"""
`VRH(VF, C, rho) -> C, rho`

Return the Voigt-Reus-Hill-averaged 6x6 Voigt stiffness matrix `C` and density `rho`
from combining n sets of elastic constants and densities, whose proportions are
listed in `VF`, where n is the number of constants and length of `VF` and `rho`.

`VF`: A vector (which need not sum to 1) containing the relative proportions of each
set of elastic constants, size n

`C`: A nx6x6 array of elastic constants

`rho`: A vector of densities of length n.
""" ->
function VRH(VF, C, rho)
    size(VF) == size(C,1) == size(rho) || error("CIJ.VRH: Arrays are not the same lengths")
    size(C)[2:3] == (6,6) || error("CIJ.VRH: C must have size (n,6,6)")
    VF = VF./sum(VF)
    voigt = zero(C)
    reuss = zero(C)
    Cave = zeros(6, 6)
    for k in 1:length(VF)
        voigt += VF[k]*reshape(C[k,:,:], 6, 6)
        reuss += VF[k]*C2S(reshape(C[k,:,:], 6, 6))
    end
    reuss = S2C(reuss)
    Cave = (voigt + reuss)/2
    rhave = VF.*rho
    Cave, rhave
end

function VoigtK(C)
    # Return the Voigt bound on a crystal's K modulus.
	# See: Hill, R., The elastic behaviour of a crystalline aggregate,
	#	   P Phys Soc Lond A (1952) vol. 65 (389) pp. 349-355
	(C[1,1] + C[2,2] + C[3,3] + 2.*(C[1,2] + C[2,3] + C[1,3]))/9.
end

function VoigtG(C)
    # Return the Voigt bound on a crystal's G modulus
	(C[1,1] + C[2,2] + C[3,3] - (C[1,2] + C[2,3] + C[1,3]) +
		3.*(C[4,4] + C[5,5] + C[6,6]))/15.
end

function ReussK(C)
    # Return the Reuss bound on a crystal's K modulus
	S = C2S(C)
	1./(S[1,1] + S[2,2] + S[3,3] + 2.*(S[1,2] + S[2,3] + S[1,3]))
end

function ReussG(C)
    # Return the Reuss bound on a crystal's G modulus
	S = C2S(C)
	15./(4.*(S[1,1] + S[2,2] + S[3,3]) - 4.*(S[1,2] + S[2,3] + S[1,3]) +
		3.*(S[4,4] + S[5,5] + S[6,6]))
end

function HillK(C)
	(VoigtK(C) + ReussK(C))/2.
end

function HillG(C)
    (VoigtG(C) + ReussG(C))/2.
end

function Au(C)
    # Return the Universal Elastic Anisotropy Index.
	# See: Ranganathan & Ostoja-Starzewksi, Universal elastic anisotropy index,
	# Phys Rev Lett (2008) vol. 101 (5) pp. 055504
	Kv, Gv, Kr, Gr = VoigtK(C), VoigtG(C), ReussK(C), ReussG(C)
	5.*(Gv/Gr) + Kv/Kr - 6.
end

function rotmatA(a)
	# Return a 3x3 array representing a rotation about the 1 axis by a degrees,
	# clockwise when looking down towards origin along axis.
	da = deg2rad(a)
	return [1.       0.       0.;
	        0.  cos(da)  sin(da);
	        0. -sin(da)  cos(da)]
end

function rotmatB(b)
	# Return a 3x3 array representing a rotation about the 1 axis by a degrees,
	# clockwise when looking down towards origin along axis.
	db = deg2rad(b)
	return [cos(db) 0. -sin(db);
	             0. 1.       0.;
	        sin(db) 0.  cos(db)]
end

function rotmatC(c)
	# Return a 3x3 array representing a rotation about the 1 axis by a degrees,
	# clockwise when looking down towards origin along axis.
	dc = deg2rad(c)
	return [ cos(dc) sin(dc) 0.;
	        -sin(dc) cos(dc) 0.;
	              0.      0. 1.]
end

function rot3(C, a, b, c)
	# Rotate a 6x6 array in turn by the x1, x2 and x3 axes by a, b and c degrees
	R = *(rotmatC(c), *(rotmatB(b), rotmatA(a)))
	return transform(C, R)
end

function transform(C, M)
	# Apply a transformation matrix M to a 6x6 matrix C
	K = zeros(6, 6)
	for i = 1:3
		i1 = (i + 1)%3
		i1 = i1 == 0 ? 3 : i1
		i2 = (i + 2)%3
		i2 = i2 == 0 ? 3 : i2
		for j = 1:3
			j1 = (j + 1)%3
			j1 = j1 == 0 ? 3 : j1
			j2 = (j + 2)%3
			j2 = j2 == 0 ? 3 : j2
			K[i,j] = M[i,j]^2
			K[i,j+3] = 2.*M[i,j1]*M[i,j2]
			K[i+3,j] = M[i1,j]*M[i2,j]
			K[i+3,j+3] = M[i1,j1]*M[i2,j2] + M[i1,j2]*M[i2,j1]
		end
	end
	return K * (C * transpose(K))
end

"""
`C2S(C) -> S`

Return the inverse of the stiffness matrix `C`, the compliance matrix `S`.
"""
C2S(C) = inv(C)

"""
`S2C(S) -> C`

Return the inverse of the compliance matrix `S`, the stiffness matrix `C`.
"""
S2C = C2S

"""
`is_iso(C) -> ::Bool`

Return `true` if `C` is isotropic.
"""
function is_iso(C)
	is_stable(C) || error("CIJ.is_iso: Not a valid elasticity matrix")
	all(abs((diag(C[2:3,2:3]) - C[1,1])) .< eps()) &&
		all(abs(diag(C[5:6,5:6]) - C[4,4]) .< eps()) &&
		all(abs([C[1,2] C[1,3] C[2,3]] - (C[1,1] - 2C[4,4])) .< eps()) &&
		return true
	return false
end

"""
`is_stable(C) -> ::Bool`

Return `false` if the input 6x6 matrix `C` is not positive definite (symmetric)
and hence not dynamically stable.
"""
function is_stable(C)
	if ! is_6x6(C)
		error("CIJ.is_stable: Matrix must be 6x6")
	end
	if ! is_symm(C)
		warn("CIJ.is_stable: Matrix is not symmetrical: Using upper half only")
	end
	# Positive definite matrices have a Cholesky decomposition
	try
		chol(C)
	catch PosDefException
		return false
	end
	return true
end

"""
`ol() -> C, rho`

Return the normalised elastic constants `C` and density `rho` for olivine,
handy for testing purposes
"""
ol() = [320.5  68.1  71.6   0.0   0.0   0.0;
         68.1 196.5  76.8   0.0   0.0   0.0;
		 71.6  76.8 233.5   0.0   0.0   0.0;
		  0.0   0.0   0.0  64.0   0.0   0.0;
		  0.0   0.0   0.0   0.0  77.0   0.0;
		  0.0   0.0   0.0   0.0   0.0  78.7]*1.e9/3355., 3355.


"""
`read(file) -> C, rho`

Return the elasticity matrix `C` and density `rho` from the ecs-format file
`file` (which usually has the extension `.ecs`).  This should be in the format:

    1 1 C11
    . . ...
    i j Cij
    . . ...
    6 6 C66
    7 7 rho

Comments are permitted at any point, preceded by a '#' character.

In this version of the code, all elastic constants must be specified, even if
they are 0.
"""
function read(file::AbstractString)
    isfile(file) || error("CIJ.read: File \"$file\" does not exist.")
    d = readdlm(file)
    size(d, 1) == 22 || error("CIJ.read: File \"$file\" is not in expected format")
    C = zeros(6,6)
    rho = 0
    for l in 1:22
        i = Int(d[l,1])
        j = Int(d[l,2])
        if i == j == 7
            rho = d[l,3]
        else
            1 <= i <= j <= 6 || error("CIJ.read: Unexpected indices '$i' and '$j' in file \"$file\"")
            C[i,j] = d[l,3]
            C[j,i] = d[l,3]
        end
    end
    C, rho
end

"""
`write(C, rho, file, comment="")`

Write the density-normalised elastic constants `C` and density `rho` to the file
`file` in the following format:

    1 1 C11
    . . ...
    i j Cij
    . . ...
    6 6 C66
    7 7 rho

If supplied, `comment` should be a string of one or more lines to write at the
end of the file to describe the constants.  By default, the following is written:

    # Saved by <user> on <hostname> on <date> using <Julia version>

Such files usually have a `.ecs` file extension.
"""
function write{T<:Number}(C::Array{T,2}, rho::Number, file::AbstractString, comment::AbstractString="")
    is_6x6(C) && is_stable(C) && is_symm(C) ||
        error("CIJ.save: Elastic constants are not in the right form.  "
              * "May be asymmetric, unstable or not a 6x6 matrix")
    isdir(dirname(file)) || error("CIJ.save: Directory \"$(dirname(file))\" does not exist")
    f = open(file, "w")
    for i = 1:6, j=i:6
        @printf(f, "%d %d %12.6e\n", i, j, C[i,j]*rho)
    end
    @printf(f, "7 7 %9.3f\n", rho)
    if comment != ""
        if comment[1] == '#' comment = "# " * comment end
        write(f, replace(chomp(comment), "\n", "\n# ") * "\n")
    end
    user = ENV["USER"]
    hostname = readchomp(`hostname`)
    write(f, "# Saved by user $user on $hostname on $(now()) using Julia $VERSION\n")
    close(f)
end

"""
`symm(C) -> C`

Return a copy of C with the upper half copied into the lower half to enforce symmetry.
"""
function symm(C)
	c = copy(C)
	CIJ.symm!(c)
	return c
end

"""
`symm!(C)`

Fill in the lower half of 6x6 matrix `C` with the upper half, making it symmetrical.
"""
symm!(C) = for i = 1:6, j = i+1:6 C[j,i] = C[i,j] end

"""
`is_symm(C) -> ::Bool`

Return `true` if `C` is symmetrical.
"""
function is_symm(C)
	# Return true if C is symmetrical
	for i = 1:6, j = i+1:6
		if C[i,j] != C[j,i] return false end
	end
	return true
end

"""
`is_6x6(C) -> ::Bool`

Return `true` if C is a 6x6 matrix.
"""
is_6x6(C) = size(C) == (6,6)

function modulo(a, b)
	# Replicate the Fortran modulo function, where the result is always
	# positive
	m = a%b
	while m < 0
		m += b
	end
	return m
end

end # module CIJ