module CIJ
# Module containing routines for dealing with elastic constants

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
	phase_vels,
	pitl,
	rot3,
	thom,
	thom_st

function cijkl(C)
	# Convert the 6x6 Voigt matrix into the 3x3x3x3 tensor
	const a = [1 6 5
	           6 2 4
			   5 4 3]
	if ! is_6x6(C)
		error("CIJ.cijkl: Input must be a 6x6 array")
	end
	c = zeros(3,3,3,3)
	for i = 1:3, j = 1:3, k = 1:3, l = 1:3
		c[i,j,k,l] = C[a[i,j],a[k,l]]
	end
	return c
end

function cij(c)
	# Convert the 3x3x3x3 elasticity tensor into the 6x6 Voigt matrix
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

function thom(vp, vs, eps, gam, del)
	# Return a 6x6 Voigt matrix defined by the weak anisotropy parameters of
	# Thomsen (1986) Weak elastic anisotropy.  Geophysics, 51, 10, 1954-1966.
	# Output is density-normalised tensor.
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

function thom_st(vp, vs, eps, gam ,delst)
	# Return a 6x6 Voigt matrix defined by the general anisotropy parameters of
	# Thomsen (1986)
	# Output has same units as input
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
	
function global_VTI(vp, vs, xi, phi, eta)
	# Return a 6x6 Voigt matrix defined by the radial anisotropy parameters
	# as used typically in global seismology
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

function pitl(d1, vp1, vs1, rho1, d2, vp2, vs2, rho2)
	# Return a tuple containing C and rho for a periodic, isotropic thin-layered
	# medium
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

function C2S(C)
	# Return the inverse of a 6x6 matrix
	inv(C)
end

# Copy of C2S for clarity when converting between C and S
S2C = C2S

function is_stable(C)
	# Return false if the input 6x6 matrix is not positive definite (symmetric)
	# and hence not dynamically stable.
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

# Return the normalised elastic constants for olivine, handy for testing purposes
ol() = [320.5  68.1  71.6   0.0   0.0   0.0;
         68.1 196.5  76.8   0.0   0.0   0.0;
		 71.6  76.8 233.5   0.0   0.0   0.0;
		  0.0   0.0   0.0  64.0   0.0   0.0;
		  0.0   0.0   0.0   0.0  77.0   0.0;
		  0.0   0.0   0.0   0.0   0.0  78.7]*1.e9/3355.

function symm(C)
	# Return a copy of C, but with the upper half copied into the lower to enforce
	# symmetry.
	c = zeros(6,6)
	for i = 1:6
		for j = i:6
			c[i,j] = C[i,j]
			if i != j; c[j,i] = C[i,j]; end
		end
	end
	return c
end

function is_symm(C)
	# Return true if C is symmetrical
	for i = 1:6
		for j = i:6
			if C[i,j] != C[j,i]
				return false
			end
		end
	end
	return true
end

function is_6x6(C)
	# Return true is C is a 6x6 matrix
	return size(C) == (6,6)
end

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