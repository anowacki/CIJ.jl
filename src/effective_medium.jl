# Functions related to creating homogenised equivalent medium tensors

"""
    pitl(d1, vp1, vs1, rho1, d2, vp2, vs2, rho2) -> C, rho

Return the effective 6x6 Voigt matrix `C` and density `rho` for a medium defined
by two periodic layers where each layer `i` is defined by:

`di`:   The proportion of the total medium

`vpi`:  P-wave velocity (m/s)

`vsi`:  S-wave velocity (m/s)

`rhoi`: Density (kg/m^3)

`C` is density-normalised.
"""
function pitl(d1, vp1, vs1, rho1, d2, vp2, vs2, rho2)
    C = zero(EC)
    # Lamé parameters from velocities
    m1 = rho1*vs1^2
    m2 = rho2*vs2^2
    l1 = rho1*vp1^2 - 2m1
    l2 = rho2*vp2^2 - 2m2
    # Time-saving terms
    l1p2m1 = l1 + 2*m1
    l2p2m2 = l2 + 2*m2
    # D term, p. 785
    D = (d1 + d2)*(d1*l2p2m2 + d2*l1p2m1)
    # Eq. (7)
    @inbounds begin
        C[1,1] = ((d1+d2)^2*l1p2m1*l2p2m2 + 4*d1*d2*(m1 - m2)*((l1 + m1) - (l2 + m2)))/D
        C[2,2] = C[1,1]
        C[1,2] = C[2,1] = ((d1 + d2)^2*l1*l2 + 2*(l1*d1 + l2*d2)*(m2*d1 + m1*d2))/D
        C[1,3] = ((d1 + d2)*(l1*d1*l2p2m2 + l2*d2*l1p2m1))/D
        C[3,1] = C[2,3] = C[1,3]
        C[3,2] = C[2,3]
        C[3,3] = ((d1 + d2)^2*l1p2m1*l2p2m2)/D
        C[4,4] = C[5,5] = (d1 + d2)*m1*m2/(d1*m2 + d2*m1)
        C[6,6] = (m1*d1 + m2*d2)/(d1 + d2)
    end
    # Normalise back to average density
    rho = (d1*rho1 + d2*rho2)/(d1 + d2)
    C ./= rho
    return C, rho
end

"""
    tandon_and_weng(vp, vs, rho, del, c, vpi, vsi, rhoi) -> C, rho

Return the effective elastic constants `C` and density `rho` of a medium with matrix
velocities `vp` and `vs` (m/s) and density `rho` (kg/m^3), and inclusion
velocities `vpi` and `vsi`, density `rhoi`.  The symmetry axis is parallel to x1.

`del` is the aspect ration of spheroidal inclusions: <1 oblate, >1 prolate

`c` is the volume fraction of inclusions (0 <= c <= 1).
"""
function tandon_and_weng(vp, vs, rho, del, c, vpi, vsi, rhoi)
    #This implementation is based on the Fortran code by Mike Kendall
    del > 0 || error("CIJ.tandon_and_weng: `del` must be > 0.")
    del == 1 && error("CIJ.tandon_and_weng: Theory not valid for `del == 1`")
    0 <= c <= 1 || error("CIJ.tandon_and_weng: `c` must be in range 0 - 1")
    vp == vpi && vs == vsi && rho == rhoi && return CIJ.iso(vp=vp, vs=vs), rho

    C = zero(EC)
    #  weighted average density
    rho_out = (1 - c)*rho + c*rhoi

    amu  = vs^2*rho
    amui = vsi^2*rhoi
    alam = vp^2*rho - 2*amu
    alami = vpi^2*rhoi - 2*amui
    bmi = alami + amui*2/3
    bmps = alam + amu
    #  Young's modulus for matrix
    E0 = amu*(3*alam + 2*amu)/(alam + amu)
    #  Poisson's ratio of the matrix.
    anu = alam/(2*(alam + amu))

    #  Some time saving terms
    t1 = del^2 - 1
    t2 = 1 - anu
    t3 = 1 - 2*anu
    t4 = 3 * del^2
    t5 = 1 - del^2

    # D1, D2 and D3 from Tandon and Weng (1984) (just before equation (18)).
    D1 = 1 + 2*(amui - amu)/(alami - alam)
    D2 = (alam + 2*amu)/(alami - alam)
    D3 = alam/(alami - alam)

    # g and g' terms (appendix of Tandon and Weng 1984). g is for prolate spheroidal
    # inclusions (del>1), whilst g' is for disc-like (oblate) inclusions (del<1).
    if (del >= 1)
        acshdel = log(del + sqrt(t1))
        g = (del*sqrt(t1) - acshdel)*del/sqrt(t1^3)
    else
        # g' below
        g = (acos(del) - del*sqrt(t5))*del/sqrt(t5^3) ;
    end

    # Eshelby's Sijkl tensor (appendix of Tandon and Weng 1984).
    s11 = (t3 + (t4-1)/t1 - (t3 + t4/t1)*g)/(2*t2)
    s22 = (t4/(t1*2) + (t3 - 9/(4*t1))*g)/(4*t2)
    s33 = s22
    s23 = (del^2/(2*t1) - (t3 + 3/(4*t1))*g)/(4*t2)
    s32 = s23
    s21 = (-2*del*del/t1 + (t4/t1 - t3)*g)/(4*t2)
    s31 = s21
    s12 = (-(t3 + 1/t1) + (t3 + 3/(2*t1))*g)/(2*t2)
    s13 = s12
    s44 = (del*del/(2*t1) + (t3 - 3/(4*t1))*g)/(4*t2)
    s66 = (t3 - (t1+2)/t1 - (t3 - 3*(t1+2)/t1)*g/2)/(4*t2)
    s55 = s66

    # Tandon and Weng's B terms (after equation 17).
    B1 = c*D1 + D2 + (1-c)*(D1*s11 + 2*s21)
    B2 = c + D3 + (1-c)*(D1*s12 + s22 + s23)
    B3 = c + D3 + (1-c)*(s11 + (1+D1)*s21)
    B4 = c*D1 + D2 + (1-c)*(s12 + D1*s22 + s23)
    B5 = c + D3 + (1-c)*(s12 + s22 + D1*s23)

    # Tandon and Weng's A terms (after equation 20).
    A1 = D1*(B4 + B5) - 2*B2
    A2 = (1 + D1)*B2 - (B4 + B5)
    A3 = B1 - D1*B3
    A4 = (1 + D1)*B1 - 2*B3
    A5 = (1 - D1)/(B4 - B5)
    A = 2*B2*B3 - B1*(B4+B5)

    # Tandon and Weng (1984) equations (25) (28) (31) (32)
    E11 = E0/(1+c*(A1+2*anu*A2)/A)
    E22 = E0/(1+c*(-2*anu*A3 + (1-anu)*A4 + (1+anu)*A5*A)/(2*A))
    amu12 = amu*(1 + c/(amu/(amui-amu) + 2*(1-c)*s66))
    amu23 = amu*(1 + c/(amu/(amui-amu) + 2*(1-c)*s44))

    # Sayers equation (36)
    anu31 = anu - c*(anu*(A1+2*anu*A2)+(A3-anu*A4)) / (A + c*(A1+2*anu*A2))

    # T&W equation (36)
    #     aK12 term; bmps=plane strain bulk modulus
    anum = (1 + anu)*(1 - 2*anu)
    denom = 1 - anu*(1+2*anu31) + c*(2*(anu31-anu)*A3 + (1-anu*(1+2*anu31))*A4)/A
    aK23 = bmps*anum/denom
    anu12tst = E11/E22 - (1/amu23 + 1/aK23)*E11/4

    # Cij - Sayers' (1992) equations (24)-(29).
    # Conversion
    C[2,2] = amu23 + aK23
    C[3,3] = C[2,2]
    C[1,1] = E11 + 4*anu12tst*aK23
    C[2,3] = -amu23 + aK23
    C[1,2] = 2*anu31*aK23
    C[1,3] = C[1,2]
    C[5,5] = amu12
    C[6,6] = C[5,5]
    C[4,4] = (C[2,2] - C[2,3])/2

    # Fill out matrix by symmetry
    CIJ.symm!(C)

    # apply density normalisation
    C ./= rho_out
    
    C, rho_out
end

"""
    hudson(vp, vs, rho, del, ϵ, vpi, vsi, rhoi) -> C, rho_bulk

Return the effective elastic constants `C` and density `rho_bulk` of a medium with
matrix isotropic velocities `vp` and `vs`, and density `rho`, and inclusion
velocities `vpi` and `vsi`, and density `rhoi`, where the crack density is `ϵ`
and the aspect ratio of ellipsoidal inclusions is `del`.

The theory is valid when `ϵ` ≪ 1, where `ϵ` is given by ``N a^2/\nu`` and ``N`` is
the number of cracks of radius ``a`` in a volume ``\nu``.

The formulation is according to Hudson (1982), as given in Crampin (1984).
"""
function hudson(vp, vs, rho, del, ϵ, vpi, vsi, rhoi)
    # error("`hudson` has not been tested yet")
    ϵ > 0.1 && warn("Theory of Hudson (1982) only valid for `ϵ` < 0.1, but using ϵ = $ϵ")
    λ, μ, κ = lame(vp, vs, rho)
    λ′, μ′, κ′ = lame(vpi, vsi, rhoi)
    K = (κ′ + 4/3*μ′)/(π*del*μ) * (λ + 2μ)/(λ + μ)
    M = 4μ′/(π*del*μ) * (λ + 2μ)/(3λ + 4μ)
    U₁ = 4/3*(λ + 2μ)/(λ + μ)/(1 + K)
    U₃ = 16/3*(λ + 2μ)/(3λ + 4μ)/(1 + M)
    q = 15*(λ/μ)^2 + 28*(λ/μ) + 28
    X = 2μ*(3λ + 8μ)/(λ + 2μ)
    c⁰ = iso(lam=λ, mu=μ)
    # First-order correction from Hudson (1981)
    c¹ = -ϵ./μ .* @SMatrix [ (λ + 2μ)^2*U₁  λ*(λ + 2μ)*U₁  λ*(λ + 2μ)*U₁ 0    0      0
                             λ*(λ + 2μ)*U₁     λ^2*U₁         λ^2*U₁     0    0      0
                             λ*(λ + 2μ)*U₁     λ^2*U₁         λ^2*U₁     0    0      0
                                   0             0              0        0    0      0
                                   0             0              0        0  μ^2*U₃   0
                                   0             0              0        0    0    μ^2*U₃]
    # Second-order correction from Hudson (1982)
    c² = ϵ.^2 ./ 15 .* @SMatrix [ (λ + 2μ)*q*U₁^2       λ*q*U₁^2             λ*q*U₁^2       0   0      0
                                     λ*q*U₁^2     λ^2*q/(λ + 2μ)*U₁^2  λ^2*q/(λ + 2μ)*U₁^2  0   0      0
                                     λ*q*U₁^2     λ^2*q/(λ + 2μ)*U₁^2  λ^2*q/(λ + 2μ)*U₁^2  0   0      0
                                        0                  0                    0           0   0      0
                                        0                  0                    0           0 X*U₃^2   0
                                        0                  0                    0           0   0    X*U₃^2]
    rho_bulk = (1 - ϵ)*rho + ϵ*rhoi
    c = EC{DEFAULT_FLOAT}((c⁰ .+ c¹ .+ c²)./rho_bulk) 
    c, rho_bulk
end

"""
    grechka_cracks!(C, ξ, ϕ=0) -> C

Add dry cracks using the theory of Grechka (2007) [1] to a tensor `C` which has VTI symmetry
about the 3-axis.  Cracks have fracture density `ξ` and strike `ϕ`°, measured from the
1-axis towards the negative 2-axis.  `C` is updated in-place.

    grechka_cracks(C, ξ, ϕ=0) -> C′

Copying version of `grechka_cracks!`.

Both forms take `ξ` and `ϕ` as both scalars and `AbstractArray`s.  In the latter case,
multiple fracture sets are added, and are assumed to be non-interacting.  This is valid for
'small' ξ only.

N.B.: No check is made on the input constants `C`.

#### References

1. Vladimir Grechka (2007). Multiple cracks in VTI rocks: Effective properties and fracture
   characterization.  Geophysics, 72, D81–D91.  doi:10.1190/1.2751500
"""
function grechka_cracks!(C, ξ, ϕ=zero(eltype(C)))
    C .= CIJ.rot3(C, 0, 0, -ϕ)
    # Excess normal and tangential compliances of crack
    Bn = 4/3*ξ*C[2,2]/C[6,6]/(C[2,2] - C[6,6])
    Bth = 16/3*ξ*C[2,2]/C[6,6]/(3C[2,2] - 2C[6,6])
    Btv = 16/3*ξ*C[3,3]/C[5,5]/(3C[3,3] - 2C[5,5])
    # Update compliance after Eschelby (1957)
    S = C2S!(C)
    S[2,2] += Bn
    S[4,4] += Btv
    S[6,6] += Bth
    C .= CIJ.rot3(S2C!(S), 0, 0, ϕ)
end
grechka_cracks(C, ξ, ϕ=zero(eltype(C))) = grechka_cracks!(deepcopy(C), ξ, ϕ)

function grechka_cracks!(C, ξ::AbstractArray, ϕ::AbstractArray=zeros(ξ, eltype(C)))
    length(ξ) == length(ϕ) || throw(ArgumentError("Lengths of ξ and ϕ must be the same"))
    for i in 1:length(ξ)
        grechka_cracks!(C, ξ[i], ϕ[i])
    end
    C
end
grechka_cracks(C, ξ::AbstractVector, ϕ::AbstractVector=zeros(ξ, eltype(C))) =
    grechka_cracks!(deepcopy(C), ξ, ϕ)

@doc (@doc grechka_cracks!) grechka_cracks
