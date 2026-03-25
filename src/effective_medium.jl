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
    tandon_and_weng(vp, vs, ρ, α, c, vpᵢ, vsᵢ, ρᵢ) -> C, ρ

Return the effective elastic constants `C` and density `ρ` of a medium with matrix
velocities `vp` and `vs` (m/s) and density `ρ` (kg/m^3), and inclusion
velocities `vpᵢ` and `vsᵢ`, density `ρᵢ`.  The symmetry axis is parallel to x1.

`α` is the aspect ration of spheroidal inclusions: <1 oblate, >1 prolate

`c` is the volume fraction of inclusions (0 <= c <= 1).
"""
function tandon_and_weng(vp, vs, ρ, α, c, vpᵢ, vsᵢ, ρᵢ)
    # This implementation is based on the Fortran code by Mike Kendall
    α > 0 || error("CIJ.tandon_and_weng: `α` must be > 0.")
    α == 1 && error("CIJ.tandon_and_weng: Theory not valid for `α == 1`")
    0 <= c <= 1 || error("CIJ.tandon_and_weng: `c` must be in range 0 - 1")
    vp == vpᵢ && vs == vsᵢ && ρ == ρᵢ && return CIJ.iso(vp=vp, vs=vs), ρ

    C = zero(EC)

    # Weighted average density
    ρ_out = (1 - c)*ρ + c*ρᵢ

    # Lamé parameters
    # Avoid integer problems by explicit conversion to float
    μ₀ = float(vs)^2*ρ
    μᵢ = float(vsᵢ)^2*ρᵢ
    λ₀ = float(vp)^2*ρ - 2μ₀
    λᵢ = float(vpᵢ)^2*ρᵢ - 2μᵢ
    Kᵢ = λᵢ + μᵢ*2/3
    # Plane-strain bulk modulus
    K̅₀ = λ₀ + μ₀ # TW below (36)
    # Young's modulus for matrix
    E0 = μ₀*(3λ₀ + 2μ₀)/(λ₀ + μ₀)
    # Poisson's ratio of the matrix
    ν₀ = λ₀/(2*(λ₀ + μ₀))

    # Some time saving terms
    t1 = α^2 - 1
    t2 = 1 - ν₀
    t3 = 1 - 2ν₀
    t4 = 3α^2
    t5 = 1 - α^2

    # D1, D2 and D3 from Tandon and Weng (1984) (just before equation (18)).
    D1 = 1 + 2*(μᵢ - μ₀)/(λᵢ - λ₀)
    D2 = (λ₀ + 2μ₀)/(λᵢ - λ₀)
    D3 = λ₀/(λᵢ - λ₀)

    # g and g' terms (appendix of Tandon and Weng 1984). g is for prolate spheroidal
    # inclusions (α>1), whilst g' is for disc-like (oblate) inclusions (α<1).
    if (α >= 1)
        acshα = log(α + sqrt(t1))
        g = α/sqrt(t1^3)*(α*sqrt(t1) - acshα)
        g = α/sqrt(t1^3)*(α*sqrt(t1) - acosh(α))
    else
        # g' below
        g = α/sqrt(t5^3)*(acos(α) - α*sqrt(t5))
    end

    # Eshelby's Sijkl tensor (appendix of Tandon and Weng 1984).
    s11 = (t3 + (t4 - 1)/t1 - (t3 + t4/t1)*g)/(2*t2)
    s22 = (t4/(t1*2) + (t3 - 9/(4*t1))*g)/(4*t2)
    s33 = s22
    s23 = (α^2/(2*t1) - (t3 + 3/(4*t1))*g)/(4*t2)
    s32 = s23
    s21 = (-2*α^2/t1 + (t4/t1 - t3)*g)/(4*t2)
    s31 = s21
    s12 = (-(t3 + 1/t1) + (t3 + 3/(2*t1))*g)/(2*t2)
    s13 = s12
    s44 = (α^2/(2*t1) + (t3 - 3/(4*t1))*g)/(4*t2)
    s66 = (t3 - (t1+2)/t1 - (t3 - 3*(t1 + 2)/t1)*g/2)/(4*t2)
    s55 = s66

    # Tandon and Weng's B terms (after equation 17).
    B1 = c*D1 + D2 + (1 - c)*(D1*s11 + 2*s21)
    B2 = c + D3 + (1 - c)*(D1*s12 + s22 + s23)
    B3 = c + D3 + (1 - c)*(s11 + (1+D1)*s21)
    B4 = c*D1 + D2 + (1 - c)*(s12 + D1*s22 + s23)
    B5 = c + D3 + (1 - c)*(s12 + s22 + D1*s23)

    # Tandon and Weng's A terms (after equation 20).
    A1 = D1*(B4 + B5) - 2*B2
    A2 = (1 + D1)*B2 - (B4 + B5)
    A3 = B1 - D1*B3
    A4 = (1 + D1)*B1 - 2*B3
    A5 = (1 - D1)/(B4 - B5)
    A = 2*B2*B3 - B1*(B4 + B5)

    # Tandon and Weng (1984) equations (25) (28) (31) (32)
    E11 = E0/(1 + c*(A1 + 2ν₀*A2)/A)
    E22 = E0/(1 + c*(-2*ν₀*A3 + (1 - ν₀)*A4 + (1 + ν₀)*A5*A)/(2*A))
    μ₁₂ = μ₀*(1 + c/(μ₀/(μᵢ - μ₀) + 2*(1 - c)*s66))
    μ₂₃ = μ₀*(1 + c/(μ₀/(μᵢ - μ₀) + 2*(1 - c)*s44))

    # Sayers equation (36)
    ν₃₁ = ν₀ - c*(ν₀*(A1 + 2*ν₀*A2)+(A3 - ν₀*A4))/(A + c*(A1 + 2*ν₀*A2))

    # T&W equation (36)
    num = (1 + ν₀)*(1 - 2ν₀)
    denom = 1 - ν₀*(1 + 2ν₃₁) + c*(2*(ν₃₁ - ν₀)*A3 + (1 - ν₀*(1 + 2ν₃₁))*A4)/A
    K₂₃ = K̅₀*num/denom
    ν₁₂² = E11/E22 - (1/μ₂₃ + 1/K₂₃)*E11/4

    # Cij - Sayers' (1992) equations (24)-(29).
    # Conversion
    C[1,1] = E11 + 4*ν₁₂²*K₂₃
    C[2,2] = μ₂₃ + K₂₃
    C[3,3] = C[2,2]
    C[1,2] = C[2,1] = 2*ν₃₁*K₂₃
    C[1,3] = C[3,1] = C[1,2]
    C[2,3] = C[3,2] = -μ₂₃ + K₂₃
    C[4,4] = (C[2,2] - C[2,3])/2
    C[5,5] = C[6,6] = μ₁₂

    # Apply density normalisation
    C = EC(C./ρ_out)

    C, ρ_out
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
    ϵ > 0.1 && @warn("Theory of Hudson (1982) only valid for `ϵ` < 0.1, but using ϵ = $ϵ")
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

"""
    thomsen_cracks(vp, vs, flu_fac, frac_den; delta_option=:weak_mensch) -> C

Create **density-normalised** elastic tensor for isotropic and porous rock
containing a single fracture set using Thomsen's (1995) formulation
assuming relatively low-frequency and low-porosity.  (Less restrictive equations
are also provided by Thomsen but not implemeted here yet.)

`vp` is the isotropic P-velocity of the porous host rock.

`vs` is the isotropic S-velocity of the porous host rock.

`flu_fac` is Thomsen's (1995) fluid factor parameter (Eq. 7a multiplied
by 1 - ratio of pore fluid-to-solid bulk moduli).
Value is confined to [0.0, 1.0]; weaker pore fluids correspond to larger
`flu_fac` (e.g., gas ≈ 0.0).

`frac_den` is the fracture density ``\\eta_c`` (Eq. 3), defined as
```math
    \\eta_c = N_v \\langle a^3 \\rangle / 8,
```
where ``N_v`` is the number density of fractures per unit volume and
``\\langle a^3\\rangle`` is the mean cubed fracture radius.

### Keyword arguments

- `delta_option = :weak_mensch`: Method to calculate C₁₃ from Thomsen's δ or
  δ* parameter, passed to [`thomsen_vti`](@ref).  See [`thomsen_vti`](@ref)
  for details.

### References

- Thomsen, L. (1995). Elastic anisotropy due to aligned cracks in porous rock.
  Geophysical Prospecting, 43(6), 805-829.
"""
function thomsen_cracks(vp, vs, flu_fac, frac_den; delta_option=:weak_mensch)
    if !(0 <= flu_fac <= 1)
        throw(ArgumentError("flu_fac must be in the range 0 to 1."))
    end

    # Define Thomsen's parameters
    vp_vs_ratio = vp/vs
    ν = 0.5/(1 - 1/vp_vs_ratio^2) - 1/(vp_vs_ratio^2 - 1)
    ϵ = (8/3)*flu_fac*frac_den # Eq. 6a; Always positive!
    γ = (8/3)*frac_den*(1 - ν)/(2 - ν) # Eq. 6b; Always positive! Note γ > ϵ for stiffer pore fluids and low aspect ratio
    δ = 2*(1 - ν)*ϵ - 2*γ*(1 - 2*ν)/(1 - ν) # Eq. 6c; Can reasonably be ±
    α = vp/sqrt(1 + 2*(((1 - ν)^2)/(1 - 2*ν))*ϵ)
    β = vs/sqrt(1 + 2*γ)

    thomsen_vti(vp, vs, ϵ, γ, δ; delta_option=delta_option)
end

"""
    thomsen_cracks(vp, vs, density, porosity, frac_aspect_ratio, frac_density, Kf; delta_option=:weak_mensch) -> C

Determine the density-normalised elastic constants using fracture density,
aspect ratio and fluid bulk modulus rather than the fluid factor.

`vp` is the isotropic P-velocity of the porous host rock.

`vs` is the isotropic S-velocity of the porous host rock.

`density` is the density of the porous host rock.

`porosity` is the background porosity of the host rock.

`frac_aspect_ratio` is the aspect ratio of the fractures (thickness/diameter).

`frac_density` is the fracture density (see Eq. 3 for definition).

`Kf` is the fluid bulk modulus (ρ_fluid*vp_fluid²)

!!! note
    The units of `vp`, `vs`, `density` and `Kf` must be consistent for `C` to
    be meaningful; e.g. velocities should be in m/s, density should be in kg/m³,
    and `Kf` should be in Pa.  No check is made on this.
"""
function thomsen_cracks(
    vp, vs, density, porosity, frac_aspect_ratio, frac_density, Kf;
    delta_option=:weak_mensch
)
    flu_fac = thomsen_fluid_factor(
        vp, vs, density, porosity, frac_aspect_ratio, frac_density, Kf
    )

    thomsen_cracks(vp, vs, flu_fac, frac_density; delta_option=delta_option)
end

"""
    thomsen_cracks(vp, vs, density, porosity, frac_aspect_ratio, frac_density, vp_fluid, density_fluid) -> C

Determine constants using fluid P-wave velocity `vp_fluid` and density
`density_fluid` rather than specifying bulk modulus directly.
"""
function thomsen_cracks(
    vp, vs, density, porosity, frac_aspect_ratio, frac_density, vp_fluid, density_fluid;
    delta_option=:weak_mensch
)
    Kf = density_fluid*vp_fluid^2
    thomsen_cracks(
        vp, vs, density, porosity, frac_aspect_ratio, frac_density, Kf;
        delta_option=delta_option
    )
end

"""
    thomsen_fluid_factor(vp, vs, porosity, density, frac_aspect_ratio, frac_density, Kf)

Compute fluid factor from Thomsen's (1995) theory of elasticity due to aligned
fractures in a porous medium.

`vp` is the isotropic P-velocity of the porous host rock.

`vs` is the isotropic S-velocity of the porous host rock.

`density` is the density of the porous host rock.

`porosity` is the background porosity of the host rock.

`frac_aspect_ratio` is the aspect ratio of the fractures (thickness/diameter).
Values less than one describe flat 'pancake' fractures, while ratios bigger
than 1 are for 'rugby ball'-shaped fractures.

`frac_density` is the fracture density ``\\eta_c``, defined by Thomsen (1995)
Eq. 3 as
```math
    \\eta_c = N_v \\langle a^3 \\rangle / 8,
```
where ``N_v`` is the number density of fractures per unit volume and
``\\langle a^3\\rangle`` is the mean cubed fracture radius.

`Kf` is the fluid bulk modulus (ρ_fluid*vp_fluid²).

!!! note
    The units of `vp`, `vs`, `density` and `Kf` must be consistent for `f` to
    be meaningful; e.g. velocities should be in m/s, density should be in kg/m³,
    and `Kf` should be in Pa.  No check is made on this.

### References

- Thomsen, L. (1995). Elastic anisotropy due to aligned cracks in porous rock.
  Geophysical Prospecting, 43(6), 805-829.
"""
function thomsen_fluid_factor(
    vp, vs, density, porosity, frac_aspect_ratio, frac_density, Kf
)
    # Derived elastic constants
    Gs = density*vs^2 # Solid shear modulus
    Ks = density*vp^2 - (4/3)*Gs # Solid bulk modulus
    νs = (3*Ks - 2*Gs)/(6*Ks + 2*Gs) # Solid poisson ratio
    frac_poro = (4/3)*pi*frac_aspect_ratio*frac_density # Fracture porosity

    # Fluid-influence factor terms
    Ac = (16/9)*(1 - νs^2)/(1 - 2*νs) # Eq. A13b
    Ap = (3/2)*(1 - νs)/(1 - 2*νs) # Eq. A30b
    Dcp = 1/(
        1 - Kf/Ks
        + (Kf/(Ks*(porosity + frac_poro)))*(Ap*porosity + Ac*frac_density)
    ) # Eq. 7a

    return (1 - Kf/Ks)*Dcp
end
