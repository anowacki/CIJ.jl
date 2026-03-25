# Create tensors with certain symmetry from medium parameters

"""
    iso(; vp=nothing, vs=nothing, lam=nothing, mu=nothing, K=nothing, G=nothing) -> C

Return an isotropic Voigt stiffness matrix `C` from a pair of the following:

`vp`, `vs` : Isotropic velocities in m/s

`lam`, `mu`: Lamé parameters **divided by density** in m^2/s^2

`K`, `G`   : Bulk and shear moduli **divided by density** in m^2/s^2
"""
function iso(; vp=nothing, vs=nothing, lam=nothing, mu=nothing, K=nothing, G=nothing)
    C = zero(EC)
    @inbounds begin
        if vp != nothing && vs != nothing
            C[1,1] = vp^2
            C[4,4] = vs^2
            C[1,2] = C[1,1] - 2C[4,4]
        elseif lam != nothing && mu != nothing
            C[1,1] = lam + 2mu
            C[4,4] = mu
            C[1,2] = lam
        elseif K != nothing && G != nothing
            C[1,1] = K + 4G/3
            C[4,4] = G
            C[1,2] = C[1,1] - 2C[4,4]
        else
            throw(ArgumentError("iso: Must request a pair of vp,vs; lam,mu; K,G"))
        end
        C[2,2] = C[3,3] = C[1,1]
        C[5,5] = C[6,6] = C[4,4]
        C[1,3] = C[2,3] = C[3,1] = C[3,2] = C[2,1] = C[1,2]
    end
    return C
end

"""
    vti(; vpv, vsv, vph, vsh, eta)

Return a density-normalised elastic tensor for a VTI medium defined using the velocities
Vpv, Vsv, Vph, Vsh and anisotropy parameter η, whereby

- Vpv is the P-wave velocity parallel to the symmetry axis
- Vsv is the S-wave velocity parallel to the symmetry axis
- Vph is the P-wave velocity perpendicular to the symmetry axis
- Vsh is the S-wave velocity for a horizontally-polarised shear wave perpendicular
  to the symmetry axis
- η describes how the velocities vary with angle from the symmetry axis
  and is C₁₃/(C₁₁ - 2C₄₄)

The constants are symmetric about the 3-axis

---
    vti(; A, C, L, N, F)

Return an unnormalised elastic tensor for a VTI medium defined using Love's (1927) notation,
whereby

- A is C₁₁, or ρ⋅Vph²
- C is C₃₃, or ρ⋅Vpv²
- L is C₄₄, or ρ⋅Vsv²
- N is C₆₆, or ρ⋅Vsh²
- F is C₁₃

### References

- Love AEH (1927) A Treatise on the Mathematical Theory of Elasticity. New York: Dover Publications.
"""
function vti(; vpv=nothing, vsv=nothing, vph=nothing, vsh=nothing, eta=nothing,
        A=nothing, C=nothing, L=nothing, N=nothing, F=nothing)
    if all(x->x!==nothing, (vpv, vsv, vph, vsh, eta))
        c = zero(EC)
        c[1,1] = c[2,2] = vph^2
        c[3,3] = vpv^2
        c[4,4] = c[5,5] = vsv^2
        c[6,6] = vsh^2
        c[1,3] = c[2,3] = c[3,1] = c[3,2] = eta*(c[1,1] - 2c[4,4])
        c[1,2] = c[2,1] = c[1,1] - 2c[6,6]
        c
    elseif all(x->x!==nothing, (A, C, L, N, F))
        EC{DEFAULT_FLOAT}((  A,  A-2N,  F,   0,   0,   0,
                           A-2N,   A,   F,   0,   0,   0,
                             F,    F,   C,   0,   0,   0,
                             0,    0,   0,   L,   0,   0,
                             0,    0,   0,   0,   L,   0,
                             0,    0,   0,   0,   0,   N))
    else
        throw(ArgumentError("Insufficient parameters defined."))
    end
end

"""
    thomsen_vti(α, β, ϵ, γ, δ; delta_option=:weak_mensch) -> C::EC

Returns density-normalised elastic constants from Thomsen's (1986)
parameters for hexagonal symmetry.  `C` has rotational symmetry about the
3-axis.

`α`, `β`, `ϵ`, `γ`, `δ` are Thomsen's (1986) parameters:

- `α`: P-wave velocity along the axis of symmetry, C₃₃
- `β`: S-wave velocity along the axis of symmetry, C₄₄
- `ϵ`: (C₁₁ - C₃₃)/2C₃₃
- `γ`: (C₆₆ - C₄₄)/2c₄₄
- `δ`: Determines how velocities change as a function of angle.  See `delta_option`
  below.

`delta_option` specifies which relationship to use in computing the
elastic constant C₁₃ from Thomsen's δ or δ* parameter. Options are:
- `:weak_thomsen`: Thomsen's (1986) weak anisotropy approximation (δ)
- `:weak_mensch`: Mensch and Rasolofosaon's (1997) formulation for weak anisotropy
  using a more accurate expansion
- `:exact`: Thomsen's (1986) exact anisotropy equations (δ*)

### References

- Thomsen, L. (1986). Weak elastic anisotropy. Geophysics, 51(10), 1954-1966.
- Mensch, T., & Rasolofosaon, P. (1997). Elastic-wave velocities in anisotropic media of arbitrary
   symmetry—generalization of Thomsen's parameters ε, δ and γ. Geophysical Journal International, 128(1), 43-64.
"""
function thomsen_vti(α, β, ϵ, γ, δ; delta_option=:weak_mensch)
    if α <= 0
        throw(ArgumentError("vp (α) must be greater than 0"))
    elseif β < 0
        throw(ArgumentError("vs (β) must be greater than or equal to 0"))
    end

    c33 = float(α)^2
    c44 = float(β)^2
    c11 = (1 + 2*ϵ)*c33
    c66 = (1 + 2*γ)*c44

    # Weak formula (Eq. 17) to be used with Thomsen's weak anisotropy equations (Eq. 16)
    if delta_option == :weak_thomsen
        sqrt_term = 2*δ*c33*(c33 - c44) + (c33 - c44)^2
        if sqrt_term < 0
            throw(ArgumentError("S velocity too high or delta too negative"))
        end
        c13 = sqrt(sqrt_term) - c44
    # Alternative (and generally more accurate) form from Mensch & Rasolofosaon
    # (1997), Eqs. 20, δₓ
    elseif delta_option == :weak_mensch
        c13 = δ*c33 + c33 - 2*c44
    # Exact formula (Eq. 8c) to be used with Thomsen's exact anisotropy equations (Eq. 10)
    elseif delta_option == :exact
        sqrt_term = δ*(c33^2) + (c33 - c44)*(c11 + c33 - 2*c44)/2
        if sqrt_term < 0
            throw(ArgumentError("S velocity too high or delta too negative"))
        end
        c13 = sqrt(sqrt_term) - c44
    else
        throw(ArgumentError("unknown delta option: ':$delta_option'"))
    end

    c12 = c11 - 2*c66

    EC((
        c11, c12, c13,   0,   0,   0,
        c12, c11, c13,   0,   0,   0,
        c13, c13, c33,   0,   0,   0,
          0,   0,   0, c44,   0,   0,
          0,   0,   0,   0, c44,   0,
          0,   0,   0,   0,   0, c66
    ))
end

"""
    thom(vp, vs, eps, gam, del) -> C

Return the 6x6 Voigt matrix `C` defined by the weak anisotropy parameters of
Thomsen (1986).

Output is density-normalised tensor.

!!! warning
    `thom` is deprecated in favour of [`thomsen_vti`](@ref) and may be removed
    in a future breaking release.

### References

Thomsen, L. (1986).  Weak elastic anisotropy.  Geophysics, 51, 10, 1954-1966.
"""
thom(vp, vs, eps, gam, del) = thomsen_vti(vp, vs, eps, gam, del; delta_option=:weak_thomsen)

"""
    thom_st(vp, vs, eps, gam, delst) -> C

Return the 6x6 Voigt matrix `C` defined by the general anisotropy parameters of
Thomsen (1986), where `delst` is Thomsen's δ*.

Output has same units as input.

!!! warning
    `thom_st` is deprecated in favour of [`thomsen_vti`](@ref) and may be removed
    in a future breaking release.

### References

Thomsen, L. (1986).  Weak elastic anisotropy.  Geophysics, 51, 10, 1954-1966.
"""
thom_st(vp, vs, eps, gam, delst) = thomsen_vti(vp, vs, eps, gam, delst; delta_option=:exact)

"""
    global_VTI(vp, vs, xi, phi, eta) -> C

Return the 6x6 Voigt matrix `C` defined by the radial anisotropy parameters
as used typically in global seismology.
"""
function global_VTI(vp, vs, xi, phi, eta)
    if vp <= 0
        error("CIJ.global_VTI: vp must be greater than 0")
    elseif vs < 0
        error("CIJ.global_VTI: vs must be greater than or equal to 0")
    elseif xi <= 0 || phi <= 0 || eta <= 0
        error("CIJ.global_VTI: xi, phi and eta must be greater then 0")
    end

    L = 15*((3*phi + 8 + 4*eta)*vs^2 - (phi + 1 - 2*eta)*vp^2) /
        ((6 + 4*eta + 5*xi)*(3*phi + 8 +4*eta) - 8*(phi + 1 -2*eta)*(1 - eta))
    A = (15*vp^2 - 8*(1 - eta)*L) / (3*phi + 8 + 4*eta)
    F = eta*(A - 2L)
    C = phi*A
    N = xi*L
    C12 = A - 2N
    return EC{DEFAULT_FLOAT}(( A,  C12, F,  0,  0,  0,
                              C12,  A,  F,  0,  0,  0,
                               F,   F,  C,  0,  0,  0,
                               0,   0,  0,  L,  0,  0,
                               0,   0,  0,  0,  L,  0,
                               0,   0,  0,  0,  0,  N))
end

"""
    Panning_VTI(vp, vs, xi, phi, eta) -> C

Return the 6x6 Voigt matrix `C` defined by the radial anisotropy parameters
defined by Panning and Romanowicz as an approximation to general VTI.
"""
function Panning_VTI(vp, vs, xi, phi)
    L = 3*vs^2/(2 + xi)
    N = xi*L
    A = 5*vp^2/(4 + phi)
    C = phi*A
    F = A - 2*L
    C12 = A - 2*N
    return EC{DEFAULT_FLOAT}(( A,  C12, F,  0,  0,  0,
                              C12,  A,  F,  0,  0,  0,
                               F,   F,  C,  0,  0,  0,
                               0,   0,  0,  L,  0,  0,
                               0,   0,  0,  0,  L,  0,
                               0,   0,  0,  0,  0,  N))
end

"""
    lame(vp, vs, ρ) -> λ, μ, κ

Return the Lamé parameters `λ`, `μ` and `κ` from the velocities `vp` and `vs`
and density `ρ`.
"""
function lame(vp, vs, ρ)
    μ = ρ*vs^2
    λ = ρ*(vp^2 - 2vs^2)
    κ = ρ*vp^2 - 4/3*μ
    λ, μ, κ
end
