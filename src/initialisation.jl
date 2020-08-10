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
Vpv, Vsv, Vph, Vsh and anisotropy parameter η.

    vti(; A, C, L, N, F)

Return an unnormalised elastic tensor for a VTI medium defined using Love's (1927) notation.

The constants are symmetric about the 3-axis.

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
    thom(vp, vs, eps, gam, del) -> C

Return the 6x6 Voigt matrix defined `C` by the weak anisotropy parameters of
Thomsen (1986).

Output is density-normalised tensor.

### References

Thomsen, L. (1986).  Weak elastic anisotropy.  Geophysics, 51, 10, 1954-1966.
"""
function thom(vp, vs, eps, gam, del)
    if vp <= 0
        error("CIJ.thom: vp must be greater than 0")
    elseif vs < 0
        error("CIJ.thom: vs must be greater than or equal to 0")
    end

    C = zero(EC)
    @inbounds begin
        C[3,3] = vp^2
        C[1,1] = C[2,2] = C[3,3]*(2*eps + 1)
        C[4,4] = C[5,5] = vs^2
        C[6,6] = C[4,4]*(2*gam + 1)

        b = 2C[4,4]
        term = C[3,3] - C[4,4]
        c = C[4,4]^2 - (2*del*C[3,3]*term + term^2)
        d = b^2 - 4*c
        if (d < 0)
            error("CIJ.thom: S-velocity too high or delta too negative")
        end
        C[1,3] = C[2,3] = C[3,1] = C[3,2] = -b/2 + sqrt(d)/2
        C[1,2] = C[2,1] = C[1,1] - 2*C[6,6]
    end
    return C
end

"""
    thom_st(vp, vs, eps, gam, delst) -> C

Return the 6x6 Voigt matrix `C` defined by the general anisotropy parameters of
Thomsen (1986).

Output has same units as input.
"""
function thom_st(vp, vs, eps, gam, delst)
    if vp <= 0
        error("CIJ.thom: vp must be greater than 0")
    elseif vs < 0
        error("CIJ.thom: vs must be greater than or equal to 0")
    end

    C = zero(EC)
    @inbounds begin
        C[3,3] = vp^2
        C[1,1] = C[2,2] = C[3,3]*(2*eps + 1)
        C[4,4] = C[5,5] = vs^2
        C[6,6] = C[4,4]*(2*gam + 1)
        a = 2.0
        b = 4*C[4,4]
        c = C[4,4]^2 - 2*delst*C[3,3]^2 - (C[3,3] - C[4,4])*(C[1,1] + C[3,3] - 2*C[4,4])
        if b^2 - 4*a*c < 0
            error("CIJ.thom_st: S velocity too high or delta too negative")
        end
        C[1,3] = C[3,1] = C[2,3] = C[3,2] = (-b + sqrt(b^2 - 4*a*c))/(2*a)
        C[1,2] = C[2,1] = C[1,1] - 2C[6,6]
    end
    return C
end

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
