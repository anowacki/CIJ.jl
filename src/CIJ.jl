# Copyright Andy Nowacki 2015-, all rights reserverd.
# See the file LICENSE.md for licence details.

"""
# CIJ

Deal with linear elastic constants, especially for geophysical applications.

## 
"""
module CIJ

import Dates.now
import LinearAlgebra
import LinearAlgebra: cross, dot, norm
import Printf.@printf

using StaticArrays
using Rotations


export
    # Types and constructors
    EC,
    # Constructors of elastic constants
    Panning_VTI,
    global_VTI,
    grechka_cracks,
    grechka_cracks!,
    iso,
    pitl,
    tandon_and_weng,
    thom,
    thom_st,
    # Properties of a tensor
    Au,
    C2S,
    C2S!,
    Hill_average,
    HillG,
    HillK,
    Reuss_average,
    ReussG,
    ReussK,
    VRH,
    Voigt_average,
    VoigtG,
    VoigtK,
    is_stable,
    phase_vels,
    # Conversion between representations
    S2C,
    S2C!,
    cij,
    cijkl,
    # Operations
    rot3,
    symm,
    symm!

"Lookup index matrix to convert between elasticity tensor and Voigt matrix"
const VOIGT_CONTRACTION_MATRIX = @SMatrix [1 6 5
                                           6 2 4
                                           5 4 3]

"Default precision of floating point numbers used in `EC` type"
const DEFAULT_FLOAT = Float64

#=
    Types
=#
"""
    EC{T} <: AbstractArray{T,2}

Wrapper type for elastic constants held as a 6×6 Voigt matrix.  The underlying storage
is a `StaticArrays.MArray`, meaning certain operations using this structure will be
quicker than using normal `Array`s.  For instance, elastic constant transformations
such as rotations, or conversion between stiffness and compliance, are each about
twice as fast.

    EC(data::AbstractArray{T}; warn=false) where T -> ec

Create a set of elastic constants `ec` from the array `data` which must have length 36.
Symmetry is enforced when converting, taking the upper half of the matrix.  No warning
is issued for non-symmetric input, unless `warn` is `true`.
"""
struct EC{T} <: AbstractArray{T,2}
    data::MArray{Tuple{6,6},T,2,36}

    function EC{T}(data; warn=false) where T
        ec = new{T}(MMatrix{6,6,T}(data))
        if warn
            atol = √eps(T)
            for i in 1:6, j in i+1:6
                if !isapprox(ec[i,j], ec[j,i], atol=atol)
                    @warn("input matrix not symmetrical: taking upper half")
                    break
                end
            end
        end
        for j in 1:6, i in j+1:6
            ec[i,j] = ec[j,i]
        end
        ec
    end
end

EC(x::AbstractArray{T}; kwargs...) where T = EC{float(T)}(MMatrix{6,6}(x); kwargs...)
EC{T}(x::NTuple{36}; kwargs...) where T = EC{T}(MArray{Tuple{6,6}}(x); kwargs...)
EC(x::NTuple{36,T}; kwargs...) where T = EC{DEFAULT_FLOAT}(MArray{Tuple{6,6},DEFAULT_FLOAT}(x); kwargs...)

Base.size(x::EC) = (6,6)
Base.length(x::EC) = 36
Base.getindex(x::EC, i...) = x.data[i...]
Base.setindex!(x::EC, val, i, j) = x.data[i,j] = x.data[j,i] = val
Base.IndexStyle(::Type{EC}) = IndexStyle(MMatrix)

Base.zero(::Type{EC{T}}) where T = EC{T}(zero(MMatrix{6,6,T}))
Base.zero(::Type{EC}) = zero(EC{DEFAULT_FLOAT})

#=
    Functions
=#
"""
    iso(; vp=nothing, vs=nothing, lam=nothing, mu=nothing, K=nothing, G=nothing) -> C

Return an isotropic Voigt stiffness matrix `C` from a pair of the following:

`vp`, `vs` : Isotropic velocities in m/s

`lam`, `mu`: Lamé parameters divided by density in m^2/s^2

`K`, `G`   : Bulk and shear moduli dvided by density in m^2/s^2
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
            error("iso: Must request a pair of vp,vs; lam,mu; K,G")
        end
        C[2,2] = C[3,3] = C[1,1]
        C[5,5] = C[6,6] = C[4,4]
        C[1,3] = C[2,3] = C[3,1] = C[3,2] = C[2,1] = C[1,2]
    end
    return C
end

"""
    cijkl(C) -> c

Convert the 6x6 Voigt matrix `C` into the 3x3x3x3 tensor `c`
"""
function cijkl(C)
    a = VOIGT_CONTRACTION_MATRIX
    if ! is_6x6(C)
        error("CIJ.cijkl: Input must be a 6x6 array")
    end
    c = MArray{Tuple{3,3,3,3},eltype(C),4,81}(undef)
    for i = 1:3, j = 1:3, k = 1:3, l = 1:3
        c[i,j,k,l] = C[a[i,j],a[k,l]]
    end
    return c
end

"""
    cij(c) -> C

Convert the 3x3x3x3 elasticity tensor `c` into the 6x6 Voigt matrix `C`
"""
function cij(c)
    a = VOIGT_CONTRACTION_MATRIX
    if size(c) != (3,3,3,3)
        error("CIJ.cij: Input must be a 3x3x3x3 tensor")
    end
    C = zero(EC{float(eltype(c))})
    for i = 1:3, j = i:3, k = 1:3, l = k:3
        C[a[i,j],a[k,l]] = c[i,j,k,l]
    end
    return C
end

"""
    thom(vp, vs, eps, gam, del) -> C

Return the 6x6 Voigt matrix defined `C` by the weak anisotropy parameters of
Thomsen (1986) Weak elastic anisotropy.  Geophysics, 51, 10, 1954-1966.

Output is density-normalised tensor.
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
    global_VTI(vp, vs, xi, phi, eta) -> C

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

function get_pol(inc, az, x, xs1)
    # Projection of fast shear wave polarisation onto wavefront plane
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
    for i = 1:3, j = 1:3, k = 1:3, l = 1:3
        m = ijkl[i,j]
        n = ijkl[k,l]
        T[i,k] = T[i,k] + C[m,n]*x[j]*x[l]
    end
    return LinearAlgebra.Hermitian(T) # Real-symmetric
end

"Return a 3-vector which is the cartesian direction corresponding to
 the azimuth `az` and inclination `inc` (degrees)"
function incaz2cart(inc, az)
    i = deg2rad(inc)
    a = deg2rad(az)
    cosi = cos(i)
    return @SVector [cos(a)*cosi, -sin(a)*cosi, sin(i)]
end

function VRH(VF1, C1, rh1, VF2, C2, rh2)
    # Return the VRH elasticity and density
    voigt = (C1.*VF1 + C2.*VF2)/(VF1 + VF2)
    reuss = (C2S(C1).*VF1 + C2S(C2).*VF2)/(VF1 + VF2)
    return ((voigt + S2C(reuss))/2, (rh1*VF1 + rh2*VF2)/(VF1 + VF2))
end

"""
    VRH(VF1, C1, rh1, VF2, C2, rh2) -> C, rho

Return the Voigt-Reuss-Hill-averaged 6x6 Voigt stiffness matrix `C` from combining
two sets of elastic constants with volume fractions `VF{1,2}`, constants `C{1,2}`
and density `rh{1,2}`.

    VRH(VF, C, rho) -> C, rho

Return the Voigt-Reus-Hill-averaged 6x6 Voigt stiffness matrix `C` and density `rho`
from combining n sets of elastic constants and densities, whose proportions are
listed in `VF`, where n is the number of constants and length of `VF` and `rho`.

`VF`: A vector (which need not sum to 1) containing the relative proportions of each
set of elastic constants, size n

`C`: A nx6x6 array of elastic constants

`rho`: A vector of densities of length n.

See: Hill, R., The elastic behaviour of a crystalline aggregate,
     P Phys Soc Lond A (1952) vol. 65 (389) pp. 349-355
"""
function VRH(VF, C::AbstractArray{T,3} where T, rho)
    length(VF) == size(C,1) == length(rho) ||
        throw(ArgumentError("arrays are not the same lengths"))
    size(C, 2) == size(C, 3) == 6 || throw(ArgumentError("C must have size (n,6,6)"))
    VRH(VF, [EC(C[i,:,:]) for i in 1:size(C, 1)], rho)
end

function VRH(VF, C::AbstractArray{<:EC}, rho)
    error("VRH hasn't been tested yet")
    length(VF) == length(C) == length(rho) || error("VF, C and rho must have same legnth")
    ΣVF = sum(VF)
    voigt = zero(EC)
    reuss = zero(EC)
    Cave = zero(EC)
    for k in 1:length(VF)
        voigt .+= VF[k]/ΣVF.*C[k]
        reuss .+= VF[k]/ΣVF.*C2S(C[k])
    end
    Cave = (voigt .+ S2C(reuss))./2
    rhave = sum(VF.*rho./ΣVF)
    Cave, rhave
end

"""
    VoigtK(C) -> K

Return the Voigt bound on a tensor `C`'s bulk modulus `K`.
"""
VoigtK(C) = (C[1,1] + C[2,2] + C[3,3] + 2*(C[1,2] + C[2,3] + C[1,3]))/9

"""
    VoigtG(C) -> G

Return the Voigt bound on a tensor `C`'s shear modulus `G`.
"""
VoigtG(C) = (C[1,1] + C[2,2] + C[3,3] - (C[1,2] + C[2,3] + C[1,3]) +
             3*(C[4,4] + C[5,5] + C[6,6]))/15

"""
    ReussK(C) -> K

Return the Reuss bound on a tensor `C`'s bulk modulus `K`.
"""
function ReussK(C)
    S = C2S(C)
    1/(S[1,1] + S[2,2] + S[3,3] + 2*(S[1,2] + S[2,3] + S[1,3]))
end

"""
    ReussG(C) -> G

Return the Reuss bound on a tensor `C`'s shear modulus `G`.
"""
function ReussG(C)
    S = C2S(C)
    15/(4*(S[1,1] + S[2,2] + S[3,3]) - 4*(S[1,2] + S[2,3] + S[1,3]) +
        3*(S[4,4] + S[5,5] + S[6,6]))
end

"""
    HillK(C) -> K

Return the Voigt-Reuss-Hill bound on a tensor `C`'s bulk modulus `K`.
"""
HillK(C) = (VoigtK(C) + ReussK(C))/2

"""
    HillG(C) -> G

Return the Voigt-Reuss-Hill bound on a tensor `C`'s shear modulus `G`.
"""
HillG(C) = (VoigtG(C) + ReussG(C))/2

"""
    Voigt_average(C) -> ⟨C⟩

Return the Voigt isotropic average of a tensor `C`.
"""
Voigt_average(C) = iso(K=VoigtK(C), G=VoigtG(C))

"""
    Reuss_average(C) -> ⟨C⟩

Return the Reuss isotropic average of a tensor `C`.
"""
Reuss_average(C) = iso(K=ReussK(C), G=ReussG(C))

"""
    Hill_average(C) -> ⟨C⟩

Return the Hill isotropic average of a tensor `C`.
"""
Hill_average(C) = iso(K=HillK(C), G=HillG(C))

"""
    Au(C) -> au

Return the Universal Elastic Anisotropy Index, `au`, of the tensor `C`.

See: Ranganathan & Ostoja-Starzewksi, Universal elastic anisotropy index,
     Phys Rev Lett (2008) vol. 101 (5) pp. 055504
"""
function Au(C)
    Kv, Gv, Kr, Gr = VoigtK(C), VoigtG(C), ReussK(C), ReussG(C)
    5*(Gv/Gr) + Kv/Kr - 6
end

"""
    rot3(C, a, b, x) -> Cr

Return a rotated tensor `Cr`, the result of rotating in turn the tensor `C`:

1. About the x1 axis by `a` degrees clockwise when looking down the x1 axis
2.   "    "  x2   "   " `b`    "        "       "     "      "   "  x2   "
3.   "    "  x3   "   " `b`    "        "       "     "      "   "  x3   "

Rotations are performed in that order,
"""
rot3(C, a, b, c) = transform(C, RotXYZ(deg2rad(a), deg2rad(b), deg2rad(c)))

"""
    transform(C, M) -> C′

Apply a transformation matrix `M` to a 6x6 matrix `C`, returning the transformed matrix `C′`.
"""
function transform(C, M)
    K = _transformation_matrix(M)
    return K * (C * transpose(K))
end

"""
    _transformation_matrix(M) -> K

Return the matrix `K` which represents the matrix with which to transform a Voigt
6x6 matrix by the 3x3 transformation matrix `M`.
"""
@inline function _transformation_matrix(M)
    K = MMatrix{6, 6, Float64}(undef)
    @inbounds for i = 1:3
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
            K[i,j+3] = 2*M[i,j1]*M[i,j2]
            K[i+3,j] = M[i1,j]*M[i2,j]
            K[i+3,j+3] = M[i1,j1]*M[i2,j2] + M[i1,j2]*M[i2,j1]
        end
    end
    K
end

"""
    C2S(C) -> S

Return the inverse of the stiffness matrix `C`, the compliance matrix `S`.
"""
C2S(C) = LinearAlgebra.inv(C)
C2S(C::EC) = LinearAlgebra.inv(C.data)

"""
    S2C(S) -> C

Return the inverse of the compliance matrix `S`, the stiffness matrix `C`.
"""
S2C(S) = C2S(S)

"""
    C2S!(C) -> S

Invert the stiffness matrix `C` in place, giving the compliance matrix `S`.
"""
C2S!(C) = C .= LinearAlgebra.inv!(LinearAlgebra.lu(C))
C2S!(C::EC) = C2S!(C.data)

"""
    S2C!(S) -> C

Invert the compliance matrix `S` in place, giving the siffness matrix `C`.
"""
S2C!(S) = C2S!(S)

"""
    is_iso(C) -> ::Bool

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
    is_stable(C) -> ::Bool

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
        LinearAlgebra.cholesky(C)
    catch PosDefException
        return false
    end
    return true
end
is_stable(C::EC) = is_stable(C.data)

"""
    ol() -> C, rho

Return the normalised elastic constants `C` and density `rho` for olivine,
handy for testing purposes
"""
ol() = EC{DEFAULT_FLOAT}((320.5,  68.1,  71.6,   0.0,   0.0,   0.0,
                           68.1, 196.5,  76.8,   0.0,   0.0,   0.0,
                           71.6,  76.8, 233.5,   0.0,   0.0,   0.0,
                            0.0,   0.0,   0.0,  64.0,   0.0,   0.0,
                            0.0,   0.0,   0.0,   0.0,  77.0,   0.0,
                            0.0,   0.0,   0.0,   0.0,   0.0,  78.7).*1.e9./3355.0), 3355.0


"""
    read(file) -> C, rho

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
    isfile(file) || error("file \"$file\" does not exist.")
    d = readdlm(file)
    size(d, 1) == 22 || error("file \"$file\" is not in expected format")
    C = zero(EC)
    rho = 0
    for l in 1:22
        i = Int(d[l,1])
        j = Int(d[l,2])
        if i == j == 7
            rho = d[l,3]
        else
            1 <= i <= j <= 6 || error("unexpected indices '$i' and '$j' in file \"$file\"")
            C[i,j] = d[l,3]
            C[j,i] = d[l,3]
        end
    end
    C, rho
end

"""
    write(C, rho, file, comment="")

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
function Base.write(C::Union{EC,Array{<:Number,2}}, rho::Number, file, comment::AbstractString="")
    is_6x6(C) || throw(ArgumentError("elastic constants must be an EC or 6x6 matrix"))
    is_stable(C) && is_symm(C) ||
        warn("elastic constants are not in the right form.  "
              * "May be asymmetric or unstable.")
    isdir(dirname(file)) || error("directory \"$(dirname(file))\" does not exist")
    open(file, "w") do f
        for i = 1:6, j=i:6
            @printf(f, "%d %d %12.6e\n", i, j, C[i,j]*rho)
        end
        @printf(f, "7 7 %9.3f\n", rho)
        if comment != ""
            if comment[1] != '#'
                comment = "# " * comment
            end
            write(f, replace(chomp(comment), "\n"=>"\n# ") * "\n")
        end
        user = ENV["USER"]
        hostname = readchomp(`hostname`)
        write(f, "# Saved by user $user on $hostname on $(now()) using CIJ.jl on Julia $VERSION\n")
    end
    nothing
end

"""
    symm(C) -> C

Return a copy of C with the upper half copied into the lower half to enforce symmetry.
"""
symm(C) = symm!(deepcopy(C))

"""
    symm!(C) -> C

Fill in the lower half of 6x6 matrix `C` with the upper half, making it symmetrical,
and returning `C`.

If `C` is an `EC`, then we assume it is synnetrical, since it is impossible to make
an asymmetrical `EC` without directly accessing its `.data` field, which is not
a supported way of manipulating `EC`s.
"""
function symm!(C)
    for i in 1:6, j in i+1:6
        C[j,i] = C[i,j]
    end
    C
end
# Because setindex! sets both upper and lower halves, assume symmetrical.
# May not be true if c.data has been accessed directly.
symm!(c::EC) = c

"""
    is_symm(C) -> ::Bool

Return `true` if `C` is symmetrical.
"""
function is_symm(C)
    for i = 1:6, j = i+1:6
        if !isapprox(C[i,j], C[j,i]) return false end
    end
    return true
end
is_symm(::EC) = true

"""
    is_6x6(C) -> ::Bool

Return `true` if C is a 6x6 matrix.
"""
is_6x6(C) = size(C) == (6,6)
is_6x6(::EC) = true

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

end # module CIJ
