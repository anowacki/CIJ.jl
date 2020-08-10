# Functions related to the averaging of several tensors

"""
    VRH(VF₁, C₁, ρ₁, VF₂, C₂, ρ₂) -> ⟨C⟩, ⟨ρ⟩

Return the Voigt-Reuss-Hill-averaged 6x6 Voigt stiffness matrix `C` from combining
two sets of elastic constants with volume fractions `VF₁,₂`, constants `C₁,₂`
and density `ρ₁,₂`.
"""
function VRH(VF1, C1, rh1, VF2, C2, rh2)
    # Return the VRH elasticity and density
    voigt = (C1.*VF1 + C2.*VF2)/(VF1 + VF2)
    reuss = (C2S(C1).*VF1 + C2S(C2).*VF2)/(VF1 + VF2)
    return ((voigt + S2C(reuss))./2, (rh1*VF1 + rh2*VF2)/(VF1 + VF2))
end

"""
    VRH(VF, C, ρ) -> ⟨C⟩, ⟨ρ⟩

Return the Voigt-Reus-Hill-averaged 6x6 Voigt stiffness matrix `⟨C⟩` and density `⟨ρ⟩`
from combining n sets of elastic constants and densities, whose proportions are
listed in `VF`, where n is the number of constants and length of `VF` and `ρ`.

    VRH(C) -> ⟨C⟩

Compute the average when all tensors have the same density and equal volume fraction.

### Arguments

`VF`: A vector (which need not sum to 1) containing the relative proportions of each
set of elastic constants, size n

`C`: A n×6×6 array of elastic constants, or an n-length vector of `EC`s.

`ρ`: A vector of densities of length n.

### References

Hill, R., The elastic behaviour of a crystalline aggregate,
     P Phys Soc Lond A (1952) vol. 65 (389) pp. 349-355
"""
function VRH(VF, C::AbstractArray{T,3} where T, rho)
    length(VF) == size(C,1) == length(rho) ||
        throw(ArgumentError("arrays are not the same lengths"))
    size(C, 2) == size(C, 3) == 6 || throw(ArgumentError("C must have size (n,6,6)"))
    VRH(VF, [EC(C[i,:,:]) for i in 1:size(C, 1)], rho)
end

function VRH(VF, C::AbstractArray{<:EC{T}}, rho) where T
    # error("VRH hasn't been tested yet")
    length(VF) == length(C) == length(rho) || error("VF, C and rho must have same length")
    ΣVF = sum(VF)
    voigt = zero(EC{T})
    reuss = zero(EC{T})
    Cave = zero(EC{T})
    for k in 1:length(VF)
        voigt .+= VF[k]/ΣVF.*C[k]
        reuss .+= VF[k]/ΣVF.*C2S(C[k])
    end
    Cave = (voigt .+ S2C(reuss))./2
    rhave = sum(VF.*rho./ΣVF)
    Cave, rhave
end

function VRH(C::AbstractArray{<:EC{T}}) where T
    n = length(C)
    Cave = zero(EC{T})
    voigt = similar(Cave)
    reuss = similar(Cave)
    S = [C2S(CC) for CC in C]
    @inbounds for k in 1:n
        for j in 1:6, i in j:6
            voigt[i,j] += C[k][i,j]
            reuss[i,j] += S[k][i,j]
        end
    end
    Cave./n
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

Return the Voigt isotropic average of a single tensor `C`.
"""
Voigt_average(C) = iso(K=VoigtK(C), G=VoigtG(C))

"""
    Reuss_average(C) -> ⟨C⟩

Return the Reuss isotropic average of a single tensor `C`.
"""
Reuss_average(C) = iso(K=ReussK(C), G=ReussG(C))

"""
    Hill_average(C) -> ⟨C⟩

Return the Hill isotropic average of a single tensor `C`.
"""
Hill_average(C) = iso(K=HillK(C), G=HillG(C))
