"""
    _cij2x(c) -> x::SVector{21}

Convert the Voigt matrix-form elastic tensor `c` into the
Browaeys & Chevrot's (2004) elastic vector 𝐗, which has 21
components.
"""
function _cij2x(c)
    @SVector[
        c[1,1], c[2,2], c[3,3], √2*c[2,3], √2*c[1,3], √2*c[1,2],
        2c[4,4], 2c[5,5], 2c[6,6], 2c[1,4], 2c[2,5], 2c[3,6],
        2c[3,4], 2c[1,5], 2c[2,6], 2c[2,4], 2c[3,5], 2c[1,6],
        2√2*c[5,6], 2√2*c[4,6], 2√2*c[4,5]
    ]
end

"""
    _x2cij(x) -> C::SMatrix{6,6}

Convert the Browaeys & Chevrot (2004) elastic vector 𝐗, which has 21
components, into the Voigt matrix-form elastic tensor `C`.
"""
function _x2cij(x)
    @SMatrix[
        x[1]     x[6]/√2  x[5]/√2  x[10]/2   x[14]/2   x[18]/2
        x[6]/√2  x[2]     x[4]/√2  x[16]/2   x[11]/2   x[15]/2
        x[5]/√2  x[4]/√2  x[3]     x[13]/2   x[17]/2   x[12]/2
        x[10]/2  x[16]/2  x[13]/2  x[7]/2    x[21]/2√2 x[20]/2√2
        x[14]/2  x[11]/2  x[17]/2  x[21]/2√2 x[8]/2    x[19]/2√2
        x[18]/2  x[15]/2  x[12]/2  x[20]/2√2 x[19]/2√2 x[9]/2
    ]
end

"""
    decompose(C) -> (Ciso, p_iso, Chex, p_hex, Ctet, p_tet, Corth, p_orth, Cmono, p_mono, Ctri, p_tri)

Decompose the elastic tensor `C` (a 6⨯6 Voigt matrix) into isotropic,
hexagonal, tetragonal, orthorhombic, monoclinic and triclinic parts
according to the method of Browaeys & Chevrot (2004).

Return a named tuple with the following properties:
- `Ciso`: Isotropic component
- `p_iso`: Proportion of the whole tensor which is isotropic
- `Chex`: Hexaganal component
- `p_hex`: Proportion of the whole tensor which is hexagonal
- `Ctet`: Tetragonal component
- `p_tet`: Proportion of the whole tensor which is tetragonal
- `Corth`: Orthorhombic component
- `p_orth`: Proportion of the whole tensor which is orthorhombic
- `Cmono`: Monoclinic component
- `p_mono`: Proportion of the whole tensor which is monoclinic
- `Ctri`: Triclinic component
- `p_tri`: Proportion of the whole tensor which is triclinic

# Example
```
julia> d = decompose(CIJ.ol()[1]);

julia> d.Ciso
6×6 EC{Float64}:
 7.08058e7  2.339e7    2.339e7    0.0        0.0        0.0
 2.339e7    7.08058e7  2.339e7    0.0        0.0        0.0
 2.339e7    2.339e7    7.08058e7  0.0        0.0        0.0
 0.0        0.0        0.0        2.37079e7  0.0        0.0
 0.0        0.0        0.0        0.0        2.37079e7  0.0
 0.0        0.0        0.0        0.0        0.0        2.37079e7

julia> d.p_iso
0.8159622595948542
```

# Reference
- Browaeys, J., Chevrot, S., 2004.
  Decomposition of the elastic tensor and geophysical applications.
  *Geophysical Journal International* **159**, 667–678.
  https://doi.org/10.1111/j.1365-246X.2004.02415.x

"""
function decompose(C)
    T = eltype(C)

    # Input elastic vector
    X = _cij2x(C)

    # Projection matrix
    M = zero(MMatrix{21,21,T,21*21})

    # Isotropic part
    _Miso!(M)
    Xiso = M*X
    Ciso = _x2cij(Xiso)
    X_minus_iso = _cij2x(C - Ciso)

    # Hexagonal part
    _Mhex!(M)
    Xhex = M*X_minus_iso
    Chex = _x2cij(Xhex)
    X_minus_hex = _cij2x(C - Ciso - Chex)

    # Tetragonal part
    _Mtet!(M)
    Xtet = M*X_minus_hex
    Ctet = _x2cij(Xtet)
    X_minus_tet = _cij2x(C - Ciso - Chex - Ctet)

    # Orthorhombic part
    _Morth!(M)
    Xorth = M*X_minus_tet
    Corth = _x2cij(Xorth)
    X_minus_orth = _cij2x(C - Ciso - Chex - Ctet - Corth)

    # Monoclinic part
    _Mmono!(M)
    Xmono = M*X_minus_orth
    Cmono = _x2cij(Xmono)
    X_minus_mono = _cij2x(C - Ciso - Chex - Ctet - Corth - Cmono)

    # Triclinic part
    Ctri = C - Ciso - Chex - Ctet - Corth - Cmono

    # Proportions
    ΣX² = sum(abs2, X)
    p_iso = 1 - sqrt(sum(abs2, X_minus_iso)/ΣX²)
    p_hex = 1 - sqrt(sum(abs2, X_minus_hex)/ΣX²) - p_iso
    p_tet = 1 - sqrt(sum(abs2, X_minus_tet)/ΣX²) - p_iso - p_hex
    p_orth = 1 - sqrt(sum(abs2, X_minus_orth)/ΣX²) - p_iso - p_hex - p_tet
    p_mono = 1 - sqrt(sum(abs2, X_minus_mono)/ΣX²) - p_iso - p_hex - p_tet - p_orth
    p_tri = 1 - p_iso - p_hex - p_tet - p_orth - p_mono

    (; Ciso, p_iso, Chex, p_hex, Ctet, p_tet, Corth, p_orth, Cmono, p_mono, Ctri, p_tri)
end

function _Miso!(M::AbstractMatrix{T}) where T
    M .= 0
    M[1:3,1:3] .= T(3)/15
    M[1:3,4:6] .= √(T(2))/15
    M[4:6,1:3] .= √(T(2))/15
    M[1:3,7:9] .= T(2)/15
    M[7:9,1:3] .= T(2)/15
    M[4:6,4:6] .= T(4)/15
    M[4:6,7:9] .= -√(T(2))/15
    M[7:9,4:6] .= -√(T(2))/15
    M[7:9,7:9] .= T(1)/5
    M    
end

function _Mhex!(M::AbstractMatrix{T}) where T
    M[1:9,1:9] .= @SMatrix[
             T(3)/8          T(3)/8    0      0      0    T(1)/4√(T(2))      0      0           T(1)/4
             T(3)/8          T(3)/8    0      0      0    T(1)/4√(T(2))      0      0           T(1)/4
                  0               0    1      0      0                0      0      0                0
                  0               0    0 T(1)/2 T(1)/2                0      0      0                0
                  0               0    0 T(1)/2 T(1)/2                0      0      0                0
    T(1)/(4√(T(2))) T(1)/(4√(T(2)))    0      0      0           T(3)/4      0      0 -T(1)/(2√(T(2)))
                  0               0    0      0      0                0 T(1)/2 T(1)/2                0
                  0               0    0      0      0                0 T(1)/2 T(1)/2                0
             T(1)/4          T(1)/4    0      0      0 -T(1)/(2√(T(2)))      0      0           T(1)/2
    ]
    M
end

function _Mtet!(M::AbstractMatrix{T}) where T
    M[1:9,1:9] .= @SMatrix[
        T(1)/2 T(1)/2      0      0      0      0      0      0      0
        T(1)/2 T(1)/2      0      0      0      0      0      0      0
             0      0      1      0      0      0      0      0      0
             0      0      0 T(1)/2 T(1)/2      0      0      0      0
             0      0      0 T(1)/2 T(1)/2      0      0      0      0
             0      0      0      0      0      1      0      0      0
             0      0      0      0      0      0 T(1)/2 T(1)/2      0
             0      0      0      0      0      0 T(1)/2 T(1)/2      0
             0      0      0      0      0      0      0      0      1
    ]
    M
end

function _Morth!(M)
    M .= 0
    for i in 1:10
        M[i,i] = 1
    end
    M
end

function _Mmono!(M)
    M .= 0
    for i in 1:21
        M[i,i] = 1
    end
    M
end
