"""
    symmetry(c; tol=1e-5) -> (; C, R, symmetry)

Determine the crystal symmetry of the 6⨯6 Voigt elastic stiffness matrix `c`
and return `C′` (`c` rotated into the optimum orientation) `R` (the
rotation matix which transforms `c` in to `C`) and `symmetry` (a `Symbol`
giving the symmetry type).  The method is that of Browaeys & Chevrot (2004)

`symmetry` takes on of the following values and axes are aligned as follows:

- `:isotropic`: Orientation is arbitrary so `C` is the same as `c`
- `:trigonal`: The symmetry axis is the new x₃
- `:hexagonal`: The symmetry axis is the new x₃
- `:tetragonal`: The symmetry axis is the new x₃
- `:orthorhombic`: x₁ is the fast P-wave direction, x₃ is the slow P-wave
  direction and x₂ is the intermediate axis
- `:monoclinic`: New axes are somewhat arbitrary
- `:triclinic`: New axes are somewhat arbitrary

!!! note
    `symmetry` is current experimental and the exact orientation of the
    transformed constants `C` are subject to change for the lower symmetry
    classes (monoclinic and triclinic).

# Reference
- Browaeys, J., Chevrot, S., 2004.
  Decomposition of the elastic tensor and geophysical applications.
  *Geophysical Journal International* **159**, 667–678.
  https://doi.org/10.1111/j.1365-246X.2004.02415.x
"""
function symmetry(C; tol=1e-5)
    is_6x6(C) || throw(ArgumentError("matrix must be 6×6"))
    is_symm(C) || @warn("matrix not symmetrical: using upper half only")

    T = eltype(C)

    # Identity transformation matrix
    I = @SMatrix[
        T(1) 0 0
        0 T(1) 0
        0 0 T(1)
    ]

    # Scaled version to avoid numerical errors in eigendecomposition
    scale = maximum(abs, C)
    C′ = SMatrix{6,6,T}(C./scale)

    # Tolerance on symmetries based on norm of tensor
    C_tol = sqrt(sum(abs2, C′))*tol

    # Dilatational stiffness tensor
    d = _dilatational_stiffness(C′)

    # Voigt stiffness tensor
    v = _voigt_stiffness(C′)

    # Eigenvectors and -values, sorted by increasing eigenvalue
    dvals, dvecs = eigen(d)
    vvals, vvecs = eigen(v)

    # Get separate eigenvectors and -values
    dvec1 = dvecs[:,1]
    dvec2 = dvecs[:,2]
    dvec3 = dvecs[:,3]
    d1, d2, d3 = dvals

    # Find number of distinct directions based on distinct eigenvalues
    # of d
    vtol = sqrt(sum(abs2, dvecs))*tol
    num_dirs = if all(≈(d1; atol=vtol), (d2, d3))
        1
    elseif all(x -> abs(x) > vtol, (d1 - d2, d2 - d3, d3 - d1))
        3
    else
        2
    end

    # Isotropic case: no particular directions and no need to rotate
    if num_dirs == 1
        symmetry = :isotropic
        Cout = C
        R = I

    # Hexagonal or tetragonal; the distinct axis is the hexad or triad
    elseif num_dirs == 2
        # Find symmetry axis `x3`
        x3 = if abs(d1 - d2) <= vtol
            # Old x3 is unique
            dvec3
        elseif abs(d2 - d3) <= vtol
            # Old x1 is unique
            dvec1
        else
            # Old x2 is unique
            dvec2
        end

        # If the new x3 axis is the same as the old, no need to rotate
        if ≈(x3[1], 0; atol=vtol) && ≈(x3[2], 0; atol=vtol)
            Cout = C
            R = I
        # Otherwise rotate with arbitrary axes x1 and x2
        else
            # Arbitrary unit vector in the x1-x2 plane
            x2 = normalize(x3 × @SVector[x3[1], x3[2], 0])
            # Other direction also in the plane and mutually orthogonal
            x1 = normalize(x3 × x2)
            R = hcat(x1, x2, x3)'
            Cout = transform(C, R)
        end

        # Decide between hexagonal, trigonal and tetragonal;
        symmetry = if abs(Cout[1,4]) > vtol
            :trigonal
        elseif ≈(Cout[6,6], (Cout[1,1] - Cout[1,2])/2; rtol=tol)
            :hexagonal
        else
            :tetragonal
        end

    # Orthorhombic or lower
    else
        # Find number of coincident eigenvectors between the dilatational
        # and Voigt stress tensors
        num_coincident = _num_coincident_columns(dvecs, vvecs, vtol)

        # Three coincident eigenvectors ⇒ orthorhombic
        # x1 is aligned with the smallest dilatational stiffness
        # eigenvector, and x3 is the largest.  This normally means
        # Vp along x1 is slowest and Vp along x3 is fastest.
        if num_coincident == 3
            imin = argmin(dvals)
            imin == 1 || error("Smallest eigenvalue is not the first")
            imax = argmax(dvals)
            imax == 3 || error("Largest eigenvalue is not the last")

            R = dvecs'
            Cout = transform(C, R')
            symmetry = :orthorhombic

        # One coincident eigenvector => monoclinic
        # Two coincident eigenvectors => triclinic
        # For these, take the bisectrices between the closest pairs
        # of the eigenvectors of d and v
        else
            # Mutable rotation matrix
            R_m = zero(MMatrix{3,3,T,9})
            for i in 1:3
                # Index of nearest v eigenvector to this d eigenvector
                # (pair with largest absolute dot product)
                _, j = findmax(ii -> abs(dot(dvecs[:,i], v[:,ii])), 1:3)
                # Bisectrix
                x_bisec = if dot(dvecs[:,i], vvecs[:,j]) < 0
                    normalize(dvecs[:,i] .- vvecs[:,j])
                else
                    normalize(dvecs[:,i] .+ vvecs[:,j])
                end
                R_m[i,:] .= x_bisec
            end

            R = SMatrix(R_m)'
            Cout = transform(C, R)

            symmetry = if num_coincident == 1
                :monoclinic
            else
                :triclinic
            end
        end
    end

    (; C=Cout, R, symmetry)
end

"""
    _dilatational_stiffness(C) -> d

Return the dilatational stiffness tensor for a Voigt stiffness matrix `C`.
"""
function _dilatational_stiffness(C)
    @SMatrix[
         C[1,1]+C[1,2]+C[1,3] C[1,6]+C[2,6]+C[3,6] C[1,5]+C[2,5]+C[3,5]
         C[1,6]+C[2,6]+C[3,6] C[1,2]+C[2,2]+C[3,2] C[1,4]+C[2,4]+C[3,4]
         C[1,5]+C[2,5]+C[3,5] C[1,4]+C[2,4]+C[3,4] C[1,3]+C[2,3]+C[3,3]
    ]
end

"""
    _voigt_stiffness(C) -> v

Return the Voigt stiffness tensor for a Voigt stiffness matrix `C`.
"""
function _voigt_stiffness(C)
    @SMatrix[
         C[1,1]+C[6,6]+C[5,5] C[1,6]+C[2,6]+C[4,5] C[1,5]+C[3,5]+C[4,6]
         C[1,6]+C[2,6]+C[4,5] C[6,6]+C[2,2]+C[4,4] C[2,4]+C[3,4]+C[5,6]
         C[1,5]+C[3,5]+C[4,6] C[2,4]+C[3,4]+C[5,6] C[5,5]+C[4,4]+C[3,3]
    ]
end

"""
    _num_coincident_columns(a, b, tol) -> count

Return the count of how many of the columns of 3×3 matrices `a` and `b`
are distinct.  Distinct vectors have an L1 norm of their
difference less than or equal to `tol`.
"""
function _num_coincident_columns(a, b, tol)
    sum(x -> norm(x) <= tol, aa - bb for aa in eachcol(a) for bb in eachcol(b))
end
