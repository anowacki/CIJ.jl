# Copyright Andy Nowacki 2015-, all rights reserved.
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
    is_iso,
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


# Types
include("types.jl")

# Functions
include("initialisation.jl")
include("inquiry.jl")
include("conversion.jl")
include("transformation.jl")
include("util.jl")
include("averaging.jl")
include("effective_medium.jl")
include("measures.jl")
include("velocities.jl")
include("io.jl")

# Data
include("sample_data.jl")

end # module CIJ
