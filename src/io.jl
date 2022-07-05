# Input/output of tensors

"""
    CIJ.read(file) -> C, rho

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
    open(file) do io
        read(io; file=file)
    end
end

function read(io::IO; file=nothing)
    d = readdlm(io; comments=true)
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
    CIJ.write(file, C, rho, comment=nothing)
    CIJ.write(io::IO, C, rho, comment=nothing)

Write the density-normalised elastic constants `C` and density `rho` to the file
`file` or `IO` `io` in the following format:

    1 1 C11
    . . ...
    i j Cij
    . . ...
    6 6 C66
    7 7 rho

If supplied, `comment` should be a string of one or more lines to write at the
end of the output to describe the constants.  By default, the following is written:

    # Saved by <user> on <hostname> on <date> using <Julia version>

Such files usually have a `.ecs` file extension.
"""
function write(file, C::EC, rho, comment=nothing)
    isdir(dirname(file)) || error("directory \"$(dirname(file))\" does not exist")
    open(file, "w") do io
        write(io, C, rho, comment)
    end
end

function write(io::IO, C::EC, rho, comment=nothing)
    is_6x6(C) || throw(ArgumentError("elastic constants must be an EC or 6x6 matrix"))
    is_stable(C) && is_symm(C) ||
        @warn("elastic constants are not in the right form.  "
              * "May be asymmetric or unstable.")

    for i = 1:6, j = i:6
        @printf(io, "%d %d %12.6e\n", i, j, C[i,j]*rho)
    end
    @printf(io, "7 7 %9.3f\n", rho)

    if comment !== nothing && !isempty(comment)
        if !occursin(r"^#", comment)
            comment = "# " * comment
        end
        Base.write(io, replace(chomp(comment), "\n"=>"\n# ") * "\n")
    elseif comment === nothing
        user = ENV["USER"]
        hostname = readchomp(`hostname`)
        Base.write(io, "# Saved by user $user on $hostname on $(now()) using CIJ.jl on Julia $VERSION\n")
    end

    nothing
end

write(file_or_io, C::Matrix, rho, comment=nothing) = 
    write(file_or_io, EC(C), rho, comment)
