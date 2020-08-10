# Input/output of tensors

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
