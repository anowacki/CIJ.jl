# Utility functions

"""
    symm(C) -> C

Return a copy of C with the upper half copied into the lower half to enforce symmetry.
"""
symm(C) = symm!(deepcopy(C))

"""
    symm!(C) -> C

Fill in the lower half of 6x6 matrix `C` with the upper half, making it symmetrical,
and returning `C`.

If `C` is an `EC`, then we assume it is symmetrical, since it is impossible to make
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
