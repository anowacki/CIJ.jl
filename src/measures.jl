# Measures of anisotropy

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
