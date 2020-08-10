# Some sample data

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
