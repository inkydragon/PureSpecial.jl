# SPDX-License-Identifier: MIT
"""
# Special functions in Pure Julia

- Elementary Functions
- Gamma Functions
- Exponential and Trigonometric Integrals
- Error Functions
- Airy Functions
- Bessel Functions
- Struve Functions
- Parabolic Cylinder Functions
- Hypergeometric Functions
- Legendre Functions
- Elliptic Integrals
- Elliptic Functions
- Zeta Functions
- Mathieu Functions
- Spheroidal Wave Functions
- Miscellaneous Functions
"""
module PureSpecial

# Submod
include("specfun/Specfun.jl")
import .Specfun

# APIs
include("elementary.jl")
include("gamma.jl")
include("exp_int.jl")
include("error.jl")
include("airy.jl")
include("bessel.jl")
include("struve.jl")
include("parabolic.jl")
include("hyper.jl")
include("legendre.jl")
include("elliptic.jl")
include("zeta.jl")
include("mathieu.jl")
include("spheroidal.jl")

end
