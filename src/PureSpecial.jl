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
include("airy.jl")

end
