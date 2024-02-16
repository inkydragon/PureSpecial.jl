# SPDX-License-Identifier: MIT
"""

## Special functions
> ref: [scipy.special - SciPy Manual](https://docs.scipy.org/doc/scipy/reference/special.html)

**TOC**
- Airy functions
- Elliptic functions
- Bessel functions
    + Zeros of Bessel functions
    + Faster Common Bessel functions
    + Integrals of Bessel functions
    + Derivatives of Bessel functions
    + Spherical Bessel functions
    + Riccati-Bessel functions
- Struve functions
- Statistical functions: Skip!
- Gamma functions
- Error function
- Fresnel integrals
- Legendre functions
- Ellipsoidal harmonics
- Orthogonal polynomials
    + evaluate values
    + roots
- Hypergeometric functions
- Parabolic Cylinder functions
- Mathieu functions
- Spheroidal Wave functions
- Kelvin functions
- Combinatorics
- Lambert W functions
- Other Special functions
"""

module Scipy4j

# Submod
include("specfun/Specfun.jl")
import .Specfun

# APIs
include("airy.jl")

end
