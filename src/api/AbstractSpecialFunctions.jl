# SPDX-License-Identifier: MIT
"""
This package provides abstract interfaces for Special Functions.

- Gamma Functions
- Exponential Integrals
- Error Functions
- Airy Functions
- Bessel Functions
- Struve Functions
- Parabolic Cylinders
- Hypergeometric Functions
- Legendre Functions
- 𝑞-Functions
- Orthogonal Polynomials
- Elliptic Integrals & Functions
- Bernoulli and Euler Polynomials
- Zeta Functions
- Combinatorial Analysis
- Number Theory Functions
- Mathieu Functions
- Spheroidal Waves
- Lamé & Heun Functions
- Coulomb Functions

# Reference

- [DLMF: NIST Digital Library of Mathematical Functions](https://dlmf.nist.gov/)
"""
module AbstractSpecialFunctions

include("gamma.jl")
include("exp_int.jl")

end # AbstractSpecialFunctions
