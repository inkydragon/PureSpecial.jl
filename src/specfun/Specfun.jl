# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Special functions in pure julia.

Translated from `specfun.f`

```fortran
C
C       COMPUTATION OF SPECIAL FUNCTIONS
C
C          Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs.
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
C
```
"""
module Specfun

include("const.jl")

## 1~4
include("bernoulli_euler.jl")
# 2. Orthogonal Polynomials
include("gamma.jl")
include("legendre.jl")

## 5~7
# 5 Bessel Functions
include("bessel_zeros.jl")
# 6 Modified Bessel Functions
# 7 Integrals of Bessel Functions

## 8~12
include("kelvin.jl")
include("airy.jl")
include("struve.jl")
include("hyper.jl")

## 13~16
include("parabolic.jl")
include("mathieu.jl")
include("spheroidal.jl")
include("error.jl")

## 17~19
include("trig_int.jl")
# Elliptic Integrals and Jacobian Elliptic Functions
include("exp.jl")

end # PureSpecial