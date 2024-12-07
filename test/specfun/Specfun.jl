# SPDX-License-Identifier: MIT
using PureSpecial.Specfun

# for Test
include("warp.jl")

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
