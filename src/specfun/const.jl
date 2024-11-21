# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Math Constants"""

"""
`+Inf` in `specfun.f`
"""
const SF_INF300 = Float64(1e300)

const SF_EPS15 = Float64(1e-15)
const SF_EPS14 = Float64(1e-14)
const SF_EPS12 = Float64(1e-12)
const SF_EPS11 = Float64(1e-11)
const SF_EPS10 = Float64(1e-10)

"""
The constant π.

Same as: Base.MathConstants.pi

```jldoctest
julia> Float64(Base.MathConstants.pi)
3.141592653589793
```
"""
const SF_PI = Float64(Base.MathConstants.pi)


"""
Euler's Constant γ

Same as: Base.MathConstants.eulergamma

```jldoctest
julia> Float64(Base.MathConstants.eulergamma)
0.5772156649015329
```
"""
const SF_EULER_GAMMA = Float64(Base.MathConstants.eulergamma)

"""
Euler's Constant γ

```jldoctest
julia> Float64(Base.MathConstants.eulergamma)
0.5772156649015329

julia> 0.57721566490153_28  # imprecise
0.5772156649015328
```
"""
const SF_EULER_GAMMA_28 = 0.5772156649015328


"""
Euler's Constant γ

```jldoctest
julia> Float64(Base.MathConstants.eulergamma)
0.5772156649015329

julia> 0.57721566490153_00  # imprecise
0.57721566490153
```
"""
const SF_EULER_GAMMA_00 = 0.57721566490153
