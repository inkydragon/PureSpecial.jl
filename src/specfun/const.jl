# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Math Constants"""

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

note: imprecise
"""
const EULER_GAMMA_28 = 0.57721566490153_28
