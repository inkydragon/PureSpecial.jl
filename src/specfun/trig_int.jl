# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""17 Cosine and Sine Integrals"""
#=
- CISIA, 647
- ✅ cisib
=#

"""
    cisib(x::Float64)

Compute cosine and sine integrals Ci(x) and Si(x) for x ≥ 0.

Parameters:
- `x`: Argument of Ci(x) and Si(x)

Returns: `(ci, si)`
- Ci(x)
- Si(x)
"""
function cisib(x::Float64)
    x2 = x * x
    ci = 0.0
    si = 0.0

    if x == 0.0
        ci = -1.0e300
        si = 0.0
    elseif x <= 1.0
        ci = ((((-3.0e-8 * x2 + 3.1e-6) * x2 - 2.3148e-4) *
                    x2 + 1.041667e-2) * x2 - 0.25) * x2 + 0.577215665 + log(x)
        si = ((((3.1e-7 * x2 - 2.834e-5) * x2 + 1.66667e-3) *
                    x2 - 5.555556e-2) * x2 + 1.0) * x
    else
        fx = ((((x2 + 38.027264) * x2 + 265.187033) * x2 + 335.67732) * x2 + 38.102495) /
             ((((x2 + 40.021433) * x2 + 322.624911) * x2 + 570.23628) * x2 + 157.105423)

        gx = ((((x2 + 42.242855) * x2 + 302.757865) * x2 + 352.018498) * x2 + 21.821899) /
             ((((x2 + 48.196927) * x2 + 482.485984) * x2 + 1114.978885) * x2 + 449.690326) /
             x

        ci = fx * sin(x) / x - gx * cos(x) / x
        si = 1.570796327 - fx * cos(x) / x - gx * sin(x) / x
    end

    return ci, si
end
