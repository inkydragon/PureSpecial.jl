# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""17 Cosine and Sine Integrals"""
#=
- ✅ cisia
- ✅ cisib
=#

"""
    cisia(x::Float64)

Compute cosine and sine integrals Ci(x) and Si(x).

Parameters:
- `x`: Argument of Ci(x) and Si(x), x ≥ 0

Returns: `(ci, si)`
- Ci(x)
- Si(x)
"""
function cisia(x::Float64)
    p2 = π / 2
    el = SF_EULER_GAMMA
    eps = 1.0e-15

    x2 = x * x
    ci = 0.0
    si = 0.0
    if x == 0.0
        ci = -1.0e300
        si = 0.0
    elseif x <= 16.0
        xr = -0.25 * x2
        ci = el + log(x) + xr
        for k in 2:40
            xr = -0.5 * xr * (k - 1) / (k^2 * (2 * k - 1)) * x2
            ci += xr
            if abs(xr) < abs(ci) * eps
                break
            end
        end

        xr = x
        si = x
        for k in 1:40
            xr = -0.5 * xr * (2 * k - 1) / (k * (4 * k^2 + 4 * k + 1)) * x2
            si += xr
            if abs(xr) < abs(si) * eps
                return (ci, si)
            end
        end
    elseif x <= 32.0
        m = trunc(Int, 47.2 + 0.82 * x)
        bj = zeros(Float64, m)
        xa1 = 0.0
        xa0 = 1.0e-100

        for k in m:-1:1
            xa = 4.0 * k * xa0 / x - xa1
            bj[k] = xa
            xa1, xa0 = xa0, xa
        end

        xs = bj[1]
        for k in 3:2:m
            xs += 2.0 * bj[k]
        end

        bj .= bj ./ xs

        xr = 1.0
        xg1 = bj[1]
        for k in 2:m
            xr = 0.25 * xr * (2 * k - 3)^2 / ((k - 1) * (2 * k - 1)^2) * x
            xg1 += bj[k] * xr
        end

        xr = 1.0
        xg2 = bj[1]
        for k in 2:m
            xr = 0.25 * xr * (2 * k - 5)^2 / ((k - 1) * (2 * k - 3)^2) * x
            xg2 += bj[k] * xr
        end

        xcs = cos(x / 2.0)
        xss = sin(x / 2.0)
        ci = el + log(x) - x * xss * xg1 + 2 * xcs * xg2 - 2 * xcs^2
        si = x * xcs * xg1 + 2 * xss * xg2 - sin(x)
    else
        xr = 1.0
        xf = 1.0
        for k in 1:9
            xr = -2.0 * xr * k * (2 * k - 1) / x2
            xf += xr
        end

        xr = 1.0 / x
        xg = xr
        for k in 1:8
            xr = -2.0 * xr * (2 * k + 1) * k / x2
            xg += xr
        end

        ci = xf * sin(x) / x - xg * cos(x) / x
        si = p2 - xf * cos(x) / x - xg * sin(x) / x
    end

    return ci, si
end

"""
    cisib(x::Float64)

Compute cosine and sine integrals Ci(x) and Si(x).

Using polynomial and rational approximations:
- `x <= 1`,         `|eps(x)| <= 1e-7`
- `1 < x <= Inf`,   `|eps_ci(x)| <= 5e-7, |eps_si(x)| <= 3e-7`

Parameters:
- `x`: Argument of Ci(x) and Si(x), x ≥ 0

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
