# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Struve Functions"""
#=
- itsh0
- itth0
- itsl0
=#

"""
Evaluate the integral of Struve function
H0(t) with respect to t from 0 and x

Input :
- x   --- Upper limit  ( x ≥ 0 )

Output:
- TH0 --- Integration of H0(t) from 0 and x
"""
function itsh0(x::Float64)
    @assert x >= 0
    _EPS = 1e-12

    if x <= 30.0
        r = 1.0
        s = 0.5
        for k in 1:100
            rd = ifelse(k == 1, 0.5, 1.0)
            r = -r * rd * k / (k + 1.0) * (x / (2.0 * k + 1.0))^2
            s += r
            if abs(r) < abs(s) * _EPS
                break
            end
        end

        th0 = 2.0 / pi * x * x * s
        return th0
    else
        r = 1.0
        s = 1.0
        for k in 1:12
            r = -r * k / (k + 1.0) * ((2.0 * k + 1.0) / x)^2
            s += r
            if abs(r) < abs(s) * _EPS
                break
            end
        end

        el = 0.57721566490153
        a = zeros(Float64, 25)
        s0 = s / (pi * x^2) + 2.0 / pi * (log(2.0 * x) + el)
        a0 = 1.0
        a1 = 5.0 / 8.0
        a[1] = a1
        for k in 1:20
            af = (1.5 * (k + 0.5) * (k + 5.0 / 6.0) * a1 - 0.5 * (k + 0.5) * (k + 0.5) * (k - 0.5) * a0) / (k + 1.0)
            a[k + 1] = af
            a0 = a1
            a1 = af
        end

        bf = 1.0
        r = 1.0
        for k in 1:10
            r = -r / (x * x)
            bf += a[2 * k] * r
        end

        bg = a[1] / x
        r = 1.0 / x
        for k in 1:10
            r = -r / (x * x)
            bg += a[2 * k + 1] * r
        end

        xp = x + 0.25 * pi
        ty = sqrt(2.0 / (pi * x)) * (bg * cos(xp) - bf * sin(xp))
        th0 = ty + s0

        return th0
    end
end
