# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
#=exp.jl
    + e1xb
=#


"""
    e1xb(x::Float64)

Compute exponential integral E1(x)

Input
x  --- Argument of E1(x)

Output
e1 --- E1(x)  ( x > 0 )
"""
function e1xb(x::Float64)
    @assert x >= 0

    EPS = 1e-15
    ga = 0.5772156649015328

    e1 = NaN
    if x == 0.0
        # Inf
        e1 = 1e300
    elseif x <= 1.0
        e1 = 1.0
        r = 1.0
        for k = 1:25
            r *= -k*x / (k+1)^2
            e1 += r
            if abs(r) <= (abs(e1)*EPS)
                break
            end
        end
        e1 = -ga - log(x) + x*e1
    else
        m = 20 + trunc(Int64, 80.0 / x)
        t0 = 0.0
        for k = m:-1:1
            t0 = k / (1.0 + k / (x + t0))
        end
        t = 1.0 / (x + t0)
        e1 = exp(-x) * t
    end

    return e1
end
