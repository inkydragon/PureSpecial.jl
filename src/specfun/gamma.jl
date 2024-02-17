# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md


const _GAM0_G = Float64[
    1.0e0, 0.5772156649015329e0, 
    -0.6558780715202538e0, -0.420026350340952e-1, 
     0.1665386113822915e0, -0.421977345555443e-1, 
    -0.96219715278770e-2, 0.72189432466630e-2, 
    -0.11651675918591e-2, -0.2152416741149e-3, 
     0.1280502823882e-3, -0.201348547807e-4, 
    -0.12504934821e-5, 0.11330272320e-5, 
    -0.2056338417e-6, 0.61160950e-8, 
     0.50020075e-8, -0.11812746e-8,
     0.1043427e-9, 0.77823e-11, 
    -0.36968e-11, 0.51e-12, 
    -0.206e-13, -0.54e-14,
    0.14e-14
]

"""
    gam0(x)

Compute gamma function Г(x) for small `x`

Input
x  --- Argument of Г(x)  ( |x| ≤ 1 )

Output
ga --- Г(x)
"""
function gam0(x::Float64)
    @assert abs(x) <= 1

    gr = _GAM0_G[end]

    for k = 24:-1:1
        gr = gr * x + _GAM0_G[k]
    end

    return 1.0 / (gr * x)
end

const _GAMMA2_G = Float64[
    _GAM0_G...,
    0.1e-15,
]

"""
    gamma2(x::Float64)

Compute gamma function Г(x)

Input
x  --- Argument of Г(x)
        ( x is not equal to 0,-1,-2,…)

Output
ga --- Г(x)
"""
function gamma2(x::Float64)
    ga = NaN
    if isinteger(x)
        if x > 0.0
            ga = 1.0
            m1 = trunc(Int64, x) - 1
            for k = 2:m1
                ga *= k
            end
        else
            # Inf
            ga = 1e300
        end
    else
        r = 1.0
        if abs(x) > 1.0
            z = abs(x)
            m = trunc(Int64, z)
            for k = 1:m
                r *= (z - k)
            end
            z -= m
        else
            z = x
        end

        gr = _GAMMA2_G[26]
        for k = 25:-1:1
            gr = gr*z + _GAMMA2_G[k]
        end
        ga = 1.0 / (gr * z)

        if abs(x) > 1.0
            ga *= r
            if x < 0.0
                ga = -pi / (x * ga * sin(pi * x))
            end
        end
    end

    return ga
end
