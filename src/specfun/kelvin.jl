# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Kelvin functions"""
#=
- ✅ klvna
- KLVNB
- ✅ klvnzo
=#

"""
    klvna(x::Float64)

Compute Kelvin functions ber x, bei x, ker x
and kei x, and their derivatives  ( x > 0 )

## Input
- `x` --- Argument of Kelvin functions

## Output
- `ber` --- ber x
- `bei` --- bei x
- `ger` --- ker x
- `gei` --- kei x
- `der` --- ber'x
- `dei` --- bei'x
- `her` --- ker'x
- `hei` --- kei'x
"""
function klvna(x::Float64)
    @assert Float64(pi) === 3.141592653589793
    el = 0.5772156649015329
    eps = 1.0e-15

    ber = 1.0
    bei = 0.0
    ger = 1.0e+300
    gei = -0.25 * pi
    der = 0.0
    dei = 0.0
    her = -1.0e+300
    hei = 0.0
    if x == 0.0
        return ber, bei, ger, gei, der, dei, her, hei
    end

    x2 = 0.25 * x * x
    x4 = x2 * x2
    if abs(x) < 10.0
        ber = 1.0
        r = 1.0
        for m = 1:60
            r *= -0.25 / (m * m) / ((2.0 * m - 1.0) * (2.0 * m - 1.0)) * x4
            ber += r
            if abs(r) < (abs(ber) * eps)
                break
            end
        end #= 10 =#

        #= 15 =#
        bei = x2
        r = x2
        for m = 1:60
            r = -0.25 * r / (m * m) / ((2.0 * m + 1.0) * (2.0 * m + 1.0)) * x4
            bei += r
            if abs(r) < (abs(bei) * eps)
                break
            end
        end #= 20 =#

        #= 25 =#
        ger = -(log(x / 2.0) + el) * ber + 0.25 * pi * bei
        r = 1.0
        gs = 0.0
        for m = 1:60
            r = -0.25 * r / (m * m) / ((2.0 * m - 1.0) * (2.0 * m - 1.0)) * x4
            gs += 1.0 / (2.0 * m - 1.0) + 1.0 / (2.0 * m)
            ger += r * gs
            if abs(r * gs) < (abs(ger) * eps)
                break
            end
        end #= 30 =#

        #= 35 =#
        gei = x2 - (log(x / 2.0) + el) * bei - 0.25 * pi * ber
        r = x2
        gs = 1.0
        for m = 1:60
            r = -0.25 * r / (m * m) / ((2.0 * m + 1.0) * (2.0 * m + 1.0)) * x4
            gs += 1.0 / (2.0 * m) + 1.0 / (2.0 * m + 1.0)
            gei += r * gs
            if abs(r * gs) < (abs(gei) * eps)
                break
            end
        end #= 40 =#

        #= 45 =#
        der = -0.25 * x * x2
        r = der
        for m = 1:60
            r = -0.25 * r / m / (m + 1.0) / ((2.0 * m + 1.0) * (2.0 * m + 1.0)) * x4
            der += r
            if abs(r) < (abs(der) * eps)
                break
            end
        end #= 50 =#

        #= 55 =#
        dei = 0.5 * x
        r = dei
        for m = 1:60
            r = -0.25 * r / (m * m) / ((2.0 * m - 1.0) * (2.0 * m + 1.0)) * x4
            dei += r
            if abs(r) < (abs(dei) * eps)
                break
            end
        end #= 60 =#

        #= 65 =#
        r = -0.25 * x * x2
        gs = 1.5
        her = 1.5 * r - ber / x - (log(x / 2.0) + el) * der + 0.25 * pi * dei
        for m = 1:60
            r *= -0.25 / m / (m + 1.0) / ((2.0 * m + 1.0) * (2.0 * m + 1.0)) * x4
            gs += 1.0 / (2.0 * m + 1.0) + 1.0 / (2.0 * m + 2.0)
            her += r * gs
            if abs(r * gs) < (abs(her) * eps)
                break
            end
        end #= 70 =#

        #= 75 =#
        r = 0.5 * x
        gs = 1.0
        hei = 0.5 * x - bei / x - (log(x / 2.0) + el) * dei - 0.25 * pi * der
        for m = 1:60
            r *= -0.25 / (m * m) / ((2 * m - 1.0) * (2 * m + 1.0)) * x4
            gs += 1.0 / (2.0 * m) + 1.0 / (2 * m + 1.0)
            hei += r * gs
            if abs(r * gs) < (abs(hei) * eps)
                return ber, bei, ger, gei, der, dei, her, hei
            end
        end #= 80 =#
    else
        @assert abs(x) >= 10.0

        pp0 = 1.0
        pn0 = 1.0
        qp0 = 0.0
        qn0 = 0.0
        r0 = 1.0
        km = 18
        if abs(x) >= 40.0
            km = 10
        end
        fac = 1.0
        for k = 1:km
            fac = -fac
            xt = 0.25 * k * pi - trunc(0.125 * k) * 2.0 * pi
            cs = cos(xt)
            ss = sin(xt)
            r0 = 0.125 * r0 * ((2.0 * k - 1.0) * (2.0 * k - 1.0)) / k / x
            rc = r0 * cs
            rs = r0 * ss
            pp0 += rc
            pn0 += fac * rc
            qp0 += rs
            qn0 += fac * rs
        end #= 85 =#

        xd = x / sqrt(2.0)
        xe1 = exp(xd)
        xe2 = exp(-xd)
        xc1 = 1.0 / sqrt(2.0 * pi * x)
        xc2 = sqrt(0.5 * pi / x)
        cp0 = cos(xd + 0.125 * pi)
        cn0 = cos(xd - 0.125 * pi)
        sp0 = sin(xd + 0.125 * pi)
        sn0 = sin(xd - 0.125 * pi)

        ger = xc2 * xe2 * (pn0 * cp0 - qn0 * sp0)
        gei = xc2 * xe2 * (-pn0 * sp0 - qn0 * cp0)
        ber = xc1 * xe1 * (pp0 * cn0 + qp0 * sn0) - gei / pi
        bei = xc1 * xe1 * (pp0 * sn0 - qp0 * cn0) + ger / pi

        pp1 = 1.0
        pn1 = 1.0
        qp1 = 0.0
        qn1 = 0.0
        r1 = 1.0
        fac = 1.0
        for k = 1:km
            fac = -fac
            xt = 0.25 * k * pi - trunc(0.125 * k) * 2.0 * pi
            cs = cos(xt)
            ss = sin(xt)
            r1 *= 0.125 * (4.0 - ((2.0 * k - 1.0) * (2.0 * k - 1.0))) / (k * x)
            rc = r1 * cs
            rs = r1 * ss
            pp1 += fac * rc
            pn1 += rc
            qp1 += fac * rs
            qn1 += rs
        end #= 90 =#
        her = xc2 * xe2 * (-pn1 * cn0 + qn1 * sn0)
        hei = xc2 * xe2 * (pn1 * sn0 + qn1 * cn0)
        der = xc1 * xe1 * (pp1 * cp0 + qp1 * sp0) - hei / pi
        dei = xc1 * xe1 * (pp1 * sp0 - qp1 * cp0) + her / pi
    end

    return ber, bei, ger, gei, der, dei, her, hei
end

"""
    klvnzo(nt::Int, kd::Int)

Compute the zeros of Kelvin functions.

## Input
- `NT`  --- Total number of zeros
- `KD`  --- Function code
    - KD=1 to 8 for ber x, bei x, ker x, kei x,
                    ber'x, bei'x, ker'x and kei'x,
                    respectively.

## Output
- ZO(M) --- the M-th zero of Kelvin function for code KD.

## Routine called
- [`Specfun.klvna`](@ref) for computing Kelvin functions and their derivatives.
"""
function klvnzo(nt::Int, kd::Int)
    @assert kd in 1:8
    EPS = 5e-10
    _rt0 = Float64[
        2.84891, 5.02622, 1.71854, 3.91467,
        6.03871, 3.77268, 2.66584, 4.93181,
    ]
    rt = _rt0[kd]

    zo = zeros(Float64, nt)
    for m in 1:nt
        while true
            ber, bei, ger, gei, der, dei, her, hei = klvna(rt)
            if kd == 1
                rt -= ber / der
            elseif kd == 2
                rt -= bei / dei
            elseif kd == 3
                rt -= ger / her
            elseif kd == 4
                rt -= gei / hei
            elseif kd == 5
                rt -= der / (-bei - der / rt)
            elseif kd == 6
                rt -= dei / (ber - dei / rt)
            elseif kd == 7
                rt -= her / (-gei - her / rt)
            else
                rt -= hei / (ger - hei / rt)
            end

            if abs(rt - _rt0[kd]) <= EPS
                break
            else
                _rt0[kd] = rt
            end
        end
        zo[m] = rt
        rt += 4.44
    end

    return zo
end
