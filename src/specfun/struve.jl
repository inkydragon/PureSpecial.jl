# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Struve Functions"""
#=
- stvh0
- stvh1
- stvhv
- stvl0
- stvl1
- stvlv

Integrals
- ✅ itsh0
- ✅ itth0
- ✅ itsl0
=#

"""
Compute Struve Function H0(x)

Input:
- x   --- Argument of H0(x) ( x ≥ 0 )

Output:
- SH0 --- H0(x)
"""
function stvh0(x::Float64)::Float64
    @assert x >= 0
    _EPS = 1.0e-12
    # TODO: isinf(x), isnan

    s = 1.0
    r = 1.0
    sh0 = 0.0
    if x <= 20.0
        # Use (11.1.4) when x ≤ 20
        a0 = 2.0 * x / pi
        for k in 1:60
            r = -r * (x / (2.0 * k + 1.0))^2
            s += r
            if abs(r / s) < _EPS
                break
            end
        end
        sh0 = a0 * s
    else
        # Use (11.1.21) when x > 20
        km = trunc(Int, 0.5 * (x + 1.0))
        if x >= 50.0
            km = 25
        end
        for k in 1:km
            r = -r * ((2.0 * k - 1.0) / x)^2
            s += r
            if abs(r / s) < _EPS
                break
            end
        end

        # Calculate Y0(x) using (5.2.18)
        t = 4.0 / x
        t2 = t^2
        p0 = (((((-0.37043e-5 * t2 + 0.173565e-4) * t2 - 0.487613e-4) *
                t2 + 0.17343e-3) * t2 - 0.1753062e-2) * t2 + 0.3989422793)
        q0 = t * ((((((0.32312e-5 * t2 - 0.142078e-4) * t2 + 0.342468e-4) *
                t2 - 0.869791e-4) * t2 + 0.4564324e-3) * t2 - 0.0124669441))
        ta0 = x - 0.25 * pi
        by0 = 2.0 / sqrt(x) * (p0 * sin(ta0) + q0 * cos(ta0))
        sh0 = 2.0 / (pi * x) * s + by0
    end

    return sh0
end

"""
Compute Struve Function H1(x)

Input:
- x   --- Argument of H1(x) ( x ≥ 0 )

Output:
- SH1 --- H1(x)
"""
function stvh1(x::Float64)::Float64
    @assert x >= 0
    _EPS = 1.0e-12

    r = 1.0
    sh1 = 0.0
    if x <= 20.0
        # <<< Use (11.1.5) when x ≤ 20
        s = 0.0
        ao = -2.0 / pi
        for k in 1:60
            r = -r * x^2 / (4.0 * k^2 - 1.0)
            s += r
            if abs(r / s) < _EPS
                break
            end
        end
        sh1 = ao * s
    else
        # <<< Use (11.1.22) when x > 20
        s = 1.0
        km = trunc(Int, 0.5 * x)
        if x > 50.0
            km = 25
        end
        for k in 1:km
            r = -r * (4.0 * k^2 - 1.0) / x^2
            s += r
            if abs(r / s) < _EPS
                break
            end
        end

        # <<< Calculate Y1(x) using (5.2.20)
        t = 4.0 / x
        t2 = t^2
        p1 = (((((0.42414e-5 * t2 - 0.20092e-4) * t2 + 0.580759e-4) *
                 t2 - 0.223203e-3) * t2 + 0.29218256e-2) * t2 + 0.3989422819)
        q1 = t * ((((((-0.36594e-5 * t2 + 0.1622e-4) * t2 - 0.398708e-4) *
                     t2 + 0.1064741e-3) * t2 - 0.63904e-3) * t2 + 0.0374008364))
        ta1 = x - 0.75 * pi
        by1 = 2.0 / sqrt(x) * (p1 * sin(ta1) + q1 * cos(ta1))
        sh1 = 2.0 / pi * (1.0 + s / x^2) + by1
    end

    return sh1
end

"""
Compute Struve Functions Hv(x) with
arbitrary order v  ( -8.0 ≤ v ≤ 12.5 )

Input  :
- v  --- Order of Hv(x)
- x  --- Argument of Hv(x) ( x ≥ 0 )

Output :
- Hv --- Hv(x)

Required:
- `gamma2` function to compute the Gamma function
"""
function stvhv(v::Float64, x::Float64)::Float64
    @assert -8.0 <= v <= 12.5
    @assert x >= 0
    _EPS = 1.0e-12

    # Special case when x == 0
    if x == 0.0
        if (v > -1.0) || (v - trunc(v) == 0.5)
            return 0.0
        elseif v < -1.0
            return (-1.0)^(trunc(0.5 - v) - 1) * 1.0e300
        elseif v == -1.0
            return 2.0 / pi
        end
    end

    hv = 0.0
    if x <= 20.0
        # Use (11.1.3) when x ≤ 20
        vo = v + 1.5
        ga = gamma2(vo)
        s = 2.0 / (sqrt(pi) * ga)
        r1 = 1.0
        for k in 1:100
            va = k + 1.5
            ga = gamma2(va)
            vb = v + k + 1.5
            gb = gamma2(vb)
            r1 *= -(0.5 * x)^2
            r2 = r1 / (ga * gb)
            s += r2
            if abs(r2 / s) < _EPS
                break
            end
        end
        hv = (0.5 * x)^(v + 1.0) * s
    else
        # Use (11.1.20) when x > 20
        sa = (0.5 * x)^(v - 1.0) / pi
        vo = v + 0.5
        ga = gamma2(vo)
        s = sqrt(pi) / ga
        r1 = 1.0
        for k in 1:12
            va = k + 0.5
            ga = gamma2(va)
            vb = -k + v + 0.5
            gb = gamma2(vb)
            r1 = r1 / (0.5 * x)^2
            s += ga / gb * r1
        end
        s0 = sa * s

        # Forward recurrence for Yv(x)
        u = abs(v)
        n = trunc(Int, u)
        u0 = u - n
        pu0, qu0, pu1, qu1 = 0.0, 0.0, 0.0, 0.0
        for l in 0:1
            vt = 4.0 * (u0 + l)^2
            r1, r2 = 1.0, 1.0
            pu1 = 1.0
            for k in 1:12
                r1 = -0.0078125 * r1 * (vt - (4.0 * k - 3.0)^2) *
                    (vt - (4.0 * k - 1.0)^2) / ((2.0 * k - 1.0) * k * x^2)
                pu1 += r1
            end
            qu1 = 1.0
            for k in 1:12
                r2 = -0.0078125 * r1 * (vt - (4.0 * k - 1.0)^2) *
                    (vt - (4.0 * k + 1.0)^2) / ((2.0 * k + 1.0) * k * x^2)
                qu1 += r2
            end
            qu1 = 0.125 * (vt - 1.0) / x * qu1
            if l == 0
                pu0, qu0 = pu1, qu1
            end
        end

        t0 = x - (0.5 * u0 + 0.25) * pi
        t1 = x - (0.5 * u0 + 0.75) * pi
        sr = sqrt(2.0 / (pi * x))
        by0 = sr * (pu0 * sin(t0) + qu0 * cos(t0))
        by1 = sr * (pu1 * sin(t1) + qu1 * cos(t1))

        bf = 0.0
        bf0, bf1 = by0, by1
        for k in 2:n
            bf = 2.0 * (k - 1.0 + u0) / x * bf1 - bf0
            bf0, bf1 = bf1, bf
        end

        byv = if n == 0
            by0
        elseif n == 1
            by1
        else
            bf
        end
        hv = byv + s0
    end

    return hv
end


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
        # CoSF (11.1.4) Power-Series Expansions
        # CoSF (11.1.15)
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
        # CoSF 11.1.2 Asymptotic Expansions
        #   When |z| -> ∞, |arg z| < π
        # CoSF (11.1.23)
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
        # CoSF 7.1.5: integrals of Y0(t)
        #   xref: CoSF SUBROUTINE ITJYA(X,TJ,TY)
        # Calculate a(k)
        a0 = 1.0
        a1 = 5.0 / 8.0
        a[1] = a1
        for k in 1:20
            af = (1.5 * (k + 0.5) * (k + 5.0 / 6.0) * a1 - 0.5 * (k + 0.5) * (k + 0.5) * (k - 0.5) * a0) / (k + 1.0)
            a[k + 1] = af
            a0 = a1
            a1 = af
        end
        # Calculate f(x)
        bf = 1.0
        r = 1.0
        for k in 1:10
            r = -r / (x * x)
            bf += a[2 * k] * r
        end
        # Calculate g(x)
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

"""
Evaluate the integral H0(t)/t with respect to t
from x to infinity

Input :
- x   --- Lower limit  ( x ≥ 0 )

Output:
- TTH --- Integration of H0(t)/t from x to infinity
"""
function itth0(x::Float64)
    @assert x >= 0
    _EPS = 1e-12

    s = 1.0
    r = 1.0
    if x < 24.5
        # CoSF (11.1.4) Power-Series Expansions
        # CoSF (11.1.16)
        for k in 1:60
            r = -r * x * x * (2.0 * k - 1.0) / (2.0 * k + 1.0)^3
            s += r
            if abs(r) < abs(s) * _EPS
                break
            end
        end

        tth = pi / 2.0 - 2.0 / pi * x * s
    else
        # CoSF 11.1.2 Asymptotic Expansions
        #   When |z| -> ∞, |arg z| < π
        # CoSF (11.1.24)
        for k in 1:10
            r = -r * (2.0 * k - 1.0)^3 / ((2.0 * k + 1.0) * x * x)
            s += r
            if abs(r) < abs(s) * _EPS
                break
            end
        end

        tth = 2.0 / (pi * x) * s
        # CoSF (7.1.14) Integrals of Y0(t)/t over the Interval (z, +Inf)
        #   xref: CoSF SUBROUTINE ITTJYB(X,TTJ,TTY)
        t = 8.0 / x
        xt = x + 0.25 * pi
        f0 = (((((0.0018118 * t - 0.0091909) * t + 0.017033) * t - 0.0009394) * t - 0.051445) * t - 0.0000011) * t + 0.7978846
        g0 = (((((-0.0023731 * t + 0.0059842) * t + 0.0024437) * t - 0.0233178) * t + 0.0000595) * t + 0.1620695) * t
        tty = (f0 * sin(xt) - g0 * cos(xt)) / sqrt(x) / x

        tth += tty
    end

    return tth
end

"""
Evaluate the integral of modified Struve function
L0(t) with respect to t from 0 to x

Input :
- x   --- Upper limit  ( x ≥ 0 )

Output:
- TL0 --- Integration of L0(t) from 0 to x
"""
function itsl0(x::Float64)
    @assert x >= 0
    _EPS = 1e-12
    el = 0.57721566490153

    r = 1.0
    if x <= 20.0
        # CoSF (11.2.1); DLMF 11.2.2: Power-Series Expansions
        # CoSF (11.2.13)
        s = 0.5
        for k in 1:100
            rd = ifelse(k == 1, 0.5, 1.0)
            r = r * rd * k / (k + 1.0) * (x / (2.0 * k + 1.0))^2
            s += r
            if abs(r / s) < _EPS
                break
            end
        end

        tl0 = 2.0 / pi * x * x * s
    else
        # CoSF 11.2.2 Asymptotic Expansions
        #   When |z| -> ∞, |arg z| < π/2
        # CoSF (11.2.19)
        s = 1.0
        for k in 1:10
            r = r * k / (k + 1.0) * ((2.0 * k + 1.0) / x)^2
            s += r
            if abs(r / s) < _EPS
                break
            end
        end
        s0 = -s / (pi * x * x) + 2.0 / pi * (log(2.0 * x) + el)

        # CoSF (7.2.4) Integrals of I0(t) over the Interval (0, z)
        #   xref: CoSF SUBROUTINE ITIKA(X,TI,TK)
        # Calculate a(k)
        a = zeros(Float64, 18)
        a0, a1 = 1.0, 5.0 / 8.0
        a[1] = a1
        for k in 1:10
            af = ((1.5 * (k + 0.5) * (k + 5.0 / 6.0) * a1 - 0.5 * (k + 0.5)^2 * (k - 0.5) * a0)) / (k+1)
            a[k+1] = af
            a0, a1 = a1, af
        end
        # CoSF (7.2.4)
        ti = 1.0
        r = 1.0
        for k in 1:11
            r /= x
            ti += + a[k] * r
        end
        ti = ti / sqrt(2 * pi * x) * exp(x)

        tl0 = ti + s0
    end

    return tl0
end
