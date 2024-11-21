# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Hypergeometric functions"""
#=
- HYGFX(GAMMA,PSI), 376
- HYGFZ(GAMMA,PSI), 380

- ✅ chgm
    - ✅ chgm_kernel
- ✅ cchg
- ✅ chgu
    - ✅ chgus
    - ✅ chgul
    - ✅ chgubi
    - ✅ chguit
=#

"""
    cchg(a::Float64, b::Float64, z::Complex{Float64})

Compute confluent hypergeometric function `M(a, b, z)`
with real parameters `a`, `b` and a complex argument `z`.

## Input
- `a` --- Parameter
- `b` --- Parameter
- `z` --- Complex argument

## Output
- M(a,b,z)

## Routine called
- [`Specfun.cgama`](@ref) for computing complex ln[Г(x)]

## Reference
- [DLMF C13: Confluent Hypergeometric Functions](https://dlmf.nist.gov/13)
- Zhang, S.J., & Jin, J.M. (1996). Computation of Special Functions. Wiley.
"""
function cchg(a::Float64, b::Float64, z::Complex{Float64})
    ci = 0.0 + 1.0im

    # Variables initialization
    a0 = a
    a1 = a
    z0 = z

    #= Check for special cases =#
    # b = 0, -1, -2, ...
    if isinteger(b) && b <= 0.0
        # TODO: warp and ret Inf
        return complex(SF_INF300)
    end
    # a = 0 OR z = 0
    if a == 0.0 || z == complex(0.0)
        # DLMF 13.6.3:  M(0,b,z) = ... = 1
        # CoSF 12.8.5:  M(a,b,0) = M(0,b,z) = 1
        return complex(1.0)
    end
    if a == -1.0
        return 1.0 - z / b
    end
    if a == b
        # DLMF 13.6.1:  M(a,a,z) = e^z
        return exp(z)
    end
    # a = b + 1
    if (a - b) == 1.0
        return (1.0 + z / b) * exp(z)
    end
    if a == 1.0 && b == 2.0
        return (exp(z) - 1.0) / z
    end

    # a is negative integer
    if isinteger(a) && a < 0.0
        # CoSF 12.8.7:  Degenerating Form.
        #   M(-m,b,z), m > 0
        m = trunc(Int64, -a)
        @assert m > 0
        cr = complex(1.0)
        chg = complex(1.0)
        for k in 1:m
            #      (-m)_k       / k! / (b)_k        * z^k      
            cr *= (a + k - 1.0) / k / (b + k - 1.0) * z
            chg += cr
        end
        return chg
    end

    x0 = real(z)
    if x0 < 0.0
        # Preparing DLMF 13.2.39:  Kummer’s Transformations
        #   M(a,b,z) = exp(z)*M(b-a,b,-z)
        a = b - a
        a0 = a
        z = -z
    end

    nl = 0
    la = 0
    if a >= 2.0
        # Preparing terms for DLMF 13.3.1
        # xref: "Applying DLMF 13.3.1"
        nl = 1
        la = trunc(Int64, a)
        a -= la + 1
    end

    cy0 = 0.0
    cy1 = 0.0
    ns = 0
    for n in 0:nl
        if a0 >= 2.0
            a += 1.0
        end
        if (abs(z) < 20.0 + abs(b)) || (a < 0.0)
            # CoSF 12.7.7
            chw = complex(0.0)
            chg = complex(1.0)
            crg = complex(1.0)
            for j in 1:500
                crg *= (a + j - 1.0) / (j * (b + j - 1.0)) * z
                chg += crg
                if abs((chg - chw) / chg) < SF_EPS15
                    break
                end

                chw = chg
            end
        else
            # When |z| --> Inf
            # CoSF 12.8.11:  Asymptotic Formulas
            y = 0.0
            cg1 = cgama(a, 0)
            cg2 = cgama(b, 0)
            cg3 = cgama(b - a, 0)
            
            cs1 = complex(1.0)
            cs2 = complex(1.0)
            cr1 = complex(1.0)
            cr2 =complex(1.0)
            for i in 1:8
                cr1 *= -(a + i - 1.0) * (a - b + i) / (z * i)
                cr2 *= (b - a + i - 1.0) * (i - a) / (z * i)
                cs1 += cr1
                cs2 += cr2
            end
            
            x = real(z)
            y = imag(z)
            phi= 0.0
            if x == 0.0 && y >= 0.0
                phi = 0.5 * pi
            elseif x == 0.0 && y <= 0.0
                phi = -0.5 * pi
            else
                phi = atan(y / x)
            end

            # -π/2 < arg(Z) < 3π/2
            if (phi > -0.5*pi) && (phi < 1.5*pi)
                ns = 1
            end
            # -3π/2 < arg(Z) <= -π/2
            if (phi > -1.5*pi) && (phi <= -0.5*pi)
                ns = -1
            end

            cfac = exp(ns * ci * pi * a)
            if y == 0.0
                cfac = cos(pi * a)
            end
            chg1 = exp(cg2 - cg3) * z ^ (-a) * cfac * cs1
            chg2 = exp(cg2 - cg1 + z) * z ^ (a - b) * cs2
            chg = chg1 + chg2
        end
        if n == 0
            cy0 = chg
        end
        if n == 1
            cy1 = chg
        end
    end

    if a0 >= 2.0
        # Applying DLMF 13.3.1:
        #   (b-a)*M(a-1,b,z) + (2a-b+z)*M(a,b,z) - a*M(a+1,b,z) = 0
        #
        #   M(a+1,b,z) = ( (2a-b+z)*cy1 + (b-a)*cy0 ) / a
        #       cy1 = M(a-1,b,z)
        #       cy0 = M(a,b,z)
        #
        for _ in 1:(la-1)
            chg = ((2.0 * a - b + z) * cy1 + (b - a) * cy0) / a
            cy0 = cy1
            cy1 = chg
            a += 1.0
        end
    end

    if x0 < 0.0
        # DLMF 13.2.39:  Kummer’s Transformations
        #   M(a,b,z) = exp(z)*M(b-a,b,-z)
        #
        # xref: "Preparing DLMF 13.2.39"
        chg *= exp(-z)
    end

    # TODO: remove this
    a = a1
    z = z0
    return chg
end

"""
    chgm(a::Float64, b::Float64, x::Float64)

Compute confluent hypergeometric function `M(a, b, x)`.

Input
- `a`  --- Parameter
- `b`  --- Parameter ( b != 0,-1,-2,... )
- `x`  --- Argument

Output
- M(a,b,x)

Routine called
- [`Specfun.cgama`](@ref) for computing complex ln[Г(x)]
"""
function chgm(a::Float64, b::Float64, x::Float64)
    # TODO: merge chgm && chgm_kernel, when removeing specfun.f tests.

    #= Check for special cases =#
    # b = 0, -1, -2, ...
    if isinteger(b) && b <= 0.0
        # TODO: return Inf
        return 1.0e300
    end
    # M(0,b,x) OR M(a,b,0)
    if a == 0.0 || x == 0.0
        # DLMF 13.6.3:  M(0,b,x) = ... = 1
        # CoSF 12.8.5:  M(a,b,0) = M(0,b,x) = 1
        return 1.0
    end
    # M(-1,b,x)
    if a == -1.0
        return 1.0 - x / b
    end
    # M(a,a,x)
    if a == b
        # DLMF 13.6.1:  M(a,a,x) = e^x
        return exp(x)
    end
    # M(b+1,b,x)
    if (a - b) == 1.0
        return (1.0 + x / b) * exp(x)
    end
    # M(1,2,x)
    if a == 1.0 && b == 2.0
        return (exp(x) - 1.0) / x
    end
    # M(-m,b,x)
    # a is negative integer
    if isinteger(a) && a < 0.0
        # CoSF 12.8.7:  Degenerating Form.
        #   M(-m,b,z), m > 0
        m = trunc(Int64, -a)
        @assert m > 0
        r = 1.0
        hg = 1.0
        for k in 1:m
            #      (-m)_k       / k! / (b)_k        * z^k      
            r *= (a + k - 1.0) / k / (b + k - 1.0) * z
            hg += r
        end
        return hg
    end

    chgm_kernel(a, b, x)
end

"""
    chgm_kernel(a::Float64, b::Float64, x::Float64)

F77 impl in scipy, without input check.
"""
function chgm_kernel(a::Float64, b::Float64, x::Float64)
    a0 = a
    x0 = x
    hg = 0.0

    # DLMF 13.2.39:  M(a,b,x) = exp(x)*M(b-a,b,-x)
    if x < 0.0
        a = b - a
        a0 = a
        x = -x
    end

    nl = 0
    la = 0
    if a >= 2.0
        # preparing terms for DLMF 13.3.1
        nl = 1
        # (la-1): Number of iterations
        la = trunc(Int64, a)
        a -= la + 1
        @assert -1.0 <= a < 0.0
    end

    y0 = 0.0
    y1 = 0.0
    @assert nl in 0:1
    for n = 0:nl
        if a0 >= 2.0
            a += 1.0
        end

        if x <= (30.0 + abs(b)) || (a < 0.0)
            hg = 1.0
            rg = 1.0
            for j = 1:500
                rg *= (a + j - 1.0) / (j * (b + j - 1.0)) * x
                hg += rg
                if hg != 0.0 && abs(rg / hg) < SF_EPS15
                    # DLMF 13.2.39:  M(a,b,z) = exp(z)*M(b-a,b,-z)
                    #   (cf. above)
                    if x0 < 0.0
                        hg *= exp(x0)
                    end
                    break
                end
            end
        else
            #
            # DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum
            #
            cta = cgama(a, 0)
            ctb = cgama(b, 0)
            xg = b - a
            ctba = cgama(xg, 0)
            sum1 = 1.0
            sum2 = 1.0
            r1 = 1.0
            r2 = 1.0
            for i = 1:8
                r1 *= -(a + i - 1.0) * (a - b + i) / (x * i)
                r2 *= -(b - a + i - 1.0) * (a - i) / (x * i)
                sum1 += r1
                sum2 += r2
            end
            if x0 >= 0.0
                hg1 = real(exp(ctb - ctba) * x^(-a) * cos(pi*a) * sum1)
                hg2 = real(exp(ctb - cta + x) * x^(a - b) * sum2)
            else
                # 
                # DLMF 13.2.39 (cf. above)
                #
                hg1 = real(exp(ctb - ctba + x0) * x^(-a) * cos(pi*a) * sum1)
                hg2 = real(exp(ctb - cta) * x^(a - b) * sum2)
            end
            hg = hg1 + hg2
        end

        if n == 0
            y0 = hg
        end
        if n == 1
            y1 = hg
        end
    end

    if a0 >= 2.0
        # DLMF 13.3.1:
        #   (b-a)*M(a-1,b,z) + (2a-b+z)*M(a,b,z) - a*M(a+1,b,z) = 0
        #
        #   M(a+1,b,z) = ( (2a-b+z)*y1 + (b-a)*y0 ) / a
        #       y0 = M(a-1,b,z)
        #       y1 = M(a,b,z)
        #
        for _ = 1:(la-1)
            hg = ((2.0 * a - b + x) * y1 + (b - a) * y0) / a
            y0 = y1
            y1 = hg
            a += 1.0
        end
    end

    return hg
end

"""
    chgul(a::Float64, b::Float64, x::Float64)

Compute the confluent hypergeometric function `U(a, b, x)`
for large argument `x`.

Input:
- `a`  --- Parameter
- `b`  --- Parameter
- `x`  --- Argument, x > 0

Output: `(hu, id)`
- `HU` --- U(a,b,x)
- `ID` --- Estimated number of significant digits
"""
function chgul(a::Float64, b::Float64, x::Float64)
    @assert x > 0
    _EPS = SF_EPS15

    id = -100
    aa = a - b + 1.0
    il1 = (a == trunc(Int, a)) && (a <= 0.0)
    il2 = (aa == trunc(Int, aa)) && (aa <= 0.0)
    nm = 0
    if il1
        nm = trunc(Int, abs(a))
    elseif il2
        nm = trunc(Int, abs(aa))
    end

    if il1 || il2
        # IL1: DLMF 13.2.7 with k = -s-a
        # IL2: DLMF 13.2.8
        r = 1.0
        hu = 1.0
        for k in 1:nm
            r = - r * (a + k - 1.0) * (a - b + k) / (k * x)
            hu += r
        end
        hu = x^(-a) * hu
        id = 10
    else
        r = 1.0
        hu = 1.0
        ra = NaN
        r0 = NaN
        # DLMF 13.7.3
        for k in 1:25
            r = - r * (a + k - 1.0) * (a - b + k) / (k * x)
            ra = abs(r)
            if ((k > 5) && (ra >= r0)) || (ra < _EPS)
                break
            end
            r0 = ra
            hu += r
        end
        hu = x^(-a) * hu
        id = trunc(Int, abs(log10(ra)))
    end

    return hu, id
end

"""
    chgus(a::Float64, b::Float64, x::Float64)

Compute confluent hypergeometric function `U(a, b, x)` for small argument `x`.

Input:
- `a`  --- Parameter
- `b`  --- Parameter ( b ≠ 0, -1, -2, ...)
- `x`  --- Argument

Output: `(hu, id)`
- `HU` --- U(a,b,x)
- `ID` --- Estimated number of significant digits

Routine called:
- [`Specfun.gamma2`](@ref) for computing gamma function
"""
function chgus(a::Float64, b::Float64, x::Float64)
    @assert !isinteger(b)
    _EPS = SF_EPS15
    # DLMF 13.2.42 with prefactors rewritten according to
    # DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2
    ga = gamma2(a)
    gb = gamma2(b)
    gab = gamma2(1.0 + a - b)
    gb2 = gamma2(2.0 - b)
    hu0 = pi / sin(pi * b)

    r1 = hu0 / (gab * gb)
    r2 = hu0 * x^(1.0 - b) / (ga * gb2)
    hu = r1 - r2
    hmax = 0.0
    hmin = SF_INF300
    h0 = 0.0
    for j in 1:150
        r1 = r1 * (a + j - 1.0) / (j * (b + j - 1.0)) * x
        r2 = r2 * (a - b + j) / (j * (1.0 - b + j)) * x
        hu = hu + r1 - r2
        hua = abs(hu)
        hmax = max(hmax, hua)
        hmin = min(hmin, hua)
        if abs(hu - h0) < abs(hu) * _EPS
            break
        end
        h0 = hu
    end

    d1 = log10(hmax)
    d2 = ifelse(hmin != 0.0, log10(hmin), 0.0)
    id = trunc(Int, 15 - abs(d1 - d2))
    return hu, id
end

"""
    chgubi(a::Float64, b::Float64, x::Float64)

Compute confluent hypergeometric function `U(a, b, x)`
with integer `b` ( b = ±1,±2,... ).

Input:
- `a`  --- Parameter
- `b`  --- Parameter
- `x`  --- Argument

Output: `(hu, id)`
- `HU` --- U(a,b,x)
- `ID` --- Estimated number of significant digits

Routines called:
- [`Specfun.gamma2`](@ref) for computing gamma function Г(x)
- [`Specfun.psi`](@ref) for computing psi function
"""
function chgubi(a::Float64, b::Float64, x::Float64)
    # (a + m - 1) > 0 && m >= 1
    @assert a > 0
    @assert isinteger(b)
    _EPS = SF_EPS15

    id = -100
    n = trunc(Int, abs(b - 1))
    rn = 1.0
    rn1 = 1.0
    for j in 1:n
        rn *= j
        if j == (n - 1)
            rn1 = rn
        end
    end

    ps = psi(a)
    ga = gamma2(a)
    ua, ub = 0.0, 0.0
    if b > 0.0
        a0 = a
        a1 = a - n
        a2 = a1
        ga1 = gamma2(a1)
        ua = (-1)^(n - 1) / (rn * ga1)
        ub = rn1 / ga * x^(-n)
    else
        a0 = a + n
        a1 = a0
        a2 = a
        ga1 = gamma2(a1)
        ua = (-1)^(n - 1) / (rn * ga) * x^n
        ub = rn1 / ga1
    end

    r = 1.0
    hm1 = 1.0
    hmax = 0.0
    hmin = SF_INF300
    h0 = 0.0
    for k in 1:150
        r = r * (a0 + k - 1) * x / ((n + k) * k)
        hm1 += r
        hu1 = abs(hm1)
        hmax = max(hmax, hu1)
        hmin = min(hmin, hu1)
        if abs(hm1 - h0) < abs(hm1) * _EPS
            break
        end
        h0 = hm1
    end

    da1 = log10(hmax)
    da2 = ifelse(hmin != 0.0, log10(hmin), 0.0)
    id = trunc(Int, 15 - abs(da1 - da2))
    hm1 *= log(x)

    s0 = 0.0
    for m in 1:n
        if b >= 0
            s0 -= 1.0 / m
        else
            s0 += (1.0 - a) / (m * (a + m - 1))
        end
    end

    hm2 = ps + 2 * SF_EULER_GAMMA + s0
    r = 1.0
    hmax = 0.0
    hmin = 1.0e+300
    for k in 1:150
        s1 = 0.0
        s2 = 0.0
        if b > 0
            for m in 1:k
                s1 -= (m + 2 * a - 2) / (m * (m + a - 1))
            end
            for m in 1:n
                s2 += 1.0 / (k + m)
            end
        else
            for m in 1:(k + n)
                s1 += (1.0 - a) / (m * (m + a - 1))
            end
            for m in 1:k
                s2 += 1.0 / m
            end
        end

        hw = 2 * SF_EULER_GAMMA + ps + s1 - s2
        r = r * (a0 + k - 1) * x / ((n + k) * k)
        hm2 += r * hw
        hu2 = abs(hm2)
        hmax = max(hmax, hu2)
        hmin = min(hmin, hu2)
        if abs(hm2 - h0) < abs(hm2) * _EPS
            break
        end
        h0 = hm2
    end

    db1 = log10(hmax)
    db2 = ifelse(hmin != 0.0, log10(hmin), 0.0)
    id1 = trunc(Int, 15 - abs(db1 - db2))
    id = min(id, id1)

    r = 1.0
    hm3 = ifelse(n == 0, 0.0, 1.0)
    for k in 1:(n - 1)
        r *= (a2 + k - 1.0) / ((k - n) * k) * x
        hm3 += r
    end

    sa = ua * (hm1 + hm2)
    sb = ub * hm3
    hu = sa + sb
    if sa != 0.0
        id1 = trunc(Int, log10(abs(sa)))
    end
    if hu != 0.0
        id2 = trunc(Int, log10(abs(hu)))
    end
    if (sa * sb) < 0.0
        id -= abs(id1 - id2)
    end

    return hu, id
end

const _CHGUIT_T = NTuple{30, Float64}((
    0.0259597723012478, 0.0778093339495366, 0.129449135396945, 0.180739964873425,
    0.231543551376029, 0.281722937423262, 0.331142848268448, 0.379670056576798,
    0.427173741583078, 0.473525841761707, 0.518601400058570, 0.562278900753945,
    0.604440597048510, 0.644972828489477, 0.683766327381356, 0.720716513355730,
    0.755723775306586, 0.788693739932264, 0.819537526162146, 0.848171984785930,
    0.874519922646898, 0.898510310810046, 0.920078476177628, 0.939166276116423,
    0.955722255839996, 0.969701788765053, 0.981067201752598, 0.989787895222222,
    0.995840525118838, 0.999210123227436
))

const _CHGUIT_W = NTuple{30, Float64}((
    0.0519078776312206, 0.0517679431749102, 0.0514884515009810, 0.0510701560698557,
    0.0505141845325094, 0.0498220356905502, 0.0489955754557568, 0.0480370318199712,
    0.0469489884891122, 0.0457343797161145, 0.0443964787957872, 0.0429388928359356,
    0.0413655512355848, 0.0396806954523808, 0.0378888675692434, 0.0359948980510845,
    0.0340038927249464, 0.0319212190192963, 0.0297524915007890, 0.0275035567499248,
    0.0251804776215213, 0.0227895169439978, 0.0203371207294572, 0.0178299010142074,
    0.0152746185967848, 0.0126781664768159, 0.0100475571822880, 0.00738993116334531,
    0.00471272992695363, 0.00202681196887362
))

"""
    chguit(a::Float64, b::Float64, x::Float64)

Compute hypergeometric function `U(a, b, x)` by
using Gaussian-Legendre integration (`n = 60`).

Input:
- `a`  --- Parameter ( a > 0 )
- `b`  --- Parameter
- `x`  --- Argument ( x > 0 )

Output: `(hu, id)`
- `HU` --- U(a,b,z)
- `ID` --- Estimated number of significant digits

Routine called:
- [`Specfun.gamma2`](@ref) for computing Г(x)
"""
function chguit(a::Float64, b::Float64, x::Float64)
    @assert a > 0
    @assert x > 0
    _EPS = SF_EPS09
    t = _CHGUIT_T
    w = _CHGUIT_W

    id = 9
    a1 = a - 1.0
    b1 = b - a - 1.0
    c = 12.0 / x
    hu0 = 0.0

    # First integration (DLMF 13.4.4)
    hu1 = 0.0
    for m in 10:5:100
        hu1 = 0.0
        g = 0.5 * c / m
        d = g
        for j in 1:m
            s = 0.0
            for k in 1:30
                t1 = d + g * t[k]
                t2 = d - g * t[k]
                f1 = exp(-x * t1) * t1^a1 * (1.0 + t1)^b1
                f2 = exp(-x * t2) * t2^a1 * (1.0 + t2)^b1
                s += w[k] * (f1 + f2)
            end
            hu1 += s * g
            d += 2.0 * g
        end
        if abs(1.0 - hu0) < abs(hu1) * _EPS
            break
        end

        hu0 = hu1
    end
    ga = gamma2(a)
    hu1 /= ga

    # Second integration (DLMF 13.4.4, substitution)
    hu2 = 0.0
    for m in 2:2:10
        hu2 = 0.0
        g = 0.5 / m
        d = g
        for j in 1:m
            s = 0.0
            for k in 1:30
                t1 = d + g * t[k]
                t2 = d - g * t[k]
                t3 = c / (1.0 - t1)
                t4 = c / (1.0 - t2)
                f1 = t3^2 / c * exp(-x * t3) * t3^a1 * (1.0 + t3)^b1
                f2 = t4^2 / c * exp(-x * t4) * t4^a1 * (1.0 + t4)^b1
                s += w[k] * (f1 + f2)
            end
            hu2 += s * g
            d += 2.0 * g
        end
        if abs(1.0 - hu0) < abs(hu2) * _EPS
            break
        end

        hu0 = hu2
    end
    hu2 /= ga

    hu = hu1 + hu2
    return hu, id
end

"""
    chgu(a::Float64, b::Float64, x::Float64)

Compute the confluent hypergeometric function `U(a, b, x)`.

Input:
- `a`  --- Parameter
- `b`  --- Parameter
- `x`  --- Argument  ( x > 0 )

Output: `(hu, md, isfer)`
- `HU` --- U(a,b,x)
- `MD` --- Method code
- `ISFER` --- Error flag

Routines called:
- [`Specfun.chgus`](@ref) for small x ( MD=1 )
- [`Specfun.chgul`](@ref) for large x ( MD=2 )
- [`Specfun.chgubi`](@ref) for integer b ( MD=3 )
- [`Specfun.chguit`](@ref) for numerical integration ( MD=4 )
"""
function chgu(a::Float64, b::Float64, x::Float64)
    aa = a - b + 1.0
    md = 0
    isfer = 0

    il1 = (a == trunc(Int, a)) && (a <= 0.0)
    il2 = (aa == trunc(Int, aa)) && (aa <= 0.0)
    il3 = abs(a * (a - b + 1.0)) / x <= 2.0
    bl1 = (x <= 5.0) || ((x <= 10.0) && (a <= 2.0))
    bl2 = (x > 5.0) && (x <= 12.5) && ((a >= 1.0) && (b >= a + 4.0))
    bl3 = (x > 12.5) && (a >= 5.0) && (b >= a + 5.0)
    bn = (b == trunc(Int, b)) && (b != 0.0)

    id = 0
    hu = 0.0
    id1 = -100
    hu1 = 0.0
    if b != trunc(Int, b)
        hu, id1 = chgus(a, b, x)
        md = 1
        if id1 >= 9
            return hu, md, isfer
        end

        hu1 = hu
    end

    if il1 || il2 || il3
        hu, id = chgul(a, b, x)
        md = 2
        if id >= 9
            return hu, md, isfer
        end

        if id1 > id
            md = 1
            id = id1
            hu = hu1
        end
    end

    if a >= 1.0
        if bn && (bl1 || bl2 || bl3)
            hu, id = chgubi(a, b, x)
            md = 3
        else
            hu, id = chguit(a, b, x)
            md = 4
        end
    else
        if b <= a
            b00 = b
            a = a - b + 1.0
            b = 2.0 - b
            hu, id = chguit(a, b, x)
            hu = x^(1.0 - b00) * hu
            md = 4
        elseif bn && (!il1)
            hu, id = chgubi(a, b, x)
            md = 3
        end
    end

    if id < 6
        isfer = 6
    end

    return hu, md, isfer
end
