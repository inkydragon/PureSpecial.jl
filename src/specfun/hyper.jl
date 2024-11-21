# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Hypergeometric functions"""
#=
- HYGFX(GAMMA,PSI), 376
- HYGFZ(GAMMA,PSI), 380

- ✅ cchg
- ✅ chgm
- ✅ chgm_kernel
- chgu
    - chgubi
    - chgus,chguit
    - chgul
=#

"""
    cchg(a::Float64, b::Float64, z::Complex{Float64})

Compute confluent hypergeometric function
`M(a,b,z)` with real parameters `a`, `b` and
a complex argument `z`.

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
        return complex(1e300)
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
                if abs((chg - chw) / chg) < 1e-15
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

Compute confluent hypergeometric function M(a,b,x)

Input
- `a`  --- Parameter
- `b`  --- Parameter ( b <> 0,-1,-2,... )
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
                if hg != 0.0 && abs(rg / hg) < 1e-15
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
Purpose: Compute the confluent hypergeometric function
U(a,b,x) for large argument x

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
    _EPS = 1e-15

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
