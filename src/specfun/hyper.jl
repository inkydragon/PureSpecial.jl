# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Hypergeometric functions

hyp2f1(a, b, c, z[, out])
    Gauss hypergeometric function 2F1(a, b; c; z)
    + hyp2f1 -> cephes_hyp2f1 @cephes/hyp2f1.c
    + hyp2f1_complex @_hyp2f1.pxd

hyp1f1(a, b, x[, out])
    Confluent hypergeometric function 1F1.
    + hyp1f1_double -> hyp1f1_wrap -> <hyp1f1_wrap> -> boost::math::hypergeometric_1F1
    +               -> hyp1f1_wrap -> specfun_chgm
    + chyp1f1_wrap -> specfun_cchg

hyperu(a, b, x[, out])
    Confluent hypergeometric function U
    + hyperu -> poch -> cephes_poch
             -> hypU_wrap -> specfun_chgu

hyp0f1(v, z[, out])
    Confluent hypergeometric limit function 0F1.
    + _hyp0f1_real -> _hyp0f1_asy@cython
    + _hyp0f1_cmplx -> cbesi_wrap @amos_wrappers.c
                    -> cbesj_wrap @amos_wrappers.c
"""

"Coefficients for the series expansion"
const _CGAMMA_A = Float64[
    8.333333333333333e-02, -2.777777777777778e-03,
    7.936507936507937e-04, -5.952380952380952e-04,
    8.417508417508418e-04, -1.917526917526918e-03,
    6.410256410256410e-03, -2.955065359477124e-02,
    1.796443723688307e-01, -1.392432216905900e+00
]

"""
    cgama(z::Complex{Float64}, kf::Int)

Compute the gamma function Г(z) or ln[Г(z)]
for a complex argument.

Input
z  --- Complex argument
kf --- Function code
       kf=0 for ln[Г(z)]
       kf=1 for Г(z)

Output
g  --- ln[Г(z)] or Г(z)
"""
function cgama(z::Complex{Float64}, kf::Int)
    @assert kf in [0, 1] "Only accept `kf=0` or `kf=1`"

    x = real(z)
    y = imag(z)

    # Check if z is a negative real integer
    if y == 0.0 && x <= 0.0 && isinteger(x)
        # Inf
        return complex(1e300)
    end

    x1 = x
    y1 = y
    if x < 0.0
        x = -x
        y = -y
        z = -z
        @assert z === complex(x, y)
    else
        x1 = x
        y1 = 0.0
    end

    x0 = x
    na = 0
    if x <= 7.0
        na = trunc(Int64, 7 - x)
        x0 += na
    end

    az0 = abs(complex(x0, y))
    th = atan(y / x0)
    gr = (x0 - 0.5)*log(az0) - th*y - x0 + 0.5*log(2.0 * pi)
    gi = th * (x0 - 0.5) + y*log(az0) - y

    for k in 1:10
        t = az0 ^ (1 - 2*k)
        gr += _CGAMMA_A[k] * t * cos((2.0*k - 1.0) * th)
        gi -= _CGAMMA_A[k] * t * sin((2.0*k - 1.0) * th)
    end

    if x <= 7.0
        gr1 = 0.0
        gi1 = 0.0
        for j in 0:(na-1)
            gr1 += 0.5 * log((x + j)^2 + y^2)
            gi1 += atan(y / (x + j))
        end
        gr -= gr1
        gi -= gi1
    end

    if x1 < 0.0
        az0 = abs(z)
        th1 = atan(y / x)
        sr = -sin(pi*x) * cosh(pi*y)
        si = -cos(pi*x) * sinh(pi*y)
        az1 = abs(complex(sr, si))
        th2 = atan(si / sr)
        if sr < 0.0
            th2 += pi
        end
        gr = log(pi / (az0 * az1)) - gr
        gi = -th1 - th2 - gi
        z = complex(x1, y1)
    end

    if kf == 1
        # Г(z)
        #= exp(gr) * complex(cos(gi), sin(gi)) =#
        return exp(gr) * cis(gi)
    else
        # ln[ Г(z) ]
        return complex(gr, gi)
    end
end
cgama(f::Float64, kf) = cgama(complex(f), kf)

"""
    cchg(a::Float64, b::Float64, z::Complex{Float64})

Compute confluent hypergeometric function
`M(a,b,z)` with real parameters `a`, `b` and
a complex argument `z`.

## Input
a --- Parameter
b --- Parameter
z --- Complex argument

## Output
CHG --- M(a,b,z)

## Routine called
CGAMA for computing complex ln[Г(x)]

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
a  --- Parameter
b  --- Parameter ( b <> 0,-1,-2,... )
x  --- Argument
         
Output
HG --- M(a,b,x)

Routine called
CGAMA for computing complex ln[Г(x)]
"""
function chgm(a::Float64, b::Float64, x::Float64)
    # TODO: merge chgm && chgm_kernel, when removeing specfun.f tests.

    #= Check for special cases =#
    # b = 0, -1, -2, ...
    if isinteger(b) && b <= 0.0
        return Inf
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
