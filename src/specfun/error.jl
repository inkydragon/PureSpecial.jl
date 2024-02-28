# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md


"""
Compute error function erf(x)
Input:  x   --- Argument of erf(x)
Output: ERR --- erf(x)
"""
function erf(x::Float64)
    EPS = 1.0e-15

    x2 = x * x
    err = NaN
    if abs(x) < 3.5
        er = 1.0
        r = 1.0
        for k = 1:50
            r *= x2 / (k + 0.5)
            er += r
            if abs(r) ≤ abs(er) * EPS
                break
            end
        end
        c0 = 2.0 / sqrt(pi) * x * exp(-x2)
        err = c0 * er
    else
        er = 1.0
        r = 1.0
        for k = 1:12
            r *= -(k - 0.5) / x2
            er += r
        end
        c0 = exp(-x2) / (abs(x) * sqrt(pi))
        err = 1.0 - c0 * er
        if x < 0.0
            err = -err
        end
    end

    return err
end

"""
Compute error function erf(z) for a complex argument (z=x+iy)
Input:  z   --- Complex argument
Output: CER --- erf(z)
"""
function erf(z::Complex{Float64})
    EPS = 1.0e-15
    sqpi = 1.7724538509055160273
    z1 = z
    if real(z) < 0.0
        z1 = -z
    end

    # Cutoff radius R = 4.36; determined by balancing rounding error
    # and asymptotic expansion error, see below.
    #
    # The resulting maximum global accuracy expected is around 1e-8
    #
    cer = complex(NaN)
    if abs(z) <= 4.36
        # Rounding error in the Taylor expansion is roughly
        # ~ R*R * EPSILON * R**(2 R**2) / (2 R**2 Gamma(R**2 + 1/2))
        cs = z1
        cr = z1
        for k = 1:121
            cr *= z1 * z1 / (k + 0.5)
            cs += cr
            if abs(cr / cs) < EPS
                break
            end
        end
        cer = 2.0 * exp(-z * z) * cs / sqpi
    else
        cl = 1.0 / z1
        cr = cl
        # Asymptotic series; maximum K must be at most ~ R^2.
        #
        # The maximum accuracy obtainable from this expansion is roughly
        #
        # ~ Gamma(2R**2 + 2) / (
        #          (2 R**2)**(R**2 + 1/2) Gamma(R**2 + 3/2) 2**(R**2 + 1/2))
        for k = 1:21
            cr *= -(k - 0.5) / (z1 * z1)
            cl += cr
            if abs(cr / cl) < EPS
                break
            end
        end
        cer = 1.0 - exp(-z * z) * cl / sqpi
    end
    # TODO: use copysign
    if real(z) < 0.0
        cer = -cer
    end

    return cer
end

"""
Compute complex Error function erf(z) & erf'(z)
    
Input
z   --- Complex argument of erf(z)

Output
CER --- erf(z)
CDER --- erf'(z)
"""
function cerf(z::Complex{Float64})
    EPS = 1.0e-12

    x = real(z)
    y = imag(z)
    x2 = x * x

    cer = complex(NaN)
    if x <= 3.5
        er = 1.0
        r = 1.0
        w = 0.0

        for k = 1:100
            r *= x2 / (k + 0.5)
            er += r
            if abs(er - w) <= EPS * abs(er)
                break
            end
            w = er
        end

        c0 = 2.0 / sqrt(pi) * x * exp(-x2)
        er0 = c0 * er
        cer = complex(er0, 0.0)
    else
        er = 1.0
        r = 1.0

        for k = 1:12
            r *= -r * (k - 0.5) / x2
            er += r
        end

        c0 = exp(-x2) / (x * sqrt(pi))
        er0 = 1.0 - c0 * er
        cer = complex(er0, 0.0)
    end

    if y == 0.0
        err = real(cer)
        cer = complex(err, 0.0)
    else
        cs = cos(2.0 * x * y)
        ss = sin(2.0 * x * y)
        er1 = exp(-x2) * (1.0 - cs) / (2.0 * pi * x)
        ei1 = exp(-x2) * ss / (2.0 * pi * x)
        er2 = 0.0
        w1 = 0.0

        for n = 1:100
            er2 += exp(-0.25 * n * n) / (n * n + 4.0 * x2) * (2.0 * x - 2.0 * x * cosh(n * y) * cs + n * sinh(n * y) * ss)
            if abs((er2 - w1) / er2) < EPS
                break
            end
            w1 = er2
        end

        c0 = 2.0 * exp(-x2) / pi
        err = real(cer) + er1 + c0 * er2
        ei2 = 0.0
        w2 = 0.0

        for n = 1:100
            ei2 += exp(-0.25 * n * n) / (n * n + 4.0 * x2) * (2.0 * x * cosh(n * y) * ss + n * sinh(n * y) * cs)
            if abs((ei2 - w2) / ei2) < EPS
                break
            end
            w2 = ei2
        end
        cer = complex(err, ei1 + c0 * ei2)
    end

    cder = 2.0 / sqrt(pi) * exp(-z*z)
    return cer, cder
end

function cerzo!(zo::Vector{Complex{Float64}}, nt::Int)
    EPS = 1.0e-11

    for nr in 1:nt
        pu = sqrt(pi * (4.0 * nr - 0.5))
        pv = pi * sqrt(2.0 * nr - 0.25)
        px = 0.5 * pu - 0.5 * log(pv) / pu
        py = 0.5 * pu + 0.5 * log(pv) / pu
        z = complex(px, py)

        it = 0
        w = 0.0
        while true
            it += 1
            zf, zd = cerf(z)
            zp = 1.0
            for i in 1:nr
                zp *= (z - zo[i])
            end
            zfd = zf / zp
            zq = 0.0
            for i in 1:nr
                zw = 1.0
                for j in 1:nr
                    if j != i
                        zw *= (z - zo[j])
                    end
                end
                zq += zw
            end
            zgd = (zd - zq * zfd) / zp
            z -= zfd / zgd
            w0 = w
            w = abs(z)
            if (it > 50) || (abs((w - w0) / w) <= EPS)
                break
            end
        end

        zo[nr] = z
    end
end



"""
Compute complex Fresnel integral C(z) and C'(z)

Input
z --- Argument of C(z)

Output
ZF --- C(z)
ZD --- C'(z)
"""
function cfc(z::Complex{Float64})
    EPS = 1.0e-14

    w0 = abs(z)
    # pi/2 * x^2
    zp = 0.5 * pi * z * z
    # (pi/2 * x^2)^2
    zp2 = zp * zp
    z0 = complex(0.0)

    c = complex(NaN, NaN)
    if z == z0
        # C(0) = 0
        c = z0
    elseif w0 <= 2.5
        # DLMF 7.6.4: Power Series
        cr = z
        c = cr
        wa0 = 0.0
        for k ∈ 1:80
            cr = -0.5 * cr * (4 * k - 3) / k / (2 * k - 1) / (4 * k + 1) * zp2
            c += cr
            wa = abs(c)
            if (abs((wa - wa0) / wa) < EPS) && (k > 10)
                break
            end

            wa0 = wa
        end
    elseif w0 < 4.5
        # CoSF 16.5.5: Expansions in Series of Bessel Functions
        # NOTE: Uses a different algorithm from [CoSF-P632]
        m = 85
        c = z0
        cf1 = z0
        cf0 = complex(1.0e-100)
        cf = complex(NaN, NaN)
        for k ∈ m:-1:0
            cf = (2.0 * k + 3.0) * cf0 / zp - cf1
            if mod(k, 2) == 0
                c += cf
            end
            cf1 = cf0
            cf0 = cf
        end
        c *= 2.0 / (pi * z) * sin(zp) / cf
    else
        # use C(z) = iC(-iz)
        #
        # Auxiliary functions f(z) and g(z) can be computed using an
        # asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
        # as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
        # Interestingly, most of the expansion code is the same across
        # the quadrants. (The forth power in Z is the equalizer here.)
        # Only one constant has to be adapted.
        d = complex(NaN)
        if (imag(z) > -real(z)) && (imag(z) <= real(z))
            # right quadrant
            d = complex(0.5)
        elseif (imag(z) > real(z)) && (imag(z) >= -real(z))
            # upper quadrant
            d = complex(0.5im)
        elseif (imag(z) < -real(z)) && (imag(z) >= real(z))
            # left quadrant
            d = complex(-0.5)
        else
            # lower quadrant
            d = complex(-0.5im)
        end

        # f(X): DLMF 7.12.2
        cr = complex(1.0)
        cf = complex(1.0)
        for k ∈ 1:20
            cr = -0.25 * cr * (4 * k - 1) * (4 * k - 3) / zp2
            cf += cr
        end

        # g(X): DLMF 7.12.3
        cr = 1.0 / (pi * z * z)
        cg = cr
        for k ∈ 1:12
            cr = -0.25 * cr * (4 * k + 1) * (4 * k - 1) / zp2
            cg += cr
        end

        # DLMF 7.5.3:  Asymptotic Expansions
        c = d + (cf * sin(zp) - cg * cos(zp)) / (pi * z)
    end

    zf = c
    # C'(z) = cos(pi/2 * z^2)
    zd = cos(0.5 * pi * z * z)
    return zf, zd
end


"""
Compute complex Fresnel Integral S(z) and S'(z)

Input
z  --- Argument of S(z)

Output
ZF --- S(z)
ZD --- S'(z)
"""
function cfs(z::Complex{Float64})
    EPS = 1.0e-14

    w0 = abs(z)
    # pi/2 * x^2
    zp = 0.5 * pi * z * z
    # (pi/2 * x^2)^2
    zp2 = zp * zp
    z0 = complex(0.0)

    s = complex(NaN, NaN)
    if z == z0
        # S(0) = 0
        s = z0
    elseif w0 <= 2.5
        # DLMF 7.6.6: Power Series
        s = z * zp / 3.0
        cr = s
        wb0 = 0.0
        for k = 1:80
            cr = -0.5 * cr * (4 * k - 1) / k / (2 * k + 1) / (4 * k + 3) * zp2
            s += cr
            wb = abs(s)
            if (abs(wb - wb0) < EPS) && (k > 10)
                break
            end

            wb0 = wb
        end
    elseif w0 < 4.5
        # CoSF 16.5.6: Expansions in Series of Bessel Functions
        # NOTE: Uses a different algorithm from [CoSF-P632]
        m = 85
        s = z0
        cf1 = z0
        cf0 = complex(1.0e-100)
        cf = complex(NaN, NaN)
        for k = m:-1:0
            cf = (2.0 * k + 3.0) * cf0 / zp - cf1
            if mod(k, 2) == 1
                s += cf
            end
            cf1 = cf0
            cf0 = cf
        end
        s = 2.0 / (pi * z) * sin(zp) / cf * s
    else
        # Auxiliary functions f(z) and g(z) can be computed using an
        # asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
        # as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
        # Interestingly, most of the expansion code is the same across
        # the quadrants. (The forth power in Z is the equalizer here.)
        # Only one constant has to be adapted.
        d = complex(NaN)
        if (imag(z) > -real(z)) && (imag(z) <= real(z))
            # right quadrant
            d = complex(0.5)
        elseif (imag(z) > real(z)) && (imag(z) >= -real(z))
            # upper quadrant
            d = complex(-0.5im)
        elseif (imag(z) < -real(z)) && (imag(z) >= real(z))
            # left quadrant
            d = complex(-0.5)
        else
            # lower quadrant
            d = complex(0.5im)
        end

        # f(X): DLMF 7.12.2
        cr = complex(1.0)
        cf = complex(1.0)
        for k = 1:20
            cr = -0.25 * cr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / zp2
            cf += cr
        end

        # g(X): DLMF 7.12.3
        cr = 1.0
        cg = 1.0
        for k = 1:12
            cr = -0.25 * cr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / zp2
            cg += cr
        end

        # DLMF 7.5.4:  Asymptotic Expansions
        cg = cg / (pi * z * z)
        s = d - (cf * cos(zp) + cg * sin(zp)) / (pi * z)
    end

    zf = s
    # S'(z) = sin(pi/2 * z^2)
    zd = sin(0.5 * pi * z * z)
    return zf, zd
end


# TODO: merge cfc/cfs, use one function (CoSF, P631)
"""
Compute Fresnel integrals C(x) and S(x)

Input
x --- Argument of C(x) and S(x)

Output
C --- C(x)
S --- S(x)
"""
function fcs(x::Float64)
    EPS = 1.0e-15

    xa = abs(x)
    px = pi * xa
    t = 0.5 * px * xa
    t2 = t * t

    c = 0.0
    s = 0.0
    if xa == 0.0
        c = 0.0
        s = 0.0
    elseif xa < 2.5
        r = xa
        c = r
        for k = 1:50
            r = -0.5 * r * (4 * k - 3) / k / (2 * k - 1) / (4 * k + 1) * t2
            c += r
            if abs(r) < abs(c) * EPS
                break
            end
        end
        s = xa * t / 3.0
        r = s
        for k = 1:50
            r = -0.5 * r * (4 * k - 1) / k / (2 * k + 1) / (4 * k + 3) * t2
            s += r
            if abs(r) < abs(s) * EPS
                break
            end
        end
    elseif xa < 4.5
        m = trunc(Int64, 42.0 + 1.75 * t)
        su = 0.0
        c = 0.0
        s = 0.0
        f1 = 0.0
        f0 = 1.0e-100
        for k = m:-1:0
            f = (2 * k + 3) * f0 / t - f1
            if mod(k, 2) == 0
                c += f
            else
                s += f
            end
            su += (2 * k + 1) * f * f
            f1 = f0
            f0 = f
        end
        q = sqrt(su)
        c *= xa / q
        s *= xa / q
    else
        r = 1.0
        f = 1.0
        for k = 1:20
            r = -0.25 * r * (4 * k - 1) * (4 * k - 3) / t2
            f += r
        end
        r = 1.0 / (px * xa)
        g = r
        for k = 1:12
            r = -0.25 * r * (4 * k + 1) * (4 * k - 1) / t2
            g += r
        end
        t0 = t - trunc(t / (2.0 * pi)) * 2.0 * pi
        c = 0.5 + (f * sin(t0) - g * cos(t0)) / px
        s = 0.5 - (f * cos(t0) + g * sin(t0)) / px
    end

    if x < 0.0
        c = -c
        s = -s
    end

    return c, s
end

"""
Compute the complex zeros of Fresnel integral C(z)
or S(z) using modified Newton's iteration method

Input
KF  --- Function code
    KF=1 for C(z) or KF=2 for S(z)
NT  --- Total number of zeros

Output
ZO(L) --- L-th zero of C(z) or S(z)

Routines called:
(1) CFC for computing Fresnel integral C(z)
(2) CFS for computing Fresnel integral S(z)
"""
function fcszo!(zo::Vector{Complex{Float64}}, kf::Int, nt::Int)
    EPS = 1.0e-12

    psq = 0.0
    w = 0.0
    for nr = 1:nt
        psq = 0.0
        if kf == 1
            psq = sqrt(4.0 * nr - 1.0)
        end
        if kf == 2
            psq = 2.0 * sqrt(nr)
        end

        px = psq - log(pi * psq) / (pi^2 * psq^3)
        py = log(pi * psq) / pi / psq
        z = Complex{Float64}(px, py)
        if kf == 2
            if nr == 2; z = complex(2.8334, 0.2443); end
            if nr == 3; z = complex(3.4674, 0.2185); end
            if nr == 4; z = complex(4.0025, 0.2008); end
        end

        it = 0
        w0 = w
        while true
            it += 1

            zf, zd = complex(0.0), complex(0.0)
            if kf == 1
                zf, zd = cfc(z)
            end
            if kf == 2
                zf, zd = cfs(z)
            end

            zp = 1.0
            for i in 1:(nr-1)
                zp *= (z - zo[i])
            end

            zfd = zf / zp
            zq = 0.0
            for i in 1:nr-1
                zw = 1.0
                for j in 1:nr-1
                    if j == i
                        continue
                    end
                    zw *= (z - zo[j])
                end
                zq += zw
            end
            zgd = (zd - zq * zfd) / zp
            z -= zfd / zgd
            w0 = w
            w = abs(z)
            
            if (it > 50) || (abs((w - w0) / w) <= EPS)
                break
            end
        end

        zo[nr] = z
    end
end

"""
Compute modified Fresnel integrals F±(x)
and K±(x)

Input
x   --- Argument of F±(x) and K±(x)
KS  --- Sign code
    KS=0 for calculating F+(x) and K+(x)
    KS=1 for calculating F_(x) and K_(x)

Output
FR  --- Re[F±(x)]
FI  --- Im[F±(x)]
FM  --- |F±(x)|
FA  --- Arg[F±(x)]  (Degs.)
GR  --- Re[K±(x)]
GI  --- Im[K±(x)]
GM  --- |K±(x)|
GA  --- Arg[K±(x)]  (Degs.)
"""
function ffk(ks::Int, x::Float64)
    srd = 57.29577951308233
    eps = 1.0e-15
    pp2 = 1.2533141373155
    p2p = 0.7978845608028654

    xa = abs(x)
    x2 = x^2
    x4 = x2 * x2

    if x == 0.0
        fr = 0.5 * sqrt(0.5 * pi)
        fi = (-1)^ks * fr
        fm = sqrt(0.25 * pi)
        fa = (-1)^ks * 45.0
        gr = 0.5
        gi = 0.0
        gm = 0.5
        ga = 0.0
        return fr, fi, fm, fa, gr, gi, gm, ga
    end

    fr, fi = 0.0, 0.0
    fm, fa = 0.0, 0.0
    gr, gi = 0.0, 0.0
    gm, ga = 0.0, 0.0
    if xa <= 2.5
        xr = p2p * xa
        c1 = xr

        for k in 1:50
            xr = -0.5 * xr * (4.0 * k - 3.0) / k / (2.0 * k - 1.0) / (4.0 * k + 1.0) * x4
            c1 += xr
            if abs(xr / c1) < eps
                break
            end
        end

        s1 = p2p * xa * xa * xa / 3.0
        xr = s1

        for k in 1:50
            xr = -0.5 * xr * (4.0 * k - 1.0) / k / (2.0 * k + 1.0) / (4.0 * k + 3.0) * x4
            s1 += xr
            if abs(xr / s1) < eps
                break
            end
        end

        fr = pp2 * (0.5 - c1)
        fi0 = pp2 * (0.5 - s1)
        fi = (-1)^ks * fi0
    elseif xa < 5.5
        m = trunc(Int64, 42 + 1.75 * x2)
        xsu = 0.0
        xc = 0.0
        xs = 0.0
        xf1 = 0.0
        xf0 = 1.0e-100

        for k in m:-1:0
            xf = (2.0 * k + 3.0) * xf0 / x2 - xf1
            if k % 2 == 0
                xc += xf
            else
                xs += xf
            end
            xsu += (2.0 * k + 1.0) * xf * xf
            xf1 = xf0
            xf0 = xf
        end

        xq = sqrt(xsu)
        xw = p2p * xa / xq
        c1 = xc * xw
        s1 = xs * xw
    else
        xr = 1.0
        xf = 1.0

        for k in 1:12
            xr = -0.25 * xr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / x4
            xf += xr
        end

        xr = 1.0 / (2.0 * xa * xa)
        xg = xr

        for k in 1:12
            xr = -0.25 * xr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / x4
            xg += xr
        end

        c1 = 0.5 + (xf * sin(x2) - xg * cos(x2)) / sqrt(2.0 * pi) / xa
        s1 = 0.5 - (xf * cos(x2) + xg * sin(x2)) / sqrt(2.0 * pi) / xa
    end

    fr = pp2 * (0.5 - c1)
    fi0 = pp2 * (0.5 - s1)
    fi = (-1)^ks * fi0
    fm = abs(fr + fi * im)
    
    if fr >= 0.0
        fa = srd * atan(fi / fr)
    elseif fi > 0.0
        fa = srd * (atan(fi / fr) + pi)
    elseif fi < 0.0
        fa = srd * (atan(fi / fr) - pi)
    end

    xp = x2 + pi / 4.0
    cs = cos(xp)
    ss = sin(xp)
    xq2 = 1.0 / sqrt(pi)

    gr = xq2 * (fr * cs + fi0 * ss)
    gi = (-1)^ks * xq2 * (fi0 * cs - fr * ss)
    gm = sqrt(gr * gr + gi * gi)

    if gr >= 0.0
        ga = srd * atan(gi / gr)
    elseif gi > 0.0
        ga = srd * (atan(gi / gr) + pi)
    elseif gi < 0.0
        ga = srd * (atan(gi / gr) - pi)
    end

    if x < 0.0
        fr = pp2 - fr
        fi = (-1)^ks * pp2 - fi
        fm = abs(fr + fi * im)
        fa = srd * atan(fi / fr)
        gr = cos(x2) - gr
        gi = -(-1)^ks * sin(x2) - gi
        gm = sqrt(gr * gr + gi * gi)
        ga = srd * atan(gi / gr)
    end
    
    fr, fi, fm, fa, gr, gi, gm, ga
end
