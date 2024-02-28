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
