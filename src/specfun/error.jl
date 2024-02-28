# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

# TODO: merge cfc/cfs, use one function (CoSF, P631)

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
