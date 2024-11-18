# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Exponential and Trigonometric Integrals"""
#=
Exponential Integrals:
- ✅ e1xa
- ✅ e1xb
- ✅ e1z
- ENXA
- ENXB
- ✅ eix
- ✅ eixz

Trigonometric Integrals:
- CISIA, 647
- CISIB, 650
=#

"""
    e1xa(x::Float64)

Compute exponential integral E1(x)

Using rational approximation:
- `0 < x <= 1`:   |eps(x)| <= 2e10-7
- `1 < x`:        |eps(x)| <= 2e10-8

Input:
- `x`  --- Argument of E1(x)

Output:
- E1(x) ( x > 0 )
"""
function e1xa(x::Float64)
    if x == 0.0
        e1 = 1.0e300
    elseif x <= 1.0
        # Use CoSF (19.2.1)
        e1 = -log(x) +
             ((((1.07857e-3 * x - 9.76004e-3) * x + 5.519968e-2) * x -
               0.24991055) * x + 0.99999193) * x - 0.57721566
    else
        # Use CoSF (19.2.7)
        es1 = (((x + 8.5733287401) * x + 18.059016973) * x +
               8.6347608925) * x + 0.2677737343
        es2 = (((x + 9.5733223454) * x + 25.6329561486) * x +
               21.0996530827) * x + 3.9584969228
        e1 = exp(-x) / x * es1 / es2
    end
    return e1
end

"""
    e1xb(x::Float64)

Compute exponential integral E1(x)

Input
- `x`  --- Argument of E1(x)

Output
- E1(x)  ( x > 0 )
"""
function e1xb(x::Float64)
    @assert x >= 0
    EPS = 1e-15

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
        e1 = -SF_EULER_GAMMA_28 - log(x) + x*e1
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


"""
    e1z(z::Complex{Float64})

Compute complex exponential integral E1(z)

Input
- `z`   --- Argument of E1(z)

Output
- E1(z)
"""
function e1z(z::Complex{Float64})
    @assert isapprox(Base.MathConstants.eulergamma, SF_EULER_GAMMA_28)
    EPS = 1e-15

    # Continued fraction converges slowly near negative real axis,
    #   so use power series in a wedge around it until radius 40.0
    a0 = abs(z)
    if a0 == 0.0
        return complex(1e300)
    end

    x = real(z)
    xt = -2.0 * abs(imag(z))
    if (a0 < 5.0) || ((x < xt) && (a0 < 40.0))
        # DLMF 6.6.2:  Power series
        ce1 = complex(1.0)
        cr = complex(1.0)
        for k = 1:500
            cr *= -z * (k / (k+1)^2)
            ce1 += cr
            if abs(cr) < (abs(ce1) * EPS)
                break
            end
        end

        if (x <= 0.0) && (imag(z) == 0.0)
            # Careful on the branch cut -- use the sign of the imaginary part
            #   to get the right sign on the factor if pi.
            ce1 = -SF_EULER_GAMMA_28 - log(-z) + z * ce1 - copysign(pi, imag(z)) * im
        else
            ce1 = -SF_EULER_GAMMA_28 - log(z) + z * ce1
        end
    else
        # DLMF 6.9.1:  Continued Fraction
        #
        #                       1     1     1     2     2     3     3
        #   E1(z) = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
        #                     Z +   1 +   Z +   1 +   Z +   1 +   Z +
        #
        zc = complex(0.0)

        zd = 1 / z
        zdc = zd
        zc += zdc
        for k = 1:500
            zd = 1.0 / (zd * k + 1.0)
            zdc *= (1.0 * zd - 1.0)
            zc += zdc

            zd = 1.0 / (zd * k + z)
            zdc *= (z * zd - 1.0)
            zc += zdc
            if (abs(zdc) <= (abs(zc)*EPS)) && (k > 20)
                break
            end
        end

        ce1 = exp(-z) * zc
        if (x <= 0.0) && (imag(z) == 0.0)
            ce1 -= pi * im
        end
    end

    return ce1
end

"""
Compute exponential integral En(x).

Parameters:
- `n` : Order of En(x)
- `x` : Argument of En(x), ( 0 < x <= 20 )

Returns:
- En(x)

Routine called:
- E1XB for computing E1(x)
"""
function enxa(n::Int, x::Float64)::Vector{Float64}
    @assert (n + 1) >= 2
    @assert 0 < x <= 20
    en = zeros(Float64, n + 1)

    en[1] = exp(-x) / x
    e1 = e1xb(x)
    en[2] = e1
    for k in 2:n
        # Use CoSF (19.1.8)
        ek = (exp(-x) - x * e1) / (k - 1)
        en[k+1] = ek
        e1 = ek
    end

    return en
end


"""
    eix(x::Float64)

Compute exponential integral Ei(x)

Input
- `x`  --- Argument of Ei(x)

Output
- Ei(x)
"""
function eix(x::Float64)
    EPS = 1e-15

    ei = NaN
    if x == 0.0
        # -Inf
        return -1.0e+300
    elseif x < 0
        # DLMF 6.2.6:  Ei(-x) = -E1(x)
        return -e1xb(-x)
    elseif abs(x) <= 40.0
        # DLMF 6.6.1:  x > 0, Power series around x=0
        ei = 1.0
        r = 1.0
        for k = 1:100
            r *= k * x / (k + 1)^2
            ei += r
            if abs(r / ei) <= EPS
                break
            end
        end
        ei = SF_EULER_GAMMA_28 + log(x) + x * ei
    else # x > 40
        # DLMF 6.12.2:  x > 0, x --> Inf, Asymptotic expansion 
        #   (the series is not convergent)
        ei = 1.0
        r = 1.0
        for k = 1:20
            r *= k / x
            ei += r
        end
        ei = exp(x) / x * ei
    end

    return ei
end

"""
    eixz(z::Complex{Float64})

Compute exponential integral Ei(x)

Input
- `x`  --- Complex argument of Ei(x)

Output
- Ei(x)
"""
function eixz(z::Complex{Float64})
    cei = -e1z(-z)

    if imag(z) > 0.0
        cei += 0.0 + pi*im
    elseif imag(z) < 0.0
        cei -= 0.0 + pi*im
    else
        if real(z) > 0.0
            cei += 0.0 + copysign(pi, imag(z))*im
        end
    end

    #= TODO: maybe opt to
        if 0 == imag(z) && real(z) <= 0
            return cei 
        end
        cei += 0.0 + copysign(pi, imag(z))*im
    =#
    return cei
end
