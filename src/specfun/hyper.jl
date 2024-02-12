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
