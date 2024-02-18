# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Error function and Fresnel integrals

## cephes,faddeeva
erf(z)
    Returns the error function of complex argument.
    + cephes_erf
    + faddeeva_erf

erfc(x)
    Complementary error function, 1 - erf(x).
    + cephes_erfc
    + faddeeva_erfc_complex

erfcx(x)
    Scaled complementary error function, exp(x**2) * erfc(x).
    + faddeeva_erfcx
    + faddeeva_erfcx_complex

erfi(z)
    Imaginary error function, -i erf(i z).
    + faddeeva_erfi
    + faddeeva_erfi_complex

erfinv(y)
    Inverse of the error function.
    + erfinv_float,erfinv_double @ boost

erfcinv(y)
    Inverse of the complementary error function.
    + cephes_erfcinv

wofz(z)
    Faddeeva function
    + faddeeva_w

dawsn(x)
    Dawson's integral.
    + faddeeva_dawsn
    + faddeeva_dawsn_complex

fresnel(z)
    Fresnel integrals.
    + cephes_fresnl
    + cfresnl_wrap -> specfun_cfs,specfun_cfc

modfresnelp(x)
    Modified Fresnel positive integrals
    + modified_fresnel_plus_wrap -> specfun_ffk(ks=0, ...)

modfresnelm(x)
    Modified Fresnel negative integrals
    + modified_fresnel_minus_wrap -> specfun_ffk(ks=1, ...)

voigt_profile(x, sigma, gamma)
    Voigt profile.
    + faddeeva_voigt_profile


## _specfun
erf_zeros(nt)
    Compute the first nt zero in the first quadrant, ordered by absolute value.
    + _specfun.cerzo -> specfun_cerzo

fresnel_zeros(nt)
    Compute nt complex zeros of sine and cosine Fresnel integrals S(z) and C(z).
    + _specfun.fcszo -> specfun_fcszo

    fresnelc_zeros(nt)
        Compute nt complex zeros of cosine Fresnel integral C(z).
        + _specfun.fcszo(1, nt)

    fresnels_zeros(nt)
        Compute nt complex zeros of sine Fresnel integral S(z).
        + _specfun.fcszo(2, nt)
"""
