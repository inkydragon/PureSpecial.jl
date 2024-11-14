# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md


"""Elementary Functions"""
#=
Lambert W

lambertw(z)
    + lambertw_scalar -> special::lambertw
        
wrightomega(z)
    + wrightomega -> wright::wrightomega
    + wrightomega_real -> wright::wrightomega_real
=#


"""Gamma Functions"""
#=

## cephes
gamma(z)
    gamma function.
    + cephes_Gamma,cgamma

gammaln(x)
    Logarithm of the absolute value of the gamma function.
    + cephes_lgam

loggamma(z)
    Principal branch of the logarithm of the gamma function.
    + loggamma_real,loggamma -> special::loggamma -> cephes::lgam

gammasgn(x)
    Sign of the gamma function.
    + cephes_gammasgn

gammainc(a, x)
    Regularized lower incomplete gamma function.
    + cephes_igam

gammaincinv(a, y)
    Inverse to the regularized lower incomplete gamma function.
    + cephes_igami

gammaincc(a, x)
    Regularized upper incomplete gamma function.
    + cephes_igamc

gammainccinv(a, y)
    Inverse of the regularized upper incomplete gamma function.
    + cephes_igamci

beta(a, b)
    Beta function.
    + cephes_beta

betaln(a, b)
    Natural logarithm of absolute value of beta function.
    + cephes_lbeta


## boost
betainc(a, b, x)
    Regularized incomplete beta function.
    + ibeta_float,ibeta_double @ boost

betaincc(a, b, x)
    Complement of the regularized incomplete beta function.
    + ibetac_float,ibetac_double @ boost

betaincinv(a, b, y)
    Inverse of the regularized incomplete beta function.
    + ibeta_inv_float,ibeta_inv_double @ boost

betainccinv(a, b, y)
    Inverse of the complemented regularized incomplete beta function.
    + ibetac_inv_float,ibetac_inv_double @ boost


psi(z)
    The digamma function.
    + digamma,cdigamma -> special::digamma

digamma(z)
    The digamma function.
    + -> psi
    
rgamma(z)
    Reciprocal of the gamma function.
    + cephes_rgamma
    + crgamma -> special::rgamma

polygamma(n, x)
    Polygamma functions.
    + gamma,zeta,psi

multigammaln(a, d)
    Returns the log of multivariate gamma, also sometimes called the generalized gamma.
    + loggam

poch(z, m)
    Pochhammer symbol.
    + cephes_poch
=#


"""Exponential and Trigonometric Integrals"""
#=
- exp1(z)
- expi(x)
=#


"""Error Functions"""
#=
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
fresnel_zeros(nt)
    fresnelc_zeros(nt)
    fresnels_zeros(nt)
=#


"""Airy Functions"""
#=
airy(z)
airye(z)
ai_zeros(nt)
bi_zeros(nt)
itairy(x)
=#


"""Bessel functions"""
#=
Bessel functions
jv(v, z)
    Bessel function of the first kind of real order and complex argument.
    + cbesj_wrap_real -> cbesj_wrap,cephes_jv
    + cbesj_wrap -> amos_besj,amos_besy
        + cbesj_wrap_e

jve(v, z)
    Exponentially scaled Bessel function of the first kind of order v.
    + cbesj_wrap_e_real -> cbesj_wrap_e -> amos_besj,amos_besy

yn(n, x)
    Bessel function of the second kind of integer order and real argument.
    + cephes_yn

yv(v, z)
    Bessel function of the second kind of real order and complex argument.
    + cbesy_wrap_real -> cbesy_wrap,cephes_yv
    + cbesy_wrap -> amos_besy,amos_besj

yve(v, z)
    Exponentially scaled Bessel function of the second kind of real order.
    + cbesy_wrap_e_real -> cbesy_wrap_e -> amos_besy,amos_besj


## amos_besk,amos_besi
kn(n, x)
    Modified Bessel function of the second kind of integer order n
    + cbesk_wrap_real_int -> cbesk_wrap_real -> cbesk_wrap -> amos_besk

kv(v, z)
    Modified Bessel function of the second kind of real order v
    + cbesk_wrap_real -> cbesk_wrap -> amos_besk

kve(v, z)
    Exponentially scaled modified Bessel function of the second kind.
    + cbesk_wrap_e_real -> cbesk_wrap_e -> amos_besk

iv(v, z)
    Modified Bessel function of the first kind of real order.
    + cbesi_wrap -> amos_besi,cbesi_wrap_e,amos_besk
        + cbesi_wrap_e -> amos_besi,amos_besk
    + cephes_iv

ive(v, z)
    Exponentially scaled modified Bessel function of the first kind.
    + cbesi_wrap_e_real -> cbesi_wrap_e -> amos_besi,amos_besk


## amos_besh
hankel1(v, z)
    Hankel function of the first kind
    + cbesh_wrap1 -> amos_besh

hankel1e(v, z)
    Exponentially scaled Hankel function of the first kind
    + cbesh_wrap1_e -> amos_besh

hankel2(v, z)
    Hankel function of the second kind
    + cbesh_wrap2 -> amos_besh

hankel2e(v, z)
    Exponentially scaled Hankel function of the second kind
    + cbesh_wrap2_e -> amos_besh

    
wright_bessel(a, b, x)
    Wright's generalized Bessel function.
    + wright_bessel_scalar@_wright_bessel.pxd


lmbda(v, x)
    Jahnke-Emden Lambda function, Lambdav(x).
    + _specfun.lamv -> specfun_lamv



# Zeros of Bessel functions

## jdzo
jnjnp_zeros(nt)

## jyzo
jnyn_zeros(n, nt)
    jn_zeros(n, nt)
    jnp_zeros(n, nt)
    yn_zeros(n, nt)
    ynp_zeros(n, nt)

## cyzo
y0_zeros(nt)
y1_zeros(nt)
y1p_zeros(nt)


# Faster Common Bessel functions

j0(x)
    Bessel function of the first kind of order 0.
    + cephes_j0

j1(x)
    Bessel function of the first kind of order 1.
    + cephes_j1

y0(x)
    Bessel function of the second kind of order 0.
    + cephes_y0

y1(x)
    Bessel function of the second kind of order 1.
    + cephes_y1

i0(x)
    Modified Bessel function of order 0.
    + cephes_i0

i0e(x)
    Exponentially scaled modified Bessel function of order 0.
    + cephes_i0e

i1(x)
    Modified Bessel function of order 1.
    + cephes_i1

i1e(x)
    Exponentially scaled modified Bessel function of order 1.
    + cephes_i1e

k0(x)
    Modified Bessel function of the second kind of order 0, K0(x)
    + cephes_k0

k0e(x)
    Exponentially scaled modified Bessel function K of order 0
    + cephes_k0e

k1(x)
    Modified Bessel function of the second kind of order 1, K1(x)
    + cephes_k1
    
k1e(x)
    Exponentially scaled modified Bessel function K of order 1
    + cephes_k1e


# Integrals of Bessel functions
itj0y0(x)
    Integrals of Bessel functions of the first kind of order 0.
    + it1j0y0_wrap -> specfun_itjya

it2j0y0(x)
    Integrals related to Bessel functions of the first kind of order 0.
    + it2j0y0_wrap -> specfun_ittjya

iti0k0(x)
    Integrals of modified Bessel functions of order 0.
    + it1i0k0_wrap -> specfun_itika

it2i0k0(x)
    Integrals related to modified Bessel functions of order 0.
    + it2i0k0_wrap -> specfun_ittika

besselpoly(a, lmb, nu)
    Weighted integral of the Bessel function of the first kind.
    + cephes_besselpoly


# Derivatives of Bessel functions

jvp(v, z)
Compute derivatives of Bessel functions of the first kind.
+ jv(v, z),_bessel_diff_formula

yvp(v, z)
Compute derivatives of Bessel functions of the second kind.
+ yv(v, z),_bessel_diff_formula

kvp(v, z)
Compute derivatives of real-order modified Bessel function Kv(z)
+ kv(v, z),_bessel_diff_formula

ivp(v, z)
Compute derivatives of modified Bessel functions of the first kind.
+ iv(v, z),_bessel_diff_formula

h1vp(v, z)
Compute derivatives of Hankel function H1v(z) with respect to z.
+ hankel1(v, z),_bessel_diff_formula

h2vp(v, z)
Compute derivatives of Hankel function H2v(z) with respect to z.
+ hankel2(v, z),_bessel_diff_formula


# Spherical Bessel functions

spherical_jn(n, z)
Spherical Bessel function of the first kind or its derivative.
+ _spherical_jn,_spherical_jn_d
    + spherical_jn_real,spherical_jn_complex -> cbesj

spherical_yn(n, z)
Spherical Bessel function of the second kind or its derivative.
+ _spherical_yn,_spherical_yn_d
    + spherical_yn_d_real -> spherical_yn_real
    + spherical_yn_d_complex -> spherical_yn_complex -> cbesy
    
spherical_in(n, z)
Modified spherical Bessel function of the first kind or its derivative.
+ _spherical_in,_spherical_in_d
    + spherical_in_real -> iv
    + spherical_in_complex -> cbesi_wrap

spherical_kn(n, z)
Modified spherical Bessel function of the second kind or its derivative.
+ _spherical_kn,_spherical_kn_d
    + spherical_yn_d_real -> spherical_yn_real
    + spherical_yn_d_complex -> spherical_yn_complex -> cbesy


# Riccati-Bessel functions
riccati_jn(n, x)
Compute Ricatti-Bessel function of the first kind and its derivative.
+ _specfun.rctj -> specfun_rctj

riccati_yn(n, x)
Compute Ricatti-Bessel function of the second kind and its derivative.
+ _specfun.rcty -> specfun_rcty


# Kelvin functions
kelvin(x)
    + kelvin_wrap -> specfun_klvna
        + ret (ber, bei, ger, gei, der, dei, her, hei)
    
    ber(x)
    bei(x)
    berp(x)
    beip(x)
    ker(x)
    kei(x)
    kerp(x)
    keip(x)

kelvin_zeros(nt)
    + _specfun.klvnzo -> specfun_klvnzo
        + ret (ber, bei, ger, gei, der, dei, her, hei)

    ber_zeros(nt)
    bei_zeros(nt)
    berp_zeros(nt)
    beip_zeros(nt)
    ker_zeros(nt)
    kei_zeros(nt)
    kerp_zeros(nt)
    keip_zeros(nt)
=#


"""Struve Functions"""
#=
struve(v, x)
Struve function.
+ cephes_struve_h
    
modstruve(v, x)
Modified Struve function.
+ cephes_struve_l
    
itstruve0(x)
Integral of the Struve function of order 0.
+ itstruve0_wrap -> specfun_itsh0
    
it2struve0(x)
Integral related to the Struve function of order 0.
+ it2struve0_wrap -> specfun_itth0
    
itmodstruve0(x)
Integral of the modified Struve function of order 0.
+ itmodstruve0_wrap -> specfun_itsl0
=#


"""Parabolic cylinder functions"""
#=
pbdv(v, x)
Parabolic cylinder function D
+ pbdv_wrap -> specfun_pbdv

pbvv(v, x)
Parabolic cylinder function V
+ pbvv_wrap -> specfun_pbvv

pbwa(a, x)
Parabolic cylinder function W.
+ pbwa_wrap -> specfun_pbwa

pbdv_seq(v, x)
Parabolic cylinder functions Dv(x) and derivatives.
+ _specfun.pbdv -> specfun_pbdv

pbvv_seq(v, x)
Parabolic cylinder functions Vv(x) and derivatives.
+ _specfun.pbvv -> specfun_pbvv

pbdn_seq(n, z)
Parabolic cylinder functions Dn(z) and derivatives.
+ _specfun.cpbdn -> specfun_cpbdn
=#


"""Hypergeometric functions"""
#=
hyp2f1(a, b, c, z)
    Gauss hypergeometric function 2F1(a, b; c; z)
    + hyp2f1 -> cephes_hyp2f1 @cephes/hyp2f1.c
    + hyp2f1_complex @_hyp2f1.pxd

hyp1f1(a, b, x)
    Confluent hypergeometric function 1F1.
    + hyp1f1_double -> hyp1f1_wrap -> <hyp1f1_wrap> -> boost::math::hypergeometric_1F1
    +               -> hyp1f1_wrap -> specfun_chgm
    + chyp1f1_wrap -> specfun_cchg

hyperu(a, b, x)
    Confluent hypergeometric function U
    + hyperu -> poch -> cephes_poch
             -> hypU_wrap -> specfun_chgu

hyp0f1(v, z)
    Confluent hypergeometric limit function 0F1.
    + _hyp0f1_real -> _hyp0f1_asy@cython
    + _hyp0f1_cmplx -> cbesi_wrap @amos_wrappers.c
                    -> cbesj_wrap @amos_wrappers.c
=#


"""Legendre functions"""
#=
sph_harm(m, n, theta, phi)
Compute spherical harmonics.
+ sph_harmonic_unsafe -> sph_harmonic -> cephes_poch


## specfun
lpmv(m, v, x)
    Associated Legendre function of integer order and real degree.
    + pmv_wrap -> specfun_lpmv

clpmn(m, n, z)
    Associated Legendre function of the first kind for complex arguments.
    + specfun_clpmn

lpn(n, z)
    Legendre function of the first kind.
    + specfun_lpn

lqn(n, z)
    Legendre function of the second kind.
    + specfun_lpn

lpmn(m, n, z)
    Sequence of associated Legendre functions of the first kind.
    + specfun_lpmn

lqmn(m, n, z)
    Sequence of associated Legendre functions of the second kind.
    + specfun_lqmn
=#


"""Elliptic Integrals and Functions"""
#=
## cephes
ellipj(u, m)
    Jacobian elliptic functions
    + cephes_ellpj

ellipk(m)
    Complete elliptic integral of the first kind.
    + ellipk -> ellpk -> cephes_ellpk

ellipkm1(p)
    Complete elliptic integral of the first kind around m = 1
    + cephes_ellpk

ellipkinc(phi, m)
    Incomplete elliptic integral of the first kind
    + cephes_ellik

ellipe(m)
    Complete elliptic integral of the second kind
    + cephes_ellpe

ellipeinc(phi, m)
    Incomplete elliptic integral of the second kind
    + cephes_ellie

    
## ellint_carlson
elliprc(x, y)
    Degenerate symmetric elliptic integral.
    + fellint_RC,cellint_RC

elliprd(x, y, z)
    Symmetric elliptic integral of the second kind.
    + fellint_RD,cellint_RD

elliprf(x, y, z)
    Completely-symmetric elliptic integral of the first kind.
    + fellint_RF,cellint_RF

elliprg(x, y, z)
    Completely-symmetric elliptic integral of the second kind.
    + fellint_RG,cellint_RG

elliprj(x, y, z, p)
    Symmetric elliptic integral of the third kind.
    + fellint_RJ,cellint_RJ
=#


"""Zeta Functions"""
#=

=#


"""Mathieu Functions"""
#=
mathieu_a(m, q)
    Characteristic value of even Mathieu functions
    + cem_cva_wrap -> specfun_cva2

mathieu_b(m, q)
    Characteristic value of odd Mathieu functions
    + sem_cva_wrap -> specfun_cva2

mathieu_even_coef(m, q)
    Fourier coefficients for even Mathieu and modified Mathieu functions.
    + _specfun.fcoef(kd=1,2) -> specfun_fcoef

mathieu_odd_coef(m, q)
    Fourier coefficients for even Mathieu and modified Mathieu functions.
    + _specfun.fcoef(kd=3,4) -> specfun_fcoef

mathieu_cem(m, q, x)
    Even Mathieu function and its derivative
    + cem_wrap -> specfun_mtu0

mathieu_sem(m, q, x)
    Odd Mathieu function and its derivative
    + sem_wrap -> specfun_mtu0

mathieu_modcem1(m, q, x)
    Even modified Mathieu function of the first kind and its derivative
    + mcm1_wrap -> specfun_mtu12

mathieu_modcem2(m, q, x)
    Even modified Mathieu function of the second kind and its derivative
    + mcm2_wrap -> specfun_mtu12

mathieu_modsem1(m, q, x)
    Odd modified Mathieu function of the first kind and its derivative
    + msm1_wrap -> specfun_mtu12

mathieu_modsem2(m, q, x)
    Odd modified Mathieu function of the second kind and its derivative
    + msm2_wrap -> specfun_mtu12
=#


"""Spheroidal Wave Functions"""
#=
pro_ang1(m, n, c, x[, out])
Prolate spheroidal angular function of the first kind and its derivative
+ prolate_aswfa_nocv_wrap -> specfun_segv,specfun_aswfa

pro_rad1(m, n, c, x[, out])
Prolate spheroidal radial function of the first kind and its derivative
+ prolate_radial1_nocv_wrap -> specfun_segv,specfun_rswfp

pro_rad2(m, n, c, x[, out])
Prolate spheroidal radial function of the second kind and its derivative
+ prolate_radial2_nocv_wrap -> specfun_segv,specfun_rswfp

obl_ang1(m, n, c, x[, out])
Oblate spheroidal angular function of the first kind and its derivative
+ oblate_aswfa_nocv_wrap -> specfun_segv,specfun_aswfa

obl_rad1(m, n, c, x[, out])
Oblate spheroidal radial function of the first kind and its derivative
+ oblate_radial1_nocv_wrap -> specfun_segv,specfun_rswfo

obl_rad2(m, n, c, x[, out])
Oblate spheroidal radial function of the second kind and its derivative.
+ oblate_radial2_nocv_wrap -> specfun_segv,specfun_rswfo

pro_cv(m, n, c[, out])
Characteristic value of prolate spheroidal function
+ prolate_segv_wrap -> specfun_segv

obl_cv(m, n, c[, out])
Characteristic value of oblate spheroidal function
+ oblate_segv_wrap -> specfun_segv


pro_cv_seq(m, n, c)
Characteristic values for prolate spheroidal wave functions.
+ _specfun.segv

obl_cv_seq(m, n, c)
Characteristic values for oblate spheroidal wave functions.
+ _specfun.segv


pro_ang1_cv(m, n, c, cv, x[, out])
Prolate spheroidal angular function pro_ang1 for precomputed characteristic value
+ prolate_aswfa_wrap -> specfun_aswfa

pro_rad1_cv(m, n, c, cv, x[, out])
Prolate spheroidal radial function pro_rad1 for precomputed characteristic value
+ prolate_radial1_wrap -> specfun_rswfp

pro_rad2_cv(m, n, c, cv, x[, out])
Prolate spheroidal radial function pro_rad2 for precomputed characteristic value
+ prolate_radial2_wrap -> specfun_rswfp

obl_ang1_cv(m, n, c, cv, x[, out])
Oblate spheroidal angular function obl_ang1 for precomputed characteristic value
+ oblate_aswfa_wrap -> specfun_aswfa

obl_rad1_cv(m, n, c, cv, x[, out])
Oblate spheroidal radial function obl_rad1 for precomputed characteristic value
+ oblate_radial1_wrap -> specfun_rswfo

obl_rad2_cv(m, n, c, cv, x[, out])
Oblate spheroidal radial function obl_rad2 for precomputed characteristic
+ oblate_radial2_wrap -> specfun_rswfo
=#


"""Miscellaneous Functions"""
#=
=#
