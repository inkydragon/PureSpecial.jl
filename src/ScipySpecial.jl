# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md


"""Elementary Functions"""
#=
lambertw(z)
    + lambertw_scalar -> special::lambertw
wrightomega(z)
    + wrightomega -> wright::wrightomega
    + wrightomega_real -> wright::wrightomega_real
=#


"""Gamma Functions"""
#=
## cephes
gamma(z)            + cephes_Gamma,cgamma
gammaln(x)          + cephes_lgam
loggamma(z)         + loggamma_real,loggamma -> special::loggamma -> cephes::lgam
gammasgn(x)         + cephes_gammasgn
gammainc(a, x)      + cephes_igam
gammaincinv(a, y)   + cephes_igami
gammaincc(a, x)     + cephes_igamc
gammainccinv(a, y)  + cephes_igamci

beta(a, b)      + cephes_beta
betaln(a, b)    + cephes_lbeta


## boost
betainc(a, b, x)        + ibeta_float,ibeta_double @ boost
betaincc(a, b, x)       + ibetac_float,ibetac_double @ boost
betaincinv(a, b, y)     + ibeta_inv_float,ibeta_inv_double @ boost
betainccinv(a, b, y)    + ibetac_inv_float,ibetac_inv_double @ boost


psi(z)      + digamma,cdigamma -> special::digamma
digamma(z)  + -> psi
rgamma(z)
    + cephes_rgamma
    + crgamma -> special::rgamma
polygamma(n, x)     + gamma,zeta,psi
multigammaln(a, d)  + loggam
poch(z, m)          + cephes_poch
=#


"""Exponential and Trigonometric Integrals"""
#=
exp1(z)
expi(x)
expn(n, x)
=#


"""Error Functions"""
#=
## cephes,faddeeva
erf(z)
    + cephes_erf
    + faddeeva_erf
erfc(x)
    + cephes_erfc
    + faddeeva_erfc_complex
erfcx(x)
    + faddeeva_erfcx
    + faddeeva_erfcx_complex
erfi(z)
    + faddeeva_erfi
    + faddeeva_erfi_complex
erfinv(y)
    + erfinv_float,erfinv_double @ boost
erfcinv(y)
    + cephes_erfcinv
wofz(z)
    + faddeeva_w
dawsn(x)
    + faddeeva_dawsn
    + faddeeva_dawsn_complex
fresnel(z)
    + cephes_fresnl
    + cfresnl_wrap -> specfun_cfs,specfun_cfc
modfresnelp(x)
    + modified_fresnel_plus_wrap -> specfun_ffk(ks=0, ...)
modfresnelm(x)
    + modified_fresnel_minus_wrap -> specfun_ffk(ks=1, ...)
voigt_profile(x, sigma, gamma)
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
    + cbesj_wrap_real -> cbesj_wrap,cephes_jv
    + cbesj_wrap -> amos_besj,amos_besy
        + cbesj_wrap_e
jve(v, z)   + cbesj_wrap_e_real -> cbesj_wrap_e -> amos_besj,amos_besy
yn(n, x)    + cephes_yn
yv(v, z)
    + cbesy_wrap_real -> cbesy_wrap,cephes_yv
    + cbesy_wrap -> amos_besy,amos_besj
yve(v, z)   + cbesy_wrap_e_real -> cbesy_wrap_e -> amos_besy,amos_besj


## amos_besk,amos_besi
kn(n, x)    + cbesk_wrap_real_int -> cbesk_wrap_real -> cbesk_wrap -> amos_besk
kv(v, z)    + cbesk_wrap_real -> cbesk_wrap -> amos_besk
kve(v, z)   + cbesk_wrap_e_real -> cbesk_wrap_e -> amos_besk
iv(v, z)
    + cbesi_wrap -> amos_besi,cbesi_wrap_e,amos_besk
        + cbesi_wrap_e -> amos_besi,amos_besk
    + cephes_iv
ive(v, z)   + cbesi_wrap_e_real -> cbesi_wrap_e -> amos_besi,amos_besk


## amos_besh
hankel1(v, z) + cbesh_wrap1 -> amos_besh
hankel1e(v, z) + cbesh_wrap1_e -> amos_besh
hankel2(v, z) + cbesh_wrap2 -> amos_besh
hankel2e(v, z) + cbesh_wrap2_e -> amos_besh

wright_bessel(a, b, x)  + wright_bessel_scalar@_wright_bessel.pxd
lmbda(v, x)     + _specfun.lamv -> specfun_lamv


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

j0(x)   + cephes_j0
j1(x)   + cephes_j1
y0(x)   + cephes_y0
y1(x)   + cephes_y1
i0(x)   + cephes_i0
i0e(x)  + cephes_i0e
i1(x)   + cephes_i1
i1e(x)  + cephes_i1e
k0(x)   + cephes_k0
k0e(x)  + cephes_k0e
k1(x)   + cephes_k1
k1e(x)  + cephes_k1e


# Integrals of Bessel functions
itj0y0(x)   + it1j0y0_wrap -> specfun_itjya
it2j0y0(x)  + it2j0y0_wrap -> specfun_ittjya
iti0k0(x)   + it1i0k0_wrap -> specfun_itika
it2i0k0(x)  + it2i0k0_wrap -> specfun_ittika

besselpoly(a, lmb, nu)  + cephes_besselpoly


# Derivatives of Bessel functions

jvp(v, z)   + jv(v, z),_bessel_diff_formula
yvp(v, z)   + yv(v, z),_bessel_diff_formula
kvp(v, z)   + kv(v, z),_bessel_diff_formula
ivp(v, z)   + iv(v, z),_bessel_diff_formula
h1vp(v, z)  + hankel1(v, z),_bessel_diff_formula
h2vp(v, z)  + hankel2(v, z),_bessel_diff_formula


# Spherical Bessel functions

spherical_jn(n, z)
+ _spherical_jn,_spherical_jn_d
    + spherical_jn_real,spherical_jn_complex -> cbesj
spherical_yn(n, z)
+ _spherical_yn,_spherical_yn_d
    + spherical_yn_d_real -> spherical_yn_real
    + spherical_yn_d_complex -> spherical_yn_complex -> cbesy
spherical_in(n, z)
+ _spherical_in,_spherical_in_d
    + spherical_in_real -> iv
    + spherical_in_complex -> cbesi_wrap
spherical_kn(n, z)
+ _spherical_kn,_spherical_kn_d
    + spherical_yn_d_real -> spherical_yn_real
    + spherical_yn_d_complex -> spherical_yn_complex -> cbesy


# Riccati-Bessel functions
riccati_jn(n, x) + _specfun.rctj -> specfun_rctj
riccati_yn(n, x) + _specfun.rcty -> specfun_rcty


# Kelvin functions
kelvin(x)
    ber(x)
    bei(x)
    berp(x)
    beip(x)
    ker(x)
    kei(x)
    kerp(x)
    keip(x)

kelvin_zeros(nt)
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
    + cephes_struve_h
modstruve(v, x)
    + cephes_struve_l

itstruve0(x)
it2struve0(x)
itmodstruve0(x)
=#


"""Parabolic cylinder functions"""
#=
pbdv(v, x)
pbvv(v, x)
pbwa(a, x)

pbdv_seq(v, x)
pbvv_seq(v, x)
pbdn_seq(n, z)
=#


"""Hypergeometric functions"""
#=
hyp2f1(a, b, c, z)
    Gauss hypergeometric function 2F1(a, b; c; z)
    + hyp2f1 -> cephes_hyp2f1 @cephes/hyp2f1.c
    + hyp2f1_complex @_hyp2f1.pxd

hyp1f1(a, b, x)
hyperu(a, b, x)

hyp0f1(v, z)
    Confluent hypergeometric limit function 0F1.
    + _hyp0f1_real -> _hyp0f1_asy@cython
    + _hyp0f1_cmplx -> cbesi_wrap @amos_wrappers.c
                    -> cbesj_wrap @amos_wrappers.c
=#


"""Legendre functions"""
#=
sph_harm(m, n, theta, phi)
    + sph_harmonic_unsafe -> sph_harmonic -> cephes_poch

## specfun
lpmv(m, v, x)
clpmn(m, n, z)
lpn(n, z)
lqn(n, z)
lpmn(m, n, z)
lqmn(m, n, z)
=#


"""Elliptic Integrals and Functions"""
#=
## cephes
ellipj(u, m)        + cephes_ellpj
ellipk(m)           + ellipk -> ellpk -> cephes_ellpk
ellipkm1(p)         + cephes_ellpk
ellipkinc(phi, m)   + cephes_ellik
ellipe(m)           + cephes_ellpe
ellipeinc(phi, m)   + cephes_ellie

## ellint_carlson
elliprc(x, y)       + fellint_RC,cellint_RC
elliprd(x, y, z)    + fellint_RD,cellint_RD
elliprf(x, y, z)    + fellint_RF,cellint_RF
elliprg(x, y, z)    + fellint_RG,cellint_RG
elliprj(x, y, z, p) + fellint_RJ,cellint_RJ
=#


"""Zeta Functions"""
#=
zeta(x[, q])    -> xsf::cephes
zetac(x)        -> xsf::cephes
=#


"""Mathieu Functions"""
#=
mathieu_a(m, q)
mathieu_b(m, q)
mathieu_even_coef(m, q)
mathieu_odd_coef(m, q)
mathieu_cem(m, q, x)
mathieu_sem(m, q, x)
mathieu_modcem1(m, q, x)
mathieu_modcem2(m, q, x)
mathieu_modsem1(m, q, x)
mathieu_modsem2(m, q, x)
=#


"""Spheroidal Wave Functions"""
#=
pro_ang1(m, n, c, x)
pro_rad1(m, n, c, x)
pro_rad2(m, n, c, x)
obl_ang1(m, n, c, x)
obl_rad1(m, n, c, x)
obl_rad2(m, n, c, x)

pro_cv(m, n, c)
obl_cv(m, n, c)
pro_cv_seq(m, n, c)
obl_cv_seq(m, n, c)

pro_ang1_cv(m, n, c, cv, x)
pro_rad1_cv(m, n, c, cv, x)
pro_rad2_cv(m, n, c, cv, x)
obl_ang1_cv(m, n, c, cv, x)
obl_rad1_cv(m, n, c, cv, x)
obl_rad2_cv(m, n, c, cv, x)
=#


"""Miscellaneous Functions"""
#=
=#
