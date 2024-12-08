# scipy.special

> ref: [scipy.special ## SciPy Manual](https://docs.scipy.org/doc/scipy/reference/special.html)  
> based commit: *846ecd4648* (2024-12-08)

- `cephes.h`: 100
- `specfun_wrappers.h`: 49
- `amos_wrappers.h`: 14
- `faddeeva.h++`: 9
- `boost_special_functions.h++`: 7

## Airy functions

- [ ] `airy(z)`: cephes_airy; amos_airy,amos_biry
- [ ] `airye(z)`: amos_airy,amos_biry
- [ ] `ai_zeros(nt)`: _specfun.airyzo -> specfun_airyzo
- [ ] `bi_zeros(nt)`: specfun_airyzo
- [ ] `itairy(x)`: specfun_itairy

## Elliptic functions

- [ ] `ellipk(m)`: cephes_ellpk
- [ ] `ellipkm1(p)`: cephes_ellpk
- [ ] `ellipkinc(phi, m)`: cephes_ellik
- [ ] `ellipe(m)`: cephes_ellpe
- [ ] `ellipeinc(phi, m)`: cephes_ellie

- [ ] `elliprf(x, y, z)`: fellint_RF,cellint_RF
- [ ] `elliprg(x, y, z)`: fellint_RG,cellint_RG
- [ ] `elliprj(x, y, z, p)`: fellint_RJ,cellint_RJ
- [ ] `elliprd(x, y, z)`: fellint_RD,cellint_RD
- [ ] `elliprc(x, y)`: fellint_RC,cellint_RC

- [ ] `ellipj(u, m)`: cephes_ellpj
  Jacobian elliptic functions

## Bessel functions

- [ ] `lmbda(v, x)`: _specfun.lamv -> specfun_lamv; _specfun.lamn -> specfun_lamn

### Zeros of Bessel functions

- [ ] `jnjnp_zeros(nt)`: _specfun.jdzo -> specfun_jdzo
- [ ] `jnyn_zeros(n, nt)`: _specfun.jyzo -> specfun_jyzo
  - `jn_zeros(n, nt)`
  - `jnp_zeros(n, nt)`
  - `yn_zeros(n, nt)`
  - `ynp_zeros(n, nt)`

- [ ] _specfun.cyzo -> specfun_cyzo:
  - `y0_zeros(nt)`
  - `y1_zeros(nt)`
  - `y1p_zeros(nt)`

### Faster Common Bessel functions

> Use `cephes`

### Integrals of Bessel functions

- [ ] `itj0y0(x)`: specfun_itjya
- [ ] `it2j0y0(x)`: specfun_ittjya
- [ ] `iti0k0(x)`: specfun_itika
- [ ] `it2i0k0(x)`: specfun_ittika
- [ ] `besselpoly(a, lmb, nu)`: cephes_besselpoly

### Derivatives of Bessel functions

> No 3rd deps.

### Spherical Bessel functions

> AMOS or pure cython.

### Riccati-Bessel functions

- [ ] `riccati_jn`: _specfun.rctj -> specfun_rctj
- [ ] `riccati_yn`: _specfun.rcty -> specfun_rcty

## Struve functions

> Use `cephes`, `specfun`

- [ ] `struve`: cephes_struve_h
- [ ] `modstruve`: cephes_struve_l

- [ ] `itstruve0`: specfun_itsh0
- [ ] `it2struve0`: specfun_itth0
- [ ] `itmodstruve0`: specfun_itsl0

## Statistical functions

> Skip!

## Information Theory functions

> Skip!


## Gamma functions
> Use `cephes`, `boost::math`

### Gamma function
- `gamma(z)`:           `cephes::Gamma(x)`, `xsf::gamma(z)`
- `gammaln(x)`:         xsf::gammaln  -> `cephes::lgam(x)`
- `loggamma(z)`:        xsf::loggamma -> `cephes::lgam(x)`
- `gammasgn(x)`:        `cephes::gammasgn(x)`
- `gammainc(a, x)`:     `cephes::igam(a, x)`
- `gammaincinv(a, y)`:  `cephes::igami(a, p)`
- `gammaincc(a, x)`:    `cephes::igamc(a, x)`
- `gammainccinv(a, y)`: `cephes::igamci(a, p)`
- `rgamma(z)`:          `cephes::rgamma(z)`
- `digamma(z)`:         `cephes::psi(z)`
- `psi(z)`:             `cephes::psi(z)`
- `poch(z, m)`:         `cephes::poch(x, m)`

In Python:
- `polygamma(n, x)`:    in python, `psi, gamma, zeta`
- `multigammaln(a, d)`: in python, `gammaln`

### Beta function
- `beta(a, b)`:     `cephes::beta(a, b)`    
- `betaln(a, b)`:   `cephes::lbeta(a, b)`

Boost::math:
- `betainc(a, b, x)`:     `boost::math::ibeta(a, b, x)`
- `betaincinv(a, b, y)`:  `boost::math::ibeta_inv(a, b, p)`
- `betaincc(a, b, x)`:    `boost::math::ibetac(a, b, x)`
- `betainccinv(a, b, y)`: `boost::math::ibetac_inv(a, b, p)`


## Error functions
> Use: faddeeva; cephes; boost::math; specfun

- `erf(z)`:   `xsf::cephes::erf(x)`;  `Faddeeva::erf(z)`
- `erfc(z)`:  `xsf::cephes::erfc(x)`; `Faddeeva::erfc(z)`
- `erfcx(z)`: `Faddeeva::erfcx(x)`;   `Faddeeva::erfcx(z)`
- `erfi(z)`:  `Faddeeva::erfi(x)`;    `Faddeeva::erfi(z)`
- `erfinv(y)`:  `boost::math::erf_inv(x)`
- `erfcinv(y)`: `xsf::cephes::erfcinv(y)`

special:
- `wofz(z)`:        `Faddeeva::w(z)`
- `dawsn(z)`:       `Faddeeva::Dawson(x)`; `Faddeeva::Dawson(z)`
- `fresnel(z)`:     `cephes::fresnl(x)`; `detail::cfs(z), detail::cfc(z)` (specfun_cfs,specfun_cfc)
- `modfresnelp(x)`: xsf::modified_fresnel_plus  -> `detail::ffk(0)` (specfun_ffk)
- `modfresnelm(x)`: xsf::modified_fresnel_minus -> `detail::ffk(1)` (specfun_ffk)
- `voigt_profile(x, sigma, gamma)`: faddeeva_voigt_profile -> `Faddeeva::w(z)`

zeros:
- `erf_zeros(nt)`:      _specfun.cerzo -> `xsf::specfun::cerzo`
- `fresnel_zeros(nt)`:  _specfun.fcszo(k, nt) -> `xsf::fcszo`
- `fresnelc_zeros(nt)`: _specfun.fcszo(k, nt) -> `xsf::fcszo`
- `fresnels_zeros(nt)`: _specfun.fcszo(k, nt) -> `xsf::fcszo`


## Legendre functions

- [ ] `sph_harm(m, n, theta, phi)`: cephes_poch

- [ ] `lpmv(m, v, x)`: specfun_lpmv
- [ ] `clpmn(m, n, z)`: _specfun.clpmn -> specfun_clpmn
- [ ] `lpn(n, z)`: _specfun.lpn -> specfun_lpn; _specfun.clpn -> specfun_clpn
- [ ] `lqn(n, z)`: _specfun.lqnb -> specfun_lqnb ; _specfun.clqn -> specfun_clqn
- [ ] `lpmn(m, n, z)`: _specfun.lpmn -> specfun_lpmn
- [ ] `lqmn(m, n, z)`: _specfun.lqmn -> specfun_lqmn; _specfun.clqmn -> specfun_clqmn

## Ellipsoidal harmonics

## Orthogonal polynomials

+ evaluate values
+ roots

## Hypergeometric functions

- [ ] `hyp2f1`: cephes_hyp2f1，hyp2f1_complex

- [ ] `hyp1f1`: boost::hyp1f1_double; chyp1f1_wrap -> specfun_cchg
- [ ] `hyperu`: cephes_poch,specfun_chgu
- [ ] `hyp0f1`: cbesi_wrap,cbesj_wrap,_hyp0f1_asy

## Parabolic Cylinder functions

- [ ] `pbdv`: specfun_pbdv
- [ ] `pbvv`: specfun_pbvv
- [ ] `pbwa`: specfun_pbwa

- [ ] `pbdv_seq`: _specfun.pbdv -> specfun_pbdv
- [ ] `pbvv_seq`: _specfun.pbvv -> specfun_pbvv
- [ ] `pbwa_seq`: _specfun.cpbdn -> specfun_cpbdn

## Mathieu functions

- [ ] `mathieu_a(m, q)`: cem_cva_wrap -> specfun_cva2
- [ ] `mathieu_b(m, q)`: sem_cva_wrap -> specfun_cva2

- [ ] _specfun.fcoef -> specfun_fcoef:
  - `mathieu_even_coef(m, q)`
  - `mathieu_odd_coef(m, q)`

- [ ] `mathieu_cem(m, q, x)`: specfun_mtu0
- [ ] `mathieu_sem(m, q, x)`: specfun_mtu0
- [ ] `mathieu_modcem1(m, q, x)`: mcm1_wrap -> specfun_mtu12
- [ ] `mathieu_modcem2(m, q, x)`: mcm2_wrap -> specfun_mtu12
- [ ] `mathieu_modsem1(m, q, x)`: msm1_wrap -> specfun_mtu12
- [ ] `mathieu_modsem2(m, q, x)`: msm2_wrap -> specfun_mtu12

## Spheroidal Wave functions

- [ ] `pro_ang1(m, n, c, x)`: prolate_aswfa_nocv_wrap -> specfun_segv,specfun_aswfa
- [ ] `pro_rad1(m, n, c, x)`: prolate_radial1_nocv_wrap -> specfun_segv(1),specfun_rswfp(1)
- [ ] `pro_rad2(m, n, c, x)`: prolate_radial2_nocv_wrap -> specfun_segv(1),specfun_rswfp(2)
- [ ] `obl_ang1(m, n, c, x)`: oblate_aswfa_nocv_wrap -> specfun_segv,specfun_aswfa
- [ ] `obl_rad1(m, n, c, x)`: oblate_radial1_nocv_wrap -> specfun_segv,specfun_rswfo
- [ ] `obl_rad2(m, n, c, x)`: oblate_radial2_nocv_wrap -> specfun_segv,specfun_rswfo
- [ ] `pro_cv(m, n, c)`: prolate_segv_wrap -> specfun_segv
- [ ] `obl_cv(m, n, c)`: oblate_segv_wrap -> specfun_segv
- [ ] `pro_cv_seq(m, n, c)`: _specfun.segv -> specfun_segv
- [ ] `obl_cv_seq(m, n, c)`: _specfun.segv -> specfun_segv

- [ ] `pro_ang1_cv(m, n, c, cv, x)`: prolate_aswfa_wrap -> specfun_aswfa
- [ ] `pro_rad1_cv(m, n, c, cv, x)`: prolate_radial1_wrap -> specfun_rswfp(1)
- [ ] `pro_rad2_cv(m, n, c, cv, x)`: prolate_radial2_wrap -> specfun_rswfp(2)
- [ ] `obl_ang1_cv(m, n, c, cv, x)`: oblate_aswfa_wrap -> specfun_aswfa
- [ ] `obl_rad1_cv(m, n, c, cv, x)`: oblate_radial1_wrap -> specfun_rswfo
- [ ] `obl_rad2_cv(m, n, c, cv, x)`: oblate_radial2_wrap -> specfun_rswfo

## Kelvin functions

- [ ] `kelvin(x)`: kelvin_wrap -> specfun_klvna
  - `ber(x)`
  - `bei(x)`
  - `berp(x)`
  - `beip(x)`
  - `ker(x)`
  - `kei(x)`
  - `kerp(x)`
  - `keip(x)`
- [ ] `kelvin_zeros(nt)`: specfun_klvnzo

## Combinatorics

## Lambert W functions

> in cpp

## Other Special functions

- [ ] `agm(a, b)`: 
- [ ] `bernoulli(n)`: _specfun.bernob -> specfun_bernob
- [ ] `binom(x, y)`: 
- [ ] `diric(x, n)`: 
- [ ] `euler(n)`: _specfun.eulerb -> specfun_eulerb
- [ ] `expn(n, x)`: 
- [ ] `exp1(z)`: exp1_wrap -> specfun_e1xb; cexp1_wrap -> specfun_e1z
- [ ] `expi(x)`: expi_wrap -> specfun_eix; cexpi_wrap -> specfun_eixz
- [ ] `factorial(n)`: 
- [ ] `factorial2(n)`: 
- [ ] `factorialk(n, k)`: 
- [ ] `shichi(x)`: 
- [ ] `sici(x)`: 
- [ ] `softmax(x)`: 
- [ ] `log_softmax(x)`: 
- [ ] `spence(z)`: 
- [ ] `zeta(x)`: 
- [ ] `zetac(x)`: 

## Convenience functions
