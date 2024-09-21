# Special functions

> ref: [scipy.special ## SciPy Manual](https://docs.scipy.org/doc/scipy/reference/special.html)  
> based commit: *2b84173*

## Airy functions

- [ ] `airy(z)`: cephes_airy; amos_airy,amos_biry
- [ ] `airye(z)`: amos_airy,amos_biry
- [ ] `ai_zeros(nt)`: specfun_airyzo
- [ ] `bi_zeros(nt)`: specfun_airyzo
- [ ] `itairy(x)`: specfun_itairy

## Elliptic functions

- [ ] `ellipj(u, m)`: cephes_ellpj
- [ ] `ellipk(m)`: cephes_ellpk
- [ ] `ellipkm1(p)`: cephes_ellpk
- [ ] `ellipkinc(phi, m)`: cephes_ellik
- [ ] `ellipe(m)`: cephes_ellpe
- [ ] `ellipeinc(phi, m)`: cephes_ellie

- [ ] `elliprc(x, y)`: fellint_RC,cellint_RC
- [ ] `elliprd(x, y, z)`: fellint_RD,cellint_RD
- [ ] `elliprf(x, y, z)`: fellint_RF,cellint_RF
- [ ] `elliprg(x, y, z)`: fellint_RG,cellint_RG
- [ ] `elliprj(x, y, z, p)`: fellint_RJ,cellint_RJ


## Bessel functions

- [ ] `lmbda(v, x)`: specfun_lamv

### Zeros of Bessel functions

- [ ] `jnjnp_zeros(nt)`: specfun_jdzo

- [ ] `jnyn_zeros(n, nt)`: specfun_jyzo
- [ ] `jn_zeros(n, nt)`: specfun_jyzo
- [ ] `jnp_zeros(n, nt)`: specfun_jyzo
- [ ] `yn_zeros(n, nt)`: specfun_jyzo
- [ ] `ynp_zeros(n, nt)`: specfun_jyzo

- [ ] `y0_zeros(nt)`: specfun_cyzo
- [ ] `y1_zeros(nt)`: specfun_cyzo
- [ ] `y1p_zeros(nt)`: specfun_cyzo

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

- [ ] `riccati_jn`: specfun_rctj
- [ ] `riccati_yn`: specfun_rcty

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

> Use `cephes`

## Error function

> faddeeva,cephes,boost_special

- [ ] `fresnel`: cephes_fresnl; cfresnl_wrap -> specfun_cfs,specfun_cfc
- [ ] `modfresnelp(x)`: modified_fresnel_plus_wrap -> specfun_ffk
- [ ] `modfresnelm(x)`: modified_fresnel_minus_wrap -> specfun_ffk

- [ ] `erf_zeros`: specfun_cerzo
- [ ] `fresnel_zeros`: specfun_fcszo
    - `fresnelc_zeros`: fresnel_zeros
    - `fresnels_zeros`: fresnel_zeros

## Legendre functions

- [ ] `sph_harm(m, n, theta, phi)`: cephes_poch

- [ ] `lpmv(m, v, x)`: specfun_lpmv
- [ ] `clpmn(m, n, z)`: specfun_clpmn
- [ ] `lpn(n, z)`: specfun_lpn
- [ ] `lqn(n, z)`: specfun_lpn
- [ ] `lpmn(m, n, z)`: specfun_lpmn
- [ ] `lqmn(m, n, z)`: specfun_lqmn

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

## Mathieu functions

- [ ] `mathieu_a(m, q)`: cem_cva_wrap -> specfun_cva2
- [ ] `mathieu_b(m, q)`: sem_cva_wrap -> specfun_cva2

- [ ] `mathieu_even_coef(m, q)`: 
- [ ] `mathieu_odd_coef(m, q)`: 

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
- [ ] `pro_cv_seq(m, n, c)`: 
- [ ] `obl_cv_seq(m, n, c)`: 

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
- [ ] `bernoulli(n)`: 
- [ ] `binom(x, y)`: 
- [ ] `diric(x, n)`: 
- [ ] `euler(n)`: 
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
