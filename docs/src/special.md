# Special functions
> ref: [scipy.special ## SciPy Manual](https://docs.scipy.org/doc/scipy/reference/special.html)

## Airy functions
- [ ] `airy(z)`: cephes_airy; amos_airy,amos_biry
- [ ] `airye(z)`: amos_airy,amos_biry
- `itairy(x)`: specfun_itairy
- `ai_zeros(nt)`: specfun_airyzo
- `bi_zeros(nt)`: specfun_airyzo


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
- [ ] `it2j0y0`: specfun_ittjya
- [ ] `iti0k0`: specfun_itika
- [ ] `it2i0k0`: specfun_ittika

- [ ] `besselpoly`: cephes_besselpoly

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

- [ ] `fresnel`: specfun_cfs,specfun_cfc
- [ ] `modfresnelp`: specfun_ffk
- [ ] `modfresnelm`: specfun_ffk

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

- [ ] `hyp1f1`: specfun_chgm,specfun_cchg
- [ ] `hyperu`: cephes_poch,specfun_chgu
- [ ] `hyp0f1`: cbesi_wrap,cbesj_wrap,_hyp0f1_asy

## Parabolic Cylinder functions
- [ ] `pbdv`: specfun_pbdv
- [ ] `pbvv`: specfun_pbvv
- [ ] `pbwa`: specfun_pbwa

## Mathieu functions
+ [ ]: specfun_cva2
+ [ ]: specfun_mtu0
+ [ ]: specfun_mtu12
+ [ ]: specfun_fcoef


## Spheroidal Wave functions
> specfun

## Kelvin functions
- `kelvin`: specfun_klvna
- `kelvin_zeros`: specfun_klvnzo

