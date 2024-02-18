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