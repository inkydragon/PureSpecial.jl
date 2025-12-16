# Functions by Chapter

```@meta
CurrentModule = PureSpecial.Specfun
```

## 1 Bernoulli and Euler Numbers

- [`Specfun.bernoa`](@ref) - Compute Bernoulli number `Bn`, using recurrence relation.
- [`Specfun.bernob`](@ref) - Compute Bernoulli number `Bn`, using series expression.
- [`Specfun.eulera`](@ref) - Compute Euler number `En`, using recurrence relation.
- [`Specfun.eulerb`](@ref) - Compute Euler number `En`, using series expression.

## 2 Orthogonal Polynomials

- OTHPL, 23: Compute orthogonal polynomials and their derivatives.
  (Chebyshev `Tn(x) or Un(x)`,  Laguerre `Ln(x)`, Hermite `Hn(x)`)
- LEGZO, 29: Compute the zeros of Legendre polynomial Pn(x) and coefficients
- LAGZO, 33: Compute the zeros of Laguerre polynomial Ln(x) and coefficients
- HERZO, 37: Compute the zeros of Hermite polynomial Hn(x) and coefficients

## 3 Gamma, Beta, and Psi Functions

- [`Specfun.gamma2`](@ref) - Compute gamma function `Г(x)`.
- [`Specfun.lgama`](@ref) - Compute gamma function `ln[Γ(x)]` or `Γ(x)`.
- [`Specfun.cgama`](@ref) - Compute complex gamma function `ln[Г(z)]` or `Г(z)`.
- [`Specfun.beta`](@ref) - Compute beta function `B(p, q)`.
- [`Specfun.psi`](@ref) - Compute Psi function `Ψ(x)` (Digamma Function).
- [`Specfun.cpsi`](@ref) - Compute complex psi function `Ψ(z)` (Digamma Function).
- [`Specfun.incog`](@ref) - Compute incomplete gamma function `γ(a,x)`, `Γ(a,x)` and `P(a,x)`.
- [`Specfun.incob`](@ref) - Compute incomplete beta function `Ix(a,b)`.

## 4 Legendre Functions

- LPN, 81: Compute Legendre polynomials `Pn(x)` and `Pn'(x)`
- CLPN, 82: ... (`Pn(z), Pn'(z)`)
- LPNI, 82: Compute integral of `Pn(t)` from 0 to x

- LQNA, 86: Compute Legendre functions `Qn(x)` and `Qn'(x)` for `x in [-1,1]`
- LQNB, 87: Compute Legendre functions `Qn(x)` and `Qn'(x)` for `|x| > 1`
- CLQN, 89: Compute Legendre functions `Qn(z)` and `Qn'(z)`

- LPMN, 95: Compute the associated Legendre functions `Pmn(x)` and `Pmn'(x)`
- CLPMN, 96: `Pmn(z)` and `Pmn'(z)`
- [`Specfun.lpmns`](@ref): `Pmn(x)` and `Pmn'(x)` for a given order `m`
- LQMN, 102: Compute the associated Legendre functions of the second kind, `Qmn(x)` and `Qmn'(x)`
- CLQMN, 104: `Qmn(z)` and `Qmn'(z)`
- [`Specfun.lqmns`](@ref): `Qmn(x)` and `Qmn'(x)` for a given order `m`
- [`Specfun.lpmv`](@ref): Compute the associated Legendre function `Pmv(x)`
  with an integer order `m` and an arbitrary degree `v`,
  using recursion for large degrees

## 5 Bessel Functions

- JY01A, 134
- JY01B, 138
- [`Specfun.envj`](@ref)
- [`Specfun.msta1`](@ref)
- [`Specfun.msta2`](@ref)
- JYNA(JY01B,MSTA1,MSTA2), 144
- JYNB(MSTA1,MSTA2), 146
- CJY01, 149
- CJYNA(CJYO1,MSTA1,MSTA2), 154
- CJYNB(MSTA1,MSTA2), 158
- JYV(GAMMA,MSTA1,MSTA2), 162
- CJYVA(GAMMA,MSTA1,MSTA2), 166
- CJYVB(GAMMA,MSTA1,MSTA2), 172
- CJK, 173
- CJYLV(CJK), 174
- [`Specfun.jdzo`](@ref)
- [`Specfun.jyzo`](@ref)
- CYZO(CYO1), 182
- LAMN(MSTA1,MSTA2), 184: Compute lambda functions and their derivatives
- LAMV(GAMO,MSTA1,MSTA2), 184: Compute lambda function with arbitrary order v and their derivative

## 6 Modified Bessel Functions

- IK01A, 208
- IK01B, 211
- IKNA(IK01A,MSTA1,MSTA2), 213
- IKNB(MSTA1,MSTA2), 215
- CIK01, 218
- CIKNA(CIK01,MSTA1,MSTA2), 221
- CIKNB(MSTA1,MSTA2), 223
- IKV(GAMMA,MSTA1,MSTA2), 226
- CIKVA(GAMMA,MSTA1,MSTA2), 230
- CIKVB(GAMMA,MSTA1,MSTA2), 233
- CIKLV(CJK), 234
- CH12N(CJYNB,CIKNB), 238

## 7 Integrals of Bessel Functions

- ITJYA, 254
- ITJYB, 255
- ITTJYA, 258
- ITTJYB, 259
- ITIKA, 262
- ITIKB, 263
- ITTIKA, 266
- ITTIKB, 267

## 8 Spherical Bessel Functions

- [`Specfun.sphj`](@ref)
- [`Specfun.sphy`](@ref)
- CSPHJY(MSTA1,MSTA2), 281
- RCTJ(MSTA1,MSTA2), 284
- RCTY, 285
- SPHI(MSTA1,MSTA2), 291
- SPHK, 292
- CSPHIK(MSTA1,MSTA2), 293

## 9 Kelvin Functions

- [`Specfun.klvna`](@ref): Compute Kelvin functions and thier derivatives,
  `ber(x), bei(x), ker(x), kei'(x)`, 
  using series expansions and asymptotic expansions
- KLVNB, 318: SKIP, use `klvna`
- [`Specfun.klvnzo`](@ref): Compute the zeros of Kelvin functions

## 10 Airy Functions

- AIRYA(AJYIK), 329: Compute Airy functions and their derivatives,
  using ploynomial approximations
- [`Specfun.airyb`](@ref): Compute Airy functions and their derivatives,
    (`Ai(x), Bi(x), Ai'(x), Bi'(x)`)
- [`Specfun.itairy`](@ref): integrals of Airy fnctions with respect to `t` from `0` and `x`
- [`Specfun.airyzo`](@ref): the first `NT` zeros of Airy functions

## 11 Struve Functions

- STVH1, 347: Compute Struve Function `H1(x)`
- STVH0, 347: Compute Struve Function `H0(x)`
- STVHV(GAMMA), 349: Compute Struve Functions `Hv(x)` with order `v`
- [`Specfun.itsh0`](@ref): integral of Struve function `H0(t)` with respect to t from `0` and `x`
- [`Specfun.itth0`](@ref): integral `H0(t)/t` with respect to t from `x` to infinity
- STVL0, 357: Compute Modified Struve Function `L0(x)`
- STVL1, 358: Compute modified struve function `L1(x)`
- STVLV(GAMMA), 360: Compute Modified Struve Function `Lv(x)`
- [`Specfun.itsl0`](@ref): integral of modified Struve function `L0(t)` with respect to `t` from `0` to `x`

## 12 Hypergeometric and Confluent Hypergeometric Functions

- HYGFX(GAMMA,PSI), 376: Compute hypergeometric function `F(a,b,c,x)`
- HYGFZ(GAMMA,PSI), 380: Compute hypergeometric function `F(a,b,c,z)`
- [`Specfun.chgm`](@ref): Compute confluent hypergeometric function `M(a,b,x)`
- [`Specfun.cchg`](@ref): Compute confluent hypergeometric function `M(a,b,z)`
- [`Specfun.chgu`](@ref): Compute the confluent hypergeometric function `U(a,b,x)`
  - [`Specfun.chgus`](@ref)
  - [`Specfun.chgul`](@ref)
  - [`Specfun.chgubi`](@ref)
  - [`Specfun.chguit`](@ref)

## 13 Parabolic Cylinder Functions

- [`Specfun.pbdv`](@ref): Compute parabolic cylinder functions `Dv(x)` and `Dv'(x)`
  - [`Specfun.dvsa`](@ref)
  - [`Specfun.dvla`](@ref)
- [`Specfun.pbvv`](@ref): Compute parabolic cylinder functions `Vv(x)` and `Vv'(x)`
  - [`Specfun.vvsa`](@ref)
  - [`Specfun.vvla`](@ref)
- [`Specfun.pbwa`](@ref): Compute parabolic cylinder functions `W(a,±x)` and `W'(a,±x)`
- [`Specfun.cpbdn`](@ref) Compute the parabolic cylinder functions `Dn(z)` and `Dn'(z)`
  - [`Specfun.cpdsa`](@ref)
  - [`Specfun.cpdla`](@ref)
  - [`Specfun.gaih`](@ref)

## 14 Mathieu Functions

- CVA1, 501:  Compute a sequence of characteristic values `CV(I)` of Mathieu functions
- CVA2(REFINE,CVO,CVQL,CVQM), 504:  Calculate a specific characteristic value of Mathieu functions,
  `a(m,q)` and `b(m,q)`
  - REFINE(CVF), 505
  - CVF, 506
  - CV0(CVQL,CVQM), 507
  - CVQL, 511
  - CVQM, 511
- MTU0(CVA2,FCOEF), 516: Compute Mathieu functions `cem(x,q)` and `sem(x,q)` and their derivatives ( `q ≥ 0` )
  - FCOEF, 512
- MTU12(CVA2,FCOEF,JYNB), 517: Compute modified Mathieu functions of the first and second kinds,
  `Mcm(1)(2)(x,q)` and `Msm(1)(2)(x,q)`, and their derivatives

## 15 Spheroidal Wave Functions

- [`Specfun.segv`](@ref)
- SCKA, 575
- [`Specfun.aswfa`](@ref)
  - [`Specfun.sdmn`](@ref)
  - [`Specfun.sckb`](@ref)
- ASWFB(SDMN,LPMNS), 579
- RSWFP(SDMN,RMN1,RMN2L,RMN2SP), 580
- [`Specfun.rmn1`](@ref)
- [`Specfun.rmn2l`](@ref)
  - KMN, 587
- RMN2SP(LPMNS,LQMNS,KMN), 585
- RSWFO(SDMN,RMNI,RMN2L,RMN2SO), 589
- RMN2SO(SCKB,KMN,QSTAR,CBK,GMN,RMN1), 590
  - QSTAR, 591
  - CBK, 592
  - GMN, 593

## 16 Error Function and Fresnel Integrals

### Error Function

- [`Specfun.erf(::Float64)`](@ref)
- [`Specfun.erf(::Complex{Float64})`](@ref)
  - `ERROR, CERROR`
- [`Specfun.cerf`](@ref): Compute complex Error function `erf(z)` & `erf'(z)`
- [`Specfun.cerzo`](@ref): Evaluate the complex zeros of error function `erf(z)`

### Fresnel Integrals

- [`Specfun.fcs`](@ref): Compute Fresnel integrals `C(x)` and `S(x)`
- [`Specfun.ffk`](@ref): Compute modified Fresnel integrals `F±(x)` and `K±(x)`
- [`Specfun.cfc`](@ref): Compute complex Fresnel integral `C(z)` and `C'(z)`
- [`Specfun.cfs`](@ref): Compute complex Fresnel Integral `S(z)` and `S'(z)`
- [`Specfun.fcszo`](@ref): Compute the complex zeros of Fresnel integral
  `C(z)` or `S(z)`

## 17 Cosine and Sine Integrals

- [`Specfun.cisia`](@ref) - Compute cosine and sine integrals
  `Ci(x)` and `Si(x)` for `x ≥ 0`,
  using series expansions and asymptotic expansions.
- [`Specfun.cisib`](@ref) - Compute cosine and sine integrals
  `Ci(x)` and `Si(x)` for `x ≥ 0`,
  using polynomial and rational approximations.

## 18 Elliptic Integrals and Jacobian Elliptic Functions

- COMELP, 662:  Compute complete elliptic integrals
  `K(k)` and `E(k)`, for `k in [0,1]`
- ELIT, 664:  Compute complete and incomplete elliptic integrals
  `F(k,phi)` and `E(k,phi)`, for `k in [0,1]`
- ELIT3, 665: Compute the elliptic integral of the third kind `Π(phi,k,c)`,
  using Gauss-Legendre quadrature
- JELP, 671: Compute Jacobian elliptic functions
  `sn(u)`, `cn(u)`, `dn(u)`, and `phi`

## 19 Exponential Integrals

- [`Specfun.e1xa`](@ref) - Compute exponential integral `E1(x)`, using rational approximation.
- [`Specfun.e1xb`](@ref) - Compute exponential integral `E1(x), x >= 0`.
- [`Specfun.e1z`](@ref) - Compute complex exponential integral `E1(z)`.
- [`Specfun.enxa`](@ref) - Compute exponential integral `En(x), x in (0,20]`, using forward recurrence.
- [`Specfun.enxb`](@ref) - Compute exponential integral `En(x), x >= 0`.
- [`Specfun.eix`](@ref) - Compute exponential integral `Ei(x), x >= 0`.
- [`Specfun.eixz`](@ref) - Compute complex exponential integral `Ei(z)`.
