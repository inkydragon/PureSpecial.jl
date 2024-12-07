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

- OTHPL, 23
- LEGZO, 29
- LAGZO, 33
- HERZO, 37

## 3 Gamma, Beta, and Psi Functions

- [`Specfun.gamma2`](@ref) - Compute gamma function `Г(x)`.
- [`Specfun.lgama`](@ref) - Compute gamma function `ln[Γ(x)]` or `Γ(x)`.
- [`Specfun.cgama`](@ref)
- [`Specfun.beta`](@ref)
- [`Specfun.psi`](@ref)
- [`Specfun.cpsi`](@ref)
- [`Specfun.incog`](@ref)
- [`Specfun.incob`](@ref)

## 4 Legendre Functions

- LPN, 81
- CLPN, 82
- LPNI, 82

- LQNA, 86
- LQNB, 87
- CLQN, 89

- LPMN, 95
- [`Specfun.lpmns`](@ref)
- CLPMN, 96
- LQMN, 102
- [`Specfun.lqmns`](@ref)
- CLQMN, 104
- [`Specfun.lpmv`](@ref)

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
- LAMN(MSTA1,MSTA2), 184
- LAMV(GAMO,MSTA1,MSTA2), 184

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

- [`Specfun.klvna`](@ref)
- KLVNB, 318: SKIP, use `klvna`
- [`Specfun.klvnzo`](@ref)

## 10 Airy Functions

- AIRYA(AJYIK), 329
- [`Specfun.airyb`](@ref)
- [`Specfun.itairy`](@ref)
- [`Specfun.airyzo`](@ref)

## 11 Struve Functions

- STVH1, 347
- STVH0, 347
- STVHV(GAMMA), 349
- [`Specfun.itsh0`](@ref)
- [`Specfun.itth0`](@ref)
- STVL0, 357
- STVL1, 358
- STVLV(GAMMA), 360
- [`Specfun.itsl0`](@ref)

## 12 Hypergeometric and Confluent Hypergeometric Functions

- HYGFX(GAMMA,PSI), 376
- HYGFZ(GAMMA,PSI), 380
- [`Specfun.chgm`](@ref)
- [`Specfun.cchg`](@ref)
- [`Specfun.chgu`](@ref)
  - [`Specfun.chgus`](@ref)
  - [`Specfun.chgul`](@ref)
  - [`Specfun.chgubi`](@ref)
  - [`Specfun.chguit`](@ref)

## 13 Parabolic Cylinder Functions

- [`Specfun.pbdv`](@ref)
  - [`Specfun.dvsa`](@ref)
  - [`Specfun.dvla`](@ref)
- [`Specfun.pbvv`](@ref)
  - [`Specfun.vvsa`](@ref)
  - [`Specfun.vvla`](@ref)
- [`Specfun.pbwa`](@ref)
- [`Specfun.cpbdn`](@ref)
  - [`Specfun.cpdsa`](@ref)
  - [`Specfun.cpdla`](@ref)
  - [`Specfun.gaih`](@ref)

## 14 Mathieu Functions

- CVA1, 501
- CVA2(REFINE,CVO,CVQL,CVQM), 504
  - REFINE(CVF), 505
  - CVF, 506
  - CV0(CVQL,CVQM), 507
  - CVQL, 511
  - CVQM, 511
- MTU0(CVA2,FCOEF), 516
  - FCOEF, 512
- MTU12(CVA2,FCOEF,JYNB), 517

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

- [`Specfun.erf(::Float64)`](@ref)
- [`Specfun.erf(::Complex{Float64})`](@ref)
- [`Specfun.fcs`](@ref)
- [`Specfun.ffk`](@ref)
- [`Specfun.cerzo`](@ref)
  - [`Specfun.cerf`](@ref)
- [`Specfun.fcszo`](@ref)
  - [`Specfun.cfc`](@ref)
  - [`Specfun.cfs`](@ref)

## 17 Cosine and Sine Integrals

- [`Specfun.cisia`](@ref)
- [`Specfun.cisib`](@ref)

## 18 Elliptic Integrals and Jacobian Elliptic Functions

- COMELP, 662
- ELIT, 664
- ELIT3, 665
- JELP, 671

## 19 Exponential Integrals

- [`Specfun.e1xa`](@ref)
- [`Specfun.e1xb`](@ref)
- [`Specfun.e1z`](@ref)
- [`Specfun.enxa`](@ref)
- [`Specfun.enxb`](@ref)
- [`Specfun.eix`](@ref)
