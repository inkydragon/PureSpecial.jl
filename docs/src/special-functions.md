# Special Function Packages

SpecialFunctions.jl consists of the following sub-packages:

- [ ] GammaFunctions.jl
- [ ] ExponentialIntegrals.jl
- [ ] ErrorFunctions.jl
- [ ] AiryFunctions.jl
- [x] JuliaMath/Bessels.jl
- (not impl) StruveFunctions.jl
- (not impl) ParabolicCylinders.jl
- [x] JuliaMath/HypergeometricFunctions.jl
- (not impl) LegendreFunctions.jl
- (not impl) QFunctions.jl
- (not impl) OrthogonalPolynomials.jl
- [ ] EllipticFunctions.jl
- (not impl) Bernoulli and Euler Polynomials:
  Maybe move to other package?
- [ ] ZetaFunctions.jl
- [x] JuliaMath/Combinatorics.jl
- (not impl) NumberTheoryFunctions.jl
- (not impl) MathieuFunctions.jl
- (not impl) SpheroidalWaves.jl
- (not impl) HeunFunctions.jl
- (not impl) CoulombFunctions.jl


## GammaFunctions.jl

[DLMF: Chapter 5 Gamma Function](https://dlmf.nist.gov/5)
[DLMF: Chapter 8 Incomplete Gamma and Related Functions](https://dlmf.nist.gov/8)

- Factorial Function:  factorial, Binomial,
- Gamma Function:  gamma, InvGamma, LogGamma, ...;
- Incomplete Gamma Function:  incomplete gamma
- Beta Function:  beta
- Pochhammer Function:  poch

## ExponentialIntegrals.jl

[DLMF: Chapter 6 Exponential, Logarithmic, Sine, and Cosine Integrals](https://dlmf.nist.gov/6)

- Exponential Integral:  expint, ...; LogIntegral
- Trigonometric Integral:  Si, Ci, Shi, Chi

## ErrorFunctions.jl

[DLMF: Chapter 7 Error Functions, Dawson’s and Fresnel Integrals](https://dlmf.nist.gov/7)

- Error Function:     erf, ...; erf_inv, ...; faddeeva;
- Dawson Integral:    dawson
- Fresnel Integral:   C(z), S(z), ...

See also:

- [kiranshila/FresnelIntegrals.jl](https://github.com/kiranshila/FresnelIntegrals.jl),
- [MartinMikkelsen/FewSpecialFunctions.jl](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl)
  for Fresnel Integrals (`FresnelC`, `FresnelS`)

## AiryFunctions.jl

[DLMF: Chapter 9 Airy and Related Functions](https://dlmf.nist.gov/9)

- Airy Function:  Ai, Bi, Ai', Bi'
- Zeros of Airy Function:  Ai_zeros, Bi_zeros
- Integral of Airy Function:  Ai_int, Bi_int
- Scorer Function:  Gi, Hi, Gi', Hi'

## Bessels.jl

[DLMF: Chapter 10 Bessel Functions](https://dlmf.nist.gov/10)

Use [JuliaMath/Bessels.jl: Bessel functions for real arguments and orders](https://github.com/JuliaMath/Bessels.jl/)

- Bessel Function:  Jν(x), Yν(x)
- Modified Bessel Function:  Iν​(x), Kν​(x)
- Hankel Function:  Hkν​(z), H1v(z) , H2v(z)
- Spherical Bessel Function:  jν(x), yν(x), iν(x), kν(x)
- Kelvin Function:  kelvin(x), kelvin_zeros(nt); ber(x), bei(x), ker(x), kei(x)

## StruveFunctions.jl

[DLMF: Chapter 11 Struve and Related Functions](https://dlmf.nist.gov/11)

- Struve Function:  Hν(z), Lν(z); integral of ...
- Lommel Function:  sμν(z), Sμν(z)
- Anger and Weber Function:  Jν(z), Eν(z), Aν(z)

See also:

- (Obsolete? last update 2022) [gwater/Struve.jl](https://github.com/gwater/Struve.jl)

## ParabolicCylinders.jl

[DLMF: Chapter 12 Parabolic Cylinder Functions](https://dlmf.nist.gov/12)

- Dν(z), V(a, z), U(a, z)
- DLMF 12.14:  W(a, z)

See also:

- [MartinMikkelsen/FewSpecialFunctions.jl](https://github.com/MartinMikkelsen/FewSpecialFunctions.jl)

## HypergeometricFunctions.jl

[DLMF: Chapter 13 Confluent Hypergeometric Functions](https://dlmf.nist.gov/13)
[DLMF: Chapter 15 Hypergeometric Function](https://dlmf.nist.gov/15)
[DLMF: Chapter 16 Generalized Hypergeometric Functions and Meijer 𝐺-Function](https://dlmf.nist.gov/16)

Use [JuliaMath/HypergeometricFunctions.jl: A Julia package for calculating hypergeometric functions](https://github.com/JuliaMath/HypergeometricFunctions.jl)

- Hypergeometric Function:  ₂F₁(a,b,c,x)
- Confluent Hypergeometric Function:  ₀F₁(a,z), ₁F₁(a,b,z), U(a,b,x), M(a,b,x)
- Kummer Functions:  `₁F₁(a,b,z)`; `M(a,b,x)` (Olver’s confluent hypergeometric function); `U(a,b,z)`
- Whittaker Functions:  Mκμ(z), Wκμ(z)
- Generalized Hypergeometric Function:  pFq(A, B, z); MeijerG

## LegendreFunctions.jl

[DLMF: Chapter 14 Legendre and Related Functions](https://dlmf.nist.gov/14)

- Legendre functions:  Pn(z), Qn(z)
- Associated Legendre Function:  Pmn(z), Qmn(z)

See also:

- [jishnub/LegendrePolynomials.jl](https://github.com/jishnub/LegendrePolynomials.jl)
- [jmert/AssociatedLegendrePolynomials.jl](https://github.com/jmert/AssociatedLegendrePolynomials.jl)

## QFunctions.jl

[DLMF: Chapter 17 𝑞-Hypergeometric and Related Functions](https://dlmf.nist.gov/17)

- TODO

## OrthogonalPolynomials.jl

[DLMF: Chapter 18 Orthogonal Polynomials](https://dlmf.nist.gov/18)

- Jacobi polynomials
- Gegenbauer polynomials (Ultraspherical)
- Chebyshev polynomials
- Legendre polynomials
- Laguerre polynomials
- Hermite polynomials
- Bessel polynomials
- classical discrete orthogonal polynomials
  - Hahn
  - Krawtchouk
  - Meixner
  - Charlier
- Other Orthogonal Polynomials

See also:

- [jverzani/SpecialPolynomials.jl](https://jverzani.github.io/SpecialPolynomials.jl/dev/#Implemented-polynomial-types)
- [JuliaApproximation/ClassicalOrthogonalPolynomials.jl](https://juliaapproximation.github.io/ClassicalOrthogonalPolynomials.jl/stable/)
- [sciml/PolyChaos.jl](https://docs.sciml.ai/PolyChaos/stable/functions/)

## EllipticFunctions.jl

[DLMF: Chapter 19 Elliptic Integrals](https://dlmf.nist.gov/19)
[DLMF: Chapter 20 Theta Functions](https://dlmf.nist.gov/20)
[DLMF: Chapter 22 Jacobian Elliptic Functions](https://dlmf.nist.gov/22)
[DLMF: Chapter 23 Weierstrass Elliptic and Modular Functions](https://dlmf.nist.gov/23)

- Legendre Integral:  K(m), E(m), Π(u, m); F(Φ, m), E(Φ, m), Π(Φ, u, m)
- Symmetric Integral:  RF(x, y, z), RG(x, y, z), RJ(x, y, z, p), RD(x, y, z), RC(x, y)
- Theta Functions:  θ(n, z, q)
- Jacobian Elliptic Functions:  sn(z, k), cn(z, k), dn(z, k)
- Weierstrass Elliptic Functions:  WeierstrassP(z), WeierstrassZeta(z), WeierstrassSigma(x)

See also:

- [EllipticFunctions.jl documentation · EllipticFunctions](https://stla.github.io/EllipticFunctions.jl/)
- [JacobiElliptic API · JacobiElliptic.jl](https://dominic-chang.com/JacobiElliptic.jl/stable/api/) for Jacobian Elliptic Functions
- (Obsolete, last update 2020) [nolta/Elliptic.jl: Elliptic integral and Jacobi elliptic special functions](https://github.com/nolta/Elliptic.jl)

## Bernoulli and Euler Polynomials

[DLMF: Chapter 24 Bernoulli and Euler Polynomials](https://dlmf.nist.gov/24)

- Bernoulli Numbers and Polynomials:  B(n), Bn(x)
- Euler Numbers and Polynomials:  En, En(x)

## ZetaFunctions.jl

[DLMF: Chapter 25 Zeta and Related Functions](https://dlmf.nist.gov/25)

- Riemann Zeta Function: ζ(s), ζ(s, a)
- Dilogarithm:  `Li₂(z)` (dilogarithm), `Liₛ(z)` (polylogarithm)
- Lerch’s Transcendent:  Φ(z, s, a)
- Dirichlet L-function:  `L(s, χ)`; `η(s)` (Dirichlet eta function)

See also:

- Polylogarithm
  - [Expander/PolyLog.jl: Implementation of polylogarithms in Julia](https://github.com/Expander/PolyLog.jl)
  - [mroughan/Polylogarithms.jl: Polylogarithm function and related special functions and sequences.](https://github.com/mroughan/Polylogarithms.jl)

## Combinatorial Analysis

[DLMF: Chapter 26 Combinatorial Analysis](https://dlmf.nist.gov/26)

Use [JuliaMath/Combinatorics.jl: A combinatorics library for Julia](https://github.com/JuliaMath/Combinatorics.jl)

## NumberTheoryFunctions.jl

[DLMF: Chapter 27 Functions of Number Theory](https://dlmf.nist.gov/27)

- TODO

## MathieuFunctions.jl

[DLMF: Chapter 28 Mathieu Functions and Hill’s Equation](https://dlmf.nist.gov/28)

- Mathieu Functions:  cem(a, q, z), sem(a, q, z), cem‘(a, q, z), sem’(b, q, z)
- Characteristic Value of Mathieu function:  mathieu_a(n, q), mathieu_b(n, q), mathieu_exp(a, q)

See also:

- [BBN-Q/MathieuFunctions.jl: Julia package for Mathieu Functions](https://github.com/BBN-Q/MathieuFunctions.jl)
- (Obsolete, last update 2020) [Integer Order Characteristic Values · Mathieu.jl](https://jebej.github.io/Mathieu.jl/dev/api/values.html)

## SpheroidalWaves.jl

[DLMF: Chapter 30 Spheroidal Wave Functions](https://dlmf.nist.gov/30)

- Angular Spheroidal Wave:  PSmn(z, γ), QSmn(z, γ), PS'mn(z, γ), QS'mn(z, γ)
- Radial Spheroidal Wave Function:  S1mn(z, γ), S2mn(z, γ), S1'mn(z, γ), S2'mn(z, γ)
- Eigenvalues:  λmn(γ)

## HeunFunctions.jl

[DLMF: Chapter 29 Lamé Functions](https://dlmf.nist.gov/29)
[DLMF: Chapter 31 Heun Functions](https://dlmf.nist.gov/31)
[Heun and Related Functions—Wolfram Documentation](https://reference.wolfram.com/language/guide/HeunAndRelatedFunctions.html)

- Lamé Functions:
- Heun Functions:

## CoulombFunctions.jl

[DLMF: Chapter 33 Coulomb Functions](https://dlmf.nist.gov/33)

- F_L(η,ρ), G_L(η,ρ), H⁺L(η,ρ), H⁻L(η,ρ)
- `σL(η)` (Coulomb phase shift)

See also:

- [Coulomb Functions · CoulombFunctions.jl](https://www.tipota.org/CoulombFunctions.jl/dev/coulomb_functions/)


## References

- Special Functions implementaion survey:  [Special Functions impl](https://gist.github.com/inkydragon/4f45468ef1a6e8498ff6ff9175175638)
    Investigated the implementation of special functions in `Mathematica`, the `Julia` ecosystem, and `scipy` (`xsf`).
- [DLMF: Software Index](https://dlmf.nist.gov/software/)
- [Special Functions — Wolfram Documentation](https://reference.wolfram.com/language/guide/SpecialFunctions.html)
- [Special Functions -- from Wolfram MathWorld](https://mathworld.wolfram.com/topics/SpecialFunctions.html)
- [SageMath Functions](https://doc.sagemath.org/html/en/reference/functions/index.html)
- [Special Functions — GSL 2.8 documentation](https://www.gnu.org/software/gsl/doc/html/specfunc.html)
