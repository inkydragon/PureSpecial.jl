# Special Functions

- [DLMF: NIST Digital Library of Mathematical Functions](https://dlmf.nist.gov/)
- [Special Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/SpecialFunctions.html)


## Gamma Functions

> - [DLMF: Chapter 5 Gamma Function](https://dlmf.nist.gov/5)
> - [#Gamma Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/GammaFunctions.html)

### Factorial Function

- `factorial(n::UInt)`:     `n!`, The factorial of non-negative integer n
- `factorial2(n::UInt)`:    `n!!`, Double factorial
- `factorialk(n::UInt, k::UInt)`: `n(!!...!)`, Multifactorial of n of order k > 0
- `binom(x::Real, y::Real)`:    `C(n, k)`, binomial coefficient

### Gamma Function

> - [DLMF: §5.2 Gamma Function](https://dlmf.nist.gov/5.2)
> - [DLMF: §5.15 Polygamma Functions](https://dlmf.nist.gov/5.15)
> - [Gamma Function - Wolfram MathWorld](https://mathworld.wolfram.com/GammaFunction.html)
> - [#Polygamma Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/PolygammaFunctions.html)

- `Γ(z)`, `≡ gamma(z::Complex)`:        gamma function
- `ln Γ(z)`, `≡ loggamma(z::Complex)`:  log gamma function
- `ln |Γ(z)|`, `≡ logabsgamma(z::Complex)`:     log abs gamma function
- `ψ(z)`, `≡ psi(z::Complex)`, `≡ digamma(z)`:  psi function
- `ψ'(z)`, `≡ trigamma(z::Complex)`:    trigamma function
- `polygamma(n::UInt, z::Complex)`:     polygamma functions

### Incomplete Gamma Function

> - [DLMF: §8.2 Incomplete Gamma and Related Functions](https://dlmf.nist.gov/8.2)
> - [Incomplete Gamma Function - Wolfram MathWorld](https://mathworld.wolfram.com/IncompleteGammaFunction.html)

- `γ(a, z)`:    (lower) incomplete gamma function
- `Γ(a, z)`:    (upper/complementary) incomplete gamma function
- `γ*(a, z)`:   Tricomi’s incomplete gamma function
- `P(a, z)`:    Normalized lower incomplete gamma function
- `Q(a, z)`:    Normalized upper incomplete gamma function

### Pochhammer Function

> - [DLMF: §5.2 Pochhammer’s Symbol](https://dlmf.nist.gov/5.2#iii)
> - [Pochhammer Symbol - Wolfram MathWorld](https://mathworld.wolfram.com/PochhammerSymbol.html)

- `poch(z, n)`, `(z)_n = Γ(z+n)/Γ(z)`, Pochhammer’s symbol (or shifted factorial)
- `poch1(z, n)`, `≡ (poch(z, n) - 1)/z`

### Beta Function

> - [DLMF: §5.12 Beta Function](https://dlmf.nist.gov/5.12)
> - [Beta Function - Wolfram MathWorld](https://mathworld.wolfram.com/BetaFunction.html)

- `beta(a, b)`:             Beta function
- `logbeta(a, b)`:          log beta function
- `logabsbeta(a, b)`:       log abs beta function
- `beta(a, b, z)`:          incomplete beta function
- `I(a, b, z)`, `betainc(a, b, z)`: regularized incomplete beta function
- `betaincinv(a, b, y)`:    Inverse of the regularized incomplete beta function


## Exponential and Trigonometric Integrals

> - [DLMF: Chapter 6 Exponential, Logarithmic, Sine, and Cosine Integrals](https://dlmf.nist.gov/6)
> - [Named Integrals - Wolfram MathWorld](https://mathworld.wolfram.com/topics/NamedIntegrals.html)

### Exponential Integral

> - [DLMF: §6.2 Exponential and Logarithmic Integrals](https://dlmf.nist.gov/6.2)
> - [DLMF: §8.19 Generalized Exponential Integral](https://dlmf.nist.gov/8.19)

- `E₁(z)`, `expint(z::Complex)`: (principal value of) exponential integral
- `Eν(z)`, `expint(ν::Complex, z::Complex)`: generalized exponential integral
- `eᶻEν(z)`, `expintx(ν::Complex, z::Complex)`: scaled (generalized) exponential integral
- `Ei(x)`, `expinti(x::Real)`: exponential integral
- `Li(x)`, `logint(x::Real)`, `{x>1}`: logarithmic integral

### Trigonometric Integral

- `Si(z)`, `sinint(z)`: sine integral function
- `Ci(z)`, `cosint(z)`: cosine integral function
- `Shi(z)`, `sinhint(z)`: hyperbolic sine integral function
- `Chi(z)`, `coshint(z)`: hyperbolic cosine integral function


## Error Functions

> - [DLMF: Chapter 7 Error Functions, Dawson’s and Fresnel Integrals](https://dlmf.nist.gov/7)
> - [Error and Exponential Integral Functions — Wolfram Language Documentation](https://reference.wolfram.com/language/guide/ErrorAndExponentialIntegralFunctions.html)

### Error Function

> - [DLMF: §7.2 Error Functions](https://dlmf.nist.gov/7.2#i)
> - [DLMF: §7.17 Inverse Error Functions](https://dlmf.nist.gov/7.17)
> - [Erf - Wolfram MathWorld](https://mathworld.wolfram.com/topics/Erf.html)

- `erf(z)`:         error function
- `erf(z0, z1)`:    generalized error function `erf(z1) - erf(z0)`
- `erfc(z)`:        complementary error function
- `erfcx(z)`:       Scaled complementary error function
- `ω(z)`, `faddeeva(z)`: Faddeeva function
- `erfi(z)`:        Imaginary error function
- `erfinv(z)`:      inverse error function
- `erfcinv(z)`:     inverse complementary error function

### Dawson Integral

> - [DLMF: §7.2 Dawson’s Integral](https://dlmf.nist.gov/7.2#ii)
> - [Dawson's Integral - Wolfram MathWorld](https://mathworld.wolfram.com/DawsonsIntegral.html)

- `dawson(z)`: Dawson integral

### Fresnel Integral

> - [DLMF: §7.2 Fresnel Integrals](https://dlmf.nist.gov/7.2#iii)
> - [Fresnel Integrals - Wolfram MathWorld](https://mathworld.wolfram.com/FresnelIntegrals.html)

- `C(z)`: Fresnel integral C(z)
- `S(z)`: Fresnel integral S(z)
- `f(z)`: Fresnel auxiliary function f(z)
- `g(z)`: Fresnel auxiliary function g(z)


## Airy Functions

> - [DLMF: §9.2 Airy’s Equation](https://dlmf.nist.gov/9.2#i)
> - [Airy Functions - Wolfram MathWorld](https://mathworld.wolfram.com/AiryFunctions.html)

- `Ai(z)`, `airyai(z)`:         Airy Ai function
- `Ai'(z)`, `airyaiprime(z)`:   derivative of Airy Ai function
- `Bi(z)`, `airybi(z)`:         Airy Bi function
- `Bi'(z)`, `airybiprime(z)`:   derivative of Airy Bi function
- `airyaix(z)`:         scaled Airy Ai function
- `airyaiprimex(z)`:    scaled derivative of Airy Ai function
- `airybix(z)`:         scaled Airy Bi function
- `airybiprimex(z)`:    scaled derivative of Airy Bi function

### Zeros of Airy Function

> - [DLMF: §9.9 Zeros ‣ Airy Functions](https://dlmf.nist.gov/9.9)

- `ai_zeros(nt)`:   Compute nt zeros and values of the Airy function Ai and its derivative
- `bi_zeros(nt)`:   Compute nt zeros and values of the Airy function Bi and its derivative

### Integral of Airy Function

> - [DLMF: §9.10 Integrals ‣ Airy Functions](https://dlmf.nist.gov/9.10)

- `airyaiint(z)`:   Integral of Airy Ai function
- `airybiint(z)`:   Integral of Airy Bi function

### Scorer Function

> - [DLMF: §9.12 Scorer Functions](https://dlmf.nist.gov/9.12)

- `Gi(z)`:  Scorer function Gi(z)
- `Gi'(z)`: derivative of the Scorer function Gi(z)
- `Hi(z)`:  Scorer function Hi(z)
- `Hi'(z)`: derivative of the Scorer function Hi(z)


## Bessel Functions

> - [DLMF: Chapter 10 Bessel Functions](https://dlmf.nist.gov/10)
> - [Bessel Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/BesselFunctions.html)

### Hankel Function

> - [DLMF: §10.2 Bessel Functions of the Third Kind (Hankel Functions)](https://dlmf.nist.gov/10.2#Px3)

### Modified Bessel Function

> - [DLMF: §10.25 Modified Bessel Functions](https://dlmf.nist.gov/10.25)

### Spherical Bessel Function

> - [DLMF: §10.47 Spherical Bessel Functions](https://dlmf.nist.gov/10.47)

### Kelvin Function

> - [DLMF: §10.61 Kelvin Functions](https://dlmf.nist.gov/10.61)
> - [Kelvin Functions - Wolfram MathWorld](https://mathworld.wolfram.com/KelvinFunctions.html)

- `kelvin(x)`:          Kelvin Function
- `kelvin_zeros(nt)`
- `ber(x), bei(x)`:     Kelvin Function of the First Kind
- `ker(x), kei(x)`:     Kelvin Function of the Second Kind
- `berp(x), beip(x)`:   Derivative of the Kelvin Function of the First Kind
- `kerp(x), keip(x)`:   Derivative of the Kelvin Function of the Second Kind


## Struve Functions

> - [DLMF: Chapter 11 Struve and Related Functions](https://dlmf.nist.gov/11)
> - [Struve Function - Wolfram MathWorld](https://mathworld.wolfram.com/StruveFunction.html)


## Parabolic Cylinder Functions

> - [DLMF: Chapter 12 Parabolic Cylinder Functions](https://dlmf.nist.gov/12)
> - [Parabolic Cylinder Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/ParabolicCylinderFunctions.html)


## Confluent Hypergeometric Functions

> - [DLMF: Chapter 13 Confluent Hypergeometric Functions](https://dlmf.nist.gov/13)
> - [Confluent Hypergeometric Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/ConfluentHypergeometricFunctions.html)


## Legendre Functions

> - [DLMF: Chapter 14 Legendre and Related Functions](https://dlmf.nist.gov/14)


## Hypergeometric Functions

> - [DLMF: Chapter 15 Hypergeometric Function](https://dlmf.nist.gov/15)
> - [Hypergeometric Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/HypergeometricFunctions.html)


## Generalized Hypergeometric Functions

> - [DLMF: Chapter 16 Generalized Hypergeometric Functions and Meijer 𝐺-Function](https://dlmf.nist.gov/16)
> - [Generalized Hypergeometric Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/GeneralizedHypergeometricFunctions.html)


## Elliptic Integrals

> - [DLMF: Chapter 19 Elliptic Integrals](https://dlmf.nist.gov/19)
> - [Elliptic Integrals - Wolfram MathWorld](https://mathworld.wolfram.com/topics/EllipticIntegrals.html)

### Legendre Integral

> - [Legendre’s Integrals - DLMF](https://dlmf.nist.gov/19.2.8)

- `K(m)`: `{m in (-Inf,1]}`, (Legendre’s) complete elliptic integral of the first kind
- `E(m)`: `{m in (-Inf,1]}`, (Legendre’s) complete elliptic integral of the second kind
- `Π(n, m)`:  (Legendre’s) complete elliptic integral of the third kind

### Symmetric Integral

> - [Symmetric Integrals - DLMF](https://dlmf.nist.gov/19.16)
> - [Carlson Elliptic Integrals - MathWorld](https://mathworld.wolfram.com/CarlsonEllipticIntegrals.html)

- `RF(x, y, z)`,    symmetric elliptic integral of first kind
- `RG(x, y, z)`,    symmetric elliptic integral of second kind
- `RJ(x, y, z, p)`, symmetric elliptic integral of third kind
- `RD(x, y, z)`,    `≡ RJ(x, y, z, z)`, elliptic integral symmetric in only two variables
- `RC(x, y)`,       `≡ RF(x, y, y)`,    Carlson’s combination of inverse circular and inverse hyperbolic functions


## Elliptic Functions

### Theta Function

> - [DLMF: §20.2 Theta Functions](https://dlmf.nist.gov/20.2)
> - [Theta Functions - Wolfram MathWorld](https://mathworld.wolfram.com/topics/ThetaFunctions.html)

- `jtheta(n, z, q)`:    (Jacobi) Theta function

### Jacobi Elliptic Function

> - [DLMF: §22.2 Jacobian Elliptic Functions](https://dlmf.nist.gov/22.2)
> - [Jacobi Elliptic Functions - Wolfram MathWorld](https://mathworld.wolfram.com/JacobiEllipticFunctions.html)

- `sn(z, k)`, Jacobi Elliptic Function sn
- `cn(z, k)`, Jacobi Elliptic Function cn
- `dn(z, k)`, Jacobi Elliptic Function dn

### Weierstrass Elliptic Function

> - [DLMF: §23.2 Weierstrass Elliptic Functions Functions](https://dlmf.nist.gov/23.2)
> - [Weierstrass Elliptic Function - Wolfram MathWorld](https://mathworld.wolfram.com/WeierstrassEllipticFunction.html)
> - [Weierstrass Zeta Function - Wolfram MathWorld](https://mathworld.wolfram.com/WeierstrassZetaFunction.html)
> - [Weierstrass Sigma Function - Wolfram MathWorld](https://mathworld.wolfram.com/WeierstrassSigmaFunction.html)

- `WeierstrassP(z)`,        Weierstrass ℘ function
- `WeierstrassZeta(z)`,     Weierstrass zeta function
- `WeierstrassSigma(x)`,    Weierstrass sigma function


## Zeta Functions

> - [DLMF: §25.2 Riemann Zeta Function](https://dlmf.nist.gov/25.2)
> - [Riemann Zeta Function - Wolfram MathWorld](https://mathworld.wolfram.com/RiemannZetaFunction.html)

### Dilogarithms

> - [DLMF: §25.12 Polylogarithms](https://dlmf.nist.gov/25.12)
> - [Spence's Function - Wolfram MathWorld](https://mathworld.wolfram.com/SpencesFunction.html)

### Dirichlet L-function

> - [DLMF: §25.15 Dirichlet 𝐿-functions](https://dlmf.nist.gov/25.15)
> - [Dirichlet L-Series - Wolfram MathWorld](https://mathworld.wolfram.com/DirichletL-Series.html)
> - [Dirichlet Eta Function - Wolfram MathWorld](https://mathworld.wolfram.com/DirichletEtaFunction.html)


## Mathieu Functions

> - [DLMF: §28.2 Mathieu Functions of Integer Order ‣ Eigenfunctions](https://dlmf.nist.gov/28.2#vi)
> - [Mathieu Function - Wolfram MathWorld](https://mathworld.wolfram.com/MathieuFunction.html)

### Eigenvalue

> - [DLMF: §28.2 Mathieu Functions of Integer Order ‣ Eigenvalues](https://dlmf.nist.gov/28.2#v)


## Spheroidal Wave Functions

> - [DLMF: Chapter 30 Spheroidal Wave Functions](https://dlmf.nist.gov/30)
> - [Prolate Spheroidal Wave Function - Wolfram MathWorld](https://mathworld.wolfram.com/ProlateSpheroidalWaveFunction.html)
> - [Oblate Spheroidal Wave Function - Wolfram MathWorld](https://mathworld.wolfram.com/OblateSpheroidalWaveFunction.html)

### Characteristic Value


## Miscellaneous Functions

> TODO
