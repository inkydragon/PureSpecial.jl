# Cephes

> CEPHES MATHEMATICAL FUNCTION LIBRARY

This computer software library is a collection of more than
400 high quality mathematical routines for scientific and
engineering applications.   All are written entirely in C
language.  Many of the functions are supplied in six different
arithmetic precisions: 32 bit single (24-bit significand), 64 bit
IEEE double (53-bit), 64 bit DEC (56-bit), 80 or 96 bit IEEE long
double (64-bit), and extended precision formats having 144-bit
and 336-bit significands.  The extended precision arithmetic is
included with the function library.

```
Copyright 1984 - 1992 by Stephen L. Moshier

Release 1.0: July, 1984
Release 1.1: March, 1985
Release 1.2: May, 1986
Release 2.0: April, 1987
Release 2.1: March, 1989
Release 2.2: July, 1992
```

## Double Precision


### Arithmetic and Algebraic
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Square root                     sqrt    2e-17   2e-16
Long integer square root        lsqrt   1       1
Cube root                       cbrt    2e-17   2e-16
Evaluate polynomial             polevl
Evaluate Chebyshev series       chbevl
Round to nearest integer value  round
Truncate upward to integer      ceil
Truncate downward to integer    floor
Extract exponent                frexp
Add integer to exponent         ldexp
Absolute value                  fabs
Rational arithmetic             euclid
Roots of a polynomial           polrt
Reversion of power series       revers
IEEE 854 arithmetic             ieee
Polynomial arithmetic (polyn.c):
  Add polynomials                 poladd
  Subtract polynomials            polsub
  Multiply polynomials            polmul
  Divide polynomials              poldiv
  Substitute polynomial variable  polsbt
  Evaluate polynomial             poleva
  Set all coefficients to zero    polclr
  Copy coefficients               polmov
  Display coefficients            polprt
 Note, polyr.c contains routines corresponding to
 the above for polynomials with rational coefficients.
Power series manipulations (polmisc.c):
  Square root of a polynomial     polsqt
  Arctangent                      polatn
  Sine                            polsin
Reversion of power series       revers
```

### Exponential and Trigonometric
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Arc cosine                      acos    3e-17   3e-16
Arc hyperbolic cosine           acosh   4e-17   5e-16
Arc hyperbolic sine             asinh   5e-17   4e-16
Arc hyperbolic tangent          atanh   3e-17   2e-16
Arcsine                         asin    6e-17   5e-16
Arctangent                      atan    4e-17   3e-16
Quadrant correct arctangent     atan2   4e-17   4e-16
Cosine                          cos     3e-17   2e-16
Cosine of arg in degrees        cosdg   4e-17   2e-16
Exponential, base e             exp     3e-17   2e-16
Exponential, base 2             exp2    2e-17   2e-16
Exponential, base 10            exp10   3e-17   2e-16
Hyperbolic cosine               cosh    3e-17   3e-16
Hyperbolic sine                 sinh    4e-17   3e-16
Hyperbolic tangent              tanh    3e-17   3e-16
Logarithm, base e               log     2e-17   2e-16
Logarithm, base 2               log2            2e-16
Logarithm, base 10              log10   3e-17   2e-16
Power                           pow     1e-15   2e-14
Integer Power                   powi		9e-14
Sine                            sin     3e-17   2e-16
Sine of arg in degrees          sindg   4e-17   2e-16
Tangent                         tan     4e-17   3e-16
Tangent of arg in degrees       tandg   3e-17   3e-16
```

### Exponential integral
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Exponential integral            expn    2e-16   2e-15
Hyperbolic cosine integral      shichi  9e-17   8e-16
Hyperbolic sine integral        shichi  9e-17   7e-16
Cosine integral                 sici    8e-17A  7e-16
Sine integral                   sici    4e-17A  4e-16
```

### Gamma
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Beta                            beta    8e-15   8e-14
Factorial                       fac     2e-17   2e-15
Gamma                           gamma   1e-16   1e-15
Logarithm of gamma function     lgam    5e-17   5e-16
Incomplete beta integral        incbet  4e-14   4e-13
Inverse beta integral           incbi   3e-13   8e-13
Incomplete gamma integral       igam    5e-15   4e-14
Complemented gamma integral     igamc   3e-15   1e-12
Inverse gamma integral          igami   9e-16   1e-14
Psi (digamma) function          psi     2e-16   1e-15
Reciprocal Gamma                rgamma  1e-16   1e-15
```

### Error function
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Error function                  erf     5e-17   4e-16
Complemented error function     erfc    5e-16   6e-14
Dawson's integral               dawsn   7e-16   7e-16
Fresnel integral (C)            fresnl  2e-16   2e-15
Fresnel integral (S)            fresnl  2e-16   2e-15
```

### Bessel
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Airy (Ai)                       airy    6e-16A  2e-15A
Airy (Ai')                      airy    6e-16A  5e-15A
Airy (Bi)                       airy    6e-16A  4e-15A
Airy (Bi')                      airy    6e-16A  5e-15A
Bessel, order 0                 j0      4e-17A  4e-16A
Bessel, order 1                 j1      4e-17A  3e-16A
Bessel, order n                 jn      7e-17A  2e-15A
Bessel, noninteger order        jv              5e-15A
Bessel, second kind, order 0    y0      7e-17A  1e-15A
Bessel, second kind, order 1    y1      9e-17A  1e-15A
Bessel, second kind, order n    yn      3e-16A  3e-15A
Bessel, noninteger order        yv      see struve.c
Modified Bessel, order 0        i0      8e-17   6e-16
Exponentially scaled i0         i0e     8e-17   5e-16
Modified Bessel, order 1        i1      1e-16   2e-15
Exponentially scaled i1         i1e     1e-16   2e-15
Modified Bessel, nonint. order  iv      3e-15   2e-14
Mod. Bessel, 3rd kind, order 0  k0      1e-16   1e-15
Exponentially scaled k0         k0e     1e-16   1e-15
Mod. Bessel, 3rd kind, order 1  k1      9e-17   1e-15
Exponentially scaled k1         k1e     9e-17   8e-16
Mod. Bessel, 3rd kind, order n  kn      1e-9    2e-8
```

### Hypergeometric
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Confluent hypergeometric        hyperg  1e-15   2e-14
Gauss hypergeometric function   hyp2f1  4e-11   9e-8
2F0                             hyp2f0f  see hyperg.c
1F2                             onef2f   see struve.c
3F0                             threef0f see struve.c
```

### Elliptic
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Complete elliptic integral (E)   ellpe  3e-17   2e-16
Incomplete elliptic integral (E) ellie  2e-16   2e-15
Complete elliptic integral (K)   ellpk  4e-17   3e-16
Incomplete elliptic integral (K) ellik  9e-17   6e-16
Jacobian elliptic function (sn)  ellpj  5e-16A  4e-15A
Jacobian elliptic function (cn)  ellpj          4e-15A
Jacobian elliptic function (dn)  ellpj          1e-12A
Jacobian elliptic function (phi) ellpj          9e-16
```

### Probability
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Binomial distribution           bdtr    4e-14   4e-13
Complemented binomial           bdtrc   4e-14   4e-13
Inverse binomial                bdtri   3e-13   8e-13
Chi square distribution         chdtr   5e-15   3e-14
Complemented Chi square         chdtrc  3e-15   2e-14
Inverse Chi square              chdtri  9e-16   6e-15
F distribution                  fdtr    4e-14   4e-13
Complemented F                  fdtrc   4e-14   4e-13
Inverse F distribution          fdtri   3e-13   8e-13
Gamma distribution              gdtr    5e-15   3e-14
Complemented gamma              gdtrc   3e-15   2e-14
Negative binomial distribution  nbdtr   4e-14   4e-13
Complemented negative binomial  nbdtrc  4e-14   4e-13
Normal distribution             ndtr    2e-15   3e-14
Inverse normal distribution     ndtri   1e-16   7e-16
Poisson distribution            pdtr    3e-15   2e-14
Complemented Poisson            pdtrc   5e-15   3e-14
Inverse Poisson distribution    pdtri   3e-15   5e-14
Student's t distribution        stdtr   2e-15   2e-14
```

### Miscellaneous
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Dilogarithm                     spence  3e-16   4e-15
Riemann Zeta function           zetac   1e-16   1e-15
Two argument zeta function      zeta
Struve function                 struve
```

### Matrix
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Fast Fourier transform          fftr
Simultaneous linear equations   simq
Simultaneous linear equations   gels (symmetric coefficient matrix)
Matrix inversion                minv
Matrix multiply                 mmmpy
Matrix times vector             mvmpy
Matrix transpose                mtransp
Eigenvectors (symmetric matrix) eigens
Levenberg-Marquardt nonlinear equations  lmdif
```

### Numerical Integration
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Simpson's rule                  simpsn
Runge-Kutta                     runge - see de118
Adams-Bashforth                 adams - see de118
```

### Complex Arithmetic
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Complex addition                cadd    1e-17   1e-16
Subtraction                     csub    1e-17   1e-16
Multiplication                  cmul    2e-17   2e-16
Division                        cdiv    5e-17   4e-16
Absolute value                  cabs    3e-17   3e-16
Square root                     csqrt   3e-17   3e-16
```

### Complex Exponential and Trigonometric
```
                                          Accuracy
Function                        Name    DEC     IEEE
--------                        ----    ----    ----
Exponential                     cexp    4e-17   3e-16
Logarithm                       clog    9e-17   5e-16A
Cosine                          ccos    5e-17   4e-16
Arc cosine                      cacos   2e-15   2e-14
Sine                            csin    5e-17   4e-16
Arc sine                        casin   2e-15   2e-14
Tangent                         ctan    7e-17   7e-16
Arc tangent                     catan   1e-16   2e-15
Cotangent                       ccot    7e-17   9e-16
```
