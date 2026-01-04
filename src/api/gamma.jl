# SPDX-License-Identifier: MIT
"""
# Gamma Functions

## Factorial Function
- factorial, doublefactorial, multifactorial, superfactorial, hyperfactorial

### Binomial Coefficient
- binomial, multinomial

## Gamma Function
## Incomplete Gamma Function
## Pochhammer Function
## Beta Function

# Reference
"""
const GammaFunctionsDoc = nothing


"""
    factorial(n::Integer)

Factorial of non-negative integer `n`, `n! = 1 * 2 * 3 * ... * n`.

# Reference
- [DLMF: factorial - Common Notations and Definitions](https://dlmf.nist.gov/front/introduction#common.t1.r15)
- [Factorial -- from Wolfram MathWorld](https://mathworld.wolfram.com/Factorial.html)
"""
function factorial end

"""
    doublefactorial(n::Integer)

Double factorial of non-negative integer `n`, `n!!`.

Also known as:  `factorial2(n)`.

# Reference
- [DLMF: double factorial - Common Notations and Definitions](https://dlmf.nist.gov/front/introduction#common.t1.r16)
- [Double Factorial -- from Wolfram MathWorld](https://mathworld.wolfram.com/DoubleFactorial.html)
"""
function doublefactorial end

"""
    multifactorial(k::Integer, n::Integer)

Multifactorial of non-negative integer `n` of order `k`, `n(!!...!)`.

Also known as:  `factorialk(k, n)`.

# Reference
- [Multifactorial -- from Wolfram MathWorld](https://mathworld.wolfram.com/Multifactorial.html)
- [Multifactorials Definitions - Wikipedia](https://en.wikipedia.org/wiki/Double_factorial#Generalizations)
"""
function multifactorial end

"""
    superfactorial(n::Integer)

Superfactorial of non-negative integer `n`, `n\$ = 1! * 2! * 3! * ... * n!`.

# Reference
- [Superfactorial -- from Wolfram MathWorld](https://mathworld.wolfram.com/Superfactorial.html)
- [Superfactorial - Wikipedia](https://en.wikipedia.org/wiki/Superfactorial)
"""
function superfactorial end

"""
    hyperfactorial(n::Integer)

Hyperfactorial of non-negative integer `n`, `H(n) = 1^1 * 2^2 * 3^3 * ... * n^n`.

# Reference
- [Hyperfactorial -- from Wolfram MathWorld](https://mathworld.wolfram.com/Hyperfactorial.html)
- [Hyperfactorial - Wikipedia](https://en.wikipedia.org/wiki/Hyperfactorial)
"""
function hyperfactorial end

"""
    binomial(n::Integer, k::Integer)
    binomial(z::Complex, k::Integer)

Binomial coefficient, `C(n, k)`.

# Reference
- [DLMF: §1.2(i) Binomial Coefficients](https://dlmf.nist.gov/1.2#i)
- [Binomial Coefficient -- from Wolfram MathWorld](https://mathworld.wolfram.com/BinomialCoefficient.html)
"""
function binomial end

"""
    multinomial(n1::Integer, ...)

Multinomial coefficient, `C(n₁, n₂, ..., nₖ)`.

# Reference
- [DLMF: §26.4(i) multinomial coefficient Definitions](https://dlmf.nist.gov/26.4#i.p1)
- [Multinomial Coefficient -- from Wolfram MathWorld](https://mathworld.wolfram.com/MultinomialCoefficient.html)
"""
function multinomial end
