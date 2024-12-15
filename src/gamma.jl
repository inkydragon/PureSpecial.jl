# SPDX-License-Identifier: MIT

"""Gamma Functions

## Factorial Function
## Gamma Function
## Incomplete Gamma Function
## Pochhammer Function
## Beta Function
"""

# Factorial Function
public factorial, binomial  # Base.factorial, Base.binomial
export factorial2, factorialk
# Gamma Function
export gamma, loggamma, logabsgamma, digamma, trigamma, polygamma


"""
    factorial(n)

Compute the factorial function `n!`.
"""
function factorial end

"""
    factorial2(n)

Compute the double factorial function `n!!`.
"""
function factorial2 end

"""
    factorialk(n, k)

Compute the multifactorial function of `n` of order `k`, `n(!!...!)`.
"""
function factorialk end

"""
    binomial(n, k)

Compute the binomial coefficient, `C(n, k)`.
"""
function binomial end


"""
    gamma(z)

Compute the gamma function for complex `z`.

# Examples
```jldoctest
julia> [ (i, gamma(i)) for i in 1:5 ]
5-element Vector{Tuple{Int64, Float64}}:
 (1, 1.0)
 (2, 1.0)
 (3, 2.0)
 (4, 6.0)
 (5, 24.0)

julia> gamma(0.5)^2 ≈ π
true

julia> gamma(4 + 1) == prod(1:4) == factorial(4)
true
```

See also: `factorial`, `loggamma`.
"""
gamma(x::Float64) = Specfun.gamma2(x)
gamma(z::ComplexF64) = Specfun.cgama(x, 1)
gamma(i::Integer) = gamma(float(i))

"""
    loggamma(z)

Compute the principal branch of the logarithm of the gamma function
`log(gamma(z))`.
"""
function loggamma end

"""
    logabsgamma(z)

Compute the logarithm of the absolute value of the gamma function
`log(abs(gamma(z)))`
"""
function logabsgamma end

"""
    digamma(z)

Compute the digamma function (psi function).
"""
function digamma end

"""
    trigamma(z)

Compute the trigamma function.
"""
function trigamma end

"""
    polygamma(n, x)

Compute the polygamma function.
"""
function polygamma end
