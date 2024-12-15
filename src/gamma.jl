# SPDX-License-Identifier: MIT

"""Gamma Functions

## Factorial Function
## Gamma Function
## Incomplete Gamma Function
## Pochhammer Function
## Beta Function
"""

export gamma

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
