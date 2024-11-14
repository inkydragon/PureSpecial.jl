# SPDX-License-Identifier: MIT

"""Airy Functions

## Zeros of Airy Function
## Integral of Airy Function
## Scorer Function
"""

export ai_zeros, bi_zeros, itairy


"""
    ai_zeros(nt)

Compute the first `nt` zeros of Airy functions `Ai(x)` and `Ai'(x)`, `a` and `a'`,
and the associated values of `Ai(a')` and `Ai'(a)`.
"""
function ai_zeros(nt)
    if nt < 0
        throw(DomainError(nt, "`nt` must be a positive integer (nt >= 0)"))
    end

    # for Ai(x) and Ai'(x)
    kf = 1
    a, b, c, d = zeros(nt), zeros(nt), zeros(nt), zeros(nt)

    Specfun.airyzo!(nt, kf, a, b, c, d)

    a, b, c, d
end

"""
    bi_zeros(nt)

Compute the first `nt` zeros of Airy functions `Bi(x)` and `Bi'(x)`, `b` and `b'`,
and the associated values of `Bi(b')` and `Bi'(b)`.
"""
function bi_zeros(nt)
    if nt < 0
        throw(DomainError(nt, "`nt` must be a positive integer (nt >= 0)"))
    end

    # for Bi(x) and Bi'(x)
    kf = 2
    a, b, c, d = zeros(nt), zeros(nt), zeros(nt), zeros(nt)

    Specfun.airyzo!(nt, kf, a, b, c, d)

    a, b, c, d
end

"""
    itairy(x)

Compute the integrals of Airy functions with respect to t from 0 and x ( x ≥ 0 )
"""
function itairy(x)
    if x < 0
        throw(DomainError(x, "`x` must be a positive number (x >= 0)"))
    end

    Specfun.itairy(x)
end
