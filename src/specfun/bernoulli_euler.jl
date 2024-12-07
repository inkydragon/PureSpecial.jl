# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Bernoulli and Euler Numbers"""
#=
- ✅ bernoa
- ✅ bernob
- ✅ eulera
- ✅ eulerb
=#

function bernoa(n::Int, bn::Vector{Float64})
    bn[1] = 1.0
    bn[2] = -0.5

    for m = 2:n
        # CoSF (1.1.5)
        s = -(1.0/(m+1) - 0.5)
        for k = 2:(m-1)
            r = 1.0
            for j = 2:k
                r *= (j + m - k) / j
            end
            s -= r * bn[k+1]
        end
        bn[m+1] = s
    end

    for m = 3:2:n
        # Set B(2n+1) = 0, n=1,2,3,...
        bn[m+1] = 0.0
    end

    return bn
end

"""
Compute Bernoulli number `Bn`, using recurrence relation.
Not work for large `n`.

Input :
- `n` --- Serial number, `n >= 2`

Output:
- `BN(n)` --- Bn, n = 0..N
"""
function bernoa(n::Int)
    bn = zeros(Float64, n+1)
    bernoa(n, bn)
    return bn
end

function bernob(n::Int, bn::Vector{Float64})
    @assert (n+1) >= 3
    @assert length(bn) >= (n+1)
    tpi = 2*pi
    @assert isequal(6.283185307179586, tpi)

    bn[1] = 1.0
    bn[2] = -0.5
    bn[3] = 1.0 / 6.0

    # DLMF 24.8.1, with x = 0
    # CoSF (1.1.16)
    r1 = (2.0 / tpi)^2
    for m in 4:2:n
        r1 = -r1 * (m - 1) * m / (tpi * tpi)
        r2 = 1.0
        for k in 2:10_000
            s = (1.0 / k)^m
            r2 += s
            if s < SF_EPS15
                break
            end
        end

        bn[m+1] = r1 * r2
    end

    return bn
end

"""
Compute Bernoulli number `Bn`,
using series expression (DLMF 24.8.1).

Input :
- `n` --- Serial number, `n >= 2`

Output:
- `BN(n)` --- Bn, n = 0..N
"""
function bernob(n::Int)
    bn = zeros(Float64, n+1)
    bernob(n, bn)
    return bn
end

function eulera(n::Int, en::Vector{Float64})
    en[1] = 1.0

    # CoSF (1.2.5)
    half_n = trunc(Int, n / 2)
    for m in 1:half_n
        s = 1.0
        for k in 1:(m-1)
            r = 1.0
            for j in 1:(2*k)
                r = r * (2.0*m - 2.0*k + j) / j
            end
            s += r * en[2*k + 1]
        end

        en[2*m + 1] = -s
    end

    return en
end

"""
Compute Euler number `En`, using recurrence relation.
Not work for large `n`.

Input :
- `n` --- Serial number, `n >= 2`

Output:
- `EN(n)` --- En, n = 0..N
"""
function eulera(n::Int)
    en = zeros(Float64, n+1)
    return eulera(n, en)
end

function eulerb(n::Int, en::Vector{Float64})
    @assert (n+1) >= 3
    @assert length(en) >= (n+1)

    en[1] = 1.0
    en[3] = -1.0

    # DLMF 24.8.4, with x = pi/2
    # CoSF (1.2.13)
    hpi = 2.0 / pi
    r1 = -4.0 * hpi^3
    for m in 4:2:n
        r1 = -r1 * (m - 1) * m * hpi * hpi
        r2 = 1.0
        isgn = 1.0
        for k in 3:2:1000
            isgn = -isgn
            s = (1.0 / k)^(m + 1)
            r2 += isgn * s
            if s < SF_EPS15
                break
            end
        end

        en[m+1] = r1 * r2
    end

    return en
end

"""
Compute Euler number `En`,
using series expression (DLMF 24.8.4).

Input :
- `n` --- Serial number, `n >= 2`

Output:
- `EN(n)` --- En, n = 0..N
"""
function eulerb(n::Int)
    en = zeros(Float64, n+1)
    return eulerb(n, en)
end
