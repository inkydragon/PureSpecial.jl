# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md
"""Bernoulli and Euler Numbers"""
#=
BERNOA, 5
BERNOB, 6
    EULERA, 8
    EULERB, 9
=#

function bernoa(n::Int, bn::Vector{Float64})
    bn[1] = 1.0
    bn[2] = -0.5

    for m = 2:n
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
        bn[m+1] = 0.0
    end

    return bn
end

"""
Compute Bernoulli number Bn

Input :
- n --- Serial number, n >= 2

Output:
- BN(n) --- Bn, n = 0..N
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
Compute Bernoulli number Bn

Input :
- n --- Serial number, n >= 2

Output:
- BN(n) --- Bn, n = 0..N
"""
function bernob(n::Int)
    bn = zeros(Float64, n+1)
    bernob(n, bn)
    return bn
end
