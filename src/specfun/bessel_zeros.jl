# SPDX-License-Identifier: MIT OR BSD-3-Clause
#   See also: src/specfun/LICENSE.md

"""Scipy: scipy.special: Zeros of Bessel functions

## jdzo
jnjnp_zeros(nt)
    Compute zeros of integer-order Bessel functions Jn and Jn'.

## jyzo
jnyn_zeros(n, nt)
    Compute nt zeros of Bessel functions Jn(x), Jn'(x), Yn(x), and Yn'(x).

    jn_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel functions Jn.

    jnp_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel function derivatives Jn'.

    yn_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel function Yn(x).

    ynp_zeros(n, nt) -> jnyn_zeros
        Compute zeros of integer-order Bessel function derivatives Yn'(x).

## cyzo
y0_zeros(nt[, complex])
    Compute nt zeros of Bessel function Y0(z), and derivative at each zero.

y1_zeros(nt[, complex])
    Compute nt zeros of Bessel function Y1(z), and derivative at each zero.

y1p_zeros(nt[, complex])
    Compute nt zeros of Bessel derivative Y1'(z), and value at each zero.
"""


"""
    bjndd(x::Float64, n::Int)

Compute Bessel functions Jn(x) and their
first and second derivatives ( 0 <= n <= 100 )

## Input
- `x` ---  Argument of Jn(x)  ( x ≥ 0 )
- `n` ---  Order of Jn(x)

## Output
- `BJ(n+1)` ---  Jn(x)
- `DJ(n+1)` ---  Jn'(x)
- `FJ(n+1)` ---  Jn"(x)
"""
function bjndd(x::Float64, n::Int)
    @assert n in 0:100
    array_szie = 101
    bj = zeros(array_szie)
    dj = zeros(array_szie)
    fj = zeros(array_szie)

    m = 1
    for nt in 1:900
        # TODO: use a table?
        mt = 0.5 * log10(6.28 * nt) - nt * log10(1.36 * abs(x) / nt)
        if mt > 20
            m = nt
            break
        end
    end #= 10 =#

    #= 15 =#
    bs = 0.0
    f = 0.0
    f0 = 0.0
    f1 = 1e-35
    for k in m:-1:0
        f = 2.0 * (k + 1.0) * f1 / x - f0
        if k <= n
            bj[k+1] = f
        end
        if k % 2 == 0
            bs += 2.0 * f
        end
        f0 = f1
        f1 = f
    end #= 20 =#

    # TODO: use @.
    for k in 1:(n + 1)
        bj[k] /= bs - f
    end #= 25 =#

    dj[1] = -bj[2]
    fj[1] = -bj[1] - dj[1] / x

    for k in 1:n
        dj[k+1] = bj[k] - k * bj[k+1] / x
        fj[k+1] = (k*k / (x*x) - 1.0) * bj[k+1] - dj[k+1] / x
    end #= 30 =#

    return bj, dj, fj
end
